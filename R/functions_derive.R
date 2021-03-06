#' Internal function for extracting intervals from pattern matches on DNA
#'
#' The \code{strand} parameter can be
#'
#' \itemize{
#'   \item 'p' search plus strand
#'   \item 'm' search minus strand
#'   \item 'u' search plus strand but report no strand
#'   \item 'b' search both strands
#' }
#'
#' @param x       DNAStringSet
#' @param pattern character regular expression
#' @param strand  which strand to search (see details)
#' @param ...     extra arguments for gregexpr
#' @return GRanges object
dnaregex <- function(x, pattern, strand=c('b'), ...){

  sinfo <- GenomeInfoDb::Seqinfo(seqnames=names(x), seqlengths=Biostrings::width(x))

  if(! strand %in% c('p', 'm', 'u', 'b')){
    stop("IllegalArgument: strand must be 'p', 'm', 'u', or 'b'")
  }

  if(strand %in% c('p', 'u', 'b')){
    f <- .dnaregex_unstranded(x, pattern, sinfo=sinfo, ...)
    if(strand != 'u'){
      GenomicRanges::strand(f) <- '+'
    }
  }

  if(strand %in% c('m', 'b')){
    r <- .dnaregex_unstranded(Biostrings::reverseComplement(x), pattern, sinfo=sinfo, ...)
    r <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(r),
      seqinfo  = sinfo,
      ranges   = IRanges::IRanges(
                  start = GenomeInfoDb::seqlengths(r)[as.character(GenomeInfoDb::seqnames(r))] - IRanges::end(r) + 1,
                  width = IRanges::width(r)
                 )
    )
    GenomicRanges::strand(r) <- '-'
  }

  switch(strand,
    p = f,
    m = r,
    u = f,
    b = c(f,r)
  )
  
}
.dnaregex_unstranded <- function(x, pattern, sinfo, ...){
  # gregexpr search for a pattern, returning the start and width on each given string

  gregexpr(pattern, as.character(x), ...) %>%
  magrittr::set_names(names(x)) %>%
  {
    data.frame(
      names = lapply(., length) %>% rep.int(names(.), times=.),
      start = unlist(.) %>% as.vector,
      width = lapply(., attr, "match.length") %>% unlist %>% as.vector
    )
  } %>%
  # gregexpr stores absence of a match as -1 in both start and width
  dplyr::filter(.data$start != -1) %>%
  {
    GenomicRanges::GRanges(
      seqnames = .$names,
      seqinfo  = sinfo, 
      ranges   = IRanges::IRanges(start=.$start, width=.$width)
    )
  }
}

#' Derive N-interval table from a genome
#'
#' @param dna DNAStringSet genome
#' @return GRanges object
#' @export
derive_nstring <- function(dna) {

  "Given a genome, find all sequences of N (unknown nucleotides)"

  dnaregex(dna, "N+", strand='u', ignore.case=TRUE) %>% synder::as_gff()

} 

extract_range <- function(dna, rng){
  d <- dna[rng]
  negative <- GenomicRanges::strand(rng) == '-'
  if(any(negative)){
    d[negative] <- Biostrings::reverseComplement(d[negative])
  }
  d
}

# This is needed since ORFik returns indices rather than names
fix_names <- function(g, dna){
  GenomeInfoDb::renameSeqlevels(
    g,
    names(dna)[as.integer(GenomeInfoDb::seqlevels(g))]
  )
}

getORFs <- function(dna, con, sense='.'){
  .findORFs <- function(x, strnd){
    orf <- ORFik::findORFs(
      x,
      startCodon    = con@orf@start,
      stopCodon     = con@orf@stop,
      minimumLength = con@orf@minlen,
      longestORF    = TRUE
    ) %>% as("GRanges")
    GenomicRanges::strand(orf) <- strnd
    fix_names(orf, dna)
  }
  plus_orfs <- if(sense == '.' || sense == '+'){
    .findORFs(dna, "+")
  } else {
    GenomicRanges::GRanges()
  }
  minus_orfs <- if(sense == '.' || sense == '-'){
    r <- .findORFs(Biostrings::reverseComplement(dna), "-")
    if(length(r) > 0){
      lengths <- Biostrings::width(dna[GenomeInfoDb::seqnames(r)])
      new_start <- lengths - GenomicRanges::end(r) + 1
      r@ranges@start <- as.integer(new_start)
    }
    r
  } else {
    GenomicRanges::GRanges()
  }
  orfs <- append(plus_orfs, minus_orfs)
  GenomicRanges::mcols(orfs) <- data.frame(
      seqid = paste0('orf_', seq_along(orfs)),
      type  = 'orf',
      stringsAsFactors=FALSE
  )
  orfs
}

#' Derive mono-exonic ORF intervals from a genome
#'
#' @param x DNAStringSet genome
#' @param con fagin config object
#' @return GenomicRanges object
#' @export
derive_genomic_ORFs <- function(dna, con) {

  "
  Given a genome, find the locations of ORFs on both strands.

  Note, here I am doing the stupidest thing possible: just finding all
  `<START>(...)+<STOP>` intervals. I should replace this with a real ORF
  finder.

  Gives ORFs unique ids (just an integer sequence for now) in meta-column `seqid`.
  "

  getORFs(dna, con, sense='.')
}

#' Derive ORFs from transcripts
#'
#' @param x DNAStringSet transcriptome
#' @param con fagin config object
#' @return GenomicRanges object
#' @export
derive_transcript_ORFs <- function(dna, con) {

  "
  Given a transcriptome, find the locations of all ORFs. Also gives ORFs unique
  ids (just an integer sequence for now) in meta-column `seqid`.
  "

  getORFs(dna, con, sense='+')
}

#' Extract DNA intervals, reverse complementing them if they are minus sense
#'
#' @param dna DNA template
#' @param gff Intervals to extract
#' @export
#' @return DNAStringSet
extract_with_complements <- function(dna, gff){

  "
  Extract a set of sequences from a genome given a GFF. Assumptions:
  
  1) All scaffolds named in the GFF are present in genome file

  2) All end positions in the GFF are less than or equal to the length of their
     respective scaffold (not currently tested).
  "

  dna <- Rsamtools::getSeq(dna, gff)

  names(dna) <- if(! is.null(gff$attr)){
    gff$attr 
  } else if(!is.null(GenomicRanges::mcols(gff)$seqid)) {
    GenomicRanges::mcols(gff)$seqid
  } else {
    NULL
  }

  dna

}
