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

  dnaregex(dna, "N+", strand='u', ignore.case=TRUE)

} 

#' Derive mono-exonic ORF intervals from a genome
#'
#' @param dna DNAStringSet genome
#' @return GenomicRanges object
#' @export
derive_orfgff <- function(dna) {

  "
  Given a genome, find the locations of ORFs on both strands.

  Note, here I am doing the stupidest thing possible: just finding all
  `<START>(...)+<STOP>` intervals. I should replace this with a real ORF
  finder.

  Gives ORFs unique ids (just an integer sequence for now) in meta-column `seqid`.

  BUG: This misses potentially longer ORFs that start within a previous ORF.
  TODO: replace with real ORF finder.
  "

  orfpat <- "ATG(AAA|AAC|AAG|AAT|ACA|ACC|ACG|ACT|AGA|AGC|AGG|AGT|ATA|ATC|ATG|ATT|CAA|CAC|CAG|CAT|CCA|CCC|CCG|CCT|CGA|CGC|CGG|CGT|CTA|CTC|CTG|CTT|GAA|GAC|GAG|GAT|GCA|GCC|GCG|GCT|GGA|GGC|GGG|GGT|GTA|GTC|GTG|GTT|TAC|TAT|TCA|TCC|TCG|TCT|TGC|TGG|TGT|TTA|TTC|TTG|TTT){49,}(TAA|TGA|TAG)"

  dnaregex(dna, orfpat, strand='b', perl=TRUE, ignore.case=TRUE) %>%
  {

    GenomicRanges::mcols(.) <- data.frame(
      seqid = paste0('orf_', seq_along(.)),
      type  = 'orf',
      stringsAsFactors=FALSE
    )

    .
  }

}

#' Extract DNA intervals, reverse complementing them if they are minus sense
#'
#' @param dna DNA template
#' @param gff Intervals to extract
#' @export
#' @return DNAStringSet
extractWithComplements <- function(dna, gff){

  "
  Extract a set of sequences from a genome given a GFF. Assumptions:
  
  1) All scaffolds named in the GFF are present in genome file

  2) All end positions in the GFF are less than or equal to the length of their
     respective scaffold (not currently tested).
  "

  dna <- Rsamtools::getSeq(dna, gff)

  if(!is.null(GenomicRanges::mcols(gff)$seqid))
    names(dna) <- GenomicRanges::mcols(gff)$seqid

  dna

}

# mergeSeqs <- function(dna, gff, tag){
#
#   "Extract all entries from `gff` that have type `tag` (for example, CDS or
#   exon). Then extract these matching intervals from `dna`. All extracted
#   sequences that share a common parent, are then merged. If they are on the
#   negative strang, the reverse complement is taken."
#
#   g <- gff %>_% {
#
#     "Assert that the GFF file has the required columns"
#
#     # TODO: need to standardize the case across columns
#
#     if(is.null(.$type) || is.null(.$Parent))
#       stop("GFF must have the meta-columns 'type' and 'Parent'")
#
#   } %>>% {
#
#     "Extract elements of type `tag` from the `gff`"
#
#     .[.$type == tag]
#
#   } %>>% {
#
#     "Sort the elements by start this is required, since they will be
#     concatenated by order"
#
#     .[order(GenomicRanges::start(.))]
#
#   }
#
#   revpar_ <- g %>>% {
#
#     "Get all unique parents on the reverse strand"
#
#     .[GenomicRanges::strand(.) == '-']$Parent %>% unique
#
#   }
#
#   forpar_ <- g %>>% { GenomicRanges::mcols(.)$Parent }
#
#   # TODO: assert no elements within a group overlap
#
#   g %>>% {
#
#     "Extract DNA sequences based on the GFF"
#
#     # this is slow -- 5s
#     dna[.]
#
#   } %>% rmonad::funnel(forpar=forpar_) %*>% {
#
#     "Collapse the sequences together, aggregating on the Parent"
#
#     q <- .
#     base::split(seq_len(length(forpar)), forpar) %>%
#     sapply(function(ids) paste(q[ids], collapse=""))
#
#   } %>>% {
#
#     "Convert list of strings back into a DNAStringSet object"
#
#     Biostrings::DNAStringSet(.)
#
#   } %>% rmonad::funnel(revpar=revpar_) %*>% {
#
#     "Take the reverse complement of all elements on the reverse strand"
#
#     .[names(.) %in% revpar] <-
#       Biostrings::reverseComplement(.[names(.) %in% revpar])
#
#     .
#
#   }
#
# }
