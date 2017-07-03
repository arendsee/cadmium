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
  # TODO: add tests for this 

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

  dnaregex(dna, "N+", strand='u', ignore.case=TRUE)

} 

#' Derive mono-exonic ORF intervals from a genome
#'
#' Note, here I am doing the stupidest thing possible: just finding all
#' '<START>(...)+<STOP>' intervals. I should replace this with a real ORF
#' finder.
#'
#' @param dna DNAStringSet genome
#' @return GenomicRanges object
derive_orfgff <- function(dna) {

  # BUG: This misses potentially longer ORFs that start within a previous ORF
  # TODO: replace with real ORF finder

  orfpat <- "ATG(AAA|AAC|AAG|AAT|ACA|ACC|ACG|ACT|AGA|AGC|AGG|AGT|ATA|ATC|ATG|ATT|CAA|CAC|CAG|CAT|CCA|CCC|CCG|CCT|CGA|CGC|CGG|CGT|CTA|CTC|CTG|CTT|GAA|GAC|GAG|GAT|GCA|GCC|GCG|GCT|GGA|GGC|GGG|GGT|GTA|GTC|GTG|GTT|TAC|TAT|TCA|TCC|TCG|TCT|TGC|TGG|TGT|TTA|TTC|TTG|TTT){49,}(TAA|TGA|TAG)"

  dnaregex(dna, orfpat, strand='b', perl=TRUE, ignore.case=TRUE)

}

#' Extract DNA intervals, reverse complementing them if they are minus sense
#'
#' @param dna DNA template
#' @param gff Intervals to extract
#' @return DNAStringSet
extractWithComplements <- function(dna, gff){
  append(
      dna[gff[GenomicRanges::strand(gff) == '+']],
      dna[gff[GenomicRanges::strand(gff) == '-']] %>% Biostrings::reverseComplement() 
    )
}

#' Derive ORF protein sequences from genome and ORF ranges
#'
#' @param dna DNAStringSet genome
#' @param orfgff GenomicRanges object
#' @return AAStringSet object
derive_orffaa <- function(dna, orfgff) {

  # TODO: get a real ORF finder and sync this accordingly

  extractWithComplements(dna, orfgff) %>% Biostrings::translate()

}

mergeSeqs <- function(fna, gff, tag){
  g <- gff[gff$type == tag]

  # sort the elements by start
  # this is required, since they will be concatenated by order
  g <- g[order(GenomicRanges::start(g))]

  revpar <- g[GenomicRanges::strand(g) == '-']$parent %>% unique

  parents <- GenomicRanges::mcols(g)$parent

  # TODO: assert no elements within a group overlap

  seqs <- fna[g]                %>%
    base::split(parents)        %>%
    lapply(paste0, collapse="") %>%
    unlist                      %>%
    Biostrings::DNAStringSet()

  seqs[names(seqs) %in% revpar] <-
    Biostrings::reverseComplement(seqs[names(seqs) %in% revpar])

  seqs
}

#' Derive protein sequence from genome and gene model GFF
#'
#' @param dna DNAStringSet genome
#' @param gff GenomicRanges object annotated with \code{type} and \code{parent}
#' @return AAStringSet object
derive_aa <- function(dna, gff) {

  mergeSeqs(dna, gff, "CDS") %>% Biostrings::translate()

}

#' Derive mRNA sequence from genome and gene model GFF
#'
#' @param dna DNAStringSet genome
#' @param gff GenomicRanges object annotated with \code{type} and \code{parent}
#' @return DNAStringSet object
derive_trans <- function(dna, gff) {

  mergeSeqs(dna, gff, "exon")

}

#' Derive ORF intervals from a transcriptome
#'
#' @param trans DNAStringSet genome
#' @return GenomicRanges object
derive_transorfgff <- function(trans) {

  derive_orfgff(trans)

}

#' Derive ORF intervals from a transcriptome
#'
#' @param trans DNAStringSet genome
#' @param transorfgff GRanges locations of ORFs on mRNAs
#' @return GenomicRanges object
derive_transorffaa <- function(trans, transorfgff) {

  extractWithComplements(trans, transorfgff) %>% Biostrings::translate()

}
