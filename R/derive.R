#' Internal function for extracting intervals from pattern matches on DNA
#'
#' @param x       DNAStringSet
#' @param pattern character regular expression
#' @param ...     extra arguments for gregexpr
#' @return GRanges object
dnaregex <- function(x, pattern, ...){
  # TODO: add tests for this 

  # gregexpr search for a pattern, returning the start and width on each given string
  lapply(x, function(s) gregexpr(pattern, as.character(s), ...)) %>%
  {
    data.frame(
      names = lapply(., '[[', 1) %>% lapply(length) %>% rep.int(names(.), times=.),
      start = unlist(.) %>% as.vector,
      width = lapply(., '[[', 1) %>% lapply(attr, "match.length") %>% unlist %>% as.vector
    )
  } %>%
  # gregexpr stores absence of a match as -1 in both start and width
  dplyr::filter(.data$start != -1) %>%
  {
    GenomicRanges::GRanges(
      seqnames = .$names,
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

  dnaregex(dna, "N+", ignore.case=TRUE)

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

  orfpat <- "ATG(...)+?(TAA|TGA|TAG)"

  f <- dnaregex(dna, orfpat, perl=TRUE, ignore.case=TRUE)
  r <- dnaregex(Biostrings::reverseComplement(dna), orfpat, perl=TRUE, ignore.case=TRUE)
  GenomicRanges::strand(f) <- '+'
  GenomicRanges::strand(r) <- '-'

  c(f,r)
}

#' Derive ORF protein sequences from genome and ORF ranges
#'
#' @param dna DNAStringSet genome
#' @param orfgff GenomicRanges object
#' @return AAStringSet object
derive_orffaa <- function(dna, orfgff) {

}

#' Derive protein sequence from genome and gene model GFF
#'
#' @param dna DNAStringSet genome
#' @param gff GenomicRanges object annotated with \code{type} and \code{parent}
#' @return AAStringSet object
derive_aa <- function(gff, dna) { }

#' Derive mRNA sequence from genome and gene model GFF
#'
#' @param dna DNAStringSet genome
#' @param gff GenomicRanges object annotated with \code{type} and \code{parent}
#' @return DNAStringSet object
derive_trans <- function(gff, dna) { }

#' Derive ORF intervals from a transcriptome
#'
#' @param trans DNAStringSet genome
#' @return GenomicRanges object
derive_transorf <- function(trans) { }
