context("derive.R")

fna <- c(
  #           start         intron   exon   stop
  paste0("A", "ATG", "CCC", "GATA", "GGG", "TAA", "T", collapse=""),
  #            stop         intron   exon   start 
  paste0("A", "TTA", "GGG", "CACA", "CCC", "CAT", "T", collapse="")
) %>%
  Biostrings::DNAStringSet() %>%
  magrittr::set_names(c("foo", "bar"))

gff <- GenomicRanges::GRanges(
    c("foo", "foo", "bar", "bar"), 
    ranges=IRanges::IRanges(start=c(2,12,2,12), width=c(6,6,6,6)),
    strand=c('+','+','-','-')
  )
GenomicRanges::mcols(gff)$type = "CDS"
GenomicRanges::mcols(gff)$parent = c("pan", "pan", "bob", "bob")

cds <- c("ATGGGGCCCTAA", "ATGCCCGGGTAA") %>%
  Biostrings::DNAStringSet() %>%
  magrittr::set_names(c("bob", "pan"))

faa <- c("MGP*", "MPG*") %>%
  Biostrings::AAStringSet() %>%
  magrittr::set_names(c("bob", "pan"))

test_that("mergeSeqs correctly merges", {
  expect_true(all( mergeSeqs(fna, gff, "CDS") == cds))
  expect_true(all( derive_aa(fna, gff) == faa))
})

test_that("dnaregex gracefully handles no matches", {
  expect_silent(dnaregex(fna, "GANDALF"))
  expect_true( length(fna[dnaregex(fna, "GANDALF")]) == 0 )
})

test_that("dnaregex works on all strand variants", {
  expect_equal(
    fna[dnaregex(fna, "ATG", strand='b')] %>% as.character %>% as.vector,
    c("ATG", "CAT")
  )
  expect_equal(
    fna[dnaregex(fna, "ATG", strand='p')] %>% as.character %>% as.vector,
    "ATG"
  )
  expect_equal(
    fna[dnaregex(fna, "ATG", strand='m')] %>% as.character %>% as.vector,
    "CAT"
  )
})
