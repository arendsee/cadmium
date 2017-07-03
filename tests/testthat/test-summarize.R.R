context("summarize.R")

dna <- load_dna("tiny/tiny.fna")
gff <- load_gff("tiny/tiny.gff")

empty_irange <- IRanges::IRanges()
empty_grange <- GenomicRanges::GRanges(seqnames=c(), ranges=empty_irange)
empty_faa <- Biostrings::AAStringSet()
empty_fna <- Biostrings::DNAStringSet()

test_that("Nothing can be summarized", {
  expect_silent(summarize_dna(empty_fna))
  expect_silent(summarize_faa(empty_faa))
  expect_silent(summarize_gff(empty_grange))
  expect_silent(summarize_granges(empty_grange))
  expect_silent(summarize_syn(data.frame()))
  expect_silent(summarize_nstring(empty_irange))
  expect_silent(summarize_nstring(empty_grange))
})
