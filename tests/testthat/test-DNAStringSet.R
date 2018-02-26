context("test `DNAStringSet -> a` functions")

cds_sample <- Biostrings::DNAStringSet(c(
  'ATGTTTTAA',
  'ATG',
  'TAA',
  'TAATTT',
  'NNN'
))

cds_sample_ambiguous <- Biostrings::DNAStringSet('CCN')

test_that("translate works", {
  expect_equal(
    translate(cds_sample), Biostrings::AAStringSet(c("MF*", "M", "*", "*F", "X"))
  )
  expect_equal(
    translate(cds_sample_ambiguous), Biostrings::AAStringSet(c("P"))
  )
  expect_warning(
    Biostrings::DNAStringSet(c(x="A")) %>>% translate(species_name='unicorn')
  )
})
