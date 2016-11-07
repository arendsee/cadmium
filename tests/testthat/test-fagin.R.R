context("fagin.R")

test_that("Data loading", {
  expect(is.null(load_gff('a', 'b')), "GFF loading function failed")
  expect(is.null(load_genome('a', 'b')), "GFF loading function failed")
  expect(is.null(load_synmap('a', 'b', 'c')), "GFF loading function failed")
})
