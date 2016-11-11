context("inputs.R")

test_that("Data loading", {
  expect(is.null(load_gff('a', 'b')), "load_gff is stub")
  expect(is.null(load_genome('a', 'b')), "load_genome is stub")
  expect(is.null(load_proteome('a', 'b')), "load_proteome is stub")
  expect(is.null(load_synmap('a', 'b', 'c')), "load_synmap is stub")
  expect(is.null(load_species_tree('a')), "load_species_tree is stub")
})
