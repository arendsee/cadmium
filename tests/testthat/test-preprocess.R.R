context("preprocess.R")

dna <- load_dna("asdf/asdf.fna")
gff <- load_gff("asdf/asdf.gff")

test_that("load_dna works", {
  expect_true(all(Biostrings::width(dna)[1] == 60))
  # assert that description was removed from the header
  expect_equal( names(dna)[1], "foo" )
})

test_that("nstrings gets the right intervals", {
  expect_equal(derive_nstring(dna) %>% Biostrings::width(), 10)
  expect_equal(derive_nstring(dna) %>% Biostrings::start(), 51)
})

test_that("load tree works", {
  expect_equal({ t <- load_tree("asdf/asdf.tree"); tail(t$tip.label, 1) }, "human")
})

test_that("protein models are correctly assembled", {
  expect_equal(
    derive_aa(dna, gff) %>% as.character %>% as.vector,
    c(
      "MFGK*",     # sample two-exon plus-sense model
      "MFFFKPPP*", # reverse complement case
      "MFFFKPPP*"  # as above, but with GFF order scrambled
     )
  ) 
})


# test_that("derive_orfgff creates valid ORFs", {
#   # Here I test only that the ouput are ORFs
#   # Currently the way I find ORFs is totally wrong
#
# })


  # derive_orfgff(dna)
  # derive_orffaa(dna, orfgff)
  #
  # derive_trans(dna, gff)
  # derive_transorfgff(trans)
  # derive_transorffaa(trans)
  #
  # summarize_faa(transorffaa),
  # summarize_fna(trans),
  # summarize_gff(transorfgff),
  # summarize_nstring(nstring)
  #
  # load_synmap(filename)
  #
  # summarize_syn(synfile)
