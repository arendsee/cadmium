context("preprocess.R")

dna <- load_dna(file.path("asdf", "asdf.fna"))
gff <- load_gff(file.path("asdf", "asdf.gff"))

test_that("load_dna works", {
  expect_true(all(Biostrings::width(dna) == 3000))
  # assert that description was removed from the header
  expect_true( grepl("NC_\\d+\\.\\d+", names(dna), perl=TRUE) %>% all )
})

test_that("nstrings gets the right intervals", {
  expect_equal(derive_nstring(dna) %>% Biostrings::width(), c(3,80,5,68,4,480,1,68))
  expect_equal(derive_nstring(dna) %>% Biostrings::start(), c(1,241,404,649,889,81,1208,2328))
})

# test_that("derive_orfgff creates valid ORFs", {
#   # Here I test only that the ouput are ORFs
#   # Currently the way I find ORFs is totally wrong
#
# })


  # derive_orfgff(dna)
  # derive_orffaa(dna, orfgff)
  #
  # derive_aa(dna, gff)
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
