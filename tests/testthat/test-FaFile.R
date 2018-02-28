context("FaFile types")

fna <- system.file("yeast", "fna", "Saccharomyces_arboricola.fna", package='fagin')
fai <- paste0(fna, ".fai")

unlink(fai)
test_that("load_dna builds DB and returns", {
  expect_silent(load_dna(fna))
  expect_true(file.exists(fai))
  expect_true(class(load_dna(fna))[1] == "FaFile")
})

test_that("summary of genomes is what it is", {
  expect_equal(
    load_dna(fna) %>% summarize_dna %>% top_class,
    "dna_summary"
  )
})

test_that("can convert FaFile to DNAStringSet", {
  expect_equal(
    load_dna(fna) %>% convert_FaFile_to_XStringSet %>% top_class,
    "DNAStringSet"
  )
})
