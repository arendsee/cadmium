context("protein extraction from genome")

fna_file <- system.file("yeast", "fna", "Saccharomyces_arboricola.fna", package='fagin')
gff_file <- system.file("yeast", "gff", "Saccharomyces_arboricola.gff", package='fagin')

fna <- load_dna(fna_file)
gff <- derive_orfgff(convert_FaFile_to_XStringSet(fna))
cds <- NULL

test_that("load_dna builds DB and returns", {
  expect_silent(cds <<- extractWithComplements(fna, gff)) 
  expect_equal(top_class(cds), "DNAStringSet")
  expect_equal(translate(cds) %>% top_class, "AAStringSet")
})
