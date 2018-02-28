context("protein extraction from genome")

fna <- load_dna(file.path("tiny", "unicorn.fna"))
gff <- derive_orfgff(convert_FaFile_to_XStringSet(fna))
cds <- NULL

test_that("load_dna builds DB and returns", {
  expect_silent(cds <<- extractWithComplements(fna, gff)) 
  expect_equal(top_class(cds), "DNAStringSet")
  expect_equal(translate(cds) %>% top_class, "AAStringSet")
})
