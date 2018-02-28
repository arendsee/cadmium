context("TxDb")

gff <- system.file("yeast", "gff", "Saccharomyces_arboricola.gff", package='fagin')
dna <- system.file("yeast", "fna", "Saccharomyces_arboricola.fna", package='fagin')

seqinfo <- GenomicFeatures::seqinfo(load_dna(dna))

# hold gene models Rmonad object 
m <- NULL

test_that("load gene models", {
  expect_silent(m <<- load_gene_models(gff, seqinfo))
  expect_true(all(rmonad::get_OK(m)))
  expect_equal(rmonad::get_value(m, m@head)[[1]] %>% top_class, "TxDb")
})

test_that("summarization doesn't explode", {
  expect_true(all(m %>>% summarize_gff %>% rmonad::get_OK())) 
})

gffsum <- m %>>% summarize_gff %>% rmonad::get_value(., .@head)
gffsum <- gffsum[[1]]

test_that("summarization is good", {
  expect_equal(top_class(gffsum), "gff_summary")
})
