context("TxDb")

gff <- file.path("tiny", "unicorn.gff")
seqinfo <- GenomicFeatures::seqinfo(load_dna(file.path("tiny", "unicorn.fna")))

# hold gene models Rmonad object 
m <- NULL

test_that("load gene models", {
  expect_silent(m <<- load_gene_models(gff, seqinfo))
  expect_true(all(rmonad::get_OK(m)))
  expect_equal(rmonad::get_value(m, m@head)[[1]] %>% top_class, "TxDb")
})

test_that("summarization doesn't explode", {
  expect_true(all(m %>>% summarize_gff %>% get_OK)) 
})

gffsum <- m %>>% summarize_gff %>% rmonad::get_value(., .@head)
gffsum <- gffsum[[1]]

test_that("summarization is good", {
  expect_equal(top_class(gffsum), "gff_summary")
})
