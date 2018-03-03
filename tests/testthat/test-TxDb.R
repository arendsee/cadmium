context("TxDb")

gff <- system.file("yeast", "gff", "Saccharomyces_arboricola.gff", package='fagin')
dna <- system.file("yeast", "fna", "Saccharomyces_arboricola.fna", package='fagin')

genomeDB <- load_dna(dna)
seqinfo <- GenomicFeatures::seqinfo(genomeDB)

# hold gene models Rmonad object 
m_gffDB <- NULL
gffDB <- NULL 

test_that("load gene models", {
  expect_silent({
    m_gffDB <<- load_gene_models(gff, seqinfo)
    gffDB <<- rmonad::get_value(m_gffDB, m_gffDB@head)[[1]]
  })
  expect_true(all(rmonad::get_OK(m_gffDB)))
  expect_equal(rmonad::get_value(m_gffDB, m_gffDB@head)[[1]] %>% top_class, "TxDb")
})

test_that("summarization doesn't explode", {
  expect_true(all(m_gffDB %>>% summarize_gff %>% rmonad::get_OK())) 
})

gffsum <- m_gffDB %>>% summarize_gff %>% rmonad::get_value(., .@head)
gffsum <- gffsum[[1]]

test_that("summarization is good", {
  expect_equal(top_class(gffsum), "gff_summary")
})

test_that("m_get_proteins works", {
  expect_silent(
    m_get_proteins(gffDB=gffDB, genomeDB=genomeDB, species_name="foo") %>% rmonad::esc()
  )
})
