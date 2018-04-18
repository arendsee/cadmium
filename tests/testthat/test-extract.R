context("protein extraction from genome")

fna_file <- system.file("yeast", "fna", "Saccharomyces_arboricola.fna", package='fagin')
gff_file <- system.file("yeast", "gff", "Saccharomyces_arboricola.gff", package='fagin')

tiny <- system.file("tiny", package='fagin')

fna <- load_dna(fna_file)
orf <- derive_orfgff(convert_FaFile_to_XStringSet(fna))
cds <- NULL

test_that("extract and translate work", {
  expect_silent(cds <<- extractWithComplements(fna, orf)) 
  expect_equal(top_class(cds), "DNAStringSet")
  expect_equal(translate(cds) %>% top_class, "AAStringSet")
})

fna_file <- file.path(tiny, "test.fna")
gff_file <- file.path(tiny, "test.gff")
fna <- load_dna(fna_file)

# TODO: MAKE IT WORK
test_that("can extract multi", {
  expect_equal({
      m_gffDB <- load_gene_models(gff_file, seqinfo_=GenomeInfoDb::seqinfo(fna))
      gffDB <- rmonad::get_value(m_gffDB, m_gffDB@head)[[1]]
      m <- m_get_proteins(gffDB=gffDB, genomeDB=fna, species_name="foo")
      seq <- rmonad::get_value(m, tag='faa')[[1]]
      names(seq) <- NULL
      seq
    },
    Biostrings::AAStringSet("MMM*")
  )
})
