context("protein extraction from genome")

fna_file <- system.file("yeast", "fna", "Saccharomyces_arboricola.fna", package='fagin')
gff_file <- system.file("yeast", "gff", "Saccharomyces_arboricola.gff", package='fagin')

tiny <- system.file("tiny", package='fagin')

fna <- load_dna(fna_file)
con <- config()
con@orf@start <- "ATG"
orf <- derive_genomic_ORFs(convert_FaFile_to_XStringSet(fna), con)
cds <- NULL

test_that("extract and translate work", {
  expect_silent(cds <<- extractWithComplements(fna, orf))
  expect_equal(top_class(cds), "DNAStringSet")
  expect_equal(fuzzy_translate(cds) %>% top_class, "AAStringSet")
})

fna_file <- file.path(tiny, "test.fna")
gff_file <- file.path(tiny, "test.gff")
fna <- load_dna(fna_file)

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


con@orf@minlen <- 3L
dna <- Biostrings::DNAStringSet(
  c("a"="ATGTTTAAATAA",    # too short
    "b"="ATGTTTAAAGGGTAA", # just right
    "c"="ATGTTATGTTTTTAAAATAA", # overlapping
    "d"="GGGTTAGGGATGAAACATAAATAA" # both strands 
  ))
test_that("ORF finding works as expected", {
  genorfs <- derive_genomic_ORFs(dna, con)
  expect_equal(as.character(GenomeInfoDb::seqnames(genorfs)), c("b", "c", "c", "d", "d"))
  expect_equal(GenomicRanges::start(genorfs), c(1,1,6,10,4))
  expect_equal(as.character(GenomicRanges::strand(genorfs)), c("+", "+", "+", "+", "-"))
  expect_equal(unname(as.character(Biostrings::translate(extract_range(dna, genorfs)))),
               c("MFKG*", "MLCF*", "MFLK*", "MKHK*", "MFHP*"))

  transorfs <- derive_transcript_ORFs(dna, con)
  expect_equal(as.character(GenomeInfoDb::seqnames(transorfs)), c("b", "c", "c", "d"))
  expect_equal(GenomicRanges::start(transorfs), c(1,1,6,10))
  expect_equal(as.character(GenomicRanges::strand(transorfs)), c("+", "+", "+", "+"))
})
