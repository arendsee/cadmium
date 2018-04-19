context("test `DNAStringSet -> a` functions")

fna_file <- system.file("yeast", "fna", "Saccharomyces_arboricola.fna", package='fagin')
gff_file <- system.file("yeast", "gff", "Saccharomyces_arboricola.gff", package='fagin')

cds_sample <- Biostrings::DNAStringSet(c(
  a='ATGTTTTAA',
  b='ATG',
  c='TAA',
  d='TAATTT',
  e='NNN'
))

cds_sample_ambiguous <- Biostrings::DNAStringSet('CCN')

test_that("translate works", {
  expect_equal(
    translate(cds_sample), Biostrings::AAStringSet(c(a="MF*", b="M", c="*", d="*F", e="X"))
  )
  expect_equal(
    translate(cds_sample_ambiguous), Biostrings::AAStringSet("P")
  )
  expect_warning(
    fagin:::translate(Biostrings::DNAStringSet(c(x="A")), label='unicorn')
  )
})

test_that("filter_with_warnings__zero_length_proteins", {
  expect_warning(
    Biostrings::AAStringSet(c("GANDALF", "")) %>%
      filter_with_warnings__zero_length_proteins
  )
  expect_equal(
    Biostrings::AAStringSet(c("GANDALF", "")) %>%
      filter_with_warnings__zero_length_proteins,
    Biostrings::AAStringSet(c("GANDALF"))
  )
})

test_that("Check nstrings", {
  expect_equal(
    cds_sample %>% derive_nstring %>% top_class,
    "GRanges"
  )
  expect_equal(cds_sample %>% derive_nstring %>% length, 1)
  expect_equal(
    cds_sample                %>%
    derive_nstring            %>%
    GenomicRanges::seqnames() %>%
    as.character, "e"
  )
  expect_equal(
    cds_sample %>% derive_nstring %>% summarize_nstring %>% unname,
    as.vector(summary(3))
  )
})

orfgff <- NULL

# TODO: make a simple test, ensure it works CORRECTLY
# Because the current implementation is dead wrong
test_that("ngORF prediction runs", {
  expect_equal(
    {
      orfgff <<- load_dna(fna_file) %>%
        convert_FaFile_to_XStringSet %>%
        derive_orfgff
      top_class(orfgff)
    },
    "GRanges"  
  ) 
})

test_that("convert_GRanges_to_SynderGFF", {
  expect_true(identical(
    GenomicRanges::ranges(convert_GRanges_to_SynderGFF(orfgff)),
    GenomicRanges::ranges(orfgff)
  ))
})

test_that("Can summarize GRanges", {
  expect_equal(
    summarize_granges(orfgff) %>% top_class,
    "granges_summary"
  )
})
