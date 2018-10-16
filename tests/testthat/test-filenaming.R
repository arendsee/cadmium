context("filenaming")

tiny <- system.file("tiny", package='fagin')

test_that("get_readable_filename", {
  expect_equal(
    get_readable_filename("unicorn", dir=tiny, ext=c("fa", "fna", "fasta")),
    file.path(tiny, "unicorn.fna")
  )
  # die if no file found
  expect_error(
    get_readable_filename("dragon", dir=tiny, ext=c("fa", "fna", "fasta"))
  )
  # die on ambiguity
  expect_error(
    get_readable_filename("dragon", dir=tiny, ext=c("fna", "gff"))
  )
  expect_equal(
    get_synmap_filename(
      focal_name  = "unicorn",
      target_name = "centaur",
      dir         = tiny
    ),
    file.path(tiny, "unicorn.vs.centaur.syn")
  )
})
