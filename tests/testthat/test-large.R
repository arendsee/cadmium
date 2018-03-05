context("whole pipeline")

# This tests entire pipelines. You may want to turn it off for most tests.

a1 <- NULL
con <- config()

test_that("analysis_1 works for yeast", {
  expect_true(
    a1 <<- load_species(species_name = con@input@focal_species, con=con)
    all(get_OK(a1))
  )
})

# TODO: cleanup
