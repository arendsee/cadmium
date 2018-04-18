context("whole pipeline")

# This tests entire pipelines. You may want to turn it off for most tests.

a1 <- NULL
con <- config()

test_that("analysis_1 works for yeast", {
  expect_true({
    a1 <<- load_species(species_name = con@input@focal_species, con=con)
    all(rmonad::get_OK(a1))
  })
})

# faafile <- system.file("yeast", "test-faa", "Saccharomyces_cerevisiae.gff", package='fagin')
# exp_faa <- Biostrings::readAAStringSet(faafile)
#
# test_that("the correct proteins were extracted", {
#   expect_true({
#       # obs_faa <- rmonad::get_value(a1, tag="faa")[[1]]
#       #
#       # obs_only <- Biostrings::setdiff(obs_faa, exp_faa)
#       #
#       # exp_only <- Biostrings::setdiff(exp_faa, obs_faa)
#
#     TRUE
#   })
# })

# TODO: cleanup
