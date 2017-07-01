context("filenaming.R")


# NOTE: this test will correctly fail on Windows.
test_that("correctly name genome files", {
  expect_equal(
    # note, spaces in the species name should be replaced with underscores, but
    # this should not be done for directories, where they are legal (though
    # frowned upon by elitists like myself).
    get_genome_filename("Pauxi unicornus", dir="Zanzibar Wisconsin"),
    "Zanzibar Wisconsin/Pauxi_unicornus.fna"
  )
})

test_that("correctly name gff files", {
  expect_equal(
    get_gff_filename("Cheesemania stellata", dir="Narnia"),
    "Narnia/Cheesemania_stellata.gff3"
  )
})

test_that("correctly name syn files", {
  expect_equal(
    get_synmap_filename("Cheesemania stellata", "Pauxi unicornus", dir="Mozambique"),
    "Mozambique/Cheesemania_stellata.vs.Pauxi_unicornus.syn"
  )
})
