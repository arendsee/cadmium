context("reading")

test_that("load tree works", {
  expect_equal({ t <- load_tree("tiny/tiny.tree"); tail(t$tip.label, 1) }, "human")
})

test_that("load_gene_list", {
  expect_equal(
    load_gene_list(file.path('tiny', 'unicorn-genelist.txt')) %>% class,
    "character"
  )
})
