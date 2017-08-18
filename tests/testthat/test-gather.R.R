context("gather.R")

test_that("load tree works", {
  expect_equal({ t <- load_tree("tiny/tiny.tree"); tail(t$tip.label, 1) }, "human")
})
