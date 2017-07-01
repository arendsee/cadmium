context("io.R")


attr_simple <- c(
  "ID=123",
  "ID=234;Parent=123"
)

d_simple <- data.frame(
  ID      = c("123", "234"),
  Parent  = c(NA_character_, "123"),
  .n_tags = c(1,2),
  stringsAsFactors=FALSE
)

d_simple_u <- data.frame(
  ID       = c("123", "234"),
  Parent   = c(NA_character_, "123"),
  Untagged = c(NA_character_, NA_character_),
  .n_tags  = c(1,2),
  stringsAsFactors=FALSE
)

d_simple_waldo <- data.frame(
  waldo    = c(NA_character_, NA_character_),
  .n_tags  = c(1,2),
  stringsAsFactors=FALSE
)

test_that("tag extraction works in the normal case", {
  expect_equal(extract_tags(attr_simple, c("ID", "Parent")), d_simple)
  expect_equal(extract_tags(attr_simple, c("ID", "Parent"), get_naked=TRUE), d_simple_u)
  expect_equal(extract_tags(attr_simple, c("waldo")), d_simple_waldo)
  expect_error(extract_tags(attr_simple, c()))
})
