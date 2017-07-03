context("io.R")



# ===================================
#    Checks for GFF error catching   
# ===================================

test_that("Repeated tags trigger error", {
  expect_silent ( extract_tags("A=1;B=2", "B") )
  expect_error  ( extract_tags("B=1;B=2", c()) )
})

test_that("Repeated tags trigger error", {
  expect_silent ( extract_tags("A=1;B=2", "B") )
  expect_error  ( extract_tags("B=1;B=2", "B") )
})

test_that("Multiple equals or OK", {
  expect_silent ( extract_tags("B=2",   "B") )
  expect_silent ( extract_tags("B=2=1", "B") )
})

test_that("Commas trigger errors", {
  expect_silent ( extract_tags("B=2", "B") )
  expect_error  ( extract_tags("B=1,2", "B") )
})



# =======================
#     Simple GFF cases   
# =======================

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
  ID      = c("123", "234"),
  Parent  = c(NA_character_, "123"),
  .U      = c(NA_character_, NA_character_),
  .n_tags = c(1,2),
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
})



# ==================================
#    Handling for AUGUSTUS output   
# ==================================

attr_augustus <- c(
  "cryptic_id",
  "ID=bill;Parent=cryptic_id"
)

d_augustus_i <- data.frame(
  ID      = c("cryptic_id", "bill"),
  Parent  = c(NA_character_, "cryptic_id"),
  .n_tags = c(1,2),
  stringsAsFactors=FALSE
)

d_augustus <- data.frame(
  ID      = c(NA_character_, "bill"),
  Parent  = c(NA_character_, "cryptic_id"),
  .n_tags = c(1,2),
  stringsAsFactors=FALSE
)

test_that("An untagged value is cast as ID if no other attributes are present", {
  expect_equal(extract_tags(attr_augustus, c("ID", "Parent"), infer_id=TRUE), d_augustus_i)
  expect_equal(extract_tags(attr_augustus, c("ID", "Parent"), infer_id=FALSE), d_augustus)
})
