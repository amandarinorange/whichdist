context("check-whichdist")
library(testthat)        # load testthat package
library(whichdist)       # load our package

# Test whether the output is a list
test_that("whichdist() returns a list", {
  output <- whichdist(df = toydata, countvar = "nb_count")
  expect_is(output, "list")
})

