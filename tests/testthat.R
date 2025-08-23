# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# To run tests, use:
#   testthat::test_local()
#   testthat::test_check("gbrsR")
#
# To run tests interactively, use:
#   testthat::test_file("tests/testthat/test-*.R")
#
# To run a single test, use:
#   testthat::test_file("tests/testthat/test-*.R", filter = "test_name")

library(testthat)
library(gbrsR)

test_check("gbrsR")
