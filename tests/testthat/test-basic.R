test_that("Package loads correctly", {
  expect_true(require(gbrsR))
})

test_that("Default constants are available", {
  expect_true(exists("DEFAULT_CONFIG"))
  expect_true(exists("DEFAULT_FOUNDERS"))
  expect_true(exists("DEFAULT_FOUNDER_COLORS"))
  expect_true(exists("DEFAULT_CHROM_LENGTHS"))
})

test_that("Default founders has correct structure", {
  expect_equal(length(DEFAULT_FOUNDERS), 8)
  expect_true(all(DEFAULT_FOUNDERS %in% c("A", "B", "C", "D", "E", "F", "G", "H")))
})

test_that("Default config has required elements", {
  required_elements <- c("bar_height", "founder_gap", "chrom_spacing", "font_size_title")
  expect_true(all(required_elements %in% names(DEFAULT_CONFIG)))
})

test_that("generate_diplotypes works", {
  diplotypes <- generate_diplotypes(c("A", "B"))
  expect_equal(length(diplotypes), 3)  # AA, AB, BB
  expect_true(all(diplotypes %in% c("AA", "AB", "BB")))
})

test_that("configure_plot works", {
  config <- configure_plot(bar_height = 0.5)
  expect_equal(config$bar_height, 0.5)
  expect_equal(config$founder_gap, DEFAULT_CONFIG$founder_gap)  # unchanged
})

test_that("create_compact_config works", {
  config <- create_compact_config()
  expect_true(config$compact_mode)
  expect_equal(config$bar_height, 0.15)
  expect_equal(config$founder_gap, 0.05)
})
