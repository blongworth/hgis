
test_that("Reading results file works", {
  expect_type(read_results_file(here::here("tests/USAMS062521R.txt")), "list")
})

test_that("Processing results file works", {
  expect_type(process_hgis_results(here::here("tests/USAMS062521R.txt")), "list")
})
