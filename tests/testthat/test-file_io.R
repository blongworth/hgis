
test_that("Reading results file works", {
  expect_type(read_results_file("USAMS062521R.txt"), "list")
})

test_that("Processing results file works", {
  expect_type(process_hgis_results("USAMS062521R.txt"), "list")
})
