test_that("Correct ratio for d13C standard returned", {
  expect_equal(calc_d13c(1.111618), -10.77, tolerance = 0.001)
})

# read reduced test data
us0625 <- readRDS(here::here("tests/us0625.rds"))

# produce new reduced data from raw file
us0625_test <- reduce_hgis(here::here("tests/USAMS062521R.txt"))

test_that("reduced data returned without error", {
  expect_type(us0625_test, "list")
  expect_equal(us0625_test, us0625)
})
