test_that("Correct ratio for d13C standard returned", {
  expect_equal(calc_d13c(1.111618), -10.77, tolerance = 0.001)
})

test_that("reduced data returned without error", {
  expect_type(reduce_hgis(here::here("data/USAMS062521R.txt")), "list")
})
