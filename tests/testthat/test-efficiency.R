test_that("Efficiency calc works", {
  expect_equal(hgis_eff(1,213), 1.000372, tolerance = 0.0001)
})
