
# Using values calculated by https://www.dolomite-microfluidics.com/support/microfluidic-calculator/ for comparison.
test_that("Pouiselle produces right numbers", {
  expect_equal(flowcalc(100, 0.000025, 0.001, 0.1), 9.2, tolerance = 0.004)
  expect_equal(flowcalc(100, 0.000025, 14.7E-6, 7), 8.95, tolerance = 0.004)
})
