context("Hagen-Poiseulle equation")

# Using values calculated by https://www.dolomite-microfluidics.com/support/microfluidic-calculator/ for comparison.
test_that("Pouiselle produces right numbers", {
  expect_equal(flowcalc(100, 0.000025, 0.001, 0.1), 9.2, tolerance = 0.004)
  expect_equal(flowcalc(100, 0.000025, 14.7E-6, 7), 8.95, tolerance = 0.004)
})


context("Vial flows")

test_that("Catches nonsense input", {
  expect_error(concCO2(0, 1, 1, 2))
})

test_that("Time 0 makes sense", {
  expect_equal(concCO2(0, 1, 1, 1), 1)
  expect_equal(concCO2(0), 1)
  expect_equal(concCO2(0, flow = 100), 100)
})

test_that("Step 1 makes sense", {
  expect_equal(concCO2(1, 1, 1), .368, tolerance = 0.001)
})

test_that("Doubling flow in capillary doubles CO2", {
  expect_equal(2 * concCO2(1, flow = 1), concCO2(1, flow = 2), tolerance = 0.001)
})

test_that("Starting at different fraction CO2 works", {
  expect_equal(concCO2(0, 1, 1, initco2 = 0.5), 0.5, tolerance = 0.001)
  expect_equal(concCO2(10, 1, 1, initco2 = 0.5), 0.44, tolerance = 0.001)
})