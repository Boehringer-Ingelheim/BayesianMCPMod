test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("dummy tests", {
  expect_error(calc_pow())
  expect_error(simulateData())
  expect_error(postShape())
  expect_error(BayesMCPMod())
  expect_error(estimateModel())
})
