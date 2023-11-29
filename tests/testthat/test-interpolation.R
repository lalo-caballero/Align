test_that("Interpolation works", {
  expect_equal(interpolation(c(1:10), c(1:10)), list(f = c(2:9), s = c(2:9), g = rep(1,8)))
})
