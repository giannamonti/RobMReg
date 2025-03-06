context("dummy-fun")

test_that("is_posdef works for a basic posdef matrix", {
  m <- cov(airquality[,3:5])
  result <- is_posdef(m)
  expect_that(result, equals(TRUE))
})

