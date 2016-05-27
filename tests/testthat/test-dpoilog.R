context("Poilog density")

test_that("dpoilog gives correct single answers", {
  expect_that(dpoilog(1, 1.0, 1.0), equals(0.208462862162023))
  expect_that(dpoilog(1, 1.0, 1.0, trunc=FALSE), equals(0.175733342729315))
  expect_that(dpoilog(0, 1.0, 1.0), equals(0.0))
})
