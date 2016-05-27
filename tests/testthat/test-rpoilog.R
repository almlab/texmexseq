context("Poilog random variates")

test_that("rpoilog gives numbers", {
  expect_that(rpoilog(10, 1.0, 1.0) %>% is.numeric(), is_true())
})
