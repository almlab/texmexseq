context("Poilog MLE fitting")

test_that("poilogMLE fits data correctly", {
  fit <- poilogMLE(c(rep_len(1, 100), rep_len(2, 50)), -1.0, 1.0)
  expect_that(fit$par['mu'], equals(c(mu=-0.5011), tolerance=1e-3))
  expect_that(fit$p, equals(1.0))
  expect_that(fit$logLval, equals(-107.4394, tolerance=1e-3))
})

test_that("texmex.fit fits data correctly", {
  fit <- texmex.fit(c(rep_len(1, 100), rep_len(2, 50)))
  expect_that(fit$par['mu'], equals(c(mu=-0.5011), tolerance=1e-3))
  expect_that(fit$p, equals(1.0))
  expect_that(fit$logLval, equals(-107.4394, tolerance=1e-3))
})
