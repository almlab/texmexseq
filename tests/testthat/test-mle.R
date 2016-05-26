context("Poilog MLE fitting")

test_that("poilogMLE fits data correctly", {
  fit <- poilogMLE(c(rep_len(1, 100), rep_len(2, 50)), -1.0, 1.0)
  expect_that(fit$par['mu'], equals(c(mu=-0.501100851401378)))
  expect_that(fit$par['sig'], equals(c(sig=0.002507109553842)))
  expect_that(fit$p, equals(1.0))
  expect_that(fit$logLval, equals(-107.439478784575))
})

test_that("texmex.fit fits data correctly", {
  fit <- texmex.fit(c(rep_len(1, 100), rep_len(2, 50)))
  expect_that(fit$par['mu'], equals(c(mu=-0.50112872990046842)))
  expect_that(fit$par['sig'], equals(c(sig=0.00132585466481472)))
  expect_that(fit$p, equals(1.0))
  expect_that(fit$logLval, equals(-107.439430897867))
})
