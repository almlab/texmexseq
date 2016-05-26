context("Table manipulation")

test_that("quad.table pulls out the right columns", {
  otu.table <- read.table('otus.dat', header=T, row.names=1)

  f.table <- f.transform.table(otu.table)
  f.quad <- quad.table(f.table, 'inoculum1.control.before', 'inoculum1.control.after', 'inoculum1.treatment.before', 'inoculum1.treatment.after')
  expect_that(colnames(f.quad), equals(c('control.before', 'control.after', 'treatment.before', 'treatment.after', 'd.control', 'd.treatment', 'otu.id')))

  z.table <- z.transform.table(otu.table)
  z.quad <- quad.table(z.table, 'inoculum1.control.before', 'inoculum1.control.after', 'inoculum1.treatment.before', 'inoculum1.treatment.after')
  expect_that(colnames(z.quad), equals(c('control.before', 'control.after', 'treatment.before', 'treatment.after', 'd.control', 'd.treatment', 'otu.id')))
})
