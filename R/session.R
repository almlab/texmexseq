library(dplyr)
library(ggplot2)

ReadOtuTable <- function(fn) read.table(fn, header=T, row.names=1, check.names=F)

z.transform.counts <- function(n) {
    fit <- texmex.fit(n)
    z <- (log(n) - fit$par['mu']) / fit$par['sig']
    return(z)
}

z.transform.table <- function(otus) {
    z.otus <- as.data.frame(apply(otus, 2, z.transform.counts))
    return(z.otus)
}

f.transform.counts <- function(n) {
    fit <- texmex.fit(n)

    # get the pdf of the poilog (parameterized by the fit we just did) for 1, 2, ..., up
    # to the maximum N value we need. that's a table of theoretical pdf values. then take
    # the cumulative sum to get a table of F values.
    pdf.table <- dpoilog(0:max(n), fit$par['mu'], fit$par['sig'], trunc=TRUE)
    f.table <- cumsum(pdf.table)

    # now look up the F value for each of our OTU counts. (the +1 is there because the
    # f.table[1] is the pdf for 0 counts, so if N=1 we want f.table[1+1], which is the
    # pdf for 1 count.
    f <- f.table[n + 1]
    return(f)
}

f.transform.table <- function(otus) {
    f.otus <- as.data.frame(apply(otus, 2, f.transform.counts))
    return(f.otus)
}

ppplot <- function(n, n.points=10) {
    # convenience function for plotting
    fit <- texmex.fit(n)

    tpdf.table <- dpoilog(0:max(n), fit$par['mu'], fit$par['sig'], trunc=TRUE)
    f.table <- cumsum(tpdf.table)

    # make the empirical pdf values
    epdf.table <- tabulate(n[n > 0] + 1) # just get a raw "table"
    epdf.table <- epdf.table / sum(epdf.table) # divide by the sum, so all values sum to 1.0

    # turn the empirical pdf into a Cdf
    ecdf.table <- cumsum(epdf.table)

    pp.data <- data.frame(empirical=ecdf.table, theoretical=f.table)

    perfect.fit <- data.frame(empirical=c(0, 1), theoretical=c(0, 1))
    p <- ggplot(pp.data, aes(x=empirical, y=theoretical)) +
      geom_point(data=pp.data[1:n.points, ]) +
      geom_line() + 
      geom_line(data=perfect.fit) +
      coord_fixed() +
      xlim(0, 1) + ylim(0, 1) +
      theme_bw()
    return(p)
}

plot.pair.counts <- function(counts) {
    if (!all(sort(names(counts)) == c("before", "after"))) {
        stop("input a data frame of counts with two columns named 'before' and 'after'")
    }

    p <- ggplot(counts, aes(x=before, y=after)) + geom_point()
    return(p)
}

plot.quad <- function(f.table) {
    if (!all(c('control.before', 'control.after', 'treatment.before', 'treatment.after') %in% names(f.table))) {
        stop("input a data frame of F values with four columns named 'control.before', 'control.after', 'treatment.before', 'treatment.after'")
    }

    diff.f.table <- mutate(f.table, control=control.after-control.before, treatment=treatment.after-treatment.before)
    p <- ggplot(diff.f.table, aes(x=control, y=treatment)) +
      geom_point() + 
      coord_fixed() +
      xlim(-1, 1) + ylim(-1, 1) +
      theme_bw()
    return(p)
}

quad.table <- function(otu, control.before, control.after, treatment.before, treatment.after) {
  if (!all(c(control.before, control.after, treatment.before, treatment.after) %in% names(otu))) {
    stop("one of the specified names is not a column name in the input table")
  }
  data.frame(control.before=otu[[control.before]],
             control.after=otu[[control.after]],
             treatment.before=otu[[treatment.before]],
             treatment.after=otu[[treatment.after]])
}
