\name{ppplot}
\alias{ppplot}

\title{PP plot for Poisson lognormal distribution}
\description{
    Generate a plot showing the quality of the fit of the Poisson lognormal to the data.
}

\usage{
    ppplot(n, n.points=10)
}

\arguments{
  \item{n}{vector of observed counts}
  \item{n.points}{number of points on the graph to highlight}
}

\details{
The function fits the Poisson lognormal to the raw read counts \code{n} and uses that
to generate theoretical percents for a PP (``percent-percent'' or ``probability-probability'')
plot. A perfect fit falls on the diagonal (marked with a dotted line).

The optional \code{n.points} plots some extra points on the graph. The first point (the lower-left-most)
represents the fraction of all OTUs that have 1 count in the empirical (x-axis) and theoretical
(y-axis) distributions. The second point represents OTUs with 1 or 2 counts. The third point
represents OTUs with up to 3 counts, and so on.
}

\value{
    a \code{ggplot} object whose default \code{data} object has columns \code{empirical} and \code{theoretical}
}

\seealso{
    \code{\link{texmex.fit}}
}

\examples{
# make up some data
n <- rpoilog(1000, 1.0, 1.0)

# plot it
p <- ppplot(n)
p

# compare to the lognormal's fit
# first, make the empirical cumulative distribution function from the data
x <- tabulate(n + 1)
empirical <- cumsum(x / sum(x))

# then, get the theoretical percents
theoretical <- plnorm(0:max(n), meanlog=mean(log(n)), sdlog=sd(log(n)))
lognormal.fit <- data.frame(empirical=empirical, theoretical=theoretical)

# add that data in a new layer
p + geom_line(data=lognormal.fit, color='red')
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{hplot}