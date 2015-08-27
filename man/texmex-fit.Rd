\name{texmex.fit}
\alias{texmext.fit}

\title{wrapper for fitting the Poisson lognormal}
\description{
    A nice wrapper for the \code{poilogMLE} function.
}

\usage{
  texmex.fit(n, start.mus=c(-2.0, -1.0, 0.0, 1.0, 2.0), start.sigs=rep(1.0, times=5), ...)
}

\arguments{
  \item{n}{vector of observed counts}
  \item{start.mus}{vector of starting \code{mu} values}
  \item{start.sigs}{vector of starting \code{sig} values}
  \item{...}{further parameters to go to \code{poilogMLE}}
}

\details{
It can be a little annoying to direclty use \code{poilogMLE}: the optimization can
fail to converge if you choose a bad starting \code{mu} value. This function makes
multiple attempts at fitting the poilog distribution, returning the first one that
works.

For the first attempt, the first \code{mu} and \code{sig} are used as starting values.
For the second attempt, the second elements of those vectors, etc.
}

\value{
    \item{par}{Maximum likelihood estimates of the parameters}
    \item{p}{Approximate fraction of OTUs revealed by the sample}
    \item{logLval}{Log likelihood of the data given the estimated parameters}
}

\seealso{
    \code{\link{poilogMLE}}
}

\examples{
# create some random data
x <- rpoilog(S=1000, mu=-2.0, sig=2.0, keep0=FALSE)

# fit that data
res <- texmex.fit(x)
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}