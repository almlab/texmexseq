\name{dpoilog}
\alias{dpoilog}
\alias{rpoilog}

\title{probability density function for the Poisson lognormal distribution}
\description{
    Probability density function for the Poisson lognormal distribution, along with
    methods for drawing random samples from that distribution and fitting experimental
    data to that distribution.
}

\usage{
    dpoilog(n, mu, sig, trunc=TRUE)
    rpoilog(S, mu, sig, condS=FALSE, keep0=FALSE)
}

\arguments{
  \item{n}{vector of observed counts}
  \item{S}{number of OTUs to sample}
  \item{mu}{mean of lognormal part of distribution}
  \item{sig}{standard deviation of lognormal part of distribution}
  \item{trunc}{remove weight from zero counts?}
  \item{condS}{are random samples conditional on \code{S}?}
  \item{keep0}{do not discard zero counts?}
}

\details{
A detailed description of the calculation of these values can be found in the
documentation for the \code{dpoilog} in the \code{poilog} package.
}

\value{
    \item{dpoilog}{returns that densities}
    \item{rpoilog}{returns the random counts}
}

\seealso{
    \code{\link{poilogMLE}}
}

\examples{
# visualize the density function
plot(dpoilog(0:100, mu=1.0, sig=1.0), type='o')

# visualize the empirical distribution of a random sample of 1000 OTUs
hist(rpoilog(1000, mu=1.0, sig=1.0, condS=TRUE))
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{univar}
\keyword{datagen}