\name{poilogMLE}
\alias{poilogMLE}

\title{maximum likelihood estimation for the Poisson lognormal}
\description{
    Fit a set of OTU counts to the Poisson lognormal distribution.
}

\usage{
  poilogMLE(n, nboot=0, trunc=TRUE, method='L-BFGS-B', start.mu=-1.0, start.sig=1.0,
  control=list(maxit=1000, fnscale=length(n), parscale=c(1.0, 1.0)))
}

\arguments{
  \item{n}{vector of observed counts}
  \item{nboot}{number of bootstrapping iterations}
  \item{trunc}{remove weight from zero counts?}
  \item{method}{optimization method for \code{optim}}
  \item{start.mu}{initial estimate for the lognormal mean}
  \item{start.sig}{initial estimate for the lognormal standard deviation}
  \item{control}{list of parameters for \code{optim}}
}

\details{
The function estimates parameters mean mu and standard deviation sig.

The parameters must be given starting values for the optimization procedure.
The default values here worked well when fitting OTUs in the referenced paper.

The function uses the optimization procedures in \code{optim} to make the
maximum likelihood estimate. The \code{method} and \code{control} arguments are
passed to \code{optim}.

A zero-truncated distribution (see \code{dpoilog}) is assumed by default.
Truncation should only be turned off if all the known OTUs of the sequenced
community are known. In most cases, this is not applicable, since it is usually
not possible to know if an OTU had zero counts because it is not present in
the environment or if it is present but in low abundance.

Parametric bootstrapping simulate new sets of observations using the estimated
parameters (see \code{rbipoilog}).
}

\value{
    \item{par}{Maximum likelihood estimates of the parameters}
    \item{p}{Approximate fraction of OTUs revealed by the sample}
    \item{logLval}{Log likelihood of the data given the estimated parameters}
    \item{gof}{Godness of fit measure obtained by checking the rank of
        \code{logLval} against \code{logLval}'s obtained during the bootstrap
        procedure, (gof<0.05) or (gof>0.95) indicates lack of fit}
    \item{boot}{A data frame containing the bootstrap replicates of parameters and \code{logLval}}
}

\seealso{
    \code{\link{dpoilog}}
    \code{\link{optim}}
}

\examples{
# create some random data
x <- rpoilog(S=1000, mu=-2.0, sig=2.0, keep0=FALSE)

# fit that data
res <- poilogMLE(x)

# the results should be fairly robust to the starting parameters
res2 <- poilogMLE(x, start.mu=2.0, start.sig=0.5)
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}