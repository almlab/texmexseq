\name{texmexseq-package}
\alias{texmexseq-package}
\alias{texmexseq}
\docType{package}

\title{Treatment Effect eXplorer for Microbial Ecology eXperiments}

\description{
Fit OTU count data to the Poisson lognormal distribution to create values
that can be compared across control and treatment units in before and
after timepoints. Plotting methods for these ``quads'' of samples.
}

\details{
\tabular{ll}{
Package: \tab texmexseq\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-03-10\cr
License: \tab GPL-3\cr
}

\code{texmexseq} is essentially a wrapper around the earlier R package
\code{poilog}. The implementations of the basic functions related to the Poisson
lognormal distribution in this package are nearly identical to the ones in
\code{poilog}. There are convenience functions, notably \code{Experiment}
deisgned to directly interact with OTU tables.
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{package}
