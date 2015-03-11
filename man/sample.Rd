\name{Sample}
\alias{Sample}
\alias{SampleFit}
\alias{PlotFitPP}

\title{make and fit samples of OTUs}
\description{
    A construct that holds OTU counts for a sample. Provides nicer interfaces
    for plotting and manipulating data.
}

\usage{
    Sample(n, name=NULL)
    SampleFit(sample, trunc=TRUE, verbose=TRUE, ...)
    PlotFitPP(f, log=FALSE, npoints=10)
}

\arguments{
  \item{n}{vector of OTU counts}
  \item{name}{a name for the sample to be used in plots}
  \item{...}{extra arguments passed to \code{poilogMLE}}
  \item{sample}{a sample as returned by \code{Sample}}
  \item{trunc}{ignore and reweight zero counts?}
  \item{verbose}{print messages about fitting?}
  \item{f}{the fit object}
  \item{log}{use log scales?}
  \item{npoints}{number of OTUs to mark on the PP plot}
}

\details{
  You can construct samples of this type bottom-up from OTU count data. They are also
  produced top-down when using the higher-level \code{Experiment} command. Samples
  can be joined into pairs.

  \code{PlotFitPP} makes a PP plot that compares the empirical distribution
  function against the one from the Poisson lognormal fit.
}

\value{
    \item{Sample}{returns a sample object}
    \item{SampleFit}{returns a fit object, like the output of \code{poilogMLE}}
}

\seealso{
    \code{\link{SamplePair}}
    \code{\link{Experiment}}
    \code{\link{poilogMLE}}
}

\examples{
# create some data
counts <- rpoilog(S=1000, mu=-2.0, sig=2.0)

# create a sample out of that
s <- Sample(counts)

# make a fit object using that sample
f <- SampleFit(s)

# visualize the quality of the fit
PlotFitPP(f)
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}
