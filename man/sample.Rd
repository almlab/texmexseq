\name{Sample}
\alias{SampleFit}

\title{make and fit samples of OTUs}
\description{
    A construct that holds OTU counts for a sample. Provides nicer interfaces
    for plotting and manipulating data.
}

\usage{
    Sample(n, name=NULL)
    SampleFit(sample, trunc=TRUE, verbose=TRUE)
}

\arguments{
  \item{n}{vector of OTU counts}
  \item{name}{a name for the sample to be used in plots}
  \item{sample}{a sample as returned by \code{Sample}}
  \item{trunc}{ignore and reweight zero counts?}
  \item{verbose}{print messages about fitting?}
}

\details{
You can construct samples of this type bottom-up from OTU count data. They are also
produced top-down when using the higher-level \code{Experiment} command. Samples
can be joined into pairs.
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

# visualize the distribution of counts in a sample
plot(s)

# make a fit object using that sample
f <- SampleFit(s)

# visualize the quality of the fit
plot(f)
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}