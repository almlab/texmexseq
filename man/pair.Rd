\name{SamplePair}
\alias{SamplePair}
\alias{PlotPair}
\alias{SampleQuad}
\alias{PlotQuad}

\title{make and fit pairs or quads of samples}
\description{
    Constructs for organizing samples that are control/experiment pairs, or
    for organzining before/after pairs of pairs (called quads).
}

\usage{
    SamplePair(sample0, sample1, name=NULL, trunc0=TRUE, trunc1=TRUE)
    PlotPair(pair, log='xy', fit=TRUE, highlight=NULL)
    SampleQuad(control, treatment, name=NULL)
    PlotQuad(quad, highlight=NULL, dF=FALSE, jitter.amount=0.0)
}

\arguments{
  \item{sample0,sample1}{samples as produced by \code{Sample}}
  \item{name}{a name for the pair or quad to be used in plots}
  \item{trunc0,trunc1}{ignore and reweight zero counts?}
  \item{pair}{pair object to be plotted}
  \item{log}{use logarithmic axes}
  \item{fit}{include a best-fit line?}
  \item{highlight}{a boolean vector of OTUs to highlight}
  \item{control,treatment}{pairs as produced by \code{SamplePair}}
  \item{quad}{quad object to be plotted}
  \item{dF}{show dF values? otherwise dz}
  \item{jitter.amount}{amount by which to scatter points in the plot}
}

\details{
  You can construct samples of this type bottom-up from OTU count data. From
  four samples (control before, control after, treatment before, treatment after),
  you can make two pairs (control, treatment) and thus one quad. These constructs
  are also mode top-down when using the higher-level \code{Experiment} command.

  \code{PlotPair} shows the counts for each OTU in the two samples in the pair.

  \code{PlotQuad} shows the dz values for each OTU in the two pairs in the quad.
}

\value{
    \item{SamplePair}{returns a pair object}
    \item{SampleQuad}{returns a quad object}
}

\seealso{
    \code{\link{Sample}}
    \code{\link{Experiment}}
}

\examples{
# create some data
random.sample <- function() Sample(rpoilog(S=1000, mu=-2.0, sig=2.0, condS=TRUE))

ct.before <- random.sample()
ct.after <- random.sample()
tx.before <- random.sample()
tx.after <- random.sample()

# join these up into pairs
ct <- SamplePair(ct.before, ct.after)
tx <- SamplePair(tx.before, tx.after)

# visualize the change in counts in a pair
PlotPair(ct)

# create a quad
q <- SampleQuad(ct, tx)

# visualize the changes in rescaled read counts
PlotQuad(q)
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}
