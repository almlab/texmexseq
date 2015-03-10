\name{SamplePair}
\alias{SampleQuad}

\title{make and fit pairs or quads of samples}
\description{
    Constructs for organizing samples that are control/experiment pairs, or
    for organzining before/after pairs of pairs (called quads).
}

\usage{
    SamplePair(sample0, sample1, name=NULL, trunc0=TRUE, trunc1=TRUE)
    SampleQuad(control, treatment, name=NULL)
}

\arguments{
  \item{sample0,sample1}{samples as produced by \code{Sample}}
  \item{name}{a name for the pair or quad to be used in plots}
  \item{control,treatment}{pairs as produced by \code{SamplePair}}
  \item{sample}{a sample as returned by \code{Sample}}
  \item{trunc0,trunc1}{ignore and reweight zero counts?}
}

\details{
You can construct samples of this type bottom-up from OTU count data. From
four samples (control before, control after, treatment before, treatment after),
you can make two pairs (control, treatment) and thus one quad. These constructs
are also mode top-down when using the higher-level \code{Experiment} command.
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
random.sample <- function() Sample(rpoilog(S=1000, mu=-2.0, sig=2.0))

ct.before <- random.sample()
ct.after <- random.sample()
tx.before <- random.sample()
tx.after <- random.sample()

# join these up into pairs
ct <- SamplePair(ct.before, ct.after)
tx <- SamplePair(tx.before, tx.after)

# visualize the change in counts in a pair
plot(ct)

# create a quad
q <- SampleQuad(ct, tx)

# visualize the changes in rescaled read counts
plot(q)

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}