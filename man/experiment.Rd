\name{Experiment}
\alias{Experiment}
\alias{ReadOtuTable}

\title{Create many samples at once}
\description{
    Constructs and organizes samples, pairs, and quads using an OTU table
    and some auxiliary tables.
}

\usage{
  Experiment(otu.table, pair.table, quad.table)
  ReadOtuTable(fn)
}

\arguments{
  \item{otu.table}{a filename or a dataframe (as produced by \code{ReadOtuTable})}
  \item{pair.table,quad.table}{a filename, as described below}
  \item{fn}{the OTU table filename}
}

\details{
  \code{ReadOtuTable} expect a table of integer counts, with a first columnas OTU
  IDs and the top row as sample names. Do not include any lineage information.

  \code{Experiment} reads an OTU table, sample table, and quad table to produce
  many samples, pairs, and quads at once in a single object. 

  The pair table should be tab-delimited and have three columns with headers
  ``name'', ``before'', and ``after''. The names can be anything, but the before
  and after values must equal sample names in the OTU table.

  The quad table should be tab-delimited and have four columns with headers
  ``name'', ``control'', and ``treatment''. The names can be anything, but the
  control and treatment values must equal names from the pair table.
}

\value{
    \item{ReadOtuTable}{returns a dataframe}
    \item{Experiment}{returns a list, including the OTU table, all the pairs,
    and all the quads}
}

\seealso{
    \code{\link{SamplePair}}
    \code{\link{SampleQuad}}
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}