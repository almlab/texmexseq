\name{Experiment}
\alias{Experiment}
\alias{ReadOtuTable}
\alias{ExtractQuadData}

\title{Create many samples at once}
\description{
    Constructs and organizes samples, pairs, and quads using an OTU table
    and some auxiliary tables.
}

\usage{
  Experiment(otu.table, pair.table, quad.table)
  ReadOtuTable(fn)
  ExtractQuadData(expt, quad.name)
}

\arguments{
  \item{otu.table}{a filename or a dataframe (as produced by \code{ReadOtuTable})}
  \item{pair.table,quad.table}{a filename, as described below}
  \item{fn}{the OTU table filename}
  \item{expt}{the experiment object, as returned by \code{Experiment}}
  \item{quad.name}{a column name in \code{expt}}
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

  \code{ExtractQuadData} creates a data frame that contains all the information
  from a single quead: the counts, relative abundances, and changes in z and F.
}

\value{
    \item{ReadOtuTable}{returns a dataframe}
    \item{Experiment}{returns a list, including the OTU table, all the pairs,
    and all the quads}
    \item{ExtractQuadData}{returns a data frame}
}

\seealso{
    \code{\link{SamplePair}}
    \code{\link{SampleQuad}}
}

\author{Scott Olesen \email{swo@mit.edu}}

\keyword{distribution}