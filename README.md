## texmexseq
* Scott Olesen <swo@mit.edu>
* project page: https://almlab.mit.edu/texmex.html
* cran page: http://cran.r-project.org/web/packages/texmexseq/index.html
* github page: https://github.com/almlab/texmexseq

Treatment Effect eXplorer for Microbial Ecology eXperiments that use sequence counts (aka
`texmexseq`) is an R package designed to normalize OTU count data and correct for community
composition changes that are common to a control and experimental unit.

## How to use it
There is a demo of the code in the `demo` folder of the source distribution. The Example
sections of the help pages for the functions in the package also have working commands.

## Installation
The package in on CRAN, so it should be as easy as `install.packages('texmexseq')`.

## New in this version
Version 0.2 is a non-backward-compatible rework of version 0.1. The underlying
theory and computational methods are the same. Rather than creating Pair and Quad
functions, the transformations are made directly to the OTU tables. Plots are
made using `ggplot2`, allowing easier manipulation and decoration. OTU tables are
manipulated using `dplyr`, allowing easier filtering for interesting OTUs.

## Acknowledgement
This package is mostly a rewording of the `poilog` package, designed to make it easier for
microbial ecologists.

## License
GPL-3
