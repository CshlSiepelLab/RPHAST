% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgc.R
\name{bgc.informative}
\alias{bgc.informative}
\title{Return features indicating regions informative for bgc}
\usage{
bgc.informative(align, foreground, tree, not.informative = FALSE)
}
\arguments{
\item{align}{An MSA object representing a multiple alignment}

\item{foreground}{A character string giving the name of a branch (or
a label given to several branches) indicating which branch should be
in the foreground.  The foreground branch is where GC-biased gene
conversion is applied, and, if using a coding model, is where a test
of positive selection can be performed.}

\item{tree}{The phylogenetic tree to be used.  Can be a newick
string describing a tree, or an object of type \code{tm}.}

\item{not.informative}{If TRUE, return the regions that are not
informative for bgc.}
}
\value{
An object of type \code{feat} indicating which regions
are informative for bgc on the named foreground branch.  If
\code{not.informative==TRUE}, it will instead return the inverse of
this, indicating which regions are not informative for bgc.  The
coordinates of the features object are in the frame of the reference
species of the alignment.
}
\description{
Return features indicating regions informative for bgc
}
\author{
Melissa J. Hubisz
}
