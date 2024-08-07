% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeModel.R
\name{unapply.bgc.sel}
\alias{unapply.bgc.sel}
\title{Unapply bgc+selection parameters from a matrix}
\usage{
unapply.bgc.sel(m, bgc = 0, sel = 0, alphabet = "ACGT")
}
\arguments{
\item{m}{A transition matrix}

\item{bgc}{The bgc parameter which was used to calculate m}

\item{sel}{The selection parameter which was used to calculate m}

\item{alphabet}{The alphabet used for nucleotide states}
}
\value{
A matrix reflecting m before bgc and sel were applied.
}
\description{
Unapply bgc+selection parameters from a matrix
}
\author{
Melissa J. Hubisz and Adam Siepel
}
