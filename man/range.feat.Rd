% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feat.R
\name{range.feat}
\alias{range.feat}
\title{Features range}
\usage{
\method{range}{feat}(..., na.rm = FALSE)
}
\arguments{
\item{...}{Objects of type \code{feat}}

\item{na.rm}{Whether to remove values of NA before calculating range.}
}
\value{
A vector of size 2 indicating minimum and maximum coord in
the features object
}
\description{
Get the range of a features object
}
\author{
Melissa J. Hubisz
}
\keyword{features}
