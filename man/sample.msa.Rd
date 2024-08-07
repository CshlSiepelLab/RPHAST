% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msa.R
\name{sample.msa}
\alias{sample.msa}
\title{Sample columns from an MSA}
\usage{
\method{sample}{msa}(x, size, replace = FALSE, prob = NULL, pointer.only = FALSE)
}
\arguments{
\item{x}{An object of type \code{msa}}

\item{size}{The number of columns to sample}

\item{replace}{Whether to sample with replacement}

\item{prob}{A vector of probability weights for sampling each column;
\code{prob=NULL} implies equal probability for all columns.  Probabilities
need not sum to one but should be non-negative and can not all be zero.}

\item{pointer.only}{If \code{TRUE}, return only a pointer to an alignment
object stored in C (useful for large objects; advanced use only).}
}
\value{
An object of type \code{msa} with columns randomly
re-sampled from the original
}
\description{
Sample columns from an MSA
}
\note{
This function is implemented using R's sample function in
conjunction with "[.msa".  It will not alter the value of x even if it
is stored as a pointer.
}
\examples{
m <- msa(seqs=c("AAAAAAAAAACCCCCGGT", "GGGGGGGGGGTTTTTCCA", "CCCCCCCCCCAAAAAGGA"),
         names=c("human", "mouse", "rat"))
sample.msa(m, 10, replace=TRUE)
sample.msa(m, 10, replace=TRUE, prob=c(rep(1, 10), rep(2, 5), rep(5, 2), 10))
}
\author{
Melissa J. Hubisz and Adam Siepel
}
\keyword{msa}
