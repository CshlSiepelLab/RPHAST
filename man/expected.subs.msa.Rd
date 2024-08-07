% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msa.R
\name{expected.subs.msa}
\alias{expected.subs.msa}
\title{Obtain expected number of substitutions on each branch and site}
\usage{
expected.subs.msa(x, tm)
}
\arguments{
\item{x}{An object of type \code{msa}}

\item{tm}{An object of type \code{tm}}
}
\value{
An array giving the expected number of substitutions on each
branch at each unique site pattern, summed across all types of
substitutions.
}
\description{
Obtain expected number of substitutions on each branch and site
}
\examples{
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, "ENr334-100k.maf")
m <- read.msa("ENr334-100k.maf")
mod <- phyloFit(m, tree="((hg18,(mm9,rn4)),canFam2)")
x <- expected.subs.msa(sub.msa(m, start.col=41447839, end.col=41448033, refseq="hg18"), mod)
dim(x)
dimnames(x)
x[,"CCCC"]
x["mm9-rn4",]
unlink("ENr334-100k.maf")
}
\author{
Melissa J. Hubisz and Adam Siepel
}
