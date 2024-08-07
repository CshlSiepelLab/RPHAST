% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{read.newick.tree}
\alias{read.newick.tree}
\title{Read a Newick Tree from a File}
\usage{
read.newick.tree(filename)
}
\arguments{
\item{filename}{The file containing the tree.}
}
\value{
a character string representing the tree in newick format
}
\description{
Read a tree from a file
}
\details{
Reads a tree in newick format
}
\examples{
cat(c("((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):0.020666,canFam2:0.193569);",
      "(human, (mouse, rat));",
      sep="\n"), file="test.nh")
read.newick.tree("test.nh")
unlink("test.nh")
}
\author{
Melissa J. Hubisz and Adam Siepel
}
\keyword{newick}
\keyword{trees}
