% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{numnodes.tree}
\alias{numnodes.tree}
\title{Number of Nodes in a Tree}
\usage{
numnodes.tree(tree)
}
\arguments{
\item{tree}{A vector of character strings, each containing a newick tree}
}
\value{
A numeric vector containing the number of nodes in each tree
}
\description{
Get the number of nodes in a tree
}
\examples{
numnodes.tree(c("((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);",
                "(human, (mouse, rat));"))
   
}
\author{
Melissa J. Hubisz and Adam Siepel
}
\keyword{trees}
