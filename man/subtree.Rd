% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{subtree}
\alias{subtree}
\title{Subtree}
\usage{
subtree(tree, node, super.tree = FALSE)
}
\arguments{
\item{tree}{A vector of character strings, each containing a newick tree}

\item{node}{A vector of character strings, each representing the name
of the node which will be the new root of the tree.  If node is shorter
than tree, values will be recycled, and a warning produced if \code{length(tree) \%\% length(node) != 0}}

\item{super.tree}{A vector of logical values.  If TRUE, then remove all
nodes which are descendants of node, rather than keeping them.}
}
\value{
A vector of trees which have been pruned, removing all nodes
which are not descendants of the given node.
}
\description{
Get a subtree
}
\examples{
trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385)
            :0.020666,canFam2:0.193569);")
trees <- name.ancestors(trees)
subtree(trees, c("hg18-panTro2", "mm9-rn4"))
}
\author{
Melissa J. Hubisz and Adam Siepel
}
\keyword{trees}
