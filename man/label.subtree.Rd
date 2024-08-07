% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees.R
\name{label.subtree}
\alias{label.subtree}
\title{Label subtree}
\usage{
label.subtree(tree, node, label, include.leading = FALSE)
}
\arguments{
\item{tree}{A vector of character strings, each containing a newick tree}

\item{node}{A character string, giving the node at the head of the subtree.}

\item{label}{A single character string giving the label to apply
to the branches in the subtree.}

\item{include.leading}{A logical value; if \code{TRUE}, include the
branch leading to the subtree in the labelled group; otherwise include
only descendants of the named node.}
}
\value{
A vector of character strings containing the modified trees;
the branches are labelled by appending a pound sign and the label
to the node name in the tree string.
}
\description{
Apply a label to a subtree
}
\examples{
trees <- c("((hg18:1.0, panTro2:2.0)hg18-panTro2:3.0, mm9:4.0);",
           "(((hg18:0.01, panTro2:0.01)hg18-panTro2:0.07,
              (mm9:0.083220,rn4:0.090564)mm9-rn4:
             0.269385)hg18-rn4:0.020666,canFam2:0.193569);")
label.subtree(trees, "hg18-panTro2", "human-chimp", include.leading=FALSE)
label.subtree(trees, "hg18-panTro2", "human-chimp", include.leading=TRUE)
}
\author{
Melissa J. Hubisz
}
\keyword{trees}
