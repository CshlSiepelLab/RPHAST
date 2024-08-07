% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeModel.R
\name{tm}
\alias{tm}
\title{Tree Models}
\usage{
tm(
  tree,
  subst.mod,
  rate.matrix = NULL,
  backgd = NULL,
  alphabet = "ACGT",
  nratecats = 1,
  alpha = 0,
  rate.consts = NULL,
  rate.weights = NULL,
  selection = NULL,
  root.leaf = NULL,
  likelihood = NULL
)
}
\arguments{
\item{tree}{A character string representing a phylogenetic tree in
newick foromat}

\item{subst.mod}{A character string giving a valid substitution mod.
See \code{\link{subst.mods}}.}

\item{rate.matrix}{A square matrix representing the rate of substitution
from one state to the next.}

\item{backgd}{A numeric vector giving the equilibrium frequencies for
each state.}

\item{alphabet}{A character vector containing all valid states, given
in the order they are represented in rate.matrix and backgd.  Defaults
to "ACGT"}

\item{nratecats}{The number of rate categories in the model.  Defaults
to 1.}

\item{alpha}{If nratecats > 1, weight for each category is computed using
a gamma distribution with shape parameter alpha.}

\item{rate.consts}{The rate for each rate category.  NULL if only
one category.}

\item{rate.weights}{Vector of numeric of length nratecats, determining
the weight of each rate category.  Must sum to 1 (will be normalized
otherwise).  May be defined implicitly by alpha.}

\item{selection}{If not NULL, then this is a numeric value giving the
selection parameter for this model.  If NULL then there is no selection
in the model.  If selection==0.0, means that selection has no effect
in the current model, but is part of the model, and by default the
selection parameter will be optimized by phyloFit.  The rate matrix
is assumed to already be scaled by the selection parameter, if provided.}

\item{root.leaf}{Usually NULL, but if set to the name of a leaf
node in the tree, the tree will be re-rooted at this leaf node.}

\item{likelihood}{an optional value giving the log likelihood of this
model for some alignment.}
}
\value{
An object of class \code{tm} representing a phylogenetic model.
}
\description{
Make a new tree model
}
\details{
Tree models represent a substitution process along a phylogenetic
tree.  They are stored as a list, with components defined by the
arguments to this function.
}
\examples{
tree <- "((human:0.01, chimp:0.01):0.03, mouse:0.3)"
subst.mod <- "JC69"
rate.mat <- matrix(runif(16), nrow=4, ncol=4)
for (i in 1:4)
  rate.mat[i,i] <- -sum(rate.mat[i,-i])
backgd <- runif(4)
backgd <- backgd/sum(backgd)
alphabet <- "ACGT"
t <- tm(tree, subst.mod, rate.mat, backgd, alphabet)
t

nratecats <- 3
alpha <- 1.5
rate.consts <- runif(nratecats, max=3.0)
root.leaf <- "human"
t <- tm(tree, subst.mod, rate.matrix=rate.mat,
        backgd=backgd, alphabet=alphabet,
        nratecats=nratecats, alpha=alpha,
        rate.consts=rate.consts, root.leaf=root.leaf)
t

}
\author{
Melissa J. Hubisz and Adam Siepel
}
\keyword{tm}
