% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hmm.R
\name{reflect.phylo.hmm}
\alias{reflect.phylo.hmm}
\title{Reflect a phylo-hmm across a strand}
\usage{
reflect.phylo.hmm(x, pivot.states, mods = NULL)
}
\arguments{
\item{x}{An object of type hmm}

\item{pivot.states}{The list of states to "reflect" across; these
should be the states that are not strand-specific.  Can be an
integer vector containing state indices, or a character vector
corresponding to state names (in \code{row.names(x$trans.mat)})}

\item{mods}{A list of objects of type \code{tm} representing phylogenetic
models corresponding to each state in the hmm.  If given, then the
models will also be reflected and the return value will be a list with
a new hmm and a new list of models.}
}
\value{
If \code{mods==NULL} then a new hmm will be returned.  Otherwise
a list containing the new hmm and the corresponding models will be
returned.
}
\description{
Reflect a phylo-hmm across a strand
}
\examples{
#state.names <- c("neutral", "conserved", "codon1", "codon2", "codon3")
#h <- hmm(t(matrix(c(0.95, 0.04, 0.01, 0, 0, 
#                    0.04, 0.95, 0.01, 0, 0,
#                    0,    0,    0,    1, 0,
#                    0,    0,    0,    0, 1,
#                    0.005, 0.005,  0.99,   0, 0), nrow=5,
#                  dimnames=list(state.names, state.names))))
#         eq.freq=c(0.6, 0.3, 0.1/3, 0.1/3, 0.1/3))
#reflect.phylo.hmm(h, c("neutral", "conserved"))

}
\author{
Melissa J. Hubisz and Adam Siepel
}
