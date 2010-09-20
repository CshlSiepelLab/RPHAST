# make barebones hmm obj
##' @nord
##' @export
.makeObj.hmm <- function() {
  hmm <- list()
  class(hmm) <- "hmm"
  hmm
}

##' @nord
##' @export
fixFreq.hmm <- function(freq, name, n) {
  if (is.null(freq)) {
    freq <- rep(1/n, n)
  } else {
    freq <- as.numeric(freq)/sum(freq)
    if (length(freq) != n)
      stop(name, " should have same length as trans.mat dimensions (", n,")")
  }
  freq
}

##' Create a new HMM object
##' @title Create an rphast HMM object
##' @param trans.mat A square matrix object of dimension n x n where n is the
##' number of states, and element [i,j] is the rate of
##' transition from state i to state j
##' @param eq.freq A vector of length n giving the equilibrium
##' frequencies of each state.  If NULL, assume a uniform distribution
##' across states.
##' @param begin.freq A vector of length n giving the initial state
##' frequencies.  If NULL, assume a uniform distribution across states.
##' @param end.freq A vector of length n giving the final state frequencies.
##' If NULL, do not condition on end frequencies.
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
hmm <- function(trans.mat, eq.freq=NULL, begin.freq=NULL,
                end.freq=NULL) {
  trans.mat <- as.matrix(trans.mat)
  n <- nrow(trans.mat)
  if (ncol(trans.mat) != n)
    stop("trans.mat should be square")
  eq.freq <- fixFreq.hmm(eq.freq, "eq.freq", n)
  begin.freq <- fixFreq.hmm(begin.freq, "begin.freq", n)
  if (! is.null(end.freq))
    end.freq <- fixFreq.hmm(end.freq, "end.freq", n)
  for (i in 1:n) trans.mat[i,] <- fixFreq.hmm(trans.mat[i,], "trans.mat", n)
  hmm <- .makeObj.hmm()
  hmm$trans.mat <- trans.mat
  hmm$eq.freq <- eq.freq
  hmm$begin.freq <- begin.freq
  if (! is.null(end.freq))
    hmm$end.freq <- end.freq
  hmm
}

##' Read an HMM object from a file
##'
##' This function uses phast's internal hmm format, which is quite
##' simple.  See \code{write.hmm} or file used in example below for
##' exaples of hmm format.
##' @param filename The file to read
##' @return An hmm object
##' @export
##' @keywords hmm
##' @author Melissa J. Hubisz and Adam Siepel
read.hmm <- function(filename) {
  h <- .makeObj.hmm()
  h$externalPtr <- .Call("rph_hmm_new_from_file", filename)
  from.pointer.hmm(h)
}


##' Write an HMM object to a file
##' @param x An object of type \code{hmm}
##' @param filename The name of the file to write to (if NULL, write
##' to terminal)
##' @param append If \code{TRUE}, append hmm to existing file, otherwise
##' overwrite.
##' @export
##' @keywords hmm
##' @author Melissa J. Hubisz and Adam Siepel
write.hmm <- function(x, filename, append=FALSE) {
  h <- as.pointer.hmm(x)
  invisible(.Call("rph_hmm_print", h$externalPtr, filename, append))
}


##' @export
stopIfNotValidHmm <- function(hmm){
  if (is.null(hmm$trans.mat) ||
      is.null(hmm$eq.freq) ||
      is.null(hmm$begin.freq) ||
      nrow(hmm$trans.mat) != ncol(hmm$trans.mat) ||
      nrow(hmm$trans.mat) != length(hmm$eq.freq) ||
      nrow(hmm$trans.mat) != length(hmm$begin.freq) ||
      ((! is.null(hmm$end.freq)) &&
       nrow(hmm$trans.mat != length(hmm$end.freq))))
    stop("invalid hmm object")
  invisible(NULL)
}

##' HMM number of states
##' @param hmm An object of type \code{hmm}
##' @return The number of states in the hidden Markov Model
##' @export
##' @keywords hmm
##' @author Melissa J. Hubisz
nstate.hmm <- function(hmm) {
  stopIfNotValidHmm(hmm)
  nrow(hmm$trans.mat)
}


##' @export
as.pointer.hmm <- function(hmm) {
  obj <- .makeObj.hmm()
  obj$externalPtr <- .Call("rph_hmm_new", hmm$trans.mat, hmm$eq.freq,
                           hmm$begin.freq, hmm$end.freq)
  obj
}

##' @nord
##' @export
from.pointer.hmm <- function(x) {
  if (is.null(x$externalPtr)) {
    stopIfNotValidHmm()
    return(x)
  }
  hmm <- .makeObj.hmm()
  hmm$trans.mat = .Call("rph_hmm_transMat", x$externalPtr)
  hmm$eq.freq = .Call("rph_hmm_eqFreq", x$externalPtr)
  hmm$begin.freq = .Call("rph_hmm_beginFreq", x$externalPtr)
  temp <- .Call("rph_hmm_endFreq", x$externalPtr)
  if (!is.null(temp))
    hmm$end.freq <- temp
  rphast.simplify.list(hmm)
}


##' Produce likelihood of an alignment given a phylo-HMM, posterior
##' probabilities of phylo-HMM states across an alignment,
##' and predict states using Viterbi algorithm
##' @title Score an alignment using a general phylo-HMM
##' @param msa An object of type \code{msa}
##' @param mod A list of tree model objects, corresponding to each state in the phylo-HMM
##' @param hmm An object of type \code{hmm} describing transitions between states,
##' equilbrium frequencies, initial frequencies, and optionally end frequencies
##' @param states A vector of characters naming
##' the states of interest in the phylo-HMM, or a vector of integers
##' corresponding to states in the transition matrix.  The post.probs will give
##' the probability of any of these states, and the viterbi regions reflect
##' regions where the state is predicted to be any of these states.
##' @param viterbi A logical value indicating whether to predict a path through the phylo-HMM
##' using the Viterbi algorithm.
##' @param ref.idx An integer value.  Use the coordinate frame of the given sequence.
##' Default is 1, indicating the first sequence in the alignment.
##' A value of 0 indicates the coordinate frame of the entire alignment.
##' @param reflect.strand Given an hmm describing
##' the forward strand, create a larger HMM that allows for features
##' on both strands by "reflecting" the original HMM about the specified
##' states.  States can be described as a vector of integers or characters
##' in the same manner as states argument (above).  The new hmm will be
##' used for prediction on both strands.
##' @param features If non-NULL, compute the likelihood of each feature
##' under the phylo-HMM.
##' @param quiet If \code{TRUE}, suppress printing of progress information.
##' @return If \code{features} is not NULL, returns a numeric vector
##' with one value per feature, giving the likelihood of the feature under
##' the phylo-HMM.
##'
##' Otherwise, returns a list with some or all of
##' the following arguments (depending on options):
##' \item{in.states}{An object of type \code{feat} which describes regions which
##' fall within the interesting states specified in the states parameter,
##' as determined by the Viterbi algorithm.}
##' \item{post.prob.wig}{A data frame giving a coordinate and posterior
##' probibility that each site falls within an interesting state.}
##' \item{likelihood}{The likelihood of the data under the estimated model.}
##' @export
##' @keywords hmm
##' @author Melissa J. Hubisz and Adam Siepel
score.hmm <- function(msa, mod, hmm, states, viterbi=TRUE, ref.idx=1,
                      reflect.strand=NULL, features=NULL,
                      quiet=(!is.null(features))) {
  if (!is.null(features)) {
    if (!is.data.frame(features)) {
      features <- as.data.frame.feat(features)
      if (!is.data.frame(features)) stop("invalid features")
    }
    rv <- numeric(nrow(features))
    if (!is.null(msa$externalPtr)) 
      msa <- as.pointer.msa(msa)
    for (i in 1:nrow(features)) {
      m <- extract.feature.msa(copy.msa(msa), features[i,], pointer.only=TRUE)
      pcResult <- phastCons.call(msa=m, mod,
                                 rho=NULL, target.coverage=NULL, expected.length=NULL, transitions=NULL,
                                 estimate.rho=FALSE, estimate.expected.length=FALSE,
                                 estimate.transitions=FALSE,
                                 estimate.trees=FALSE, viterbi=FALSE, score.viterbi=FALSE,
                                 gc=NULL, nrates=NULL, compute.lnl=TRUE, suppress.probs=FALSE,
                                 ref.idx=ref.idx, hmm=hmm, states=states,
                                 reflect.strand=reflect.strand, quiet=quiet)
      rv[i] <- pcResult$likelihood
    }
  } else {
    rv <- phastCons.call(msa, mod,
                         rho=NULL, target.coverage=NULL, expected.length=NULL, transitions=NULL,
                         estimate.rho=FALSE, estimate.expected.length=FALSE, estimate.transitions=FALSE,
                         estimate.trees=FALSE,
                         viterbi=viterbi, score.viterbi=viterbi,
                         gc=NULL, nrates=NULL, compute.lnl=TRUE, suppress.probs=FALSE,
                         ref.idx=ref.idx,
                         hmm=hmm,
                         states=states,
                         reflect.strand=reflect.strand, quiet=quiet)
    if (!is.null(rv$most.conserved)) {
      w <- which(names(rv) == "most.conserved")
      names(rv)[w] <- "in.states"
    }
  }
  rv
}
