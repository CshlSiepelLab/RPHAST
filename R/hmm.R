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
nstate.hmm <- function(hmm) {
  stopIfNotValidHmm(hmm)
  nrow(hmm$trans.mat)
}


##' @export
as.pointer.hmm <- function(hmm) {
  obj <- list()
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
  hmm <- .makeObj.tm()
  hmm$trans.mat = .Call("rph_hmm_transMat", x$externalPtr)
  hmm$eq.freq = .Call("rph_hmm_eqFreq", x$externalPtr)
  hmm$begin.freq = .Call("rph_hmm_beginFreq", x$externalPtr)
  temp <- .Call("rph_hmm_endFreq", x$externalPtr)
  if (!is.null(temp))
    hmm$end.freq <- temp
  hmm
}
