# don't export in the future, 
##'
##'
##' @export
phastCons.call <- function(msa,
                           mod,
                           rho=0.3,
                           estimate.trees=FALSE,
                           estimate.rho=FALSE,
                           gc=NULL,
                           nrates=NULL,
                           transitions=NULL,
                           init.transitions=NULL,
                           target.coverage=NULL,
                           expected.length=NULL,
                           init.expected.length=NULL,
                           viterbi=FALSE,
                           score.viterbi=FALSE,
                           compute.lnl=FALSE,
                           suppress.probs=FALSE,
                           ref.idx=1) {
  # check parameters
  check.arg(rho, "rho", "numeric", null.OK=FALSE)
  if (rho < 0 || rho > 1) stop("rho should be in range [0,1]")
  check.arg(estimate.trees, "estimate.trees", "logical", null.OK=FALSE)
  check.arg(estimate.rho, "estimate.rho", "logical", null.OK=FALSE)
  check.arg(gc, "gc", "numeric", null.OK=TRUE)
  if (!is.null(gc) && (gc < 0 || gc > 1)) stop("gc should be in range [0,1]")
  check.arg(nrates, "nrates", "integer", null.OK=TRUE)
  if (!is.null(nrates) && nrates <=0) stop("nrates should be >=1")
  check.arg(transitions, "transitions", "numeric", min.length=1L, max.length=2L,
            null.OK=TRUE)
  if (!is.null(transitions) && (sum(transitions > 1)>0 ||
                                sum(transitions < 0)>0))
    stop("transitions values should be in range [0,1]")
  if (!is.null(transitions) && length(transitions)==1L)
    transitions <- rep(transitions, 2)
  check.arg(init.transitions, "init.transitions", "numeric", min.length=1L,
            max.length=2L, null.OK=TRUE)
  if (!is.null(init.transitions) && (sum(init.transitions > 1)>0 ||
                                     sum(init.transitions < 0)>0))
    stop("init.transitions values should be in range [0,1]")
  if (!is.null(init.transitions) && length(init.transitions)==1L)
    init.transitions <- rep(init.transitions, 2)
  check.arg(target.coverage, "target.coverage", "numeric", null.OK=TRUE)
  if (!is.null(target.coverage) && (target.coverage <= 0 ||
                                    target.coverage >=1))
    stop("target.coverage should be in range (0,1)")
  check.arg(expected.length, "expected.length", "numeric", null.OK=TRUE)
  if (!is.null(expected.length) && expected.length <= 0)
    stop("expected.length should be greater than 0")
  check.arg(init.expected.length, "init.expected.length", "numeric",
            null.OK=TRUE)
  if (!is.null(init.expected.length) && init.expected.length <= 0)
    stop("init.expected.length should be greater than 0")
  check.arg(viterbi, "viterbi", "logical", null.OK=FALSE)
  check.arg(score.viterbi, "score.viterbi", "logical", null.OK=FALSE)
  check.arg(compute.lnl, "compute.lnl", "logical", null.OK=FALSE)
  check.arg(suppress.probs, "suppress.probs", "logical", null.OK=FALSE)
  check.arg(ref.idx, "ref.idx", "integer", null.OK=FALSE)
  if (ref.idx < 0) stop("ref.idx must be >=0")
  if (ref.idx > nrow.msa(msa)) stop("ref.idx must be <= nrow.msa(msa)")

  # check for bad param combinations
  if (estimate.trees && estimate.rho)
    stop("cannot specify both estimate.trees and estimate.rho")
  if (!is.null(gc) && !estimate.trees && !estimate.rho)
    stop("can only specify gc if estimate.trees or estimate.rho")
  if (!is.null(transitions) && !is.null(init.transitions))
    stop("cannot specify both transitions and init.transitions")
  if ((!is.null(transitions) || !is.null(init.transitions)) &&
      (!is.null(target.coverage) || !is.null(expected.length) ||
       !is.null(init.expected.length)))
    stop("cannot specify transitions or init.transitions with target.coverage, expected.length, or init.expected.length")
  if (score.viterbi && !viterbi)
    stop("cannot specify score.viterbi without viterbi")

  
  if (is.null(mod$tree)) {
    l <- length(mod$tree)
    modList <- list()
    for (i in 1:l) {
      if (is.null(mod[[l]]$tree))
        stop("invalid tree model")
      modList[[i]] <- as.pointer.tm(mod[[l]])$externalPtr
    }
  } else {
    modList <- as.pointer.tm(mod)
  }

  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  
  result <- .Call("rph_phastCons", msa$externalPtr, modList,
                  rho, estimate.trees, estimate.rho,
                  gc, nrates, transitions, init.transitions,
                  target.coverage, expected.length,
                  init.expected.length, viterbi,
                  score.viterbi, compute.lnl, suppress.probs,
                  ref.idx)
  rphast.simplify.list(result)
}



##' Produce conservation scores and identify conserved elements,
##' given a multiple alignment and a phylo-HMM.
##'
##' A phylo-HMM consisting of two states is assumed: a "conserved"
##' state and a "non-conserved" state.  If two phylogenetic models
##' are given, the first is the conserved state, and the second
##' is the non-conserved state.  If only one model is given, then
##' this is used as the non-conserved state, and the conserved state
##' is obtained by multiplying the branch lengths by the parameter
##' rho.
##'
##' @param msa An object of type \code{msa} representing the multiple
##' alignment to be scored.
##' @param mod Either a single object of type \code{tm}, or a list
##' containing two \code{tm} objects.  If two objects are given, they
##' represent the conserved and non-conserved phylogenetic models.  If
##' one is given, then this represents the non-conserved model, and the
##' conserved model is obtained by scaling the branches by a parameter
##' rho.
##' @param rho The scaling parameter
##' @export
phastCons <- function(msa,
                      mod,
                      rho=0.3,
                      estimate.trees=FALSE,
                      estimate.rho=FALSE,
                      gc=NULL,
                      nrates=NULL,
                      transitions=NULL,
                      init.transitions=NULL,
                      target.coverage=NULL,
                      expected.length=NULL,
                      init.expected.length=NULL,
                      viterbi=FALSE,
                      ref.idx=1) {
  score.viterbi <- viterbi;
  compute.lnl <- TRUE
  suppress.probs <- FALSE
  phastCons.call(msa,mod, rho=rho, estimate.trees=estimate.trees,
                 estimate.rho=estimate.rho, gc=gc, nrates=nrates,
                 transitions=transitions,
                 init.transitions=init.transitions,
                 target.coverage=target.coverage,
                 expected.length=expected.length,
                 init.expected.length=init.expected.length,
                 viterbi=viterbi, score.viterbi=score.viterbi,
                 compute.lnl=compute.lnl,
                 suppress.probs=suppress.probs,
                 ref.idx=ref.idx)
}
