# don't export in the future, 
##' @nord
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
                           ref.idx=1,
                           hmm=NULL, states=NULL, reflect.strand=NULL) {
  # check parameters
  check.arg(rho, "rho", "numeric", null.OK=FALSE)
  if (rho < 0 || rho > 1) stop("rho should be in range [0,1]")
  check.arg(estimate.trees, "estimate.trees", "logical", null.OK=FALSE)
  check.arg(estimate.rho, "estimate.rho", "logical", null.OK=FALSE)
  check.arg(gc, "gc", "numeric", null.OK=TRUE)
  if (!is.null(gc) && (gc < 0 || gc > 1)) stop("gc should be in range [0,1]")
  check.arg(nrates, "nrates", "integer", null.OK=TRUE, min.length=1L, max.length=2L)
  if (!is.null(nrates) || sum(nrates <= 0) >= 1L) stop("nrates should be >=1")
  check.arg(transitions, "transitions", "numeric", min.length=1L, max.length=2L,
            null.OK=TRUE)
  if (!is.null(transitions) && (sum(transitions > 1)>0L ||
                                sum(transitions < 0)>0L))
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
  if (!is.null(states)) {
    if (is.null(hmm)) {
      warning("states is not used when hmm is NULL")
    } else {
      if (is.character(states)) {
        orig <- states
        states <- sapply(states, function(x, mat) {
          which(row.names(mat) == x)}, hmm$trans.mat)
        if (length(states) != length(orig))
          warning("some states not found in hmm")
      } else {
        nstate <- nrow(hmm$trans.mat)
        check.arg(states, "states", "integer", null.OK=FALSE,
                  min.length=1L, max.length=nstate)
        if (sum(states <= 0 | states > max.length) > 0L)
          stop("invalid integers in states")
      }
    }
  }
  if (!is.null(reflect.strand)) {
    if (is.null(hmm)) {
      warning("reflect.strand not used when hmm is NULL")
    } else {
      if (is.character(reflect.strand)) {
        orig <- reflect.strand
        reflect.strand <- sapply(reflect.strand, function(x, mat) {
          which(row.names(mat) == x)}, hmm$trans.mat)
        if (length(reflect.strand) != length(orig))
          warning("some reflect.strand elements not found in hmm")
      } else {
        nstate <- nrow(hmm$trans.mat)
        check.arg(reflect.strand, "reflect.strand", "integer", null.OK=FALSE,
                  min.length=1L, max.length=nstate)
        if (sum(reflect.strand <= 0 | reflect.strand > max.length) > 0L)
          stop("invalid integers in reflect.strand")
      }
    }
  }
  if (!is.null(hmm)) hmm <- as.pointer.hmm(hmm)
  if (is.null(states) && (!is.null(hmm)) && (!suppress || viterbi))
    warning("no states given; scoring state 2")

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
    l <- length(mod)
    modList <- list()
    for (i in 1:l) {
      if (is.null(mod[[i]]$tree))
        stop("invalid tree model")
      modList[[i]] <- as.pointer.tm(mod[[i]])$externalPtr
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
                  ref.idx,
                  if(is.null(hmm)) NULL else hmm$externalPtr,
                  states, reflect.strand)
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
##' @param rho Set the scale (overall evolutionary rate) of the model for the
##' conserved state to be <rho> times tht of the model for the non-conserved state
##' ( 0 < rho < 1).  If used with estimate.trees or estimate.rho, the specified
##' value will be used for initialization only, and rho will be estimated.  This
##' argument is ignored if mod contains two tree model objects.
##' @param estimate.trees A logical value.  If \code{TRUE}, estimate free
##' parameters of tree models for conserved and non-conserved state.
##' Estimated models are stored in return list.
##' @param estimate.rho A logical value.  If \code{TRUE}, Estimate the parameter
##' rho (as described above).  Estimated value is reported in return list.
##' @param gc A single numeric value given the fraction of bases that are G or C, for
##' optional use with estimate.trees or estimate.rho.  This overrides the default
##' behavior of estimating the base composition empirically from the data.
##' @param nrates An integer vector of length one or two, for optional use with
##' estimate.trees and a discrete-gamma model.  Assume the specified number of
##' rate categories, rather than the number given in the input tree model(s).  If
##' two values are given they apply to the conserved and nonconserved models,
##' respectively.
##' @param transitions A numeric vector of length one or two, representing the
##' transition probabilities for the two-state HMM.  If not provided transition rates
##' will be estimated by maximum likelihood.  The first value represents mu, the
##' transition rate from the conserved to the non-conserved state, and the second
##' value is nu, the rate from non-conserved to conserved.  If only one value is
##' provided then mu=nu.  The rate of self-transitions are then 1-mu and 1-nu, and
##' the expected lengths of conserved and non-conserved elements are 1/mu and 1/nu,
##' respectively.
##' @param init.transitions Like the transitions argument (described above), but
##' only specifies initial values for mu, and nu.  Transitions probabilities will
##' then be estimated by maximum likelihood.
##' @param target.coverage A single numeric value, representing the parameter gamma,
##' which describes the fraction of sites in conserved elements.  This argument is
##' an alternative to the transitions argument.  This argument sets a prior
##' expectation rather than a posterior and assumes stationarity of the
##' state-transition process.  Adding this constraint causes the ratio of
##' between-state transitions to be fixed at (1-gamma)/gamma.  If used with the
##' expected.length argument, the transition probailities will be completely fixed,
##' otherwise the expected-length parameter will be estimated by maximum likelihood.
##' @param expected.length A single numeric value, representing the parameter omega,
##' which describes the expected length of conserved elements.  This is an
##' alternative to the transitions argument.  If provided with target.coverage, than
##' transition rates are fully determined, otherwise the target-coverage parameter
##' will be estimated by maximum likelihood.
##' @param init.expected.length Like expected.length above, but only sets the initial
##' value of the expected length parameter, which will be estimated by maximum
##' likelihood.
##' @param viterbi A logical value.  If \code{TRUE}, produce discrete elements
##' using the Viterbi algorithm.
##' @param ref.idx An integer value.  Use the coordinate frame of the given sequence.
##' Default is 1, indicating the first sequence in the alignment.
##' A value of 0 indicates the coordinate frame of the entire alignment.
##' @param hmm An object of type \code{hmm} describing a custom HMM
##' @param states (For use with hmm) A vector of characters naming
##' the states of interest in the phylo-HMM, or a vector of integers
##' corresponding to states in the transition matrix.  The post.probs will give
##' the probability of any of these states, and the viterbi regions reflect
##' regions where the state is predicted to be any of these states.
##' @param reflect.strand (For use with hmm) Given an hmm describing
##' the forward strand, create a larger HMM that allows for features
##' on both strands by "reflecting" the original HMM about the specified
##' states.  States can be described as a vector of integers or characters
##' in the same manner as states argument (above).  The new hmm will be
##' used for prediction on both strands.
##' @return A list containing parameter estimates.  The list may have any of the
##' following elements, depending on the arguments:
##' \item{transition.rates}{A numeric vector of length two giving the rates from the
##' conserved to the non-conserved state, and from the non-conserved to the conserved
##' state.}
##' \item{rho}{The relative evolutionary rate of the conserved state compared to the
##' non-conserved state.}
##' \item{tree.models}{Tree model objects describing the evolutionary process in the
##' conserved and non-conserved states.}
##' \item{most.conserved}{An object of type gff which describes conserved elements
##' detected by the Viterbi algorithm.}
##' \item{post.prob.wig}{A data frame giving a coordinate and score for individual
##' bases in the alignment}
##' \item{likelihood}{The likelihood of the data under the estimated model.}
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
                      ref.idx=1, hmm=NULL,
                      states=NULL,
                      reflect.strand=NULL) {
  score.viterbi <- viterbi
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
                 ref.idx=ref.idx, hmm=hmm, states, reflect.strand)
}
