

##' Fit a Phylogenetic model to an alignment
##' @param msa An alignment object.  May be altered if passed in as a pointer to
##' C memory (see Note).
##' @param tree A character string containing a Newick formatted tree
##' defining the topology.  Required if the number of species > 3, unles
##' init.mod is specified.  The topology must be rooted, although the
##' root is ignored if the substitution model is reversible.
##' @param subst.mod The substitution model to use.  See
##' \code{\link{subst.mods}}.
##' @param init.mod An object of class \code{tm} used to initialize the model
##' @param no.freqs (Only applies when init.mod provided). If \code{TRUE}, do
##' not estimate equilibrium frequencies; just use the ones from init.mod.
##' @param no.rates (Only applies when init.mod provided). If \code{TRUE},
##' do not estimate transition rate parameters; just use the transition
##' matrix in init.mod.
##' @param features An object of type \code{feat}.  If given, a separate model will be
##' estimated for each feature type.
##' @param scale.only A logical value. If \code{TRUE}, estimate only the
##' scale of the tree.  Branches will be held at initial values.  Useful in
##' conjunction with init.mod.
##' @param scale.subtree A character string giving the name of a node in
##' a tree.  This option implies scale.only=TRUE.  If given, estimate
##' separate scale factors for subtree beneath identified node and the rest
##' of the tree.  The branch leading to the subtree is included in the subtree.
##' @param nrates An integer.  The number of rate categories to use.
##' Specifying a value greater than one causes the discrete gamma model for
##' rate variation to be used,  unless rate constants are specified.
##' @param alpha A numeric value > 0, for use with "nrates".  Initial value
##' for alpha, the shape parameter of the gamma distribution.
##' @param rate.constants A numeric vector.  Implies nrates =
##' length(rate.constants).  Also implies EM=TRUE.  Uses a non-parametric
##' mixture model for rates, instead of a gamma distribution.  The weight
##' associated with each rate will be estimated.  alpha may still be used to
##' initialize these weights.
##' @param init.random A logical value.  If \code{TRUE}, parameters will be
##' initialized randomly.
##' @param init.parsimony A logical value.  If \code{TRUE}, branch lengths
##' will be estimated based on parsimony counts for the alignments.
##' Currently only works for models of order 0.
##' @param clock A logical value.  If \code{TRUE}, assume a molecular clock
##' in estimation.
##' @param EM A logical value.  If \code{TRUE}, the model is fit using EM
##' rather than the default BFGS quasi-Newton algorithm.  Not available
##' for all models/options.
##' @param precision A character vector, one of "HIGH", "MED", or "LOW",
##' denoting the level of precision to use in estimating model parameters.
##' Affects convergence criteria for iterative algorithms: higher precision
##' means more iterations and longer execution time.
##' @param ninf.sites An integer.  Require at least this many "informative"
##' sites in order to estimate a model.  An informative site as an alignment
##' column with at least two non-gap and non-missing-data characers.
##' @param quiet A logical value.  If \code{TRUE}, do not report progress
##' to screen.
##' @return An object of class \code{tm} (tree model), or (if several models
##' are computed, as is possible with the features or windows options), a list of
##' objects of class \code{tm}.
##' @note If msa or features object are passed in as pointers to C memory,
##' they may be altered by this function!  Use \code{copy.msa(msa)} or
##' \code{copy.feat(features)} to avoid this behavior!
##' @keywords msa tm features trees
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
phyloFit <- function(msa,
                     tree=NULL,
                     subst.mod="REV",
                     init.mod=NULL,
                     no.freqs=FALSE,
                     no.rates=FALSE,
                     features=NULL,
#                     do.cats=NULL,
#                     reverse.groups=NULL,
                     scale.only=FALSE,
                     scale.subtree=NULL,
                     nrates=1,
                     alpha=1,
                     rate.constants=NULL,
                     init.random=FALSE,
                     init.parsimony=FALSE,
                     clock=FALSE,
                     EM=FALSE,
                     precision="HIGH",
                     ninf.sites=50,
                     quiet=FALSE
#                     no.opt=NULL,
#                     init.freqs="from.data",
#                     symmetric.freqs=FALSE,
#                     report.error=FALSE,
#                     ancestor=NULL,
                     ) {
  if (nrow.msa(msa) > 2 && is.null(tree) && is.null(init.mod))
    stop("need tree if MSA has more than two sequences")
  if (!isSubstMod.tm(subst.mod))
    stop("invalid substitution model ", subst.mod)
  check.arg(tree, "tree", "character", null.OK=TRUE)
  check.arg(scale.only, "scale.only", "logical", null.OK=FALSE)
  check.arg(scale.subtree, "scale.subtree", "character", null.OK=TRUE)
  check.arg(nrates, "nrates", "integer", null.OK=FALSE)
  check.arg(alpha, "alpha", "numeric", null.OK=FALSE)
  if (!is.null(alpha) && alpha <= 0)
    stop("alpha must be > 0")
  check.arg(rate.constants, "rate.constants", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  if (!is.null(init.mod))
    init.mod <- as.pointer.tm(init.mod)
  check.arg(no.freqs, "no.freqs", "logical", null.OK=FALSE)
  check.arg(no.rates, "no.rates", "logical", null.OK=FALSE)
  if (no.freqs && is.null(init.mod)) {
    warning("no.freqs only applies when init.mod is specified; ignoring")
    no.freqs <- FALSE
  }
  if (no.rates && is.null(init.mod)) {
    warning("no.rates only applies when init.mod is specified; ignoring")
    no.rates <- FALSE
  }
  check.arg(init.random, "init.random", "logical", null.OK=FALSE)
  check.arg(init.parsimony, "init.parsimony", "logical", null.OK=FALSE)
  check.arg(clock, "clock", "logical", null.OK=FALSE)
  check.arg(EM, "EM", "logical", null.OK=FALSE)
  check.arg(precision, "precision", "character", null.OK=FALSE)

  check.arg(ninf.sites, "ninf.sites", "integer", null.OK=TRUE)
  check.arg(quiet, "quiet", "logical", null.OK=TRUE)
  # phyloFit will complain if precision is an invalid string, no need to
  # check here.

  if (!is.null(rate.constants)) {
    EM <- TRUE
    nrates <- length(rate.constants)
  }

  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  if (!is.null(features)) {
    if (is.null(features$externalPtr))
      features <- as.pointer.feat(features)
  }

  result <- list()
  result$externalPtr <- .Call("rph_phyloFit",
                              msa$externalPtr,
                              tree,
                              subst.mod,
                              scale.only,
                              scale.subtree,
                              nrates,
                              alpha,
                              rate.constants,
                              init.mod$externalPtr,
                              no.freqs,
                              no.rates,
                              init.random,
                              init.parsimony,
                              clock,
                              EM,
                              precision,
                              features$externalPtr,
                              ninf.sites,
                              quiet)

  #need to parse result to make a list in R
  numModels <- .Call("rph_phyloFit_result_num_models", result$externalPtr)
  if (numModels > 1L)
    resultList <- list()
  for (i in 1:numModels) {
    tm <- list()
    tm$externalPtr <- .Call("rph_phyloFit_result_get_model",
                            result$externalPtr, i)
    tm <- from.pointer.tm(tm)

    # need to incorporate error here if option was given

    if (numModels==1L) return(tm)  # no need for labels if result has only
                                   # one model
    currName <- .Call("rph_phyloFit_result_get_name",
                      result$externalPtr, i)
    if (is.null(currName)) currName <- i
    
    resultList[[currName]] <- tm
  }
  resultList
}



