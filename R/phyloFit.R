

#  Note: the man page for this function has become too complex for
# roxygen to currently handle; it is created manually and copied
# to RPHAST/man after roxygen is called.
##' @nord
##' @export
phyloFit <- function(msa,
                     tree=NULL,
                     subst.mod="REV",
                     init.mod=NULL,
                     no.opt=c("backgd"),
                     init.backgd.from.data=ifelse(is.null(init.mod),TRUE,FALSE),
                     features=NULL,
#                     do.cats=NULL,
#                     reverse.groups=NULL,
                     scale.only=FALSE,
                     scale.subtree=NULL,
                     nrates=1,
                     alpha=1,
                     rate.constants=NULL,
                     selection=NULL,
                     init.random=FALSE,
                     init.parsimony=FALSE,
                     clock=FALSE,
                     EM=FALSE,
                     precision="HIGH",
                     ninf.sites=50,
                     quiet=FALSE,
                     bound=NULL,
                     log.file=FALSE
#                     alt.mod=NULL  # don't use alt.mod here anymore, instead require init.mod with alt.mod defined for alt.mod
#                     symmetric.freqs=FALSE,
#                     report.error=FALSE,
#                     ancestor=NULL,
                     ) {
  if (nrow.msa(msa) > 2 && is.null(tree) && is.null(init.mod))
    stop("need tree if MSA has more than two sequences")
  if (!is.subst.mod.tm(subst.mod))
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
  check.arg(selection, "selection", "numeric", null.OK=TRUE)
  check.arg(no.opt, "no.opt", "character", null.OK=TRUE, min.length=1L,
            max.length=NULL)
  check.arg(init.random, "init.random", "logical", null.OK=FALSE)
  check.arg(init.parsimony, "init.parsimony", "logical", null.OK=FALSE)
  check.arg(init.backgd.from.data, "init.backgd.from.data",
            "logical", null.OK=FALSE,
            min.length=1L, max.length=1L)
  if (init.backgd.from.data == FALSE && is.null(init.mod))
    stop("init.backgd.from.data cannot be FALSE unless init.mod is provided")
  check.arg(clock, "clock", "logical", null.OK=FALSE)
  check.arg(EM, "EM", "logical", null.OK=FALSE)
  check.arg(precision, "precision", "character", null.OK=FALSE)

  check.arg(ninf.sites, "ninf.sites", "integer", null.OK=TRUE)
  check.arg(quiet, "quiet", "logical", null.OK=TRUE)
  # phyloFit will complain if precision is an invalid string, no need to
  # check here.
  check.arg(bound, "bound", "character", null.OK=TRUE)

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
  if (!is.null(init.mod))
    init.mod <- as.pointer.tm(init.mod)

  rphast.simplify.list(.Call.rphast("rph_phyloFit",
                                    msa$externalPtr,
                                    tree,
                                    subst.mod,
                                    scale.only,
                                    scale.subtree,
                                    nrates,
                                    alpha,
                                    rate.constants,
                                    if (is.null(init.mod)) NULL else init.mod$externalPtr,
                                    init.backgd.from.data,
                                    init.random,
                                    init.parsimony,
                                    clock,
                                    EM,
                                    precision,
                                    features$externalPtr,
                                    ninf.sites,
                                    quiet,
                                    no.opt,
                                    bound,
                                    log.file,
                                    selection))
}
