#' @nord
tm.makeObj <- function() {
  tm <- list()
  class(tm) <- "tm"
  tm
}


# this needs a special function to avoid having to make a matrix object
# in R from within C.  The C code returns a vector which is coerced
# into a matrix here.  Watch out for row/column order!
#' @nord
tm.rateMatrix.from.pointer <- function(x) {
  if (is.null(x$externalPtr))
    stop("tm.rateMatrix.from.pointer expects list with externalPtr")
  m <- .Call("rph_tm_rateMatrix", x$externalPtr)
  if (is.null(m)) return(NULL)
  matrix(m, nrow=sqrt(length(m)), byrow=TRUE)
}


#' @nord
tm.rootLeaf.from.pointer <- function(x, tree) {
  if (is.null(x$externalPtr))
    stop("tm.rootLeaf.from.pointer expects list with externalPtr")
  id <- .Call("rph_tm_rootLeaf", x$externalPtr)
  if (is.null(id) || id == -1) return(NULL)
  .Call("rph_tree_nodeName", tree, id)
}


#' @nord
tm.from.pointer <- function(x) {
  cat("tm.from.pointer\n")
  if (is.null(x$externalPtr))
    stop("tm.from.pointer expects list with externalPtr")
  tm <- tm.makeObj()
  tm$alphabet <- .Call("rph_tm_alphabet", x$externalPtr)
  tm$backgd <- .Call("rph_tm_backgd", x$externalPtr)
  tm$rate.matrix <- tm.rateMatrix.from.pointer(x)
  tm$subst.mod <- .Call("rph_tm_substMod", x$externalPtr)
  tm$alpha <- .Call("rph_tm_alpha", x$externalPtr)
  tm$nratecats <- .Call("rph_tm_nratecats", x$externalPtr)
  tm$rate.consts <- .Call("rph_tm_rK", x$externalPtr)
  tm$rate.weights <- .Call("rph_tm_freqK", x$externalPtr)
  tm$tree <- .Call("rph_tm_tree", x$externalPtr)
  tm$root.leaf <- tm.rootLeaf.from.pointer(x, tm$tree)
  tm
}


#' @nord
tm.free <- function(x) {
  invisible(.Call("rph_tm_free", x))
}


#' @nord
tm.to.pointer <- function(tm) {
  x <- list()
  x$externalPtr <- .Call("rph_tm_new",
                         tm$tree,
                         tm$alphabet,
                         tm$backgd,
                         tm$rate.matrix,
                         tm$subst.mod,
                         tm$alpha,
                         tm$nratecats,
                         tm$rate.consts,
                         tm$rate.weights,
                         tm$root.leaf)
  reg.finalizer(x$externalPtr, tm.free)
  x
}

##' Read a tree model from a file
##' @param filename The file containing a tree model
##' @title Read a Tree Model
##' @return An object of class "tm"
##' @seealso \code{\link{tm.new}}
##' @export
read.tm <- function(filename) {
  cat("read.tm\n")
  check.arg(filename, "filename", "character", null.OK=FALSE)
  x <- list()
  x$externalPtr <- .Call("rph_tm_read", filename)
  reg.finalizer(x$externalPtr, tm.free)
  tm <- tm.from.pointer(x)
  tm
}


##' Write a tree model to a file (or to the terminal)
##' @title Wrting Tree Models
##' @param tm An object of class "tm"
##' @param filename The filename to write to (use NULL for output to terminal)
##' @param append Whether to append the tree to the end of the file
##' (if FALSE, overwrites file).  Not used if filename is \code{NULL}
##' @seealso \code{\link{tm.new}}
##' @export
tm.write <- function(tm, filename=NULL, append=FALSE) {
  check.arg(filename, "filename", "character", null.OK=TRUE)
  check.arg(append, "append", "logical", null.OK=TRUE)
  tm <- tm.to.pointer(tm)
  invisible(.Call("rph_tm_print", tm$externalPtr, filename, append))
}


##' Tree model summary
##' @title Tree Model Summary
##' @param object An object of class tm
##' @param ... Parameters to be passed to/from other functions
##' @seealso \code{\link{tm.new}}
##' @export
summary.tm <- function(object, ...) {
  tm.write(object, NULL)
}


##' Coerce a tree model into a list
##' @title Tree Model to List
##' @param x an object of class tm
##' @param ... arguments to be passed to/from other functions
##' @seealso \code{\link{tm.new}}
##' @export
as.list.tm <- function(x, ...) {
  class(x) <- "list"
  x
}

##' Print a tree model
##' @title Printing Tree Models
##' @param x An object of class "tm"
##' @param ... arguments to be passed to/from other functions
##' @seealso \code{\link{tm.new}}
##' @export
print.tm <- function(x, ...) {
  print(as.list(x), ...)
}


##' Check if a string represents a phast substitution model
##' @title Check Substitution Model Strings
##' @param mod A vector of character strings representing substitution
##' model names to be tested
##' @return A vector of logical values indicating whether each string
##' represents a defined substitution model
##' @export
tm.isValidSubstMod <- function(mod) {
  result <- logical(length(mod))
  for (i in 1:length(mod)) 
    result[i] <- .Call("rph_subst_mods_is_valid_string", mod[i])
  result
}
    


##' Make a new tree model
##'
##' Tree models represent a substitution process along a phylogenetic
##' tree.  They are stored as a list, with components defined by the
##' arguments to this function.
##' @title Tree Models
##' @param tree A character string representing a phylogenetic tree in
##' newick foromat
##' @param subst.mod A character string giving a valid substitution mod.
##' See \code{\link{tm.isValidSubstMod}}.
##' @param rate.matrix A square matrix representing the rate of substitution
##' from one state to the next.
##' @param backgd A numeric vector giving the equilibrium frequencies for
##' each state.
##' @param alphabet A character vector containing all valid states, given
##' in the order they are represented in rate.matrix and backgd.  Defaults
##' to "ACGT"
##' @param nratecats The number of rate categories in the model.  Defaults
##' to 1.
##' @param alpha  If nratecats > 0, weight for each category is given
##' by a gamma distribution with shape parameter alpha
##' @param rate.consts The rate for each rate category.  NULL if only
##' one category
##' @param rate.weights Vector of numeric of length nratecats, determining
##' the weight of each rate category.  Must sum to 1 (will be normalized
##' otherwise).  May be defined implicitly by alpha.
##' @param root.leaf Usually NULL, but if set to the name of a leaf
##' node in the tree, the tree will be re-rooted at this leaf node.
##' @return An object of class tm
##' @export
tm.new <- function(tree, subst.mod, rate.matrix=NULL, backgd=NULL,
                   alphabet="ACGT", nratecats=1, alpha=0.0,
                   rate.consts=NULL, rate.weights=NULL,
                   root.leaf=NULL) {
  check.arg(tree, "tree", "character", null.OK=FALSE)
  check.arg(subst.mod, "subst.mod", "character", null.OK=FALSE)
  check.arg(rate.matrix, "rate.matrix", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL)  #will check later
  check.arg(backgd, "backgd", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL) # will check later
  check.arg(alphabet, "alphabet", "character", null.OK=FALSE)
  check.arg(nratecats, "nratecats", "integer", null.OK=FALSE)
  check.arg(alpha, "alpha", "numeric", null.OK=FALSE)
  check.arg(rate.consts, "rate.consts", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(root.leaf, "root.leaf", "character", null.OK=TRUE)

  if (!tm.isValidSubstMod(subst.mod))
    stop("invalid subst mod ", subst.mod)
  matsize <- NULL
  if (!is.null(rate.matrix)) {
    matsize <- as.integer(sqrt(length(rate.matrix)))
    if (matsize*matsize != length(rate.matrix))
      stop("rate.matrix should be square matrix")
    if (!is.matrix(rate.matrix))
      rate.matrix <- matrix(rate.matrix, nrow=matsize, ncol=matsize, byrow=TRUE)
  }
  if (!is.null(backgd)) {
    if (!is.null(matsize)) {
      if (matsize != length(backgd))
        stop("length of backgd should match dimensions of rate.matrix")
    }
  }
  if (!is.null(rate.consts)) {
    if (nratecats != length(rate.consts))
      stop("rate.consts should have length equal to nratecats")
  }
  if (!is.null(root.leaf)) {
    cat("calling .Call\n")
    if ( ! (.Call("rph_tree_isNode", tree, root.leaf))) {
      stop("tree has no node named ", root.leaf)
    }
    cat("done\n")
  }

  tm <- tm.makeObj()
  tm$alphabet <- alphabet
  tm$backgd <- backgd
  tm$rate.matrix <- rate.matrix
  tm$subst.mod <- subst.mod
  tm$alpha <- alpha
  tm$nratecats <- nratecats
  tm$rate.consts <- rate.consts
  tm$tree <- tree
  tm$root.leaf <- root.leaf
  tm
}
