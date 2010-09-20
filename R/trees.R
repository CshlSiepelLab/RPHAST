##' Add a semi-colon to end of tree string
##'
##' Check if tree string ends in semi-colon and if not add one.  This
##' is mostly done for compatibility with ape, which requires them.
##' @param x A character string or vector of character strings each
##' representing a tree in Newick format.
##' @return The same value, but with a semi-colon added to the end
##' of any strings which did not already end in semi-colons.
##' @export
##' @author Melissa J. Hubisz
fix.semicolon.tree <- function(x) {
  n <- nchar(x)
  lastCh <- substr(x, n, n)
  f <- lastCh != ";"
  if (sum(f) > 0L)
    x[f] <- paste(x[f], ";", sep="")
  x
}
    


##' Read a tree from a file
##'
##' Reads a tree in newick format
##' @title Read a Newick Tree from a File
##' @param filename The file containing the tree.
##' @return a character string representing the tree in newick format
##' @keywords trees newick
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
read.newick.tree <- function(filename) {
  check.arg(filename, "filename", "character", null.OK=FALSE)
  tr <- .Call("rph_tree_read", filename)
  fix.semicolon.tree(tr)
}


##' Get the number of nodes in a tree
##' @title Number of Nodes in a Tree
##' @param tree A vector of character strings, each containing a newick tree
##' @return A numeric vector containing the number of nodes in each tree
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
numnodes.tree <- function(tree) {
  check.arg(tree, "tree", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  result <- integer(length(tree))
  for (i in 1:length(tree)) {
    result[i] <- .Call("rph_tree_numnodes", tree[i])
  }
  result
}


##' Get the total length of the edges of a tree
##' @param tree A vector of character strings, each containing a newick tree
##' @return A numeric vector containing the total branchlength of each tree
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
branchlen.tree <- function(tree) {
  check.arg(tree, "tree", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  result <- numeric(length(tree))
  for (i in 1:length(tree))
    result[i] <- .Call("rph_tree_branchlen", tree[i])
  result
}


##' Get the distance from a node to the root of a tree
##' @param tree A vector of character strings, each containing a newick tree
##' @param node A vector of character strings, giving the node name to
##' use for each tree.  Will be recycled to the length of the first argument.
##' @return A numeric vector containing the distance from each given
##' node to the root of the corresponding tree.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
depth.tree <- function(tree, node) {
  check.arg(tree, "tree", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  node <- rep(node, length.out = length(tree))
  result <- numeric(length(tree))
  for (i in 1:length(tree))
    result[i] <- .Call("rph_tree_depth", tree[i], node[i])
  result
}



##' Prune sequences from a file
##' @title Prune a Tree
##' @param tree A vector of character strings, each containing a newick tree
##' @param seqs The sequences to prune from the trees
##' @param all.but A logical value.  If false, prunes all the named sequences
##' from the tree.  If TRUE, prunes all sequences except the ones named.
##' @return a vector of character strings representing the pruned trees.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
prune.tree <- function(tree, seqs, all.but=FALSE) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(seqs, "seqs", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  check.arg(all.but, "all.but", "logical", null.OK=FALSE)
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    result[i] <- .Call("rph_tree_prune", tree[i], seqs, all.but)
  }
  fix.semicolon.tree(result)
}


##' Name ancestors of a tree
##' @title Name Ancestral Nodes
##' @param tree A vector of character strings, each containing a newick tree
##' @return A vector of character strings containing newick trees with all
##' ancestors named.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
name.ancestors <- function(tree) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    result[i] <- .Call("rph_tree_name_ancestors", tree[i])
  }
  fix.semicolon.tree(result)
}


##' Get a subtree
##' @title Subtree
##' @param tree A vector of character strings, each containing a newick tree
##' @param node A vector of character strings, each representing the name
##' of the node which will be the new root of the tree.  If node is shorter
##' than tree, values will be recycled, and a warning produced if \code{length(tree) \%\% length(node) != 0}
##' @param super.tree A vector of logical values.  If TRUE, then remove all
##' nodes which are descendants of node, rather than keeping them.
##' @return A vector of trees which have been pruned, removing all nodes
##' which are not descendants of the given node.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
subtree <- function(tree, node, super.tree=FALSE) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(node, "node", "character", null.OK=FALSE,
            min.length=1, max.length=length(tree))
  check.arg(super.tree, "super.tree", "logical", null.OK=FALSE,
            min.length=1, max.length=length(tree))
  if (length(tree) %% length(node) != 0)
    warning("number of trees is not multiple of number of given nodes")
  node <- rep(node, length.out = length(tree))
  super.tree <- rep(super.tree, length.out = length(tree))
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    if (super.tree[i]) {
      result[i] <- .Call("rph_tree_supertree", tree[i], node[i])
    } else result[i] <- .Call("rph_tree_subtree", tree[i], node[i])
  }
  fix.semicolon.tree(result)
}



##' Rescale a tree
##' @title Scale a Tree or Subtree
##' @param tree A vector of character strings, each containing a newick tree
##' @param scale A vector of scale factors for each tree (will be recycled
##' as necessary if shorter than trees)
##' @param subtree If not NULL, scaling will be on subtree defined by the
##' named node.  Subtrees will be recycled as necessary if shorter than trees.
##' @return A vector of trees whose branches have been scaled
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
rescale.tree <- function(tree, scale, subtree=NULL) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(scale, "scale", "numeric", null.OK=FALSE,
            min.length=1, max.length=length(tree))
  check.arg(subtree, "subtree", "character", null.OK=TRUE,
            min.length=1, max.length=length(tree))
  if (length(tree) %% length(scale) != 0)
    warning("number of trees is not multiple of number of given scales")
  if ((!is.null(subtree)) && length(tree) %% length(subtree) != 0)
    warning("number of trees is not multiple of number of given subtrees")
  if (is.null(subtree)) subtreeVal <-  NULL else subtreeIdx <- 1
  scaleIdx <- 1
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    if (!is.null(subtree)) {
      subtreeVal <- subtree[subtreeIdx]
      subtreeIdx <- subtreeIdx + 1
      if (subtreeIdx > length(subtree)) subtreeIdx <- 1
    }
    result[i] <- .Call("rph_tree_scale", tree[i], scale[scaleIdx], subtreeVal)
    scaleIdx <- scaleIdx+1
    if (scaleIdx > length(scale)) scaleIdx <- 1
  }
  fix.semicolon.tree(result)
}

##' Rename nodes of trees
##' @title Tree Node Renaming
##' @param tree A vector of character strings, each containing a newick tree
##' @param old.names A vector of current names to be substituted
##' @param new.names A vector of equal length to old.names giving the
##' substitutions
##' @return A vector of character strings, in which all nodes with names
##' given in old.names are replaced with values from new.names
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
rename.tree <- function(tree, old.names, new.names) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(old.names, "old.names", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(new.names, "new.names", "character", null.OK=FALSE,
            min.length=length(old.names), max.length=length(old.names))
  fix.semicolon.tree(.Call("rph_tree_rename", tree, old.names, new.names))
}
