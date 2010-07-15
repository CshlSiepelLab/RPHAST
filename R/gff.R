#' @nord
.makeObj.gff <- function() {
  gff <- list()
  attr(gff, "class") <- "gff"
  gff
}


##' Read a GFF object from a file
##'
##' The function will guess the format of the input file automatically.
##'
##' @title Read a Feature File (GFF, BED, or GenePred)
##' @param filename the name of the file
##' @param pointer.only Whether to store object by reference instead of a
##' data.frame
##' @return If pointer.only==FALSE, a data.frame with columns corresponding
##' to the GFF specification.  Otherwise, an object which is a pointer to
##' an object stored in C.
##' @seealso \code{\link{gff.new}} for more description of GFF objects.
##'
##' \code{\link{msa.new}} for more explanation of the pointer.only option.
##'
##' \url{http://www.sanger.ac.uk/resources/software/gff/spec.html}
##' for a detailed description of GFF file format.
##' 
##' \url{http://genome.ucsc.edu/FAQ/FAQformat} for descriptions
##' of BED and GenePred formats.
##' @keywords GFF
##' @keywords Genepred
##' @keywords BED
##' @export
read.gff <- function(filename, pointer.only=FALSE) {
  gff <- .makeObj.gff()
  gff$externalPtr <- .Call("rph_gff_read", filename)
  if (!pointer.only) {
    gff <- as.data.frame.gff(gff)
  }
  gff
}


##' Create a new GFF object
##'
##' See \url{http://www.sanger.ac.uk/resources/software/gff/spec.html}
##' for more detailed description of each parameter.
##' 
##' All arguments which are provided should be vectors of equal length.
##' 
##' If pointer.only==FALSE, the new GFF object is a data frame, with
##' columns mirroring the GFF Specification
##' Otherwise, it is a list containing a single element, which is
##' a pointer to an object stored in C.
##' @title GFF Objects
##' @param seqname a character vector containing the name of the sequence
##' @param src The source of the feature
##' @param feature The feature type name
##' @param start The start of the feature.  Sequence numbering begins at 1.
##' @param end The end of the feature.  This is the last coordinate included
##' in the feature.
##' @param score The feature score, or NA if there is no score.
##' @param strand A character string which is either "+", "-", or "." (if
##' strand is not available or relevant).
##' @param frame A 0, 1, or 2, which specifies whether the feature is in frame.
##' @param attribute A feature attribute (character string).
##' @param pointer.only Whether to store object as a pointer to an object
##' in C, rather than as a data.frame in R.
##' @return If pointer.only==FALSE, returns a data.frame whose format
##' mirrors the GFF specification.  Otherwise, returns a list with a single
##' object, which is a external pointer to a C structure representing a
##' GFF file.
##' @seealso \code{\link{read.gff}}
##'
##' \code{\link{msa.new}} for more details on the pointer.only option.
##' @keywords GFF feature
##' @export
gff.new <- function(seqname, src, feature, start, end, score=NULL,
                strand=NULL, frame=NULL, attribute=NULL,
                pointer.only=FALSE) {
  check.arg(seqname, "seqname", "character", null.OK=FALSE, min.length=NULL,
             max.length=NULL)
  len <- length(seqname)
  check.arg(src, "src", "character", null.OK=FALSE, min.length=len,
             max.length=len)
  check.arg(feature, "feature", "character", null.OK=FALSE, min.length=len,
             max.length=len)
  check.arg(start, "start", "integer", null.OK=FALSE, min.length=len,
             max.length=len)
  check.arg(end, "end", "integer", null.OK=FALSE, min.length=len,
             max.length=len)
  check.arg(score, "score", "numeric", null.OK=TRUE, min.length=len,
             max.length=len)
  check.arg(strand, "strand", "character", null.OK=TRUE, min.length=len,
             max.length=len)
  check.arg(frame, "frame", "integer", null.OK=TRUE, min.length=len,
             max.length=len)
  check.arg(attribute, "attribute", "character", null.OK=TRUE,
             min.length=len, max.length=len)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  if (pointer.only) {
    if (!is.null(score)) score <- as.numeric(score)
    if (!is.null(strand)) strand <- as.character(strand)
    if (!is.null(frame)) frame <- as.integer(frame)
    if (!is.null(attribute)) attribute <- as.character(attribute)
    ptr <- .Call("rph_gff_new", as.character(seqname),
                 as.character(src), as.character(feature),
                 as.integer(start), as.integer(end),
                 score, strand, frame, attribute)
    gff <- .makeObj.gff()
    gff$externalPtr <- ptr
  } else {
    gff <- data.frame(seqname=seqname, src=src, feature=feature,
                      start=start, end=end)
    if (!is.null(score)) gff <- cbind(gff, score)
    if (!is.null(strand)) gff <- cbind(gff, strand)
    if (!is.null(frame)) gff <- cbind(gff, frame)
    if (!is.null(attribute)) gff <- cbind(gff, attribute)
  }
  gff
}


##' Take a GFF stored in R and return one stored by reference
##' @title GFF To Pointer
##' @param gff a GFF object stored by value in R
##' @return a GFF object stored by reference as a pointer to an
##' object created in C.
##' @seealso \code{\link{gff.new}} for more details on GFF storage
##' options.
##' @export
as.pointer.gff <- function(gff) {
  if (!is.null(gff$externalPtr))
    return(gff)
  gff.new(gff$seqname, gff$src, gff$feature,
          gff$start,   gff$end,    gff$score,
          gff$strand,  gff$frame,  gff$attribute,
          pointer.only=TRUE)
}


##' Prints a GFF object.
##' @title Printing a GFF Object
##' @param x a GFF object
##' @param ... further arguments to be passed to or from other methods
##' @keywords GFF
##' @seealso \code{\link{write.gff}}
##' @export
##' @S3method print gff
print.gff <- function(x, ...) {
  cat(paste("GFF object\n"))
  write.gff(NULL, x)
}


##' Write a GFF object to a file.
##' @title Writing a GFF Object
##' @param filename The name of the file to write to (will be overwritten)
##' @param gff a GFF object
##' @keywords GFF
##' @export
write.gff <- function(filename, gff) {
  check.arg(filename, "filename", "character", null.OK=TRUE)
  if (is.null(gff$externalPtr))
    gff <- as.pointer.gff(gff)
  invisible(.Call("rph_gff_print", filename, gff$externalPtr))
}


##' Get the number of rows in a GFF object
##' @title GFF Number of Rows
##' @param gff A gff object
##' @return An integer containing the number of rows in each gff object
##' @export
##' @S3method nrow gff
nrow.gff <- function(gff) {
  if (is.null(gff$externalPtr))
    return(dim(gff)[1])
  .Call("rph_gff_numrow", gff$externalPtr)
}


##' Get the number of columns in a GFF object
##' @title GFF Number of Columns
##' @param gff A gff object
##' @return An integer containing the number of columns in the gff object
##' @note If the GFF object is stored as a pointer in C, the number
##' of columns is always 9.
##' @export
##' @S3method ncol gff
ncol.gff <- function(gff) {
  if (is.null(gff$externalPtr))
    return(dim(gff)[2])
  9 # gff objects stored in C always have 9 columns
}


##' Prints a brief summary of a GFF object.
##' @title GFF Summary
##' @param object a GFF object
##' @param ... further arguments to be passed to or from other methods
##' @keywords GFF
##' @export
##' @S3method summary gff
summary.gff <- function(object, ...) {
  if (is.null(object$externalPtr)) {
    as <- "stored as data frame"
    object <- as.pointer.gff(object)
  } else as <- "stored as a pointer to a C structure"
  nrow <- nrow.gff(object)
  cat(paste("GFF object with", nrow, "rows", as))
  cat("\n")
}


##' Convert a GFF object to a data frame
##' @title GFF to Data Frame
##' @param x a GFF object
##' @param row.names optional names for each feature
##' @param optional logical, if \code{TRUE}, setting row names and 
##' converting column names (to syntactic names: see
##' \code{\link{make.names}} is optional.
##' @param ... additional arguments to be passed to other methods
##' @return a data frame containing the GFF data
##' @seealso \code{\link{gff.new}} for a description of GFF data frames,
##' and \code{\link{as.pointer.gff}} for conversion in the other
##' direction.
##' @export
##' @S3method as.data.frame gff
as.data.frame.gff <- function(x, row.names=NULL, optional=FALSE, ...) {
  if (is.data.frame(x)) return(x)
  if (!is.null(x$externalPtr)) {
    x <- .Call("rph_gff_dataframe", x$externalPtr)
  }
  attr(x, "class") <- "list"
  as.data.frame(x, row.names, optional, ...)
}


##' Get the dimensions of a GFF feature object
##' @title GFF dimensions
##' @param x a GFF object
##' @return An integer vector of length two containing the number of
##' rows and number of columns in the GFF object.
##' @export
##' @S3method dim gff
dim.gff <- function(x) {
  c(nrow.gff(x), ncol.gff(x))
}
