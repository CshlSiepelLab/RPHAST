#' @nord
gff.makeObj <- function() {
  gff <- list()
  attr(gff, "class") <- "gff"
  gff
}


#' @nord
# don't call explicitly; this is the registered finalizer for GFF
# external pointers
gff.free <- function(extGffPtr) {
  .Call("rph_gff_free", extGffPtr)
}

read.gff <- function(filename, pointer.only=FALSE) {
  obj <- .Call("rph_gff_read", filename)
  if (pointer.only) {
    gff <- gff.makeObj()
    gff$externalPtr <- obj
    reg.finalizer(gff$externalPtr, gff.free)
    return(gff)
  } 
  return(as.data.frame(.Call("rph_gff_dataframe", obj)))
}


#' Create a new GFF object
#'
#' If pointer.only==FALSE, the new GFF object is a data frame, with
#' columns mirroring the GFF Specification \url{http://www.sanger.ac.uk/resources/software/gff/spec.html}.
#' Otherwise, it is a list containing a single element, which is
#' a pointer to an object stored in C.
#'
#' @param seqname a character vector
#' @seealso \code{\link{gff.read}}
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
    reg.finalizer(ptr, gff.free)
    gff <- gff.makeObj()
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


gff.to.pointer <- function(gff) {
  if (!is.null(gff$externalPtr))
    return(gff$externalPtr)
  gff.new(gff$seqname, gff$src, gff$feature,
          gff$start,   gff$end,    gff$score,
          gff$strand,  gff$frame,  gff$attribute,
          pointer.only=TRUE)
}


# todo: print meta-data?
print.gff <- function(gff) {
  cat(paste("GFF object\n"))
  if (is.null(gff$externalPtr))
    gff <- gff.to.pointer(gff)
  invisible(.Call("rph_gff_print", NULL, gff$externalPtr))
}

write.gff <- function(filename, gff) {
  check.arg(filename, "filename", "character", null.OK=TRUE)
  print.gff(filename, gff)
}


summary.gff <- function(gff) {
  if (is.null(gff$externalPtr)) {
    as <- "stored as data frame"
    gff <- gff.to.pointer(gff)
  } else as <- "stored as a pointer to a C structure"
  nrow <- .Call("rph_gff_numrow", gff$externalPtr)
  cat(paste("GFF object with", nrow, "rows", as))
  cat("\n")
}

as.data.frame.gff <- function(gff) {
  attr(gff, "class") <- "data.frame"
  gff
}
