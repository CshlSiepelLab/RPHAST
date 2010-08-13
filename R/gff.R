#' @nord
#' @export
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
##' @seealso \code{\link{gff}} for more description of GFF objects.
##'
##' \code{\link{msa}} for more explanation of the pointer.only option.
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
##' \code{\link{msa}} for more details on the pointer.only option.
##' @keywords GFF feature
##' @export
gff <- function(seqname, src, feature, start, end, score=NULL,
                strand=NULL, frame=NULL, attribute=NULL,
                pointer.only=FALSE) {
  
  check.arg(start, "start", "integer", null.OK=FALSE, min.length=NULL,
             max.length=NULL)
  len <- length(start)
  seqname <- rep(seqname, length.out = len)
  check.arg(seqname, "seqname", "character", null.OK=FALSE, min.length=len,
            max.length=len)
  src <- rep(src, length.out = len)
  check.arg(src, "src", "character", null.OK=FALSE, min.length=len,
             max.length=len)
  feature <- rep(feature, length.out = len)
  check.arg(feature, "feature", "character", null.OK=FALSE, min.length=len,
             max.length=len)
  check.arg(end, "end", "integer", null.OK=FALSE, min.length=len,
             max.length=len)
  if (!is.null(score)) score <- rep(score, length.out=len)
  check.arg(score, "score", "numeric", null.OK=TRUE, min.length=len,
             max.length=len)
  if (!is.null(strand)) strand <- rep(strand, length.out=len)
  check.arg(strand, "strand", "character", null.OK=TRUE, min.length=len,
             max.length=len)
  if (!is.null(frame)) frame <- rep(frame, length.out = len)
  check.arg(frame, "frame", "integer", null.OK=TRUE, min.length=len,
             max.length=len)
  if (!is.null(attribute)) attribute <- rep(attribute, length.out = len)
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
##' @seealso \code{\link{gff}} for more details on GFF storage
##' options.
##' @export
as.pointer.gff <- function(gff) {
  if (!is.null(gff$externalPtr))
    return(gff)
  gff(gff$seqname, gff$src, gff$feature,
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
##' @seealso \code{\link{gff}} for a description of GFF data frames,
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


##' Get the range of a GFF feature object
##' @title GFF range
##' @param ... Objects of type \code{gff}
##' @param na.rm Whether to remove values of NA before calculating range.
##' @return A vector of size 2 indicating minimum and maximum coord in
##' the GFF object
##' @S3method range gff
##' @export
range.gff <- function(..., na.rm=FALSE) {
  gffs <- list(...)
  mins <- numeric(length(gffs))
  maxs <- numeric(length(gffs))
  for (i in 1:length(gffs)) {
    gff <- gffs[[i]]
    if (is.data.frame(gff)) {
      r <- range(c(gff$start, gff$end), na.rm=na.rm)
    } else {
      r <- c(.Call("rph_gff_minCoord", gff$externalPtr),
             .Call("rph_gff_maxCoord", gff$externalPtr))
    }
    mins[i] <- r[1]
    maxs[i] <- r[2]
  }
  c(min(mins, na.rm=na.rm),
    max(maxs, na.rm=na.rm))
}


##' plot a GFF
##' @title GFF plot
##' @param x an object of type GFF
##' @param y the location of the plot on the y axis
##' @param width the width of the boxes
##' @param plottype either "r" for rectangles or "a" for arrows
##' @param density the density of the shading lines, in lines per inch.  The
##' default of \code{NULL} implies no lines are drawn.  A zero value of
##' density \code{density} means no shading lines whereas negative values
##' (and \code{NA}) suppress shading (and so allow color filling).
##' @param angle angle (in degrees) of the shading lines.
##' @param col color(s) to fill or shade the boxes with.  The default
##' \code{NA} (or also \code{NULL} means do not fill, i.e., draw
##' transparent rectangles, unless \code{density} is specified.
##' @param border color for rectangle border(s).  The default means
##' \code{par("fg")}.  Use \code{border = NA} to omit borders.  If there
##' are shading lines, \code{border = TRUE} means use the same color for
##' the border as for the shading lines
##' @param lty line type for borders and shading; defaults to \code{"solid"}.
##' @param lwd line width for borders and shading.
##' @param add if \code{TRUE}, add to existing plot
## @param labels whether to label the boxes.  If a vector of strings gives
## the lables for each element.  If \code{TRUE}, use gff$feature for labels.
##' @param xlim A numerical vector of length 2 giving the range for the x-axis.
##' @param ylim A numerical vector of length 2 giving the range for the y-axis.
##' @S3method plot gff
##' @param ... graphical parameters to be passed to \code{plot}.
##' @export
plot.gff <- function(x, y=0, width=1, plottype="r",
                     density=NULL, angle=45,
                     col=NA, border=NULL,
                     lty=par("lty"), lwd=par("lwd"),
                     add=FALSE,
#                     labels=FALSE,
                     xlim=range.gff(x),
                     ylim=c(y-width*3/4, y+width*3/4), ...) {

  if (is.data.frame(x)) {
    starts <- x$start
    ends <- x$end
    if (is.logical(labels) && labels==TRUE) 
      labels <- x$feature
  } else if (!is.null(x$externalPtr)) {
    starts <- .Call("rph_gff_starts", x$externalPtr)
    ends <- .Call("rph_gff_ends", x$externalPtr)
    if (is.logical(labels) && labels==TRUE)
      labels <- .Call("rph_gff_feature", x$externalPtr)
  } else {
    stop("invalid gff object")
  }

  if (!add) {
    coord <- c(0); y<-c(0)
    plot(coord, y, type="n", xlim=xlim, ylim=ylim, ...)
  }
  if (plottype == "a") {  # arrows plot- just do lines if no strand
    stop("gff plottype a is not implemented yet")
    strand <- x$strand
    if (is.null(x$strand)) code <- 0 else code <- ifelse(x$strand=="+", 2, 1)
    arrows(starts, y, ends, y, col=col, code=ifelse(x$strand)) # TODO!
  } else {
    rect(starts, y-width/2, ends, y+width/2,
         density=density, angle=angle, col=col, border=border, lty=lty,
         lwd=lwd)
    # how to do strand in rect plot?  Could use shading lines, might be very slow.
  }
#  if (labels) {  #TODO
    
#  }
  invisible(NULL)
}


##' GFF kernel density
##' @param x An object of type \code{gff}
##' @param type a character string, denoting the value to compute
##' the density for.  Currently the only valid types are "length"
##' and "score"
##' @param ... additional arguments to be passed to \code{density}
##' @return A kernel density object as defined by \link{\code{density}}
##' @export
density.gff <- function(x, type="length", ...) {
  if (type == "length") {
    if (!is.null(x$externalPtr)) {
      vals <- .Call("rphast_gff_lengths", x$externalPtr)
    } else {
      vals <- x$end - x$start
    }
  } else if (type == "score") {
    if (!is.null(x$externalPtr)) {
      vals <- .Call("rphast_gff_getScores", x$externalPtr)
    } else {
      vals <- x$score
    }
  } else stop("unknown type (should be \"length\" or \"score\")")
  density(vals, ...)
}

##' plot histogram of feature lengths
##' @param x an object of type \code{gff}
##' @param type a character string, denoting the value to make the histogram with.
##' Currently the only valid types are "length" or "score"
##' @param ... additional arguments to be passed to \code{hist}
##' @S3method hist gff
##' @export
hist.gff <- function(x, type="length", ...) {
  if (type == "length") {
    if (!is.null(x$externalPtr)) {
      vals <- .Call("rphast_gff_lengths", x$externalPtr)
    } else {
      vals <- x$end - x$start
    }
  } else if (type == "score") {
    if (!is.null(x$externalPtr)) {
      vals <- .Call("rphast_gff_getScores", x$externalPtr)
    } else {
      vals <- x$score
    }
  } else stop("unknown type (should be \"length\" or \"score\")")
  hist(vals, ...)
}

##' Feature overlap
##' @param gff An object of type \code{GFF} containing features to select
##' @param filterGff An object of type \code{GFF} which determines which elements of gff to select
##' @param numbase The number of bases of overlap between gff and filterGff required to choose
##' a record.  Use NULL to ignore (but then min.percent must be defined)
##' @param min.percent The minimum percent that a record must overlap with the combined records in filterGff
##' in order to be chosen
##' @param overlapping If \code{FALSE}, choose records with less than numbase overlapping bases,
##' and less than min.percent fraction overlap if min.percent is not \code{NULL}
##' @param get.fragments If \code{FALSE}, entire records are selected from gff based on whether
##' they meet selection criteria.   
##' If \code{TRUE}, return only the fragments of \code{GFF} that overlap
##' with filterGff.  In this case, the same fragments may be output multiple times, if they are
##' selected by multiple entries in filterGff.  numbase and min.percent apply in either case.
##' @return an object of type GFF containing the selected entries from gff.
##' @export
overlap.gff <- function(gff, filterGff, numbase=1, min.percent=NULL,
                        overlapping=TRUE, get.fragments=FALSE) {
  check.arg(numbase, "numbase", "integer", null.OK=TRUE)
  check.arg(min.percent, "min.percent", "numeric", null.OK=TRUE)
  if (is.null(numbase) && is.null(min.percent)) stop("one of numbase or min.percent should not be NULL")
  if (!is.null(numbase) && numbase < 0) stop("numbase should be at least 1 (if it isn't NULL)")
  if (!is.null(min.percent) && (min.percent < 0 || min.percent > 1))
    stop("min.percent should be NULL or in the range (0,1)")
  check.arg(overlapping, "overlapping", "logical", null.OK=FALSE)
  check.arg(get.fragments, "get.fragments", "logical", null.OK=FALSE)

  if (is.null(gff$externalPtr)) 
    gff <- as.pointer.gff(gff)
  if (is.null(filterGff$externalPtr))
    filterGff <- as.pointer.gff(filterGff)

  rv <- .makeObj.gff()
  rv$externalPtr <- .Call("rph_gff_overlapSelect", gff$externalPtr, filterGff$externalPtr,
                          numbase, min.percent, !overlapping, get.fragments)
  if (!is.null(rv)) {
    rv <- as.data.frame.gff(rv)
  }
  if (nrow.gff(rv) == 0) return(NULL)
  rv
}


##' Take the inverse of a GFF
##' @param gff An object of type \code{gff}
##' @param region.bounds An object of type \code{gff} which defines
##' the boundaries of all relevant chromosomes in the first argument
##' @return An object of type \code{gff} which contains all regions in
##' region.bounds that are not in the first argument
##' @export
inverse.gff <- function(gff, region.bounds) {
  if (is.null(gff$externalPtr))
    gff <- as.pointer.gff(gff)
  if (is.null(region.bounds$externalPtr))
    region.bounds <- as.pointer.gff(region.bounds)
  rv <- .makeObj.gff()
  rv$externalPtr <- .Call("rph_gff_inverse",
                          gff$externalPtr,
                          region.bounds$externalPtr)
  as.data.frame.gff(rv)
}


##' GFF coverage
##' @param ... objects of type \code{gff}
##' @param or if \code{TRUE}, get the coverage of union of gff arguments.
##' or is \code{FALSE} by default, which takes the intersection.
##' @param get.feats if \code{TRUE}, return an object of type \code{gff}
##' representing the intersection (or union of \code{or==TRUE}) of the
##' gff arguments
##' @param not If not \code{NULL}, a vector of logicals the same length
##' as the number of gffs
##' provided.  For each value which is \code{TRUE}, the inverse of the gff
##' will be used.  If any element of \code{not} is \code{TRUE}, then
##' the region.bounds arg must also be provided.  If \code{NULL}, do not
##' take the inverse of any features.
##' @param region.bounds An object of type \code{gff} which defines the
##' start and end coordinates of all relevant chromosomes/regions in the
##' provided gff objects.  Used for taking inverses of the gff objects as
##' required by the argument \code{not}.  If \code{not==NULL} or
##' \code{not==FALSE} for all gff objects, then this argument is not used.
##' @return The number of bases covered by the gff arguments, or the
##' combined gff object if \code{get.feats==TRUE}.
##' @export
coverage.gff <- function(..., or=FALSE, get.feats=FALSE,
                         not=NULL,
                         region.bounds=NULL) {
  check.arg(or, "or", "logical", null.OK=FALSE)
  check.arg(get.feats, "get.feats", "logical", null.OK=FALSE)
  gff <- .makeObj.gff()
  gfflist <- list(...)
  check.arg(not, "not", null.OK=TRUE, min.length=length(gfflist),
            max.length=length(gfflist))
  if (is.null(not)) not <- rep(FALSE, length(gfflist))
  if (sum(not==TRUE) > 0L) {
    if (is.null(region.bounds))
      stop("region.bounds must be provided to take inverse gffs")
  }
  for (i in 1:length(gfflist)) {
    currgff <- gfflist[[i]]
    if (not[i]) currgff <- inverse.gff(currgff, region.bounds)
    if (is.null(currgff$externalPtr)) 
      currgff <- as.pointer.gff(currgff)
    gfflist[[i]] <- currgff$externalPtr
  }
  if (get.feats) {
    rv <- .makeObj.gff()
    rv$externalPtr <- .Call("rph_gff_featureBits", gfflist, or, get.feats)
    return(as.pointer.gff(rv))
  }
  .Call("rph_gff_featureBits", gfflist, or, get.feats)
}


##' Add UTRs to GFF
##' @param gff An object of type \code{gff}.  CDS regions must be present with type "CDS", and
##' the transcript_id must be indicated in the attribute field.
##' @return An object of type \code{gff}, with all the entries of the original GFF, but
##' also with UTR annotations.
##' @export
addUTRs.gff <- function(gff) {
  if (is.null(gff$externalPtr))
    gff <- as.pointer.gff(gff)
  rv <- .makeObj.gff()
  rv$externalPtr <- .Call("rph_gff_add_UTRs", gff$externalPtr)
  as.data.frame.gff(rv)
}

##' Add introns to GFF
##' @param gff An object of type \code{gff}.  CDS regions must be present with type "CDS", and
##' the transcript_id must be indicated in the attribute field.
##' @return An object of type \code{gff}, with all the entries of the original GFF, but
##' also with intron annotations.
##' @export
addIntrons.gff <- function(gff) {
  if (is.null(gff$externalPtr))
    gff <- as.pointer.gff(gff)
  rv <- .makeObj.gff()
  rv$externalPtr <- .Call("rph_gff_add_introns", gff$externalPtr)
  as.data.frame.gff(rv)
}


##' concatenate GFF objects
##' @param ... gff objects to be combined into a single object
##' @return An object of type \code{gff} containing entries from all
##' given gffs
##' @export
rbind.gff <- function(...) {
  gff <- .makeObj.gff()
  gfflist <- list(...)
  for (i in 1:length(gfflist)) {
    currgff <- gfflist[[i]]
    if (is.null(currgff$externalPtr)) 
      currgff <- as.pointer.gff(currgff)
    gfflist[[i]] <- currgff$externalPtr
  }
  gff$externalPtr <- .Call("rph_gff_append", gfflist)
  as.data.frame.gff(gff)
}

##' Split the elements of a GFF
##' @param gff An object of type GFF
##' @param max.length the maximum length of features in new object.  Can
##' be a vector giving a different length for each row, or a single numeric
##' value.  Values will be recycled to the same length
##' as \code{nrow.gff(gff)}.
##' @return An object of type GFF with the same features as gff but
##' with all features of length > max.length broken into segments
##' (starting from the first position in feature).  The last piece
##' of each split segment may be smaller than max.length
##' @export
split.gff <- function(gff, max.length) {
  gffSize <- nrow.gff(gff)
  check.arg(max.length, "max.length", "integer", null.OK=FALSE,
            min.length = 1, max.length = NULL)
  if (is.null(gff$externalPtr))
    gff <- as.pointer.gff(gff)
  splitGff <- .makeObj.gff()
  splitGff$externalPtr <- .Call("rph_gff_split", gff$externalPtr,
                                max.length)
  as.data.frame.gff(splitGff)
}


##' Sort a GFF
##' @param x An object of type \code{gff}
##' @param decreasing Set to TRUE to sort from highest to lowest coordinates
##' @param ... Currently not used
##' @return An object of type \code{gff} sorted primarily by
##' start position, then by end position.
##' @export
##' @S3method sort gff
sort.gff <- function(x, decreasing = FALSE, ...) {
  gff <- x
  if (is.null(gff$externalPtr))
    gff <- as.pointer.gff(gff)
  rv <- .makeObj.gff()
  rv$externalPtr <- .Call("rph_gff_sort", gff$externalPtr)
  rv <- as.data.frame.gff(rv)
  if (decreasing) 
    rv <- rv[dim(rv)[1]:1,]
  rv
}
  

##' Composition of features with respect to annotations
##' @param features An object of type \code{gff}.
##' @param annotations An object of type \code{gff} containing
##' some annotations.
##' @return A data frame with two columns and a row for
##' each type of element in the annotations.  The second
##' column gives the fraction of 
##' features which fall in the corresponding annotation type.
##' Assuming non-overlapping annotations which cover the
##' entire range of interest, the second column should sum
##' to 1.
##' @export
composition.gff <- function(features, annotations) {
  if (!is.null(annotations$externalPtr))
    annotations <- as.data.frame.gff(annotations)
  if (!is.null(features$externalPtr))
    features <- as.data.frame.gff(features)
  annTypes <- unique(annotations$feature)
  rv <- list()
  for (anntype in annTypes) {
    annfeat <- annotations[annotations$feature == anntype,]
    rv[[anntype]] <- coverage.gff(features, annfeat)/coverage.gff(features)
  }
  data.frame(type=names(rv), composition=as.numeric(rv), stringsAsFactors=TRUE)
}


##' Enrichment of features with respect to annotation types
##' @param features An object of type \code{gff}
##' @param annotations An object of type \code{gff} containing
##' some annotations.
##' @param region.bounds An object of type \code{gff} representing
##' the boundary coordinates of the regions of interest (such as chromosome
##' boundaries).  All elements from the first two arguments should fall
##' entirely within region.bounds.
##' @return A data frame with two columns and a row for each type of
##' element in the annotations.  The second column gives the 
##' fold-enrichment
##' of the features across the corresponding annotation type,
##' which is equal to the
##' fraction of the features which fall within the annotation type,
##' divided by the fraction of the entire region covered by the
##' annotation type.
##' @export
enrichment.gff <- function(features, annotations, region.bounds) {
  if (!is.null(annotations$externalPtr))
    annotations <- as.data.frame.gff(annotations)
  if (!is.null(features$externalPtr))
    features <- as.data.frame.gff(feature)
  annTypes <- unique(annotations$feature)
  rv <- list()
  totalNumBase <- coverage.gff(region.bounds)
  for (anntype in annTypes) {
    annfeat <- annotations[annotations$feature == anntype,]
    cat(dim(annfeat),"\n")
    rv[[anntype]] <- coverage.gff(features, annfeat)*totalNumBase/(coverage.gff(annfeat)*coverage.gff(features))
  }
  data.frame(type=names(rv), enrichment=as.numeric(rv), stringsAsFactors=TRUE)
}
