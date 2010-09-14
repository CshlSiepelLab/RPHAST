# make barebones msa obj
#' @nord
#' @export
.makeObj.msa <- function() {
  msa <- list()
  class(msa) <- "msa"
  msa
}

##' Creates a copy of an MSA sequence
##'
##' If m is stored in R (as it is by default), then m2 <- copy.msa(m1)
##' is no different than m2 <- m1.  But if it is stored as a pointer
##' to a C structure, this is the only way to make an explicit copy
##' of the MSA object.
##' @title MSA copy
##' @param msa an MSA object
##' @return an MSA object which can be modified independently from the
##' original object
##' @export
copy.msa <- function(msa) {
  if (is.null(msa$externalPtr)) return(msa)
  result <- .makeObj.msa()
  result$externalPtr <- .Call("rph_msa_copy", msa$externalPtr)
  result
}


##' Check an MSA object
##' @param msa An object to tests
##' @return A logical indicating whether object is of type \code{msa}
##' @export
##' @keywords msa
is.msa <- function(msa) {
  attr(msa, "class")=="msa"
}


##' Creates a new MSA object given sequences.
##'
##' Make a new multiple sequence alignment (MSA) object given a vector of
##' character strings.  They can be optionally annotated with sample names.
##'
##' Each character string in seqs must be the same length, and number of
##' elements in names (if provided) must match the number of elements in
##' seqs.
##'
##' Alphabet generally does not have to be specified if working with
##' DNA alignments.
##'
##' About storing objects as pointers:
##' If \code{pointer.only==FALSE}, the MSA object will be stored in R and can be
##' viewed and modified by base R code as well as RPHAST functions.
##' Setting \code{pointer.only=TRUE} will cause the object to be stored by
##' reference, as an external pointer to an object created by C code.  This
##' may be necessary to improve performance, but the object can then only
##' be viewed/manipulated via RPHAST functions.  Furthermore, if an object
##' is stored as a pointer, then its value is liable to be changed when
##' passed as an argument to a function.  All RPHAST functions which change
##' the value of an external pointer make a note of this in the help pages
##' for that function.  For example, extract.feature.msa will alter an
##' alignment if it is passed in as an external pointer (the argument will
##' be changed into the return value).  If this is undesireable, the copy.msa
##' function can be used: extract.feature.msa(copy.msa(align)) will preserve
##' the original alignment.  Simple copying, ie, \code{align2->align1} of
##' objects stored as pointer will not behave like normal R objects: both objects
##' will point to the same C structure, and both will be changed if either one
##' is altered.  Instead \code{align2 <- copy.msa(align1)} should be used.
##' @title MSA Objects
##' @param seqs a character vector containing sequences, one per sample
##' @param names a character vector identifying the sample name for each
##' sequence
##' @param alphabet a character string containing valid non-missing character
##' states
##' @param is.ordered a logical indicating whether the alignment columns
##' are stored in order.  If NULL, assume columns are ordered.
##' @param offset an integer giving the offset of coordinates for the
##' reference sequence from the beginning of the chromosome.  The reference
##' sequence is assumed to be the first sequence.  Not used
##' if is.ordered==FALSE.
##' @param pointer.only a boolean indicating whether returned alignment object
##' should be stored by reference (see Details)
##' @useDynLib rphast
##' @keywords msa
##' @export msa
msa <- function(seqs, names = NULL, alphabet="ACGT", is.ordered=TRUE,
                offset=NULL, pointer.only=FALSE) {
  #checks
  seqlen <-  unique(sapply(seqs, nchar))
  if (length(seqlen) > 1L) {
    stop("sequences should all have same length")
  }
  # check number of names
  if (!is.null(names) && length(names) != length(seqs)) {
    stop("number of names needs to match number of sequences")
  }
  check.arg(alphabet, "alphabet", "character")
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  check.arg(is.ordered, "is.ordered", "logical", null.OK=TRUE)
  check.arg(offset, "offset", "integer", null.OK=TRUE)
  if (is.null(is.ordered)) is.ordered <- TRUE
  if ( (!is.ordered) && (!is.null(offset)) && offset!=0) 
    offset <- NULL

  # TODO: if alphabet non-null, check that seqs only contains those chars?

  msa <- .makeObj.msa()

  if (pointer.only) {
    msa$externalPtr <- .Call("rph_msa_new",
                             seqsP=seqs,
                             namesP=names,
                             nseqsP=length(seqs),
                             lengthP=seqlen,
                             alphabetP=alphabet,
                             orderedP=is.ordered,
                             offsetP=offset)
  } else {
    msa$seqs <- seqs
    if (! is.null(names)) msa$names <- names
    if (! is.null(alphabet)) msa$alphabet <- alphabet
    msa$is.ordered <- is.ordered
    if (!is.null(offset)) msa$offset <- offset
  }
  msa
}

##' Returns the length of sequence in an MSA alignment.
##' @title MSA Sequence Length.
##' @param x an MSA object
##' @param refseq character vector giving name(s) of sequence whose
##' length to return.  The default \code{NULL} implies the frame of
##' reference of the entire alignment.
##' @return an integer vector containing the length of the named sequences.
##' If refseq is NULL, returns the number of columns in the alignment.
##' @keywords msa
##' @seealso \code{\link{msa}}
##' @export
ncol.msa <- function(x, refseq=NULL) {
  msa <- x
  if (is.null(msa$externalPtr) && is.null(refseq)) {
    return(nchar(msa$seqs[1]))
  }
  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  if (is.null(refseq))
    return(.Call("rph_msa_seqlen", msa$externalPtr, NULL))
  result <- integer(length(refseq))
  for (i in 1:length(refseq)) 
    result[i] <- .Call("rph_msa_seqlen", msa$externalPtr, refseq[i])
  result
}


##' Returns the dimensions of an msa object as (# of species, # of columns)
##' @param x An object of type \code{msa}
##' @return An integer vector of length two giving number of species and
##' number of columns in the alignment
##' @S3method dim msa
##' @export
##' @keywords msa
dim.msa <- function(x) {
  c(nrow.msa(x), ncol.msa(x, NULL))
}


##' The number of informative columns in an alignment
##' @param msa An object of type \code{msa}
##' @return The number of "informative" columns in the msa.  An informative
##' column has at least two non-missing and non-gap characters.
##' @export
##' @keywords msa
ninf.msa <- function(msa) {
  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  .Call("rph_msa_ninformative_sites", msa$externalPtr)
}

##' Returns the number of sequence in an MSA alignment.
##' @title MSA Number of Sequences
##' @param msa an MSA object
##' @return an integer containing the number of sequences in an alignment.
##' @keywords msa
##' @seealso \code{\link{msa}}
##' @export
nrow.msa <- function(msa) {
  if (is.null(msa$externalPtr)) {
    return(length(msa$seqs))
  }
  .Call("rph_msa_nseq", msa$externalPtr)
}

##' Returns TRUE if the argument is a valid string describing a
##' multiple sequence alignment (MSA) format.
##'
##' Valid formats include "FASTA", "PHYLIP", "SS" (Sufficient statistics
##' format used by PHAST), "MPM" (format used by MultiPipMaker),
##' "LAV" (used by blastz), or "MAF" (Multiple Alignment Format used by
##' MULTIZ and TBA.
##' @title Check an MSA Format String
##' @param format a character vector of strings to test
##' @return a logical vector indicating whether each element of the
##' input parameter is a valid format string.
##' @keywords msa
##' @export
validFormatStr.msa <- function(format) {
  if (is.null(format)) return(NULL)
  result <- logical(length(format))
  for (i in 1:length(format)) {
    result[i] <- .Call("rph_msa_valid_fmt_str", format[i]);
  }
  result
}

##' Returns the offset of the first position in an alignment from
##' some reference sequence.
##' @title MSA Index Offset
##' @param msa an MSA object
##' @return The difference between the first position in an alignment
##' from the beginning of a chromosome.
##' @keywords msa
##' @export
offset.msa <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_idxOffset", msaP=msa$externalPtr))
  if (is.null(msa$offset)) return (0)
  msa$offset
}

##' Returns the alphabet used by an MSA object.
##' @title MSA Alphabet
##' @param msa an MSA object
##' @return the valid non-missing-data characters for an MSA object.
##' @keywords msa
##' @export
alphabet.msa <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_alphabet", msaP=msa$externalPtr))
  msa$alphabet
}

##' Determines if an MSA object represents an ordered alignment.
##' @title MSA is Ordered?
##' @param msa an MSA object
##' @return a boolean indicating whether the columns are in order
##' @keywords msa
##' @export
is.ordered.msa <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_isOrdered", msaP=msa$externalPtr))

  if (is.null(msa$is.ordered)) return (TRUE)
  msa$is.ordered
}


##' Returns the sequence names for an MSA object.
##' @title MSA Sequence Names
##' @param x an MSA object
##' @return a character vector giving the names of the sequences, or
##' NULL if they are not defined
##' @keywords msa
##' @export
##' @S3method names msa
names.msa <- function(x) {
  if (!is.null(x$externalPtr))
    return(.Call("rph_msa_seqNames", msaP=x$externalPtr))
  x$names
}


##' Take an MSA stored by reference and return one stored in R
##' @title MSA From Pointer
##' @param src an MSA object stored by reference
##' @return an MSA object stored in R.  If src is already stored in R,
##' returns a copy of the object.
##' @seealso \code{\link{msa}} for details on MSA storage options.
##' @keywords msa
##' @export
from.pointer.msa <- function(src) {
  if (is.null(src$externalPtr)) return(src)
  seqs <- .Call("rph_msa_seqs", src$externalPtr)
  
  names <- .Call("rph_msa_seqNames", src$externalPtr)
  alphabet <- .Call("rph_msa_alphabet", src$externalPtr)
  ordered <- .Call("rph_msa_isOrdered", src$externalPtr)
  offset <- .Call("rph_msa_idxOffset", src$externalPtr)
  msa(seqs, names, alphabet, is.ordered=ordered,
          offset=offset, pointer.only=FALSE)
}


##' Take an MSA stored in R and return one stored by reference
##' @title MSA To Pointer
##' @param src an MSA object stored by value in R
##' @return an MSA object stored by reference as a pointer to an object
##' created in C.
##' @seealso \code{\link{msa}} for details on MSA storage options.
##' @keywords msa
##' @export
as.pointer.msa <- function(src) {
  if (!is.null(src$externalPtr)) return(src)
  msa(seqs=src$seq,
      names=src$names,
      alphabet=src$alphabet,
      is.ordered=src$is.ordered,
      offset=src$offset,
      pointer.only=TRUE)
}


##' Guess the format of an MSA file by looking at the filename extension.
##' @title MSA Format From Filename Extension
##' @param filename The name of an MSA file
##' @return A string describing an MSA file (one of "MAF", "FASTA",
##' "LAV", "SS", "MPM", "PHYLIP"), or NULL if a guess can't be made
##' based on the file extension.
##' @seealso validFormatStr.msa
##' @keywords msa
##' @export
guess.format.msa  <- function(filename) {
  if (is.null(filename)) return(NULL)
  x <- strsplit(filename, ".", fixed=TRUE)[[1]]
  x <- x[length(x)]
  if (x=="MAF" || x=="maf") return("MAF")
  if (x=="FA"  || x=="fa") return("FASTA")
  if (x=="LAV" || x=="lav") return("LAV")
  if (x=="SS"  || x=="ss") return("SS")
  NULL
}


##' Writes a multiple sequence alignment (MSA) object to a file
##' in one of several formats.
##' @title Writing MSA Objects to Files
##' @param msa an object of class msa
##' @param filename File to write (will be overwritten).  If NULL, output
##' goes to terminal.
##' @param format format to write MSA object.  Valid values are "FASTA",
##' "PHYLIP", "MPM", or "SS".
##' @param pretty.print Whether to pretty-print alignment (turning
##' bases which match the first base in the same column to ".").
##' @note pretty.print does not work if format="SS".
##' @keywords msa FASTA PHYLIP MPM SS
##' @export
write.msa <- function(msa, filename=NULL,
                      format=c(guess.format.msa(filename), "FASTA")[1],
                      pretty.print=FALSE) {
  #checks
  check.arg(filename, "filename", "character", null.OK=TRUE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(pretty.print, "pretty.print", "logical", null.OK=FALSE)
  if (! validFormatStr.msa(format)) {
    stop(paste("invalid MSA FORMAT \"", format, "\"", sep=""))
  }
  if (is.null(msa$externalPtr)) {
    printMsa <- as.pointer.msa(msa)
  } else {
    printMsa <- msa
  }
  invisible(.Call("rph_msa_printSeq",
                  msaP=printMsa$externalPtr,
                  filenameP=filename,
                  formatP=format,
                  pretty.printP=pretty.print))
}




##' Prints a short description of an MSA (multiple sequence alignment)
##' object.
##'
##' @title MSA Summary
##' @param object an MSA object
##' @param ... additional arguments sent to \code{print}
##' @param print.seq whether to supress printing of the alignment
##' @param format to print sequence in if printing alignment
##' @param pretty.print whether to pretty.print pretty-print sequence if printing alignment
##' @keywords msa
##' @seealso \code{\link{print.msa}}
##' @export
##' @S3method summary msa
summary.msa <- function(object, ...,
                        print.seq=ifelse(ncol.msa(object)*nrow.msa(object) < 500, TRUE, FALSE),
                        format="FASTA",
                        pretty.print=FALSE) {
  msa <- object
  check.arg(print.seq, "print.seq", "logical", null.OK=FALSE)
  # format and pretty.print are checked in write.msa

  cat(paste("msa object with", nrow.msa(msa), "sequences and",
            ncol.msa(msa),"columns, stored"))
  if (is.null(msa$externalPtr)) cat(" in R") else cat(" as a pointer to a C structure")
  cat("\n")

  if (!is.null(format) || pretty.print) {
    print.seq=TRUE
  }
  
  pointer <- msa$externalPtr
  names <- names.msa(msa)
  alphabet <- alphabet.msa(msa)
  is.ordered <- is.ordered.msa(msa)
  offset <- offset.msa(msa)

  printMsa <- list()
  printed <- FALSE
  if (!is.null(names)) printMsa$names <- names
  if (!is.null(alphabet)) printMsa$alphabet <- alphabet
  if (!is.null(is.ordered)) printMsa$is.ordered <- is.ordered
  if (!is.null(offset) && offset!=0) printMsa$offset <- offset
  if (print.seq && is.null(pointer) && is.null(format) && !pretty.print) {
    printMsa$seq <- msa$seq
    printed <- TRUE
  }
  
  print(printMsa, ...)

  if (!printed && print.seq) {
    cat("$seq\n")
    if (is.null(format)) format <- "FASTA"
      write.msa(msa, filename=NULL, format, pretty.print)
    cat("\n")
  }
}


##' Prints an MSA (multiple sequence alignment) object.
##'
##' Valid formats for printing are "FASTA", "PHYLIP", "MPM", and "SS".
##' See \code{\link{validFormatStr.msa}} for details on these formats.
##' If format is specified, the alignment is printed regardless of
##' print.seq.
##'
##' Pretty-printing will cause all characters in a column which match
##' the value in the first row to be printed as ".".  It only works for
##' FASTA, PHYLIP, or MPM formats.
##'
##' If print.seq==TRUE, then the default printing format depends on whether
##' the sequence is stored by value (the default storage mode), or by 
##' reference.  If the MSA is stored by value, the default format is
##' as a R character vector.  Otherwise, the default format is FASTA.
##'
##' @title Printing MSA objects
##' @param x an object of class msa
##' @param ... additional arguments sent to \code{print}
##' @param print.seq whether to supress printing of the alignment
##' @param format to print sequence in if printing alignment
##' @param pretty.print whether to pretty.print pretty-print sequence if printing alignment
##' @keywords msa
##' @export
##' @S3method print msa
print.msa <- function(x, ..., print.seq=ifelse(ncol.msa(x)*nrow.msa(x) < 500, TRUE, FALSE),
                      format=NULL, pretty.print=FALSE) {
  check.arg(print.seq, "print.seq", "logical", null.OK=FALSE)
  # format and pretty.print are checked in write.msa

  summary.msa(x, print.seq=print.seq, format=format, pretty.print=pretty.print, ...)

  if (!print.seq && is.null(format)) {
    cat("(alignment output suppressed)\n")
    cat("\n")
  }
}



##' Reads an MSA from a file.
##' @title Reading an MSA Object
##' @param filename The name of the input file containing an alignment.
##' @param format input file format: one of "FASTA", "MAF", "SS", "PHYLIP",
##' "MPM", must be correctly specified.
##' @param alphabet the alphabet of non-missing-data chraracters in the
##' alignment.  Determined automatically from the alignment if not given.
##' @param features An object of type \code{feat}.  If provided, the return
##' value will only
##' contain portions of the alignment which fall within a feature.
##' The alignment will not be ordered.
##' The loaded regions can be further constrained with the do.4d or
##' do.cats options.  Note that if this object is passed as a pointer to a
##' structure stored in C, the values will be altered by this function!
##' @param do.4d Logical.  If \code{TRUE}, the return value will contain only
##' the columns corresponding to four-fold degenerate sties.  Requires
##' features to be specified.
##' @param ordered Logical.  If \code{FALSE}, the MSA object may not retain
##' the original column order.
##' @param tuple.size Integer.  If given, and if pointer.only is \code{TRUE},
##' MSA will be stored in sufficient statistics format, where each tuple
##' contains tuple.size consecutive columns of the alignment.
##' @param do.cats Character vector.  If given, and if features is specified,
##' then only the types of features named here will be represented in the
##' returned alignment.
##' @param refseq Character string specifying a FASTA format file with a
##' reference sequence.  If given, the reference sequence will be
##' "filled in" whereever missing from the alignment.
##' @param offset An integer giving offset of reference sequence from
##' beginning of chromosome.  Not used for MAF format.
##' @param seqnames A character vector.  If provided, discard any sequence
##' in the msa that is not named here.  This is only implemented efficiently
##' for MAF input files, but in this case, the reference sequence must be
##' named.
##' @param discard.seqnames A character vector.  If provided, discard
##' sequenced named here.  This is only implemented efficiently for MAF
##' input files, but in this case, the reference sequenced must NOT be
##' discarded.
##' @param pointer.only If \code{TRUE}, MSA will be stored by reference as
##' an external pointer to an object created by C code, rather than
##' directly in R memory.  This improves performance and may be necessary
##' for large alignments, but reduces functionality.  See
##' \code{\link{msa}} for more details on MSA object storage options.
##' @note If the input is in "MAF" format and features is specified, the
##' resulting alignment will be stripped of gaps in the reference (1st)
##' sequence.
##' @note If the features argument is an object stored in C, its values will
##' be changed by this function!
##' @return an MSA object.
##' @keywords msa FASTA MAF PHYLIP SS
##' @seealso \code{\link{msa}}, \code{\link{read.feat}}
##' @export
read.msa <- function(filename,
                     format=c(guess.format.msa(filename), "FASTA")[1],
                     alphabet=NULL,                     
                     features=NULL,
                     do.4d=FALSE,
                     ordered=ifelse(do.4d || !is.null(features), FALSE, TRUE),
                     tuple.size=(if(do.4d) 3 else NULL),
                     do.cats=NULL,
                     refseq=NULL,
                     offset=0,
                     seqnames=NULL, discard.seqnames=NULL,
                     pointer.only=FALSE) {

  check.arg(filename, "filename", "character", null.OK=FALSE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(alphabet, "alphabet", "character", null.OK=TRUE)
  check.arg(features, "feat", null.OK=TRUE, min.length=NULL, max.length=NULL)
  check.arg(do.4d, "do.4d", "logical", null.OK=FALSE)
  check.arg(ordered, "ordered", "logical", null.OK=FALSE)
  check.arg(tuple.size, "tuple.size", "integer", null.OK=TRUE)
  check.arg(do.cats, "do.cats", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  check.arg(offset, "offset", "integer", null.OK=TRUE)
  check.arg(seqnames, "seqnames", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(discard.seqnames, "discard.seqnames", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  if (!validFormatStr.msa(format))
    stop(paste("invalid format string", format))
  if (!is.null(tuple.size) && tuple.size <= 0)
    stop("tuple.size should be integer >= 1")
  
  if (do.4d) {
    if (!is.null(do.cats))
      stop("should not specify do.cats if do.4d==TRUE")
    if (is.null(features))
      stop("features needs to be specified with do.4d")
    if (tuple.size != 3)
      stop("tuple.size must be 3 if do.4d==TRUE")
  }

  if (!is.null(do.cats) && is.null(features))
    stop("features required with do.cats")

  if (!is.null(features)) {
    if (is.null(features$externalPtr))
      features <- as.pointer.feat(features)
  }

  msa <- .makeObj.msa()
  msa$externalPtr <- .Call("rph_msa_read", filename, format,
                           features$externalPtr, do.4d, alphabet,
                           tuple.size, refseq, ordered, do.cats,
                           offset, seqnames, discard.seqnames)
  if (format != "MAF") {
    if (!is.null(seqnames))
      msa <- sub.msa(msa, seqnames, keep=TRUE)
    if (!is.null(discard.seqnames))
      msa <- sub.msa(msa, discard.seqnames, keep=FALSE)
  }
  if (pointer.only == FALSE) 
    msa <- from.pointer.msa(msa)
  msa
}


##' DNA complement
##'
##' @title complement
##' @param x A character vector with DNA sequences to be complemented
##' @return The complement of the given sequence(s).  Characters other
##' than A,C,G,T,a,c,g,t are unchanged.
##' @keywords msa
##' @export
complement <- function(x) {
  chartr("ACGTacgt", "TGCAtgca", x)
}


##' Reverse complement a multiple sequence alignment
##' @param x An object of type \code{msa}.
##' @return The reverse complement of msa.
##' @note If x is stored as a pointer to an object in C, x will be changed to
##' its reverse complement.  Use reverse.complement(copy.msa(x)) to avoid this
##' behavior.  The return value will be a pointer if the input value was stored
##' as a pointer.
##' @keywords msa
##' @export
reverse.complement <- function(x) {
  if (is.null(x$externalPtr)) {
    pointer.only <- FALSE
    x <- as.pointer.msa(x)
  } else pointer.only <- TRUE
  rv <- .makeObj.msa()
  rv$externalPtr <- .Call("rph_msa_reverse_complement", x$externalPtr)
  if (!pointer.only)
    rv <- from.pointer.msa(rv)
  rv
}

       
##' Get a subset of an alignment
##'
##' @title MSA Subset
##' @param x An object of type \code{msa}
##' @param seqs The sequence names to keep (or to remove if keep is
##' \code{FALSE})
##' @param keep Whether to keep the named sequences or remove them
##' @param start.col the first column to keep (columns indices start at 1)
##' @param end.col the last column to keep (inclusive)
##' @param refseq the sequence in the alignment which determines the
##' coordinates for start.col and end.col.  If NULL, start.col and
##' end.col are column indices in the multiple alignment.
##' @param pointer.only If \code{TRUE}, return an msa object which is only
##' a pointer to a C structure (advanced use only).
##' @return A new MSA object containing a subset of the original MSA.
##' @S3method sub msa
##' @note This function will not modify x even if it is stored as a pointer.
##' @export
##' @keywords msa
sub.msa <- function(x, seqs=NULL, keep=TRUE, start.col=NULL, end.col=NULL,
                    refseq=NULL, pointer.only=FALSE) {
  check.arg(keep, "keep", "logical", null.OK=FALSE)
  check.arg(seqs, "seqs", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(start.col, "start.col", "integer", null.OK=TRUE)
  check.arg(end.col, "end.col", "integer", null.OK=TRUE)
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  
  result <- .makeObj.msa()
  if (is.null(x$externalPtr)) {
    x <- as.pointer.msa(x)
  }
  result$externalPtr <- .Call("rph_msa_sub_alignment",
                              x$externalPtr, seqs, keep,
                              start.col, end.col, refseq)

  if (!pointer.only)
    result <- from.pointer.msa(result)
  result
}


##' Strip gaps from an alignment.
##'
##' If strip.mode can be a vector of integers or a vector of character
##' strings.  If it is a vector of integers, these are the indices of
##' the sequences from which to strip gaps.
##' If strip.mode is vector of character strings, each string names a
##' sequence from which to strip gaps.
##'
##' strip.mode can also be the string "all.gaps" or "any.gaps".  The former
##' will strip columns containing only gaps, whereas the latter strips
##' columns containing even a single gap.
##'
##' @title MSA Strip Gaps
##' @param x MSA object
##' @param strip.mode Determines which gaps to strip.  See Details
##' @return an MSA object, with gaps stripped according to strip.mode.
##' @note If x is passed as a pointer to a C structure (ie,
##' it was created with pointer.only=TRUE), then this function will directly
##' modify x.  Use strip.gaps.msa(copy.msa(x)) to avoid this behavior.  Also,
##' the return value will be stored as a pointer if x is stored as a pointer;
##' otherwise the return value will be stored in R.
##' @keywords msa
##' @export
strip.gaps.msa <- function(x, strip.mode=1) {
  names <- NULL
  nseq <- NULL
  if (is.null(x$externalPtr)) {
    names <- names.msa(x)
    nseq <- nrow.msa(x)
    x <- as.pointer.msa(x)
    copy.to.R <- TRUE
  } else copy.to.R <- FALSE
  for (s in strip.mode) {
    if (s=="all.gaps" || s=="any.gaps")
      x$externalPtr <- .Call("rph_msa_strip_gaps", x$externalPtr, 0, s)
    else {
      if (!is.character(s)) {
        if (is.null(nseq)) nseq <- nrow.msa(x)
        if (as.integer(s) != s || s <=0 || s>nseq)
          stop(cat("invalid sequence index", s))
        w <- s
      } else {
        if (is.null(names))
          names <- names.msa(x)
        w <- which(names==s)
        if (is.null(w))
          stop(cat("no sequence with name", s))
      }
      x$externalPtr <- .Call("rph_msa_strip_gaps", x$externalPtr, w, NULL)
    }
  }
  if (copy.to.R) 
    x <- from.pointer.msa(x)
  x
}



##' Extract, replace, reorder MSA
##'
##' Treat multiple sequence alignment as a matrix where each row
##' corresponds to a sequence for one species, and each column
##' is one position aligned across all species.
##'
##' The bracket notation can return a subset of the alignment,
##' or re-order rows and columns.
##' @param x An object of type \code{msa}
##' @param rows A numeric vector of sequence indices,
##' character vector (containing sequence name), or
##' logical vector (containing sequences to keep).  If logical vector it
##' will be recycled as necessary to the same length as \code{nrow.msa(x)}.
##' @param cols A numeric vector of alignment columns, or a logical vector
##' containing columns to keep.  If logical vector it will be recycled as
##' necessary to the same length as \code{ncol.msa(x)}.  Note that these are
##' in coordinates with respect to the entire alignment.  x$idx.offset
##' is ignored here.
##' @param pointer.only If \code{TRUE}, return an object which is only
##' a pointer to a structure stored in C (useful for large alignments;
##' advanced use only).  In certain cases when the original alignment
##' is stored in R, it may be more efficient return an object in R, in which
##' case this argument will be ignored.
##' @seealso \code{\link{sub.msa}} which can subset columns based on genomic
##' coordinates, and \code{\link{extract.feature.msa}} which can subset based
##' on genomic coordinates denoted in a features object.
##' @usage \method{[}{msa}(x, rows, cols, pointer.only)
##' @S3method "[" msa
##' @note This function will not alter the value of x even if it is stored as
##' a pointer to a C structure.
##' @keywords msa
##' @export "[.msa"
"[.msa" <- function(x, rows, cols, pointer.only=FALSE) {
  if (!missing(rows)) {
    if (is.null(rows)) stop("rows cannot be empty")
  } else rows=NULL
  if (!missing(cols)) {
    if (is.null(cols)) stop("cols cannot be empty")
  } else cols=NULL
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  if (!is.null(rows)) {
    # if rows are given by names, convert to integer
    if (is.character(rows)) {# names are given
      names <- names.msa(x)
      rows <- as.numeric(sapply(rows, function(x) {which(x ==  names)}))
      if (sum(is.na(rows)) > 0L)
        stop("unknown names in first dimension subset")
    }
  }

  # check if arguments are given as logicals.
  if (is.logical(rows)) {
    cat("is.logical rows")
    rows = which(rep(rows, length.out = nrow.msa(x)))
  }
  if (is.logical(cols)) 
    cols = which(rep(cols, length.out = ncol.msa(x)))

  # if x is stored in R, sampling rows is easier and more efficient to do here
  if (!is.null(rows) && is.null(x$externalPtr)) {
    x$names <- x$names[rows]
    x$seqs <- x$seqs[rows]
    rows <- NULL
  }
  if (is.null(rows) && is.null(cols)) return(x)
  check.arg(rows, "rows", "integer", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(cols, "cols", "integer", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)
  rv <- .makeObj.msa()
  rv$externalPtr <- .Call("rph_msa_square_brackets", x$externalPtr,
                          rows, cols)
  if (!pointer.only)
    rv <- from.pointer.msa(rv)
  rv
}




##' Likelihood of an alignment given a tree model
##' @title MSA Likelihood
##' @param x An object of class \code{msa} representing the multiple alignment
##' @param tm An object of class \code{tm} representing the tree and model of
##' substitution
##' @param features A features object.  If non-null, compute likelihoods
##' for each feature rather than the whole alignment.
##' @param by.column Logical indicating whether to get likelihoods for
##' each alignment column.  If FALSE, returns total likelihood.  Ignored
##' if features is not NULL.
##' @return Either the likelihood of the entire alignment (if
##' \code{by.column==FALSE && is.null(features)},
##' or a numeric vector giving the likelihood of each feature
##' (if \code{!is.null(features)}), or a numeric vector giving the likelihood
##' of each column (if \code{by.column==TRUE}).
##' @seealso \code{phyloFit}, \code{tm}
##' @keywords msa tm features
##' @export
likelihood.msa <- function(x, tm, features=NULL, by.column=FALSE) {
  if (is.null(features))
    check.arg(by.column, "by.column", "logical", null.OK=FALSE)
  else {
    if (is.null(features$externalPtr))
      features <- as.pointer.feat(features)
    else features <- copy.feat(features)  # if we don't make a copy features gets
                                          # destroyed by re-mapping coordinates to msa
    if (by.column) {
      warning("by.column ignored when features is NULL")
      by.column <- FALSE
    }
  }
  if (is.null(x$externalPtr)) 
    x <- as.pointer.msa(x)
  tm <- as.pointer.tm(tm)
  if (by.column && !is.ordered.msa(x))
    warning("by.column may not be a sensible option for unordered MSA")
  .Call("rph_msa_likelihood", x$externalPtr, tm$externalPtr,
        features$externalPtr,
        by.column)
}

##' Simulate a MSA given a tree model and HMM.
##'
##' Simulates a multiple sequence alignment of specified length.  Deals
##' with base-substitution only, not indels.  If one tree model is given,
##' simply simulates a sequence from this model.  If an HMM is provided,
##' then the mod parameter should be a list of tree models with the same
##' length as the number of states in the HMM.
##' @param object An object of type \code{tm} (or a list of these objects)
##' describing the phylogenetic model from which to simulate.  If it
##' is a list of tree models then an HMM should be provided to describe
##' transition rates between models.
##' @param nsim The number of columns in the simulated alignment.
##' @param seed A random number seed.  Either \code{NULL} (the default;
##' do not re-seed random  number generator), or an integer to be sent to
##' set.seed.
##' @param hmm an object of type HMM describing transitions between the
##' tree models across the columns of the alignment.
##' @param get.features (For use with hmm).  If \code{TRUE}, return object will
##' be a list of length two.  The first element will be the alignment, and the
##' second will be an object of type \code{feat} describing the path through
##' the phylo-hmm in the simulated alignment.  
##' @param pointer.only (Advanced use only). If TRUE, return only a pointer
##' to the simulated alignment.  Possibly useful for very (very) large
##' alignments.  Cannot be used if \code{get.features==TRUE}.
##' @param ... Currently not used (for S3 compatibility)
##' @return An object of type MSA containing the simulated alignment.
##' @keywords msa hmm
##' @export
simulate.msa <- function(object, nsim, seed=NULL, hmm=NULL, get.features=FALSE,
                         pointer.only=FALSE, ...) {
  nsites <- nsim
  mod <- object
  check.arg(get.features, "get.features", "logical", null.OK=FALSE, min.length=1L,
            max.length=1L)
  if (!is.null(seed)) set.seed(seed)
  if (get.features && (!is.null(hmm)) && pointer.only) 
    warning("pointer.only cannot be TRUE with get.features==TRUE; using pointer.only=FALSE")

  nstate <- 1L
  if (!is.null(hmm)) {
    nstate <- nstate.hmm(hmm)
    hmm <- (as.pointer.hmm(hmm))$externalPtr
  }
  if (is.tm(mod)) {
    tmlist <- list(mod)
  } else tmlist <- mod
  nmod <- length(tmlist)
  for (i in 1:nmod) {
    if (!is.tm(tmlist[[i]]))
      stop("mod should be a list of tree models (one for every state of HMM)")
    tmlist[[i]] <- (as.pointer.tm(tmlist[[i]]))$externalPtr
  }
  if (nstate != nmod) 
    stop("number of states in HMM (", nstate, ") does not match number of models (", nmod,")")

  if ((!is.null(hmm)) && get.features) {
    temp <- list()
    temp[[1]] <- .Call("rph_msa_base_evolve", tmlist, nsites, hmm, get.features)
    result <- list()
    msa <- .makeObj.msa()
    msa$externalPtr <- .Call("rph_msa_base_evolve_struct_get_msa", temp[[1]])
    result$msa <- from.pointer.msa(msa)
    features <- .makeObj.feat()
    features$externalPtr <- .Call("rph_msa_base_evolve_struct_get_labels", temp[[1]], nsites)
    result$feats <- as.data.frame.feat(features)
    return(result)
  }
  msa <- .makeObj.msa()
  msa$externalPtr <- .Call("rph_msa_base_evolve", tmlist, nsites, hmm, get.features)
  if (pointer.only == FALSE) msa <- from.pointer.msa(msa)
  msa
}


##' Sample columns from an MSA
##' @param x An object of type \code{msa}
##' @param size The number of columns to sample
##' @param replace Whether to sample with replacement
##' @param prob A vector of probability weights for sampling each column;
##' \code{prob=NULL} implies equal probability for all columns.  Probabilities
##' need not sum to one but should be non-negative and can not all be zero.
##' @param pointer.only If \code{TRUE}, return only a pointer to an alignment
##' object stored in C (useful for large objects; advanced use only).
##' @return An object of type \code{msa} with columns randomly
##' re-sampled from the original
##' @note This function is implemented using R's sample function in
##' conjunction with "[.msa".  It will not alter the value of x even if it
##' is stored as a pointer.
##' @S3method sample msa
##' @keywords msa
##' @export
sample.msa <- function(x, size, replace=FALSE, prob=NULL, pointer.only=FALSE) {
  check.arg(size, "size", "integer", null.OK=FALSE)
  check.arg(replace, "replace", "logical", null.OK=FALSE)
  if (!is.null(prob)) prob <- rep(prob, length.out=ncol.msa(x))
  if (size > ncol.msa(x) && replace==FALSE)
    stop("cannot sample more columns than in msa unless replace=TRUE")
  x[,sample(1:ncol.msa(x), size, replace=replace, prob=prob),
    pointer.only=pointer.only]
}

##' Extract fourfold degenerate sites from an MSA object
##' @param x An object of type \code{msa}
##' @param features an object of type \code{feat}.  Should have defined coding regions
##' with feature type "CDS"
##' @return an unordered msa object containing only the sites which are
##' fourfold degenerate
##' @note If x is stored as a pointer, it will be 
##' reduced to four-fold degenerate sites, so the original alignment will be lost.
##' Use get4d.msma(copy.msa(x), features) to avoid this behavior.
##' @note For very large MSA objects it is more efficient to use the do.4d option
##' in the read.msa function instead.
##' @export
get4d.msa <- function(x, features) {
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)
  if (is.null(features$externalPtr) && sum(features$feature=="CDS")==0L) 
    stop("features has no features labelled \"CDS\"... cannot extract 4d sites")
  if (is.null(features$externalPtr))
    features <- as.pointer.feat(features)
  x$externalPtr <- .Call("rph_msa_reduce_to_4d",
                         x$externalPtr,
                         features$externalPtr)
  x <- from.pointer.msa(x)
}


##' Extract features from an MSA object
##'
##' Returns the subset of the MSA which appears in the features object.
##' @param x An object of type MSA
##' @param features An object of type \code{features} denoting the regions
##' of the alignment to extract.
##' @param do4d If \code{TRUE}, then some elements of features must have type "CDS", and only
##' fourfold-degenerate sites will be extracted.
##' @param pointer.only If \code{TRUE}, return only a pointer to an object
##' stored in C (useful for large alignments; advanced use only)
##' @return An msa object containing only the regions of x
##' appearing in the features object.
##' @note If x was loaded with \code{pointer.only==TRUE}, then x
##' will be modified to the return value of the function.
##' Use \code{extract.feature.msa(copy.msa(x),...)} if you don't want this behavior!
##' @seealso \code{sub.msa}, \code{[.msa}
##' @keywords msa features
##' @export
extract.feature.msa <- function(x, features, do4d=FALSE, pointer.only=FALSE) {
  if (!is.ordered.msa(x))
    stop("extract.feature.msa requires ordered alignment")
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)

  if (do4d) {
    if (sum(features$feature=="CDS")==0L) 
      stop("features has no elements of type \"CDS\"... cannot extract 4d sites")
    rv <- get4d.msa(x, features)
  } else {
    
    if (is.null(features$externalPtr)) 
      features <- as.pointer.feat(features)
    rv <- .makeObj.msa()
    rv$externalPtr <- .Call("rph_msa_extract_feature",
                            x$externalPtr,
                            features$externalPtr)
  }
  if (!is.null(rv$externalPtr) && !pointer.only)
    rv <- from.pointer.msa(rv)
  rv$is.ordered=FALSE
  rv
}


##' Concatenate msa objects
##'
##' If the MSAs do not contain the same set of sequences, the sequences
##' will be added to each MSA and filled with missing data.  The order
##' of sequences is taken from the first MSA, and sequences are added to
##' this as necessary.
##' @param msas A list of MSA objects to concatenate together.
##' @param ordered If FALSE, disregard the order of columns in the combined
##' MSA.
##' @param pointer.only (Advanced use only, for very large MSA objects) If
##' TRUE, return object will be a pointer to an object stored in C.
##' @return An object of type MSA
##' @note None of the msas passed to this function will be altered, even if
##' they are stored as pointers to objects in C.
##' @keywords msa
##' @export
concat.msa <- function(msas, ordered=FALSE, pointer.only=FALSE) {
  # have to do a little dance to make sure this behaves OK if
  # some msas are empty
  isZero <- logical(length(msas))
  for (i in 1:length(msas)) {
    if (is.null(msas[[i]]$externalPtr)) 
      msas[[i]] <- as.pointer.msa(msas[[i]])
    isZero[i] <- (ncol.msa(msas[[i]]) == 0L)
  }
  if (sum(isZero) > 0L) 
    msas[isZero] <- NULL
  if (length(msas) == 0L) return(NULL)

  aggMsa <- copy.msa(msas[[1]])
  if (is.null(aggMsa$externalPtr))
    aggMsa <- as.pointer.msa(aggMsa)
  if (length(msas) >= 2L) {
    for (i in 2:length(msas)) {
      currMsa <- msas[[i]]
      if (is.null(currMsa$externalPtr))
        currMsa <- as.pointer.msa(currMsa)
      aggMsa$externalPtr <- .Call("rph_msa_concat",
                                  aggMsa$externalPtr,
                                  currMsa$externalPtr)
    }
  }
  if (pointer.only == FALSE) 
    aggMsa <- from.pointer.msa(aggMsa)
  aggMsa
}


##' Split an MSA by feature
##' @param x An object of type \code{msa}
##' @param f An object of type \code{feat}
##' @param drop Not currently used
##' @param pointer.only If \code{TRUE}, returned list elements are pointers to
##' objects stored in C (advanced use only).
##' @param ... Not currently used
##' @return A list of msa objects, representing the sub-alignments for
##' each element in f
##' @note If f is stored as a pointer to an object in C, its values will
##' be altered by this function.  Use
##' \code{split.by.feature.msa(x, copy.feat(f), ...)} to avoid this behavior!
##' @note x will not be altered even if it is stored as a pointer to an
##' object in C.
##' @keywords msa features
##' @export
split.by.feature.msa <- function(x, f, drop=FALSE, pointer.only=FALSE, ...) {
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  
  if (is.null(x$externalPtr)) x <- as.pointer.msa(x)
  if (is.null(f$externalPtr)) f <- as.pointer.feat(f)
  l <- list()
  l$externalPtr <- .Call("rph_msa_split_by_gff", x$externalPtr,
                         f$externalPtr)
  lst.size <- .Call("rph_list_len", l$externalPtr)
  rv <- list()
  for (i in 1:lst.size) {
    cat(i, "\n")
    m <- .makeObj.msa()
    m$externalPtr <- .Call("rph_msaList_get", l$externalPtr, i)
    if (pointer.only) {
      rv[[i]] <- m
    } else  {
      rv[[i]] <- from.pointer.msa(m)
    }
  }
  rv
}


##' Get informative regions of an alignment
##' @param x An object of type \code{msa}.
##' @param min.numspec The minimum number of species with non-missing data
##' required for an alignment column to be considered informative.
##' @param spec A character vector of species names, or an integer vector
##' of species indices.  Only data in
##' the named species count towards deciding if a site is informative.  The
##' default value of \code{NULL} implies use all species in the alignment.
##' @param refseq Defines the frame of reference for the return value.  A
##' value of \code{1} means use the first sequence in the alignment.  A value
##' of \code{0} or \code{NULL} means use the reference frame of the entire
##' alignment.  Valid values are integers from \code{seq(0..ncol.msa(x))},
##' \code{NULL}, or the name of a species in the alignment.
##' @param gaps.inf Logical value indicating whether a gap should be considered
##' informative.  The default value of \code{FALSE} indicates that gaps as
##' well as missing data are not counted as informative.
##' @return An object of type \code{feat} indicating the regions of the
##' alignment which meet the informative criteria.  Note that unless
##' \code{refseq==0 || refseq==NULL}, columns with gaps in the reference
##' sequence will be ignored, and will fall in "informative" or "uninformative"
##' features based on the informativeness of neighboring columns.
##' @note If the msa object has an idx.offset, it is assumed to be a coordinate
##' offset for the first species in the alignment.  So the idx.offset will
##' be added to the coordinates in the returned features object only if
##' \code{refseq==1}.
##' @note This function will not alter the value of x even if it is stored as
##' a pointer.
##' @keywords msa
##' @export
informative.regions.msa <- function(x, min.numspec, spec=NULL, refseq=1,
                                    gaps.inf=FALSE) {
  numspec <- nrow.msa(x)
  check.arg(min.numspec, "min.numspec", "integer", null.OK=FALSE)
  check.arg(gaps.inf, "gaps.inf", "logical", null.OK=FALSE)
  if (min.numspec <= 0L || min.numspec > numspec)
    stop("min.numspec expected to be between 1 and ", numspec)
  if (!is.null(spec)) {
    if (is.integer(spec)) {
      check.arg(spec, "spec", "integer", null.OK=TRUE, min.length=1L,
                max.length=numspec)
      if (sum(spec <= 0 | spec > numspec) > 0L)
          stop("expected spec values between 1 and ", numspec)
    } else {
      check.arg(spec, "spec", "character", null.OK=TRUE, min.length=1L,
                max.length=nrow.msa(x))
      intspec <- as.integer(sapply(spec, function(s) {which(s==names(x))}))
      if (sum(is.na(intspec)) > 0L)
        stop("don't know species names ", spec[is.na(intspec)])
      spec <- intspec
    }
  }
  if (is.character(refseq)) {
    check.arg(refseq, "refseq", "character", null.OK=FALSE)
    intrefseq <- which(refseq==names(msa))
    if (length(intrefseq) == 0L)
      stop("don't know refseq name ", refseq)
    refseq <- intrefseq
  } else if (!is.null(refseq)) {
    check.arg(refseq, "refseq", "integer", null.OK=FALSE)
    if (refseq < 0L || refseq > numspec)
      stop("expected refseq between 0 and ", numspec)
  } else refseq <- 0

  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)

  feats <- .makeObj.feat()
  feats$externalPtr <- .Call("rph_msa_informative_feats",
                             x$externalPtr, min.numspec, spec, refseq,
                             gaps.inf)
  as.data.frame.feat(feats)
}

