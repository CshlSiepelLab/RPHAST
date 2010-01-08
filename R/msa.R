# don't call explicitly; this is the registered finalizer for MSA external
# pointers
#' @nord
msa.free <- function(extMsaPtr) {
  .Call("rph_msa_free", extMsaPtr)
}


# make barebones msa obj
#' @nord
msa.makeObj <- function() {
  msa <- list()
  class(msa) <- "msa"
  msa
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
##' If pointer.only==FALSE, the MSA object will be stored in R and can be
##' viewed and modified by base R code as well as RPHAST functions.
##' Setting pointer.only=TRUE will cause the MSA object to be stored by
##' reference, as an external pointer to an object created by C code.  This
##' may be necessary to improve performance, but the object can then only
##' be viewed/manipulated via RPHAST functions.  Furthermore, if an object
##' is stored as a pointer, the object can only be copied with msa.copy().
##' See examples below.
##' @title MSA Objects
##' @param seqs a character vector containing sequences, one per sample
##' @param names a character vector identifying the sample name for each
##' sequence
##' @param alphabet a character string containing valid non-missing character
##' states
##' @param pointer.only a boolean indicating whether MSA should be stored by
##' reference (see Details)
##' @useDynLib rphast
##' @export
msa.new <- function(seqs, names=NULL, alphabet="ACGT",
                    pointer.only=FALSE) {
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

  # TODO: if alphabet non-null, check that seqs only contains those chars?

  msa <- msa.makeObj()

  if (pointer.only) {
    msa$externalPtr <- .Call("rph_msa_new",
                             seqsP=seqs,
                             namesP=names,
                             nseqsP=length(seqs),
                             lengthP=seqlen,
                             alphabetP=alphabet)
    reg.finalizer(msa$externalPtr, msa.free)
  } else {
    msa$seqs <- seqs
    if (! is.null(names)) msa$names <- names
    if (! is.null(alphabet)) msa$alphabet <- alphabet
  }
  msa
}

##' Returns the length of sequence in an MSA alignment.
##' @title MSA Sequence Length.
##' @param msa an MSA object
##' @param refseq character vector giving name(s) of sequence whose
##' length to return.
##' @return an integer vector containing the length of the named sequences.
##' If refseq is NULL, returns the number of columns in the alignment.
##' @keywords msa
##' @seealso \code{\link{msa.new}}
##' @export
msa.seqlen <- function(msa, refseq=NULL) {
  if (is.null(msa$externalPtr) && is.null(refseq)) {
    return(nchar(msa$seqs[1]))
  }
  if (is.null(msa$externalPtr))
    msa <- msa.to.pointer(msa)
  if (is.null(refseq))
    return(.Call("rph_msa_seqlen", msa$externalPtr, NULL))
  result <- integer(length(refseq))
  for (i in 1:length(refseq)) 
    result[i] <- .Call("rph_msa_seqlen", msa$externalPtr, refseq[i])
  result
}

##' Returns the number of sequence in an MSA alignment.
##' @title MSA Number of Sequences
##' @param msa an MSA object
##' @return an integer containing the number of sequences in an alignment.
##' @keywords msa
##' @seealso \code{\link{msa.new}}
##' @export
msa.nseq <- function(msa) {
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
##' @keywords format
##' @export
msa.validFormatStr <- function(format) {
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
##' @export
msa.idx.offset <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_idxOffset", msaP=msa$externalPtr, tagP=msa))
  msa$idx.offset
}

##' Returns the alphabet used by an MSA object.
##' @title MSA Alphabet
##' @param msa an MSA object
##' @return the valid non-missing-data characters for an MSA object.
##' @export
msa.alphabet <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_alphabet", msaP=msa$externalPtr))
  msa$alphabet
}

##' Determines if an MSA object represents an ordered alignment.
##' @title MSA is Ordered?
##' @param msa an MSA object
##' @return a boolean indicating whether the columns are in order
##' @export
msa.is.ordered <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_isOrdered", msaP=msa$externalPtr))

  if (is.null(msa$is.ordered)) return (TRUE)
  msa$is.ordered
}


##' Returns the sequence names for an MSA object.
##' @title MSA Sequence Names
##' @param msa an MSA object
##' @return a character vector giving the names of the sequences, or
##' NULL if they are not defined
##' @export
msa.names <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_seqNames", msaP=msa$externalPtr))
  msa$names
}


##' Take an MSA stored by reference and return one stored in R
##' @title MSA From Pointer
##' @param src an MSA object stored by reference
##' @return an MSA object stored in R.  If src is already stored in R,
##' returns a copy of the object.
##' @seealso \code{\link{msa.new}} for details on MSA storage options.
##' @export
msa.from.pointer <- function(src) {
  if (is.null(src$externalPtr)) return(src)
  seqs <- .Call("rph_msa_seqs", src$externalPtr)
  names <- .Call("rph_msa_seqNames", src$externalPtr)
  alphabet <- .Call("rph_msa_alphabet", src$externalPtr)
  msa.new(seqs, names, alphabet, pointer.only=FALSE)
}


##' Take an MSA stored in R and return one stored by reference
##' @title MSA To Pointer
##' @param src an MSA object stored by value in R
##' @return an MSA object stored by reference as a pointer to an object
##' created in C.
##' @seealso \code{\link{msa.new}} for details on MSA storage options.
##' @export
msa.to.pointer <- function(src) {
  if (!is.null(src$externalPtr)) return(src)
  msa.new(src$seq, src$names, src$alphabet, pointer.only=TRUE)
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
##' @keywords write
##' @keywords msa
##' @export
msa.write <- function(msa, filename=NULL, format="FASTA", pretty.print=FALSE) {
  #checks
  check.arg(filename, "filename", "character", null.OK=TRUE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(pretty.print, "pretty.print", "logical", null.OK=FALSE)
  if (! msa.validFormatStr(format)) {
    stop(paste("invalid MSA FORMAT \"", format, "\"", sep=""))
  }
  
  if (is.null(msa$externalPtr)) {
    printMsa <- msa.to.pointer(msa)
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
##' @keywords summary
##' @seealso \code{\link{print.msa}}
##' @export
summary.msa <- function(object, ..., print.seq=FALSE, format="FASTA",
                        pretty.print=FALSE) {
  msa <- object
  check.arg(print.seq, "print.seq", "logical", null.OK=FALSE)
  # format and pretty.print are checked in msa.write

  cat(paste("msa object with", msa.nseq(msa), "sequences and",
            msa.seqlen(msa),"columns"))
  cat("\n")

  if (!is.null(format)) {
    print.seq=TRUE
  }
  
  pointer <- msa$externalPtr
  names <- msa.names(msa)
  alphabet <- msa.alphabet(msa)
  is.ordered <- msa.is.ordered(msa)
  idx.offset <- msa.idx.offset(msa)

  printMsa <- list()
  printed <- FALSE
  if (!is.null(names)) printMsa$names <- names
  if (!is.null(alphabet)) printMsa$alphabet <- alphabet
  if (!is.null(is.ordered)) printMsa$is.ordered <- is.ordered
  if (!is.null(idx.offset) && idx.offset!=0) printMsa$idx.offset <- idx.offset
  if (print.seq && is.null(pointer) && is.null(format)) {
    printMsa$seq <- msa$seq
    printed <- TRUE
  }
  
  print(printMsa, ...)

  if (!printed && print.seq) {
    cat("$seq\n")
    if (is.null(format)) format <- "FASTA"
      msa.write(msa, filename=NULL, format, pretty.print)
    cat("\n")
  }
}


##' Prints an MSA (multiple sequence alignment) object.
##'
##' Valid formats for printing are "FASTA", "PHYLIP", "MPM", and "SS".
##' See \code{\link{msa.validFormatStr}} for details on these formats.
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
##' @keywords alignment
##' @export
print.msa <- function(x, ..., print.seq=FALSE, format=NULL, pretty.print=FALSE) {
  check.arg(print.seq, "print.seq", "logical", null.OK=FALSE)
  # format and pretty.print are checked in msa.write

  summary.msa(x, print.seq=print.seq, format=format, pretty.print=pretty.print, ...)

  if (!print.seq && is.null(format)) {
    cat("(alignment output suppressed)\n")
    cat("\n")
  }
}



##' Reads an MSA from a file.
##' @title Reading an MSA Object
##' @param filename
##' @param format input file format: one of "FASTA", "MAF", "SS", "PHYLIP",
##' "MPM", must be correctly specified.
##' @param alphabet the alphabet of non-missing-data chraracters in the
##' alignment.  Determined automatically from the alignment if not given.
##' @param gff.only a GFF object.  If provided, the return value will only
##' contain portions of the alignment which fall within a feature in the GFF
##' @param do4d Logical.  If \code{TRUE}, the return value will contain only
##' the columns corresponding to four-fold degenerate sties.  Requires
##' gff.only to be specified.
##' @param ordered Logical.  If \code{FALSE}, the MSA object may not retain
##' the original column order.
##' @param tuple.size Integer.  If given, and if pointer.only is \code{TRUE},
##' MSA will be stored in sufficient statistics format, where each tuple
##' contains tuple.size consecutive columns of the alignment.
##' @param do.cats Character vector.  If given, and if gff.only is specified,
##' then only the types of features named here will be represented in the
##' returned alignment.
##' @param refseq Character string specifying a FASTA format file with a
##' reference sequence.  If given, the reference sequence will be
##' "filled in" whereever missing from the alignment.
##' @param pointer.only If \code{TRUE}, MSA will be stored by reference as
##' an external pointer to an object created by C code, rather than
##' directly in R memory.  This improves performance and may be necessary
##' for large alignments, but reduces functionality.  See
##' \code{\link{msa.new}} for more details on MSA object storage options.
##' @return an MSA object.  
##' @seealso \code{\link{msa.new}}, \code{\link{read.gff}}
##' @export
read.msa <- function(filename, format="FASTA",
                     alphabet=NULL,                     
                     gff.only=NULL,
                     do4d=FALSE,
                     ordered=ifelse(do4d || !is.null(gff.only), FALSE, TRUE),
                     tuple.size=(if(do4d) 3 else NULL),
                     do.cats=NULL,
                     refseq=NULL,
                     pointer.only=FALSE) {

  check.arg(filename, "filename", "character", null.OK=FALSE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(alphabet, "alphabet", "character", null.OK=TRUE)
  check.arg(gff.only, "gff.only", null.OK=TRUE, min.length=NULL, max.length=NULL)
  check.arg(do4d, "do4d", "logical", null.OK=FALSE)
  check.arg(ordered, "ordered", "logical", null.OK=FALSE)
  check.arg(tuple.size, "tuple.size", "integer", null.OK=TRUE)
  check.arg(do.cats, "do.cats", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  if (!msa.validFormatStr(format))
    stop(paste("invalid format string", format))
  if (!is.null(tuple.size) && tuple.size <= 0)
    stop("tuple.size should be integer >= 1")
  
  if (do4d) {
    if (!is.null(do.cats))
      stop("should not specify do.cats if do4d==TRUE")
    if (is.null(gff.only))
      stop("gff.only needs to be specified with do4d")
    if (tuple.size != 3)
      stop("tuple.size must be 3 if do4d==TRUE")
  }

  if (!is.null(do.cats) && is.null(gff.only))
    stop("gff.only required with do.cats")

  gff <- NULL
  (!is.null(gff.only)) 
    gff <- gff.only
  
  if (!is.null(gff)) {
    if (is.null(gff$externalPtr))
      gff <- gff.to.pointer(gff)
  }

  msa <- msa.makeObj()
  msa$externalPtr <- .Call("rph_msa_read", filename, format,
                           gff$externalPtr, do4d, alphabet,
                           tuple.size, refseq, ordered, do.cats)
  if (pointer.only == FALSE) {
    msa <- msa.from.pointer(msa)
    if (ordered==FALSE)
      msa$is.ordered = FALSE
  }
  msa
}


msa.compl.char <- function(ch) {
  b <- c("A", "C", "G", "T", "a", "c", "g", "t")
  c <- c("T", "G", "C", "A", "t", "g", "c", "a")
  w <- which(ch==b)
  if (is.null(w)) return(ch)
  c[w]
}

       
##' Get a subset of an alignment
##'
##' @title MSA Subset
##' @param msa MSA object
##' @param seqs The sequence names to keep (or to remove if keep is
##' \code{FALSE})
##' @param keep Whether to keep the named sequences or remove them
##' @param start.col the first column to keep (columns indices start at 1)
##' @param end.col the last column to keep (inclusive)
##' @param refseq the sequence in the alignment which determines the
##' coordinates for start.col and end.col.  If NULL, start.col and
##' end.col are column indices in the multiple alignment.
##' @return A new MSA object containing a subset of the original MSA.
##' @export
subMsa <- function(msa, seqs=NULL, keep=TRUE, start.col=NULL, end.col=NULL,
                   refseq=NULL) {
  check.arg(keep, "keep", "logical", null.OK=FALSE)
  check.arg(seqs, "seqs", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(start.col, "start.col", "integer", null.OK=TRUE)
  check.arg(end.col, "end.col", "integer", null.OK=TRUE)
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  
  result <- msa.makeObj()
  if (is.null(msa$externalPtr)) {
    msa <- msa.to.pointer(msa)
    copy.to.R<-TRUE
  } else copy.to.R <- FALSE
  result$externalPtr <- .Call("rph_msa_sub_alignment",
                              msa$externalPtr, seqs, keep,
                              start.col, end.col, refseq)
  if (copy.to.R) 
    result <- msa.from.pointer(result)
  result
}
