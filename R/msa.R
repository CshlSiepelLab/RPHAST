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
  result$externalPtr <- .Call("rph_msa_create_copy", msa$externalPtr)
  result
}


##' @export
is.msa <- function(msa) {
  if (is.null(msa$externalPtr)) {
    
  }
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
##' is stored as a pointer, the object can only be copied with copy.msa().
##' See examples below.
##' @title MSA Objects
##' @param seqs a character vector containing sequences, one per sample
##' @param names a character vector identifying the sample name for each
##' sequence
##' @param alphabet a character string containing valid non-missing character
##' states
##' @param is.ordered a logical indicating whether the alignment columns
##' are stored in order.  If NULL, assume columns are ordered.
##' @param offset an integer giving the offset of coordinates for the
##' reference sequence from the beginning of the chromosome.  Not used
##' if is.ordered==FALSE.
##' @param pointer.only a boolean indicating whether MSA should be stored by
##' reference (see Details)
##' @useDynLib rphast
##' @aliases is.msa
##' @export
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
##' length to return.
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

##' The number of informative columns in an alignment
##' @param msa An object of type \code{msa}
##' @return The number of "informative" columns in the msa.  An informative
##' column has at least two non-missing and non-gap characters.
##' @export
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
##' @keywords format
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
##' @keywords write
##' @keywords msa
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
##' @keywords summary
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
##' @keywords alignment
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
##' @param filename
##' @param format input file format: one of "FASTA", "MAF", "SS", "PHYLIP",
##' "MPM", must be correctly specified.
##' @param alphabet the alphabet of non-missing-data chraracters in the
##' alignment.  Determined automatically from the alignment if not given.
##' @param gff a GFF object.  If provided, the return value will only
##' contain portions of the alignment which fall within a feature in the GFF.
##' The alignment will not be ordered.
##' The loaded regions can be further constrained with the do.4d or
##' do.cats options.
##' @param do.4d Logical.  If \code{TRUE}, the return value will contain only
##' the columns corresponding to four-fold degenerate sties.  Requires
##' gff to be specified.
##' @param ordered Logical.  If \code{FALSE}, the MSA object may not retain
##' the original column order.
##' @param tuple.size Integer.  If given, and if pointer.only is \code{TRUE},
##' MSA will be stored in sufficient statistics format, where each tuple
##' contains tuple.size consecutive columns of the alignment.
##' @param do.cats Character vector.  If given, and if gff is specified,
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
##' @note If the input is in "MAF" format and a gff is given, the
##' resulting alignment will be stripped of gaps in the reference (1st)
##' sequence.
##' @return an MSA object.  
##' @seealso \code{\link{msa}}, \code{\link{read.gff}}
##' @export
read.msa <- function(filename,
                     format=c(guess.format.msa(filename), "FASTA")[1],
                     alphabet=NULL,                     
                     gff=NULL,
                     do.4d=FALSE,
                     ordered=ifelse(do.4d || !is.null(gff), FALSE, TRUE),
                     tuple.size=(if(do.4d) 3 else NULL),
                     do.cats=NULL,
                     refseq=NULL,
                     offset=0,
                     seqnames=NULL, discard.seqnames=NULL,
                     pointer.only=FALSE) {

  check.arg(filename, "filename", "character", null.OK=FALSE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(alphabet, "alphabet", "character", null.OK=TRUE)
  check.arg(gff, "gff", null.OK=TRUE, min.length=NULL, max.length=NULL)
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
    if (is.null(gff))
      stop("gff needs to be specified with do.4d")
    if (tuple.size != 3)
      stop("tuple.size must be 3 if do.4d==TRUE")
  }

  if (!is.null(do.cats) && is.null(gff))
    stop("gff required with do.cats")

  if (!is.null(gff)) {
    if (is.null(gff$externalPtr))
      gff <- as.pointer.gff(gff)
  }

  msa <- .makeObj.msa()
  msa$externalPtr <- .Call("rph_msa_read", filename, format,
                           gff$externalPtr, do.4d, alphabet,
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
##' @export
complement <- function(x) {
  chartr("ACGTacgt", "TGCAtgca", x)
}


##' Reverse complement a multiple sequence alignment
##' @param msa An object of type \code{msa}.
##' @return The reverse complement of msa.
##' @export
reverse.complement <- function(msa) {
  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  rv <- .makeObj.msa()
  rv$externalPtr <- .Call("rph_msa_reverse_complement", msa$externalPtr)
  from.pointer.msa(rv)
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
##' @S3method sub msa
##' @export
sub.msa <- function(msa, seqs=NULL, keep=TRUE, start.col=NULL, end.col=NULL,
                   refseq=NULL) {
  check.arg(keep, "keep", "logical", null.OK=FALSE)
  check.arg(seqs, "seqs", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(start.col, "start.col", "integer", null.OK=TRUE)
  check.arg(end.col, "end.col", "integer", null.OK=TRUE)
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  
  result <- .makeObj.msa()
  if (is.null(msa$externalPtr)) {
    msa <- as.pointer.msa(msa)
    copy.to.R<-TRUE
  } else copy.to.R <- FALSE
  result$externalPtr <- .Call("rph_msa_sub_alignment",
                              msa$externalPtr, seqs, keep,
                              start.col, end.col, refseq)

  if (copy.to.R) 
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
##' @param msa MSA object
##' @param strip.mode Determines which gaps to strip.  See Details
##' @note if MSA is stored as a pointer, then changes will occur to the
##' input alignment
##' @return an MSA object, with gaps stripped according to strip.mode
##' @export
strip.gaps.msa <- function(msa, strip.mode=1) {
  names <- NULL
  nseq <- NULL
  if (is.null(msa$externalPtr)) {
    names <- names.msa(msa)
    nseq <- nrow.msa(msa)
    msa <- as.pointer.msa(msa)
    copy.to.R <- TRUE
  } else copy.to.R <- FALSE
  for (s in strip.mode) {
    if (s=="all.gaps" || s=="any.gaps")
      msa$externalPtr <- .Call("rph_msa_strip_gaps", msa$externalPtr, 0, s)
    else {
      if (!is.character(s)) {
        if (is.null(nseq)) nseq <- nrow.msa(msa)
        if (as.integer(s) != s || s <=0 || s>nseq)
          stop(cat("invalid sequence index", s))
        w <- s
      } else {
        if (is.null(names))
          names <- names.msa(msa)
        w <- which(names==s)
        if (is.null(w))
          stop(cat("no sequence with name", s))
      }
      msa$externalPtr <- .Call("rph_msa_strip_gaps", msa$externalPtr, w, NULL)
    }
  }
  if (copy.to.R) 
    msa <- from.pointer.msa(msa)
  msa
}



##' Extract, replace, reorder MSA
##'
##' Treat multiple sequence alignment as a matrix where each row
##' corresponds to a sequence for one species, and each column
##' is one position aligned across all species.
##'
##' The bracket notation can return a subset of the alignment,
##' or re-order rows and columns.
##'
##' @param msa An object of type \code{msa}
##' @param rows A numeric vector of sequence indices,
##' character vector (containing sequence name), or
##' logical vector (containing sequences to keep).  If logical vector it
##' will be recycled as necessary to the same length as \code{nrow.msa(msa)}.
##' @param cols A numeric vector of alignment columns, or a logical vector
##' containing columns to keep.  If logical vector it will be recycled as
##' necessary to the same lenth as \code{ncol.msa}.  Note that these are
##' in coordinates with respect to the entire alignment.  msa$idx.offset
##' is ignored here.
##' @seealso \code\{link{sub.msa}} which can subset columns based on genomic
##' coordinates.
##' @seealso \code{link{extract.feature.msa} which can subset based on
##' genomic coordinates denoted in a features object.
##' @S3method "[" msa
##' @export "["
"[.msa" <- function(msa, rows, cols) {
  if (!missing(rows)) {
    if (is.null(rows)) stop("rows cannot be empty")
  } else rows=NULL
  if (!missing(cols)) {
    if (is.null(cols)) stop("cols cannot be empty")
  } else cols=NULL

  if (!is.null(rows)) {
    # if rows are given by names, convert to integer
    if (is.character(rows)) {# names are given
      names <- names.msa(msa)
      rows <- as.numeric(sapply(rows, function(x) {which(x ==  names)}))
      if (sum(is.na(rows)) > 0L)
        stop("unknown names in first dimension subset")
    }
  }

  # check if arguments are given as logicals.
  if (is.logical(rows)) {
    cat("is.logical rows")
    rows = which(rep(rows, length.out = nrow.msa(msa)))
  }
  if (is.logical(cols)) 
    cols = which(rep(cols, length.out = ncol.msa(msa)))

  # if msa is stored in R, sampling rows is easier and more efficient to do here
  if (!is.null(rows) && is.null(msa$externalPtr)) {
    msa$names <- msa$names[rows]
    msa$seqs <- msa$seqs[rows]
    rows=NULL
  }
  if (is.null(rows) && is.null(cols)) return(msa)
  check.arg(rows, "rows", "integer", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(cols, "cols", "integer", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  rv <- .makeObj.msa()
  rv$externalPtr <- .Call("rph_msa_square_brackets", msa$externalPtr,
                          rows, cols)
  from.pointer.msa(rv)
}




##' Likelihood of an alignment given a tree model
##' @title MSA Likelihood
##' @param msa An object of class msa representing the multiple alignment
##' @param tm An object of class TM representing the tree and model of
##' substitution
##' @param by.column Logical indicating whether to get likelihoods for
##' each alignment column.  If FALSE, returns total likelihood
##' @return Either the likelihood of the entire alignment (if by.column is
##' \code{FALSE}, or a numeric vecotr giving the likelihood of each
##' column in the alignment
##' @export
likelihood.msa <- function(msa, tm, by.column=FALSE) {
  check.arg(by.column, "by.column", "logical", null.OK=FALSE)
  if (is.null(msa$externalPtr)) 
    msa <- as.pointer.msa(msa)
  tm <- as.pointer.tm(tm)
  if (by.column && !is.ordered.msa(msa))
    warning("by.column may not be a sensible option for unordered MSA")
  .Call("rph_msa_likelihood", msa$externalPtr, tm$externalPtr, by.column)
}

##' Simulate a MSA given a tree model and HMM.
##'
##' Simulates a multiple sequence alignment of specified length.  Deals
##' with base-substitution only, not indels.  If one tree model is given,
##' simply simulates a sequence from this model.  If an HMM is provided,
##' then the mod parameter should be a list of tree models with the same
##' length as the number of states in the HMM.
##' @param mod A tree model or a list of tree models from which to simulate.
##' @param nsites The number of columns in the simulated alignment.
##' @param hmm an object of type HMM describing transitions between the
##' tree models across the columns of the alignment.
##' @param pointer.only (Advanced use only). If TRUE, return only a pointer
##' to the simulated alignment.  Possibly useful for very (very) large
##' alignments.
##' @return An object of type MSA containing the simulated alignment.
##' @export
simulate.msa <- function(mod, nsites, hmm=NULL, pointer.only=FALSE) {
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
  msa <- .makeObj.msa()
  msa$externalPtr <- .Call("rph_msa_base_evolve", tmlist, nsites, hmm)
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
##' @return An object of type \code{msa} with columns randomly
##' re-sampled from the original
##' @S3method sample msa
##' @export
sample.msa <- function(x, size, replace=FALSE, prob=NULL) {
  check.arg(size, "size", "integer", null.OK=FALSE)
  check.arg(replace, "replace", "logical", null.OK=FALSE)
  if (!is.null(prob)) prob <- rep(prob, length.out=ncol.msa(x))
  if (size > ncol.msa(x) && replace==FALSE)
    stop("cannot sample more columns than in msa unless replace=TRUE")
  x[,sample(1:ncol.msa(x), replace=replace, prob=prob)]
}

##' Extract fourfold degenerate sites from an MSA object
##' @param msa An object of type MSA
##' @param gff an object of type GFF.  Should have defined coding regions
##' with feature type "CDS"
##' @return an unordered msa object containing only the sites which are
##' fourfold degenerate
##' @note if original msa is stored as a pointer it will be destroyed.  For
##' very large MSA objects it is more efficient to use the do.4d option
##' in the read.msa function instead.
##' @export
get4d.msa <- function(msa, gff) {
  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)
  if (is.null(gff$externalPtr) && sum(gff$feature=="CDS")==0L) 
    stop("gff has no features labelled \"CDS\"... cannot extract 4d sites")
  if (is.null(gff$externalPtr))
    gff <- as.pointer.gff(gff)
  msa$externalPtr <- .Call("rph_msa_reduce_to_4d",
                           msa$externalPtr,
                           gff$externalPtr)
  msa <- from.pointer.msa(msa)
}


##' Extract features from an MSA object
##'
##' Returns the subset of the MSA which appears in the features (GFF) object.
##' @param msa An object of type MSA
##' @param gff A GFF object denoting the regions of the alignment to extract.
##' @param do4d If TRUE, then gff must have features of type "CDS", and only
##' fourfold-degenerate sites will be extracted.
##' @return An msa object containing only the regions of the msa
##' appearing in the GFF object.
##' @note if input msa is stored as a pointer it will be destroyed.
##' @export
extract.feature.msa <- function(msa, gff, do4d=FALSE) {
  if (!is.ordered.msa(msa))
    stop("extract.feature.msa requires ordered alignment")
  if (is.null(msa$externalPtr))
    msa <- as.pointer.msa(msa)

  if (do4d) {
    if (sum(gff$feature=="CDS")==0L) 
      stop("gff has no features labelled \"CDS\"... cannot extract 4d sites")
    rv <- get4d.msa(msa, gff)
  } else {
    
    if (is.null(gff$externalPtr)) 
      gff <- as.pointer.gff(gff)
    rv <- .makeObj.msa()
    rv$externalPtr <- .Call("rph_msa_extract_feature",
                            msa$externalPtr,
                            gff$externalPtr)
  }
  if (!is.null(rv$externalPtr))
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
##' @param msa A list of MSA objects to concatenate together.
##' @param ordered If FALSE, disregard the order of columns in the combined
##' MSA.
##' @param pointer.only (Advanced use only, for very large MSA objects) If
##' TRUE, return object will be a pointer to an object stored in C.
##' @param return An object of type MSA
##' @export
concat.msa <- function(msas, ordered=FALSE, pointer.only=FALSE) {
  aggMsa <- copy.msa(msas[[1]])
  if (is.null(aggMsa$externalPtr))
    aggMsa <- as.pointer.msa(aggMsa)
  for (i in 2:length(msas)) {
    currMsa <- msas[[i]]
    if (is.null(currMsa$externalPtr))
      currMsa <- as.pointer.msa(currMsa)
    aggMsa$externalPtr <- .Call("rph_msa_concat",
                                aggMsa$externalPtr,
                                currMsa$externalPtr)
  }
  if (pointer.only == FALSE) 
    aggMsa <- from.pointer.msa(aggMsa)
  aggMsa
}


##' Split an MSA by feature
##' @param msa An object of type \code{msa}
##' @param gff A feature object
##' @return A list of msa objects, representing the sub-alignments for
##' each feature in the gff object
##' @export
split.by.feature.msa <- function(msa, gff) {
  if (is.null(msa$externalPtr)) msa <- as.pointer.msa(msa)
  if (is.null(gff$externalPtr)) gff <- as.pointer.gff(gff)
  l <- list()
  l$externalPtr <- .Call("rph_msa_split_by_gff", msa$externalPtr,
                         gff$externalPtr)
  lst.size <- .Call("rph_list_len", l$externalPtr)
  rv <- list()
  for (i in 1:lst.size) {
    cat(i, "\n")
    m <- .makeObj.msa()
    m$externalPtr <- .Call("rph_msaList_get", l$externalPtr, i)
    rv[[i]] <- from.pointer.msa(m)
  }
  rv
}
