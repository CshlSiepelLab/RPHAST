# don't call explicitly; this is the registered finalizer for MSA external pointers
msa.free <- function(extMsaPtr) {
  .Call("rph_msa_free", extMsaPtr)
}


# make barebones msa obj
msa.makeObj <- function() {
  msa <- list()
  class(msa) <- "msa"
  msa
}


# make new MSA object given sequences and other information stored in R memory
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


msa.seqlen <- function(msa) {
  if (is.null(msa$externalPtr)) {
    return(nchar(msa$seqs[1]))
  }
  .Call("rph_msa_seqlen", msa$externalPtr)
}


msa.nseq <- function(msa) {
  if (is.null(msa$externalPtr)) {
    return(length(msa$seqs))
  }
  .Call("rph_msa_nseq", msa$externalPtr)
}


msa.validFormatStr <- function(format) {
  if (is.null(format)) return(NULL)
  result <- logical(length(format))
  for (i in 1:length(format)) {
    result[i] <- .Call("rph_msa_valid_fmt_str", format[i]);
  }
  result
}


msa.idx.offset <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_idxOffset", msaP=msa$externalPtr, tagP=msa))
  msa$idx.offset
}


msa.alphabet <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_alphabet", msaP=msa$externalPtr))
  msa$alphabet
}


msa.is.ordered <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_isOrdered", msaP=msa$externalPtr))

  if (is.null(msa$is.ordered)) return (TRUE)
  msa$is.ordered
}


msa.names <- function(msa) {
  if (!is.null(msa$externalPtr))
    return(.Call("rph_msa_seqNames", msaP=msa$externalPtr))
  msa$names
}


# take an MSA stored as a pointer and return one stored in R
msa.from.pointer <- function(src) {
  if (is.null(src$externalPtr)) return(src)
  seqs <- .Call("rph_msa_seqs", src)
  names <- .Call("rph_msa_seqNames", src)
  alphabet <- .Call("rph_msa_alphabet", src)
  msa.new(seqs, names, alphabet, pointer.only=FALSE)
}


# take an MSA stored in R and return one stored as a pointer
msa.to.pointer <- function(src) {
  if (!is.null(src$externalPtr)) return(src)
  msa.new(src$seq, src$names, src$alphabet, pointer.only=TRUE)
}


# write sequence to a file (or stdout if filename==NULL)
# If msa is not stored as pointer, then need to make temporary phast
# object with contents of msa, print, then delete.
# possible TODO: print directly from R in this case?
# Would be easy except for formatting.
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
  .Call("rph_msa_printSeq",
        msaP=printMsa$externalPtr,
        filenameP=filename,
        formatP=format,
        pretty.printP=pretty.print)
}


# print the msa.  Only print alignment if printSeq==TRUE.
print.msa <- function(msa, ..., printSeq=FALSE, format=NULL, pretty.print=FALSE) {
  check.arg(printSeq, "printSeq", "logical", null.OK=FALSE)
  # format and pretty.print are checked in msa.write

  cat(paste("msa object with", msa.nseq(msa), "sequences and ",
            msa.seqlen(msa),"columns"))
  cat("\n")

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
  if (printSeq && is.null(pointer) && is.null(format)) {
    printMsa$seq <- msa$seq
    printed <- TRUE
  }
  
  print(printMsa, ...)

  if (!printed && printSeq) {
    cat("$seq\n")
    if (is.null(format)) format <- "FASTA"
      msa.write(msa, filename=NULL, format, pretty.print)
  }

  if (!printSeq)
    cat("(alignment output suppressed)\n")
  cat("\n")
}


# same as print.msa but no option to print sequence
summary.msa <- function(msa, ...) {
  print.msa(msa, ...)
}



msa.read <- function(filename, format="FASTA",
                     alphabet=NULL,                     
                     gff.label=NULL, gff.only=NULL,
                     do4d=FALSE,
                     ordered=ifelse(do4d || !is.null(gff.only), FALSE, TRUE),
                     tuple.size=ifelse(do4d, 3, NULL),
                     category.map=NULL, doCats=NULL,
                     pointer.only=FALSE) {

  check.arg(filename, "filename", "character", null.OK=FALSE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(gff.label, "gff.label", null.OK=TRUE)
  check.arg(gff.only, "gff.only", null.OK=TRUE)
  check.arg(do4d, "do4d", "logical", null.OK=FALSE)
  check.arg(tuple.size, "tuple.size", "integer", null.OK=TRUE)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  check.arg(alphabet, "alphabet", "character", null.OK=TRUE)
  check.arg(category.map, "category.map", "cm", null.OK=TRUE)
  check.arg(ordered, "ordered", "logical", null.OK=FALSE)

  if (!msa.validFormatStr(format))
    stop(paste("invalid format string", format))
  if (tuple.size <= 0)
    stop("tuple.size should be integer >= 1")
  
  if (do4d) {
    if (!is.null(category.map)) stop("category.map should not be specified with do4d")
    if (is.null(gff.only) && is.null(gff.label))
      stop("one of gff.only or gff.label needs to be specified with do4d")
    category.map <- cm.new("NCATS=3; CDS 1-3")
    reverse.groups.tag <- "transcript_id"
  } else reverse.groups.tag <- NULL

  gff <- NULL
  if (!is.null(gff.only) && !is.null(gff.label)) {
    if (gff.only!=gff.label)
      stop("gff.only and gff.label cannot be different GFF objects")
    gff <- gff.only
  } else if (!is.null(gff.only)) {
    gff <- gff.only
  } else if (!is.null(gff.label))
    gff <- gff.label
  
  if (!is.null(gff)) {
    gffPtr <- .Call("rph_gff_get_ptr", gff)
    if (is.null(category.map)) {
      cmPtr <- .Call("rph_cm_new_from_gff", gffPtr)
    } else cmPtr <- .Call("rph_cm_new_from_R", cmPtr)
  }

  msa <- msa.makeObj()
  
  if (format != "MAF") {
    msa$externalPtr <- .Call("rph_msa_new_from_file", filename, format, alphabet)
    if (gff != NULL) {
      idx.offset <- msa.idx.offset(msa)
      if (!is.null(idx.offset) && idx.offset != 0) {
        gff$start <- gff$start - idx.offset
        gff$end <- gff$end - idx.offset
      }
      .Call("rph_msa_map_gff_coords", msa$externalPtr, gffPtr, -1, 0, 0)
      if (do4d) {
        .Call("rph_gff_group", gffPtr, reverse.groups.tag)
        .Call("rph_msa_reverse_compl_feats", msa$externalPtr, gffPtr, NULL)
      }
      .Call("rph_msa_label_categories", msa$externalPtr, gffPtr, cmPtr)
      # in any of these cases we need to compute sufficient statistics
      # (it seems that phast can only do sub-alignment based on category through SS
      if (!is.null(tuple.size) ||
          !is.null(gff.only) ||
          do4d ||
          (!ordered && pointer.only)) 
        .Call("rph_ss_from_msas", msaPtr, tuple.size, ordered, doCats, NULL, NULL, -1)
    }
  } else {
    msa$externalPtr <- .Call("rph_maf_read_cats", filename, NULL, tuple.size,
                             alphabet, gffPtr, cmPtr, -1, ordered,
                             reverse.groups.tag, NULL, FALSE, doCats)
  }

  if (do4d) .Call("rph_msa_reduce_to_4d", msa$externalPtr, cmPtr)

  if (pointer.only == FALSE)
    msa <- msa.extract.from.C(msa)
}


msa.compl.char <- function(ch) {
  b <- c("A", "C", "G", "T", "a", "c", "g", "t")
  c <- c("T", "G", "C", "A", "t", "g", "c", "a")
  w <- which(ch==b)
  if (is.null(w)) return(ch)
  c[w]
}

#msa.compl.string <- function(text) {
#  for (i in 1:length(text)) {
#    s <- text[i]
#    if (length(s)) 
#      for (j in 1:length(
       

