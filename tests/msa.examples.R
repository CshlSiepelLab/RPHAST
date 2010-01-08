library(rphast)

#' print.msa
# read in an MSA stored in R
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"))
print(m)
print(m, format="FASTA")
print(m, format="PHYLIP", pretty.print=TRUE)
#'
# read in an MSA stored by reference in C
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"),
             pointer.only=TRUE)
print(m)


#' summary.msa
# read in an MSA stored in R
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"))
summary(m)
#'
# read in an MSA stored by reference in C
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"),
             pointer.only=TRUE)
summary(m)

#' msa.validFormatStr
msa.validFormatStr(c("MAF", "SS", "PHYLIP", "MPM", "LAV", "FASTA",
                     "BAD_FORMAT_STRING"))

#' msa.write
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"))
msa.write(m, "foo.ss", format="SS")
msa.write(m, "foo.fasta", format="FASTA", pretty.print=TRUE)
msa.write(m, NULL, format="PHYLIP", pretty.print=TRUE)
#clean up
unlink("foo.ss")
unlink("foo.fasta")


#' msa.new
# Here is an MSA object stored in the default mode
m1 <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
              names=c("human", "mouse", "rat"))
m2 <- m1
# NOTE seqs would not be directly accessible if stored by reference
m2$seqs[3] <- "AAAAAA"
print(m1)
print(m1, printSeq=TRUE)
print(m2, printSeq=TRUE)


#' subMsa
m1 <- msa.new(seqs=c("ACGT---AT", "AGGTAGTAA", "AGGAAGTAG"),
              names=c("human", "mouse", "rat"), pointer.only=TRUE)
print(subMsa(m1, c("human", "rat"), start.col=3, end.col=6),
      print.seq=TRUE)
print(subMsa(m1, c("mouse"), keep=FALSE, refseq="human", start.col=3, end.col=4),
      print.seq=TRUE)

#' msa.names
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"))
msa.names(m)
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"),
             pointer.only=TRUE)
msa.names(m)
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"))
msa.names(m)

#' read.msa
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, c("ENr334.maf", "gencode.ENr334.gff"))
m <- read.msa("ENr334.maf", format="MAF")
m
msa.seqlen(m) # total number of columns
msa.seqlen(m, c("hg18", "canFam2"))  # individual sequence lengths in alignment
#'
# this doesn't work yet!
m <- read.msa("ENr334.maf", format="MAF", gff.only=read.gff("gencode.ENr334.gff"))
m
msa.seqlen(m, c("hg18", "canFam2"))
unlink(c("ENr334.maf", "gencode.ENr334.gff")) #clean up

