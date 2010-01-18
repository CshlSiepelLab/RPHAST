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

###############################################

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

################################################

#' msa.validFormatStr
msa.validFormatStr(c("MAF", "SS", "PHYLIP", "MPM", "LAV", "FASTA",
                     "BAD_FORMAT_STRING"))

################################################

#' msa.format.from.ext
msa.format.from.ext("file/asdf.fa")
msa.format.from.ext("blah\blah.maf")
msa.format.from.ext("abcd/efg/hij.ss")
msa.format.from.ext("maf.fasta.lav")
msa.format.from.ext("helloWorld")

################################################

#' msa.seqlen
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"))
msa.seqlen(m)
msa.seqlen(m, msa.names(m))

## non-example tests
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
msa.seqlen(m)
msa.seqlen(m, msa.names(m))
print(m, print.seq=TRUE)

################################################

#' msa.nseq
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"))
msa.nseq(m)

## non-example tests
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
msa.nseq(m)


################################################

#' msa.offset
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"))
msa.offset(m)
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), offset=500000)
msa.offset(m)

## non-example tests
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
msa.offset(m)
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"),
             offset=500000, pointer.only=TRUE)
msa.offset(m)


################################################

#' msa.alphabet
m <- msa.new(seqs=c("a--acgtaa", "NN-nnnTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"))
msa.alphabet(m)

## non-example tests
m <- msa.new(seqs=c("a--acgtaa", "NN-nnnTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
msa.alphabet(m)

################################################

#' msa.is.ordered
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"))
msa.is.ordered(m)
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), is.ordered=FALSE)
msa.is.ordered(m)

# non-example tests
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
msa.is.ordered(m)
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), is.ordered=FALSE,
             pointer.only=TRUE)
msa.is.ordered(m)

################################################


#' msa.from.pointer
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
m
m <- msa.from.pointer(m)
m

# non-example tests
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"))
m
m <- msa.from.pointer(m)
m

################################################

#' msa.to.pointer
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"))
m
m <- msa.to.pointer(m)
m

# non-example tests
m <- msa.new(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
m
m <- msa.to.pointer(m)
m


################################################

#' msa.write
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"))
msa.write(m, "foo.ss")
msa.write(m, "foo.fa", pretty.print=TRUE)
msa.write(m, NULL, format="PHYLIP", pretty.print=TRUE)
#'
#clean up
unlink("foo.ss")
unlink("foo.fasta")

################################################

#' msa.new
# Here is an MSA object stored in the default mode
m1 <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
              names=c("human", "mouse", "rat"))
m2 <- m1
# NOTE seqs would not be directly accessible if stored by reference
m2$seqs[3] <- "AAAAAA"
print(m1)
print(m1, print.seq=TRUE)
print(m2, print.seq=TRUE)

################################################

#' msa.sub
m <- msa.new(seqs=c("ACGT---AT", "AGGTAGTAA", "AGGAAGTAG"),
              names=c("human", "mouse", "rat"))
print(msa.sub(m, c("human", "rat"), start.col=3, end.col=6),
      print.seq=TRUE)
print(msa.sub(m, c("mouse"), keep=FALSE, refseq="human", start.col=3, end.col=4),
      print.seq=TRUE)

# non-example tests
print(msa.sub(m, c("human", "rat")), print.seq=TRUE)
print(msa.sub(m, c("human", "rat"), end.col=6), print.seq=TRUE)
print(msa.sub(m, c("human", "rat"), start.col=3), print.seq=TRUE)
              
m <- msa.new(seqs=c("ACGT---AT", "AGGTAGTAA", "AGGAAGTAG"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
print(msa.sub(m, c("human", "rat"), start.col=3, end.col=6),
      print.seq=TRUE)
print(msa.sub(m, c("mouse"), keep=FALSE, refseq="human", start.col=3, end.col=4),
      print.seq=TRUE)

################################################

#' msa.names
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"))
msa.names(m)
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"))
msa.names(m)

# non-example tests
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
             names=c("human", "mouse", "rat"),
             pointer.only=TRUE)
msa.names(m)
m <- msa.new(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"), pointer.only=TRUE)
msa.names(m)


################################################

#' msa.strip.gaps
m <- msa.new(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
             names=c("human", "mouse", "rat"))
print(msa.strip.gaps(m, c("human", "mouse")), print.seq=TRUE)
print(msa.strip.gaps(m, strip.mode="any.gaps"), print.seq=TRUE)
print(msa.strip.gaps(m, strip.mode="all.gaps"), print.seq=TRUE)
print(m, print.seq=TRUE)
#' NOTE if msa stored as pointer, original object is changed
m <- msa.to.pointer(m)
temp <- msa.strip.gaps(m, "any.gaps")
print(m, print.seq=TRUE)

# non-example tests
m <- msa.new(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
             names=c("human", "mouse", "rat"), pointer.only=TRUE)
print(msa.strip.gaps(m, c("human", "mouse")), print.seq=TRUE)
print(msa.strip.gaps(m, strip.mode="any.gaps"), print.seq=TRUE)
print(msa.strip.gaps(m, strip.mode="all.gaps"), print.seq=TRUE)
print(m, print.seq=TRUE)
################################################

#' read.msa
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334.maf", "ENr334.fa", "gencode.ENr334.gff")
unzip(exampleArchive, files)
#'
# Read a fasta file, ENr334.fa
# this file represents a 4-way alignment of the encode region
# ENr334 starting from hg18 chr6 position 41405894
idx.offset <- 41405894
m1 <- read.msa("ENr334.fa", offset=idx.offset)
m1
#'
# Now read in only a subset represented in a feature file
g <- read.gff("gencode.ENr334.gff")
g$seqname <- "hg18"  # need to tweak source name to match name in alignment
m1 <- read.msa("ENr334.fa", gff=g, offset=idx.offset)
#'
# Can also subset on certain features
do.cats <- c("CDS", "5'flank", "3'flank")
m1 <- read.msa("ENr334.fa", gff=g, offset=idx.offset,
              do.cats=do.cats)
#'
# Can read MAFs similarly, but don't need offset because
# MAF file is annotated with coordinates
m2 <- read.msa("ENr334.maf", gff=g, do.cats=do.cats)
# Also, note that when a GFF is given and the file is
# in MAF format, the first sequence is automatically
# stripped of gaps
msa.seqlen(m1)
msa.seqlen(m2)
msa.seqlen(m1, "hg18")
#'
unlink(files) # clean up

# non-example tests
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, files)
#'
# Read a fasta file, ENr334.fa
# this file represents a 4-way alignment of the encode region
# ENr334 starting from hg18 chr6 position 41405894
idx.offset <- 41405894
m1 <- read.msa("ENr334.fa", offset=idx.offset, pointer.only=TRUE)
m1
#'
# Now read in only a subset represented in a feature file
g <- read.gff("gencode.ENr334.gff")
g$seqname <- "hg18"  # need to tweak source name to match name in alignment
m1 <- read.msa("ENr334.fa", gff=g, offset=idx.offset, pointer.only=TRUE)
#'
# Can also subset on certain features
do.cats <- c("CDS", "5'flank", "3'flank")
m1 <- read.msa("ENr334.fa", gff=g, offset=idx.offset,
              do.cats=do.cats, pointer.only=TRUE)
#'
# Can read MAFs similarly, but don't need offset because
# MAF file is annotated with coordinates
m2 <- read.msa("ENr334.maf", gff=g, do.cats=do.cats, pointer.only=TRUE)
# Also, note that when a GFF is given and the file is
# in MAF format, the first sequence is automatically
# stripped of gaps
msa.seqlen(m1)
msa.seqlen(m2)
msa.seqlen(m1, "hg18")
#'
unlink(files) # clean up

################################################

#' msa.likelhood
files <- c("rev.mod", "ENr334.maf", "ENr334.fa")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, files)
msa <- read.msa("ENr334.fa")
tm <- read.tm("rev.mod")
l <- msa.likelihood(msa, tm)
msa.likelihood(msa, tm)
like1 <- msa.likelihood(msa, tm, by.column=TRUE)
like2 <- msa.likelihood(msa, tm, by.column=TRUE)
msa <- read.msa("ENr334.maf")
msa.likelihood(msa, tm)
msa.likelihood(msa, tm)
like3 <- msa.likelihood(msa, tm, by.column=TRUE)
like4 <- msa.likelihood(msa, tm, by.column=TRUE)
tm$subst.mod <- "JC69"
msa.likelihood(msa, tm)
