library(rphast)

#' print.msa
# read in an MSA stored in R
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
print(m)
print(m, format="FASTA")
print(m, format="PHYLIP", pretty.print=TRUE)
#'
# read in an MSA stored by reference in C
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"),
         pointer.only=TRUE)
print(m)

###############################################

#' summary.msa
# read in an MSA stored in R
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
summary(m)
#'
# read in an MSA stored by reference in C
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"),
         pointer.only=TRUE)
summary(m)

################################################

#' is.format.msa
is.format.msa(c("MAF", "SS", "PHYLIP", "MPM", "LAV", "FASTA",
                "BAD_FORMAT_STRING"))

################################################

#' guess.format.msa
guess.format.msa("file/asdf.fa")
guess.format.msa("blah\blah.maf")
guess.format.msa("abcd/efg/hij.ss")
guess.format.msa("maf.fasta.lav")
guess.format.msa("helloWorld")

################################################

#' ncol.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
ncol.msa(m)
ncol.msa(m, names.msa(m))

## non-example tests
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
ncol.msa(m)
ncol.msa(m, names.msa(m))
print(m, print.seq=TRUE)

################################################

#' nrow.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
nrow.msa(m)

## non-example tests
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
nrow.msa(m)


################################################

#' offset.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
offset.msa(m)
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), offset=500000)
offset.msa(m)

## non-example tests
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
offset.msa(m)
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"),
         offset=500000, pointer.only=TRUE)
offset.msa(m)


################################################

#' alphabet.msa
m <- msa(seqs=c("a--acgtaa", "NN-nnnTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
alphabet.msa(m)

## non-example tests
m <- msa(seqs=c("a--acgtaa", "NN-nnnTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
alphabet.msa(m)

################################################

#' is.ordered.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
is.ordered.msa(m)
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), is.ordered=FALSE)
is.ordered.msa(m)

# non-example tests
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
is.ordered.msa(m)
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), is.ordered=FALSE,
         pointer.only=TRUE)
is.ordered.msa(m)

################################################


#' from.pointer.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
m
m <- from.pointer.msa(m)
m

# non-example tests
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
m
m <- from.pointer.msa(m)
m

################################################

#' as.pointer.msa
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
m
m <- as.pointer.msa(m)
m

# non-example tests
m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
m
m <- as.pointer.msa(m)
m


################################################

#' write.msa
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
write.msa(m, "foo.ss")
write.msa(m, "foo.fa", pretty.print=TRUE)
write.msa(m, NULL, format="PHYLIP", pretty.print=TRUE)
#'
#clean up
unlink("foo.ss")
unlink("foo.fa")

################################################

#' msa
# Here is an MSA object stored in the default mode
m1 <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
          names=c("human", "mouse", "rat"))
m2 <- m1
# NOTE seqs would not be directly accessible if stored by reference
m2$seqs[3] <- "AAAAAA"
print(m1)
print(m1, print.seq=TRUE)
print(m2, print.seq=TRUE)

################################################

#' sub.msa
m <- msa(seqs=c("ACGT---AT", "AGGTAGTAA", "AGGAAGTAG"),
         names=c("human", "mouse", "rat"))
print(sub.msa(m, c("human", "rat"), start.col=3, end.col=6),
      print.seq=TRUE)
print(sub.msa(m, c("mouse"), keep=FALSE, refseq="human",
              start.col=3, end.col=4),
      print.seq=TRUE)

# non-example tests
print(sub.msa(m, c("human", "rat")), print.seq=TRUE)
print(sub.msa(m, c("human", "rat"), end.col=6), print.seq=TRUE)
print(sub.msa(m, c("human", "rat"), start.col=3), print.seq=TRUE)
              
m <- msa(seqs=c("ACGT---AT", "AGGTAGTAA", "AGGAAGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
print(sub.msa(m, c("human", "rat"), start.col=3, end.col=6),
      print.seq=TRUE)
print(sub.msa(m, c("mouse"), keep=FALSE, refseq="human",
              start.col=3, end.col=4),
      print.seq=TRUE)
################################################

#' [.msa
require("rphast")
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
print(m[c("rat", "rat", "human"), ], print.seq=TRUE)
print(m[c(3,3,1),], print.seq=TRUE)
print(m[c(TRUE, FALSE, TRUE),], print.seq=TRUE)
print(m[TRUE,], print.seq=TRUE)
print("[.msa"(m, "mouse",c(1,6,3,5)), print.seq=TRUE)


################################################

#' sample.msa
require("rphast")
m <- msa(seqs=c("AAAAAAAAAACCCCCGGT", "GGGGGGGGGGTTTTTCCA", "CCCCCCCCCCAAAAAGGA"),
         names=c("human", "mouse", "rat"))
sample.msa(m, 10, replace=TRUE)
sample.msa(m, 10, replace=TRUE, prob=c(rep(1, 10), rep(2, 5), rep(5, 2), 10))

################################################

#' names.msa
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
names.msa(m)
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"))
names.msa(m)

# non-example tests
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"),
         pointer.only=TRUE)
names.msa(m)
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"), pointer.only=TRUE)
names.msa(m)


################################################

#' strip.gaps.msa
m <- msa(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
         names=c("human", "mouse", "rat"))
print(strip.gaps.msa(m, c("human", "mouse")), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="any.gaps"), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="all.gaps"), print.seq=TRUE)
print(m, print.seq=TRUE)
#' NOTE if msa stored as pointer, original object is changed
m <- as.pointer.msa(m)
temp <- strip.gaps.msa(m, "any.gaps")
print(m, print.seq=TRUE)

# non-example tests
m <- msa(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
print(strip.gaps.msa(m, c("human", "mouse")), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="any.gaps"), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="all.gaps"), print.seq=TRUE)
print(m, print.seq=TRUE)
################################################

#' read.msa
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
f <- read.feat("gencode.ENr334.gff")
f$seqname <- "hg18"  # need to tweak source name to match name in alignment
m1 <- read.msa("ENr334.fa", features=f, offset=idx.offset)
#'
# Can also subset on certain features
do.cats <- c("CDS", "5'flank", "3'flank")
m1 <- read.msa("ENr334.fa", features=f, offset=idx.offset,
               do.cats=do.cats)
#'
# Can read MAFs similarly, but don't need offset because
# MAF file is annotated with coordinates
m2 <- read.msa("ENr334.maf", features=f, do.cats=do.cats)
# Also, note that when features is given and the file is
# in MAF format, the first sequence is automatically
# stripped of gaps
ncol.msa(m1)
ncol.msa(m2)
ncol.msa(m1, "hg18")
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
f <- read.feat("gencode.ENr334.gff")
f$seqname <- "hg18"  # need to tweak source name to match name in alignment
m1 <- read.msa("ENr334.fa", features=f, offset=idx.offset, pointer.only=TRUE)
#'
# Can also subset on certain features
do.cats <- c("CDS", "5'flank", "3'flank")
m1 <- read.msa("ENr334.fa", features=f, offset=idx.offset,
               do.cats=do.cats, pointer.only=TRUE)
#'
# Can read MAFs similarly, but don't need offset because
# MAF file is annotated with coordinates
m2 <- read.msa("ENr334.maf", features=f, do.cats=do.cats, pointer.only=TRUE)
# Also, note that when features is given and the file is
# in MAF format, the first sequence is automatically
# stripped of gaps
ncol.msa(m1)
ncol.msa(m2)
ncol.msa(m1, "hg18")
#'
unlink(files) # clean up

################################################

#' likelihood.msa
require("rphast")
files <- c("rev.mod", "ENr334.maf", "ENr334.fa", "small.gff")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, files)
msa <- read.msa("ENr334.fa")
mod <- read.tm("rev.mod")
likelihood.msa(msa, mod)
like1 <- likelihood.msa(msa, mod, by.column=TRUE)
length(like1)==ncol.msa(msa)
sum(like1)
msa <- read.msa("ENr334.maf")
likelihood.msa(msa, mod)
like2 <- likelihood.msa(msa, mod, by.column=TRUE)
sum(like2)
mod$subst.mod <- "JC69"
likelihood.msa(msa, mod)
#'
# can also get likelihood by feature
features <- read.feat("small.gff")
features$seqname <- names(msa)[1]
likelihood.msa(msa, mod, features=features)
unlink(files)


####################################################

#' simulate.msa
filename <- "rev.mod"
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, filename)
m <- matrix(nrow=3, ncol=3)
m[1,] <- c(1,2,3)
m[2,] <- c(1,5,10)
m[3,] <- c(10,4,2)
eq.freq <- c(1,2,3)
h <- hmm(m, eq.freq)
mod <- read.tm(filename)
mod2 <- mod
mod2$backgd <- rep(0.25, 4)
mod3 <- mod
mod3$backgd <- c(0.6, 0.1, 0.2, 0.1)
m <- simulate.msa(mod, 20)
m <- simulate.msa(list(mod, mod2, mod3), 20, hmm=h)
m <- matrix(1, nrow=9, ncol=9)
h <- hmm(m)
m <- simulate.msa(list(mod, mod2, mod3, mod, mod2, mod3, mod, mod2, mod3),
                  20, hmm=h)
unlink(filename)

#' get4d.msa
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334.maf", "ENr334.fa", "gencode.ENr334.gff")
unzip(exampleArchive, files)
f <- read.feat("gencode.ENr334.gff")
f$seqname <- "hg18.chr6"
m1 <- read.msa("ENr334.maf", features=f, do.4d=TRUE)
m2 <- read.msa("ENr334.maf")
m3 <- get4d.msa(m2, features=f)
m4 <- get4d.msa(read.msa("ENr334.maf"), features=f)
m5 <- get4d.msa(read.msa("ENr334.fa", offset=41405894), features=f)
unlink(files)

#' informative.regions.msa
require("rphast")
m <- msa(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
         names=c("human", "mouse", "rat"))
informative.regions.msa(m, 1, refseq=NULL)
informative.regions.msa(m, 3, refseq=NULL)
informative.regions.msa(m, 3, refseq="mouse", spec=c("mouse", "rat"))


#' postprobs.msa
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, "ENr334.maf")
m <- read.msa("ENr334.maf")
mod <- phyloFit(m, tree="((human,(mm9,rn4)),canFam2)")
x <- postprob.msa(sub.msa(m, start.col=41447839, end.col=41448033, refseq="hg18"), mod)
dim(x)
dimnames(x)
x[,,"CCCC"]

#' expected.subs.msa
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, "ENr334.maf")
m <- read.msa("ENr334.maf")
mod <- phyloFit(m, tree="((human,(mm9,rn4)),canFam2)")
x <- expected.subs.msa(sub.msa(m, start.col=41447839, end.col=41448033, refseq="hg18"), mod)
dim(x)
dimnames(x)
x[,"CCCC"]
x["mm9-rn4",]

#' total.expected.subs.msa
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, "ENr334.maf")
m <- read.msa("ENr334.maf")
mod <- phyloFit(m, tree="((human,(mm9,rn4)),canFam2)")
x <- total.expected.subs.msa(sub.msa(m, start.col=41447839, end.col=41448033, refseq="hg18"), mod)
dim(x)
dimnames(x)
x["mm9-rn4",,]


#' plot.msa
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, "ENr334.maf")
m <- read.msa("ENr334.maf")
plot.msa(m)
plot.msa(m[, 1:2])
plot.msa(m[,1:20])
plot.msa(m[1:3,1:40])
plot.msa(m[,1:100])
plot.msa(m[,1:50], refseq=NULL)
rm(m)
unlink("ENr334.maf")

rm(list = ls())
gc()
