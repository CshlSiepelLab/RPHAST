library(rphast)

#' read.gff
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
gffFile <- "gencode.ENr334.gff"
unzip(exampleArchive, gffFile)
g <- read.gff(gffFile)
dim(g)
g[1:10,]
unlink(gffFile)

# non-example tests
unzip(exampleArchive, gffFile)
g <- read.gff(gffFile, pointer.only=TRUE)
dim.gff(g)
unlink(gffFile)

################################################

#' plot.gff
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
gffFile <- "gencode.ENr334.gff"
unzip(exampleArchive, gffFile)
g <- read.gff(gffFile)
# note that plot(g) does not work because GFFs are stored as data.frames
plot.gff(g[g$feature=="CDS",])
unlink(gffFile)


################################################

#' gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff(seq, src, feature, start, end)
dim(g)
dim.gff(g)
g <- gff(seq, src, feature, start, end, pointer.only=TRUE)
dim.gff(g)
       

################################################

#' as.pointer.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g1 <- gff(seq, src, feature, start, end)
g2 <- as.pointer.gff(g1)
g1
g2


# non-example test
g <- as.pointer.gff(g2)
summary(g)

################################################

#' write.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff(seq, src, feature, start, end)
write.gff("test.gff", g)
#'
unlink("test.gff") # clean up

# non-example tests
g <- gff(seq, src, feature, start, end)
write.gff("test.gff", g)
unlink("test.gff")

################################################

#' nrow.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff(seq, src, feature, start, end)
nrow.gff(g)

# non-example tests
nrow.gff(as.pointer.gff(g))

################################################

#' ncol.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff(seq, src, feature, start, end)
ncol.gff(g)
ncol.gff(as.pointer.gff(g))


################################################

#' summary.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff(seq, src, feature, start, end)
summary(g)  # this calls summary.data.frame
summary(as.pointer.gff(g))


################################################

#' as.data.frame.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g1 <- gff(seq, src, feature, start, end, pointer.only=TRUE)
summary(g1)
g2 <- as.data.frame(g1)
summary(g2)
dim(g2)


################################################

#' dim.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g1 <- gff(seq, src, feature, start, end)
dim(g1)
dim.gff(g1)
g2 <- as.pointer.gff(g1)
dim(g2)
dim.gff(g2)

###############################################

#' overlap.gff
require("rphast")
gff1 <- gff(seqname=c(rep("chr1", 3), rep("chr2", 2)),
            src="example",
            feature=".",
            start=c(1, 5, 100, 10, 20),
            end=c(7, 10, 105, 15, 30))
gff2 <- gff(seqname=c("chr1","chr2"), src="test", feature=".", start=c(1,1), end=c(5,10))
#'
overlap.gff(gff1, gff1)
overlap.gff(gff1, gff2, min.percent=0.25)
overlap.gff(gff1, gff2, min.percent=0.25, overlapping=FALSE)
overlap.gff(gff1, gff2, get.fragments=TRUE)
overlap.gff(gff1, gff2, get.fragments=TRUE)
rm(gff1, gff2)


#' coverage.gff
q("no")
require("rphast")
gff1 <- gff(seqname=c(rep("chr1", 3), rep("chr2", 2)),
            src="example",
            feature=".",
            start=c(1, 5, 100, 10, 20),
            end=c(7, 10, 105, 15, 30))
gff2 <- gff(seqname=c("chr1","chr2"), src="test", feature=".", start=c(1,1), end=c(5,10))
coverage.gff(gff1, gff2, or=FALSE)
coverage.gff(gff1, gff2, or=TRUE)
coverage.gff(gff1, gff2, get.feats=TRUE, or=TRUE)
coverage.gff(gff1, gff2, or=TRUE)

q("no")
require("rphast")
gff1 <- read.gff("/home/melissa/rphast/RPHAST/inst/extdata/neutral-loci.bed")
coverage.gff(gff1)
gff2 <- read.gff("/home/melissa/rphast/RPHAST/inst/extdata/kb1.homoDif.bed")
coverage.gff(gff2)
coverage.gff(gff1, gff2)
coverage.gff(gff1, gff2, or=TRUE)
gff3 <- read.gff("/home/melissa/rphast/RPHAST/inst/extdata/watson.homoDif.bed")
coverage.gff(gff1, gff2, gff3)
coverage.gff(gff1, gff2, gff3, or=TRUE)
