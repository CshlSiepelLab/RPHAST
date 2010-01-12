q("no")
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

#' gff.new
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff.new(seq, src, feature, start, end)
dim(g)
dim.gff(g)
g <- gff.new(seq, src, feature, start, end, pointer.only=TRUE)
dim.gff(g)
       

################################################

#' gff.to.pointer
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g1 <- gff.new(seq, src, feature, start, end)
g2 <- gff.to.pointer(g1)
g1
g2


# non-example test
g <- gff.to.pointer(g2)
summary(g)

################################################

#' write.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff.new(seq, src, feature, start, end)
write.gff("test.gff", g)
#'
unlink("test.gff") # clean up

# non-example tests
g <- gff.new(seq, src, feature, start, end)
write.gff("test.gff", g)
unlink("test.gff")

################################################

#' gff.numrow
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff.new(seq, src, feature, start, end)
gff.numrow(g)

# non-example tests
gff.numrow(gff.to.pointer(g))

################################################

#' gff.numcol
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff.new(seq, src, feature, start, end)
gff.numcol(g)
gff.numcol(gff.to.pointer(g))


################################################

#' summary.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g <- gff.new(seq, src, feature, start, end)
summary(g)  # this calls summary.data.frame
summary(gff.to.pointer(g))


################################################

#' as.data.frame.gff
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
g1 <- gff.new(seq, src, feature, start, end, pointer.only=TRUE)
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
g1 <- gff.new(seq, src, feature, start, end)
dim(g1)
dim.gff(g1)
g2 <- gff.to.pointer(g1)
dim(g2)
dim.gff(g2)

