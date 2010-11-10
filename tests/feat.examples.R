library(rphast)

#' read.feat
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gff"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
dim(f)
f[1:10,]
unlink(featFile)

# non-example tests
unzip(exampleArchive, featFile)
f <- read.feat(featFile, pointer.only=TRUE)
dim.feat(f)
unlink(featFile)

################################################

#' plot.feat
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gff"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
# note that plot(f) does not work because features are stored as data.frames
plot.feat(f[f$feature=="CDS",])
unlink(featFile)


################################################

#' feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f <- feat(seq, src, feature, start, end)
dim(f)
dim.feat(f)
f <- feat(seq, src, feature, start, end, pointer.only=TRUE)
dim.feat(f)
       

################################################

#' as.pointer.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f1 <- feat(seq, src, feature, start, end)
f2 <- as.pointer.feat(f1)
f1
f2


# non-example test
f <- as.pointer.feat(f2)
summary(f)

################################################

#' write.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f <- feat(seq, src, feature, start, end)
write.feat(f, "test.gff")
#'
unlink("test.gff") # clean up

# non-example tests
f <- feat(seq, src, feature, start, end)
write.feat(f, "test.gff")
unlink("test.gff")

################################################

#' nrow.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f <- feat(seq, src, feature, start, end)
nrow.feat(f)

# non-example tests
nrow.feat(as.pointer.feat(f))

################################################

#' ncol.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f <- feat(seq, src, feature, start, end)
ncol.feat(f)
ncol.feat(as.pointer.feat(f))


################################################

#' summary.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f <- feat(seq, src, feature, start, end)
summary(f)  # this calls summary.data.frame
summary(as.pointer.feat(f))


################################################

#' as.data.frame.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f1 <- feat(seq, src, feature, start, end, pointer.only=TRUE)
summary(f1)
f2 <- as.data.frame(f1)
summary(f2)
dim(f2)


################################################

#' dim.feat
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f1 <- feat(seq, src, feature, start, end)
dim(f1)
dim.feat(f1)
f2 <- as.pointer.feat(f1)
dim(f2)
dim.feat(f2)

###############################################

#' overlap.feat
feat1 <- feat(seqname=c(rep("chr1", 3), rep("chr2", 2)),
              start=c(1, 5, 100, 10, 20),
              end=c(7, 10, 105, 15, 30))
feat2 <- feat(seqname=c("chr1","chr2"),
              start=c(1,1),
              end=c(5,10))
#'
overlap.feat(feat1, feat1)
overlap.feat(feat1, feat2, min.percent=0.25)
overlap.feat(feat1, feat2, min.percent=0.25, overlapping=FALSE)
overlap.feat(feat1, feat2, get.fragments=TRUE)
overlap.feat(feat1, feat2, get.fragments=TRUE)
rm(feat1, feat2)


#' coverage.feat
require("rphast")
feat1 <- feat(seqname=c(rep("chr1", 3), rep("chr2", 2)),
              start=c(1, 5, 100, 10, 20),
              end=c(7, 10, 105, 15, 30))
feat2 <- feat(seqname=c("chr1","chr2"),
              start=c(1,1), end=c(5,10))
coverage.feat(feat1, feat2, or=FALSE)
coverage.feat(feat1, feat2, or=TRUE)
coverage.feat(feat1, feat2, get.feats=TRUE, or=TRUE)
coverage.feat(feat1, feat2, or=TRUE)
rm(feat1, feat2)

#' tagval
require("rphast")
tags <- c("tag1 \"val 1a\"; tag2 \"val 2a\" \"val2a.1\" 123; tag3 \"val3a\"",
          "tag1 \"val 1b\"; tag2 \"val 2b\"; tag4 \"val4b\"",
          "tag3 \"val3a\" 1; tag4 2;")
tagval(tags, "tag1")
tagval(tags, "tag2")
tagval(tags, "tag3")
tagval(tags, "tag4")
tagval(tags, "notag")
rm(tags)

#' tagval.feat
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "sol1.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
geneName <- tagval.feat(f, "transcript_id")
geneName[1:10]
length(unique(geneName)) # number of unique genes
unlink(featFile)
rm(f, geneName)


#' add.UTRs.feat
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
table(f$feature)
coverage.feat(f[f$feature=="CDS",])
coverage.feat(f[f$feature=="exon",])
f <- add.UTRs.feat(f)
table(f$feature)
coverage.feat(f[f$feature=="3'UTR",])
coverage.feat(f[f$feature=="5'UTR",])

#' add.introns.feat
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
table(f$feature)
coverage.feat(f[f$feature=="CDS",])
coverage.feat(f[f$feature=="exon",])
f <- add.introns.feat(f)
table(f$feature)
coverage.feat(f[f$feature=="intron",])


#' add.signals.feat
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
table(f$feature)
coverage.feat(f[f$feature=="CDS",])
coverage.feat(f[f$feature=="exon",])
f <- add.signals.feat(f)
table(f$feature)
coverage.feat(f[f$feature=="3'splice",])
coverage.feat(f[f$feature=="5'splice",])
coverage.feat(f[f$feature=="start_codon",])
coverage.feat(f[f$feature=="stop_codon",])


#' fix.start.stop.feat
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
f <- add.signals.feat(f)
# let's just look at one gene
geneNames <- tagval.feat(f, "transcript_id")
f <- f[geneNames==geneNames[1],]
#'
# This features file already is correct, so let's mess it up to see
# how fix.start.stop can fix it:
#'
#modify first CDS to not include start
startCodon <- f[f$feature=="start_codon",]
firstCds <- which(f$feature=="CDS" & f$start==startCodon$start)
f[firstCds,]$start <- startCodon$end+1
#modify last CDS to include stop
stopCodon <- f[f$feature=="stop_codon",]
lastCds <- which(f$feature=="CDS" & f$end+1==stopCodon$start)
f[lastCds,]$end <- stopCodon$end
# now call fix.start.stop to fix
f.fixed <- fix.start.stop.feat(f)
#'
# first CDS has been fixed to include start codon
f[firstCds,]
f.fixed[firstCds,]
# last CDS has been fixed to not include stop codon
f[lastCds,]
f.fixed[lastCds,]

rm(list = ls())
gc()
