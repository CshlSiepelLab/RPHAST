exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gff"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
# note that plot(f) does not work because features are stored as data.frames
plot.feat(f[f$feature=="CDS",])
unlink(featFile)
