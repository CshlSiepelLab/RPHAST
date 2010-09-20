require("rphast")

#' feat.track
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "sol1.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
featTrack <- feat.track(f, "basic feature track")
f <- addIntrons.feat(f)
geneTrack <- gene.track(f, "gene track")
plot.track(list(featTrack, geneTrack))
plot.track(list(featTrack, geneTrack, geneTrack, geneTrack, geneTrack), xlim=c(14800, 16000))

#' plot.gene
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "sol1.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
plot.gene(f)
plot.gene(f, xlim=c(0, 10000))  #zoom in
