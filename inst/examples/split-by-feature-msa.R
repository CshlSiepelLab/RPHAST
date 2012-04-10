require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, c("ENr334.maf", "gencode.ENr334.gff"))      
m <- read.msa("ENr334.maf")
feats <- read.feat("gencode.ENr334.gff")
feats$seqname <- "hg18"
cdsAlign <- split.by.feature.msa(m, feats[feats$feature=="CDS",])
