require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, c("ENr334-100k.maf", "gencode.ENr334-100k.gff"))      
m <- read.msa("ENr334-100k.maf")
feats <- read.feat("gencode.ENr334-100k.gff")
feats$seqname <- "hg18"
cdsAlign <- split.by.feature.msa(m, feats[feats$feature=="CDS",])
