#' phastCons
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334.fa", "rev.mod")
unzip(exampleArchive, files)
mod <- read.tm("rev.mod")
msa <- read.msa("ENr334.fa")
rv <- phastCons(msa, mod)
names(rv)
rv2 <- phastCons(msa, mod, estimate.trees=TRUE)
names(rv2)
rv2$tree.models
unlink(files)
