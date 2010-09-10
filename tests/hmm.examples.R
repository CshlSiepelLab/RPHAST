library(rphast)

#' hmm
h <- hmm(matrix(1, nrow=4, ncol=4))
h


#' score.hmm
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334.maf", "rev.mod", "gencode.ENr334.gff")
unzip(exampleArchive, files)
# make "conserved" and "neutral" models and a phylo-HMM that describes
# transitions between them, and predict conserved elements (this
# is the same thing that phastCons does, but can be extended to general
# phylo-HMMs)
align <- read.msa("ENr334.maf")
neutralMod <- read.tm("rev.mod")
#'
# create a conserved model
conservedMod <- neutralMod
conservedMod$tree <- rescale.tree(neutralMod$tree, 0.3)
#'
# create a simple phylo-HMM
state.names <- c("neutral", "conserved")
h <- hmm(matrix(c(0.99, 0.01, 0.01, 0.99), nrow=2, dimnames=list(state.names, state.names)),
                eq.freq=c(neutral=0.9, conserved=0.1))
scores <- score.hmm(align, mod=list(neutral=neutralMod, conserved=conservedMod),
                    hmm=h, states="conserved")
# try an alternate approach of comparing likelihoods of genes 
feats <- read.feat("gencode.ENr334.gff")
# plot in a region with some genes
plot.track(list(feat.track(scores$in.states, name="hmmScores"),
                feat.track(feats[feats$feature=="CDS",], name="genes")),
           xlim=c(41650000, 41680000))
unlink(files)

#' read.hmm
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
file <- "coding.hmm"
unzip(exampleArchive, file)
# this is a 5-state hmm with states representing
# intergenic, intron, first, second, and third codon positions.
h <- read.hmm(file)
h
unlink(file)


#' write.hmm
state.names <- c("neutral", "conserved")
h <- hmm(matrix(c(0.99, 0.01, 0.01, 0.99), nrow=2,
                dimnames=list(state.names, state.names)),
                eq.freq=c(neutral=0.9, conserved=0.1))
filename <- tempfile()
write.hmm(h, filename)
unlink(filename)

gc()
