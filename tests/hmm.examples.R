library(rphast)

#'hmm
m <- matrix(1, nrow=4,ncol=4)
h <- hmm(m)
h2 <- as.pointer.hmm(h)
rm(h2)
gc()
