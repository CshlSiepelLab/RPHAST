require("rphast")

#' print.tm
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
tm
print(tm, aslist=TRUE)
unlink(filename)

#' read.tm
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
tm
unlink(filename)

#' write.tm
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
tm
write.tm(tm, NULL)
write.tm(tm, "test.mod")
unlink(c(filename, "test.mod"))

#' summary.tm
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
read.tm(filename)
unlink(c(filename, "test.mod"))


#' as.list.tm
tm <- tm(tree="((human:0.01, chimp:0.01):0.03, mouse:0.3)",
         subst.mod="JC69")
is.list(tm)
is.list(as.list(tm))

#' subst.mods
subst.mods()

#' is.subst.mod.tm
is.subst.mod.tm(c("JC69", "K80", "F81", "HKY85", "HKY85+Gap",
                  "REV", "SSREV", "REV+GC", "UNREST", "R2", "U2", "R2S",
                  "U2S", "R3", "R3S", "U3", "U3S", "GC", "HB",
                  "bad.model"))


#' tm
tree <- "((human:0.01, chimp:0.01):0.03, mouse:0.3)"
subst.mod <- "JC69"
rate.mat <- matrix(runif(16), nrow=4, ncol=4)
for (i in 1:4)
  rate.mat[i,i] <- -sum(rate.mat[i,-i])
backgd <- runif(4)
backgd <- backgd/sum(backgd)
alphabet <- "ACGT"
t <- tm(tree, subst.mod, rate.mat, backgd, alphabet)
t
#'
nratecats <- 3
alpha <- 1.5
rate.consts <- runif(nratecats, max=3.0)
root.leaf <- "human"
t <- tm(tree, subst.mod, rate.matrix=rate.mat,
        backgd=backgd, alphabet=alphabet,
        nratecats=nratecats, alpha=alpha,
        rate.consts=rate.consts, root.leaf=root.leaf)
t


#' plot.tm
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
plot(tm)
plot(tm, show.eq.freq=FALSE)
plot(tm, max.cex=20, eq.freq.max.cex=1, col=matrix(1:16, nrow=4), eq.freq.col=c("red", "green"), filled=TRUE, add=TRUE)
plot.rate.matrix(tm[["rate.matrix"]], eq.freq=tm[["backgd"]], filled=FALSE)
plot.rate.matrix(tm[["rate.matrix"]], eq.freq=tm[["backgd"]], filled=TRUE, add=TRUE)
unlink(filename)


#' plot.altmodel.tm
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
tm <- add.alt.mod(tm, branch="mm9", subst.mod="HKY85")
plot.altmodel.tm(tm, 1)
tm$alt.model$backgd <- c(0.9, 0.05, 0.03, 0.02)
plot.altmodel.tm(tm, 1)
plot.rate.matrix(tm[["rate.matrix"]], eq.freq=tm[["backgd"]], filled=FALSE, alphabet=tm[["alphabet"]])
unlink(filename)


#' mod.backgd.tm
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
filename <- "rev.mod"
unzip(exampleArchive, filename)
tm <- read.tm(filename)
#'
# change background frequencies to new value, adjusting rate matrix
mod.backgd.tm(tm, c(0.25, 0.25, 0.25, 0.25))
#'
# change background frequencies so that GC content is 0.6
mod.backgd.tm(tm, gc=0.6)

rm(list=ls())
gc()

