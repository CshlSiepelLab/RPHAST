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

#' tm.summary
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

#' isSubstMod.tm
isSubstMod.tm(c("JC69", "K80", "F81", "HKY85", "HKY85+Gap",
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
