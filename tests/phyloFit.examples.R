require("rphast")

#' phyloFit
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334.maf", "ENr334.fa", "gencode.ENr334.gff", "rev.mod")
unzip(exampleArchive, files)
m <- read.msa("ENr334.maf")
tm <- phyloFit(m, tree="((hg18, (mm9, rn4)), canFam2)")
tm
phyloFit(m, init.mod=tm)
likelihood.msa(m, tm)
tm$likelihood
print(tm$likelihood, digits=10)
f <- read.feat("gencode.ENr334.gff")
t <- phyloFit(m, tree="((hg18, (mm9, rn4)), canFam2)",
              features=f, quiet=TRUE)
names(t)
t$other
t[["5'flank"]]
phyloFit(m, init.mod=t$AR, nrates=3, alpha=4.0)
phyloFit(m, init.mod=t$AR, rate.constants=c(10, 5, 1))
