require("rphast")

#' phyloFit
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334.maf", "ENr334.fa", "gencode.ENr334.gff", "rev.mod")
unzip(exampleArchive, files)
m <- read.msa("ENr334.maf")
mod <- phyloFit(m, tree="((hg18, (mm9, rn4)), canFam2)")
mod
phyloFit(m, init.mod=mod)
likelihood.msa(m, mod)
mod$likelihood
print(mod$likelihood, digits=10)
f <- read.feat("gencode.ENr334.gff")
mod <- phyloFit(m, tree="((hg18, (mm9, rn4)), canFam2)",
                features=f, quiet=TRUE)
names(mod)
mod$other
mod[["5'flank"]]
phyloFit(m, init.mod=mod$AR, nrates=3, alpha=4.0)
phyloFit(m, init.mod=mod$AR, rate.constants=c(10, 5, 1))
