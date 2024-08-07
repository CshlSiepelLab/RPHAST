% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feat.R
\name{fix.start.stop.feat}
\alias{fix.start.stop.feat}
\title{Fix start and stop signals}
\usage{
fix.start.stop.feat(x)
}
\arguments{
\item{x}{An object of type \code{feat}.  CDS regions must be present with type
"CDS", and the transcript_id must be indicated in the attribute field.
Start and stop codons should have feature type "start_codon" and "stop_codon"
(as produced by addSignals.feat).}
}
\value{
An object of type \code{feat}, in which CDS regions are ensured to
include start codons and exclude stop codons, as required by the GTF2 standard.
}
\description{
Fix start and stop signals
}
\note{
\itemize{
\item{If x is stored as a pointer to an object stored in C, signals will be
added to x.}
\item{Assumes at most one start_codon and at most one stop_codon per transcript.}
}
}
\examples{
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
f <- add.signals.feat(f)
# let's just look at one gene
geneNames <- tagval.feat(f, "transcript_id")
f <- f[geneNames==geneNames[1],]

# This features file already is correct, so let's mess it up to see
# how fix.start.stop can fix it:

#modify first CDS to not include start
startCodon <- f[f$feature=="start_codon",]
firstCds <- which(f$feature=="CDS" & f$start==startCodon$start)
f[firstCds,]$start <- startCodon$end+1
#modify last CDS to include stop
stopCodon <- f[f$feature=="stop_codon",]
lastCds <- which(f$feature=="CDS" & f$end+1==stopCodon$start)
f[lastCds,]$end <- stopCodon$end
# now call fix.start.stop to fix
f.fixed <- fix.start.stop.feat(f)

# first CDS has been fixed to include start codon
f[firstCds,]
f.fixed[firstCds,]
# last CDS has been fixed to not include stop codon
f[lastCds,]
f.fixed[lastCds,]

unlink(featFile)
}
\author{
Melissa J. Hubisz and Adam Siepel
}
\keyword{features}
