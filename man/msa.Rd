\name{msa.new}
\alias{msa.new}
\alias{msa}
\title{MSA ojbects}
\description{
  Creates a new MSA (multiple sequence alignment) object.
}
\usage{
msa.new(seqs, names=NULL, alphabet="ACGT", pointerOnly=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
    a character vector of length nseqs, each string of length seqlen,
    containing the alignment.
    }
  \item{names}{
    a character vector of length nseqs with the name of each sample,
    in the same order as given in the alignment.
}
  \item{alphabet}{
    a character string containing all valid data characters, not
    including missing data (assumed to be "N") or gaps (assumed to be
    "-").  The alphabet is not case-sensitive.
}
\item{pointerOnly}{
  a logical value indicating whether to store alignment in R, or only
  as a pointer to an object in C.  See details below for important
  warnings about this option.  Setting this to true saves memory and
  improves performance, but for a loss of functionality.  See also
  important warning about pointer objects in RPHAST.
}}
\details{
  Creates a new alignment object with given sequences.

  Setting pointerOnly to TRUE saves memory and improves performance, but
  then the msa object can only be viewed or modified through calls to
  RPHAST functions.  Furthermore, pointer objects will not always behave as
  expected.  See warning below.
}
\value{
  Returns a new msa object.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\section{Warning }{ Unlike most R objects, msa objects do not "copy by
  value".  Due to the potentially large size of alignments, the
  alignment is not copied unless msa.copy() is called explicitly.  An
  assignment msa2 <- msa1 will result in both MSA objects pointing to
  the same alignment in memory.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{msa.copy}}, \code{\link{print.msa}}, \code{\link{write.msa}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

m <- msa.new(nseqs=3, seqlen=6, names=c("human", "mouse", "rat"))
m
print(m, printSeq=TRUE)

m <- msa.new(names=c("human", "mouse", "rat"),
             seqs=c("ACGTAT", "AGGTAA", "AGGTAG"))
m
print(m, printSeq=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
