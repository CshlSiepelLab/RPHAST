% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feat.R
\name{feat}
\alias{feat}
\title{Features Objects}
\usage{
feat(
  seqname = "default",
  src = ".",
  feature = ".",
  start,
  end,
  score = NULL,
  strand = NULL,
  frame = NULL,
  attribute = NULL,
  pointer.only = FALSE
)
}
\arguments{
\item{seqname}{a character vector containing the name of the sequence.  If
the features correspond to regions of an alignment, then seqname
should be the name of the sequence in the alignment that is used
as the frame of reference in the features.  To use the entire alignment
as a frame of reference, set seqname to "MSA".}

\item{src}{The source of the feature}

\item{feature}{The feature type name}

\item{start}{The start of the feature.  Sequence numbering begins at 1.}

\item{end}{The end of the feature.  This is the last coordinate included
in the feature.}

\item{score}{The feature score, or NA if there is no score.}

\item{strand}{A character string which is either "+", "-", or "." (if
strand is not available or relevant).}

\item{frame}{A 0, 1, or 2, which specifies whether the feature is in frame.}

\item{attribute}{A feature attribute (character string).}

\item{pointer.only}{Whether to store object as a pointer to an object
in C, rather than as a data.frame in R.}
}
\value{
If pointer.only==FALSE, returns a data.frame whose format
mirrors the GFF specification.  Otherwise, returns a list with a single
object, which is a external pointer to a C structure representing a
features object.
}
\description{
Create a new features object
}
\details{
See \url{http://www.sanger.ac.uk/resources/software/gff/spec.html}
for more detailed description of each parameter.

All arguments which are provided should be vectors of equal length.

If pointer.only==FALSE, the new object is a data frame, with
columns mirroring the GFF Specification
Otherwise, it is a list containing a single element, which is
a pointer to an object stored in C.
}
\examples{
seq <- rep("hg18.chr6", 10)
src <- rep("fake_example", 10)
feature <- rep("CDS", 10)
start <- seq(1, 100, by=10)
end <- seq(10, 100, by=10)
f <- feat(seq, src, feature, start, end)
dim(f)
dim.feat(f)
f <- feat(seq, src, feature, start, end, pointer.only=TRUE)
dim.feat(f)
}
\seealso{
\code{\link{read.feat}}

\code{\link{msa}} for more details on the pointer.only option.
}
\author{
Melissa J. Hubisz and Adam Siepel
}
\keyword{features}
