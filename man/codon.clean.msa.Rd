% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msa.R
\name{codon.clean.msa}
\alias{codon.clean.msa}
\title{Clean an alignment for codon analysis}
\usage{
codon.clean.msa(x, refseq = NULL, strand = "+")
}
\arguments{
\item{x}{An object of type \code{msa}}

\item{refseq}{The name of the reference sequence to be used.  If given,
strip all columns which contain gaps in refseq.  Once this is done,
alignment should be in frame.  If \code{refseq==NULL} then alignment
should be in frame as it is sent in (no gaps are stripped).}

\item{strand}{Either "+" or "-".  If "-", reverse complement the
alignment.}
}
\value{
An object of type \code{msa}.  It will be the same as the
original msa, with the following modifications:
\itemize{
\item If refseq is not NULL, columns with gaps in refseq will be stripped.
\item If strand is "-", the new msa will be the reverse complement of
the original.
\item After the gap stripping and reverse complementing steps, each
sequence is searched for stop codons.  If encountered, the stop codon
and the rest of the sequence to follow is converted to missing data.  The
resulting msa has a length equal to the longest remaining sequence (end
columns with all missing data are removed).
}
}
\description{
Clean an alignment for codon analysis
}
\note{
If the input msa (x) is stored as a pointer, its value will be
changed to the return value of the function.
}
\author{
Melissa J. Hubisz
}
