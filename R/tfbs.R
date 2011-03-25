##' Reads Multiple Sequences from a file.
##' @title Reading an MS Object
##' @param filename The name of the input file containing the sequences.
##' @param alphabet the alphabet of non-missing-data chraracters in the
##' sequences.  Determined automatically from the alignment if not given.
##' @param pointer.only If \code{TRUE}, MS will be stored by reference as
##' an external pointer to an object created by C code, rather than
##' directly in R memory.  This improves performance and may be necessary
##' for large alignments, but reduces functionality.  See
##' \code{\link{ms}} for more details on MS object storage options.
##' @return an MS object.
##' @keywords ms
##' @seealso \code{\link{ms}}
##' @export
read.ms <- function(filename, alphabet=NULL, pointer.only=FALSE) {

  cats.cycle <- NULL
  check.arg(filename, "filename", "character", null.OK=FALSE)
  check.arg(alphabet, "alphabet", "character", null.OK=TRUE)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  x <- list()
  x$externalPtr <- .Call.rphast("rph_ms_read", filename, alphabet)
 # if (!pointer.only) x$seqs <- .Call.rphast("rph_ms_seqs", src$externalPtr)
  x
}


##' Reads names of Multiple Sequences from a file.
##' @title Reading names of an MS Object
##' @param msListP List of Multiple Sequence objects 
##' @return List of names
##' @keywords ms, names
##' @seealso \code{\link{ms}}
##' @export
names.ms <- function(msListP) {
  if (!is.null(msListP$externalPtr)) {
    return(.Call.rphast("rph_ms_seqNames", msListP=msListP$externalPtr))
  }
  msListP$names
}


##' Clip sequences to supplied windows and group by GC content.
##' @title Clipping sequences to supplied windows and grouping by GC content
##' @param[in,out] msP Multiple sequences
##' @param[in] windowsP Set of windows 
##' @param[in] ngroupsP Number of GC content groups
##' @return List of windows, each window containing a Multiple Sequences object
##' @export
split.ms <- function(msListP, windowsP) {
  if (!is.null(msListP$externalPtr)) {
  ##' @note check windowsP pointer
    x <- list()
    x$externalPtr <-.Call.rphast("rph_ms_clip_to_supplied_windows", msListP=msListP$externalPtr, windowsP=as.pointer.feat(windowsP)$externalPtr)
    x
  }
}


##' Group by GC content
##' @param[in]
##' @param[in] ngroupsP Number of groups of GC content
##' @return Sequences in Groups of GC content
##' @export
group.ms <- function(msListP, ngroupsP) {
   if (!is.null(msListP$externalPtr)) {
  ##' @note check windowsP pointer
    x <- list()
    x$externalPtr <-.Call.rphast("rph_ms_group_by_gc_content", msListP=msListP$externalPtr, ngroupsP=ngroupsP)
    x
  }
}

##' Build Markov Model
##' @param[in] windowsP Single group of sequences, split by window
##' @param[in] norderP Order of N-order Markov Model
##' @return Markov Matrix
##' @export
build.mm <- function(msListP, norderP) {
  if(!is.null(msListP$externalPtr)) {
     x <- list()
     x$externalPtr <- .Call.rphast("rph_build_mm", msListP=msListP$externalPtr, norderP=norderP)
     x$order = norderP
     x$alph_size <- .Call.rphast("rph_alph_size", msListP=msListP$externalPtr);
     x
  }
}

##' Read PWM from a MEME formated file
##' @param[in] filenameP Filename of the MEME formated file containing PWM as motif
##' @return Position Weight Matrix from file
##' @export
read.pwm <- function(filenameP) {
  x$externalPtr <- .Call.rphast("rph_pwm_read", filenameP)
  x
}

##' Compute scores for a give PWM 
##' @todo fill in params
##' @return scores
##' @export
compute.scores <- function(msGroupListP, pwmP, mmListP) {
  
  x <- .Call.rphast("rph_compute_scores", msGroupListP=msGroupListP$externalPtr, pwmP$externalPtr[[1]], mmListP=mmListP$externalPtr, mmListP$order);
  x
}

##' Generate simulated sequence based on Markov Model
##' @todo fill in params
##' @return simulated sequence
##' @export
generate.seq <- function(mmP, lengthP) {
  x$externalPtr <- .Call.rphast("rph_mm_simulate_seq", mmP=mmP$externalPtr, mmP$order, mmP$alph_size, lengthP);
  x
}


#calcFDR <- function(groupScoresReal, simulatedScores, threshold)
#{
#  numAboveThreshReal <- 0;
#  flatScoresReal <- unlist(groupScoresReal);
#  lenRealSeqs <- length(flatScoresReal);
#  
#  for(baseNum in 1:lenRealSeqs)
#  {
#    if(flatScoresReal[baseNum] >= threshold)
#      numAboveThreshReal = numAboveThreshReal +1;
#  }
#
#  numAboveThreshSim <- 0;
#  flatScoresSim <- unlist(simulatedScores);
#  lenSimSeqs <- length(flatScoresSim);
#  
#  for(baseNum in 1:lenSimSeqs)
#  {
#    if(flatScoresSim[baseNum] >= threshold)
#      numAboveThreshSim <- numAboveThreshSim +1;
#  }
#  
#  #fdrVal <- ((numAboveThreshReal / lenRealSeqs) * (lenSimSeqs / numAboveThreshSim));
#  fdrVal <-  ((numAboveThreshSim / lenSimSeqs) * (lenRealSeqs / numAboveThreshReal));
#  fdrVal
#} 

##' FDR helper function
##' @todo fill in params
##' @export
calcFDR <- function(groupScoresReal, simulatedScores, realScoreIndex)
{
  if (length(simulatedScores[[1]]) > 0)
  {
    numAboveThreshReal <- 0;
    threshold <- groupScoresReal[[realScoreIndex]][[1]];
    lenRealSeqs <- groupScoresReal[[realScoreIndex]][[5]];
    
    for(scoreNum in 1:length(groupScoresReal))
    {
      if(groupScoresReal[[scoreNum]][[1]] >= threshold)
        numAboveThreshReal = numAboveThreshReal + 1;
    }
    
    numAboveThreshSim <- simulatedScores[[1]][[1]];
    lenSimSeqs <- simulatedScores[[1]][[5]];
    
    fdrVal <- ((numAboveThreshSim / lenSimSeqs) * (lenRealSeqs / numAboveThreshReal));
    fdrVal;
  } else {
    0; #If we don't have at least score for the simulated sequence then the FDR function will eval to 0
  }
}

##' Plot FDR scores 
##' @todo fill in params
##' @export
plot.fdr <- function(realSeqsScores, simSeqsScores) {
   xToPlot <- list();
   yToPlot <- list();
   numGroups <- length(realSeqsScores);
   #For each group
   for(groupNum in 1:length(realSeqsScores))
   {
     xaxis <- list();
     yaxis <- list();
 
     realGroup = realSeqsScores[[groupNum]];
     simGroup = simSeqsScores[[groupNum]];

     for(scoreNum in 1:length(realGroup))
     {
       if (length(realGroup[[scoreNum]]) > 0)
       { 
         xaxis <- c(xaxis, realGroup[[scoreNum]][[1]]); #Add score to xaxis list
         yaxis <- c(yaxis, calcFDR(realGroup, simGroup, scoreNum));
       }
     } 
     xToPlot[[groupNum]] <- xaxis;
     yToPlot[[groupNum]] <- yaxis;
     
   }
   cc <- unlist(xToPlot);
   xmin <- Inf; xmax <- -Inf; for (i in 1:length(cc)) { if (abs(cc[[i]]) != Inf) { if (cc[i] > xmax) xmax <- cc[i];  if (cc[i] < xmin) xmin <- cc[i] } }
   cc <- unlist(yToPlot);
   ymin <- Inf; ymax <- -Inf; for (i in 1:length(cc)) { if (abs(cc[[i]]) != Inf) { if (cc[i] > ymax) ymax <- cc[i];  if (cc[i] < ymin) ymin <- cc[i] } }
   
  # get the range for the x and y axis
  xrange <- list(xmin, xmax)
  yrange <- list(ymin, ymax)

  # set up the plot
  plot(xrange, yrange, type="n", xlab="score",
     ylab="FDR" )
  colors <- rainbow(numGroups)
  linetype <- c(1:numGroups)
  plotchar <- seq(18,18+numGroups,1)

  # add lines
  for (i in 1:numGroups) {
    lines(xToPlot[[i]], yToPlot[[i]], type="b", lwd=2,
      lty=linetype[i], col=colors[i], pch=plotchar[i])  
  }

  # add a title and subtitle
  title("False Discovery Rate", "")

  # add a legend
  legend(xrange[1], yrange[2], 1:numGroups, cex=0.8, col=colors,
     pch=plotchar, lty=linetype, title="GC Content Group")
     
}


