#basic idea: take any of the following:
# phyloP results
# phastCons Results
# gff

# and plot all in one plot, vertically stacking each set of results with a
# single x axis and local y-axis for wig-type plots.  Should be an argument
# which optionally scales relative sizes of each plot.  Should also be able
# to plot generic wig with "coord" and "score" or generic bed with
# "start" "end"

# question: what do we do about results, like phyloP results, that have
# a coord column and several statistics (from base-by-base type output)
# always plot the last statistic?


##' Guess the contents of a list containing RPHAST results
##' @param l a list containing plottable results, such as wig or gff
##' @return a character vector of plot types, the same length as the input
##' list, each element is either "wig" for lists that have a coord element,
##' "gff" for lists that have start and end elements, or "unknown" for others
##' @export
guessResultType.rphast <- function(l) {
  lapply(l, function(x) {
    if (!is.null(x$coord)) {
      return("wig")
    } else if ((!is.null(x$start)) && (!is.null(x$end))) {
      return("gff")
    }
    "unknown"
  }
         )
}


##' Is this a track?
##' @param x An object to test
##' @param ... ignored
##' @return A logical indicating whether x is an object of type track
##' @export
is.track <- function(x, ...) {
  if (is.null(attr(x, "class"))) return(FALSE)
  attr(x, "class") == "track"
}


##' Get the coordinate range of a list of RPHAST results
##' @param ... a list of tracks
##' @param resultType a vector of character strings, each should be either
##' "wig" or "gff"
##' @return a numeric vector of length two giving the minimum and maximum
##' coordinates in any element of the list
##' @export
range.track <- function(..., na.rm=TRUE) {
  l <- list(...)
  if (length(l)==1 && !is.track(l[[1]])) 
    l <- l[[1]]
  currrange <- c()
  for (i in 1:length(l)) {
    if (!is.track(l[[i]])) stop("element ", i, " is not a track")
    if (l[[i]]$type == "wig") {
      currrange <- c(currrange, range(l[[i]]$data$coord))
    } else {
      currrange <- c(currrange, range.gff(l[[i]]$data))
    }
  }
  range(currrange)
}



##' Smooth a wig plot in rphast
##' @param coord The x coordinates of un-smoothed plot
##' @param score The scores cooresponding to the x coordinates (should be same length as coord)
##' @param numpoints The number of points to use in the new plot
##' @return A data frame with numpoints rows and columns "coord" and "score"
##' with smoothed values.  If \code{length(coord) <= numpoints}, it will contain the
##' original data
##' @export
smooth.wig <- function(coord, score, numpoints=300) {
  if (length(coord) != length(score)) stop("smooth.wig expects length(coord) == length(score)")
  if (length(coord) <= numpoints) return(data.frame(coord=coord, score=score))

  xlim <- range(coord)
  windowSize <- (xlim[2]-xlim[1])/(numpoints-1)
  newcoord <- seq(from=xlim[1], to=xlim[2], by=windowSize)

  newscore <- sapply(newcoord, function(x) {
    f <- (coord >= (x-windowSize) &
          coord <= (x+windowSize))
    if (sum(f)==0) return(NA)
    sum(score[f])/sum(f)
  })
  f <- !is.na(newscore)
  data.frame(coord=newcoord[f], score=newscore[f])
}


##' Make browser-like plot in rphast
##' @param x a list of tracks, created by wigTrack or gffTrack.
##' @param doLabels Logical.  Whether to plot the label above each plot.  Will be
##' recycled to the length of x.  Does not affect printing of shortLabels.
##' @param labels Labels to appear directly above each plot.
##' @param cex.labels The character expansion factor for the labels
##' @param shortLabels If not NULL, a character vector for each plot to print in the left margin.
##' @param cex.shortLabels The character expansion factor for the shortLabels
##' @param relWigSize The relative size of wig plots compared to feature plots
##' @param xlim The range of the x coordinate to be plotted.  If \code{NULL} (the default), will
##' use the entire range represented in the resultList.
##' @param xlab The label for the x axis
##' @param ylab The label for the y axis
##' @param blankSpace The amount of vertical blank space between each plot.  This should be a single numeric
##' value between 0 and 1, representing the total fraction of the plot occupied by blank space.
##' @param axisDigits The number of digits to use on the y-axis for wig plots.
##' @param labelSpace The total fraction of vertical space given to plot labels.
##' @param belowLabelSpace The amount of space between a label and the plot it corresponds to,
##' in fractions of a character width.
##' @param lmar The size of the left margin (in number of lines)
##' @param ... Other options to be passed to \code{plot}.  See \link{par}.
##' @seealso \code{plotPhast}, which may be easier to use but less flexible
##' @export
plot.track <- function(x,
                       doLabels=TRUE,
                       cex.labels=0.75,
                       shortLabels=NULL,
                       cex.shortLabels=0.5,
                       relWigSize=5,
                       xlim=NULL,
                       xlab="coord", ylab="", blankSpace=0.25, axisDigits=3,
                       labelSpace=min(length(x)*0.05, 0.25),
                       belowLabelSpace=0.2, lmar=4, ...) {
  tracks <- x
  if (length(tracks) == 1 && is.track(tracks)) tracks <- list(tracks)
  numresult <- length(tracks)

  doLabels <- rep(doLabels, length.out=numresult)
  check.arg(doLabels, "doLabels", "logical", null.OK=FALSE,
            min.length=numresult, max.length=numresult)
  numlabels <- sum(doLabels)
  
  check.arg(relWigSize, "relWigSize", "numeric", null.OK=FALSE)
  check.arg(xlim, "xlim", "numeric", null.OK=TRUE, min.length=2L, max.length=2L)
  check.arg(xlab, "xlab", "character", null.OK=FALSE)
  check.arg(ylab, "ylab", "character", null.OK=FALSE)
  check.arg(blankSpace, "blankSpace", "numeric", null.OK=FALSE)
  if (blankSpace < 0 || blankSpace > 1) stop("blankSpace should be between 0 and 1")
  check.arg(axisDigits, "axisDigits", "integer", null.OK=FALSE)

  resultType <- character(length(tracks))
  for (i in 1:length(tracks)) resultType[i] <- tracks[[i]]$type
  plotScale <- as.numeric(ifelse(resultType=="wig", relWigSize, 1))
  plotScale <- plotScale/sum(plotScale)

  if (is.null(xlim))
    xlim <- range.track(tracks)

  plot(c(0), c(0), type="n", xlim=xlim, ylim=c(0, 1), xlab=xlab,
       ylab=ylab, yaxt="n", bty="n", ...)
  
  maxy <- 1
  if (numlabels > 0L) {
    check.arg(belowLabelSpace, "belowLabelSpace", "numeric", null.OK=FALSE)
    check.arg(labelSpace, "labelSpace", "numeric", null.OK=FALSE)
    cex.labels <- rep(cex.labels, length.out=numresult)
    check.arg(cex.labels, "cex.labels", "numeric", null.OK=FALSE, min.length=NULL, max.length=NULL)
    if (labelSpace < 0 || labelSpace > 1) stop("labelSpace should be between 0 and 1")
    labelSize <- sum(plotScale)*labelSpace/numlabels
    plotScale <- plotScale*(1.0-labelSpace)
  }
  check.arg(shortLabels, "shortLabels", "character", null.OK=TRUE, 
            min.length=numresult, max.length=numresult)
  if (!is.null(shortLabels)) {
    cex.shortLabels <- rep(cex.shortLabels, length.out=numresult)
    check.arg(cex.shortLabels, "cex.shortLabels", "numeric", null.OK=FALSE,
              min.length=NULL, max.length=NULL)
  }
  plotScale <- plotScale*(1-blankSpace)
  blankSpace <- blankSpace/numresult
  
  for (i in 1:numresult) {
    el <- tracks[[i]]
    if (doLabels[i]) {
      miny <- maxy - labelSize
      text(x=mean(xlim), y=miny, el$name, pos=3, offset=belowLabelSpace, cex=cex.labels[i])
      maxy <- miny
    }
    miny <- maxy - plotScale[i]
    yrange <- c(miny, maxy)

    if (resultType[[i]] == "wig") {
      coord <- el$data$coord
      f <- (coord >= xlim[1] & coord <= xlim[2])
      coord <- coord[f]
      score <- el$data$score[f]
      if (el$smooth) {
        if (is.null(el$numpoints)) smoothData <- smooth.wig(coord, score)
        else smoothData <- smooth.wig(coord, score, el$numpoints)
        coord <- smoothData$coord
        score <- smoothData$score
      }
      if (is.null(el$ylim))
        oldrange <- range(score)
      else oldrange <- el$ylim
      newscore <- (score - oldrange[1])*(maxy - yrange[1])/(oldrange[2]-oldrange[1]) + yrange[1]
      lines(coord, newscore, col=el$col)
      axis(side=2, at=yrange, labels=FALSE)
      mtext(format(oldrange[1], digits=axisDigits), side=2, line=0.5, at=yrange[1], las=1, cex=0.75)
      mtext(format(oldrange[2], digits=axisDigits), side=2, line=0.5, at=yrange[2], las=1, cex=0.75)
      if (!is.null(el$horiz.line)) {
        horiz.line <- (el$horiz.line - oldrange[1])*(maxy - yrange[1])/(oldrange[2]-oldrange[1]) + yrange[1]
        for (i in 1:length(el$horiz.line)) 
          abline(h=horiz.line, lty=el$horiz.lty[i], col=el$horiz.col[i])
      }
    } else {
      plot.gff(el$data, doStrand=doStrand[i], y=mean(yrange), width=(yrange[2]-yrange[1]),
               add=TRUE, col=el$col)
    }
    if (!is.null(shortLabels))
      mtext(shortLabels[i], side=2, line=0.5, at=mean(yrange), las=1, cex=cex.shortLabels[i])
    maxy <- miny-blankSpace
  }

}



# recursively "flatten" a list of lists, keeping only elements which are
# wig-type (have coord/score) or are gff-type.
##' Get a list suitable for sending to plot.rphast
##' @param resultList a possibly recursive list of \code{GFF} objects, phastCons results, and/or
##' phyloP results.  Can also contain generic data frames or lists.  If one column has name "coord",
##' will use this column as the x-axis and the last column whose name isn't coord as the y-axis.
##' If columns exist named "start" and "end", will plot as a \code{GFF}.  Other elements of the
##' recursive list are ignored.
##' @param prefix A character string to append to the names of all elements.  Usually \code{NULL}
##' on the initial call (but used for the recursion).
##' @export
getPlottableResults <- function(resultList, prefix=NULL) {
  rv.list <- list()
  nextIdx <- 1
  for (i in 1:length(resultList)) {
    curr <- resultList[[i]]
    if (is.null(names(resultList))) {
      currName <- as.character(i)
    } else currName <- names(resultList)[i]
    if (!is.null(prefix)) 
      currName <- paste(prefix, currName, sep=" ")
    if (is.list(curr)) {
      if (((!is.null(curr$coord)) && length(curr) >= 2L) ||
          ((!is.null(curr$start)) && (!is.null(curr$end)))) {
        rv.list[[nextIdx]] <- curr
        names(rv.list)[nextIdx] <- currName
        nextIdx <- nextIdx + 1
      } else {
        temp <- getPlottableResults(curr,prefix=currName)
        if (length(temp) >= 1L) {
          for (j in 1:length(temp)) {
            rv.list[[nextIdx]] <- temp[[j]]
            names(rv.list)[nextIdx] <- names(temp)[j]
            nextIdx <- nextIdx + 1
          }
        }
      }
    }
  }
  rv.list
}

##' Create a wig track
##' @param coord A numeric vector of coordinates (to be used for x-axis)
##' @param score A numeric vector of scores (y-axis coords), should be same length as coord
##' @param name The name of the track (a character string)
##' @param col The color to be used to plot this track.
##' @param ylim The limits to be used on the y-axis.  If NULL use entire range
##' of score.
##' @param smooth A logical value indicating whether to perform smoothing when plotting
##' this track
##' @param numpoints (Only used if \code{smooth==TRUE}).  An integer value indicating how many
##' points to display in the smoothed wig.
##' @param horiz.line If non-NULL, draw horizontal lines on the display at the given y coordinates
##' @param horiz.lty If horiz.line is defined, use this line type.
##' @param horiz.col If horiz.line is defined, use this color
##' @return An object of type \code{track} which can be plotted with the plot.track
##' function
##' @export
wigTrack <- function(coord, score, name, col="black", ylim=NULL, smooth=FALSE, numpoints=250,
                     horiz.line=NULL, horiz.lty=2, horiz.col="black") {
  rv <- list()
  attr(rv, "class") <- "track"
  rv$data <- data.frame(coord=coord, score=score)
  rv$name <- name
  rv$col <- col
  if (!is.null(ylim))
    rv$ylim <- ylim
  rv$smooth <- smooth
  if (rv$smooth)
    rv$numpoints <- numpoints
  if (!is.null(horiz.line)) {
    rv$horiz.line <- horiz.line
    rv$horiz.lty <- rep(horiz.lty, length.out=length(horiz.line))
    rv$horiz.col <- rep(horiz.col, length.out=length(horiz.line))
  }
  rv$type="wig"
  rv
}


##' Create a features track
##' @param gff An object of type \code{gff}
##' @param name The name of the track (a character string)
##' @param col The color to use plotting this track (can be a single
##' color or a color for each element)
##' @return An object of type \code{track} which can be plotted with plot.track
##' function
##' @export
gffTrack <- function(gff, name, col="black") {
  rv <- list()
  attr(rv, "class") <- "track"
  rv$data <- gff
  rv$name <- name
  rv$col <- col
  rv$type="gff"
  rv
}


