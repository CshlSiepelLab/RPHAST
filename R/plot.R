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


##' Get the coordinate range of a list of RPHAST results
##' @param resultList a list of RPHAST results
##' @param resultType a vector of character strings, each should be either
##' "wig" or "gff"
##' @return a numeric vector of length two giving the minimum and maximum
##' coordinates in any element of the list
##' @export
range.rphast <- function(resultList, resultType) {
  xmin <- min(mapply(function(el,ty) {
    if (nrow(el)==0) return(NA)
    if (ty=="wig") return(min(el$coord))
    return(min(el$start))},
                     resultList, resultType), na.rm=TRUE)
  xmax <- max(mapply(function(el, ty) {
    if (nrow(el)==0) return(NA)
    if (ty=="wig") return(max(el$coord))
    return(max(el$end))},
                     resultList, resultType), na.rm=TRUE)
  c(xmin, xmax)
}


##' Plot RPHAST results
##' @param resultList a list of plottable features or RPHAST results.  Possible elemets are
##' any object of type \code{gff}, lists with two elements: coord and a score (as output by phyloP)
##' @param doLabels Logical.  Whether to plot the label above each plot.  Will be
##' recycled to the length of resultList.  Does not affect printing of shortLabels.
##' @param labels Labels to appear directly above each plot.
##' @param cex.labels The character expansion factor for the labels
##' @param shortLabels If not NULL, a character vector for each plot to print in the left margin.
##' @param cex.shortLabels The character expansion factor for the shortLabels
##' @param col Colors to use for each plot.  Will be truncated/recycled to achieve the length of resultList.
##' @param doStrand Logical.  Whether to plot strand info for \code{GFF} objects.  Can be single value or
##' vector of \code{length(resultList)}.  If it is a vector, elements that do not correspond to feature
##' plots will be ignored.
##' @param relWigSize The relative size of wig plots compared to feature plots
##' @param xlim The range of the x coordinate to be plotted.  If \code{NULL} (the default), will
##' use the entire range represented in the resultList.
##' @param xlab The label for the x axis
##' @param  ylab The label for the y axis
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
plot.rphast <- function(resultList,
                        doLabels=TRUE,
                        labels=names(resultList),
                        cex.labels=0.75,
                        shortLabels=NULL,
                        cex.shortLabels=0.5,
                        col=c("red", "green", "blue"),
                        doStrand=TRUE,
                        relWigSize=5,
                        xlim=NULL,
                        xlab="coord", ylab="", blankSpace=0.25, axisDigits=3,
                        labelSpace=min(length(resultList)*0.05, 0.25),
                        belowLabelSpace=0.2, lmar=4, ...) {
  numresult <- length(resultList)

  doLabels <- rep(doLabels, length.out=numresult)
  check.arg(doLabels, "doLabels", "logical", null.OK=FALSE,
            min.length=numresult, max.length=numresult)
  numlabels <- sum(doLabels)
  if (numlabels > 0L) 
    check.arg(labels, "labels", "character", null.OK=FALSE,
              min.length=numresult, max.length=numresult)
  
  check.arg(col, "col", null.OK=FALSE, min.length=NULL, max.length=NULL)
  col <- rep(col, length.out=numresult)

  check.arg(relWigSize, "relWigSize", "numeric", null.OK=FALSE)
  check.arg(xlim, "xlim", "numeric", null.OK=TRUE, min.length=2L, max.length=2L)
  check.arg(xlab, "xlab", "character", null.OK=FALSE)
  check.arg(ylab, "ylab", "character", null.OK=FALSE)
  check.arg(blankSpace, "blankSpace", "numeric", null.OK=FALSE)
  if (blankSpace < 0 || blankSpace > 1) stop("blankSpace should be between 0 and 1")
  check.arg(axisDigits, "axisDigits", "integer", null.OK=FALSE)

  resultType <- guessResultType.rphast(resultList)
  if (sum(resultType=="gff") > 0L) {
    doStrand <- rep(doStrand, length.out=numresult)
    check.arg(doStrand, "doStrand", "logical", null.OK=FALSE, min.length=NULL, max.length=NULL)
  }
  plotScale <- as.numeric(ifelse(resultType=="wig", relWigSize, 1))
  plotScale <- plotScale/sum(plotScale)

  if (is.null(xlim))
    xlim <- range.rphast(resultList, resultType)

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
    el <- resultList[[i]]
    if (doLabels[i]) {
      miny <- maxy - labelSize
      text(x=mean(xlim), y=miny, labels[i], pos=3, offset=belowLabelSpace, cex=cex.labels[i])
      maxy <- miny
    }
    miny <- maxy - plotScale[i]
    yrange <- c(miny, maxy)

    if (nrow(el) > 0) {
    if (resultType[[i]] == "wig") {
      coord <- el$coord
      othercols <- which(names(el)!="coord")
      if (length(othercols)!=1L)
        warning("element ", names(resultList)[i],
                " has several statistics, using last item ", names(el)[othercols[length(othercols)]])
      stat <- el[[othercols[length(othercols)]]]
      oldrange <- range(stat[coord >= xlim[1] & coord <= xlim[2]])
      newstat <- (stat - oldrange[1])*(maxy - yrange[1])/(oldrange[2]-oldrange[1]) + yrange[1]
      lines(el$coord, newstat, col=col[i])
      axis(side=2, at=yrange, labels=FALSE)
      mtext(format(oldrange[1], digits=axisDigits), side=2, line=0.5, at=yrange[1], las=1, cex=0.75)
      mtext(format(oldrange[2], digits=axisDigits), side=2, line=0.5, at=yrange[2], las=1, cex=0.75)
    } else {
      plot.gff(el, doStrand=doStrand[i], y=mean(yrange), width=(yrange[2]-yrange[1]),
               add=TRUE, col=col[i])
    }
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
