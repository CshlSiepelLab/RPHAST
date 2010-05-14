# (may not want to export in the long run but good to have visible for testing)
##' @export
rphast.simplify.list <- function(lol) {
  if (!is.list(lol)) return(lol)
  if (length(lol) == 1) 
    return(rphast.simplify.list(lol[[1]]))
  currClass <- attr(lol, "class")
  attr(lol, "class") <- NULL
  # this is a little ugly but I can't think of a better way to deal with special
  # conversion issues.  No way to assign "NA" in C so if frames are undefined
  # they are -1, set to NA here.
  isGff <- (!is.null(currClass) && currClass=="gff")
  if (isGff) {
    if (!is.null(lol$frame)) {
      lol$frame[lol$frame < 0] <- NA
    }
    isDataFrame <- TRUE;
    isMatrix <- FALSE
  } else {
    isMatrix <- (!is.null(currClass) && currClass=="matrix")
    isDataFrame <- isGff || (!is.null(currClass) && currClass=="data.frame")
  }
  if (isMatrix || isDataFrame) {
    if (!is.null(lol$row.names)) {
      rowNames <- lol$row.names
      lol$row.names <- NULL
    } else rowNames <- NULL
    if (isMatrix) lol <- as.matrix(as.data.frame(lol), check.names=FALSE)
    if (isDataFrame) lol <- as.data.frame(lol, check.names=FALSE)
    if (!is.null(rowNames)) row.names(lol) <- rowNames
  } else {
    for (i in 1:length(lol))
      lol[[i]] <- rphast.simplify.list(lol[[i]])
  }
  lol
}
