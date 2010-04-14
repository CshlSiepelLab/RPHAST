##' @export
rphast.simplify.list <- function(lol) {
  if (!is.list(lol)) return(lol)
  if (length(lol) == 1) 
    return(rphast.simplify.list(lol[[1]]))
  ats <- attributes(lol)
  isMatrix <- (!is.null(ats$isMatrix) && ats$isMatrix==TRUE)
  isDataFrame <- (!is.null(ats$isDataFrame) && ats$isDataFrame==TRUE)
  if (isMatrix || isDataFrame) {
    if (!is.null(lol$row.names)) {
      rowNames <- lol$row.names
      lol$row.names <- NULL
    } else rowNames <- NULL
    if (isMatrix) lol <- as.matrix(as.data.frame(lol))
    if (isDataFrame) lol <- as.data.frame(lol)
    if (!is.null(rowNames)) row.names(lol) <- rowNames
  } else {
    for (i in 1:length(lol))
      lol[[i]] <- rphast.simplify.list(lol[[i]])
  }
  lol
}
