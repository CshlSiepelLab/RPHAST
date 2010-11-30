##' @export
##' @nord
freeall.rphast <- function() {
  invisible(.Call("rph_free_all"))
}

##' @export
##' @nord
.Call.rphast <- function(func, ...) {
  on.exit(freeall.rphast())
  .Call(func, ...)
}
