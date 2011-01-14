##' @export
##' @nord
freeall.rphast <- function() {
  invisible(.Call("rph_free_all"))
}

##' @export
##' @nord
.Call.rphast <- function(func, ...) {
  .Call("rph_new_mem_handler")
  on.exit(freeall.rphast())
  .Call(func, ...)
}
