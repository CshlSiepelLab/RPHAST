##' @export
##' @nord
freeall.rphast <- function() {
  invisible(.Call("rph_free_all"))
}
