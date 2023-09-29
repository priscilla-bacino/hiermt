#' @export
#'
#'
print.hiermt <- function(x, ...) {
  cat("This is an object of class 'hiermt'.\n")
  cat("Call:", deparse(x$Call), "\n")
}
