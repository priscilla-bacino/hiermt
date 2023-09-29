#' Summary of Hierarchical Model Testing
#'
#' This function produces the summary of hierarchical model testing.
#'
#' @param object An object of class 'hiermt'.
#' @param ... Additional arguments.
#'
#' @returns `summary` summary of results
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' df <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x = sample(rep(0:1, c(n / 2, n / 2))))
#' hmt <- hiermt(formula = cbind(y1, y2) ~ x, data = df, global_test = "bonferroni", alpha = 0.05)
#' plot(hmt)
#'
summary.hiermt <- function(object, ...) {

  object$summary <- "Summary of Hierarchical Model Testing"
  print(object$summary)
}
