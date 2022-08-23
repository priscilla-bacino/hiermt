#' Title: Hierarchical multiple testing
#'
#' Description: Calculate and visualize hierarchically bonferroni-adjusted
#' p-values.
#'
#' @param formula A formula specifying the model.
#' @param data A data frame in containing the variables named in `formula`.
#' @param global_test Global test to use when testing the intersection of hypotheses.
#' @param alpha Probability of Type I error.
#' @export
#'
#' @importFrom data.table data.table
#' @import MASS
#' @import tidyverse
#' @import dendextend
#' @import collapse
#'
#' @examples
#' set.seed(1)
#' library(MASS)
#' df <- data.frame(mvrnorm(n=80,mu=rep(0,100),Sigma=diag(100)), label = rep(c(0,1),40))
#' hiermt(formula=.~label, data = df, global_test="bonferroni", alpha = 0.05)


hiermt <- function(formula,data,global_test,alpha){
  m = 100
  return(m)
}

