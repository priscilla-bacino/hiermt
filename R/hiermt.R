#' Hierarchical multiple testing
#'
#' Calculate and visualize correlation-based hierarchically adjusted p-values.
#'
#' @param pvalues A numeric vector of p-values of length Q.
#' @param test_stats  A numeric vector of test statistics of length Q.
#' @param cor_mat A Q*Q square matrix of the correlations between the p-values or test statistics.
#' @param global_test Global test to use when testing the intersection of hypotheses. This should be one of "bonferonni", "ghc", "gbj".
#' @param alpha Probability of Type I error. Default is set to 1.
#'
#' @returns `dend` A dendrogram with the hierarchically-adjusted p-values.
#' @returns `leaves` A dataframe containing the adjusted p-values for the individual hypotheses.
#' @returns `all` A dataframe containing all attributes of each node in the hierarchy.
#' @export
#'
#' @import stats
#' @import MASS
#' @import ggplot2
#' @import ggdendro
#' @import GBJ
#' @importFrom clusterGeneration rcorrmatrix
#' @importFrom magrittr %>%
#' @importFrom data.table data.table :=
#' @importFrom purrr map map2 map_dbl map2_dbl map_lgl
#' @importFrom dendextend nnodes partition_leaves which_node which_leaf get_nodes_xy
#' @importFrom collapse fsubset rapply2d
#' @importFrom dplyr mutate nth
#'
#' @examples
#' set.seed(1)
#' Q <- 50
#' test_stats <- rnorm(Q, mean = c(rep(0, Q/2), rep(3, Q/2)))
#' pvalues <- 2*pnorm(sort(-abs(test_stats)))
#' cor_mat <- clusterGeneration::rcorrmatrix(Q)
#' hiermt(pvalues, test_stats, cor_mat, global_test="bonferroni", alpha = 0.05)



hiermt <- function(pvalues, test_stats, cor_mat, global_test = "ghc", alpha=1){
  global_tests <- c("bonferroni", "ghc", "gbj")
  if (!global_test %in% global_tests)stop("Global test arguement should be one of \"bonferonni\", \"ghc\", \"gbj\".")

  Q <- length(pvalues)
  if(length(test_stats)!= Q)stop("P-value vector and test statistic vector have differing lengths.")
  if(any(dim(cor_mat) != c(Q,Q)))stop("Pairwise correlation matrix is of the wrong size.")

  hc <- hclust(as.dist(sqrt(2*(1-abs(cor_mat)))),method="ward.D2")
  dend <- as.dendrogram(hc)
  nn <- nnodes(dend)

  S <- nodes <- isnodeleaf <- descendents <- ancestors <-
    parent <- sibling <- issiblingleaf <- adjustment <-
    node_pvalue <- adj_pvalue <- hadj_pvalue <- X1 <- X2 <-
    color_label <- ggendplot <- line_color <- show_label <-
    x <- xend <- y <- yend <- NULL

  calculate_global_pvalue <- function(p, t, corr, test){
    global_pvalue <- 1
    if(test == "bonferroni"){
      global_pvalue = min(p*length(p))
    }
    else if (length(p) == 1){
      global_pvalue = p
    }
    else if (test ==  "ghc"){
      global_pvalue = GHC(t, corr)$GHC_pvalue
    }
    else {
      global_pvalue = GBJ(t, corr)$GBJ_pvalue
    }
    return(global_pvalue)
  }

  attr <- data.table(S = c(1:nn),
                     nodes = partition_leaves(dend),
                     isnodeleaf = as.integer(
                       which_leaf(dend)))%>%mutate(
                       descendants = map(nodes,function(y)map_lgl(nodes,function(x)all(x%in%y))%>%which()),
                       ancestors = map(nodes, function(x)which_node(dend, x, max_id=FALSE)))%>%mutate(
                         parent = map(ancestors, function(y)nth(y,-2)))%>%mutate(
                           sibling = map2(parent, S, function(x,y)setdiff(which(parent==x),y)))%>%mutate(
                             issiblingleaf = map(sibling, function(x)isnodeleaf[x]))%>%mutate(
                               adjustment = map2_dbl(nodes, issiblingleaf, function(x,y) max((Q/(length(x)+y)),1)))%>%mutate(
                                 node_pvalue = map_dbl(nodes, function(x)min(calculate_global_pvalue(pvalues[x],test_stats[x],fsubset(cor_mat,x,x),global_test),1)))%>%mutate(
                                   adj_pvalue = map2_dbl(node_pvalue, adjustment, function(x,y) min(x*y,1)))%>%mutate(
                                     hadj_pvalue = map_dbl(ancestors, function(x) max(adj_pvalue[x])))
  attr_leaf <- attr[attr$isnodeleaf==1]
  attr_leaf$nodes = unlist(attr_leaf$nodes)
  branch_df <- segment(dendro_data(dend))%>%
    mutate(line_color = rep(ifelse(attr$hadj_pvalue<=alpha,"true","false")[-1], each=2))
  node_labels_df <- label(dendro_data(dend))%>%
    mutate(color_label = ifelse(attr_leaf$hadj_pvalue<=alpha,"true","false"))
  leaf_labels_df <- data.frame(get_nodes_xy(dend))%>%
    mutate(show_label = ifelse(attr$hadj_pvalue<=alpha, round(attr$hadj_pvalue,3),NA))
  ggdend_plot <- ggplot()+
    geom_segment(data = branch_df, aes(x = x, y = y, xend = xend, yend = yend, color= line_color), size=0.3)+
    geom_text(data = node_labels_df,
              aes(x = x, y = y, label = label, color = color_label),
              hjust ="outward",
              nudge_x = 0.15,
              nudge_y = -0.01,
              size = ifelse(Q<=10,3,ifelse(Q<=50,2.5,ifelse(Q<=100,2,1))))+
    geom_label(data = leaf_labels_df,
               aes(x=X1, y=X2, label= show_label),
               size = ifelse(Q<=10,2,ifelse(Q<=50,1.5,1)),
               na.rm=TRUE)+
    scale_color_manual(values = c("#bfbfbf",
                                  "#000000"))+
    #scale_y_continuous(limits = c(-0.8,4.5))+
    theme_void()+
    coord_flip()+
    theme(legend.position="none")
  leaves <- data.frame(index = attr_leaf$nodes, h_adj_pvalue = attr_leaf$hadj_pvalue)
  return(list(dend = ggdend_plot, leaves = leaves, all = attr))
}


