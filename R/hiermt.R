#' Hierarchical multiple testing
#'
#' Calculate and visualize correlation-based hierarchically adjusted p-values.
#'
#' @param formula A model formula object. Left hand side contains response variable(s) and right hand side contains explanatory variable(s).
#' @param data  A data frame containing variables listed in the formula.
#' @param global_test Global test to use when testing the hypotheses. This should be one of "bonferonni", "ghc", "gbj".
#' @param alpha Probability of Type I error. Default is set to 1.
#' @param show_plot Logical. Default set to `FALSE`, if `TRUE` a dendrogram showing the results is returned.
#'
#' @returns `dend` A dendrogram with the hierarchically-adjusted p-values.
#' @returns `leaves` A dataframe containing the adjusted p-values for the terminal node hypotheses.
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
#' @importFrom data.table data.table := fifelse fcase
#' @importFrom purrr map map2 map_dbl map2_dbl map_lgl transpose
#' @importFrom dendextend nnodes partition_leaves which_node which_leaf get_nodes_xy
#' @importFrom collapse fsubset rapply2d
#' @importFrom dplyr mutate nth bind_cols rename arrange
#' @importFrom tidyr unnest_wider
#' @importFrom car Anova
#' @importFrom emmeans emmeans
#' @importFrom multcomp cld
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' df <- data.frame(y1=rnorm(n), y2=rnorm(n),x=sample(rep(0:1,c(n/2,n/2))))
#' hiermt(formula=cbind(y1,y2)~x, data=df, global_test="bonferroni", alpha = 0.05)
#' hiermt(formula=.~x, data=df, global_test="bonferroni", alpha = 0.05)


hiermt <- function(formula,
                   data = NULL,
                   global_test = "ghc",
                   alpha = 1L,
                   show_plot = FALSE){

  # Check if global test is one of three: "bonferroni", "ghc", "gbj".
  if (!global_test %in% c("bonferroni", "ghc", "gbj"))stop("Global test
                                          arguement should be one of
                                          \"bonferonni\", \"ghc\",
                                          \"gbj\".")
  # Assign lhs and rhs of from formula.
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  # Check if lhs and rhs are the same.
  if (lhs == rhs)stop(
    "The left-hand side and right-hand side of the formula can not be
    the same."
  )

  # Define terms and termlables for rhs of the formula.
  Terms <- terms(formula, data = data)
  termlabels <- attr(Terms, "term.labels")

  # Check lhs contains >1 response. If lhs = ".", assign all responses except
  # grouping/treatment variable.
  if (lhs == "."){
    if (missing(data))stop("\'.\' in formula but no 'data'
                        argument is supplied")
    predictors <- as.character(attr(Terms, "variables"))[-c(1L,2L)]
    response_names <- as.list(setdiff(colnames(data), predictors))
    response_names <- map(response_names,as.name)
  } else {
    if(length(lhs) == 1L)stop("Must have two or more responses to cluster.")
    response_names <- lhs
    response_names[[1L]] <- NULL
  }

  # Fit models for each response variable and obtain t-stats and p-values.
  Q <- length(response_names)
  response_names_c <- as.character(response_names)
  formulas <- map(response_names,
                  function(y)reformulate(termlabels = termlabels,
                                         response = y))
  modelframes<-map(formulas,
                   function(y)model.frame(formula = y,
                                          data = data))
  Y <- matrix(unlist(map(modelframes,
                         model.response)),
              ncol = Q)
  colnames(Y) <- response_names_c
  cor_mat <- as.matrix(cor(Y, method = "spearman"))
  models <- map(formulas,
                function(y)lm(y,data=data))
  names(models) <- response_names_c
  anovas <- map(models,Anova)
  names(anovas) <- response_names_c
  anovas <- anovas%>%transpose()%>%map(bind_cols)
  test_stats <- t(anovas$`F value`)[,1L]
  pvalues <- t(anovas$`Pr(>F)`)[,1L]

  # Multiple comparisons
  emmeans_formula <- as.formula(paste("pairwise~",deparse(rhs)))
  multiple_comparisons <- function(model){
    emmeanss <- emmeans(model,specs = emmeans_formula)
    multcomp_pvalues <- data.frame(emmeanss$contrasts)[,c(1,6)]
    return(multcomp_pvalues)
  }

  hc <- hclust(as.dist(sqrt(2L*(1L-abs(cor_mat)))),method="ward.D2")
  dend <- as.dendrogram(hc)
  nn <- as.integer(nnodes(dend))

  node_counter <- node <- is_node_leaf <- descendents <- ancestors <-
    parent <- sibling <- is_sibling_leaf <- adjustment <-
    node_pvalue <- adj_pvalue <- h_adj_pvalue <- detected <- labels <-
    V1 <- V2 <- label_color <- ggendplot <- line_color <- show_label <-
    x <- xend <- y <- yend <- show_cld <- cld_values <- pairwise_comp <-
    pairwise_pvalues <-  pairwise_comp_1 <- pairwise_comp_2 <- p.value <- NULL

  calculate_global_pvalue <- function(p, t, corr, test){
    global_pvalue <- fcase(
      length(p) == 1, p[1L],
      test == "bonferroni", min(p*length(p)),
      test == "ghc", GHC(t, corr)$GHC_pvalue,
      test == "gbj", GBJ(t, corr)$GBJ_pvalue,
      default = 1
      )
    return(global_pvalue)
  }

  attr <- data.table(node_counter = c(1L:nn),
                     node = partition_leaves(dend),
                     is_node_leaf = as.integer(which_leaf(dend)))%>%
    mutate(descendants = map(node,
                             function(y)map_lgl(node,
                                                function(x)all(x%in%y))%>%
                               which()),
           ancestors = map(node,
                           function(x)which_node(dend, x, max_id=FALSE)))%>%
    mutate(parent = map(ancestors,
                        function(y)nth(y,-2L)))%>%
    mutate(sibling = map2(parent,
                          node_counter,
                          function(x,y)setdiff(which(parent==x),y)))%>%
    mutate(is_sibling_leaf = map(sibling,
                                 function(x)is_node_leaf[x]))%>%
    mutate(adjustment = map2_dbl(node,
                                 is_sibling_leaf,
                                 function(x,y) max((Q/(length(x)+y)),1)))%>%
    mutate(node_pvalue = map_dbl(node,
                                 function(x)min(calculate_global_pvalue(
                                   pvalues[x],
                                   test_stats[x],
                                   fsubset(cor_mat,x,x),
                                   global_test),1L)))%>%
    mutate(adj_pvalue = map2_dbl(node_pvalue,
                                 adjustment,
                                 function(x,y) min(x*y,1L)))%>%
    mutate(h_adj_pvalue = map_dbl(ancestors,
                                  function(x) max(adj_pvalue[x])),
           detected = fifelse(h_adj_pvalue<=alpha,
                             TRUE, FALSE),
           labels = map_dbl(h_adj_pvalue,
                            function(x) fifelse(x <= alpha,
                                                round(x,3L),
                                                NA_integer_)))



  attr_leaf <- attr[is_node_leaf==1,]%>%
    mutate(node = unlist(node),
           multcomp_pvalues = fifelse(detected==TRUE,
             map2(node,
                  is_sibling_leaf,
                  function(x,y)multiple_comparisons(models[[x]])%>%
                    mutate(p.value = map_dbl(p.value,
                                             function(z)round(
                                               min((z*Q/(y+1)),1),4)))
                  ),
             list(NA)))

  ggdend_plot <- NULL

  if (show_plot){
    dend_data <- dendro_data(dend)

    branch_df <- data.table(
    segment(dend_data),
    line_color = rep(attr$detected[-1],each=2))

    node_labels_df <- data.table(
    label(dend_data),
    label_color = attr_leaf$detected,
    show_cld = attr_leaf$multcomp_pvalues)

    leaf_labels_df <-data.table(
      get_nodes_xy(dend),
      show_label = attr$labels)

    ggdend_plot <- ggplot()+
      geom_segment(data = branch_df,
                   aes(x = x,
                       y = y,
                       xend = xend,
                       yend = yend,
                       color = line_color),
                   size=0.3)+
      geom_text(data = node_labels_df,
                aes(x = x,
                    y = y,
                    label = label,
                    color = label_color),
                hjust ="outward",
                nudge_x = 0.15,
                nudge_y = -0.01,
                size = fifelse(Q<=10,
                              3,
                              fifelse(Q<=50,
                                     2.5,
                                     fifelse(Q<=100,
                                            2,
                                            1))))+
      geom_text(data = node_labels_df,
                aes(x = x,
                    y = y,
                    label = show_cld,
                    color = label_color),
                hjust ="outward",
                nudge_x = 0.15,
                nudge_y = -0.10,
                size = fifelse(Q<=10,
                              3,
                              fifelse(Q<=50,
                                     2.5,
                                     fifelse(Q<=100,
                                            2,
                                            1))))+
      geom_label(data = leaf_labels_df,
                 aes(x=V1,
                     y=V2,
                     label= show_label),
                 size = fifelse(Q<=10,
                               2,
                               fifelse(Q<=50,
                                      1.5,
                                      1)),
                 na.rm=TRUE)+
      scale_color_manual(values = c("#bfbfbf",
                                    "#000000"))+
      #scale_y_continuous(limits = c(-0.6,5))+
      theme_void()+
      coord_flip()+
      theme(legend.position="none")
  }
  #leaves <- attr_leaf[,list(node,h_adj_pvalue,detected,pairwise_pvalues)]

  #data.table(index = attr_leaf$node,h_adj_pvalue = attr_leaf$h_adj_pvalue)

  return(list(dend = ggdend_plot,
              leaves = attr_leaf,
              all = attr))
}


