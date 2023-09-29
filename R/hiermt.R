#' Fit linear models and hierarchically adjust p-values
#'
#' hiermt is used to fit multiple univariate linear models and perform hierarchical-based multiple testing adjustment to the resulting p-values.
#'
#' @param formula A model formula object. Left hand side contains response variable(s) and right hand side contains explanatory variable(s).
#' @param data  A data frame containing variables listed in the formula.
#' @param global_test Global test to use when testing the hypotheses. This should be one of "bonferroni", "ghc", "gbj".
#' @param alpha Probability of Type I error. Default is set to 1.
#' @param mult_comp Logical: TRUE or FALSE. Whether to compute multiple comparisons.
#'
#' @returns `dend` A dendrogram with the hierarchically-adjusted p-values.
#' @returns `leaves` A dataframe containing the adjusted p-values for the terminal node hypotheses.
#' @returns `all` A dataframe containing all attributes of each node in the hierarchy.
#' @export
#'
#' @import stats
#' @importFrom GBJ GHC GBJ
#' @importFrom data.table data.table := fcase setkey
#' @importFrom dendextend nnodes partition_leaves which_node which_leaf
#' @importFrom collapse fsubset
#' @importFrom car Anova
#' @importFrom emmeans emmeans
#'
#' @examples
#' set.seed(1)
#' c
#' df <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x = sample(rep(0:1, c(n / 2, n / 2))))
#' hiermt(formula = cbind(y1, y2) ~ x, data = df, global_test = "bonferroni", alpha = 0.05)
#' hiermt(formula = . ~ x, data = df, global_test = "bonferroni", alpha = 0.05)
hiermt <- function(formula,
                   data = NULL,
                   global_test = "ghc",
                   alpha = 1L,
                   mult_comp = FALSE) {

  # Fit linear models

  lhs <- formula[[2L]]

  rhs <- formula[[3L]]

  if (lhs == rhs) {
    stop("The left-hand side and right-hand side of the
         formula can not be the same.")
  }

  Terms <- terms(formula, data = data)

  termlabels <- attr(Terms, "term.labels")

  if (lhs == ".") {
    if (missing(data)) {
      stop("\'.\' in formula but no data argument supplied.")
    }
    predictors <- as.character(attr(Terms, "variables"))[-c(1L, 2L)]
    response_names <- lapply(as.list(setdiff(colnames(data),
                                             predictors)),
                             as.name)
  } else {
    if (is.symbol(lhs)) {
      stop("Must have two or more responses to cluster.")
    }
    response_names <- lhs
    response_names[[1L]] <- NULL
  }

  model_attr <- data.table(
    response_names = as.character(response_names)
  )

  setkey(model_attr, response_names)

  model_attr[, formulas := lapply(
    response_names,
    function(x) {
      reformulate(termlabels = termlabels,
                  response = x)
    }
  )]

  model_attr[, modelframes := lapply(
    formulas,
    function(x) {
      model.frame(formula = x,
                  data = data)
    }
  )]

  model_attr[, models := lapply(
    formulas,
    function(x) {
      lm(x, data = data)
    }
  )]

  model_attr[, responses := lapply(
    modelframes, model.response
  )]

  model_attr[, model_factors := vapply(
    models,
    function(x) {
      sum(lengths(x$xlevels))
    },
    1
  )]

  model_attr[, anovas := lapply(
    models, Anova
  )]

  model_attr[, test_stats := vapply(
    anovas,
    function(x) x$`F value`[1],
    1
  )]

  model_attr[, pvalues := vapply(
    anovas,
    function(x) x$`Pr(>F)`[1],
    1
  )]

  # Hierarchical testing

  Q <- nrow(model_attr)

  Y <- matrix(
    unlist(
      model_attr[, responses]
    ),
    ncol = Q
  )
  colnames(Y) <- model_attr[, response_names]

  cor_mat <- as.matrix(
    cor(
      Y, method = "spearman"
    )
  )

  hc <- hclust(
    as.dist(
      sqrt(
        2L * (1L - abs(cor_mat))
      )
    )
  )

  dend <- as.dendrogram(hc)

  nn <- as.integer(nnodes(dend))


  global_test <- match.arg(
    global_test,
    c("bonferroni", "ghc", "gbj")
  )

  get_global_pvalue <- function(ind_pvalue,
                                ind_test_stat,
                                corr,
                                test) {
    ind_pvalue_len <- length(ind_pvalue)
    global_pvalue <- fcase(
      ind_pvalue_len == 1L, ind_pvalue[1],
      test == "bonferroni", min(ind_pvalue * ind_pvalue_len),
      test == "ghc", GHC(ind_test_stat, corr)$GHC_pvalue,
      test == "gbj", GBJ(ind_test_stat, corr)$GBJ_pvalue,
      default = 1
    )
    global_pvalue <- min(global_pvalue, 1)
    return(global_pvalue)
  }

  hier_attr <- data.table(
    node_counter = 1L:nn,
    node = partition_leaves(dend),
    is_node_leaf = as.integer(which_leaf(dend))
  )

  hier_attr[, descendants := lapply(
    node,
    function(x){
      which(sapply(node, function(z) all(z %in% x)))
    }
  )]

  hier_attr[, ancestors := lapply(
    node,
    function(x){
      which_node(dend, x, max_id = FALSE)
    }
  )]

  hier_attr[, parent := lapply(
    ancestors,
    function(x) {
      x[length(x) - 1L]
    }
  )]

  hier_attr[, sibling := mapply(
    function(x,z) {
      setdiff(which(hier_attr[, parent] == x), z)
    },
    hier_attr[, parent],
    hier_attr[, node_counter]
  )]

  hier_attr[, is_sibling_leaf := sapply(
    sibling,
    function(x) {
      hier_attr[x, is_node_leaf]
    }
  )]

  hier_attr[, adjustment := mapply(
    function(x, z) {
      max((Q / (length(x) + z)), 1L)
    },
    hier_attr[, node],
    hier_attr[, is_sibling_leaf]
  )]

  hier_attr[, node_pvalue := sapply(
    node,
    function(x) {
      get_global_pvalue(
        model_attr[response_names %in% x, pvalues],
        model_attr[response_names %in% x, test_stats],
        fsubset(cor_mat, x, x),
        global_test
      )
    }
  )]

  hier_attr[, adj_pvalue := pmin(
    node_pvalue * adjustment,
    rep(1L, nn)
  )]

  hier_attr[, h_adj_pvalue := sapply(
    ancestors,
    function(x){
      max(hier_attr[x, adj_pvalue])
    }
  )]

  hier_attr[, detected := h_adj_pvalue <= alpha]

  hier_attr[detected == TRUE,
            labels := round(h_adj_pvalue, 3)]

  emmeans_formula <- as.formula(
    paste0(
      "pairwise~",
      deparse(rhs),
      collapse = " "
    )
  )

  leaf_attr <- hier_attr[is_node_leaf == 1L,]

  leaf_attr[, node := unlist(node)]

  setkey(leaf_attr, node)

  leaf_attr <- model_attr[leaf_attr]

  if (mult_comp) {
    if (any(model_attr[, model_factors] < 3)) {
      stop("Need more than two levels in group to perform multiple comparisons.")
    }

    leaf_attr[detected == TRUE,
              emmeanss := lapply(
                models,
                function(x){
                  summary(emmeans(x, specs = emmeans_formula))$contrasts
                }
              )]

    leaf_attr[detected == TRUE,
              mult_comp_pvalues := lapply(
                emmeanss,
                function(x){
                  vapply(x$p.value,
                         function(z) {min(round(z,4),1L)},
                         1)
                }
              )]

    leaf_attr[detected == TRUE,
              mult_comp_contrasts := lapply(
                emmeanss,
                function(x) x$contrast
              )]
  }

  structure(
    list(
      hier_attr = hier_attr,
      leaf_attr = leaf_attr,
      dend = dend,
      mult_comp = mult_comp,
      alpha = alpha,
      Call = match.call()
    ),
    class = "hiermt"
  )


}
