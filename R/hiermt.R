#' Fit linear models and hierarchically adjust p-values
#'
#' hiermt is used to fit multiple univariate linear models and perform hierarchical-based multiple testing adjustment to the resulting p-values.
#'
#' @param formula An object of class \link{formula} or one that can be coerced to that class. Left hand side of the formula object contains response variable(s)and right hand side contains explanatory variable(s).
#' @param data  An optional data frame containing variables listed in the formula. If data argument is not supplied, the variables should be explicitly referenced in the formula.
#' @param global_test Global test to use when combining p-values. This should be one of "bonferroni", "ghc", "gbj".
#' @param alpha Probability of Type I error. Default is set to 1.
#' @param linkage Agglomeration method used for hierarchical clustering in \link{hclust}. Defaults to Ward's method.
#' @param mult_comp Logical: TRUE or FALSE. Whether to compute multiple comparisons.
#'
#' @returns An object of class hiermt which describes the hierachcical testing process. The object is a list with components:
#' @returns `hier_attr` Attributes of the hierarchical structure.
#' @returns `grid_attr` A subset of hier_attr showing attributes of the terminal nodes only.
#' @export
#'
#' @import stats
#' @importFrom GBJ GHC GBJ
#' @importFrom data.table data.table := fcase setkey
#' @importFrom dendextend nnodes partition_leaves which_leaf
#' @importFrom collapse fsubset
#' @importFrom car Anova
#' @importFrom emmeans emmeans
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' df <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x = sample(rep(0:1, c(n / 2, n / 2))))
#' hiermt(formula = cbind(y1, y2) ~ x, data = df, global_test = "bonferroni", alpha = 0.05)
#' hiermt(formula = . ~ x, data = df, global_test = "bonferroni", alpha = 0.05)
#'
hiermt <- function(formula,
                   data = NULL,
                   global_test = "ghc",
                   alpha = 0.05,
                   linkage = "ward.D2",
                   mult_comp = FALSE) {

  global_test <- match.arg(
    global_test,
    c("bonferroni",
      "ghc",
      "gbj"
    )
  )

  linkage <- match.arg(
    linkage,
    c("ward.D2",
      "ward.D",
      "single",
      "complete",
      "average",
      "mcquitty",
      "median",
      "centroid"
    )
  )

  lhs <- formula[[2L]]

  rhs <- formula[[3L]]

  if (lhs == rhs) {
    stop("The left-hand side and right-hand side of the
         formula can not be the same.")
  }

  Terms <- terms(formula, data = data)

  termlabels <- attr(Terms, "term.labels")

  factor_names <- as.character(attr(Terms, "variables"))[-c(1L, 2L)]

  if (lhs == ".") {
    if (missing(data)) {
      stop(
        "\'.\' in formula but no data argument supplied."
      )
    }

    response_names <- lapply(as.list(setdiff(colnames(data),
                                             factor_names)),
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
    ), method = linkage
  )

  dend <- as.dendrogram(hc)

  nn <- as.integer(nnodes(dend))

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
      which(sapply(node, function(z) all(x %in% z)))
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
    parent,
    node_counter
  )]

  hier_attr[, is_sibling_leaf := sapply(
    sibling,
    function(x) {
      is_node_leaf[x]
    }
  )]

  get_adjustment <- function(Q, pairwise, node, is_sibling_leaf){
    ifelse(
      pairwise,
      Q / length(node),
      Q / (length(node) + sum(is_sibling_leaf))
    )
  }


  hier_attr[, adjustment := mapply(
    function(x, z) {
      get_adjustment(Q, mult_comp, x, z)
    },
    node,
    is_sibling_leaf
  )]

  hier_attr[, select_pvalues := lapply(
    node,
    function(x){
      model_attr[.(x), pvalues]
    }
  )]

  hier_attr[, select_teststats := lapply(
    node,
    function(x){
      model_attr[.(x), test_stats]
    }
  )]

  hier_attr[, select_cor := lapply(
    node,
    function(x){
      fsubset(cor_mat, x, x)
    }
  )]

  get_global_pvalue <- function(pvalues, test_stats, cor, test){

    if (length(pvalues) == 1L) {

      return(pvalues[1])

    }

    switch(
      test,
      "bonferroni" = min(pvalues) * length(pvalues),
      "ghc" = GHC(test_stats, cor)$GHC_pvalue,
      "gbj" = GBJ(test_stats, cor)$GBJ_pvalue,
      1
    )

  }

  hier_attr[, node_pvalue := mapply(
    get_global_pvalue,
    select_pvalues,
    select_teststats,
    select_cor,
    MoreArgs = list(test = global_test)
  )]

  hier_attr[, adj_pvalue := pmin(
    node_pvalue * adjustment,
    rep(1L, nn)
  )]

  hier_attr[, h_adj_pvalue := sapply(
    ancestors,
    function(x){
      max(adj_pvalue[x])
    }
  )]

  grid_attr <- "No grid attributes, multiple comparisons were not performed."

  if (mult_comp) {
    if (any(model_attr[, model_factors] < 3)) {
      stop("Need more than two levels for a factor to perform multiple comparisons.")
    }

    grid_attr <- hier_attr[is_node_leaf == 1L,]

    grid_attr[, node := unlist(node)]

    setkey(grid_attr, node)

    grid_attr <- model_attr[grid_attr]

    emmeans_formula <- as.formula(
      paste(
        "pairwise", "~", deparse(rhs)
      )
    )

    grid_attr[, emmeanss := lapply(
      models,
      function(x){
        emmeans(x, specs = emmeans_formula, adjust = "none")
      }
    )]

    grid_attr[, mult_comp_pvalues := lapply(
      emmeanss,
      function(x){
        summary(x)$contrasts$p.value
      }
    )]

    grid_attr[, adj_mult_comp_pvalues := mapply(
      function(x, z) {
        pmin((x * Q * (z - 1) * (z - 2) / 2), 1L)
      },
      mult_comp_pvalues,
      model_factors,
      SIMPLIFY = FALSE
    )]

    grid_attr[, h_adj_mult_comp_pvalues := mapply(
      function(x, z) {
        pmax(x, z)
      },
      adj_mult_comp_pvalues,
      h_adj_pvalue,
      SIMPLIFY = FALSE
    )]

    grid_attr[, mult_comp_contrasts := lapply(
      emmeanss,
      function(x) summary(x)$contrasts$contrast
    )]

    grid_attr <- grid_attr[, .(
    response_names,
    emmeanss,
    mult_comp_pvalues,
    adj_mult_comp_pvalues,
    h_adj_mult_comp_pvalues,
    mult_comp_contrasts
  )]

  }

  hier_attr[, c("adjustment",
                "select_pvalues",
                "select_teststats",
                "select_cor"
              ) := NULL]

  structure(
    list(
      hier_attr = hier_attr,
      grid_attr = grid_attr,
      dend = dend,
      mult_comp = mult_comp,
      alpha = alpha,
      Call = match.call()
    ),
    class = "hiermt"
  )


}
