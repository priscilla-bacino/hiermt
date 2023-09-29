#' Plot Hierarchical Model Testing Results
#'
#' This function plots the results of hierarchical model testing.
#'
#' @param x An object of class 'hiermt'.
#' @param ... Additional arguments.
#'
#' @importFrom ggdendro dendro_data label segment
#' @importFrom data.table fcase fifelse
#' @importFrom dendextend get_nodes_xy
#' @importFrom grDevices rainbow
#' @import ggplot2
#' @import patchwork
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' df <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x = sample(rep(0:1, c(n / 2, n / 2))))
#' hmt <- hiermt(formula = cbind(y1, y2) ~ x, data = df, global_test = "bonferroni", alpha = 0.05)
#' plot(hmt)

plot.hiermt <- function(x, ...) {

  hier_attr <- x$hier_attr

  leaf_attr <- x$leaf_attr

  dend <- x$dend

  mult_comp = x$mult_comp

  alpha = x$alpha

  Q <- nrow(leaf_attr)

  label_size <- function(Q) {
    fcase(
      Q <= 10, 3,
      Q <= 50, 2.5,
      Q <= 100, 2,
      default = 1
    )
  }

  dend_data <- dendro_data(dend)

  branch_dt <- data.table(
    segment(dend_data),
    line_color = rep(hier_attr[-1, detected], each = 2)
  )

  node_labels_dt <- data.table(
    get_nodes_xy(dend),
    show_label = hier_attr[, labels]
  )

  leaf_labels_dt <- data.table(
    label(dend_data)
  )
  setkey(leaf_labels_dt,label)

  leaf_attr <- leaf_attr[leaf_labels_dt]

  base_dend <- ggplot() +
    geom_segment(
      data = branch_dt,
      aes(x = x, y = y, xend = xend, yend = yend, color = line_color),
      linewidth = 0.3
    ) +
    scale_color_manual(
      values = c("#bfbfbf", "#000000")
    ) +
    theme_void() +
    #coord_flip() +
    theme(legend.position = "none")

  top_tier_dend <- base_dend +
    geom_label(
      data = node_labels_dt,
      aes(x = V1, y = V2, label = show_label),
      size = label_size(Q),
      na.rm = TRUE
    ) +
    geom_text(
      data = leaf_attr,
      aes(x = x, y = y, label = response_names, color = detected),
      angle = 90,
      hjust = "outward",
      size = label_size(Q),
      na.rm = TRUE,
      nudge_y = -0.05,
    )


  if (mult_comp) {

    contrasts_no <- max(
      sapply(leaf_attr[, mult_comp_contrasts],
             length)
    )

    leaf_attr <- leaf_attr[is.null(mult_comp_pvalues),
                           mult_comp_pvalues := NA]

    leaf_attr[detected == TRUE,
              mult_comp_detected := lapply(
                mult_comp_pvalues,
                function(x) {
                  sapply(x, function(z) {fifelse(z <= alpha, 0, NA_integer_)})
                }
              )]

    mult_comp_points <- data.table(
      plot_no = 1:contrasts_no,
      point_color = rainbow(contrasts_no),
      contrasts = leaf_attr[
        sapply(
          leaf_attr[, mult_comp_contrasts],
          function(x) !is.null(x)
        ), mult_comp_contrasts
      ][[1]]
    )

    mult_comp_points[, plots := mapply(
      function(z,v) {
        geom_point(
          data = leaf_attr,
          aes(
            x = x,
            y = sapply(mult_comp_detected, function(x) {
              if (is.null(x)) NA else x[z]
            }) - (z/50)
          ),
          na.rm = TRUE,
          color = v
        )
      },
      plot_no,
      point_color
    )]

    followup_tier_dend <- base_dend + ggtitle('Follow-up tier')

    for (geom in mult_comp_points[, plots]) {
      followup_tier_dend <- followup_tier_dend + geom
    }

    top_tier_dend <- top_tier_dend + ggtitle('Top tier')

    hier_plot <- top_tier_dend / followup_tier_dend

  } else {

    hier_plot <- top_tier_dend

  }

  return(hier_plot)

}
