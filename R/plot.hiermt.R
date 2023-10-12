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

  dend <- x$dend

  mult_comp = x$mult_comp

  alpha = x$alpha

  nrow_hier_attr <- nrow(hier_attr)

  label_size <- function(n_row) {
    fcase(
      n_row <= 20, 3,
      n_row <= 100, 2.5,
      n_row <= 200, 2,
      default = 1
    )
  }

  hier_attr[, detected := h_adj_pvalue < alpha]

  hier_attr[detected == TRUE, labels := round(h_adj_pvalue, 3)]

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

  setkey(leaf_labels_dt, label)

  leaf_attr <- data.table(
    label = unlist(hier_attr[is_node_leaf == 1, node]),
    detected = hier_attr[is_node_leaf == 1, h_adj_pvalue] < alpha
  )

  setkey(leaf_attr, label)

  leaf_attr <- leaf_attr[leaf_labels_dt]

  base_plot <- ggplot() +
    geom_segment(
      data = branch_dt,
      aes(x = x, y = y, xend = xend, yend = yend, color = line_color),
      linewidth = 0.3
    ) +
    geom_label(
      data = node_labels_dt,
      aes(x = V1, y = V2, label = show_label),
      size = label_size(nrow_hier_attr),
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = c("#bfbfbf", "#000000")
    ) +
    coord_cartesian() +
    theme_void() +
    theme(legend.position = "none")

  if (mult_comp) {

    grid_attr <- x$grid_attr

    contrast_labels <- grid_attr[, mult_comp_contrasts][[1]]

    contrast_n <- length(contrast_labels)

    y_seq <- -seq(from = 0.05, by = 0.03, length.out = contrast_n)

    contrast_pvalues <- data.table(
      response_names = rep(grid_attr[, response_names], each = contrast_n),
      detected = unlist(grid_attr[, h_adj_mult_comp_pvalues]) < alpha
    )

    setkey(contrast_pvalues, response_names)

    contrast_pvalues <- contrast_pvalues[leaf_labels_dt]

    contrast_pvalues[, y := y + y_seq]

    hier_plot <- base_plot +
      geom_tile(
        data = contrast_pvalues,
        aes(x = x,
            y = y,
            fill = detected
        ),
        color = "white",
        lwd = 1,
        linetype = 1
      ) +
      scale_fill_manual(
        values = c("#efefef", "#000000")
      ) + geom_text(
        data = leaf_attr,
        aes(x = x,
            y = min(contrast_pvalues$y),
            label = label,
            color = detected
        ),
        angle = 90,
        hjust = "outward",
        size = label_size(nrow_hier_attr),
        na.rm = TRUE,
        nudge_y = -0.05
      ) +
      geom_text(
        aes(x = 0,
            y = y_seq,
            label = contrast_labels
        ),
        hjust = "outward",
        size = label_size(nrow_hier_attr),
        na.rm = TRUE
      )

  } else {

    hier_plot <- base_plot +
      geom_text(
        data = leaf_attr,
        aes(x = x, y = y, label = label, color = detected),
        angle = 90,
        hjust = "outward",
        size = label_size(nrow_hier_attr),
        na.rm = TRUE,
        nudge_y = -0.05
      )

  }

return(hier_plot)
}
