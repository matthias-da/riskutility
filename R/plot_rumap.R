#' Plot method for rumap objects
#'
#' Provides multiple visualization approaches for multivariate Risk-Utility evaluation
#' as described in "Beyond the Trade-off Curve" (Thees, Müller, Templ 2026).
#'
#' @param x An object of class "rumap"
#' @param y Not used
#' @param which Integer vector specifying which plots to produce:
#'   \itemize{
#'     \item 1: Composite scatterplot (R-U map with Pareto front)
#'     \item 2: Heatmap (all measures x SDGs)
#'     \item 3: Dot plot (risk and utility facets)
#'     \item 4: Parallel coordinates plot
#'     \item 5: Radial/radar chart
#'     \item 6: PCA biplot (joint PCA)
#'     \item 7: Blockwise PCA scatterplot
#'   }
#' @param ... Additional arguments passed to plotting functions
#' @param show_pareto Logical, whether to highlight Pareto-optimal SDGs. Default TRUE.
#' @param show_labels Logical, whether to show SDG labels. Default TRUE.
#' @param col_pareto Color for Pareto-optimal points. Default "steelblue".
#' @param col_dominated Color for dominated points. Default "gray60".
#'
#' @details
#' The visualization approaches are designed to reveal different aspects of the
#' Risk-Utility trade-off:
#'
#' \strong{Plot 1 - Composite Scatterplot:} Classic R-U map showing composite
#' risk vs utility scores. Pareto front is highlighted. Good for overall comparison.
#'
#' \strong{Plot 2 - Heatmap:} Shows all measures for all SDGs in a matrix format.
#' Useful for identifying patterns across measures.
#'
#' \strong{Plot 3 - Dot Plot:} Separates risk and utility measures into facets.
#' Shows individual measure values for detailed comparison.
#'
#' \strong{Plot 4 - Parallel Coordinates:} Each SDG is a line crossing all measures.
#' Reveals trade-offs and patterns across dimensions.
#'
#' \strong{Plot 5 - Radial/Radar Chart:} Polygonal profile for each SDG.
#' Provides gestalt view of multivariate performance.
#'
#' \strong{Plot 6 - PCA Biplot:} Joint PCA of all measures. Shows measure
#' correlations (arrows) and SDG positions simultaneously.
#'
#' \strong{Plot 7 - Blockwise PCA:} Separate PCA for risk and utility blocks.
#' X-axis = PC1(utility), Y-axis = PC1(risk).
#'
#' @seealso \code{\link{rumap}} for creating rumap objects
#'
#' @return The \code{rumap} object \code{x}, invisibly.
#' @family utility
#' @author Matthias Templ
#' @export
#' @importFrom graphics plot points lines polygon segments text legend axis mtext abline arrows
#' @importFrom graphics par layout image
#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom stats prcomp
plot.rumap <- function(x, y = NULL, which = 1, ...,
                       show_pareto = TRUE,
                       show_labels = TRUE,
                       col_pareto = "steelblue",
                       col_dominated = "gray60") {

  # Validate which
  valid_plots <- 1:7
  which <- intersect(which, valid_plots)
  if (length(which) == 0) {
    stop("'which' must contain values from 1 to 7")
  }

  # Setup for multiple plots
  n_plots <- length(which)
  if (n_plots > 1) {
    if (n_plots == 2) {
      op <- par(mfrow = c(1, 2))
    } else if (n_plots <= 4) {
      op <- par(mfrow = c(2, 2))
    } else {
      op <- par(mfrow = c(3, 3))
    }
    on.exit(par(op))
  }

  # Generate requested plots
  for (w in which) {
    switch(w,
           `1` = plot_rumap_composite(x, show_pareto, show_labels, col_pareto, col_dominated, ...),
           `2` = plot_rumap_heatmap(x, show_pareto, ...),
           `3` = plot_rumap_dotplot(x, show_pareto, col_pareto, col_dominated, ...),
           `4` = plot_rumap_parallel(x, show_pareto, col_pareto, col_dominated, ...),
           `5` = plot_rumap_radar(x, show_pareto, col_pareto, ...),
           `6` = plot_rumap_pca_biplot(x, show_pareto, show_labels, col_pareto, col_dominated, ...),
           `7` = plot_rumap_blockwise_pca(x, show_pareto, show_labels, col_pareto, col_dominated, ...)
    )
  }

  invisible(x)
}


#' Plot 1: Composite Scatterplot (R-U Map)
#' @keywords internal
plot_rumap_composite <- function(x, show_pareto, show_labels, col_pareto, col_dominated, ...) {
  df <- x$composites
  is_pareto <- x$pareto

  # Colors and shapes
  cols <- ifelse(is_pareto, col_pareto, col_dominated)
  pchs <- ifelse(is_pareto, 19, 1)

  # Plot
  plot(df$utility_mean, df$risk_mean,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Composite Utility (higher = better)",
       ylab = "Composite Risk (lower = better)",
       main = "Risk-Utility Map",
       pch = pchs, col = cols, cex = 1.5, ...)

  # Add error bars if we have SD
  if ("risk_sd" %in% names(df) && "utility_sd" %in% names(df)) {
    for (i in 1:nrow(df)) {
      # Vertical (risk) error bars
      segments(df$utility_mean[i], df$risk_mean[i] - df$risk_sd[i],
               df$utility_mean[i], df$risk_mean[i] + df$risk_sd[i],
               col = cols[i], lwd = 0.5)
      # Horizontal (utility) error bars
      segments(df$utility_mean[i] - df$utility_sd[i], df$risk_mean[i],
               df$utility_mean[i] + df$utility_sd[i], df$risk_mean[i],
               col = cols[i], lwd = 0.5)
    }
  }

  # Draw Pareto front
  if (show_pareto && sum(is_pareto) > 1) {
    pareto_df <- df[is_pareto, ]
    pareto_df <- pareto_df[order(pareto_df$utility_mean), ]
    lines(pareto_df$utility_mean, pareto_df$risk_mean,
          col = col_pareto, lty = 2, lwd = 2)
  }

  # Labels
  if (show_labels) {
    text(df$utility_mean, df$risk_mean, labels = df$sdg,
         pos = 3, cex = 0.7, col = cols)
  }

  # Reference lines
  abline(h = 0.5, v = 0.5, col = "gray80", lty = 3)

  # Quadrant labels
  text(0.75, 0.25, "Best", col = "darkgreen", cex = 0.8, font = 2)
  text(0.25, 0.75, "Worst", col = "darkred", cex = 0.8, font = 2)

  # Legend
  legend("topright",
         legend = c("Pareto-optimal", "Dominated"),
         pch = c(19, 1),
         col = c(col_pareto, col_dominated),
         cex = 0.8, bg = "white")
}


#' Plot 2: Heatmap
#' @keywords internal
plot_rumap_heatmap <- function(x, show_pareto, ...) {
  if (is.null(x$normalized)) {
    warning("Heatmap requires normalized data")
    return(invisible(NULL))
  }

  # Prepare data matrix
  risk_cols <- x$risk_measures
  utility_cols <- x$utility_measures
  all_cols <- c(risk_cols, utility_cols)

  mat <- as.matrix(x$normalized[, all_cols, drop = FALSE])
  rownames(mat) <- x$sdg_names
  if (show_pareto) {
    rownames(mat)[x$pareto] <- paste0(rownames(mat)[x$pareto], " *")
  }

  # Color palette: risk (red), utility (blue)
  n_risk <- length(risk_cols)
  n_utility <- length(utility_cols)

  # Create color matrix
  col_palette <- colorRampPalette(c("white", "steelblue"))(100)

  # Transpose for image (rows become columns)
  mat_t <- t(mat[nrow(mat):1, , drop = FALSE])

  # Plot
  image(1:ncol(mat), 1:nrow(mat), mat_t,
        col = col_palette,
        xlab = "", ylab = "",
        main = "Risk-Utility Heatmap",
        axes = FALSE, ...)

  # Add axes
  axis(1, at = 1:ncol(mat), labels = colnames(mat), las = 2, cex.axis = 0.7)
  axis(2, at = 1:nrow(mat), labels = rev(rownames(mat)), las = 1, cex.axis = 0.7)

  # Add values as text
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      text(j, nrow(mat) - i + 1, sprintf("%.2f", mat[i, j]),
           cex = 0.6, col = ifelse(mat[i, j] > 0.5, "white", "black"))
    }
  }

  # Add separator between risk and utility
  abline(v = n_risk + 0.5, lwd = 2)

  # Labels
  mtext("Risk", side = 1, at = n_risk/2, line = 3, col = "firebrick", font = 2)
  mtext("Utility", side = 1, at = n_risk + n_utility/2, line = 3, col = "steelblue", font = 2)

  box()
}


#' Plot 3: Dot Plot
#' @keywords internal
plot_rumap_dotplot <- function(x, show_pareto, col_pareto, col_dominated, ...) {
  if (is.null(x$normalized)) {
    warning("Dot plot requires normalized data")
    return(invisible(NULL))
  }

  risk_cols <- x$risk_measures
  utility_cols <- x$utility_measures
  n_sdgs <- x$n_sdgs

  # Setup layout: risk on left, utility on right
  op <- par(mfrow = c(1, 2), mar = c(5, 8, 4, 2))
  on.exit(par(op))

  # Colors
  cols <- ifelse(x$pareto, col_pareto, col_dominated)

  # Risk panel
  if (length(risk_cols) > 0) {
    plot(NULL, xlim = c(0, 1), ylim = c(0.5, n_sdgs + 0.5),
         xlab = "Normalized Value", ylab = "",
         main = "Risk Measures", yaxt = "n", ...)
    axis(2, at = 1:n_sdgs, labels = x$sdg_names, las = 1, cex.axis = 0.8)

    for (i in 1:n_sdgs) {
      for (j in seq_along(risk_cols)) {
        val <- x$normalized[i, risk_cols[j]]
        pch_offset <- (j - (length(risk_cols) + 1) / 2) * 0.15
        points(val, i + pch_offset, pch = 14 + j, col = cols[i], cex = 1.2)
      }
    }

    legend("topright", legend = risk_cols, pch = 14 + seq_along(risk_cols),
           cex = 0.7, bg = "white")
    abline(v = 0.5, lty = 2, col = "gray50")
  }

  # Utility panel
  if (length(utility_cols) > 0) {
    plot(NULL, xlim = c(0, 1), ylim = c(0.5, n_sdgs + 0.5),
         xlab = "Normalized Value", ylab = "",
         main = "Utility Measures", yaxt = "n", ...)
    axis(2, at = 1:n_sdgs, labels = x$sdg_names, las = 1, cex.axis = 0.8)

    for (i in 1:n_sdgs) {
      for (j in seq_along(utility_cols)) {
        val <- x$normalized[i, utility_cols[j]]
        pch_offset <- (j - (length(utility_cols) + 1) / 2) * 0.15
        points(val, i + pch_offset, pch = 14 + j, col = cols[i], cex = 1.2)
      }
    }

    legend("topright", legend = utility_cols, pch = 14 + seq_along(utility_cols),
           cex = 0.7, bg = "white")
    abline(v = 0.5, lty = 2, col = "gray50")
  }
}


#' Plot 4: Parallel Coordinates
#' @keywords internal
plot_rumap_parallel <- function(x, show_pareto, col_pareto, col_dominated, ...) {
  if (is.null(x$normalized)) {
    warning("Parallel coordinates requires normalized data")
    return(invisible(NULL))
  }

  risk_cols <- x$risk_measures
  utility_cols <- x$utility_measures
  all_cols <- c(risk_cols, utility_cols)
  n_cols <- length(all_cols)
  n_sdgs <- x$n_sdgs

  # Setup plot
  plot(NULL, xlim = c(1, n_cols), ylim = c(0, 1),
       xlab = "", ylab = "Normalized Value",
       main = "Parallel Coordinates Plot",
       xaxt = "n", ...)

  axis(1, at = 1:n_cols, labels = all_cols, las = 2, cex.axis = 0.7)

  # Draw vertical axes
  for (i in 1:n_cols) {
    segments(i, 0, i, 1, col = "gray80")
  }

  # Separator between risk and utility
  abline(v = length(risk_cols) + 0.5, lwd = 2, col = "gray40")

  # Colors and line widths
  cols <- ifelse(x$pareto, col_pareto, col_dominated)
  lwds <- ifelse(x$pareto, 2, 1)
  ltys <- ifelse(x$pareto, 1, 2)

  # Draw lines for each SDG
  for (i in 1:n_sdgs) {
    vals <- as.numeric(x$normalized[i, all_cols])
    lines(1:n_cols, vals, col = cols[i], lwd = lwds[i], lty = ltys[i])
    points(1:n_cols, vals, col = cols[i], pch = 19, cex = 0.8)
  }

  # Legend
  legend("topright",
         legend = x$sdg_names,
         col = cols, lwd = lwds, lty = ltys,
         cex = 0.7, bg = "white")

  # Axis labels
  mtext("Risk", side = 1, at = length(risk_cols)/2, line = 3.5, col = "firebrick", font = 2)
  mtext("Utility", side = 1, at = length(risk_cols) + length(utility_cols)/2,
        line = 3.5, col = "steelblue", font = 2)
}


#' Plot 5: Radar/Spider Chart
#' @keywords internal
plot_rumap_radar <- function(x, show_pareto, col_pareto, ...) {
  if (is.null(x$normalized)) {
    warning("Radar chart requires normalized data")
    return(invisible(NULL))
  }

  # For radar chart, invert risk measures so higher = better for all
  risk_cols <- x$risk_measures
  utility_cols <- x$utility_measures
  all_cols <- c(risk_cols, utility_cols)
  n_vars <- length(all_cols)
  n_sdgs <- x$n_sdgs

  if (n_vars < 3) {
    warning("Radar chart requires at least 3 measures")
    return(invisible(NULL))
  }

  # Create inverted data (1 - risk so all are "higher = better")
  data <- x$normalized[, all_cols, drop = FALSE]
  for (col in risk_cols) {
    data[[col]] <- 1 - data[[col]]
  }

  # Angles for each axis
  angles <- seq(0, 2 * pi, length.out = n_vars + 1)[1:n_vars]

  # Setup plot
  plot(NULL, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3),
       xlab = "", ylab = "", asp = 1,
       main = "Radar Chart (higher = better)",
       axes = FALSE, ...)

  # Draw circular grid
  for (r in c(0.25, 0.5, 0.75, 1)) {
    theta <- seq(0, 2 * pi, length.out = 100)
    lines(r * cos(theta), r * sin(theta), col = "gray80", lty = 3)
  }

  # Draw axis lines
  for (i in 1:n_vars) {
    segments(0, 0, cos(angles[i]), sin(angles[i]), col = "gray60")
    # Labels
    label_col <- ifelse(all_cols[i] %in% risk_cols, "firebrick", "steelblue")
    text(1.15 * cos(angles[i]), 1.15 * sin(angles[i]),
         all_cols[i], cex = 0.7, col = label_col)
  }

  # Draw polygons for each SDG
  colors <- ifelse(x$pareto, col_pareto, "gray50")
  alphas <- ifelse(x$pareto, 0.3, 0.1)

  for (i in 1:n_sdgs) {
    vals <- as.numeric(data[i, ])
    # Convert to coordinates
    x_coords <- vals * cos(angles)
    y_coords <- vals * sin(angles)
    # Close polygon
    x_coords <- c(x_coords, x_coords[1])
    y_coords <- c(y_coords, y_coords[1])

    # Fill
    polygon(x_coords, y_coords,
            col = adjustcolor(colors[i], alpha.f = alphas[i]),
            border = colors[i], lwd = ifelse(x$pareto[i], 2, 1))
  }

  # Legend
  legend("topright",
         legend = x$sdg_names[order(x$pareto, decreasing = TRUE)],
         fill = adjustcolor(colors[order(x$pareto, decreasing = TRUE)], alpha.f = 0.3),
         border = colors[order(x$pareto, decreasing = TRUE)],
         cex = 0.7, bg = "white")
}


#' Plot 6: PCA Biplot (Joint)
#' @keywords internal
plot_rumap_pca_biplot <- function(x, show_pareto, show_labels, col_pareto, col_dominated, ...) {
  if (is.null(x$normalized)) {
    warning("PCA biplot requires normalized data")
    return(invisible(NULL))
  }

  risk_cols <- x$risk_measures
  utility_cols <- x$utility_measures
  all_cols <- c(risk_cols, utility_cols)

  if (length(all_cols) < 2) {
    warning("PCA requires at least 2 measures")
    return(invisible(NULL))
  }

  # Prepare data matrix
  data <- x$normalized[, all_cols, drop = FALSE]
  data <- data[complete.cases(data), , drop = FALSE]

  if (nrow(data) < 2) {
    warning("Not enough complete cases for PCA")
    return(invisible(NULL))
  }

  # Remove constant columns (variance = 0) to avoid prcomp error
  col_vars <- apply(data, 2, var, na.rm = TRUE)
  non_const_cols <- names(col_vars)[col_vars > 1e-10]
  if (length(non_const_cols) < 2) {
    warning("Not enough non-constant columns for PCA biplot")
    return(invisible(NULL))
  }
  data <- data[, non_const_cols, drop = FALSE]
  all_cols <- non_const_cols

  # Perform PCA
  pca <- prcomp(data, scale. = TRUE, center = TRUE)

  # Extract scores and loadings
  scores <- pca$x[, 1:2]
  loadings <- pca$rotation[, 1:2]

  # Variance explained
  var_exp <- summary(pca)$importance[2, 1:2] * 100

  # Colors
  cols <- ifelse(x$pareto, col_pareto, col_dominated)
  pchs <- ifelse(x$pareto, 19, 1)

  # Scaling for biplot
  scale_factor <- max(abs(scores)) / max(abs(loadings)) * 0.8

  # Setup plot
  plot(scores[, 1], scores[, 2],
       xlab = sprintf("PC1 (%.1f%%)", var_exp[1]),
       ylab = sprintf("PC2 (%.1f%%)", var_exp[2]),
       main = "PCA Biplot",
       pch = pchs, col = cols, cex = 1.5, ...)

  # Add origin
  abline(h = 0, v = 0, col = "gray80", lty = 2)

  # Labels for points
  if (show_labels) {
    text(scores[, 1], scores[, 2], labels = x$sdg_names,
         pos = 3, cex = 0.7, col = cols)
  }

  # Draw loading arrows
  for (i in 1:nrow(loadings)) {
    arrow_col <- ifelse(rownames(loadings)[i] %in% risk_cols, "firebrick", "steelblue")
    arrows(0, 0,
           loadings[i, 1] * scale_factor,
           loadings[i, 2] * scale_factor,
           col = arrow_col, lwd = 1.5, length = 0.1)
    text(loadings[i, 1] * scale_factor * 1.1,
         loadings[i, 2] * scale_factor * 1.1,
         rownames(loadings)[i],
         col = arrow_col, cex = 0.7)
  }

  # Legend
  legend("topright",
         legend = c("Pareto-optimal", "Dominated", "Risk measure", "Utility measure"),
         pch = c(19, 1, NA, NA),
         lty = c(NA, NA, 1, 1),
         col = c(col_pareto, col_dominated, "firebrick", "steelblue"),
         cex = 0.7, bg = "white")
}


#' Plot 7: Blockwise PCA
#' @keywords internal
plot_rumap_blockwise_pca <- function(x, show_pareto, show_labels, col_pareto, col_dominated, ...) {
  if (is.null(x$normalized)) {
    warning("Blockwise PCA requires normalized data")
    return(invisible(NULL))
  }

  risk_cols <- x$risk_measures
  utility_cols <- x$utility_measures

  if (length(risk_cols) < 2 || length(utility_cols) < 2) {
    warning("Blockwise PCA requires at least 2 measures in each block")
    return(invisible(NULL))
  }

  # Prepare data
  risk_data <- x$normalized[, risk_cols, drop = FALSE]
  utility_data <- x$normalized[, utility_cols, drop = FALSE]

  # Remove NAs
  complete_idx <- complete.cases(risk_data) & complete.cases(utility_data)
  if (sum(complete_idx) < 2) {
    warning("Not enough complete cases for blockwise PCA")
    return(invisible(NULL))
  }

  risk_data <- risk_data[complete_idx, , drop = FALSE]
  utility_data <- utility_data[complete_idx, , drop = FALSE]
  is_pareto_sub <- x$pareto[complete_idx]
  sdg_names_sub <- x$sdg_names[complete_idx]

  # Remove constant columns (variance = 0) to avoid prcomp error
  risk_vars <- apply(risk_data, 2, var, na.rm = TRUE)
  risk_non_const <- names(risk_vars)[risk_vars > 1e-10]
  utility_vars <- apply(utility_data, 2, var, na.rm = TRUE)
  utility_non_const <- names(utility_vars)[utility_vars > 1e-10]

  if (length(risk_non_const) < 2 || length(utility_non_const) < 2) {
    warning("Blockwise PCA requires at least 2 non-constant measures in each block")
    return(invisible(NULL))
  }

  risk_data <- risk_data[, risk_non_const, drop = FALSE]
  utility_data <- utility_data[, utility_non_const, drop = FALSE]
  risk_cols <- risk_non_const
  utility_cols <- utility_non_const

  # PCA on each block
  pca_risk <- prcomp(risk_data, scale. = TRUE, center = TRUE)
  pca_utility <- prcomp(utility_data, scale. = TRUE, center = TRUE)

  # PC1 scores
  pc1_risk <- pca_risk$x[, 1]
  pc1_utility <- pca_utility$x[, 1]

  # Flip signs if needed for interpretability
  # Risk: higher PC1 should mean higher risk
  if (cor(pc1_risk, rowMeans(risk_data)) < 0) {
    pc1_risk <- -pc1_risk
  }
  # Utility: higher PC1 should mean higher utility
  if (cor(pc1_utility, rowMeans(utility_data)) < 0) {
    pc1_utility <- -pc1_utility
  }

  # Variance explained
  var_risk <- summary(pca_risk)$importance[2, 1] * 100
  var_utility <- summary(pca_utility)$importance[2, 1] * 100

  # Colors
  cols <- ifelse(is_pareto_sub, col_pareto, col_dominated)
  pchs <- ifelse(is_pareto_sub, 19, 1)

  # Plot
  plot(pc1_utility, pc1_risk,
       xlab = sprintf("Utility PC1 (%.1f%% var)", var_utility),
       ylab = sprintf("Risk PC1 (%.1f%% var)", var_risk),
       main = "Blockwise PCA",
       pch = pchs, col = cols, cex = 1.5, ...)

  # Reference lines at mean
  abline(h = 0, v = 0, col = "gray80", lty = 2)

  # Labels
  if (show_labels) {
    text(pc1_utility, pc1_risk, labels = sdg_names_sub,
         pos = 3, cex = 0.7, col = cols)
  }

  # Draw Pareto front
  if (show_pareto && sum(is_pareto_sub) > 1) {
    pareto_order <- order(pc1_utility[is_pareto_sub])
    lines(pc1_utility[is_pareto_sub][pareto_order],
          pc1_risk[is_pareto_sub][pareto_order],
          col = col_pareto, lty = 2, lwd = 2)
  }

  # Quadrant labels
  usr <- par("usr")
  text(usr[2] * 0.8, usr[3] * 0.8, "High U\nLow R", col = "darkgreen", cex = 0.7)
  text(usr[1] * 0.8, usr[4] * 0.8, "Low U\nHigh R", col = "darkred", cex = 0.7)

  # Legend
  legend("topright",
         legend = c("Pareto-optimal", "Dominated"),
         pch = c(19, 1),
         col = c(col_pareto, col_dominated),
         cex = 0.7, bg = "white")

  # Add loading contribution bars in margin
  # (simplified version - just add text)
  mtext(paste("Risk loadings:", paste(risk_cols, collapse = ", ")),
        side = 1, line = 4, cex = 0.6, col = "firebrick")
  mtext(paste("Utility loadings:", paste(utility_cols, collapse = ", ")),
        side = 3, line = 1, cex = 0.6, col = "steelblue")
}
