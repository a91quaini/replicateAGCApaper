# Reporting figures for the portfolio AGCA empirical section.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else {
  "."
}

find_compendium_root <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        grepl("^Package: replicateAGCApaper", readLines(file.path(path, "DESCRIPTION"), n = 1L))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the replicateAGCApaper root.", call. = FALSE)
    }
    path <- parent
  }
}

repo_dir <- find_compendium_root(script_dir)
source(file.path(repo_dir, "R", "GeodesicExtreme.R"))

output_dir <- file.path(repo_dir, "inst", "empirics", "results",
                        "portfolio_agca_reporting")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
manuscript_output_dir <- file.path(dirname(repo_dir), "full_paper", "R", "empirics",
                                   "empirics_output",
                                   "portfolio_agca_reporting")

stage_manuscript_outputs <- function() {
  manuscript_dir <- dirname(dirname(dirname(dirname(manuscript_output_dir))))
  if (!dir.exists(manuscript_dir)) {
    return(invisible(FALSE))
  }
  dir.create(manuscript_output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(output_dir, pattern = "\\.pdf$", full.names = TRUE)
  if (length(files) > 0L) {
    file.copy(files, manuscript_output_dir, overwrite = TRUE)
  }
  invisible(TRUE)
}

sample_start <- as.Date("1973-07-01")
tail_fraction <- 0.05
anchor_grid_fractions <- c(0.005, 0.01, 0.015, 0.025, 0.05, 0.075, 0.10)
bootstrap_reps <- 499L
set.seed(20260707)

dataset_specs <- list(
  ff_2x3_daily = list(
    label = "Fama-French size-based 2 x 3 sorts",
    short_label = "FF 2 x 3",
    data_file = file.path(repo_dir, "data", "empirics", "ff",
                          "ff_2x3_sorts_daily.rds"),
    group_column = "sort",
    group_label_column = "sort_label",
    max_rank = 15L,
    loading_ranks = 5L
  ),
  osap_daily_quintile_vw = list(
    label = "OSAP liquidity/trading quintile portfolios",
    short_label = "OSAP",
    data_file = file.path(repo_dir, "data", "empirics", "osap",
                          "osap_daily_quintile_vw.rds"),
    group_column = "signal",
    group_label_column = "signal_label",
    max_rank = 15L,
    loading_ranks = 5L
  )
)

rank_to_pareto <- function(x) {
  n <- length(x)
  r <- rank(x, ties.method = "average", na.last = "keep")
  (n + 1) / (n + 1 - r)
}

pareto_transform <- function(x) {
  transformed <- apply(x, 2L, rank_to_pareto)
  transformed <- as.matrix(transformed)
  colnames(transformed) <- colnames(x)
  rownames(transformed) <- rownames(x)
  transformed
}

select_top_k <- function(radius, k) {
  order(-radius, seq_along(radius))[seq_len(k)]
}

principal_anchor <- function(g) {
  moment <- crossprod(normalize_rows(g)) / nrow(g)
  eig <- eigen(moment, symmetric = TRUE)
  mu <- eig$vectors[, which.max(eig$values)]
  if (sum(mu) < 0) {
    mu <- -mu
  }
  if (any(mu <= 0)) {
    mu <- abs(mu)
  }
  unit_vector(mu, "principal anchor")
}

mean_anchor <- function(g) {
  unit_vector(colMeans(normalize_rows(g)), "mean anchor")
}

spherical_distance <- function(x, y) {
  acos(pmin(1, pmax(-1, sum(x * y))))
}

save_pdf <- function(file, plot_fun, width = 7.2, height = 4.8) {
  pdf(file, width = width, height = height, bg = "white")
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

load_dataset <- function(spec) {
  prepared <- readRDS(spec$data_file)
  loss_data <- prepared$complete_losses
  loss_data <- loss_data[loss_data$date >= sample_start, ]
  portfolio_columns <- prepared$portfolio_columns
  loss_matrix <- as.matrix(loss_data[, portfolio_columns, drop = FALSE])
  storage.mode(loss_matrix) <- "double"
  rownames(loss_matrix) <- format(loss_data$date, "%Y-%m-%d")

  pareto <- pareto_transform(loss_matrix)
  g <- normalize_rows(pareto)
  radius <- row_norms(pareto)
  k <- as.integer(round(tail_fraction * nrow(g)))
  index <- select_top_k(radius, k)

  list(
    prepared = prepared,
    loss_data = loss_data,
    g = g,
    radius = radius,
    index = index,
    g_extreme = g[index, , drop = FALSE],
    k = k
  )
}

bootstrap_variation <- function(g_extreme, mu, max_rank, reps) {
  n <- nrow(g_extreme)
  out <- matrix(NA_real_, reps, max_rank + 1L)
  colnames(out) <- paste0("rank_", 0L:max_rank)
  for (b in seq_len(reps)) {
    sample_index <- sample.int(n, n, replace = TRUE)
    fit_b <- agca_fit(g_extreme[sample_index, , drop = FALSE], mu = mu)
    rs_b <- agca_rank_summary(fit_b)
    out[b, ] <- rs_b$variation_explained[seq_len(max_rank + 1L)]
  }
  out
}

bootstrap_loading <- function(g_extreme, mu, target_loading, reps,
                              component = 1L) {
  n <- nrow(g_extreme)
  out <- matrix(NA_real_, reps, length(target_loading))
  for (b in seq_len(reps)) {
    sample_index <- sample.int(n, n, replace = TRUE)
    fit_b <- agca_fit(g_extreme[sample_index, , drop = FALSE], mu = mu)
    loading_b <- fit_b$loadings[, component]
    if (sum(loading_b * target_loading) < 0) {
      loading_b <- -loading_b
    }
    out[b, ] <- loading_b
  }
  out
}

summarize_bootstrap <- function(boot) {
  data.frame(
    rank = seq_len(ncol(boot)) - 1L,
    lower = apply(boot, 2L, quantile, probs = 0.025, na.rm = TRUE),
    median = apply(boot, 2L, quantile, probs = 0.5, na.rm = TRUE),
    upper = apply(boot, 2L, quantile, probs = 0.975, na.rm = TRUE)
  )
}

fit_summaries <- lapply(names(dataset_specs), function(dataset_name) {
  spec <- dataset_specs[[dataset_name]]
  data <- load_dataset(spec)
  mu <- canonical_anchor(ncol(data$g_extreme))
  fit <- agca_fit(data$g_extreme, mu = mu)
  rank_summary <- agca_rank_summary(fit)
  boot <- bootstrap_variation(data$g_extreme, mu, max_rank = spec$max_rank,
                              reps = bootstrap_reps)
  boot_summary <- summarize_bootstrap(boot)
  loading_boot <- bootstrap_loading(data$g_extreme, mu, fit$loadings[, 1L],
                                    reps = bootstrap_reps)

  anchors <- list(
    canonical = mu,
    principal = principal_anchor(data$g_extreme),
    mean = mean_anchor(data$g_extreme)
  )
  anchor_summary <- do.call(rbind, lapply(names(anchors), function(anchor) {
    rs <- agca_rank_summary(agca_fit(data$g_extreme, mu = anchors[[anchor]]))
    data.frame(
      dataset = dataset_name,
      anchor = anchor,
      rank = rs$rank,
      residual_risk = rs$residual_risk,
      variation_explained = rs$variation_explained
    )
  }))

  k_values <- unique(as.integer(round(anchor_grid_fractions * nrow(data$g))))
  k_values <- k_values[k_values > ncol(data$g) + 2L & k_values < nrow(data$g)]
  anchor_distance <- do.call(rbind, lapply(k_values, function(k) {
    idx <- select_top_k(data$radius, k)
    g_k <- data$g[idx, , drop = FALSE]
    mu_k <- canonical_anchor(ncol(g_k))
    anchors_k <- list(
      principal = principal_anchor(g_k),
      mean = mean_anchor(g_k)
    )
    do.call(rbind, lapply(names(anchors_k), function(anchor) {
      data.frame(
        dataset = dataset_name,
        k = k,
        tail_fraction = k / nrow(data$g),
        anchor = anchor,
        distance_to_canonical = spherical_distance(mu_k, anchors_k[[anchor]]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  threshold_summary <- do.call(rbind, lapply(k_values, function(k) {
    idx <- select_top_k(data$radius, k)
    rs <- agca_rank_summary(agca_fit(data$g[idx, , drop = FALSE], mu = mu))
    data.frame(
      dataset = dataset_name,
      k = k,
      tail_fraction = k / nrow(data$g),
      rank = rs$rank,
      variation_explained = rs$variation_explained
    )
  }))

  list(
    data = data,
    fit = fit,
    rank_summary = data.frame(
      dataset = dataset_name,
      rank_summary
    ),
    boot_summary = data.frame(
      dataset = dataset_name,
      boot_summary
    ),
    loading_boot = loading_boot,
    anchor_summary = anchor_summary,
    anchor_distance = anchor_distance,
    threshold_summary = threshold_summary
  )
})
names(fit_summaries) <- names(dataset_specs)

rank_summary_all <- do.call(rbind, lapply(fit_summaries, `[[`, "rank_summary"))
boot_summary_all <- do.call(rbind, lapply(fit_summaries, `[[`, "boot_summary"))
anchor_summary_all <- do.call(rbind, lapply(fit_summaries, `[[`,
                                            "anchor_summary"))
anchor_distance_all <- do.call(rbind, lapply(fit_summaries, `[[`,
                                             "anchor_distance"))
threshold_summary_all <- do.call(rbind, lapply(fit_summaries, `[[`,
                                               "threshold_summary"))

write.csv(rank_summary_all, file.path(output_dir, "rank_summary_main.csv"),
          row.names = FALSE)
write.csv(boot_summary_all, file.path(output_dir,
                                      "rank_summary_bootstrap.csv"),
          row.names = FALSE)
write.csv(anchor_summary_all, file.path(output_dir,
                                        "anchor_sensitivity.csv"),
          row.names = FALSE)
write.csv(anchor_distance_all, file.path(output_dir,
                                         "anchor_distance_by_k.csv"),
          row.names = FALSE)
write.csv(threshold_summary_all, file.path(output_dir,
                                           "threshold_sensitivity.csv"),
          row.names = FALSE)

plot_bootstrap_spectra <- function(file) {
  save_pdf(file, function() {
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    for (dataset_name in names(dataset_specs)) {
      spec <- dataset_specs[[dataset_name]]
      rs <- rank_summary_all[rank_summary_all$dataset == dataset_name &
                               rank_summary_all$rank <= spec$max_rank, ]
      bs <- boot_summary_all[boot_summary_all$dataset == dataset_name &
                               boot_summary_all$rank <= spec$max_rank, ]
      plot(
        rs$rank, rs$variation_explained,
        type = "n",
        ylim = c(0, 1),
        xlab = "AGCA rank p",
        ylab = "Cumulative AVE",
        main = spec$short_label
      )
      polygon(c(bs$rank, rev(bs$rank)), c(bs$lower, rev(bs$upper)),
              border = NA, col = grDevices::adjustcolor("#1B6CA8", 0.18))
      lines(rs$rank, rs$variation_explained, type = "b", pch = 16,
            col = "#1B6CA8", lwd = 1.6)
      grid(col = "gray90")
    }
  }, width = 7.8, height = 2.55)
}

plot_ff_spectrum_functionals <- function(file) {
  tail_dir <- file.path(repo_dir, "inst", "empirics", "results",
                        "portfolio_tail_functionals_agca")
  functional_summary <- read.csv(file.path(tail_dir,
                                           "functional_rank_summary.csv"))
  functional_boot <- read.csv(file.path(tail_dir,
                                        "ff_functional_bootstrap_summary.csv"))

  dataset_name <- "ff_2x3_daily"
  max_rank <- dataset_specs[[dataset_name]]$max_rank
  rs <- rank_summary_all[rank_summary_all$dataset == dataset_name &
                           rank_summary_all$rank >= 1L &
                           rank_summary_all$rank <= max_rank, ]
  bs <- boot_summary_all[boot_summary_all$dataset == dataset_name &
                           boot_summary_all$rank >= 1L &
                           boot_summary_all$rank <= max_rank, ]
  fs <- functional_summary[functional_summary$dataset == dataset_name &
                             functional_summary$rank >= 1L &
                             functional_summary$rank <= max_rank, ]
  fb <- functional_boot[functional_boot$dataset == dataset_name &
                          functional_boot$rank >= 1L &
                          functional_boot$rank <= max_rank, ]

  save_pdf(file, function() {
    par(mfrow = c(1L, 2L), mar = c(3.7, 4.3, 2.0, 0.8),
        mgp = c(2.2, 0.65, 0), tcl = -0.25)

    plot(
      rs$rank, 100 * rs$variation_explained,
      type = "n",
      ylim = c(0, 100),
      xlab = "AGCA rank p",
      ylab = "Cumulative AVE (%)",
      main = "Anchored variation"
    )
    grid(col = "gray90")
    spectrum_col <- "#1B6CA8"
    spectrum_band_col <- "#BFD9EE"
    spectrum_band_border <- grDevices::adjustcolor(spectrum_col, 0.70)
    polygon(c(bs$rank, rev(bs$rank)), 100 * c(bs$lower, rev(bs$upper)),
            border = NA, col = spectrum_band_col)
    lines(bs$rank, 100 * bs$lower, col = spectrum_band_border, lwd = 0.65)
    lines(bs$rank, 100 * bs$upper, col = spectrum_band_border, lwd = 0.65)
    abline(h = c(80, 90), col = "gray70", lty = 3, lwd = 0.9)
    lines(rs$rank, 100 * rs$variation_explained, type = "b", pch = 16,
          col = spectrum_col, lwd = 1.7)
    box()

    err_max <- max(100 * c(fb$capped_upper, fb$var_upper,
                           fs$capped_relative_error_agca_mean,
                           fs$var_relative_error_agca_mean),
                   finite = TRUE)
    plot(
      fs$rank, 100 * fs$capped_relative_error_agca_mean,
      type = "n",
      ylim = c(0, 1.08 * err_max),
      xlab = "AGCA rank p",
      ylab = "Mean relative error (%)",
      main = "Portfolio tail summaries"
    )
    grid(col = "gray90")
    capped_col <- "#B23A48"
    var_col <- "#3B8C5A"
    capped_band_col <- "#EFA9B3"
    var_band_col <- "#A9D8BA"
    capped_band_border <- grDevices::adjustcolor(capped_col, 0.80)
    var_band_border <- grDevices::adjustcolor(var_col, 0.80)
    polygon(c(fb$rank, rev(fb$rank)),
            100 * c(fb$capped_lower, rev(fb$capped_upper)),
            border = NA, col = capped_band_col)
    polygon(c(fb$rank, rev(fb$rank)),
            100 * c(fb$var_lower, rev(fb$var_upper)),
            border = NA, col = var_band_col)
    lines(fb$rank, 100 * fb$capped_lower, col = capped_band_border, lwd = 0.65)
    lines(fb$rank, 100 * fb$capped_upper, col = capped_band_border, lwd = 0.65)
    lines(fb$rank, 100 * fb$var_lower, col = var_band_border, lwd = 0.65)
    lines(fb$rank, 100 * fb$var_upper, col = var_band_border, lwd = 0.65)
    lines(fs$rank, 100 * fs$capped_relative_error_agca_mean, type = "b",
          pch = 16, col = capped_col, lwd = 1.6)
    lines(fs$rank, 100 * fs$var_relative_error_agca_mean, type = "b",
          pch = 17, col = var_col, lwd = 1.6, lty = 2)
    legend("topright", legend = c("Capped portfolio loss", "Normalized VaR"),
           col = c(capped_col, var_col), pch = c(16, 17), lty = c(1, 2),
           lwd = 1.6, bty = "n", cex = 0.78)
    box()
  }, width = 7.8, height = 2.75)
}

plot_anchor_sensitivity <- function(file) {
  save_pdf(file, function() {
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    cols <- c(canonical = "#1B6CA8", principal = "#B23A48",
              mean = "#3B8C5A")
    for (dataset_name in names(dataset_specs)) {
      spec <- dataset_specs[[dataset_name]]
      z <- anchor_summary_all[anchor_summary_all$dataset == dataset_name &
                                anchor_summary_all$rank <= spec$max_rank, ]
      plot(
        NA,
        xlim = c(0, spec$max_rank),
        ylim = c(0, 1),
        xlab = "AGCA rank p",
        ylab = "Cumulative AVE",
        main = spec$short_label
      )
      grid(col = "gray90")
      for (anchor in names(cols)) {
        zz <- z[z$anchor == anchor, ]
        lines(zz$rank, zz$variation_explained, type = "b", pch = 16,
              col = cols[anchor], lwd = 1.4)
      }
      legend("bottomright", legend = names(cols), col = cols, lty = 1,
             pch = 16, bty = "n", cex = 0.75)
    }
  }, width = 7.8, height = 2.55)
}

plot_anchor_residual_risk <- function(file) {
  save_pdf(file, function() {
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    cols <- c(canonical = "#1B6CA8", principal = "#B23A48",
              mean = "#3B8C5A")
    for (dataset_name in names(dataset_specs)) {
      spec <- dataset_specs[[dataset_name]]
      z <- anchor_summary_all[anchor_summary_all$dataset == dataset_name &
                                anchor_summary_all$rank <= spec$max_rank, ]
      ylim <- range(z$residual_risk, finite = TRUE)
      plot(
        NA,
        xlim = c(0, spec$max_rank),
        ylim = ylim,
        xlab = "AGCA rank p",
        ylab = "Residual risk",
        main = spec$short_label
      )
      grid(col = "gray90")
      for (anchor in names(cols)) {
        zz <- z[z$anchor == anchor, ]
        lines(zz$rank, zz$residual_risk, type = "b", pch = 16,
              col = cols[anchor], lwd = 1.4)
      }
      legend("topright", legend = names(cols), col = cols, lty = 1,
             pch = 16, bty = "n", cex = 0.75)
    }
  }, width = 7.8, height = 2.55)
}

plot_anchor_distance <- function(file) {
  save_pdf(file, function() {
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    cols <- c(principal = "#B23A48", mean = "#3B8C5A")
    ylim <- c(0, max(anchor_distance_all$distance_to_canonical,
                     na.rm = TRUE))
    for (dataset_name in names(dataset_specs)) {
      spec <- dataset_specs[[dataset_name]]
      z <- anchor_distance_all[anchor_distance_all$dataset == dataset_name, ]
      plot(
        NA,
        xlim = range(z$k),
        ylim = ylim,
        xlab = "Selected extremes k",
        ylab = "Distance to canonical (radians)",
        main = spec$short_label
      )
      grid(col = "gray90")
      for (anchor in names(cols)) {
        zz <- z[z$anchor == anchor, ]
        zz <- zz[order(zz$k), ]
        lines(zz$k, zz$distance_to_canonical, type = "b", pch = 16,
              col = cols[anchor], lwd = 1.4)
      }
      legend("topright", legend = names(cols), col = cols, lty = 1,
             pch = 16, bty = "n", cex = 0.75)
    }
  }, width = 7.8, height = 2.55)
}

plot_threshold_sensitivity <- function(file) {
  save_pdf(file, function() {
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    for (dataset_name in names(dataset_specs)) {
      spec <- dataset_specs[[dataset_name]]
      z <- threshold_summary_all[
        threshold_summary_all$dataset == dataset_name &
          threshold_summary_all$rank <= spec$max_rank,
      ]
      k_values <- sort(unique(z$k))
      cols <- grDevices::hcl.colors(length(k_values), "Dark 3")
      plot(
        NA,
        xlim = c(0, spec$max_rank),
        ylim = c(0, 1),
        xlab = "AGCA rank p",
        ylab = "Cumulative AVE",
        main = spec$short_label
      )
      grid(col = "gray90")
      for (i in seq_along(k_values)) {
        zz <- z[z$k == k_values[i], ]
        lines(zz$rank, zz$variation_explained, type = "b", pch = 16,
              col = cols[i], lwd = 1.25)
      }
      legend("bottomright", legend = paste0("k=", k_values), col = cols,
             lty = 1, pch = 16, bty = "n", cex = 0.65)
    }
  }, width = 7.8, height = 2.55)
}

asset_labels <- function(metadata, portfolio_columns, group_column) {
  z <- metadata[match(portfolio_columns, metadata$variable), ]
  if (group_column == "sort") {
    short_group <- c(
      size_bm = "B/M",
      size_op = "OP",
      size_inv = "Inv",
      size_mom = "Mom"
    )
    paste0(short_group[z$sort], "\n", gsub("SMALL ", "S ",
                                           gsub("BIG ", "B ", z$portfolio)))
  } else {
    paste0(z[[group_column]], "\n", gsub("port0", "Q", z$portfolio))
  }
}

group_boundaries <- function(metadata, portfolio_columns, group_column) {
  z <- metadata[match(portfolio_columns, metadata$variable), ]
  rle_groups <- rle(z[[group_column]])
  cumsum(rle_groups$lengths) + 0.5
}

loading_boot_summary_all <- do.call(rbind, lapply(names(dataset_specs),
                                                  function(dataset_name) {
  spec <- dataset_specs[[dataset_name]]
  data <- fit_summaries[[dataset_name]]$data
  fit <- fit_summaries[[dataset_name]]$fit
  metadata <- data$prepared$portfolios
  portfolio_columns <- data$prepared$portfolio_columns
  labels <- asset_labels(metadata, portfolio_columns, spec$group_column)
  boot <- fit_summaries[[dataset_name]]$loading_boot
  data.frame(
    dataset = dataset_name,
    asset_index = seq_along(portfolio_columns),
    asset = portfolio_columns,
    asset_label = labels,
    loading = fit$loadings[, 1L],
    lower = apply(boot, 2L, quantile, probs = 0.025, na.rm = TRUE),
    median = apply(boot, 2L, quantile, probs = 0.5, na.rm = TRUE),
    upper = apply(boot, 2L, quantile, probs = 0.975, na.rm = TRUE)
  )
}))

write.csv(loading_boot_summary_all,
          file.path(output_dir, "loading_agc1_bootstrap.csv"),
          row.names = FALSE)

plot_loading_heatmap <- function(dataset_name, file) {
  spec <- dataset_specs[[dataset_name]]
  data <- fit_summaries[[dataset_name]]$data
  fit <- fit_summaries[[dataset_name]]$fit
  metadata <- data$prepared$portfolios
  portfolio_columns <- data$prepared$portfolio_columns
  p <- spec$loading_ranks
  loadings <- t(fit$loadings[, seq_len(p), drop = FALSE])
  rownames(loadings) <- paste0("AGC", seq_len(p))
  colnames(loadings) <- asset_labels(metadata, portfolio_columns,
                                     spec$group_column)
  max_abs <- max(abs(loadings))
  palette <- colorRampPalette(c("#B23A48", "white", "#1B6CA8"))(101)
  breaks <- seq(-max_abs, max_abs, length.out = length(palette) + 1L)

  save_pdf(file, function() {
    par(mar = c(4.1, 4.6, 2.6, 4.8))
    image(
      x = seq_len(ncol(loadings)),
      y = seq_len(nrow(loadings)),
      z = t(loadings[nrow(loadings):1L, , drop = FALSE]),
      col = palette,
      breaks = breaks,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = spec$short_label
    )
    axis(1, at = seq_len(ncol(loadings)), labels = FALSE)
    par(xpd = NA)
    text(seq_len(ncol(loadings)), par("usr")[3] - 0.10,
         labels = gsub("\n", " ", colnames(loadings)), srt = 45,
         adj = c(1, 1), cex = 0.52)
    par(xpd = FALSE)
    axis(2, at = seq_len(nrow(loadings)),
         labels = rev(rownames(loadings)), las = 2, cex.axis = 0.85)
    boundary <- group_boundaries(metadata, portfolio_columns,
                                 spec$group_column)
    boundary <- boundary[boundary < ncol(loadings) + 0.5]
    abline(v = boundary, col = "gray35", lwd = 0.8)
    box()
    par(xpd = TRUE)
    legend_x <- ncol(loadings) + 1.8
    legend_y <- seq(1, nrow(loadings), length.out = length(palette))
    rect(legend_x, legend_y[-length(legend_y)], legend_x + 0.55,
         legend_y[-1], col = palette, border = NA)
    text(legend_x + 0.72, min(legend_y), sprintf("%.2f", -max_abs),
         adj = 0, cex = 0.65)
    text(legend_x + 0.72, mean(range(legend_y)), "0", adj = 0, cex = 0.65)
    text(legend_x + 0.72, max(legend_y), sprintf("%.2f", max_abs),
         adj = 0, cex = 0.65)
    par(xpd = FALSE)
  }, width = 7.4, height = 2.8)
}

plot_loading_bootstrap <- function(file) {
  save_pdf(file, function() {
    par(mfrow = c(2L, 1L), mar = c(4.7, 4.5, 2.2, 0.8))
    for (dataset_name in names(dataset_specs)) {
      spec <- dataset_specs[[dataset_name]]
      data <- fit_summaries[[dataset_name]]$data
      metadata <- data$prepared$portfolios
      portfolio_columns <- data$prepared$portfolio_columns
      z <- loading_boot_summary_all[
        loading_boot_summary_all$dataset == dataset_name,
      ]
      z <- z[order(z$asset_index), ]
      ylim <- range(c(z$lower, z$upper), finite = TRUE)
      ylim <- ylim + c(-0.08, 0.08) * diff(ylim)
      cols <- ifelse(z$loading >= 0, "#1B6CA8", "#B23A48")
      mids <- barplot(
        z$loading,
        names.arg = rep("", nrow(z)),
        col = grDevices::adjustcolor(cols, 0.82),
        border = "gray35",
        ylim = ylim,
        main = spec$short_label,
        ylab = "AGC1 loading"
      )
      arrows(mids, z$lower, mids, z$upper, angle = 90, code = 3,
             length = 0.025, col = "gray20", lwd = 0.8)
      abline(h = 0, col = "gray30", lwd = 0.8)
      boundary <- group_boundaries(metadata, portfolio_columns,
                                   spec$group_column)
      boundary <- boundary[boundary < nrow(z) + 0.5]
      if (length(boundary) > 0L) {
        left_index <- as.integer(boundary - 0.5)
        boundary_x <- (mids[left_index] + mids[left_index + 1L]) / 2
        abline(v = boundary_x, col = "gray55", lwd = 0.7)
      }
      par(xpd = NA)
      text(mids, par("usr")[3] - 0.025 * diff(par("usr")[3:4]),
           labels = gsub("\n", " ", z$asset_label), srt = 45,
           adj = c(1, 1), cex = 0.48)
      par(xpd = FALSE)
    }
  }, width = 7.8, height = 5.4)
}

plot_bootstrap_spectra(file.path(output_dir,
                                 "variation_explained_bootstrap.pdf"))
plot_ff_spectrum_functionals(file.path(output_dir,
                                       "ff_spectrum_functional_summary.pdf"))
plot_anchor_sensitivity(file.path(output_dir,
                                  "variation_explained_anchor_sensitivity.pdf"))
plot_anchor_residual_risk(file.path(output_dir,
                                    "residual_risk_anchor_sensitivity.pdf"))
plot_anchor_distance(file.path(output_dir, "anchor_distance_by_k.pdf"))
plot_threshold_sensitivity(file.path(output_dir,
                                     "variation_explained_threshold_sensitivity.pdf"))
plot_loading_heatmap("ff_2x3_daily",
                     file.path(output_dir, "ff_loading_heatmap.pdf"))
plot_loading_heatmap("osap_daily_quintile_vw",
                     file.path(output_dir, "osap_loading_heatmap.pdf"))
plot_loading_bootstrap(file.path(output_dir,
                                 "loading_agc1_bootstrap_intervals.pdf"))
stage_manuscript_outputs()

cat("\nPortfolio AGCA reporting figures written to:\n  ", output_dir, "\n",
    sep = "")
