# Coexceedance comparison for both portfolio panels, plus Fama--French diagnostics.
# Use --comparison-only to regenerate the two-panel comparison without rerunning
# the existing Fama--French sensitivity analyses or overwriting their outputs.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else {
  "."
}

find_compendium_root <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  repeat {
    description <- file.path(path, "DESCRIPTION")
    if (file.exists(description) &&
        grepl("^Package: replicateAGCApaper",
              readLines(description, n = 1L))) {
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
                        "portfolio_coexceedances_agca")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paper_figure_dir <- file.path(dirname(repo_dir), "jrssb_paper", "figures")

stage_paper_figures <- function() {
  if (!dir.exists(dirname(paper_figure_dir))) {
    return(invisible(FALSE))
  }
  dir.create(paper_figure_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(output_dir, pattern = "^fig_emp_.*\\.pdf$",
                      full.names = TRUE)
  # Stage only assets referenced by the paper; retain other diagnostics in the
  # result directory without restoring retired figures to the manuscript.
  tex_files <- list.files(dirname(paper_figure_dir), pattern = "\\.tex$",
                          full.names = TRUE)
  tex_source <- paste(unlist(lapply(tex_files, readLines, warn = FALSE)),
                      collapse = "\n")
  files <- files[vapply(basename(files), grepl, logical(1L),
                        x = tex_source, fixed = TRUE)]
  if (length(files) > 0L) {
    file.copy(files, paper_figure_dir, overwrite = TRUE)
  }
  invisible(TRUE)
}

sample_start <- as.Date("1973-07-01")
tail_fraction <- 0.05
ave_target <- 0.90
max_plot_rank <- 15L
threshold_fractions <- c(0.005, 0.01, 0.015, 0.025, 0.05, 0.075, 0.10)
severity_multipliers <- c(1, 2, 4)
bootstrap_reps <- 1000L
bootstrap_seed <- 20260707L

rank_to_pareto <- function(x) {
  n <- length(x)
  r <- rank(x, ties.method = "average", na.last = "keep")
  (n + 1) / (n + 1 - r)
}

pareto_transform <- function(x) {
  transformed <- as.matrix(apply(x, 2L, rank_to_pareto))
  colnames(transformed) <- colnames(x)
  rownames(transformed) <- rownames(x)
  transformed
}

select_top_k <- function(radius, k) {
  if (length(k) != 1L || k < 1L || k > length(radius) ||
      k != as.integer(k)) {
    stop("k must be an integer between 1 and the sample size.",
         call. = FALSE)
  }
  order(-radius, seq_along(radius))[seq_len(k)]
}

positive_post_project <- function(g) {
  normalize_rows(pmax(g, 0))
}

row_order_statistics <- function(g) {
  g <- as.matrix(g)
  t(vapply(seq_len(nrow(g)), function(i) {
    sort(g[i, ], decreasing = TRUE)
  }, numeric(ncol(g))))
}

fitted_coexceedance_curve <- function(g, tail_mass, x = 1,
                                      positive_part = FALSE) {
  order_statistics <- row_order_statistics(g)
  if (positive_part) {
    order_statistics <- pmax(order_statistics, 0)
  }
  tail_mass * colMeans(order_statistics) / x
}

empirical_marginal_quantile <- function(loss_matrix, q) {
  vapply(seq_len(ncol(loss_matrix)), function(j) {
    # Type 1 is the inverse empirical distribution function. Strict
    # inequalities below exclude observations tied at the quantile threshold.
    quantile(loss_matrix[, j], probs = q, type = 1, names = FALSE)
  }, numeric(1L))
}

observed_coexceedances <- function(loss_matrix, q) {
  thresholds <- empirical_marginal_quantile(loss_matrix, q)
  breaches <- sweep(loss_matrix, 2L, thresholds, FUN = ">")
  counts <- rowSums(breaches)
  curve <- vapply(seq_len(ncol(loss_matrix)), function(m) {
    mean(counts >= m)
  }, numeric(1L))
  list(
    thresholds = thresholds,
    breaches = breaches,
    counts = counts,
    curve = curve
  )
}

curve_errors <- function(estimate, target) {
  difference <- estimate - target
  c(
    MAE = mean(abs(difference)),
    MAX = max(abs(difference))
  )
}

analyze_matrix <- function(loss_matrix, dates, fraction = tail_fraction) {
  x <- pareto_transform(loss_matrix)
  radius <- row_norms(x)
  n <- nrow(x)
  k <- as.integer(round(fraction * n))
  index <- select_top_k(radius, k)
  radial_threshold <- min(radius[index])
  q_star <- 1 - 1 / radial_threshold
  g_extreme <- normalize_rows(x[index, , drop = FALSE])

  fit <- agca_fit(g_extreme, mu = canonical_anchor(ncol(g_extreme)))
  geometric_rank_summary <- agca_rank_summary(fit)
  p_star <- min(geometric_rank_summary$rank[
    geometric_rank_summary$variation_explained >= ave_target
  ])
  tail_mass <- k / n
  full_curve <- fitted_coexceedance_curve(g_extreme, tail_mass)

  rank_rows <- lapply(geometric_rank_summary$rank, function(p) {
    raw <- agca_reconstruct(fit, p = p)
    corrected <- positive_post_project(raw)
    raw_curve <- fitted_coexceedance_curve(
      raw, tail_mass, positive_part = TRUE
    )
    corrected_curve <- fitted_coexceedance_curve(corrected, tail_mass)
    raw_error <- curve_errors(raw_curve, full_curve)
    corrected_error <- curve_errors(corrected_curve, full_curve)
    needs_correction <- rowSums(raw < 0) > 0
    correction_distance <- row_norms(raw - corrected)

    data.frame(
      rank = p,
      variation_explained = geometric_rank_summary$variation_explained[
        geometric_rank_summary$rank == p
      ],
      residual_risk = geometric_rank_summary$residual_risk[
        geometric_rank_summary$rank == p
      ],
      mae_raw = unname(raw_error["MAE"]),
      max_raw = unname(raw_error["MAX"]),
      mae_corrected = unname(corrected_error["MAE"]),
      max_corrected = unname(corrected_error["MAX"]),
      correction_proportion = mean(needs_correction),
      mean_negative_coordinates = mean(rowSums(raw < 0)),
      mean_correction_distance = mean(correction_distance),
      mean_squared_correction_distance = mean(correction_distance^2),
      stringsAsFactors = FALSE
    )
  })
  rank_summary <- do.call(rbind, rank_rows)

  raw_star <- agca_reconstruct(fit, p = p_star)
  corrected_star <- positive_post_project(raw_star)
  raw_curve <- fitted_coexceedance_curve(
    raw_star, tail_mass, positive_part = TRUE
  )
  corrected_curve <- fitted_coexceedance_curve(corrected_star, tail_mass)
  observed <- observed_coexceedances(loss_matrix, q_star)

  list(
    n = n,
    d = ncol(x),
    k = k,
    dates = range(dates),
    x = x,
    radius = radius,
    index = index,
    radial_threshold = radial_threshold,
    q_star = q_star,
    g_extreme = g_extreme,
    fit = fit,
    rank_summary = rank_summary,
    p_star = p_star,
    tail_mass = tail_mass,
    observed_thresholds = observed$thresholds,
    observed_breaches = observed$breaches,
    observed_counts = observed$counts,
    observed_curve = observed$curve,
    full_curve = full_curve,
    raw_curve = raw_curve,
    corrected_curve = corrected_curve,
    raw_star = raw_star,
    corrected_star = corrected_star
  )
}

bootstrap_ave <- function(g, mu, max_rank, reps, seed) {
  set.seed(seed)
  n <- nrow(g)
  estimates <- matrix(NA_real_, nrow = reps, ncol = max_rank)
  for (b in seq_len(reps)) {
    index <- sample.int(n, n, replace = TRUE)
    fit_b <- agca_fit(g[index, , drop = FALSE], mu = mu)
    estimates[b, ] <- agca_variation_explained(fit_b)[seq_len(max_rank)]
  }
  data.frame(
    rank = seq_len(max_rank),
    lower = apply(estimates, 2L, quantile, probs = 0.025),
    median = apply(estimates, 2L, quantile, probs = 0.5),
    upper = apply(estimates, 2L, quantile, probs = 0.975)
  )
}

bootstrap_coexceedance_curves <- function(g, mu, tail_mass, p, reps, seed) {
  set.seed(seed)
  n <- nrow(g)
  d <- ncol(g)
  full_estimates <- matrix(NA_real_, nrow = reps, ncol = d)
  agca_estimates <- matrix(NA_real_, nrow = reps, ncol = d)

  for (b in seq_len(reps)) {
    index <- sample.int(n, n, replace = TRUE)
    g_b <- g[index, , drop = FALSE]
    fit_b <- agca_fit(g_b, mu = mu)
    reconstructed_b <- positive_post_project(agca_reconstruct(fit_b, p = p))
    full_estimates[b, ] <- fitted_coexceedance_curve(g_b, tail_mass)
    agca_estimates[b, ] <- fitted_coexceedance_curve(
      reconstructed_b, tail_mass
    )
  }

  summarize <- function(estimates, method) {
    data.frame(
      method = method,
      m = seq_len(d),
      lower = apply(estimates, 2L, quantile, probs = 0.025),
      median = apply(estimates, 2L, quantile, probs = 0.5),
      upper = apply(estimates, 2L, quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    summarize(full_estimates, "full"),
    summarize(agca_estimates, "agca"),
    summarize(agca_estimates - full_estimates, "agca_minus_full")
  )
}

survival_to_pmf <- function(survival) {
  c(
    1 - survival[1L],
    survival[-length(survival)] - survival[-1L],
    survival[length(survival)]
  )
}

save_pdf <- function(file, plot_fun, width, height, pointsize = 12) {
  pdf(file, width = width, height = height, bg = "white", version = "1.4",
      useDingbats = FALSE, pointsize = pointsize)
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

curve_colors <- c(
  observed = "#222222",
  full = "#1B6CA8",
  corrected = "#B23A48"
)
curve_band_colors <- c(
  full = grDevices::adjustcolor(curve_colors["full"], 0.25),
  agca = grDevices::adjustcolor(curve_colors["corrected"], 0.25)
)
curve_lty <- c(observed = 3L, full = 1L, corrected = 2L)
curve_pch <- c(observed = 1L, full = 16L, corrected = 17L)

draw_coexceedance_panel <- function(z, main = "", ylim = NULL,
                                    show_legend = TRUE, intervals = NULL,
                                    reference_label = "Full", count_ticks = NULL) {
  if (is.null(ylim)) {
    upper <- c(z$observed, z$full, z$agca)
    if (!is.null(intervals)) {
      upper <- c(upper, intervals$upper)
    }
    ylim <- c(0, 1.04 * max(upper))
  }
  plot(z$m, z$observed, type = "n",
       ylim = ylim, xlab = "Number of exceedances m",
       ylab = expression(hat(J)[m]), main = main, xaxt = "n")
  if (is.null(count_ticks)) {
    count_ticks <- unique(c(1L, pretty(z$m)))
  }
  axis(1, at = count_ticks)
  grid(col = "gray90")
  if (!is.null(intervals)) {
    for (method in c("full", "agca")) {
      band <- intervals[intervals$method == method, , drop = FALSE]
      band <- band[order(band$m), , drop = FALSE]
      polygon(
        c(band$m, rev(band$m)),
        c(band$lower, rev(band$upper)),
        border = NA,
        col = curve_band_colors[method]
      )
      band_line_color <- if (method == "full") {
        grDevices::adjustcolor(curve_colors["full"], 0.55)
      } else {
        grDevices::adjustcolor(curve_colors["corrected"], 0.55)
      }
      lines(band$m, band$lower, col = band_line_color, lwd = 0.55)
      lines(band$m, band$upper, col = band_line_color, lwd = 0.55)
      segments(
        band$m, band$lower, band$m, band$upper,
        col = band_line_color, lwd = 0.9
      )
      segments(
        band$m - 0.055, band$lower, band$m + 0.055, band$lower,
        col = band_line_color, lwd = 0.7
      )
      segments(
        band$m - 0.055, band$upper, band$m + 0.055, band$upper,
        col = band_line_color, lwd = 0.7
      )
    }
  }
  lines(z$m, z$observed, type = "b", pch = curve_pch["observed"],
        lty = curve_lty["observed"], col = curve_colors["observed"],
        cex = 0.75, lwd = 1.2)
  lines(z$m, z$full, type = "b", pch = curve_pch["full"],
        lty = curve_lty["full"], col = curve_colors["full"],
        cex = 0.75, lwd = 1.5)
  lines(z$m, z$agca, type = "b", pch = curve_pch["corrected"],
        lty = curve_lty["corrected"], col = curve_colors["corrected"],
        cex = 0.75, lwd = 1.5)
  box()
  if (show_legend) {
    legend("topright", legend = c("Observed", reference_label, "AGCA"),
           col = curve_colors, lty = curve_lty, pch = curve_pch,
           lwd = c(1.2, 1.5, 1.5), bty = "n", cex = 0.74)
  }
}

if (!"--comparison-only" %in% commandArgs(TRUE)) {

data_file <- file.path(repo_dir, "data", "empirics", "ff",
                       "ff_2x3_sorts_daily.rds")
prepared <- readRDS(data_file)
data <- prepared$complete_losses
data <- data[data$date >= sample_start, ]
portfolio_columns <- prepared$portfolio_columns
loss_matrix <- as.matrix(data[, portfolio_columns, drop = FALSE])
storage.mode(loss_matrix) <- "double"
rownames(loss_matrix) <- format(data$date, "%Y-%m-%d")

main_fit <- analyze_matrix(loss_matrix, data$date)
p_star_row <- main_fit$rank_summary[
  main_fit$rank_summary$rank == main_fit$p_star,
]

main_curves <- data.frame(
  m = seq_len(main_fit$d),
  observed = main_fit$observed_curve,
  full = main_fit$full_curve,
  agca = main_fit$corrected_curve,
  agca_raw = main_fit$raw_curve,
  observed_minus_full = main_fit$observed_curve - main_fit$full_curve,
  agca_minus_full = main_fit$corrected_curve - main_fit$full_curve,
  agca_raw_minus_full = main_fit$raw_curve - main_fit$full_curve
)

main_error_summary <- rbind(
  data.frame(comparison = "observed_vs_full",
             t(curve_errors(main_fit$observed_curve,
                            main_fit$full_curve))),
  data.frame(comparison = "raw_agca_vs_full",
             t(curve_errors(main_fit$raw_curve,
                            main_fit$full_curve))),
  data.frame(comparison = "corrected_agca_vs_full",
             t(curve_errors(main_fit$corrected_curve,
                            main_fit$full_curve)))
)
names(main_error_summary)[2:3] <- c("mae", "max")

needs_correction <- rowSums(main_fit$raw_star < 0) > 0
correction_distance <- row_norms(
  main_fit$raw_star - main_fit$corrected_star
)
marginal_breach_totals <- colSums(main_fit$observed_breaches)
positivity_summary <- data.frame(
  rank = main_fit$p_star,
  n_directions = main_fit$k,
  n_corrected = sum(needs_correction),
  correction_proportion = mean(needs_correction),
  mean_negative_coordinates = mean(rowSums(main_fit$raw_star < 0)),
  mean_negative_coordinates_if_corrected =
    mean(rowSums(main_fit$raw_star[needs_correction, , drop = FALSE] < 0)),
  mean_correction_distance = mean(correction_distance),
  mean_correction_distance_if_corrected =
    mean(correction_distance[needs_correction]),
  mean_squared_correction_distance = mean(correction_distance^2),
  rms_correction_distance = sqrt(mean(correction_distance^2)),
  curve_mae_raw_vs_corrected =
    mean(abs(main_fit$raw_curve - main_fit$corrected_curve)),
  curve_max_raw_vs_corrected =
    max(abs(main_fit$raw_curve - main_fit$corrected_curve)),
  stringsAsFactors = FALSE
)

main_summary <- data.frame(
  n = main_fit$n,
  d = main_fit$d,
  k = main_fit$k,
  start_date = main_fit$dates[1L],
  end_date = main_fit$dates[2L],
  radial_threshold = main_fit$radial_threshold,
  q_star = main_fit$q_star,
  p_star = main_fit$p_star,
  ave_p_star = p_star_row$variation_explained,
  mae_p_star = p_star_row$mae_corrected,
  max_p_star = p_star_row$max_corrected,
  observed_full_mae = main_error_summary$mae[
    main_error_summary$comparison == "observed_vs_full"
  ],
  observed_full_max = main_error_summary$max[
    main_error_summary$comparison == "observed_vs_full"
  ],
  min_marginal_breaches = min(marginal_breach_totals),
  max_marginal_breaches = max(marginal_breach_totals),
  j6_observed = main_fit$observed_curve[6L],
  j6_full = main_fit$full_curve[6L],
  j6_agca = main_fit$corrected_curve[6L],
  j12_observed = main_fit$observed_curve[12L],
  j12_full = main_fit$full_curve[12L],
  j12_agca = main_fit$corrected_curve[12L],
  stringsAsFactors = FALSE
)

severity_curves <- do.call(rbind, lapply(severity_multipliers, function(x) {
  q_x <- 1 - 1 / (main_fit$radial_threshold * x)
  data.frame(
    x = x,
    q_x = q_x,
    m = seq_len(main_fit$d),
    observed = observed_coexceedances(loss_matrix, q_x)$curve,
    full = main_fit$full_curve / x,
    agca = main_fit$corrected_curve / x
  )
}))

severity_error_summary <- do.call(rbind, lapply(
  split(severity_curves, severity_curves$x), function(z) {
    observed_error <- curve_errors(z$observed, z$full)
    agca_error <- curve_errors(z$agca, z$full)
    data.frame(
      x = z$x[1L],
      q_x = z$q_x[1L],
      observed_full_mae = unname(observed_error["MAE"]),
      observed_full_max = unname(observed_error["MAX"]),
      agca_full_mae = unname(agca_error["MAE"]),
      agca_full_max = unname(agca_error["MAX"])
    )
  }
))

threshold_analyses <- lapply(threshold_fractions, function(fraction) {
  analyze_matrix(loss_matrix, data$date, fraction = fraction)
})

threshold_rank_summary <- do.call(rbind, Map(
  function(fraction, analysis) {
    z <- analysis$rank_summary
    data.frame(
      target_tail_fraction = fraction,
      k = analysis$k,
      tail_mass = analysis$tail_mass,
      radial_threshold = analysis$radial_threshold,
      q_star = analysis$q_star,
      p_star = analysis$p_star,
      z,
      scaled_mae_corrected = z$mae_corrected / analysis$tail_mass,
      scaled_max_corrected = z$max_corrected / analysis$tail_mass
    )
  },
  threshold_fractions, threshold_analyses
))

threshold_selected_summary <- do.call(rbind, lapply(
  split(threshold_rank_summary,
        threshold_rank_summary$target_tail_fraction), function(z) {
    z[z$rank == z$p_star[1L], , drop = FALSE]
  }
))
rownames(threshold_selected_summary) <- NULL

block_ids <- unique(prepared$portfolios$sort)
block_short_labels <- c(
  size_bm = "Book-to-market",
  size_op = "Profitability",
  size_inv = "Investment",
  size_mom = "Momentum"
)

block_curves <- do.call(rbind, lapply(block_ids, function(block) {
  block_columns <- prepared$portfolios$variable[
    prepared$portfolios$sort == block
  ]
  column_index <- match(block_columns, portfolio_columns)
  observed_counts <- rowSums(
    main_fit$observed_breaches[, column_index, drop = FALSE]
  )
  block_d <- length(column_index)
  data.frame(
    block = block,
    block_label = unname(block_short_labels[block]),
    m = seq_len(block_d),
    observed = vapply(seq_len(block_d), function(m) {
      mean(observed_counts >= m)
    }, numeric(1L)),
    full = fitted_coexceedance_curve(
      main_fit$g_extreme[, column_index, drop = FALSE], main_fit$tail_mass
    ),
    agca = fitted_coexceedance_curve(
      main_fit$corrected_star[, column_index, drop = FALSE],
      main_fit$tail_mass
    )
  )
}))

block_pmf <- do.call(rbind, lapply(split(block_curves,
                                         block_curves$block), function(z) {
  data.frame(
    block = z$block[1L],
    block_label = z$block_label[1L],
    count = 0L:nrow(z),
    observed = survival_to_pmf(z$observed),
    full = survival_to_pmf(z$full),
    agca = survival_to_pmf(z$agca)
  )
}))

block_error_summary <- do.call(rbind, lapply(
  split(block_curves, block_curves$block), function(z) {
    observed_error <- curve_errors(z$observed, z$full)
    agca_error <- curve_errors(z$agca, z$full)
    data.frame(
      block = z$block[1L],
      block_label = z$block_label[1L],
      observed_full_mae = unname(observed_error["MAE"]),
      observed_full_max = unname(observed_error["MAX"]),
      agca_full_mae = unname(agca_error["MAE"]),
      agca_full_max = unname(agca_error["MAX"])
    )
  }
))

split_point <- floor(nrow(loss_matrix) / 2)
split_indices <- list(
  early = seq_len(split_point),
  late = (split_point + 1L):nrow(loss_matrix)
)
split_analyses <- lapply(split_indices, function(index) {
  analyze_matrix(loss_matrix[index, , drop = FALSE], data$date[index])
})

split_summary <- do.call(rbind, lapply(names(split_analyses), function(name) {
  analysis <- split_analyses[[name]]
  selected <- analysis$rank_summary[
    analysis$rank_summary$rank == analysis$p_star,
  ]
  observed_error <- curve_errors(
    analysis$observed_curve, analysis$full_curve
  )
  data.frame(
    sample = name,
    n = analysis$n,
    k = analysis$k,
    start_date = analysis$dates[1L],
    end_date = analysis$dates[2L],
    radial_threshold = analysis$radial_threshold,
    q_star = analysis$q_star,
    p_star = analysis$p_star,
    ave_p_star = selected$variation_explained,
    mae_p_star = selected$mae_corrected,
    max_p_star = selected$max_corrected,
    observed_full_mae = unname(observed_error["MAE"]),
    observed_full_max = unname(observed_error["MAX"]),
    j6_observed = analysis$observed_curve[6L],
    j6_full = analysis$full_curve[6L],
    j6_agca = analysis$corrected_curve[6L],
    j12_observed = analysis$observed_curve[12L],
    j12_full = analysis$full_curve[12L],
    j12_agca = analysis$corrected_curve[12L]
  )
}))
rownames(split_summary) <- NULL

split_curves <- do.call(rbind, lapply(names(split_analyses), function(name) {
  analysis <- split_analyses[[name]]
  data.frame(
    sample = name,
    m = seq_len(analysis$d),
    observed = analysis$observed_curve,
    full = analysis$full_curve,
    agca = analysis$corrected_curve
  )
}))

ave_bootstrap <- bootstrap_ave(
  main_fit$g_extreme,
  main_fit$fit$mu,
  max_rank = max_plot_rank,
  reps = bootstrap_reps,
  seed = bootstrap_seed
)

coexceedance_bootstrap <- bootstrap_coexceedance_curves(
  main_fit$g_extreme,
  main_fit$fit$mu,
  tail_mass = main_fit$tail_mass,
  p = main_fit$p_star,
  reps = bootstrap_reps,
  seed = bootstrap_seed + 1L
)

write.csv(main_summary, file.path(output_dir, "main_summary.csv"),
          row.names = FALSE)
write.csv(main_curves, file.path(output_dir, "coexceedance_curves.csv"),
          row.names = FALSE)
write.csv(main_error_summary,
          file.path(output_dir, "main_error_summary.csv"), row.names = FALSE)
write.csv(main_fit$rank_summary,
          file.path(output_dir, "coexceedance_rank_errors.csv"),
          row.names = FALSE)
write.csv(ave_bootstrap, file.path(output_dir, "ave_bootstrap.csv"),
          row.names = FALSE)
write.csv(coexceedance_bootstrap,
          file.path(output_dir, "coexceedance_bootstrap.csv"),
          row.names = FALSE)
write.csv(severity_curves,
          file.path(output_dir, "coexceedance_severity_curves.csv"),
          row.names = FALSE)
write.csv(severity_error_summary,
          file.path(output_dir, "coexceedance_severity_summary.csv"),
          row.names = FALSE)
write.csv(threshold_rank_summary,
          file.path(output_dir, "coexceedance_threshold_rank_errors.csv"),
          row.names = FALSE)
write.csv(threshold_selected_summary,
          file.path(output_dir, "coexceedance_threshold_summary.csv"),
          row.names = FALSE)
write.csv(block_curves,
          file.path(output_dir, "coexceedance_block_survival.csv"),
          row.names = FALSE)
write.csv(block_pmf, file.path(output_dir, "coexceedance_block_pmf.csv"),
          row.names = FALSE)
write.csv(block_error_summary,
          file.path(output_dir, "coexceedance_block_summary.csv"),
          row.names = FALSE)
write.csv(positivity_summary,
          file.path(output_dir, "coexceedance_positivity_summary.csv"),
          row.names = FALSE)
write.csv(split_summary,
          file.path(output_dir, "coexceedance_chronological_summary.csv"),
          row.names = FALSE)
write.csv(split_curves,
          file.path(output_dir, "coexceedance_chronological_curves.csv"),
          row.names = FALSE)

save_pdf(
  file.path(output_dir, "fig_emp_ff_ave_rank_selection.pdf"),
  function() {
    par(mar = c(3.7, 4.3, 2.0, 0.8),
        mgp = c(2.2, 0.65, 0), tcl = -0.25)

    rank_data <- main_fit$rank_summary[
      main_fit$rank_summary$rank >= 1L &
        main_fit$rank_summary$rank <= max_plot_rank,
    ]
    plot(rank_data$rank, 100 * rank_data$variation_explained,
         type = "n", ylim = c(0, 100), xlab = "AGCA rank p",
         ylab = "Cumulative AVE (%)", main = "Anchored variation")
    grid(col = "gray90")
    polygon(
      c(ave_bootstrap$rank, rev(ave_bootstrap$rank)),
      100 * c(ave_bootstrap$lower, rev(ave_bootstrap$upper)),
      border = NA, col = grDevices::adjustcolor("#1B6CA8", 0.18)
    )
    abline(h = c(80, 90), col = "gray70", lty = 3, lwd = 0.9)
    abline(v = main_fit$p_star, col = "gray45", lty = 3, lwd = 0.9)
    lines(rank_data$rank, 100 * rank_data$variation_explained,
          type = "b", pch = 16, col = "#1B6CA8", lwd = 1.7)
    box()
  },
  width = 4.7, height = 3.2
)

save_pdf(
  file.path(output_dir, "fig_emp_ff_coexceedance_fidelity.pdf"),
  function() {
    par(mfrow = c(1L, 2L), mar = c(3.7, 4.2, 2.0, 0.7),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)

    rank_data <- main_fit$rank_summary[
      main_fit$rank_summary$rank >= 1L &
        main_fit$rank_summary$rank <= max_plot_rank,
    ]
    error_max <- 1.08 * max(10000 * rank_data$max_corrected)
    plot(rank_data$rank, 10000 * rank_data$mae_corrected,
         type = "n", ylim = c(0, error_max), xlab = "AGCA rank p",
         ylab = "Absolute probability error (bp)",
         main = "Compression error")
    grid(col = "gray90")
    abline(v = main_fit$p_star, col = "gray45", lty = 3, lwd = 0.9)
    lines(rank_data$rank, 10000 * rank_data$mae_corrected,
          type = "b", pch = 16, col = "#B23A48", lwd = 1.6)
    lines(rank_data$rank, 10000 * rank_data$max_corrected,
          type = "b", pch = 17, col = "#3B8C5A", lty = 2, lwd = 1.6)
    legend("topright", legend = c("MAE", "MAX"),
           col = c("#B23A48", "#3B8C5A"), pch = c(16, 17),
           lty = c(1, 2), lwd = 1.6, bty = "n", cex = 0.78)
    box()

    curve_intervals <- coexceedance_bootstrap[
      coexceedance_bootstrap$method %in% c("full", "agca"),
    ]
    draw_coexceedance_panel(
      main_curves,
      main = "Coexceedance probabilities",
      show_legend = TRUE,
      intervals = curve_intervals,
      reference_label = "Unreduced"
    )
  },
  width = 7.8, height = 2.75, pointsize = 12
)

save_pdf(
  file.path(output_dir, "fig_emp_ff_coexceedance_severity.pdf"),
  function() {
    par(mfrow = c(1L, 3L), mar = c(3.7, 4.0, 2.2, 0.5),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    common_ylim <- c(0, 1.04 * max(severity_curves$observed,
                                   severity_curves$full,
                                   severity_curves$agca))
    for (i in seq_along(severity_multipliers)) {
      x <- severity_multipliers[i]
      z <- severity_curves[severity_curves$x == x, ]
      draw_coexceedance_panel(
        z,
        main = sprintf("x=%d, q=%.4f", x, z$q_x[1L]),
        ylim = common_ylim,
        show_legend = i == length(severity_multipliers)
      )
    }
  },
  width = 8.2, height = 2.8
)

save_pdf(
  file.path(output_dir, "fig_emp_ff_coexceedance_threshold_sensitivity.pdf"),
  function() {
    par(mfrow = c(1L, 2L), mar = c(3.7, 4.3, 2.0, 0.8),
        mgp = c(2.2, 0.65, 0), tcl = -0.25)
    z <- threshold_rank_summary[
      threshold_rank_summary$rank >= 1L &
        threshold_rank_summary$rank <= max_plot_rank,
    ]
    colors <- grDevices::hcl.colors(length(threshold_fractions), "Dark 3")
    plot(NA, xlim = c(1, max_plot_rank),
         ylim = c(0, 1.04 * 10000 * max(z$mae_corrected)),
         xlab = "AGCA rank p", ylab = "Absolute probability error (bp)",
         main = "MAE")
    grid(col = "gray90")
    for (i in seq_along(threshold_fractions)) {
      zz <- z[z$target_tail_fraction == threshold_fractions[i], ]
      lines(zz$rank, 10000 * zz$mae_corrected, type = "b", pch = 16,
            col = colors[i], lwd = 1.2, cex = 0.62)
    }
    box()

    plot(NA, xlim = c(1, max_plot_rank),
         ylim = c(0, 1.04 * 10000 * max(z$max_corrected)),
         xlab = "AGCA rank p", ylab = "Absolute probability error (bp)",
         main = "MAX")
    grid(col = "gray90")
    for (i in seq_along(threshold_fractions)) {
      zz <- z[z$target_tail_fraction == threshold_fractions[i], ]
      lines(zz$rank, 10000 * zz$max_corrected, type = "b", pch = 16,
            col = colors[i], lwd = 1.2, cex = 0.62)
    }
    legend("topright", legend = sprintf("%.1f%%", 100 * threshold_fractions),
           title = "tail fraction", col = colors, lty = 1, pch = 16,
           lwd = 1.2, bty = "n", cex = 0.62, title.cex = 0.68)
    box()
  },
  width = 7.8, height = 2.75
)

save_pdf(
  file.path(output_dir, "fig_emp_ff_coexceedance_blocks.pdf"),
  function() {
    par(mfrow = c(2L, 2L), mar = c(3.5, 4.0, 2.0, 0.6),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    common_ylim <- c(0, 1.04 * max(block_curves$observed,
                                   block_curves$full,
                                   block_curves$agca))
    for (i in seq_along(block_ids)) {
      block <- block_ids[i]
      z <- block_curves[block_curves$block == block, ]
      draw_coexceedance_panel(
        z,
        main = z$block_label[1L],
        ylim = common_ylim,
        show_legend = i == 1L
      )
    }
  },
  width = 7.8, height = 5.0
)

save_pdf(
  file.path(output_dir, "fig_emp_ff_coexceedance_positivity.pdf"),
  function() {
    par(mfrow = c(1L, 2L), mar = c(3.7, 4.3, 3.2, 0.8),
        mgp = c(2.2, 0.65, 0), tcl = -0.25)
    raw_difference <- 10000 * (main_fit$raw_curve - main_fit$full_curve)
    corrected_difference <- 10000 * (
      main_fit$corrected_curve - main_fit$full_curve
    )
    difference_ylim <- range(raw_difference, corrected_difference)
    plot(seq_len(main_fit$d), raw_difference, type = "b", pch = 1,
         col = "#6A51A3", lwd = 1.5, ylim = difference_ylim,
         xlab = "Number of breaches m", ylab = "Difference from full (bp)",
         main = "Before and after correction")
    grid(col = "gray90")
    abline(h = 0, col = "gray55", lty = 3)
    lines(seq_len(main_fit$d), corrected_difference, type = "b", pch = 16,
          col = "#B23A48", lwd = 1.5, lty = 2)
    legend("top", inset = c(0, -0.17), horiz = TRUE, xpd = NA,
           legend = c("Raw", "Corrected"),
           col = c("#6A51A3", "#B23A48"), pch = c(1, 16),
           lty = c(1, 2), lwd = 1.5, bty = "n", cex = 0.76)
    box()

    correction_effect <- 10000 * (
      main_fit$corrected_curve - main_fit$raw_curve
    )
    plot(seq_len(main_fit$d), correction_effect, type = "b", pch = 16,
         col = "#3B8C5A", lwd = 1.5,
         xlab = "Number of breaches m", ylab = "Corrected minus raw (bp)",
         main = "Effect of positivity correction")
    grid(col = "gray90")
    abline(h = 0, col = "gray55", lty = 3)
    box()
  },
  width = 7.8, height = 3.0
)

save_pdf(
  file.path(output_dir, "fig_emp_ff_coexceedance_chronological.pdf"),
  function() {
    par(mfrow = c(1L, 2L), mar = c(3.7, 4.2, 2.2, 0.6),
        mgp = c(2.2, 0.65, 0), tcl = -0.25)
    common_ylim <- c(0, 1.04 * max(split_curves$observed,
                                   split_curves$full,
                                   split_curves$agca))
    for (i in seq_along(split_analyses)) {
      name <- names(split_analyses)[i]
      z <- split_curves[split_curves$sample == name, ]
      row <- split_summary[split_summary$sample == name, ]
      draw_coexceedance_panel(
        z,
        main = sprintf("%s--%s", format(row$start_date, "%Y"),
                       format(row$end_date, "%Y")),
        ylim = common_ylim,
        show_legend = i == 2L
      )
    }
  },
  width = 7.8, height = 2.75
)

if (any(diff(main_fit$full_curve) > 1e-12) ||
    any(diff(main_fit$corrected_curve) > 1e-12) ||
    any(diff(main_fit$observed_curve) > 1e-12)) {
  stop("A coexceedance survival curve is not monotone.", call. = FALSE)
}

stage_paper_figures()

cat("\nAGCA coexceedance diagnostics\n")
cat("  sample: ", main_fit$n, " days; k=", main_fit$k, "; d=",
    main_fit$d, "\n", sep = "")
cat(sprintf("  r_k=%.9f; q_star=%.8f\n",
            main_fit$radial_threshold, main_fit$q_star))
cat(sprintf("  p_star=%d; AVE=%.6f; MAE=%.8f; MAX=%.8f\n",
            main_fit$p_star, p_star_row$variation_explained,
            p_star_row$mae_corrected, p_star_row$max_corrected))
cat("  results: ", output_dir, "\n", sep = "")

}

run_panel_comparison <- function() {
  panel_specs <- list(
    ff_2x3_daily = list(label = "Fama--French", directory = "ff",
                        file = "ff_2x3_sorts_daily.rds"),
    osap_daily_quintile_vw = list(label = "OSAP", directory = "osap",
                                 file = "osap_daily_quintile_vw.rds")
  )
  factor_curve_file <- file.path(
    repo_dir, "inst", "empirics", "results", "portfolio_factor_attribution_agca",
    "factor_attribution_coexceedance_curves.csv"
  )
  cached_factor_curves <- if (file.exists(factor_curve_file)) {
    read.csv(factor_curve_file)
  } else {
    NULL
  }

  panels <- lapply(seq_along(panel_specs), function(panel_index) {
    dataset <- names(panel_specs)[panel_index]
    spec <- panel_specs[[panel_index]]
    prepared_panel <- readRDS(file.path(
      repo_dir, "data", "empirics", spec$directory, spec$file
    ))
    panel_data <- prepared_panel$complete_losses
    panel_data <- panel_data[panel_data$date >= sample_start, ]
    panel_data <- panel_data[order(panel_data$date), ]
    losses <- as.matrix(panel_data[, prepared_panel$portfolio_columns,
                                   drop = FALSE])
    analysis <- analyze_matrix(losses, panel_data$date)
    curves <- data.frame(
      dataset = dataset, m = seq_len(analysis$d),
      observed = analysis$observed_curve,
      full = analysis$full_curve, agca = analysis$corrected_curve,
      observed_minus_full = analysis$observed_curve - analysis$full_curve,
      agca_minus_full = analysis$corrected_curve - analysis$full_curve
    )
    for (method in c("observed", "full", "agca")) {
      values <- curves[[method]]
      stopifnot(all(is.finite(values)), all(values >= 0 & values <= 1),
                all(diff(values) <= 1e-12))
    }
    if (!is.null(cached_factor_curves)) {
      reference <- cached_factor_curves[
        cached_factor_curves$dataset == dataset, ]
      reference <- reference[order(reference$breach_count), ]
      stopifnot(nrow(reference) == analysis$d,
                max(abs(curves$full - reference$full)) < 1e-12,
                max(abs(curves$agca - reference$agca_selected)) < 1e-12)
    }
    # The Fama--French seed and bootstrap procedure are unchanged. OSAP uses
    # its own deterministic stream with the same paired, fixed-rank design.
    message("Bootstrapping coexceedance curves: ", spec$label)
    panel_seed <- bootstrap_seed + panel_index
    intervals <- bootstrap_coexceedance_curves(
      analysis$g_extreme, analysis$fit$mu, analysis$tail_mass,
      p = analysis$p_star, reps = bootstrap_reps, seed = panel_seed
    )
    intervals <- data.frame(dataset = dataset, intervals)
    stopifnot(all(is.finite(intervals$lower)),
              all(is.finite(intervals$upper)),
              all(intervals$lower <= intervals$median),
              all(intervals$median <= intervals$upper))
    selected <- analysis$rank_summary[
      analysis$rank_summary$rank == analysis$p_star, ]
    observed_error <- curve_errors(curves$observed, curves$full)
    marginal_counts <- colSums(analysis$observed_breaches)
    marginal_full <- analysis$tail_mass * colMeans(analysis$g_extreme)
    summary <- data.frame(
      dataset = dataset, n = analysis$n, d = analysis$d, k = analysis$k,
      start_date = analysis$dates[1L], end_date = analysis$dates[2L],
      radial_threshold = analysis$radial_threshold, q_star = analysis$q_star,
      p_star = analysis$p_star, ave_p_star = selected$variation_explained,
      agca_full_mae_bp = 1e4 * selected$mae_corrected,
      agca_full_max_bp = 1e4 * selected$max_corrected,
      observed_full_mae_bp = 1e4 * unname(observed_error["MAE"]),
      observed_full_max_bp = 1e4 * unname(observed_error["MAX"]),
      min_marginal_breaches = min(marginal_counts),
      max_marginal_breaches = max(marginal_counts),
      nominal_marginal_breaches = analysis$n / analysis$radial_threshold,
      min_fitted_marginal_probability = min(marginal_full),
      max_fitted_marginal_probability = max(marginal_full),
      bootstrap_reps = bootstrap_reps, bootstrap_seed = panel_seed
    )
    list(summary = summary, curves = curves, intervals = intervals,
         rank_errors = data.frame(dataset = dataset, analysis$rank_summary))
  })
  names(panels) <- names(panel_specs)

  # Check the original Fama--French display and bootstrap numerically, without
  # overwriting their saved files.
  ff_curve_file <- file.path(output_dir, "coexceedance_curves.csv")
  if (file.exists(ff_curve_file)) {
    reference <- read.csv(ff_curve_file)
    for (method in c("observed", "full", "agca")) {
      stopifnot(max(abs(panels[[1L]]$curves[[method]] -
                         reference[[method]])) < 1e-12)
    }
  }
  ff_band_file <- file.path(output_dir, "coexceedance_bootstrap.csv")
  if (file.exists(ff_band_file)) {
    reference <- read.csv(ff_band_file)
    stopifnot(isTRUE(all.equal(
      panels[[1L]]$intervals[, names(reference)], reference,
      tolerance = 1e-12, check.attributes = FALSE
    )))
  }

  combined <- lapply(c("summary", "curves", "intervals", "rank_errors"),
                     function(field) do.call(rbind, lapply(panels, `[[`, field)))
  names(combined) <- c("summary", "curves", "bootstrap", "rank_errors")
  for (field in names(combined)) {
    write.csv(combined[[field]],
              file.path(output_dir, paste0("panel_comparison_", field, ".csv")),
              row.names = FALSE)
  }

  panel_heading <- function(index) {
    result <- panels[[index]]$summary
    title(main = panel_specs[[index]]$label, line = 1.65, cex.main = 1.0)
    mtext(sprintf("p* = %d, q* = %.2f%%", result$p_star, 100 * result$q_star),
          side = 3, line = 0.4, cex = 0.82)
  }
  save_pdf(
    file.path(output_dir, "fig_emp_portfolio_coexceedance_survival.pdf"),
    function() {
      layout(matrix(c(1L, 1L, 2L, 3L), nrow = 2L, byrow = TRUE),
             heights = c(0.35, 2.85))
      par(mar = rep(0, 4), cex = 1)
      plot.new()
      legend("center", legend = c("Observed", "Full", "AGCA"),
             horiz = TRUE, bty = "n", cex = 0.82,
             col = curve_colors, lty = curve_lty, pch = curve_pch,
             lwd = c(1.2, 1.5, 1.5))
      par(mar = c(4.2, 4.6, 2.9, 0.8), las = 1,
          mgp = c(3.0, 0.65, 0), tcl = -0.25)
      curve_bands <- combined$bootstrap[
        combined$bootstrap$method %in% c("full", "agca"), ]
      common_ylim <- c(0, 1.04 * max(combined$curves$observed,
                                    combined$curves$full, combined$curves$agca,
                                    curve_bands$upper))
      for (i in seq_along(panels)) {
        d <- panels[[i]]$summary$d
        ticks <- pretty(c(1, d))
        ticks <- unique(c(1L, ticks[ticks > 1 & ticks < d], d))
        draw_coexceedance_panel(
          panels[[i]]$curves, ylim = common_ylim, show_legend = FALSE,
          intervals = panels[[i]]$intervals, count_ticks = ticks
        )
        panel_heading(i)
      }
    }, width = 7.8, height = 3.2
  )
  save_pdf(
    file.path(output_dir, "fig_emp_portfolio_coexceedance_difference.pdf"),
    function() {
      par(mfrow = c(1L, 2L), mar = c(3.7, 4.3, 2.9, 0.8),
          las = 1, mgp = c(2.2, 0.65, 0), tcl = -0.25)
      difference_bands <- combined$bootstrap[
        combined$bootstrap$method == "agca_minus_full", ]
      common_ylim <- 1e4 * range(difference_bands$lower,
                                 difference_bands$upper,
                                 combined$curves$agca_minus_full, 0)
      for (i in seq_along(panels)) {
        z <- panels[[i]]$curves
        band <- panels[[i]]$intervals
        band <- band[band$method == "agca_minus_full", ]
        plot(z$m, 1e4 * z$agca_minus_full, type = "n", ylim = common_ylim,
             xlab = "Number of exceedances m", ylab = "AGCA minus Full (bp)",
             xaxt = "n")
        ticks <- pretty(z$m)
        axis(1, at = unique(c(1L, ticks[ticks > 1 & ticks < max(z$m)], max(z$m))))
        grid(col = "gray90")
        polygon(c(band$m, rev(band$m)),
                1e4 * c(band$lower, rev(band$upper)), border = NA,
                col = grDevices::adjustcolor(curve_colors["corrected"], 0.25))
        abline(h = 0, col = "gray55", lty = 3)
        lines(z$m, 1e4 * z$agca_minus_full, type = "b", pch = 17,
              col = curve_colors["corrected"], lwd = 1.5, cex = 0.75)
        panel_heading(i)
        box()
      }
    }, width = 7.8, height = 3.0
  )
  stage_paper_figures()
  print(combined$summary, row.names = FALSE)
  invisible(combined)
}

run_panel_comparison()
