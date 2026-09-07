# Attribution of portfolio tail geometry to conventional factor exposures.

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

integer_argument <- function(name, default) {
  prefix <- paste0("--", name, "=")
  value <- commandArgs(TRUE)[startsWith(commandArgs(TRUE), prefix)]
  if (length(value) == 0L) {
    return(as.integer(default))
  }
  parsed <- suppressWarnings(as.integer(sub(prefix, "", value[1L],
                                             fixed = TRUE)))
  if (is.na(parsed) || parsed < 0L) {
    stop("--", name, " must be a nonnegative integer.", call. = FALSE)
  }
  parsed
}

repo_dir <- find_compendium_root(script_dir)
source(file.path(repo_dir, "R", "GeodesicExtreme.R"))

output_dir <- file.path(repo_dir, "inst", "empirics", "results",
                        "portfolio_factor_attribution_agca")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(dirname(repo_dir), "jrssb_paper", "figures")

sample_start <- as.Date("1973-07-01")
tail_fraction <- 0.05
threshold_fractions <- c(0.005, 0.01, 0.015, 0.025, 0.05, 0.075, 0.10)
candidate_factors <- c("SMB", "HML", "RMW", "CMA", "MOM")
nuisance_factors <- "MKT_RF"
bootstrap_reps <- integer_argument("factor-bootstrap-reps", 1000L)
block_length <- integer_argument("factor-block-length", 20L)
random_seed <- integer_argument("factor-seed", 20260826L)
skip_bootstrap <- "--skip-factor-bootstrap" %in% commandArgs(TRUE)
set.seed(random_seed)

dataset_specs <- list(
  ff_2x3_daily = list(
    label = "Fama--French 2 x 3",
    data_file = file.path(repo_dir, "data", "empirics", "ff",
                          "ff_2x3_sorts_daily.rds")
  ),
  osap_daily_quintile_vw = list(
    label = "OSAP anomalies",
    data_file = file.path(repo_dir, "data", "empirics", "osap",
                          "osap_daily_quintile_vw.rds")
  )
)
factor_file <- file.path(repo_dir, "data", "empirics", "ff",
                         "ff_six_factors_daily.rds")

if (!file.exists(factor_file)) {
  stop("Prepared factor data not found: ", factor_file,
       ". Run scripts/run_data.R first.", call. = FALSE)
}
factor_data <- readRDS(factor_file)$factors

rank_to_pareto <- function(x) {
  n <- length(x)
  ranks <- rank(x, ties.method = "average", na.last = "keep")
  (n + 1) / (n + 1 - ranks)
}

pareto_transform <- function(x) {
  transformed <- apply(x, 2L, rank_to_pareto)
  transformed <- as.matrix(transformed)
  colnames(transformed) <- colnames(x)
  transformed
}

select_top_k <- function(radius, k) {
  order(-radius, seq_along(radius))[seq_len(k)]
}

load_panel <- function(spec) {
  prepared <- readRDS(spec$data_file)
  returns <- prepared$complete_returns
  returns <- returns[returns$date >= sample_start, ]
  columns <- prepared$portfolio_columns
  merged <- merge(returns[, c("date", columns)], factor_data,
                  by = "date", all = FALSE)
  merged <- merged[complete.cases(merged), ]
  merged <- merged[order(merged$date), ]
  return_matrix <- as.matrix(merged[, columns, drop = FALSE])
  storage.mode(return_matrix) <- "double"
  factor_matrix <- as.matrix(
    merged[, c(nuisance_factors, candidate_factors, "RF"), drop = FALSE]
  )
  storage.mode(factor_matrix) <- "double"
  list(
    date = merged$date,
    returns = return_matrix,
    losses = -return_matrix,
    factors = factor_matrix,
    portfolio_columns = columns,
    metadata = prepared$portfolios
  )
}

subset_panel <- function(panel, index) {
  list(
    date = panel$date[index],
    returns = panel$returns[index, , drop = FALSE],
    losses = panel$losses[index, , drop = FALSE],
    factors = panel$factors[index, , drop = FALSE],
    portfolio_columns = panel$portfolio_columns,
    metadata = panel$metadata
  )
}

estimate_factor_exposures <- function(panel) {
  excess_returns <- sweep(panel$returns, 1L, panel$factors[, "RF"], "-")
  design <- cbind(
    Intercept = 1,
    panel$factors[, c(nuisance_factors, candidate_factors), drop = FALSE]
  )
  design_qr <- qr(design)
  if (design_qr$rank < ncol(design)) {
    stop("Factor-regression design is rank deficient.", call. = FALSE)
  }
  coefficients <- qr.coef(design_qr, excess_returns)
  rownames(coefficients) <- colnames(design)
  candidate <- t(coefficients[candidate_factors, , drop = FALSE])
  all_slopes <- t(coefficients[c(nuisance_factors, candidate_factors),
                               , drop = FALSE])
  colnames(candidate) <- candidate_factors
  colnames(all_slopes) <- c(nuisance_factors, candidate_factors)
  rownames(candidate) <- panel$portfolio_columns
  rownames(all_slopes) <- panel$portfolio_columns
  list(candidate = candidate, all = all_slopes)
}

factor_basis <- function(exposures, mu, tolerance = 1e-10) {
  exposures <- as.matrix(exposures)
  centered <- exposures - tcrossprod(mu, drop(crossprod(mu, exposures)))
  decomposition <- svd(centered, nu = ncol(centered), nv = 0L)
  threshold <- max(decomposition$d) * tolerance
  rank <- sum(decomposition$d > threshold)
  if (rank < ncol(centered)) {
    stop("Centered factor-exposure matrix is rank deficient.", call. = FALSE)
  }
  q <- decomposition$u[, seq_len(rank), drop = FALSE]
  projector <- tcrossprod(q)
  list(
    centered = centered,
    q = q,
    projector = projector,
    rank = rank,
    singular_values = decomposition$d,
    condition_number = max(decomposition$d) / min(decomposition$d)
  )
}

tail_fit <- function(losses, fraction = tail_fraction) {
  pareto <- pareto_transform(losses)
  directions <- normalize_rows(pareto)
  radius <- row_norms(pareto)
  k <- as.integer(round(fraction * nrow(directions)))
  index <- select_top_k(radius, k)
  selected <- directions[index, , drop = FALSE]
  mu <- canonical_anchor(ncol(selected))
  fit <- agca_fit(selected, mu = mu)
  rank_summary <- agca_rank_summary(fit)
  p_star <- min(rank_summary$rank[
    rank_summary$variation_explained >= 0.90
  ])
  list(
    fit = fit,
    mu = mu,
    index = index,
    radius = radius,
    threshold_radius = min(radius[index]),
    k = k,
    n = nrow(directions),
    p_star = p_star,
    rank_summary = rank_summary
  )
}

captured_fraction <- function(projector, sigma) {
  sum(projector * sigma) / sum(diag(sigma))
}

positive_reconstruction <- function(g) {
  normalize_rows(pmax(g, 0))
}

subspace_reconstruction <- function(fit, basis) {
  tangent <- fit$u %*% basis %*% t(basis)
  numerator <- tcrossprod(fit$anchor_coordinate, fit$mu) + tangent
  positive_reconstruction(normalize_rows(numerator))
}

coexceedance_curve <- function(g, tail_mass) {
  ordered <- t(apply(g, 1L, sort, decreasing = TRUE))
  tail_mass * colMeans(ordered)
}

curve_error <- function(estimate, target) {
  difference <- estimate - target
  c(
    mae_bps = 1e4 * mean(abs(difference)),
    max_bps = 1e4 * max(abs(difference))
  )
}

hybrid_factor_path <- function(fit, factor_geometry) {
  d <- ncol(fit$g)
  m <- factor_geometry$rank
  augmented <- cbind(fit$mu, factor_geometry$q)
  complete_qr <- qr.Q(qr(augmented), complete = TRUE)
  complement <- complete_qr[, (m + 2L):d, drop = FALSE]
  residual_sigma <- crossprod(complement, fit$sigma %*% complement)
  eig <- eigen((residual_sigma + t(residual_sigma)) / 2,
               symmetric = TRUE)
  residual_loadings <- complement %*% eig$vectors
  path <- lapply(0:ncol(residual_loadings), function(r) {
    basis <- if (r == 0L) {
      factor_geometry$q
    } else {
      cbind(factor_geometry$q,
            residual_loadings[, seq_len(r), drop = FALSE])
    }
    projector <- tcrossprod(basis)
    total_rank <- ncol(basis)
    optimal <- sum(fit$eigenvalues[seq_len(total_rank)]) /
      sum(fit$eigenvalues)
    data.frame(
      residual_components = r,
      total_rank = total_rank,
      factor_hybrid_ave = captured_fraction(projector, fit$sigma),
      optimal_agca_ave = optimal,
      efficiency = captured_fraction(projector, fit$sigma) / optimal
    )
  })
  list(path = do.call(rbind, path), loadings = residual_loadings)
}

analyze_tail_with_exposures <- function(tail, exposures, dataset_name) {
  fit <- tail$fit
  sigma <- fit$sigma
  total_variation <- sum(diag(sigma))
  geometry <- factor_basis(exposures, tail$mu)
  factor_ave <- captured_fraction(geometry$projector, sigma)
  factor_rank <- geometry$rank
  optimal_same_rank <- sum(fit$eigenvalues[seq_len(factor_rank)]) /
    sum(fit$eigenvalues)
  p_star <- tail$p_star
  p_star_projector <- tcrossprod(
    fit$loadings[, seq_len(p_star), drop = FALSE]
  )

  by_factor <- do.call(rbind, lapply(candidate_factors, function(factor) {
    single_geometry <- factor_basis(exposures[, factor, drop = FALSE],
                                    tail$mu)
    reduced <- exposures[, setdiff(candidate_factors, factor), drop = FALSE]
    reduced_geometry <- factor_basis(reduced, tail$mu)
    data.frame(
      dataset = dataset_name,
      factor = factor,
      single_factor_ave = captured_fraction(single_geometry$projector, sigma),
      alignment_with_selected_agca =
        sum(single_geometry$projector * p_star_projector),
      leave_one_out_increment =
        factor_ave - captured_fraction(reduced_geometry$projector, sigma),
      stringsAsFactors = FALSE
    )
  }))

  hybrid <- hybrid_factor_path(fit, geometry)
  r_90 <- min(hybrid$path$residual_components[
    hybrid$path$factor_hybrid_ave >= 0.90
  ])
  hybrid_basis <- if (r_90 == 0L) {
    geometry$q
  } else {
    cbind(geometry$q,
          hybrid$loadings[, seq_len(r_90), drop = FALSE])
  }
  hybrid_row <- hybrid$path[hybrid$path$residual_components == r_90, ]

  full_curve <- coexceedance_curve(fit$g, tail$k / tail$n)
  factor_curve <- coexceedance_curve(
    subspace_reconstruction(fit, geometry$q), tail$k / tail$n
  )
  agca_same_curve <- coexceedance_curve(
    positive_reconstruction(agca_reconstruct(fit, factor_rank)),
    tail$k / tail$n
  )
  agca_selected_curve <- coexceedance_curve(
    positive_reconstruction(agca_reconstruct(fit, p_star)),
    tail$k / tail$n
  )
  hybrid_curve <- coexceedance_curve(
    subspace_reconstruction(fit, hybrid_basis), tail$k / tail$n
  )
  errors <- rbind(
    factor = curve_error(factor_curve, full_curve),
    agca_same_rank = curve_error(agca_same_curve, full_curve),
    agca_selected = curve_error(agca_selected_curve, full_curve),
    factor_hybrid_90 = curve_error(hybrid_curve, full_curve)
  )

  tangent_residual <- fit$u - fit$u %*% geometry$projector
  asset_mse <- colMeans(tangent_residual^2)
  asset_residual <- data.frame(
    dataset = dataset_name,
    asset = colnames(fit$g),
    residual_mse = asset_mse,
    share_of_factor_residual_risk = asset_mse / sum(asset_mse),
    stringsAsFactors = FALSE
  )

  curves <- data.frame(
    dataset = dataset_name,
    breach_count = seq_along(full_curve),
    full = full_curve,
    factor = factor_curve,
    agca_same_rank = agca_same_curve,
    agca_selected = agca_selected_curve,
    factor_hybrid_90 = hybrid_curve
  )

  summary <- data.frame(
    dataset = dataset_name,
    n = tail$n,
    k = tail$k,
    dimension = ncol(fit$g),
    p_star = p_star,
    ave_p_star = tail$rank_summary$variation_explained[
      tail$rank_summary$rank == p_star
    ],
    factor_rank = factor_rank,
    factor_ave = factor_ave,
    optimal_same_rank_ave = optimal_same_rank,
    factor_efficiency = factor_ave / optimal_same_rank,
    exposure_condition_number = geometry$condition_number,
    residual_components_to_90 = r_90,
    hybrid_total_rank_90 = factor_rank + r_90,
    hybrid_ave_90 = hybrid_row$factor_hybrid_ave,
    factor_mae_bps = errors["factor", "mae_bps"],
    factor_max_bps = errors["factor", "max_bps"],
    agca_same_rank_mae_bps = errors["agca_same_rank", "mae_bps"],
    agca_same_rank_max_bps = errors["agca_same_rank", "max_bps"],
    agca_selected_mae_bps = errors["agca_selected", "mae_bps"],
    agca_selected_max_bps = errors["agca_selected", "max_bps"],
    hybrid_mae_bps = errors["factor_hybrid_90", "mae_bps"],
    hybrid_max_bps = errors["factor_hybrid_90", "max_bps"],
    total_anchored_variation = total_variation,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    by_factor = by_factor,
    hybrid_path = data.frame(dataset = dataset_name, hybrid$path),
    curves = curves,
    asset_residual = asset_residual,
    factor_geometry = geometry
  )
}

analyze_panel <- function(panel, dataset_name, exposures = NULL,
                          fraction = tail_fraction) {
  if (is.null(exposures)) {
    exposures <- estimate_factor_exposures(panel)$candidate
  }
  tail <- tail_fit(panel$losses, fraction = fraction)
  analysis <- analyze_tail_with_exposures(tail, exposures, dataset_name)
  analysis$tail <- tail
  analysis$exposures <- exposures
  analysis
}

circular_block_indices <- function(n, length) {
  blocks <- ceiling(n / length)
  starts <- sample.int(n, blocks, replace = TRUE)
  indices <- unlist(lapply(starts, function(start) {
    ((start - 1L + seq_len(length) - 1L) %% n) + 1L
  }), use.names = FALSE)
  indices[seq_len(n)]
}

metrics_to_long <- function(analysis, replicate) {
  global_names <- c(
    "factor_ave", "factor_efficiency", "exposure_condition_number",
    "residual_components_to_90", "hybrid_total_rank_90",
    "hybrid_ave_90", "factor_mae_bps", "factor_max_bps",
    "agca_same_rank_mae_bps", "agca_selected_mae_bps",
    "hybrid_mae_bps", "hybrid_max_bps"
  )
  global <- data.frame(
    dataset = analysis$summary$dataset,
    replicate = replicate,
    statistic = global_names,
    value = as.numeric(analysis$summary[1L, global_names]),
    stringsAsFactors = FALSE
  )
  factor_rows <- do.call(rbind, lapply(seq_len(nrow(analysis$by_factor)),
                                       function(i) {
    row <- analysis$by_factor[i, ]
    data.frame(
      dataset = row$dataset,
      replicate = replicate,
      statistic = paste(
        c("single_factor_ave", "alignment", "leave_one_out_increment"),
        row$factor,
        sep = "__"
      ),
      value = c(row$single_factor_ave,
                row$alignment_with_selected_agca,
                row$leave_one_out_increment),
      stringsAsFactors = FALSE
    )
  }))
  rbind(global, factor_rows)
}

point_metrics_to_long <- function(analysis) {
  metrics_to_long(analysis, replicate = 0L)
}

bootstrap_panel <- function(panel, dataset_name, reps, block_length) {
  if (reps == 0L) {
    return(data.frame())
  }
  output <- vector("list", reps)
  for (b in seq_len(reps)) {
    index <- circular_block_indices(nrow(panel$returns), block_length)
    boot_panel <- subset_panel(panel, index)
    output[[b]] <- tryCatch({
      fit <- analyze_panel(boot_panel, dataset_name)
      metrics_to_long(fit, replicate = b)
    }, error = function(e) {
      warning("Bootstrap replicate ", b, " failed for ", dataset_name,
              ": ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (b %% 50L == 0L || b == reps) {
      message("Factor bootstrap ", dataset_name, ": ", b, "/", reps)
    }
  }
  output <- Filter(Negate(is.null), output)
  if (length(output) == 0L) {
    stop("All factor bootstrap replicates failed for ", dataset_name, ".",
         call. = FALSE)
  }
  do.call(rbind, output)
}

summarize_bootstrap <- function(points, bootstrap) {
  if (nrow(bootstrap) == 0L) {
    return(data.frame())
  }
  keys <- split(bootstrap, paste(bootstrap$dataset, bootstrap$statistic,
                                 sep = "\r"))
  summaries <- do.call(rbind, lapply(keys, function(z) {
    point <- points$value[
      points$dataset == z$dataset[1L] &
        points$statistic == z$statistic[1L]
    ][1L]
    data.frame(
      dataset = z$dataset[1L],
      statistic = z$statistic[1L],
      estimate = point,
      lower = unname(quantile(z$value, 0.025, na.rm = TRUE)),
      median = unname(quantile(z$value, 0.5, na.rm = TRUE)),
      upper = unname(quantile(z$value, 0.975, na.rm = TRUE)),
      bootstrap_sd = sd(z$value, na.rm = TRUE),
      successful_replicates = sum(is.finite(z$value)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summaries) <- NULL
  summaries[order(summaries$dataset, summaries$statistic), ]
}

chronological_analysis <- function(panel, dataset_name) {
  split_point <- floor(nrow(panel$returns) / 2L)
  halves <- list(
    early = subset_panel(panel, seq_len(split_point)),
    late = subset_panel(panel, (split_point + 1L):nrow(panel$returns))
  )
  exposures <- lapply(halves, function(z) {
    estimate_factor_exposures(z)$candidate
  })
  geometry <- lapply(exposures, function(b) {
    factor_basis(b, canonical_anchor(ncol(panel$returns)))
  })
  projector_distance <- sqrt(max(
    0,
    length(candidate_factors) -
      sum(geometry$early$projector * geometry$late$projector)
  ))
  rows <- list()
  counter <- 1L
  for (evaluation in names(halves)) {
    for (exposure_period in names(exposures)) {
      result <- analyze_panel(
        halves[[evaluation]], dataset_name,
        exposures = exposures[[exposure_period]]
      )
      rows[[counter]] <- data.frame(
        dataset = dataset_name,
        evaluation_period = evaluation,
        exposure_period = exposure_period,
        evaluation_start = min(halves[[evaluation]]$date),
        evaluation_end = max(halves[[evaluation]]$date),
        n_evaluation = nrow(halves[[evaluation]]$returns),
        factor_ave = result$summary$factor_ave,
        optimal_same_rank_ave = result$summary$optimal_same_rank_ave,
        factor_efficiency = result$summary$factor_efficiency,
        factor_mae_bps = result$summary$factor_mae_bps,
        exposure_projector_distance = projector_distance,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  do.call(rbind, rows)
}

threshold_analysis <- function(panel, dataset_name, exposures) {
  do.call(rbind, lapply(threshold_fractions, function(fraction) {
    result <- analyze_panel(panel, dataset_name, exposures = exposures,
                            fraction = fraction)
    data.frame(
      dataset = dataset_name,
      tail_fraction = fraction,
      k = result$summary$k,
      p_star = result$summary$p_star,
      factor_ave = result$summary$factor_ave,
      optimal_same_rank_ave = result$summary$optimal_same_rank_ave,
      factor_efficiency = result$summary$factor_efficiency,
      residual_components_to_90 =
        result$summary$residual_components_to_90,
      hybrid_total_rank_90 = result$summary$hybrid_total_rank_90,
      stringsAsFactors = FALSE
    )
  }))
}

save_pdf <- function(file, plot_function, width, height) {
  pdf(file, width = width, height = height, bg = "white")
  on.exit(dev.off(), add = TRUE)
  plot_function()
}

plot_factor_attribution <- function(file, bootstrap_summary) {
  save_pdf(file, function() {
    par(mfrow = c(2L, 2L), mar = c(4.2, 4.3, 2.0, 0.8),
        mgp = c(2.25, 0.65, 0), tcl = -0.25)
    for (metric in c("single_factor_ave", "leave_one_out_increment")) {
      for (dataset_name in names(dataset_specs)) {
        statistics <- paste(metric, candidate_factors, sep = "__")
        z <- bootstrap_summary[
          bootstrap_summary$dataset == dataset_name &
            bootstrap_summary$statistic %in% statistics,
        ]
        z <- z[match(statistics, z$statistic), ]
        values <- 100 * z$estimate
        lower <- 100 * z$lower
        upper <- 100 * z$upper
        ylim <- range(c(0, lower, upper), finite = TRUE)
        padding <- max(0.5, 0.08 * diff(ylim))
        plot(
          seq_along(candidate_factors), values,
          type = "n", xaxt = "n", xlab = "", ylim = ylim + c(-padding, padding),
          ylab = if (metric == "single_factor_ave") {
            "Single-axis attributed variation (%)"
          } else {
            "Leave-one-out increment (%)"
          },
          main = dataset_specs[[dataset_name]]$label
        )
        grid(nx = NA, ny = NULL, col = "gray90")
        arrows(seq_along(values), lower, seq_along(values), upper,
               angle = 90, code = 3, length = 0.045,
               col = "#4D4D4D", lwd = 1.1)
        points(seq_along(values), values, pch = 16, cex = 1.05,
               col = "#1B6CA8")
        axis(1, at = seq_along(candidate_factors),
             labels = candidate_factors, las = 2)
        abline(h = 0, col = "gray60", lty = 3)
        box()
      }
    }
  }, width = 8.2, height = 6.2)
}

plot_hybrid_path <- function(file, hybrid_path) {
  save_pdf(file, function() {
    par(mfrow = c(1L, 2L), mar = c(3.8, 4.3, 2.0, 0.8),
        mgp = c(2.25, 0.65, 0), tcl = -0.25)
    for (dataset_name in names(dataset_specs)) {
      z <- hybrid_path[hybrid_path$dataset == dataset_name, ]
      plot(
        z$total_rank, 100 * z$optimal_agca_ave,
        type = "n", ylim = c(min(50, 100 * min(z$factor_hybrid_ave)), 100),
        xlab = "Total subspace dimension",
        ylab = "Anchored variation captured (%)",
        main = dataset_specs[[dataset_name]]$label
      )
      grid(col = "gray90")
      abline(h = 90, lty = 3, col = "gray55")
      lines(z$total_rank, 100 * z$optimal_agca_ave, type = "b",
            pch = 1, col = "#1B6CA8", lwd = 1.5)
      lines(z$total_rank, 100 * z$factor_hybrid_ave, type = "b",
            pch = 16, col = "#B23A48", lwd = 1.5)
      legend("bottomright", c("Optimal AGCA", "Factors + residual AGCs"),
             col = c("#1B6CA8", "#B23A48"), pch = c(1, 16),
             lty = 1, lwd = 1.5, bty = "n", cex = 0.78)
      box()
    }
  }, width = 8.0, height = 3.1)
}

panels <- lapply(names(dataset_specs), function(dataset_name) {
  message("Analyzing factor attribution for ", dataset_name, ".")
  panel <- load_panel(dataset_specs[[dataset_name]])
  exposure_fit <- estimate_factor_exposures(panel)
  main <- analyze_panel(panel, dataset_name,
                        exposures = exposure_fit$candidate)
  exposures <- data.frame(
    dataset = dataset_name,
    asset = rownames(exposure_fit$all),
    exposure_fit$all,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(
    panel = panel,
    main = main,
    exposures = exposures,
    chronology = chronological_analysis(panel, dataset_name),
    threshold = threshold_analysis(panel, dataset_name,
                                   exposure_fit$candidate)
  )
})
names(panels) <- names(dataset_specs)

main_summary <- do.call(rbind, lapply(panels, function(z) z$main$summary))
by_factor <- do.call(rbind, lapply(panels, function(z) z$main$by_factor))
hybrid_path <- do.call(rbind, lapply(panels, function(z) z$main$hybrid_path))
coexceedance_curves <- do.call(rbind, lapply(panels,
                                             function(z) z$main$curves))
asset_residuals <- do.call(rbind, lapply(panels,
                                         function(z) z$main$asset_residual))
exposures <- do.call(rbind, lapply(panels, `[[`, "exposures"))
chronology <- do.call(rbind, lapply(panels, `[[`, "chronology"))
threshold <- do.call(rbind, lapply(panels, `[[`, "threshold"))
point_metrics <- do.call(rbind, lapply(panels, function(z) {
  point_metrics_to_long(z$main)
}))

bootstrap <- if (skip_bootstrap || bootstrap_reps == 0L) {
  data.frame()
} else {
  do.call(rbind, lapply(names(panels), function(dataset_name) {
    bootstrap_panel(
      panels[[dataset_name]]$panel,
      dataset_name,
      reps = bootstrap_reps,
      block_length = block_length
    )
  }))
}
bootstrap_summary <- summarize_bootstrap(point_metrics, bootstrap)

write.csv(main_summary, file.path(output_dir, "factor_attribution_summary.csv"),
          row.names = FALSE)
write.csv(by_factor, file.path(output_dir, "factor_attribution_by_factor.csv"),
          row.names = FALSE)
write.csv(exposures, file.path(output_dir, "factor_exposures.csv"),
          row.names = FALSE)
write.csv(asset_residuals,
          file.path(output_dir, "factor_attribution_asset_residuals.csv"),
          row.names = FALSE)
write.csv(hybrid_path,
          file.path(output_dir, "factor_attribution_hybrid_path.csv"),
          row.names = FALSE)
write.csv(coexceedance_curves,
          file.path(output_dir, "factor_attribution_coexceedance_curves.csv"),
          row.names = FALSE)
write.csv(chronology,
          file.path(output_dir, "factor_attribution_chronology.csv"),
          row.names = FALSE)
write.csv(threshold,
          file.path(output_dir, "factor_attribution_threshold_sensitivity.csv"),
          row.names = FALSE)
write.csv(bootstrap,
          file.path(output_dir, "factor_attribution_block_bootstrap.csv"),
          row.names = FALSE)
write.csv(bootstrap_summary,
          file.path(output_dir, "factor_attribution_block_bootstrap_summary.csv"),
          row.names = FALSE)

if (nrow(bootstrap_summary) > 0L) {
  attribution_figure <- file.path(output_dir,
                                  "fig_emp_factor_attribution.pdf")
  plot_factor_attribution(attribution_figure, bootstrap_summary)
} else {
  attribution_figure <- character(0L)
}
hybrid_figure <- file.path(output_dir, "fig_emp_factor_residual_path.pdf")
plot_hybrid_path(hybrid_figure, hybrid_path)

if (dir.exists(dirname(figure_dir))) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  staged <- c(attribution_figure, hybrid_figure)
  if (length(staged) > 0L) {
    file.copy(staged, figure_dir, overwrite = TRUE)
  }
}

cat("\nAGCA factor-attribution analysis\n")
print(main_summary, row.names = FALSE, digits = 4)
cat("\nFactor-specific attribution\n")
print(by_factor, row.names = FALSE, digits = 4)
cat("\nChronological exposure validation\n")
print(chronology, row.names = FALSE, digits = 4)
cat("\nOutputs written to:\n  ", output_dir, "\n", sep = "")
