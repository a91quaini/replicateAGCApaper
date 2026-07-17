# Time-dependence robustness checks for the portfolio AGCA empirical section.

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
        grepl("^Package: replicateAGCApaper", readLines(description, n = 1L))) {
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

sample_start <- as.Date("1973-07-01")
tail_fraction <- 0.05
max_rank <- 15L
selected_ranks <- c(3L, 5L, 8L, 10L, 12L)
decluster_run_length <- 3L

ff_data_file <- file.path(repo_dir, "data", "empirics", "ff",
                          "ff_2x3_sorts_daily.rds")

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

softmax_garch_parameters <- function(theta) {
  if (!all(is.finite(theta)) || theta[2L] < -50 || theta[2L] > 50) {
    return(NULL)
  }
  m <- max(0, theta[3L], theta[4L])
  e0 <- exp(-m)
  ea <- exp(theta[3L] - m)
  eb <- exp(theta[4L] - m)
  denom <- e0 + ea + eb
  c(
    mu = theta[1L],
    omega = exp(theta[2L]),
    alpha = ea / denom,
    beta = eb / denom
  )
}

garch_filter <- function(x, pars) {
  eps <- x - pars[["mu"]]
  n <- length(eps)
  sigma2 <- numeric(n)
  unconditional <- pars[["omega"]] /
    max(1 - pars[["alpha"]] - pars[["beta"]], 1e-8)
  sigma2[1L] <- max(unconditional, stats::var(eps), 1e-8)
  if (!is.finite(sigma2[1L])) {
    sigma2[1L] <- max(stats::var(eps), 1e-8)
  }
  for (i in 2L:n) {
    sigma2[i] <- pars[["omega"]] + pars[["alpha"]] * eps[i - 1L]^2 +
      pars[["beta"]] * sigma2[i - 1L]
    if (!is.finite(sigma2[i]) || sigma2[i] <= 0) {
      sigma2[i] <- NA_real_
      break
    }
  }
  list(
    eps = eps,
    sigma2 = sigma2,
    std_residual = eps / sqrt(sigma2)
  )
}

garch_objective <- function(theta, x) {
  pars <- softmax_garch_parameters(theta)
  if (is.null(pars)) {
    return(1e12)
  }
  filtered <- garch_filter(x, pars)
  sigma2 <- filtered$sigma2
  if (any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
    return(1e12)
  }
  eps <- filtered$eps
  value <- 0.5 * sum(log(sigma2) + eps^2 / sigma2)
  if (is.finite(value)) value else 1e12
}

fit_garch_11 <- function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("GARCH input contains non-finite values.", call. = FALSE)
  }
  x_var <- max(stats::var(x), 1e-8)
  alpha0 <- 0.05
  beta0 <- 0.90
  gap0 <- 1 - alpha0 - beta0
  start <- c(
    mean(x),
    log(x_var * gap0),
    log(alpha0 / gap0),
    log(beta0 / gap0)
  )
  fit <- try(
    stats::optim(
      start, garch_objective, x = x, method = "BFGS",
      control = list(maxit = 1000L, reltol = 1e-8)
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error") || !is.finite(fit$value)) {
    fit <- stats::optim(
      start, garch_objective, x = x, method = "Nelder-Mead",
      control = list(maxit = 3000L, reltol = 1e-8)
    )
  }
  pars <- softmax_garch_parameters(fit$par)
  filtered <- garch_filter(x, pars)
  list(
    convergence = fit$convergence,
    objective = fit$value,
    par = pars,
    sigma = sqrt(filtered$sigma2),
    std_residual = filtered$std_residual
  )
}

lag1_cor <- function(x) {
  x0 <- x[-length(x)]
  x1 <- x[-1L]
  ok <- is.finite(x0) & is.finite(x1)
  if (sum(ok) < 3L || stats::sd(x0[ok]) == 0 || stats::sd(x1[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor(x0[ok], x1[ok])
}

serial_diagnostic_detail <- function(series_matrix, sample_label) {
  out <- vector("list", ncol(series_matrix))
  for (j in seq_len(ncol(series_matrix))) {
    x <- series_matrix[, j]
    loss <- -x
    threshold <- stats::quantile(loss, probs = 0.95, type = 8, na.rm = TRUE)
    tail_indicator <- as.numeric(loss >= threshold)
    out[[j]] <- data.frame(
      sample = sample_label,
      asset = colnames(series_matrix)[j],
      return = lag1_cor(x),
      absolute_return = lag1_cor(abs(x)),
      squared_return = lag1_cor(x^2),
      loss_tail_indicator = lag1_cor(tail_indicator),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

summarize_serial_diagnostics <- function(detail) {
  statistic_columns <- c("return", "absolute_return", "squared_return",
                         "loss_tail_indicator")
  out <- list()
  index <- 1L
  for (sample_label in unique(detail$sample)) {
    z <- detail[detail$sample == sample_label, , drop = FALSE]
    for (statistic in statistic_columns) {
      values <- z[[statistic]]
      out[[index]] <- data.frame(
        sample = sample_label,
        statistic = statistic,
        min = min(values, na.rm = TRUE),
        q25 = unname(stats::quantile(values, 0.25, na.rm = TRUE)),
        median = stats::median(values, na.rm = TRUE),
        q75 = unname(stats::quantile(values, 0.75, na.rm = TRUE)),
        max = max(values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  do.call(rbind, out)
}

decluster_exceedances <- function(radius, exceedance_index, run_length) {
  exceedance_index <- sort(exceedance_index)
  if (length(exceedance_index) == 0L) {
    stop("Declustering requires at least one exceedance.", call. = FALSE)
  }
  cluster_id <- integer(length(exceedance_index))
  current_cluster <- 1L
  cluster_id[1L] <- current_cluster
  if (length(exceedance_index) > 1L) {
    for (i in seq.int(2L, length(exceedance_index))) {
      non_exceedance_gap <- exceedance_index[i] - exceedance_index[i - 1L] - 1L
      if (non_exceedance_gap >= run_length) {
        current_cluster <- current_cluster + 1L
      }
      cluster_id[i] <- current_cluster
    }
  }

  clusters <- split(exceedance_index, cluster_id)
  do.call(rbind, lapply(seq_along(clusters), function(i) {
    cluster <- clusters[[i]]
    representative <- cluster[which.max(radius[cluster])]
    data.frame(
      cluster = i,
      start_index = min(cluster),
      end_index = max(cluster),
      representative_index = representative,
      start_date = names(radius)[min(cluster)],
      end_date = names(radius)[max(cluster)],
      representative_date = names(radius)[representative],
      cluster_size = length(cluster),
      representative_radius = radius[representative],
      stringsAsFactors = FALSE
    )
  }))
}

agca_spectrum <- function(loss_matrix, sample_label, index = NULL) {
  pareto <- pareto_transform(loss_matrix)
  g <- normalize_rows(pareto)
  radius <- row_norms(pareto)
  if (is.null(index)) {
    k <- as.integer(round(tail_fraction * nrow(g)))
    index <- select_top_k(radius, k)
  } else {
    k <- length(index)
  }
  g_extreme <- g[index, , drop = FALSE]
  fit <- agca_fit(g_extreme, mu = canonical_anchor(ncol(g_extreme)))
  rank_summary <- agca_rank_summary(fit)
  rank_summary <- rank_summary[rank_summary$rank <= max_rank, , drop = FALSE]
  data.frame(
    sample = sample_label,
    n = nrow(loss_matrix),
    d = ncol(loss_matrix),
    k = k,
    tail_fraction = k / nrow(loss_matrix),
    rank_summary,
    stringsAsFactors = FALSE
  )
}

summarize_selected_ranks <- function(spectrum) {
  out <- list()
  for (sample_label in unique(spectrum$sample)) {
    z <- spectrum[spectrum$sample == sample_label, , drop = FALSE]
    row <- data.frame(
      sample = sample_label,
      n = z$n[1L],
      d = z$d[1L],
      k = z$k[1L],
      rank_80 = min(z$rank[z$variation_explained >= 0.80]),
      rank_90 = min(z$rank[z$variation_explained >= 0.90]),
      stringsAsFactors = FALSE
    )
    for (rank in selected_ranks) {
      value <- z$variation_explained[z$rank == rank]
      row[[paste0("ave_rank_", rank)]] <- if (length(value) == 1L) value else NA_real_
    }
    out[[sample_label]] <- row
  }
  do.call(rbind, out)
}

save_pdf <- function(file, plot_fun, width = 6.5, height = 4.2) {
  pdf(file, width = width, height = height, bg = "white")
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

plot_garch_spectrum <- function(file, spectrum) {
  save_pdf(file, function() {
    raw <- spectrum[spectrum$sample == "Raw returns", , drop = FALSE]
    filtered <- spectrum[spectrum$sample == "GARCH standardized innovations", ,
                         drop = FALSE]
    declustered <- spectrum[spectrum$sample == "Declustered raw exceedances",
                            , drop = FALSE]
    plot(
      raw$rank, 100 * raw$variation_explained,
      type = "n",
      xlim = c(0, max_rank),
      ylim = c(0, 100),
      xlab = "AGCA rank p",
      ylab = "Cumulative AVE (%)",
      main = "Fama-French canonical anchor"
    )
    grid(col = "gray90")
    abline(h = c(80, 90), col = "gray70", lty = 3, lwd = 0.9)
    lines(raw$rank, 100 * raw$variation_explained, type = "b",
          pch = 16, col = "#1B6CA8", lwd = 1.7)
    lines(filtered$rank, 100 * filtered$variation_explained, type = "b",
          pch = 17, col = "#B23A48", lwd = 1.7, lty = 2)
    lines(declustered$rank, 100 * declustered$variation_explained, type = "b",
          pch = 15, col = "#3B8C5A", lwd = 1.7, lty = 3)
    legend("bottomright",
           legend = c("Raw returns", "GARCH filtered",
                      "Declustered raw (m=3)"),
           col = c("#1B6CA8", "#B23A48", "#3B8C5A"),
           pch = c(16, 17, 15), lty = c(1, 2, 3), lwd = 1.7,
           bty = "n", cex = 0.76)
    box()
  }, width = 6.4, height = 4.0)
}

stage_manuscript_outputs <- function(pdf_file) {
  figure_dir <- file.path(dirname(repo_dir), "jrssb_paper", "figures")
  if (!dir.exists(figure_dir)) {
    return(invisible(FALSE))
  }
  file.copy(pdf_file,
            file.path(figure_dir, "fig_emp_garch_canonical_spectrum.pdf"),
            overwrite = TRUE)
  invisible(TRUE)
}

prepared <- readRDS(ff_data_file)
return_data <- prepared$complete_returns
return_data <- return_data[return_data$date >= sample_start, , drop = FALSE]
portfolio_columns <- prepared$portfolio_columns
return_matrix <- as.matrix(return_data[, portfolio_columns, drop = FALSE])
storage.mode(return_matrix) <- "double"
rownames(return_matrix) <- format(return_data$date, "%Y-%m-%d")

garch_fits <- vector("list", ncol(return_matrix))
names(garch_fits) <- colnames(return_matrix)
std_residuals <- matrix(NA_real_, nrow(return_matrix), ncol(return_matrix))
colnames(std_residuals) <- colnames(return_matrix)
rownames(std_residuals) <- rownames(return_matrix)

for (j in seq_len(ncol(return_matrix))) {
  asset <- colnames(return_matrix)[j]
  message("Fitting GARCH(1,1) for ", asset)
  garch_fits[[j]] <- fit_garch_11(return_matrix[, j])
  std_residuals[, j] <- garch_fits[[j]]$std_residual
}

garch_fit_summary <- do.call(rbind, lapply(seq_along(garch_fits), function(j) {
  fit <- garch_fits[[j]]
  data.frame(
    asset = names(garch_fits)[j],
    convergence = fit$convergence,
    objective = fit$objective,
    mu = fit$par[["mu"]],
    omega = fit$par[["omega"]],
    alpha = fit$par[["alpha"]],
    beta = fit$par[["beta"]],
    persistence = fit$par[["alpha"]] + fit$par[["beta"]],
    stringsAsFactors = FALSE
  )
}))

diagnostic_detail <- rbind(
  serial_diagnostic_detail(return_matrix, "Raw returns"),
  serial_diagnostic_detail(std_residuals, "GARCH standardized innovations")
)
diagnostic_summary <- summarize_serial_diagnostics(diagnostic_detail)

raw_loss_matrix <- -return_matrix
filtered_loss_matrix <- -std_residuals
raw_pareto <- pareto_transform(raw_loss_matrix)
raw_radius <- row_norms(raw_pareto)
names(raw_radius) <- rownames(raw_loss_matrix)
raw_k <- as.integer(round(tail_fraction * length(raw_radius)))
raw_exceedance_index <- select_top_k(raw_radius, raw_k)
declustered_exceedances <- decluster_exceedances(
  raw_radius, raw_exceedance_index, decluster_run_length
)
declustered_index <- declustered_exceedances$representative_index
spectrum <- rbind(
  agca_spectrum(raw_loss_matrix, "Raw returns"),
  agca_spectrum(filtered_loss_matrix, "GARCH standardized innovations"),
  agca_spectrum(raw_loss_matrix, "Declustered raw exceedances",
                index = declustered_index)
)
spectrum_summary <- summarize_selected_ranks(spectrum)

write.csv(garch_fit_summary,
          file.path(output_dir, "garch_fit_summary.csv"),
          row.names = FALSE)
write.csv(diagnostic_detail,
          file.path(output_dir, "garch_serial_dependence_detail.csv"),
          row.names = FALSE)
write.csv(diagnostic_summary,
          file.path(output_dir, "garch_serial_dependence_summary.csv"),
          row.names = FALSE)
write.csv(declustered_exceedances,
          file.path(output_dir, "declustered_raw_exceedances_m3.csv"),
          row.names = FALSE)
write.csv(spectrum,
          file.path(output_dir, "garch_canonical_spectrum.csv"),
          row.names = FALSE)
write.csv(spectrum_summary,
          file.path(output_dir, "garch_canonical_spectrum_summary.csv"),
          row.names = FALSE)

garch_spectrum_file <- file.path(output_dir, "garch_canonical_spectrum.pdf")
plot_garch_spectrum(garch_spectrum_file, spectrum)
stage_manuscript_outputs(garch_spectrum_file)

cat("\nPortfolio time-dependence robustness outputs written to:\n  ",
    output_dir, "\n",
    sep = "")
