# Portfolio tail-function diagnostics for fitted AGCA reconstructions.

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
                        "portfolio_tail_functionals_agca")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
manuscript_output_dir <- file.path(dirname(repo_dir), "full_paper", "R", "empirics",
                                   "empirics_output",
                                   "portfolio_tail_functionals_agca")

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
alpha <- 1
random_portfolios <- 250L
leverage_limit <- 1.5
functional_bootstrap_reps <- 99L
set.seed(20260707)

dataset_specs <- list(
  ff_2x3_daily = list(
    label = "FF daily 2x3 size-based sorts",
    data_file = file.path(repo_dir, "data", "empirics", "ff",
                          "ff_2x3_sorts_daily.rds"),
    group_column = "sort",
    group_label_column = "sort_label",
    max_plot_rank = 15L
  ),
  osap_daily_quintile_vw = list(
    label = "OSAP daily VW quintile liquidity/trading sorts",
    data_file = file.path(repo_dir, "data", "empirics", "osap",
                          "osap_daily_quintile_vw.rds"),
    group_column = "signal",
    group_label_column = "signal_label",
    max_plot_rank = 15L
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
  if (length(k) != 1L || k < 1L || k > length(radius) || k != as.integer(k)) {
    stop("k must be an integer between 1 and the number of observations.",
         call. = FALSE)
  }
  order(radius, decreasing = TRUE)[seq_len(k)]
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

positive_post_project <- function(g) {
  gp <- pmax(g, 0)
  normalize_rows(gp)
}

normalize_columns_to_one <- function(x) {
  sweep(x, 2L, colSums(x), "/")
}

random_long_only_weights <- function(d, n) {
  weights <- matrix(rexp(d * n), d, n)
  normalize_columns_to_one(weights)
}

random_limited_leverage_weights <- function(d, n, leverage_limit) {
  short_gross <- runif(n, min = 0, max = (leverage_limit - 1) / 2)
  long_gross <- 1 + short_gross
  long_leg <- normalize_columns_to_one(matrix(rexp(d * n), d, n))
  short_leg <- normalize_columns_to_one(matrix(rexp(d * n), d, n))
  sweep(long_leg, 2L, long_gross, "*") -
    sweep(short_leg, 2L, short_gross, "*")
}

make_portfolios <- function(portfolio_columns, metadata, group_column,
                            group_label_column, n_random, leverage_limit) {
  d <- length(portfolio_columns)
  weights <- list()
  labels <- list()

  add_weight <- function(id, class, label, w) {
    weights[[id]] <<- as.numeric(w)
    labels[[id]] <<- data.frame(
      portfolio_id = id,
      portfolio_class = class,
      portfolio_label = label,
      gross_exposure = sum(abs(w)),
      l2_norm = sqrt(sum(w^2)),
      stringsAsFactors = FALSE
    )
  }

  add_weight("equal_all", "equal", "Equal weight, all assets", rep(1 / d, d))

  groups <- unique(metadata[[group_column]])
  for (group in groups) {
    variables <- metadata$variable[metadata[[group_column]] == group]
    index <- match(variables, portfolio_columns)
    index <- index[!is.na(index)]
    if (length(index) > 0L) {
      w <- rep(0, d)
      w[index] <- 1 / length(index)
      group_label <- unique(metadata[[group_label_column]][
        metadata[[group_column]] == group
      ])[1L]
      add_weight(
        paste0("block_", group),
        "block_equal",
        paste0("Equal weight, ", group_label),
        w
      )
    }
  }

  long_only <- random_long_only_weights(d, n_random)
  for (j in seq_len(ncol(long_only))) {
    add_weight(
      sprintf("long_only_%03d", j),
      "random_long_only",
      sprintf("Random long-only %03d", j),
      long_only[, j]
    )
  }

  leveraged <- random_limited_leverage_weights(d, n_random, leverage_limit)
  for (j in seq_len(ncol(leveraged))) {
    add_weight(
      sprintf("limited_leverage_%03d", j),
      "random_limited_leverage",
      sprintf("Random limited-leverage %03d", j),
      leveraged[, j]
    )
  }

  weight_matrix <- do.call(cbind, weights)
  rownames(weight_matrix) <- portfolio_columns
  portfolio_info <- do.call(rbind, labels)
  rownames(portfolio_info) <- NULL

  list(weights = weight_matrix, info = portfolio_info)
}

capped_pareto_expectation_one <- function(a, trigger, cap, alpha) {
  out <- numeric(length(a))
  keep <- is.finite(a) & a > 0
  if (!any(keep)) {
    return(out)
  }

  ak <- a[keep]
  p0 <- trigger / ak
  p1 <- (trigger + cap) / ak
  full <- p1 <= 1
  values <- numeric(length(ak))
  values[full] <- 1

  partial <- !full
  if (any(partial)) {
    lower <- pmax(1, p0[partial])
    upper <- p1[partial]
    if (abs(alpha - 1) < sqrt(.Machine$double.eps)) {
      first_moment <- log(upper / lower)
    } else {
      first_moment <- alpha *
        (upper^(1 - alpha) - lower^(1 - alpha)) / (1 - alpha)
    }
    probability <- lower^(-alpha) - upper^(-alpha)
    integral <- (ak[partial] / cap) * first_moment -
      (trigger / cap) * probability
    values[partial] <- pmax(0, pmin(1, integral + upper^(-alpha)))
  }

  out[keep] <- values
  out
}

capped_pareto_mean <- function(score_matrix, trigger, cap, alpha) {
  vapply(seq_len(ncol(score_matrix)), function(j) {
    mean(capped_pareto_expectation_one(score_matrix[, j], trigger[j],
                                       cap[j], alpha))
  }, numeric(1L))
}

tail_constant <- function(score_matrix, alpha) {
  colMeans(pmax(score_matrix, 0)^alpha)
}

relative_error <- function(estimate, target) {
  out <- rep(NA_real_, length(target))
  keep <- is.finite(estimate) & is.finite(target) & target > 0
  out[keep] <- abs(estimate[keep] - target[keep]) / target[keep]
  out
}

relative_var_error <- function(log_error) {
  exp(abs(log_error)) - 1
}

summarize_by_rank <- function(x) {
  keys <- paste(x$dataset, x$rank, sep = "\r")
  split_x <- split(x, keys)
  out <- do.call(rbind, lapply(split_x, function(z) {
    data.frame(
      dataset = z$dataset[1L],
      rank = z$rank[1L],
      variation_explained = z$variation_explained[1L],
      residual_risk = z$residual_risk[1L],
      positive_post_projection_delta = z$positive_post_projection_delta[1L],
      capped_error_agca_mean = mean(abs(z$capped_error_agca)),
      capped_error_agca_q90 = unname(quantile(abs(z$capped_error_agca), 0.9)),
      capped_relative_error_agca_mean =
        mean(z$capped_relative_error_agca, na.rm = TRUE),
      capped_relative_error_agca_q90 =
        unname(quantile(z$capped_relative_error_agca, 0.9, na.rm = TRUE)),
      capped_error_post_mean = mean(abs(z$capped_error_post)),
      capped_error_post_q90 = unname(quantile(abs(z$capped_error_post), 0.9)),
      capped_relative_error_post_mean =
        mean(z$capped_relative_error_post, na.rm = TRUE),
      capped_relative_error_post_q90 =
        unname(quantile(z$capped_relative_error_post, 0.9, na.rm = TRUE)),
      var_log_error_agca_mean = mean(abs(z$var_log_error_agca),
                                     na.rm = TRUE),
      var_log_error_agca_q90 = unname(quantile(abs(z$var_log_error_agca),
                                               0.9, na.rm = TRUE)),
      var_relative_error_agca_mean =
        mean(z$var_relative_error_agca, na.rm = TRUE),
      var_relative_error_agca_q90 =
        unname(quantile(z$var_relative_error_agca, 0.9, na.rm = TRUE)),
      var_log_error_post_mean = mean(abs(z$var_log_error_post),
                                     na.rm = TRUE),
      var_log_error_post_q90 = unname(quantile(abs(z$var_log_error_post),
                                               0.9, na.rm = TRUE)),
      var_relative_error_post_mean =
        mean(z$var_relative_error_post, na.rm = TRUE),
      var_relative_error_post_q90 =
        unname(quantile(z$var_relative_error_post, 0.9, na.rm = TRUE))
    )
  }))
  out <- out[order(out$dataset, out$rank), ]
  rownames(out) <- NULL
  out
}

summarize_by_class <- function(x) {
  keys <- paste(x$dataset, x$rank, x$portfolio_class, sep = "\r")
  split_x <- split(x, keys)
  out <- do.call(rbind, lapply(split_x, function(z) {
    data.frame(
      dataset = z$dataset[1L],
      rank = z$rank[1L],
      portfolio_class = z$portfolio_class[1L],
      n_portfolios = nrow(z),
      capped_error_agca_mean = mean(abs(z$capped_error_agca)),
      capped_error_agca_q90 = unname(quantile(abs(z$capped_error_agca), 0.9)),
      capped_relative_error_agca_mean =
        mean(z$capped_relative_error_agca, na.rm = TRUE),
      capped_relative_error_agca_q90 =
        unname(quantile(z$capped_relative_error_agca, 0.9, na.rm = TRUE)),
      capped_error_post_mean = mean(abs(z$capped_error_post)),
      capped_error_post_q90 = unname(quantile(abs(z$capped_error_post), 0.9)),
      capped_relative_error_post_mean =
        mean(z$capped_relative_error_post, na.rm = TRUE),
      capped_relative_error_post_q90 =
        unname(quantile(z$capped_relative_error_post, 0.9, na.rm = TRUE)),
      var_log_error_agca_mean = mean(abs(z$var_log_error_agca),
                                     na.rm = TRUE),
      var_log_error_agca_q90 = unname(quantile(abs(z$var_log_error_agca),
                                               0.9, na.rm = TRUE)),
      var_relative_error_agca_mean =
        mean(z$var_relative_error_agca, na.rm = TRUE),
      var_relative_error_agca_q90 =
        unname(quantile(z$var_relative_error_agca, 0.9, na.rm = TRUE)),
      var_log_error_post_mean = mean(abs(z$var_log_error_post),
                                     na.rm = TRUE),
      var_log_error_post_q90 = unname(quantile(abs(z$var_log_error_post),
                                               0.9, na.rm = TRUE)),
      var_relative_error_post_mean =
        mean(z$var_relative_error_post, na.rm = TRUE),
      var_relative_error_post_q90 =
        unname(quantile(z$var_relative_error_post, 0.9, na.rm = TRUE))
    )
  }))
  out <- out[order(out$dataset, out$rank, out$portfolio_class), ]
  rownames(out) <- NULL
  out
}

write_weight_table <- function(dataset_name, portfolio_columns, portfolio_set) {
  weights <- portfolio_set$weights
  rows <- lapply(seq_len(ncol(weights)), function(j) {
    data.frame(
      dataset = dataset_name,
      portfolio_id = colnames(weights)[j],
      asset = portfolio_columns,
      weight = weights[, j],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

analyze_dataset <- function(dataset_name, spec) {
  if (!file.exists(spec$data_file)) {
    stop("Prepared data not found for ", dataset_name, ": ", spec$data_file,
         call. = FALSE)
  }

  prepared <- readRDS(spec$data_file)
  portfolio_columns <- prepared$portfolio_columns
  loss_data <- prepared$complete_losses
  loss_data <- loss_data[loss_data$date >= sample_start, ]
  loss_matrix <- as.matrix(loss_data[, portfolio_columns, drop = FALSE])
  storage.mode(loss_matrix) <- "double"
  rownames(loss_matrix) <- format(loss_data$date, "%Y-%m-%d")

  pareto_matrix <- pareto_transform(loss_matrix)
  g <- normalize_rows(pareto_matrix)
  radius <- row_norms(pareto_matrix)
  k <- as.integer(round(tail_fraction * nrow(g)))
  index <- select_top_k(radius, k)
  g_extreme <- g[index, , drop = FALSE]

  fit <- agca_fit(g_extreme, mu = canonical_anchor(ncol(g_extreme)))
  rank_summary <- agca_rank_summary(fit)
  max_rank <- max(rank_summary$rank)

  portfolio_set <- make_portfolios(
    portfolio_columns = portfolio_columns,
    metadata = prepared$portfolios,
    group_column = spec$group_column,
    group_label_column = spec$group_label_column,
    n_random = random_portfolios,
    leverage_limit = leverage_limit
  )
  weights <- portfolio_set$weights
  portfolio_info <- portfolio_set$info

  scores_true <- g_extreme %*% weights
  true_constant <- tail_constant(scores_true, alpha)
  l2_norm <- portfolio_info$l2_norm
  trigger <- 2 * l2_norm
  cap <- l2_norm
  true_capped <- capped_pareto_mean(scores_true, trigger, cap, alpha)

  rows <- vector("list", max_rank + 1L)
  for (p in 0L:max_rank) {
    g_hat <- agca_reconstruct(fit, p = p)
    g_post <- positive_post_project(g_hat)
    scores_hat <- g_hat %*% weights
    scores_post <- g_post %*% weights

    constant_hat <- tail_constant(scores_hat, alpha)
    constant_post <- tail_constant(scores_post, alpha)
    capped_hat <- capped_pareto_mean(scores_hat, trigger, cap, alpha)
    capped_post <- capped_pareto_mean(scores_post, trigger, cap, alpha)

    var_log_error_hat <- rep(NA_real_, length(true_constant))
    var_log_error_post <- rep(NA_real_, length(true_constant))
    keep_hat <- true_constant > 0 & constant_hat > 0
    keep_post <- true_constant > 0 & constant_post > 0
    var_log_error_hat[keep_hat] <-
      (log(constant_hat[keep_hat]) - log(true_constant[keep_hat])) / alpha
    var_log_error_post[keep_post] <-
      (log(constant_post[keep_post]) - log(true_constant[keep_post])) / alpha

    capped_relative_error_hat <- relative_error(capped_hat, true_capped)
    capped_relative_error_post <- relative_error(capped_post, true_capped)
    var_relative_error_hat <- relative_var_error(var_log_error_hat)
    var_relative_error_post <- relative_var_error(var_log_error_post)

    rows[[p + 1L]] <- data.frame(
      dataset = dataset_name,
      dataset_label = spec$label,
      rank = p,
      portfolio_info,
      true_capped_excess = true_capped,
      agca_capped_excess = capped_hat,
      post_capped_excess = capped_post,
      capped_error_agca = capped_hat - true_capped,
      capped_error_post = capped_post - true_capped,
      capped_relative_error_agca = capped_relative_error_hat,
      capped_relative_error_post = capped_relative_error_post,
      true_tail_constant = true_constant,
      agca_tail_constant = constant_hat,
      post_tail_constant = constant_post,
      var_log_error_agca = var_log_error_hat,
      var_log_error_post = var_log_error_post,
      var_relative_error_agca = var_relative_error_hat,
      var_relative_error_post = var_relative_error_post,
      variation_explained =
        rank_summary$variation_explained[rank_summary$rank == p],
      residual_risk = rank_summary$residual_risk[rank_summary$rank == p],
      positive_post_projection_delta = mean(rowSums((g_hat - g_post)^2)),
      stringsAsFactors = FALSE
    )
  }

  list(
    errors = do.call(rbind, rows),
    weights = write_weight_table(dataset_name, portfolio_columns, portfolio_set),
    dataset_summary = data.frame(
      dataset = dataset_name,
      dataset_label = spec$label,
      n_complete_days = nrow(g),
      n_extreme_days = length(index),
      n_assets = ncol(g),
      n_portfolios = ncol(weights),
      sample_start = sample_start,
      start_date = min(loss_data$date),
      end_date = max(loss_data$date),
      tail_fraction = tail_fraction,
      alpha = alpha,
      leverage_limit = leverage_limit,
      stringsAsFactors = FALSE
    ),
    g_extreme = g_extreme,
    portfolio_set = portfolio_set,
    trigger = trigger,
    cap = cap
  )
}

analyze_anchor_functionals <- function(dataset_name, result, max_rank) {
  g_extreme <- result$g_extreme
  weights <- result$portfolio_set$weights
  scores_true <- g_extreme %*% weights
  true_constant <- tail_constant(scores_true, alpha)
  true_capped <- capped_pareto_mean(scores_true, result$trigger, result$cap,
                                    alpha)
  anchors <- list(
    canonical = canonical_anchor(ncol(g_extreme)),
    principal = principal_anchor(g_extreme),
    mean = mean_anchor(g_extreme)
  )

  do.call(rbind, lapply(names(anchors), function(anchor) {
    fit <- agca_fit(g_extreme, mu = anchors[[anchor]])
    rank_summary <- agca_rank_summary(fit)
    ranks <- 0L:min(max_rank, max(rank_summary$rank))
    do.call(rbind, lapply(ranks, function(p) {
      g_hat <- agca_reconstruct(fit, p = p)
      scores_hat <- g_hat %*% weights
      constant_hat <- tail_constant(scores_hat, alpha)
      capped_hat <- capped_pareto_mean(scores_hat, result$trigger,
                                       result$cap, alpha)
      var_log_error_hat <- rep(NA_real_, length(true_constant))
      keep <- true_constant > 0 & constant_hat > 0
      var_log_error_hat[keep] <-
        (log(constant_hat[keep]) - log(true_constant[keep])) / alpha
      data.frame(
        dataset = dataset_name,
        anchor = anchor,
        rank = p,
        capped_relative_error_mean =
          mean(relative_error(capped_hat, true_capped), na.rm = TRUE),
        var_relative_error_mean =
          mean(relative_var_error(var_log_error_hat), na.rm = TRUE),
        variation_explained =
          rank_summary$variation_explained[rank_summary$rank == p],
        stringsAsFactors = FALSE
      )
    }))
  }))
}

save_pdf <- function(file, plot_fun, width = 7.2, height = 4.8) {
  pdf(file, width = width, height = height, bg = "white")
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

class_label <- function(x) {
  labels <- c(
    equal = "Equal",
    block_equal = "Block equal",
    random_long_only = "Random long-only",
    random_limited_leverage = "Limited leverage"
  )
  unname(ifelse(x %in% names(labels), labels[x], x))
}

bootstrap_functional_errors <- function(result, dataset_name, max_rank, reps) {
  g_extreme <- result$g_extreme
  weights <- result$portfolio_set$weights
  trigger <- result$trigger
  cap <- result$cap
  n <- nrow(g_extreme)
  rows <- vector("list", reps * (max_rank + 1L))
  row_id <- 1L

  for (b in seq_len(reps)) {
    sample_index <- sample.int(n, n, replace = TRUE)
    g_b <- g_extreme[sample_index, , drop = FALSE]
    fit_b <- agca_fit(g_b, mu = canonical_anchor(ncol(g_b)))
    scores_true <- g_b %*% weights
    true_constant <- tail_constant(scores_true, alpha)
    true_capped <- capped_pareto_mean(scores_true, trigger, cap, alpha)

    for (p in 0L:max_rank) {
      g_hat <- agca_reconstruct(fit_b, p = p)
      scores_hat <- g_hat %*% weights
      constant_hat <- tail_constant(scores_hat, alpha)
      capped_hat <- capped_pareto_mean(scores_hat, trigger, cap, alpha)

      var_log_error <- rep(NA_real_, length(true_constant))
      keep <- true_constant > 0 & constant_hat > 0
      var_log_error[keep] <-
        (log(constant_hat[keep]) - log(true_constant[keep])) / alpha

      rows[[row_id]] <- data.frame(
        dataset = dataset_name,
        bootstrap = b,
        rank = p,
        capped_relative_error_agca_mean =
          mean(relative_error(capped_hat, true_capped), na.rm = TRUE),
        var_relative_error_agca_mean =
          mean(relative_var_error(var_log_error), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }

  do.call(rbind, rows)
}

summarize_functional_bootstrap <- function(x) {
  keys <- paste(x$dataset, x$rank, sep = "\r")
  split_x <- split(x, keys)
  out <- do.call(rbind, lapply(split_x, function(z) {
    data.frame(
      dataset = z$dataset[1L],
      rank = z$rank[1L],
      capped_lower = unname(quantile(z$capped_relative_error_agca_mean,
                                     0.025, na.rm = TRUE)),
      capped_median = unname(quantile(z$capped_relative_error_agca_mean,
                                      0.5, na.rm = TRUE)),
      capped_upper = unname(quantile(z$capped_relative_error_agca_mean,
                                     0.975, na.rm = TRUE)),
      var_lower = unname(quantile(z$var_relative_error_agca_mean,
                                  0.025, na.rm = TRUE)),
      var_median = unname(quantile(z$var_relative_error_agca_mean,
                                   0.5, na.rm = TRUE)),
      var_upper = unname(quantile(z$var_relative_error_agca_mean,
                                  0.975, na.rm = TRUE))
    )
  }))
  out <- out[order(out$dataset, out$rank), ]
  rownames(out) <- NULL
  out
}

plot_rank_errors <- function(rank_summary, y_column_agca, y_column_post,
                             ylab, file) {
  save_pdf(file, function() {
    datasets <- unique(rank_summary$dataset)
    par(mfrow = c(length(datasets), 1L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)
    for (dataset_name in datasets) {
      z <- rank_summary[rank_summary$dataset == dataset_name, ]
      z <- z[z$rank <= dataset_specs[[dataset_name]]$max_plot_rank, ]
      y_max <- max(z[[y_column_agca]], z[[y_column_post]], na.rm = TRUE)
      plot(
        z$rank, z[[y_column_agca]],
        type = "b", pch = 16, col = "#1B6CA8",
        ylim = c(0, y_max),
        xlab = "AGCA rank p",
        ylab = ylab,
        main = dataset_specs[[dataset_name]]$label
      )
      lines(z$rank, z[[y_column_post]], type = "b", pch = 17,
            col = "#B23A48")
      grid(col = "gray90")
      legend("topright", legend = c("AGCA", "Positive post-projection"),
             col = c("#1B6CA8", "#B23A48"), pch = c(16, 17), lty = 1,
             bty = "n", cex = 0.8)
    }
  }, height = 2.8 * length(unique(rank_summary$dataset)))
}

plot_class_errors <- function(class_summary, y_column, ylab, file,
                              selected_ranks = c(3L, 5L, 8L, 10L, 12L)) {
  save_pdf(file, function() {
    datasets <- unique(class_summary$dataset)
    par(mfrow = c(length(datasets), 1L), mar = c(4.2, 4.8, 2.6, 1.0))
    cols <- c("#1B6CA8", "#3B8C5A", "#9C6B1D", "#B23A48")
    for (dataset_name in datasets) {
      z <- class_summary[class_summary$dataset == dataset_name &
                           class_summary$rank %in% selected_ranks, ]
      classes <- unique(z$portfolio_class)
      plot(
        NA,
        xlim = range(selected_ranks),
        ylim = c(0, max(z[[y_column]], na.rm = TRUE)),
        xlab = "AGCA rank p",
        ylab = ylab,
        main = dataset_specs[[dataset_name]]$label
      )
      grid(col = "gray90")
      for (i in seq_along(classes)) {
        zz <- z[z$portfolio_class == classes[i], ]
        zz <- zz[order(zz$rank), ]
        lines(zz$rank, zz[[y_column]], type = "b", pch = 15 + i,
              col = cols[((i - 1L) %% length(cols)) + 1L])
      }
      legend("topright", legend = classes,
             col = cols[seq_along(classes)], pch = 16 + seq_along(classes),
             lty = 1, bty = "n", cex = 0.75)
    }
  }, height = 4.2 * length(unique(class_summary$dataset)))
}

plot_ff_relative_rank_bootstrap <- function(rank_summary, boot_summary, file) {
  save_pdf(file, function() {
    dataset_name <- "ff_2x3_daily"
    max_rank <- dataset_specs[[dataset_name]]$max_plot_rank
    z <- rank_summary[rank_summary$dataset == dataset_name &
                        rank_summary$rank <= max_rank, ]
    b <- boot_summary[boot_summary$dataset == dataset_name &
                        boot_summary$rank <= max_rank, ]
    par(mfrow = c(1L, 2L), mar = c(2.2, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)

    plot_one <- function(y, lower, upper, ylab, main) {
      y <- 100 * y
      lower <- 100 * lower
      upper <- 100 * upper
      ylim <- range(c(0, y, lower, upper), finite = TRUE)
      plot(
        z$rank, y,
        type = "n",
        ylim = ylim,
        xlab = "",
        ylab = ylab,
        main = main
      )
      polygon(c(b$rank, rev(b$rank)), c(lower, rev(upper)),
              border = NA, col = grDevices::adjustcolor("#1B6CA8", 0.18))
      lines(z$rank, y, type = "b", pch = 16, col = "#1B6CA8", lwd = 1.5)
      grid(col = "gray90")
    }

    plot_one(
      z$capped_relative_error_agca_mean,
      b$capped_lower,
      b$capped_upper,
      "Mean relative error (%)",
      "Capped excess"
    )
    plot_one(
      z$var_relative_error_agca_mean,
      b$var_lower,
      b$var_upper,
      "Mean relative error (%)",
      "Normalized VaR"
    )
  }, width = 7.8, height = 2.55)
}

plot_ff_relative_class_errors <- function(class_summary, file) {
  save_pdf(file, function() {
    dataset_name <- "ff_2x3_daily"
    max_rank <- dataset_specs[[dataset_name]]$max_plot_rank
    z <- class_summary[class_summary$dataset == dataset_name &
                         class_summary$rank <= max_rank, ]
    classes <- unique(z$portfolio_class)
    cols <- c("#1B6CA8", "#3B8C5A", "#9C6B1D", "#B23A48")
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 0.8, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)

    plot_one <- function(y_column, ylab, main) {
      ylim <- c(0, 100 * max(z[[y_column]], na.rm = TRUE))
      plot(
        NA,
        xlim = c(0, max_rank),
        ylim = ylim,
        xlab = "AGCA rank p",
        ylab = ylab,
        main = ""
      )
      grid(col = "gray90")
      for (i in seq_along(classes)) {
        zz <- z[z$portfolio_class == classes[i], ]
        zz <- zz[order(zz$rank), ]
        lines(zz$rank, 100 * zz[[y_column]], type = "b", pch = 15 + i,
              col = cols[((i - 1L) %% length(cols)) + 1L])
      }
      legend("topright", legend = class_label(classes),
             col = cols[seq_along(classes)], pch = 16 + seq_along(classes),
             lty = 1, bty = "n", cex = 0.72)
    }

    plot_one("capped_relative_error_agca_mean",
             "Mean relative error (%)", "Capped excess")
    plot_one("var_relative_error_agca_mean",
             "Mean relative error (%)", "Normalized VaR")
  }, width = 7.8, height = 2.55)
}

plot_osap_relative_rank <- function(rank_summary, file) {
  save_pdf(file, function() {
    dataset_name <- "osap_daily_quintile_vw"
    z <- rank_summary[rank_summary$dataset == dataset_name &
                        rank_summary$rank <=
                        dataset_specs[[dataset_name]]$max_plot_rank, ]
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)

    plot_one <- function(y_column, main) {
      plot(
        z$rank, 100 * z[[y_column]],
        type = "b", pch = 16, col = "#1B6CA8",
        ylim = c(0, 100 * max(z[[y_column]], na.rm = TRUE)),
        xlab = "AGCA rank p",
        ylab = "Mean relative error (%)",
        main = main
      )
      grid(col = "gray90")
    }

    plot_one("capped_relative_error_agca_mean", "Capped excess")
    plot_one("var_relative_error_agca_mean", "Normalized VaR")
  }, width = 7.8, height = 2.55)
}

plot_osap_relative_class <- function(class_summary, file,
                                     selected_ranks = c(3L, 5L, 8L, 10L, 12L)) {
  save_pdf(file, function() {
    dataset_name <- "osap_daily_quintile_vw"
    z <- class_summary[class_summary$dataset == dataset_name &
                         class_summary$rank %in% selected_ranks, ]
    classes <- unique(z$portfolio_class)
    cols <- c("#1B6CA8", "#3B8C5A", "#9C6B1D", "#B23A48")
    par(mfrow = c(1L, 2L), mar = c(3.5, 4.2, 2.0, 0.8),
        mgp = c(2.1, 0.65, 0), tcl = -0.25)

    plot_one <- function(y_column, main) {
      plot(
        NA,
        xlim = range(selected_ranks),
        ylim = c(0, 100 * max(z[[y_column]], na.rm = TRUE)),
        xlab = "AGCA rank p",
        ylab = "Mean relative error (%)",
        main = main
      )
      grid(col = "gray90")
      for (i in seq_along(classes)) {
        zz <- z[z$portfolio_class == classes[i], ]
        zz <- zz[order(zz$rank), ]
        lines(zz$rank, 100 * zz[[y_column]], type = "b", pch = 15 + i,
              col = cols[((i - 1L) %% length(cols)) + 1L])
      }
      legend("topright", legend = class_label(classes),
             col = cols[seq_along(classes)], pch = 16 + seq_along(classes),
             lty = 1, bty = "n", cex = 0.72)
    }

    plot_one("capped_relative_error_agca_mean", "Capped excess")
    plot_one("var_relative_error_agca_mean", "Normalized VaR")
  }, width = 7.8, height = 2.55)
}

plot_anchor_functional_sensitivity <- function(anchor_summary, file) {
  save_pdf(file, function() {
    par(mfrow = c(2L, 2L), mar = c(3.3, 4.2, 1.9, 0.8),
        mgp = c(2.0, 0.62, 0), tcl = -0.25)
    cols <- c(canonical = "#1B6CA8", principal = "#B23A48",
              mean = "#3B8C5A")

    plot_one <- function(dataset_name, y_column, main) {
      max_rank <- dataset_specs[[dataset_name]]$max_plot_rank
      z <- anchor_summary[anchor_summary$dataset == dataset_name &
                            anchor_summary$rank <= max_rank, ]
      plot(
        NA,
        xlim = c(0, max_rank),
        ylim = c(0, 100 * max(z[[y_column]], na.rm = TRUE)),
        xlab = "AGCA rank p",
        ylab = "Mean relative error (%)",
        main = main
      )
      grid(col = "gray90")
      for (anchor in names(cols)) {
        zz <- z[z$anchor == anchor, ]
        zz <- zz[order(zz$rank), ]
        lines(zz$rank, 100 * zz[[y_column]], type = "b", pch = 16,
              col = cols[anchor], lwd = 1.35)
      }
      legend("topright", legend = names(cols), col = cols, lty = 1,
             pch = 16, bty = "n", cex = 0.68)
    }

    plot_one("ff_2x3_daily", "capped_relative_error_mean",
             "FF: capped excess")
    plot_one("ff_2x3_daily", "var_relative_error_mean",
             "FF: normalized VaR")
    plot_one("osap_daily_quintile_vw", "capped_relative_error_mean",
             "OSAP: capped excess")
    plot_one("osap_daily_quintile_vw", "var_relative_error_mean",
             "OSAP: normalized VaR")
  }, width = 7.8, height = 5.1)
}

results <- lapply(names(dataset_specs), function(dataset_name) {
  analyze_dataset(dataset_name, dataset_specs[[dataset_name]])
})
names(results) <- names(dataset_specs)

all_errors <- do.call(rbind, lapply(results, `[[`, "errors"))
all_weights <- do.call(rbind, lapply(results, `[[`, "weights"))
dataset_summary <- do.call(rbind, lapply(results, `[[`, "dataset_summary"))
rank_summary <- summarize_by_rank(all_errors)
class_summary <- summarize_by_class(all_errors)
functional_anchor_summary <- do.call(rbind, lapply(names(results),
                                                   function(dataset_name) {
  analyze_anchor_functionals(
    dataset_name,
    results[[dataset_name]],
    max_rank = dataset_specs[[dataset_name]]$max_plot_rank
  )
}))
ff_functional_bootstrap <- bootstrap_functional_errors(
  results$ff_2x3_daily,
  dataset_name = "ff_2x3_daily",
  max_rank = dataset_specs$ff_2x3_daily$max_plot_rank,
  reps = functional_bootstrap_reps
)
ff_functional_bootstrap_summary <-
  summarize_functional_bootstrap(ff_functional_bootstrap)

write.csv(dataset_summary, file.path(output_dir, "dataset_summary.csv"),
          row.names = FALSE)
write.csv(all_weights, file.path(output_dir, "portfolio_weights.csv"),
          row.names = FALSE)
write.csv(all_errors, file.path(output_dir, "portfolio_functional_errors.csv"),
          row.names = FALSE)
write.csv(rank_summary, file.path(output_dir, "functional_rank_summary.csv"),
          row.names = FALSE)
write.csv(class_summary, file.path(output_dir, "functional_class_summary.csv"),
          row.names = FALSE)
write.csv(functional_anchor_summary,
          file.path(output_dir, "functional_anchor_sensitivity.csv"),
          row.names = FALSE)
write.csv(ff_functional_bootstrap,
          file.path(output_dir, "ff_functional_bootstrap.csv"),
          row.names = FALSE)
write.csv(ff_functional_bootstrap_summary,
          file.path(output_dir, "ff_functional_bootstrap_summary.csv"),
          row.names = FALSE)

plot_ff_relative_rank_bootstrap(
  rank_summary,
  ff_functional_bootstrap_summary,
  file.path(output_dir, "ff_relative_error_by_rank_bootstrap.pdf")
)
plot_ff_relative_class_errors(
  class_summary,
  file.path(output_dir, "ff_relative_error_by_class.pdf")
)
plot_osap_relative_rank(
  rank_summary,
  file.path(output_dir, "osap_relative_error_by_rank.pdf")
)
plot_osap_relative_class(
  class_summary,
  file.path(output_dir, "osap_relative_error_by_class.pdf")
)
plot_anchor_functional_sensitivity(
  functional_anchor_summary,
  file.path(output_dir, "functional_anchor_sensitivity.pdf")
)

plot_rank_errors(
  rank_summary,
  "capped_error_agca_mean",
  "capped_error_post_mean",
  "Mean absolute capped-excess error",
  file.path(output_dir, "capped_excess_error_by_rank.pdf")
)
plot_rank_errors(
  rank_summary,
  "var_log_error_agca_mean",
  "var_log_error_post_mean",
  "Mean absolute log-VaR error",
  file.path(output_dir, "var_log_error_by_rank.pdf")
)
plot_class_errors(
  class_summary,
  "capped_error_agca_mean",
  "Mean absolute capped-excess error",
  file.path(output_dir, "capped_excess_error_by_class.pdf")
)
plot_class_errors(
  class_summary,
  "var_log_error_agca_mean",
  "Mean absolute log-VaR error",
  file.path(output_dir, "var_log_error_by_class.pdf")
)
stage_manuscript_outputs()

cat("\nAGCA portfolio tail-function diagnostics\n")
print(dataset_summary, row.names = FALSE)
cat("\nRank summaries, selected ranks:\n")
print(
  rank_summary[rank_summary$rank %in% c(0L, 3L, 5L, 8L, 10L, 12L, 15L),
               c("dataset", "rank", "variation_explained",
                 "capped_error_agca_mean", "capped_error_post_mean",
                 "capped_relative_error_agca_mean",
                 "var_log_error_agca_mean",
                 "var_relative_error_agca_mean",
                 "var_log_error_post_mean")],
  row.names = FALSE
)
cat("\nOutput directory:\n  ", output_dir, "\n", sep = "")
