# Standard 10-variable EVT generator simulation for a full AGCA workflow.
#
# The design combines a shared low-dimensional extremal mechanism with
# variable-specific asymptotically independent extremes. Variables X1-X8 are
# generated from a logistic-block Pareto source whose spectral rays lie in a
# two-dimensional anchored angular subspace. Variables X9-X10 are independent
# Pareto sources, so their extremes create near-axis regimes outside the shared
# low-rank structure.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else if (file.exists(file.path("R", "simulations", "standard_simulations.R"))) {
  file.path("R", "simulations")
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

project_root <- find_compendium_root(script_dir)
source(file.path(project_root, "R", "GeodesicExtreme.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[1L])
}

arg_flag <- function(name) {
  any(args %in% paste0("--", name))
}

parse_integer_list <- function(value, name) {
  parts <- strsplit(value, ",", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  if (length(parts) == 0L || any(!nzchar(parts))) {
    stop(name, " must be a comma-separated list of integers.", call. = FALSE)
  }
  numeric_values <- suppressWarnings(as.numeric(parts))
  values <- as.integer(numeric_values)
  if (any(is.na(values)) || any(values != numeric_values)) {
    stop(name, " must be a comma-separated list of integers.", call. = FALSE)
  }
  unique(values)
}

seed <- as.integer(get_arg("seed", "20260627"))
main_n <- as.integer(get_arg("n", "10000"))
main_k <- as.integer(get_arg("k", "500"))
p_rank <- as.integer(get_arg("rank", "2"))
p_reconstruction <- as.integer(get_arg("reconstruction-rank", "4"))
bootstrap_reps <- as.integer(get_arg("bootstrap-reps", "300"))
population_n <- as.integer(get_arg("population-n", "200000"))
population_seed <- as.integer(get_arg(
  "population-seed",
  as.character(seed + 100000L)
))
coverage_seed <- as.integer(get_arg("coverage-seed", "20260709"))
coverage_reps <- as.integer(get_arg("coverage-reps", "5000"))
coverage_population_n <- as.integer(get_arg("coverage-population-n", "500000"))
coverage_population_seed <- as.integer(get_arg(
  "coverage-population-seed",
  as.character(coverage_seed + 200000L)
))
coverage_calibration_n <- as.integer(get_arg("coverage-calibration-n", "500000"))
coverage_calibration_seed <- as.integer(get_arg(
  "coverage-calibration-seed",
  as.character(coverage_seed + 100000L)
))
coverage_calibration_tail_k_arg <- get_arg("coverage-calibration-tail-k", NA_character_)
coverage_calibration_tail_k <- if (is.na(coverage_calibration_tail_k_arg)) {
  NA_integer_
} else {
  as.integer(coverage_calibration_tail_k_arg)
}
coverage_ranks <- parse_integer_list(get_arg("coverage-ranks", "2"), "--coverage-ranks")
coverage_conf_level <- as.numeric(get_arg("coverage-conf-level", "0.95"))
coverage_checkpoint_every <- as.integer(get_arg("coverage-checkpoint-every", "50"))
coverage_progress_every <- as.integer(get_arg("coverage-progress-every", "25"))
logistic_theta <- as.numeric(get_arg("theta", "0.45"))
finite_tau <- as.numeric(get_arg("tau", "0.25"))
axis9_scale <- as.numeric(get_arg("axis9-scale", "1.00"))
axis10_scale <- as.numeric(get_arg("axis10-scale", "1.00"))
skip_bootstrap <- arg_flag("skip-bootstrap")
skip_coverage <- arg_flag("skip-coverage")
use_rank_transform <- !arg_flag("raw-margins")

if (!is.finite(seed)) {
  stop("--seed must be an integer.", call. = FALSE)
}
if (!is.finite(main_n) || main_n < 100L) {
  stop("--n must be an integer at least 100.", call. = FALSE)
}
if (!is.finite(main_k) || main_k < 5L || main_k >= main_n) {
  stop("--k must be an integer between 5 and n - 1.", call. = FALSE)
}
if (!is.finite(p_rank) || p_rank < 1L || p_rank > 9L) {
  stop("--rank must be an integer between 1 and 9.", call. = FALSE)
}
if (!is.finite(p_reconstruction) || p_reconstruction < p_rank || p_reconstruction > 9L) {
  stop("--reconstruction-rank must be an integer between --rank and 9.", call. = FALSE)
}
if (!is.finite(bootstrap_reps) || bootstrap_reps < 1L) {
  stop("--bootstrap-reps must be a positive integer.", call. = FALSE)
}
if (!is.finite(population_n) || population_n < 100L) {
  stop("--population-n must be an integer at least 100.", call. = FALSE)
}
if (!is.finite(population_seed)) {
  stop("--population-seed must be an integer.", call. = FALSE)
}
if (!is.finite(coverage_seed)) {
  stop("--coverage-seed must be an integer.", call. = FALSE)
}
if (!is.finite(coverage_reps) || coverage_reps < 1L) {
  stop("--coverage-reps must be a positive integer.", call. = FALSE)
}
if (!is.finite(coverage_population_n) || coverage_population_n < 100L) {
  stop("--coverage-population-n must be an integer at least 100.", call. = FALSE)
}
if (!is.finite(coverage_population_seed)) {
  stop("--coverage-population-seed must be an integer.", call. = FALSE)
}
if (!is.finite(coverage_calibration_n) || coverage_calibration_n < 100L) {
  stop("--coverage-calibration-n must be an integer at least 100.", call. = FALSE)
}
if (!is.na(coverage_calibration_tail_k) &&
    (!is.finite(coverage_calibration_tail_k) ||
     coverage_calibration_tail_k < 10L ||
     coverage_calibration_tail_k >= coverage_calibration_n)) {
  stop(
    "--coverage-calibration-tail-k must be between 10 and ",
    "coverage-calibration-n - 1.",
    call. = FALSE
  )
}
if (!is.finite(coverage_calibration_seed)) {
  stop("--coverage-calibration-seed must be an integer.", call. = FALSE)
}
if (any(coverage_ranks < 1L) || any(coverage_ranks > 8L)) {
  stop("--coverage-ranks must contain integers between 1 and 8.", call. = FALSE)
}
if (!is.finite(coverage_conf_level) ||
    coverage_conf_level <= 0 ||
    coverage_conf_level >= 1) {
  stop("--coverage-conf-level must lie in (0, 1).", call. = FALSE)
}
if (!is.finite(coverage_checkpoint_every) || coverage_checkpoint_every < 0L) {
  stop("--coverage-checkpoint-every must be a nonnegative integer.", call. = FALSE)
}
if (!is.finite(coverage_progress_every) || coverage_progress_every < 0L) {
  stop("--coverage-progress-every must be a nonnegative integer.", call. = FALSE)
}
if (!is.finite(logistic_theta) || logistic_theta <= 0 || logistic_theta >= 1) {
  stop("--theta must lie in (0, 1).", call. = FALSE)
}
if (!is.finite(finite_tau) || finite_tau < 0) {
  stop("--tau must be nonnegative.", call. = FALSE)
}
if (!is.finite(axis9_scale) || axis9_scale <= 0 ||
    !is.finite(axis10_scale) || axis10_scale <= 0) {
  stop("--axis9-scale and --axis10-scale must be positive.", call. = FALSE)
}

set.seed(seed)

output_dir <- file.path(
  project_root,
  "inst",
  "simulations",
  "results",
  "standard_simulation_output"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

variable_labels <- paste0("X", seq_len(10L))
regime_levels <- c("shared_low_rank", "axis_9", "axis_10")
x9_color <- "#5B3A29"
x10_color <- "#7570B3"
population_loading_color <- "#D95F02"
regime_cols <- c(
  shared_low_rank = "#1B9E77",
  axis_9 = x9_color,
  axis_10 = x10_color
)

threshold_k_values <- unique(pmin(
  main_n - 1L,
  pmax(5L, as.integer(round(main_k * c(0.5, 0.7, 1, 1.5, 2))))
))

rpareto <- function(n, shape = 1, scale = 1) {
  scale * runif(n)^(-1 / shape)
}

rpositive_stable <- function(n, alpha) {
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0, 1).", call. = FALSE)
  }
  u <- runif(n, min = 0, max = pi)
  w <- rexp(n)
  sin(alpha * u) / (sin(u)^(1 / alpha)) *
    (sin((1 - alpha) * u) / w)^((1 - alpha) / alpha)
}

rlogistic_frechet <- function(n, d, theta) {
  stable <- rpositive_stable(n, theta)
  expo <- matrix(rexp(n * d), n, d)
  sweep(1 / expo, 1L, stable, "*")^theta
}

frechet_to_pareto <- function(z) {
  1 / (-expm1(-1 / z))
}

rlogistic_pareto <- function(n, d, theta) {
  frechet_to_pareto(rlogistic_frechet(n, d, theta))
}

rank_pareto_transform <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  out <- apply(x, 2L, function(z) {
    r <- rank(z, ties.method = "average")
    (n + 1) / (n + 1 - r)
  })
  colnames(out) <- colnames(x)
  out
}

shared_contrast_basis <- function() {
  raw <- cbind(
    gradient = seq(-1, 1, length.out = 8L),
    block = c(rep(1, 4L), rep(-1, 4L))
  )
  raw <- scale(raw, center = TRUE, scale = FALSE)
  orthonormalize_columns(raw)
}

shared_rays <- function() {
  mu8 <- canonical_anchor(8L)
  basis <- shared_contrast_basis()
  scores <- rbind(
    c(-0.42, -0.26),
    c(-0.30, 0.30),
    c(0.00, -0.36),
    c(0.28, 0.24),
    c(0.48, -0.02)
  )
  rays <- matrix(mu8, nrow(scores), 8L, byrow = TRUE) +
    scores %*% t(basis)
  if (any(rays <= 0)) {
    stop("Shared low-rank rays left the positive orthant.", call. = FALSE)
  }
  normalize_rows(rays)
}

add_finite_threshold_noise <- function(x, tau = finite_tau) {
  x <- as.matrix(x)
  if (tau == 0) {
    return(x)
  }
  x + tau * matrix(rexp(length(x)), nrow(x), ncol(x))
}

simulate_standard_10d <- function(n) {
  rays8 <- shared_rays()
  z_shared <- rlogistic_pareto(n, d = nrow(rays8), theta = logistic_theta)
  shared_signal <- z_shared %*% rays8
  axis9_signal <- axis9_scale * rpareto(n)
  axis10_signal <- axis10_scale * rpareto(n)

  x_signal <- cbind(shared_signal, axis9_signal, axis10_signal)
  x <- add_finite_threshold_noise(x_signal)
  colnames(x) <- variable_labels
  colnames(x_signal) <- variable_labels

  regime_score <- cbind(
    shared_low_rank = row_norms(shared_signal),
    axis_9 = axis9_signal,
    axis_10 = axis10_signal
  )
  label <- colnames(regime_score)[max.col(regime_score, ties.method = "first")]

  list(
    x = x,
    signal = x_signal,
    latent_shared = z_shared,
    latent_axis9 = axis9_signal,
    latent_axis10 = axis10_signal,
    label = factor(label, levels = regime_levels),
    regime_score = regime_score
  )
}

axis_contrast <- function(d, axis, mu) {
  e <- rep(0, d)
  e[axis] <- 1
  u <- drop(project_to_tangent(matrix(e, nrow = 1L), mu))
  unit_vector(u, "axis contrast")
}

true_contrasts_full <- function(mu) {
  shared_basis <- shared_contrast_basis()
  embedded <- rbind(shared_basis, matrix(0, 2L, ncol(shared_basis)))
  out <- target_contrasts_for_anchor(mu, embedded)
  colnames(out) <- c("shared AGC 1", "shared AGC 2")
  cbind(
    out,
    axis9 = axis_contrast(10L, axis = 9L, mu = mu),
    axis10 = axis_contrast(10L, axis = 10L, mu = mu)
  )
}

sphere_log_map <- function(mu, g) {
  mu <- unit_vector(mu, "mu")
  g <- normalize_rows(g)
  inner <- pmax(-1, pmin(1, drop(g %*% mu)))
  theta <- acos(inner)
  residual <- g - tcrossprod(inner, mu)
  factor <- rep(1, length(theta))
  regular <- theta > 1e-10
  factor[regular] <- theta[regular] / sin(theta[regular])
  sweep(residual, 1L, factor, "*")
}

sphere_exp_map <- function(mu, v) {
  mu <- unit_vector(mu, "mu")
  v <- as.numeric(v)
  nrm <- sqrt(sum(v^2))
  if (!is.finite(nrm) || nrm <= 1e-12) {
    return(mu)
  }
  cos(nrm) * mu + sin(nrm) * v / nrm
}

spherical_frechet_anchor <- function(g, max_iter = 200L, tol = 1e-10) {
  g <- normalize_rows(g)
  mu <- unit_vector(colMeans(g), "initial Frechet anchor")
  for (iter in seq_len(max_iter)) {
    update <- colMeans(sphere_log_map(mu, g))
    update_norm <- sqrt(sum(update^2))
    if (!is.finite(update_norm) || update_norm <= tol) {
      break
    }
    step <- 1
    repeat {
      candidate <- sphere_exp_map(mu, step * update)
      if (all(candidate > 0) || step < 1e-6) {
        mu <- unit_vector(pmax(candidate, 1e-12), "Frechet anchor")
        break
      }
      step <- step / 2
    }
  }
  unit_vector(mu, "Frechet anchor")
}

principal_anchor <- function(g) {
  g <- normalize_rows(g)
  moment <- crossprod(g) / nrow(g)
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

target_contrasts_for_anchor <- function(mu, contrasts) {
  contrasts <- as.matrix(contrasts)
  projected <- apply(contrasts, 2L, function(x) {
    unit_vector(
      project_to_tangent(matrix(x, nrow = 1L), mu),
      "projected contrast"
    )
  })
  if (is.null(dim(projected))) {
    projected <- matrix(projected, ncol = 1L)
  }
  colnames(projected) <- colnames(contrasts)
  projected
}

with_seed <- function(seed_value, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  set.seed(seed_value)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

align_loading_signs <- function(loadings, reference) {
  loadings <- as.matrix(loadings)
  reference <- as.matrix(reference)
  n_components <- min(ncol(loadings), ncol(reference))
  for (j in seq_len(n_components)) {
    if (drop(crossprod(loadings[, j], reference[, j])) < 0) {
      loadings[, j] <- -loadings[, j]
    }
  }
  loadings
}

fit_at_threshold <- function(x, labels, k, mu, p, reference_space = NULL,
                             target_contrasts = NULL) {
  thr <- threshold_directions(x, k = k)
  fit <- agca_fit(thr$g, mu = mu, p = p)
  variation <- agca_variation_explained(fit)
  selected_labels <- labels[thr$index]
  label_table <- prop.table(table(factor(selected_labels, levels = regime_levels)))

  projector_distance <- if (is.null(reference_space)) {
    NA_real_
  } else {
    subspace_distance(fit$loadings[, seq_len(p), drop = FALSE], reference_space)
  }
  target_alignment <- if (is.null(target_contrasts)) {
    NA_real_
  } else {
    mean(principal_angle_cosines(
      fit$loadings[, seq_len(p), drop = FALSE],
      target_contrasts[, seq_len(p), drop = FALSE]
    ))
  }

  eigen_cols <- as.data.frame(as.list(fit$eigenvalues))
  names(eigen_cols) <- paste0("eig", seq_along(fit$eigenvalues))

  summary <- cbind(
    data.frame(k = k, threshold = thr$threshold),
    eigen_cols,
    data.frame(
      variation_explained_rank = variation[p],
      residual_risk_rank = agca_mean_residual(fit, p = p),
      projector_distance_to_reference = projector_distance,
      target_space_mean_cosine = target_alignment,
      shared_fraction = unname(label_table[["shared_low_rank"]]),
      axis9_fraction = unname(label_table[["axis_9"]]),
      axis10_fraction = unname(label_table[["axis_10"]])
    )
  )

  list(
    fit = fit,
    threshold = thr,
    selected_labels = selected_labels,
    summary = summary
  )
}

threshold_path <- function(x, labels, k_values, mu, p, reference_fit,
                           target_contrasts = NULL) {
  reference_space <- reference_fit$loadings[, seq_len(p), drop = FALSE]
  rows <- lapply(k_values, function(k) {
    fit_at_threshold(
      x = x,
      labels = labels,
      k = k,
      mu = mu,
      p = p,
      reference_space = reference_space,
      target_contrasts = target_contrasts
    )$summary
  })
  do.call(rbind, rows)
}

eigen_summary <- function(fit) {
  data.frame(
    component = seq_along(fit$eigenvalues),
    eigenvalue = fit$eigenvalues,
    cumulative_variation = agca_variation_explained(fit)
  )
}

anchor_sensitivity <- function(g_selected, canonical_fit, p) {
  canonical_mu <- canonical_fit$mu
  canonical_model <- cbind(
    canonical_mu,
    canonical_fit$loadings[, seq_len(p), drop = FALSE]
  )
  anchors <- list(
    canonical = canonical_mu,
    Frechet = spherical_frechet_anchor(g_selected),
    principal = principal_anchor(g_selected)
  )

  rows <- lapply(names(anchors), function(anchor_name) {
    anchor <- anchors[[anchor_name]]
    current_fit <- agca_fit(g_selected, mu = anchor, p = p)
    current_model <- cbind(
      anchor,
      current_fit$loadings[, seq_len(p), drop = FALSE]
    )
    data.frame(
      anchor = anchor_name,
      anchor_distance_to_canonical = sphere_geodesic_distance(
        matrix(anchor, nrow = 1L),
        matrix(canonical_mu, nrow = 1L)
      ),
      total_anchored_variation = sum(current_fit$eigenvalues),
      variation_explained_rank = agca_variation_explained(current_fit)[p],
      residual_risk_rank = agca_mean_residual(current_fit, p = p),
      model_projector_distance_to_canonical = subspace_distance(
        current_model,
        canonical_model
      )
    )
  })
  do.call(rbind, rows)
}

bootstrap_agca_stability <- function(g_selected, fit, p, labels,
                                     b = bootstrap_reps,
                                     loading_components = p_reconstruction,
                                     seed_value = seed + 200000L) {
  with_seed(seed_value, {
    n <- nrow(g_selected)
    d <- ncol(g_selected)
    loading_components <- min(loading_components, ncol(fit$loadings))
    reference_space <- fit$loadings[, seq_len(p), drop = FALSE]
    boot_loadings <- array(
      NA_real_,
      dim = c(b, d, loading_components),
      dimnames = list(NULL, variable_labels, paste0("component", seq_len(loading_components)))
    )
    variation_curves <- matrix(
      NA_real_,
      nrow = b,
      ncol = d,
      dimnames = list(NULL, paste0("rank", 0:(d - 1L)))
    )
    rows <- data.frame(
      iteration = seq_len(b),
      variation_explained_rank = NA_real_,
      residual_risk_rank = NA_real_,
      projector_distance_to_main = NA_real_,
      loading1_abs_alignment = NA_real_,
      loading2_abs_alignment = NA_real_,
      shared_fraction = NA_real_,
      axis9_fraction = NA_real_,
      axis10_fraction = NA_real_
    )

    for (iter in seq_len(b)) {
      sample_index <- sample.int(n, n, replace = TRUE)
      boot_fit <- agca_fit(g_selected[sample_index, , drop = FALSE],
                           mu = fit$mu, p = p)
      rows$variation_explained_rank[iter] <- agca_variation_explained(boot_fit)[p]
      rows$residual_risk_rank[iter] <- agca_mean_residual(boot_fit, p = p)
      variation_curves[iter, ] <- c(0, agca_variation_explained(boot_fit))
      rows$projector_distance_to_main[iter] <- subspace_distance(
        boot_fit$loadings[, seq_len(p), drop = FALSE],
        reference_space
      )
      boot_labels <- labels[sample_index]
      label_table <- prop.table(table(factor(boot_labels, levels = regime_levels)))
      rows$shared_fraction[iter] <- unname(label_table[["shared_low_rank"]])
      rows$axis9_fraction[iter] <- unname(label_table[["axis_9"]])
      rows$axis10_fraction[iter] <- unname(label_table[["axis_10"]])

      for (component in seq_len(loading_components)) {
        loading <- boot_fit$loadings[, component]
        reference <- fit$loadings[, component]
        alignment <- sum(loading * reference)
        if (alignment < 0) {
          loading <- -loading
          alignment <- -alignment
        }
        boot_loadings[iter, , component] <- loading
        if (component == 1L) {
          rows$loading1_abs_alignment[iter] <- alignment
        } else if (component == 2L) {
          rows$loading2_abs_alignment[iter] <- alignment
        }
      }
    }

    loading_rows <- do.call(
      rbind,
      lapply(seq_len(loading_components), function(component) {
        values <- boot_loadings[, , component, drop = TRUE]
        data.frame(
          component = component,
          variable = variable_labels,
          main_loading = fit$loadings[, component],
          mean = colMeans(values, na.rm = TRUE),
          sd = apply(values, 2L, sd, na.rm = TRUE),
          q025 = apply(values, 2L, quantile, probs = 0.025, na.rm = TRUE),
          q50 = apply(values, 2L, quantile, probs = 0.50, na.rm = TRUE),
          q975 = apply(values, 2L, quantile, probs = 0.975, na.rm = TRUE),
          row.names = NULL
        )
      })
    )
    variation_intervals <- data.frame(
      rank = 0:(d - 1L),
      main_variation = c(0, agca_variation_explained(fit)),
      mean = colMeans(variation_curves, na.rm = TRUE),
      sd = apply(variation_curves, 2L, sd, na.rm = TRUE),
      q025 = apply(variation_curves, 2L, quantile, probs = 0.025, na.rm = TRUE),
      q50 = apply(variation_curves, 2L, quantile, probs = 0.50, na.rm = TRUE),
      q975 = apply(variation_curves, 2L, quantile, probs = 0.975, na.rm = TRUE),
      row.names = NULL
    )

    list(
      iterations = rows,
      loading_intervals = loading_rows,
      variation_intervals = variation_intervals
    )
  })
}

bootstrap_summary <- function(boot_iterations, main_values) {
  rows <- lapply(names(main_values), function(metric) {
    values <- boot_iterations[[metric]]
    data.frame(
      metric = metric,
      main_value = main_values[[metric]],
      mean = mean(values, na.rm = TRUE),
      sd = sd(values, na.rm = TRUE),
      q025 = quantile(values, probs = 0.025, na.rm = TRUE),
      q50 = quantile(values, probs = 0.50, na.rm = TRUE),
      q975 = quantile(values, probs = 0.975, na.rm = TRUE),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

oracle_alignment <- function(fit, target_contrasts, n_components = 6L) {
  n_components <- min(n_components, ncol(fit$loadings))
  align <- abs(t(fit$loadings[, seq_len(n_components), drop = FALSE]) %*%
                 target_contrasts)
  rownames(align) <- as.character(seq_len(n_components))
  out <- as.data.frame(as.table(align))
  names(out) <- c("component", "target", "abs_alignment")
  out$component <- as.integer(as.character(out$component))
  out
}

estimate_population_agca <- function(mu, p, n = population_n,
                                     seed_value = population_seed,
                                     tail_fraction = main_k / main_n) {
  with_seed(seed_value, {
    observations <- simulate_standard_10d(n)
    x_analysis <- if (use_rank_transform) {
      rank_pareto_transform(observations$x)
    } else {
      observations$x
    }
    colnames(x_analysis) <- variable_labels
    k <- pmin(n - 1L, pmax(5L, as.integer(round(n * tail_fraction))))
    threshold <- threshold_directions(x_analysis, k = k)
    fit <- agca_fit(threshold$g, mu = mu, p = p)
    list(
      fit = fit,
      threshold = threshold,
      n = n,
      k = k,
      tail_fraction = k / n
    )
  })
}

population_loading_table <- function(population_fit) {
  do.call(
    rbind,
    lapply(seq_len(ncol(population_fit$loadings)), function(component) {
      data.frame(
        component = component,
        variable = variable_labels,
        loading = population_fit$loadings[, component],
        eigenvalue = population_fit$eigenvalues[component],
        cumulative_variation = agca_variation_explained(population_fit)[component],
        row.names = NULL
      )
    })
  )
}

matched_population_agca_loadings <- function(fit, population_fit,
                                             components = seq_len(4L)) {
  components <- components[components <= ncol(fit$loadings)]
  population_components <- seq_len(ncol(population_fit$loadings))
  align <- abs(t(fit$loadings[, components, drop = FALSE]) %*%
                 population_fit$loadings[, population_components, drop = FALSE])

  selected <- integer(0)
  population_index <- integer(length(components))
  for (row in seq_along(components)) {
    order_j <- order(align[row, ], decreasing = TRUE)
    choice <- order_j[!(order_j %in% selected)][1L]
    if (is.na(choice)) {
      choice <- order_j[1L]
    }
    population_index[row] <- choice
    selected <- c(selected, choice)
  }

  out <- population_fit$loadings[, population_index, drop = FALSE]
  for (j in seq_along(components)) {
    if (drop(crossprod(fit$loadings[, components[j]], out[, j])) < 0) {
      out[, j] <- -out[, j]
    }
  }
  colnames(out) <- as.character(components)

  list(
    loadings = out,
    match_table = data.frame(
      component = components,
      population_component = population_components[population_index],
      abs_alignment = align[cbind(seq_along(components), population_index)],
      row.names = NULL
    )
  )
}

shared_block_oracle_fit <- function(x_analysis, labels, selected_index) {
  shared_index <- selected_index[labels[selected_index] == "shared_low_rank"]
  if (length(shared_index) < 10L) {
    return(NULL)
  }
  g8 <- normalize_rows(x_analysis[shared_index, seq_len(8L), drop = FALSE])
  mu8 <- canonical_anchor(8L)
  fit8 <- agca_fit(g8, mu = mu8, p = 2L)
  target8 <- shared_contrast_basis()
  list(
    fit = fit8,
    index = shared_index,
    summary = data.frame(
      n_shared_selected = length(shared_index),
      variation_explained_rank2 = agca_variation_explained(fit8)[2L],
      residual_risk_rank2 = agca_mean_residual(fit8, p = 2L),
      projector_distance_to_truth = subspace_distance(
        fit8$loadings[, seq_len(2L), drop = FALSE],
        target8
      )
    )
  )
}

empirical_sd <- function(x) {
  x <- as.numeric(x)
  sqrt(mean((x - mean(x))^2))
}

build_oracle_margin_transform <- function(calibration_x, tail_k = NA_integer_) {
  calibration_x <- as.matrix(calibration_x)
  n <- nrow(calibration_x)
  d <- ncol(calibration_x)
  if (is.na(tail_k)) {
    tail_k <- min(n - 1L, max(1000L, as.integer(round(0.01 * n))))
  }
  if (tail_k < 10L || tail_k >= n) {
    stop("tail_k must be between 10 and n - 1.", call. = FALSE)
  }

  sorted <- lapply(seq_len(d), function(j) sort(calibration_x[, j]))
  threshold <- vapply(sorted, function(z) z[n - tail_k + 1L], numeric(1L))
  tail_probability <- tail_k / (n + 1)
  tail_constant <- threshold * tail_probability

  transform <- function(x) {
    x <- as.matrix(x)
    if (ncol(x) != d) {
      stop("x has incompatible dimension for oracle margin transform.", call. = FALSE)
    }
    out <- matrix(NA_real_, nrow(x), d)
    for (j in seq_len(d)) {
      rank_leq <- findInterval(x[, j], sorted[[j]], rightmost.closed = TRUE)
      empirical_survival <- (n + 1 - rank_leq) / (n + 1)
      tail_survival <- tail_constant[j] / pmax(x[, j], .Machine$double.xmin)
      survival <- empirical_survival
      use_tail <- x[, j] > threshold[j]
      survival[use_tail] <- tail_survival[use_tail]
      survival <- pmin(1, pmax(survival, .Machine$double.xmin))
      out[, j] <- 1 / survival
    }
    colnames(out) <- colnames(x)
    out
  }

  list(
    transform = transform,
    summary = data.frame(
      variable = colnames(calibration_x),
      calibration_n = n,
      tail_k = tail_k,
      threshold = threshold,
      tail_probability = tail_probability,
      tail_constant = tail_constant,
      row.names = NULL
    )
  )
}

coverage_point_estimates <- function(fit, ranks = coverage_ranks) {
  max_rank <- max(ranks)
  eigen_components <- seq_len(max_rank)
  eigen_values <- fit$eigenvalues[eigen_components]
  names(eigen_values) <- paste0("lambda", eigen_components)
  ave_values <- vapply(
    ranks,
    function(p) agca_variation_explained(fit)[p],
    numeric(1L)
  )
  names(ave_values) <- paste0("AVE", ranks)
  c(
    tau = sum(fit$eigenvalues),
    eigen_values,
    ave_values
  )
}

coverage_fit_margin_sample <- function(x, k, mu, p) {
  threshold <- threshold_directions(x, k = k)
  fit <- agca_fit(threshold$g, mu = mu, p = p)
  list(threshold = threshold, fit = fit)
}

coverage_interval_rows <- function(fit, rep_id, margin, targets, z_crit,
                                   ranks = coverage_ranks) {
  u <- fit$u
  k <- nrow(u)
  norm2 <- rowSums(u^2)
  tau_hat <- sum(fit$eigenvalues)
  variation_explained <- agca_variation_explained(fit)
  max_rank <- max(ranks)

  make_row <- function(statistic, estimate, influence_values) {
    sigma_hat <- empirical_sd(influence_values)
    se <- sigma_hat / sqrt(k)
    lower <- estimate - z_crit * se
    upper <- estimate + z_crit * se
    target <- unname(targets[[statistic]])
    data.frame(
      replicate = rep_id,
      margin = margin,
      statistic = statistic,
      k = k,
      estimate = estimate,
      target = target,
      se = se,
      lower = lower,
      upper = upper,
      interval_length = upper - lower,
      covered = lower <= target && target <= upper,
      studentized = if (is.finite(se) && se > 0) {
        (estimate - target) / se
      } else {
        NA_real_
      },
      row.names = NULL
    )
  }

  rows <- list(make_row("tau", tau_hat, norm2))
  for (j in seq_len(max_rank)) {
    scores_j <- drop(u %*% fit$loadings[, j])
    rows[[length(rows) + 1L]] <- make_row(
      paste0("lambda", j),
      fit$eigenvalues[j],
      scores_j^2
    )
  }

  for (p in ranks) {
    ave_hat <- variation_explained[p]
    projected_norm2 <- rowSums((u %*% fit$loadings[, seq_len(p), drop = FALSE])^2)
    psi <- (projected_norm2 - ave_hat * norm2) / tau_hat
    rows[[length(rows) + 1L]] <- make_row(paste0("AVE", p), ave_hat, psi)
  }

  do.call(rbind, rows)
}

coverage_summary <- function(intervals) {
  pieces <- split(intervals, list(intervals$margin, intervals$statistic), drop = TRUE)
  rows <- lapply(pieces, function(z) {
    valid <- is.finite(z$lower) & is.finite(z$upper) & is.finite(z$target)
    covered <- z$covered[valid]
    estimates <- z$estimate[valid]
    ses <- z$se[valid]
    studentized <- z$studentized[valid]
    coverage <- mean(covered)
    data.frame(
      margin = z$margin[1L],
      statistic = z$statistic[1L],
      target = z$target[1L],
      reps = nrow(z),
      valid_reps = sum(valid),
      coverage = coverage,
      coverage_mcse = sqrt(coverage * (1 - coverage) / sum(valid)),
      mean_estimate = mean(estimates),
      bias = mean(estimates - z$target[valid]),
      empirical_sd = sd(estimates),
      mean_se = mean(ses),
      se_to_empirical_sd = mean(ses) / sd(estimates),
      mean_interval_length = mean(z$interval_length[valid]),
      studentized_mean = mean(studentized, na.rm = TRUE),
      studentized_sd = sd(studentized, na.rm = TRUE),
      row.names = NULL
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$statistic, out$margin), , drop = FALSE]
}

run_oracle_coverage_simulation <- function(file_dir) {
  z_crit <- qnorm((1 + coverage_conf_level) / 2)
  tail_fraction <- main_k / main_n
  max_coverage_rank <- max(coverage_ranks)
  coverage_population_k <- pmin(
    coverage_population_n - 1L,
    pmax(5L, as.integer(round(coverage_population_n * tail_fraction)))
  )

  files <- c(
    metadata = file.path(file_dir, "oracle_coverage_metadata.csv"),
    calibration = file.path(file_dir, "oracle_coverage_margin_calibration.csv"),
    population_targets = file.path(file_dir, "oracle_coverage_population_targets.csv"),
    intervals = file.path(file_dir, "oracle_coverage_replicate_intervals.csv"),
    summary = file.path(file_dir, "oracle_coverage_summary.csv")
  )

  metadata <- data.frame(
    coverage_seed = coverage_seed,
    coverage_reps = coverage_reps,
    n = main_n,
    k = main_k,
    main_rank = p_rank,
    coverage_ranks = paste(coverage_ranks, collapse = ","),
    tail_fraction = tail_fraction,
    population_n = coverage_population_n,
    population_k = coverage_population_k,
    population_seed = coverage_population_seed,
    calibration_n = coverage_calibration_n,
    calibration_seed = coverage_calibration_seed,
    calibration_tail_k = if (is.na(coverage_calibration_tail_k)) {
      NA_integer_
    } else {
      coverage_calibration_tail_k
    },
    logistic_theta = logistic_theta,
    finite_tau = finite_tau,
    axis9_scale = axis9_scale,
    axis10_scale = axis10_scale,
    conf_level = coverage_conf_level,
    z_crit = z_crit,
    oracle_margin_method = paste(
      "independent calibration empirical CDF with Pareto tail extrapolation"
    )
  )
  write.csv(metadata, files[["metadata"]], row.names = FALSE)

  message("Calibrating oracle margins with n = ", coverage_calibration_n)
  calibration <- with_seed(coverage_calibration_seed, simulate_standard_10d(coverage_calibration_n))
  oracle_transform <- build_oracle_margin_transform(
    calibration$x,
    tail_k = coverage_calibration_tail_k
  )
  write.csv(oracle_transform$summary, files[["calibration"]], row.names = FALSE)

  mu <- canonical_anchor(10L)
  message(
    "Estimating oracle coverage population target with n = ",
    coverage_population_n,
    " and k = ",
    coverage_population_k
  )
  population <- with_seed(coverage_population_seed, simulate_standard_10d(coverage_population_n))
  population_oracle_x <- oracle_transform$transform(population$x)
  population_fit <- coverage_fit_margin_sample(
    population_oracle_x,
    k = coverage_population_k,
    mu = mu,
    p = max_coverage_rank
  )$fit
  targets <- coverage_point_estimates(population_fit, ranks = coverage_ranks)
  write.csv(
    data.frame(statistic = names(targets), target = unname(targets), row.names = NULL),
    files[["population_targets"]],
    row.names = FALSE
  )

  message(
    "Running ",
    coverage_reps,
    " oracle coverage replicates with n = ",
    main_n,
    " and k = ",
    main_k
  )
  set.seed(coverage_seed)
  interval_list <- vector("list", coverage_reps)
  for (rep_id in seq_len(coverage_reps)) {
    observations <- simulate_standard_10d(main_n)
    oracle_x <- oracle_transform$transform(observations$x)
    rank_x <- rank_pareto_transform(observations$x)

    oracle_fit <- coverage_fit_margin_sample(
      oracle_x,
      k = main_k,
      mu = mu,
      p = max_coverage_rank
    )$fit
    rank_fit <- coverage_fit_margin_sample(
      rank_x,
      k = main_k,
      mu = mu,
      p = max_coverage_rank
    )$fit

    interval_list[[rep_id]] <- rbind(
      coverage_interval_rows(
        oracle_fit,
        rep_id = rep_id,
        margin = "oracle",
        targets = targets,
        z_crit = z_crit,
        ranks = coverage_ranks
      ),
      coverage_interval_rows(
        rank_fit,
        rep_id = rep_id,
        margin = "rank",
        targets = targets,
        z_crit = z_crit,
        ranks = coverage_ranks
      )
    )

    if (coverage_progress_every > 0L && rep_id %% coverage_progress_every == 0L) {
      message("Completed oracle coverage replicate ", rep_id, " / ", coverage_reps)
    }
    if (coverage_checkpoint_every > 0L && rep_id %% coverage_checkpoint_every == 0L) {
      write.csv(
        do.call(rbind, interval_list[seq_len(rep_id)]),
        files[["intervals"]],
        row.names = FALSE
      )
    }
  }

  intervals <- do.call(rbind, interval_list)
  summary <- coverage_summary(intervals)
  write.csv(intervals, files[["intervals"]], row.names = FALSE)
  write.csv(summary, files[["summary"]], row.names = FALSE)

  list(
    targets = targets,
    intervals = intervals,
    summary = summary,
    files = files
  )
}

save_pdf <- function(file, plot_fun, width = 6.5, height = 5.5) {
  pdf(file, width = width, height = height, bg = "white")
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

plot_text_cex <- list(axis = 1.05, lab = 1.15, main = 1.05, legend = 1.00)

plot_eigenvalues <- function(fit, file) {
  save_pdf(file, function() {
    par(mar = c(4.8, 5.0, 3.2, 1.2), cex.axis = plot_text_cex$axis,
        cex.lab = plot_text_cex$lab, cex.main = plot_text_cex$main)
    barplot(
      fit$eigenvalues,
      names.arg = seq_along(fit$eigenvalues),
      xlab = "Component",
      ylab = "Eigenvalue",
      main = "AGCA eigenvalues",
      col = "gray70",
      border = "gray35"
    )
  })
}

plot_variation_explained <- function(fits, file, text_scale = 1,
                                     title_scale = text_scale,
                                     label_scale = text_scale) {
  save_pdf(file, function() {
    par(mar = c(4.8, 5.0, 3.4, 1.2),
        cex.axis = plot_text_cex$axis * text_scale,
        cex.lab = plot_text_cex$lab * label_scale,
        cex.main = plot_text_cex$main * title_scale)
    if (inherits(fits, "agca_fit")) {
      fits <- list(canonical = fits)
    }
    max_rank <- max(vapply(fits, function(fit) length(fit$eigenvalues), integer(1L)))
    rank <- 0L:max_rank
    plot(
      rank,
      rep(NA_real_, length(rank)),
      type = "n",
      ylim = c(0, 1),
      xaxt = "n",
      xlab = "Rank p",
      ylab = "Cumulative anchored variation",
      main = "Variation explained by anchor"
    )
    axis(1, at = rank, labels = rank)
    abline(h = c(0.8, 0.9, 0.95), col = "gray85", lty = 2)
    anchor_cols <- c(canonical = "black", Frechet = "#2166AC", principal = "#B2182B")
    anchor_lty <- c(canonical = 1, Frechet = 2, principal = 3)
    anchor_pch <- c(canonical = 16, Frechet = 17, principal = 15)
    for (anchor_name in names(fits)) {
      variation <- c(0, agca_variation_explained(fits[[anchor_name]]))
      lines(
        seq_along(variation) - 1L,
        variation,
        type = "b",
        pch = anchor_pch[[anchor_name]],
        lty = anchor_lty[[anchor_name]],
        col = anchor_cols[[anchor_name]],
        lwd = 2
      )
    }
    legend(
      "bottomright",
      legend = names(fits),
      col = anchor_cols[names(fits)],
      lty = anchor_lty[names(fits)],
      pch = anchor_pch[names(fits)],
      lwd = 2,
      bty = "n",
      cex = plot_text_cex$legend * text_scale
    )
  })
}

plot_canonical_variation_ci <- function(variation_intervals, file) {
  save_pdf(file, function() {
    par(mar = c(4.8, 5.0, 3.2, 1.2), cex.axis = plot_text_cex$axis,
        cex.lab = plot_text_cex$lab, cex.main = plot_text_cex$main)
    plot(
      variation_intervals$rank,
      variation_intervals$main_variation,
      type = "n",
      ylim = c(0, 1),
      xlim = c(min(variation_intervals$rank), max(variation_intervals$rank) + 0.25),
      xaxt = "n",
      xlab = "Rank p",
      ylab = "Cumulative anchored variation",
      main = "Canonical variation explained"
    )
    axis(1, at = variation_intervals$rank, labels = variation_intervals$rank)
    abline(h = c(0.8, 0.9, 0.95), col = "gray85", lty = 2)
    polygon(
      c(variation_intervals$rank, rev(variation_intervals$rank)),
      c(variation_intervals$q025, rev(variation_intervals$q975)),
      col = "gray88",
      border = NA
    )
    ci_x <- variation_intervals$rank + 0.10
    segments(
      ci_x,
      variation_intervals$q025,
      ci_x,
      variation_intervals$q975,
      col = "gray45",
      lwd = 2.0
    )
    segments(
      ci_x - 0.055,
      variation_intervals$q025,
      ci_x + 0.055,
      variation_intervals$q025,
      col = "gray45",
      lwd = 2.0
    )
    segments(
      ci_x - 0.055,
      variation_intervals$q975,
      ci_x + 0.055,
      variation_intervals$q975,
      col = "gray45",
      lwd = 2.0
    )
    lines(
      variation_intervals$rank,
      variation_intervals$main_variation,
      type = "b",
      pch = 16,
      lty = 1,
      col = "black",
      lwd = 2
    )
    legend(
      "bottomright",
      legend = c("canonical", "95% bootstrap CI"),
      col = c("black", "gray45"),
      lty = 1,
      pch = c(16, NA),
      lwd = c(2, 4),
      bty = "n",
      cex = plot_text_cex$legend
    )
  }, width = 6.2, height = 5.0)
}

plot_scores <- function(scores, file, title = "First two AGCA scores") {
  save_pdf(file, function() {
    par(mar = c(4.8, 5.0, 3.2, 1.2), cex.axis = plot_text_cex$axis,
        cex.lab = plot_text_cex$lab, cex.main = plot_text_cex$main)
    cols <- adjustcolor(regime_cols[as.character(scores$label)], alpha.f = 0.70)
    xlim <- extendrange(c(scores$score1, 0), f = 0.18)
    ylim <- extendrange(c(scores$score2, 0), f = 0.18)
    plot(
      scores$score1,
      scores$score2,
      pch = 16,
      cex = 0.65,
      col = cols,
      xlim = xlim,
      ylim = ylim,
      xlab = "AGCA score 1",
      ylab = "AGCA score 2",
      main = title
    )
    abline(h = 0, v = 0, col = "gray80")
    points(0, 0, pch = 4, lwd = 1.4, cex = 1.1)
    legend(
      "topright",
      legend = names(regime_cols),
      col = regime_cols,
      pch = 16,
      bty = "n",
      cex = plot_text_cex$legend
    )
  })
}

plot_loading <- function(fit, file, component, loading_intervals = NULL,
                         population_loading = NULL,
                         text_scale = 1,
                         title_scale = text_scale) {
  save_pdf(file, function() {
    par(mar = c(5.8, 4.8, 3.5, 1.0),
        cex.axis = 0.95 * text_scale,
        cex.lab = 1.05 * text_scale,
        cex.main = 1.00 * title_scale)
    cols <- rep("gray70", length(variable_labels))
    cols[9:10] <- c(x9_color, x10_color)
    current_intervals <- if (is.null(loading_intervals)) {
      NULL
    } else {
      loading_intervals[loading_intervals$component == component, , drop = FALSE]
    }

    loading <- fit$loadings[, component]
    display_sign <- 1
    if (!is.null(population_loading) &&
        drop(crossprod(loading, population_loading)) < 0) {
      display_sign <- -1
    }
    loading <- display_sign * loading

    if (!is.null(current_intervals) && nrow(current_intervals) > 0L) {
      current_intervals$main_loading <- display_sign * current_intervals$main_loading
      if (display_sign < 0) {
        q025 <- -current_intervals$q975
        q975 <- -current_intervals$q025
        current_intervals$q025 <- q025
        current_intervals$q975 <- q975
      } else {
        current_intervals$q025 <- display_sign * current_intervals$q025
        current_intervals$q975 <- display_sign * current_intervals$q975
      }
    }

    ylim_values <- c(loading, 0)
    if (!is.null(current_intervals) && nrow(current_intervals) > 0L) {
      ylim_values <- c(ylim_values, current_intervals$q025, current_intervals$q975)
    }
    if (!is.null(population_loading)) {
      ylim_values <- c(ylim_values, population_loading)
    }
    ylim <- extendrange(ylim_values, f = 0.25)

    centers <- barplot(
      loading,
      names.arg = variable_labels,
      las = 2,
      ylim = ylim,
      col = cols,
      border = "gray35",
      ylab = "Loading",
      main = paste0("AGC ", component)
    )
    abline(h = 0, col = "gray65")

    legend_items <- character(0)
    legend_lty <- numeric(0)
    legend_lwd <- numeric(0)
    legend_pch <- numeric(0)
    legend_col <- character(0)

    if (!is.null(current_intervals) && nrow(current_intervals) == length(variable_labels)) {
      segments(
        centers,
        current_intervals$q025,
        centers,
        current_intervals$q975,
        lwd = 1.2,
        col = "black"
      )
      segments(
        centers - 0.055,
        current_intervals$q025,
        centers + 0.055,
        current_intervals$q025,
        lwd = 1.2,
        col = "black"
      )
      segments(
        centers - 0.055,
        current_intervals$q975,
        centers + 0.055,
        current_intervals$q975,
        lwd = 1.2,
        col = "black"
      )
      points(centers, current_intervals$main_loading, pch = 16, cex = 0.55)
      legend_items <- c(legend_items, "bootstrap 95% CI")
      legend_lty <- c(legend_lty, 1)
      legend_lwd <- c(legend_lwd, 1.2)
      legend_pch <- c(legend_pch, NA)
      legend_col <- c(legend_col, "black")
    }

    if (!is.null(population_loading)) {
      points(centers, population_loading, pch = 16, cex = 0.85,
             col = population_loading_color)
      legend_items <- c(legend_items, "population loading")
      legend_lty <- c(legend_lty, NA)
      legend_lwd <- c(legend_lwd, NA)
      legend_pch <- c(legend_pch, 16)
      legend_col <- c(legend_col, population_loading_color)
    }

    if (length(legend_items) > 0L) {
      legend(
        "top",
        legend = legend_items,
        lty = legend_lty,
        lwd = legend_lwd,
        pch = legend_pch,
        col = legend_col,
        bty = "n",
        horiz = TRUE,
        x.intersp = 0.60,
        cex = 0.78 * text_scale
      )
    }
  }, width = 6.5, height = 5.0)
}

plot_loading_heatmap <- function(fit, file, n_components = 6L) {
  save_pdf(file, function() {
    n_components <- min(n_components, ncol(fit$loadings))
    par(mar = c(5.8, 5.0, 5.4, 1.0), cex.axis = 1.0,
        cex.lab = 1.1, cex.main = 1.0)
    z <- fit$loadings[, seq_len(n_components), drop = FALSE]
    max_abs <- max(abs(z))
    breaks <- seq(-max_abs, max_abs, length.out = 101L)
    cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100L)
    image(
      x = seq_len(nrow(z)),
      y = seq_len(ncol(z)),
      z = z,
      breaks = breaks,
      col = cols,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "Component",
      main = ""
    )
    axis(1, at = seq_len(nrow(z)), labels = variable_labels, las = 2)
    axis(2, at = seq_len(ncol(z)), labels = paste0("AGC ", seq_len(ncol(z))))
    mtext("AGCA loading heatmap", side = 3, line = 3.2, font = 2)
    usr <- par("usr")
    y_top <- usr[4L] + 0.14 * diff(usr[3:4])
    y_bottom <- usr[4L] + 0.085 * diff(usr[3:4])
    x_left <- usr[1L] + 0.22 * diff(usr[1:2])
    x_right <- usr[1L] + 0.78 * diff(usr[1:2])
    legend_x <- seq(x_left, x_right, length.out = length(cols) + 1L)
    for (i in seq_along(cols)) {
      rect(
        legend_x[i],
        y_bottom,
        legend_x[i + 1L],
        y_top,
        col = cols[i],
        border = NA,
        xpd = NA
      )
    }
    text(x_left, y_top + 0.018 * diff(usr[3:4]), "- loading", adj = c(0, 0),
         cex = 0.78, xpd = NA)
    text((x_left + x_right) / 2, y_top + 0.018 * diff(usr[3:4]), "0",
         cex = 0.78, xpd = NA)
    text(x_right, y_top + 0.018 * diff(usr[3:4]), "+ loading", adj = c(1, 0),
         cex = 0.78, xpd = NA)
  }, width = 8.0, height = 5.8)
}

plot_threshold_sensitivity <- function(path, file) {
  save_pdf(file, function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4.8, 5.0, 3.1, 1.0),
        cex.axis = 0.95, cex.lab = 1.05, cex.main = 0.98)
    plot(path$k, path$variation_explained_rank, type = "b", pch = 16,
         ylim = c(0, 1), xlab = "Number of extremes k",
         ylab = paste0("Rank-", p_rank, " variation"),
         main = "Explained variation")
    plot(path$k, path$residual_risk_rank, type = "b", pch = 16,
         xlab = "Number of extremes k", ylab = "Residual risk",
         main = "Reconstruction risk")
    plot(path$k, path$projector_distance_to_reference, type = "b", pch = 16,
         ylim = c(0, max(path$projector_distance_to_reference, na.rm = TRUE) * 1.1),
         xlab = "Number of extremes k", ylab = "Projector distance",
         main = "Subspace stability")
    matplot(
      path$k,
      as.matrix(path[, c("shared_fraction", "axis9_fraction", "axis10_fraction")]),
      type = "b",
      pch = 16,
      lty = 1,
      col = regime_cols,
      ylim = c(0, 1),
      xlab = "Number of extremes k",
      ylab = "Selected fraction",
      main = "Selected regimes"
    )
    legend("topright", legend = names(regime_cols), col = regime_cols,
           lty = 1, pch = 16, bty = "n", cex = 0.75)
  }, width = 9.0, height = 7.2)
}

plot_anchor_sensitivity_metric <- function(anchor_summary, metric, file,
                                           ylab, main, ylim = NULL,
                                           bar_width = 0.65,
                                           bar_space = 0,
                                           width = 3.6,
                                           height = 4.8,
                                           text_scale = 1,
                                           title_scale = text_scale,
                                           label_scale = text_scale,
                                           x_label_angle = 90) {
  save_pdf(file, function() {
    bottom_margin <- if (x_label_angle == 90) 5.8 else 4.9
    par(mar = c(bottom_margin, 4.8, 3.4, 1.0),
        cex.axis = 0.95 * text_scale,
        cex.lab = 1.05 * label_scale,
        cex.main = 1.00 * title_scale,
        xpd = NA)
    values <- anchor_summary[[metric]]
    if (is.null(ylim)) {
      ylim <- extendrange(c(0, values), f = 0.10)
    }
    centers <- barplot(
      values,
      names.arg = if (x_label_angle == 90) anchor_summary$anchor else FALSE,
      las = 2,
      ylim = ylim,
      width = bar_width,
      space = bar_space,
      col = c("gray70", "#2166AC", "#B2182B")[seq_len(nrow(anchor_summary))],
      border = "gray35",
      ylab = ylab,
      main = main
    )
    if (x_label_angle != 90) {
      axis(1, at = centers, labels = FALSE)
      text(
        centers,
        par("usr")[3L] - 0.035 * diff(par("usr")[3:4]),
        labels = anchor_summary$anchor,
        srt = x_label_angle,
        adj = 1,
        cex = 0.95 * text_scale
      )
    }
  }, width = width, height = height)
}

plot_oracle_alignment <- function(alignment, file) {
  save_pdf(file, function() {
    par(mar = c(5.8, 8.4, 3.2, 1.0), cex.axis = 0.82,
        cex.lab = 1.05, cex.main = 1.00)
    components <- sort(unique(alignment$component))
    targets <- unique(alignment$target)
    z <- xtabs(abs_alignment ~ component + target, data = alignment)
    z <- z[as.character(components), targets, drop = FALSE]
    rownames(z) <- paste0("AGC ", components)
    image(
      x = seq_len(nrow(z)),
      y = seq_len(ncol(z)),
      z = z,
      col = colorRampPalette(c("white", "#2166AC"))(100L),
      zlim = c(0, 1),
      xaxt = "n",
      yaxt = "n",
      xlab = "Estimated component",
      ylab = "",
      main = "Oracle loading alignment"
    )
    axis(1, at = seq_len(nrow(z)), labels = rownames(z))
    axis(2, at = seq_len(ncol(z)), labels = colnames(z), las = 1,
         cex.axis = 0.78)
    mtext("Population contrast", side = 2, line = 6.4, cex = 1.05)
    for (i in seq_len(nrow(z))) {
      for (j in seq_len(ncol(z))) {
        text(i, j, labels = sprintf("%.2f", z[i, j]), cex = 0.75)
      }
    }
  }, width = 8.0, height = 5.8)
}

write_design_metadata <- function(file) {
  metadata <- data.frame(
    seed = seed,
    n = main_n,
    k = main_k,
    rank = p_rank,
    reconstruction_rank = p_reconstruction,
    threshold_path_k = paste(threshold_k_values, collapse = ","),
    bootstrap_reps = if (skip_bootstrap) 0L else bootstrap_reps,
    population_n = population_n,
    population_seed = population_seed,
    population_tail_fraction = main_k / main_n,
    coverage_reps = if (skip_coverage) 0L else coverage_reps,
    coverage_seed = coverage_seed,
    coverage_ranks = paste(coverage_ranks, collapse = ","),
    logistic_theta = logistic_theta,
    finite_tau = finite_tau,
    axis9_scale = axis9_scale,
    axis10_scale = axis10_scale,
    rank_pareto_transform = use_rank_transform
  )
  write.csv(metadata, file, row.names = FALSE)
}

run_standard_10d_simulation <- function() {
  metadata_file <- file.path(output_dir, "design_metadata.csv")
  write_design_metadata(metadata_file)

  observations <- simulate_standard_10d(main_n)
  x_analysis <- if (use_rank_transform) {
    rank_pareto_transform(observations$x)
  } else {
    observations$x
  }
  colnames(x_analysis) <- variable_labels

  mu <- canonical_anchor(10L)
  target_contrasts <- true_contrasts_full(mu)

  main <- fit_at_threshold(
    x = x_analysis,
    labels = observations$label,
    k = main_k,
    mu = mu,
    p = p_rank,
    target_contrasts = target_contrasts
  )
  fit <- main$fit
  selected_labels <- main$selected_labels
  selected_table <- table(factor(selected_labels, levels = regime_levels))

  print_fit_summary <- function(title, fit, p) {
    eigen_table <- data.frame(
      component = seq_along(fit$eigenvalues),
      eigenvalue = fit$eigenvalues,
      cumulative_variation = agca_variation_explained(fit)
    )
    cat("\n", title, "\n", sep = "")
    print(round(eigen_table, 4), row.names = FALSE)
    cat("\nResidual risk by reconstruction rank:\n")
    print(round(agca_rank_summary(fit), 4), row.names = FALSE)
    cat("\nRank ", p, " variation explained: ",
        round(agca_variation_explained(fit)[p], 4), "\n", sep = "")
  }

  print_fit_summary("Standard 10D simulation: full AGCA fit", fit, p_rank)
  cat("Selected regimes:\n")
  print(selected_table)

  threshold <- main$threshold
  path <- threshold_path(
    x = x_analysis,
    labels = observations$label,
    k_values = threshold_k_values,
    mu = mu,
    p = p_rank,
    reference_fit = fit,
    target_contrasts = target_contrasts
  )

  anchor_summary <- anchor_sensitivity(
    g_selected = threshold$g,
    canonical_fit = fit,
    p = p_rank
  )

  anchors <- list(
    canonical = mu,
    sample_frechet = spherical_frechet_anchor(threshold$g),
    sample_principal = principal_anchor(threshold$g)
  )
  anchor_coordinates <- do.call(
    rbind,
    lapply(names(anchors), function(anchor_name) {
      data.frame(
        anchor = anchor_name,
        coordinate = seq_along(anchors[[anchor_name]]),
        variable = variable_labels,
        value = anchors[[anchor_name]]
      )
    })
  )
  anchor_fits <- list(
    canonical = fit,
    Frechet = agca_fit(threshold$g, mu = anchors$sample_frechet, p = p_rank),
    principal = agca_fit(threshold$g, mu = anchors$sample_principal, p = p_rank)
  )

  scores <- as.data.frame(fit$scores)
  names(scores) <- paste0("score", seq_len(ncol(scores)))
  scores$index <- threshold$index
  scores$radius <- threshold$radius
  scores$label <- selected_labels
  scores <- scores[, c("index", "radius", "label",
                       paste0("score", seq_len(ncol(fit$scores))))]

  rank_summary <- agca_rank_summary(fit)
  eigen_table <- eigen_summary(fit)
  oracle <- oracle_alignment(fit, target_contrasts, n_components = 6L)
  message(
    "Estimating Monte Carlo population AGCA with n = ",
    population_n,
    " and tail fraction = ",
    signif(main_k / main_n, 4)
  )
  population_agca <- estimate_population_agca(mu = mu, p = p_rank)
  population_loading_matches <- matched_population_agca_loadings(
    fit,
    population_agca$fit,
    components = seq_len(min(4L, ncol(fit$loadings)))
  )
  shared_oracle <- shared_block_oracle_fit(
    x_analysis = x_analysis,
    labels = observations$label,
    selected_index = threshold$index
  )

  bootstrap <- NULL
  bootstrap_metrics <- NULL
  if (!skip_bootstrap) {
    message("Running bootstrap stability with ", bootstrap_reps, " replicates")
    bootstrap <- bootstrap_agca_stability(
      g_selected = threshold$g,
      fit = fit,
      p = p_rank,
      labels = selected_labels,
      b = bootstrap_reps
    )
    bootstrap_metrics <- bootstrap_summary(
      bootstrap$iterations,
      main_values = c(
        variation_explained_rank = agca_variation_explained(fit)[p_rank],
        residual_risk_rank = agca_mean_residual(fit, p = p_rank),
        projector_distance_to_main = 0,
        loading1_abs_alignment = 1,
        loading2_abs_alignment = if (p_rank >= 2L) 1 else NA_real_
      )
    )
  }

  coverage <- NULL
  if (!skip_coverage) {
    coverage <- run_oracle_coverage_simulation(output_dir)
  }

  files <- c(
    metadata = metadata_file,
    raw_observations = file.path(output_dir, "raw_observations.csv"),
    analyzed_observations = file.path(output_dir, "analyzed_observations.csv"),
    main_summary = file.path(output_dir, "main_threshold_summary.csv"),
    eigen_summary = file.path(output_dir, "eigen_summary.csv"),
    rank_summary = file.path(output_dir, "rank_summary.csv"),
    scores = file.path(output_dir, "scores.csv"),
    selected_regime_counts = file.path(output_dir, "selected_regime_counts.csv"),
    threshold_sensitivity = file.path(output_dir, "threshold_sensitivity.csv"),
    anchor_sensitivity = file.path(output_dir, "anchor_sensitivity.csv"),
    anchor_coordinates = file.path(output_dir, "anchor_coordinates.csv"),
    population_threshold_summary = file.path(output_dir, "population_threshold_summary.csv"),
    population_eigen_summary = file.path(output_dir, "population_eigen_summary.csv"),
    population_loadings = file.path(output_dir, "population_loadings.csv"),
    population_loading_matches = file.path(output_dir, "population_loading_matches.csv"),
    oracle_alignment = file.path(output_dir, "oracle_alignment.csv"),
    shared_block_oracle = file.path(output_dir, "shared_block_oracle.csv")
  )

  write.csv(observations$x, files[["raw_observations"]], row.names = FALSE)
  write.csv(x_analysis, files[["analyzed_observations"]], row.names = FALSE)
  write.csv(main$summary, files[["main_summary"]], row.names = FALSE)
  write.csv(eigen_table, files[["eigen_summary"]], row.names = FALSE)
  write.csv(rank_summary, files[["rank_summary"]], row.names = FALSE)
  write.csv(scores, files[["scores"]], row.names = FALSE)
  write.csv(as.data.frame(selected_table),
            files[["selected_regime_counts"]], row.names = FALSE)
  write.csv(path, files[["threshold_sensitivity"]], row.names = FALSE)
  write.csv(anchor_summary, files[["anchor_sensitivity"]], row.names = FALSE)
  write.csv(anchor_coordinates, files[["anchor_coordinates"]], row.names = FALSE)
  write.csv(
    data.frame(
      n = population_agca$n,
      k = population_agca$k,
      tail_fraction = population_agca$tail_fraction,
      threshold = population_agca$threshold$threshold
    ),
    files[["population_threshold_summary"]],
    row.names = FALSE
  )
  write.csv(
    eigen_summary(population_agca$fit),
    files[["population_eigen_summary"]],
    row.names = FALSE
  )
  write.csv(
    population_loading_table(population_agca$fit),
    files[["population_loadings"]],
    row.names = FALSE
  )
  write.csv(
    population_loading_matches$match_table,
    files[["population_loading_matches"]],
    row.names = FALSE
  )
  write.csv(oracle, files[["oracle_alignment"]], row.names = FALSE)
  if (is.null(shared_oracle)) {
    write.csv(data.frame(), files[["shared_block_oracle"]], row.names = FALSE)
  } else {
    write.csv(shared_oracle$summary, files[["shared_block_oracle"]],
              row.names = FALSE)
    write.csv(eigen_summary(shared_oracle$fit),
              file.path(output_dir, "shared_block_eigen_summary.csv"),
              row.names = FALSE)
    write.csv(agca_rank_summary(shared_oracle$fit),
              file.path(output_dir, "shared_block_rank_summary.csv"),
              row.names = FALSE)
  }

  if (!is.null(bootstrap)) {
    files <- c(
      files,
      bootstrap_stability = file.path(output_dir, "bootstrap_stability.csv"),
      bootstrap_summary = file.path(output_dir, "bootstrap_summary.csv"),
      bootstrap_loading_intervals = file.path(output_dir, "bootstrap_loading_intervals.csv")
    )
    write.csv(bootstrap$iterations, files[["bootstrap_stability"]],
              row.names = FALSE)
    write.csv(bootstrap_metrics, files[["bootstrap_summary"]],
              row.names = FALSE)
    write.csv(bootstrap$loading_intervals,
              files[["bootstrap_loading_intervals"]], row.names = FALSE)
  }
  if (!is.null(coverage)) {
    files <- c(
      files,
      oracle_coverage_metadata = coverage$files[["metadata"]],
      oracle_coverage_margin_calibration = coverage$files[["calibration"]],
      oracle_coverage_population_targets = coverage$files[["population_targets"]],
      oracle_coverage_replicate_intervals = coverage$files[["intervals"]],
      oracle_coverage_summary = coverage$files[["summary"]]
    )
  }

  plot_files <- c(
    eigenvalues = file.path(output_dir, "eigenvalues.pdf"),
    variation_explained = file.path(output_dir, "anchored_variation_explained.pdf"),
    scores = file.path(output_dir, "scores_rank2.pdf"),
    loading_agc1 = file.path(output_dir, "loading_agc1.pdf"),
    loading_agc2 = file.path(output_dir, "loading_agc2.pdf"),
    loading_agc3 = file.path(output_dir, "loading_agc3.pdf"),
    loading_agc4 = file.path(output_dir, "loading_agc4.pdf"),
    loading_heatmap = file.path(output_dir, "loadings_heatmap.pdf"),
    threshold_sensitivity = file.path(output_dir, "threshold_sensitivity.pdf"),
    anchor_distance = file.path(output_dir, "anchor_sensitivity_distance.pdf"),
    oracle_alignment = file.path(output_dir, "oracle_alignment.pdf")
  )
  if (!is.null(bootstrap)) {
    plot_files <- c(
      plot_files,
      canonical_variation_ci = file.path(output_dir, "canonical_variation_explained_ci.pdf")
    )
  }

  plot_eigenvalues(fit, plot_files[["eigenvalues"]])
  plot_variation_explained(
    anchor_fits,
    plot_files[["variation_explained"]],
    text_scale = 1.18,
    title_scale = 1.32,
    label_scale = 1.28
  )
  if (!is.null(bootstrap)) {
    plot_canonical_variation_ci(
      bootstrap$variation_intervals,
      plot_files[["canonical_variation_ci"]]
    )
  }
  plot_scores(scores, plot_files[["scores"]])
  for (component in seq_len(min(4L, ncol(fit$loadings)))) {
    plot_loading(
      fit,
      plot_files[[paste0("loading_agc", component)]],
      component = component,
      loading_intervals = if (is.null(bootstrap)) {
        NULL
      } else {
        bootstrap$loading_intervals
      },
      population_loading = if (as.character(component) %in%
                               colnames(population_loading_matches$loadings)) {
        population_loading_matches$loadings[, as.character(component)]
      } else {
        NULL
      },
      text_scale = if (component <= 2L) 1.22 else 1,
      title_scale = if (component <= 2L) 1.45 else 1
    )
  }
  plot_loading_heatmap(fit, plot_files[["loading_heatmap"]])
  plot_threshold_sensitivity(path, plot_files[["threshold_sensitivity"]])
  plot_anchor_sensitivity_metric(
    anchor_summary,
    metric = "anchor_distance_to_canonical",
    file = plot_files[["anchor_distance"]],
    ylab = "Geodesic distance",
    main = "Anchor distance",
    text_scale = 1.18,
    title_scale = 1.32,
    label_scale = 1.28,
    x_label_angle = 35
  )
  plot_oracle_alignment(oracle, plot_files[["oracle_alignment"]])

  cat("\nStandard 10D simulation output directory:\n")
  cat("  ", output_dir, "\n", sep = "")
  cat("\nKey outputs:\n")
  cat("  ", files[["main_summary"]], "\n", sep = "")
  cat("  ", files[["threshold_sensitivity"]], "\n", sep = "")
  cat("  ", files[["anchor_sensitivity"]], "\n", sep = "")
  if (!is.null(bootstrap)) {
    cat("  ", files[["bootstrap_summary"]], "\n", sep = "")
  }
  if (!is.null(coverage)) {
    cat("  ", files[["oracle_coverage_summary"]], "\n", sep = "")
  }
  cat("  ", plot_files[["scores"]], "\n", sep = "")
  cat("  ", plot_files[["oracle_alignment"]], "\n", sep = "")

  invisible(list(
    observations = observations,
    analyzed = x_analysis,
    main = main,
    threshold_path = path,
    anchor_sensitivity = anchor_summary,
    anchor_coordinates = anchor_coordinates,
    oracle_alignment = oracle,
    shared_block_oracle = shared_oracle,
    bootstrap = bootstrap,
    coverage = coverage,
    files = files,
    plot_files = plot_files
  ))
}

standard_simulation_results <- run_standard_10d_simulation()
