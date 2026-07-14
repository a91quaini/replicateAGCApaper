# Standard EVT generator simulations for three-variable AGCA diagnostics.
#
# This script mirrors the output structure of simulations_3d.R and
# simulations_3d_n_sensitivity.R, but replaces the hand-built angular laws by
# logistic-block heavy-tailed generators motivated by Gnecco et al., Example 3.
#
# Model 1: a three-dimensional logistic tail block is embedded so that two
#          dominant rays generate span(mu0, b1), while a weaker third ray bends
#          the selected angular cloud through the b2 direction.
# Model 2: a bivariate logistic X1-X2 tail block is combined with an independent
#          Pareto X3 source, giving one asymptotically independent variable.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else if (file.exists(file.path("R", "simulations", "standard_simulations_3d.R"))) {
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

n_rep <- as.integer(get_arg("reps", "300"))
bootstrap_reps <- as.integer(get_arg("bootstrap-reps", "300"))
seed <- as.integer(get_arg("seed", "20260627"))
skip_sensitivity <- arg_flag("skip-sensitivity")
skip_bootstrap <- arg_flag("skip-bootstrap")

if (!is.finite(n_rep) || n_rep < 1L) {
  stop("--reps must be a positive integer.", call. = FALSE)
}
if (!is.finite(bootstrap_reps) || bootstrap_reps < 1L) {
  stop("--bootstrap-reps must be a positive integer.", call. = FALSE)
}
if (!is.finite(seed)) {
  stop("--seed must be an integer.", call. = FALSE)
}

set.seed(seed)

output_dir <- file.path(
  project_root,
  "inst",
  "simulations",
  "results",
  "standard_simulation_output_3d"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theta_arg <- get_arg("theta", NA_character_)
logistic_theta_model1 <- as.numeric(get_arg(
  "theta-model1",
  if (is.na(theta_arg)) "0.50" else theta_arg
))
logistic_theta_model2 <- as.numeric(get_arg(
  "theta-model2",
  if (is.na(theta_arg)) "0.25" else theta_arg
))
embedding_rho <- as.numeric(get_arg("rho", "0.25"))
model1_third_eta <- as.numeric(get_arg("model1-third-eta", "0.45"))
model1_third_scale <- as.numeric(get_arg("model1-third-scale", "0.45"))
finite_tau <- as.numeric(get_arg("tau", "0.30"))
axis_scale <- as.numeric(get_arg("axis-scale", "1.00"))
adaptive_anchor_type <- get_arg("adaptive-anchor", "frechet")
main_n <- as.integer(get_arg("main-n", "2400"))
main_k <- as.integer(get_arg("main-k", "120"))
population_anchor_n <- as.integer(get_arg("population-anchor-n", "200000"))
population_anchor_seed <- as.integer(get_arg(
  "population-anchor-seed",
  as.character(seed + 100000L)
))

if (!is.finite(logistic_theta_model1) ||
    logistic_theta_model1 <= 0 ||
    logistic_theta_model1 >= 1) {
  stop("--theta-model1 must lie in (0, 1).", call. = FALSE)
}
if (!is.finite(logistic_theta_model2) ||
    logistic_theta_model2 <= 0 ||
    logistic_theta_model2 >= 1) {
  stop("--theta-model2 must lie in (0, 1).", call. = FALSE)
}
if (!is.finite(embedding_rho) || embedding_rho <= 0 || embedding_rho >= 1) {
  stop("--rho must lie in (0, 1).", call. = FALSE)
}
if (!is.finite(model1_third_eta) || model1_third_eta <= 0) {
  stop("--model1-third-eta must be positive.", call. = FALSE)
}
if (!is.finite(model1_third_scale) || model1_third_scale <= 0) {
  stop("--model1-third-scale must be positive.", call. = FALSE)
}
if (!is.finite(finite_tau) || finite_tau < 0) {
  stop("--tau must be nonnegative.", call. = FALSE)
}
if (!is.finite(axis_scale) || axis_scale <= 0) {
  stop("--axis-scale must be positive.", call. = FALSE)
}
if (!(adaptive_anchor_type %in% c("frechet", "principal"))) {
  stop("--adaptive-anchor must be either 'frechet' or 'principal'.", call. = FALSE)
}
if (!is.finite(main_n) || main_n < 2L) {
  stop("--main-n must be an integer at least 2.", call. = FALSE)
}
if (!is.finite(main_k) || main_k < 1L || main_k >= main_n) {
  stop("--main-k must be an integer between 1 and --main-n - 1.", call. = FALSE)
}
if (!is.finite(population_anchor_n) || population_anchor_n < 2L) {
  stop("--population-anchor-n must be an integer at least 2.", call. = FALSE)
}
if (!is.finite(population_anchor_seed)) {
  stop("--population-anchor-seed must be an integer.", call. = FALSE)
}

main_k_values <- unique(pmin(
  main_n - 1L,
  pmax(1L, as.integer(round(main_k * c(0.25, 0.5, 1, 1.5, 2.5))))
))

contrast_basis_d3 <- function() {
  cbind(
    b1 = c(1, -1, 0) / sqrt(2),
    b2 = c(1, 1, -2) / sqrt(6)
  )
}

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

standard_embedding_matrix <- function(rho = embedding_rho,
                                      third_eta = model1_third_eta,
                                      third_scale = model1_third_scale) {
  gamma <- (1 + rho) / 2
  mu <- canonical_anchor(3L)
  b2 <- contrast_basis_d3()[, "b2"]
  a3 <- mu - third_eta * b2
  if (any(a3 <= 0)) {
    stop(
      "The Model 1 weak third ray left the positive orthant. ",
      "Reduce --model1-third-eta.",
      call. = FALSE
    )
  }

  cbind(
    a1 = c(1, rho, gamma),
    a2 = c(rho, 1, gamma),
    a3 = third_scale * a3
  )
}

add_finite_threshold_noise <- function(x, tau = finite_tau) {
  x <- as.matrix(x)
  if (tau == 0) {
    return(x)
  }
  x + tau * matrix(rexp(length(x)), nrow(x), ncol(x))
}

simulate_standard_model1 <- function(n) {
  z <- rlogistic_pareto(n, d = 3L, theta = logistic_theta_model1)
  x_signal <- z %*% t(standard_embedding_matrix())
  x <- add_finite_threshold_noise(x_signal)
  list(
    x = x,
    signal = x_signal,
    latent = z,
    label = rep("logistic_low_dim", n)
  )
}

simulate_standard_model2 <- function(n) {
  z_pair <- rlogistic_pareto(n, d = 2L, theta = logistic_theta_model2)
  z_axis <- rpareto(n)
  x_signal <- cbind(z_pair[, 1L], z_pair[, 2L], axis_scale * z_axis)
  x <- add_finite_threshold_noise(x_signal)

  pair_radius <- sqrt(rowSums(z_pair^2))
  axis_radius <- axis_scale * z_axis
  label <- ifelse(axis_radius > pair_radius, "axis_3", "shared_12")

  list(
    x = x,
    signal = x_signal,
    latent_pair = z_pair,
    latent_axis = z_axis,
    label = label
  )
}

axis_contrast <- function(d, axis, mu) {
  e <- rep(0, d)
  e[axis] <- 1
  u <- drop(project_to_tangent(matrix(e, nrow = 1L), mu))
  unit_vector(u, "axis contrast")
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

adaptive_anchor <- function(g, type = adaptive_anchor_type) {
  switch(
    type,
    frechet = spherical_frechet_anchor(g),
    principal = principal_anchor(g),
    stop("Unknown adaptive anchor type.", call. = FALSE)
  )
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

estimate_population_frechet_anchor <- function(simulate_fun,
                                               n = population_anchor_n,
                                               tail_fraction = main_k / main_n,
                                               seed_value = population_anchor_seed) {
  with_seed(seed_value, {
    observations <- simulate_fun(n)
    k <- pmin(n - 1L, pmax(1L, as.integer(round(n * tail_fraction))))
    thr <- threshold_directions(observations$x, k = k)
    list(
      anchor = spherical_frechet_anchor(thr$g),
      n = n,
      k = k,
      tail_fraction = k / n
    )
  })
}

anchor_table <- function(anchor_type, sample_anchor,
                         population_anchor_fit = NULL) {
  rows <- data.frame(
    anchor_type = anchor_type,
    anchor_role = "sample",
    coordinate = seq_along(sample_anchor),
    value = sample_anchor,
    population_anchor_n = NA_integer_,
    population_anchor_k = NA_integer_,
    tail_fraction = NA_real_
  )

  if (!is.null(population_anchor_fit)) {
    rows <- rbind(
      rows,
      data.frame(
        anchor_type = "frechet",
        anchor_role = "population_mc",
        coordinate = seq_along(population_anchor_fit$anchor),
        value = population_anchor_fit$anchor,
        population_anchor_n = population_anchor_fit$n,
        population_anchor_k = population_anchor_fit$k,
        tail_fraction = population_anchor_fit$tail_fraction
      )
    )
  }

  rows
}

target_contrasts_for_anchor <- function(mu, contrasts) {
  contrasts <- as.matrix(contrasts)
  projected <- apply(contrasts, 2L, function(x) {
    unit_vector(project_to_tangent(matrix(x, nrow = 1L), mu),
                "projected population contrast")
  })
  if (is.null(dim(projected))) {
    projected <- matrix(projected, ncol = 1L)
  }
  colnames(projected) <- colnames(contrasts)
  projected
}

fit_at_threshold <- function(x, k, mu, p, true_space = NULL,
                             target_contrast = NULL, labels = NULL) {
  thr <- threshold_directions(x, k = k)
  fit <- agca_fit(thr$g, mu = mu, p = p)
  variation <- agca_variation_explained(fit)

  selected_labels <- if (is.null(labels)) {
    NULL
  } else {
    labels[thr$index]
  }

  subspace_err <- if (is.null(true_space)) {
    NA_real_
  } else {
    subspace_distance(fit$loadings[, seq_len(p), drop = FALSE], true_space)
  }

  loading1_alignment <- if (is.null(target_contrast)) {
    NA_real_
  } else {
    abs(drop(crossprod(fit$loadings[, 1L], target_contrast)))
  }

  axis_fraction <- if (is.null(selected_labels)) {
    NA_real_
  } else {
    mean(selected_labels == "axis_3")
  }

  eigen_cols <- as.data.frame(as.list(fit$eigenvalues))
  names(eigen_cols) <- paste0("eig", seq_along(fit$eigenvalues))

  summary <- cbind(
    data.frame(k = k, threshold = thr$threshold),
    eigen_cols,
    data.frame(
      variation_explained_p = variation[p],
      residual_risk_p = agca_mean_residual(fit, p = p),
      subspace_distance = subspace_err,
      loading1_target_alignment = loading1_alignment,
      axis_fraction = axis_fraction
    )
  )

  list(
    fit = fit,
    threshold = thr,
    selected_labels = selected_labels,
    summary = summary
  )
}

threshold_path <- function(x, k_values, mu, p, true_space = NULL,
                           target_contrast = NULL, labels = NULL) {
  rows <- lapply(
    k_values,
    function(k) {
      fit_at_threshold(
        x = x,
        k = k,
        mu = mu,
        p = p,
        true_space = true_space,
        target_contrast = target_contrast,
        labels = labels
      )$summary
    }
  )
  do.call(rbind, rows)
}

eigen_summary <- function(fit) {
  data.frame(
    component = seq_along(fit$eigenvalues),
    eigenvalue = fit$eigenvalues,
    cumulative_variation = agca_variation_explained(fit)
  )
}

anchor_fits_for_directions <- function(g, p) {
  anchors <- list(
    canonical = canonical_anchor(ncol(g)),
    Frechet = spherical_frechet_anchor(g),
    principal = principal_anchor(g)
  )
  lapply(anchors, function(anchor) {
    agca_fit(g, mu = anchor, p = p)
  })
}

anchor_distance_summary <- function(anchor_fits) {
  canonical_mu <- anchor_fits[["canonical"]]$mu
  data.frame(
    anchor = names(anchor_fits),
    anchor_distance_to_canonical = vapply(
      anchor_fits,
      function(fit) {
        sphere_geodesic_distance(
          matrix(fit$mu, nrow = 1L),
          matrix(canonical_mu, nrow = 1L)
        )
      },
      numeric(1L)
    ),
    row.names = NULL
  )
}

bootstrap_agca_loadings <- function(g, fit, loading_components = 2L,
                                    b = bootstrap_reps,
                                    seed_value = seed + 200000L) {
  with_seed(seed_value, {
    n <- nrow(g)
    d <- ncol(g)
    loading_components <- min(loading_components, ncol(fit$loadings))
    boot_loadings <- array(
      NA_real_,
      dim = c(b, d, loading_components),
      dimnames = list(NULL, paste0("X", seq_len(d)),
                      paste0("component", seq_len(loading_components)))
    )

    for (iter in seq_len(b)) {
      sample_index <- sample.int(n, n, replace = TRUE)
      boot_fit <- agca_fit(g[sample_index, , drop = FALSE],
                           mu = fit$mu, p = fit$p)
      for (component in seq_len(loading_components)) {
        loading <- boot_fit$loadings[, component]
        if (drop(crossprod(loading, fit$loadings[, component])) < 0) {
          loading <- -loading
        }
        boot_loadings[iter, , component] <- loading
      }
    }

    do.call(
      rbind,
      lapply(seq_len(loading_components), function(component) {
        values <- boot_loadings[, , component, drop = TRUE]
        data.frame(
          component = component,
          variable = paste0("X", seq_len(d)),
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
  })
}

oracle_alignment <- function(fit, target_contrasts, target_labels = NULL,
                             n_components = 2L) {
  target_contrasts <- as.matrix(target_contrasts)
  n_components <- min(n_components, ncol(fit$loadings))
  n_targets <- min(ncol(target_contrasts), n_components)
  target_contrasts <- target_contrasts[, seq_len(n_targets), drop = FALSE]
  if (is.null(target_labels)) {
    target_labels <- colnames(target_contrasts)
  }
  if (is.null(target_labels)) {
    target_labels <- paste0("population AGC ", seq_len(n_targets))
  }
  align <- abs(t(fit$loadings[, seq_len(n_components), drop = FALSE]) %*%
                 target_contrasts)
  rownames(align) <- as.character(seq_len(n_components))
  colnames(align) <- target_labels
  out <- as.data.frame(as.table(align))
  names(out) <- c("component", "target", "abs_alignment")
  out$component <- as.integer(as.character(out$component))
  out
}

write_common_results <- function(out_dir, fit, main_summary, path, scores) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- c(
    threshold_path = file.path(out_dir, "threshold_path.csv"),
    main_summary = file.path(out_dir, "main_threshold_summary.csv"),
    eigen_summary = file.path(out_dir, "eigen_summary.csv"),
    rank_summary = file.path(out_dir, "rank_summary.csv"),
    scores = file.path(out_dir, "scores.csv")
  )

  write.csv(path, files[["threshold_path"]], row.names = FALSE)
  write.csv(main_summary, files[["main_summary"]], row.names = FALSE)
  write.csv(eigen_summary(fit), files[["eigen_summary"]], row.names = FALSE)
  write.csv(agca_rank_summary(fit), files[["rank_summary"]], row.names = FALSE)
  write.csv(scores, files[["scores"]], row.names = FALSE)

  files
}

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

save_pdf <- function(file, plot_fun, width = 6.5, height = 5.5) {
  pdf(file, width = width, height = height, bg = "white")
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

plot_text_cex <- list(
  axis = 1.25,
  lab = 1.35,
  main = 1.20,
  legend = 1.30
)

set_plot_text <- function(mar = c(5.2, 5.5, 3.4, 1.2), xpd = FALSE,
                          text_scale = 1) {
  par(
    mar = mar,
    cex.axis = plot_text_cex$axis * text_scale,
    cex.lab = plot_text_cex$lab * text_scale,
    cex.main = plot_text_cex$main * text_scale,
    xpd = xpd
  )
}

draw_score_anchor <- function(xlim, ylim, score_x = NULL, score_y = NULL) {
  label <- expression(mu)
  label_cex <- plot_text_cex$lab
  xspan <- diff(xlim)
  yspan <- diff(ylim)
  candidates <- rbind(
    c(0.055, 0.060),
    c(-0.055, 0.060),
    c(0.055, -0.060),
    c(-0.055, -0.060)
  )
  candidates[, 1L] <- candidates[, 1L] * xspan
  candidates[, 2L] <- candidates[, 2L] * yspan
  choice <- 1L
  if (!is.null(score_x) && !is.null(score_y)) {
    local_counts <- apply(candidates, 1L, function(z) {
      sum(abs(score_x - z[1L]) < 0.10 * xspan &
            abs(score_y - z[2L]) < 0.10 * yspan)
    })
    choice <- which.min(local_counts)
  }
  x <- candidates[choice, 1L]
  y <- candidates[choice, 2L]
  label_width <- strwidth(label, cex = label_cex)
  label_height <- strheight(label, cex = label_cex)
  points(0, 0, pch = 4, lwd = 1.3, cex = 1.0, col = "gray35")
  rect(x - 0.55 * label_width, y - 0.60 * label_height,
       x + 0.55 * label_width, y + 0.60 * label_height,
       col = "white", border = "gray75")
  text(x, y, labels = label, cex = label_cex, col = "black")
}

plot_eigenvalues <- function(fit, file, title) {
  save_pdf(file, function() {
    set_plot_text()
    barplot(
      fit$eigenvalues,
      names.arg = seq_along(fit$eigenvalues),
      xlab = "Component",
      ylab = "Eigenvalue",
      main = title,
      col = "gray70",
      border = "gray35"
    )
  })
}

plot_variation_explained <- function(fits, file, title, text_scale = 1) {
  save_pdf(file, function() {
    set_plot_text(text_scale = text_scale)
    if (inherits(fits, "agca_fit")) {
      fits <- list(canonical = fits)
    }
    max_rank <- max(vapply(fits, function(fit) {
      length(fit$eigenvalues)
    }, integer(1L)))
    rank <- 0L:max_rank

    plot(
      rank,
      rep(NA_real_, length(rank)),
      type = "n",
      ylim = c(0, 1),
      xaxt = "n",
      xlab = "Rank p",
      ylab = "Cumulative proportion explained",
      main = title
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

plot_reconstruction_risk <- function(fit, file, title) {
  save_pdf(file, function() {
    set_plot_text()
    rank_summary <- agca_rank_summary(fit)

    plot(
      rank_summary$rank,
      rank_summary$residual_risk,
      type = "b",
      pch = 16,
      xlab = "Rank p",
      ylab = "Mean residual risk",
      main = title
    )
  })
}

plot_scores <- function(fit, file, labels = NULL) {
  save_pdf(file, function() {
    set_plot_text()

    if (is.null(labels)) {
      point_cols <- adjustcolor("#2C6AA0", alpha.f = 0.35)
    } else {
      label_cols <- c(
        shared_12 = adjustcolor("#1B9E77", alpha.f = 0.35),
        axis_3 = adjustcolor("#D95F02", alpha.f = 0.55)
      )
      point_cols <- label_cols[labels]
    }

    xlim <- extendrange(c(fit$scores[, 1L], 0), f = 0.16)
    ylim <- extendrange(c(fit$scores[, 2L], 0), f = 0.16)
    plot(
      fit$scores[, 1L], fit$scores[, 2L],
      pch = 16,
      cex = 0.55,
      col = point_cols,
      xlim = xlim,
      ylim = ylim,
      xlab = "AGCA score 1",
      ylab = "AGCA score 2",
      main = "Scores of thresholded directions"
    )
    abline(h = 0, v = 0, col = "gray80")
    draw_score_anchor(xlim, ylim, fit$scores[, 1L], fit$scores[, 2L])

    if (!is.null(labels)) {
      legend(
        "topright",
        legend = c("shared_12", "axis_3"),
        col = c("#1B9E77", "#D95F02"),
        pch = 16,
        bty = "n",
        cex = plot_text_cex$legend
      )
    }
  })
}

plot_subspace_path <- function(path, file, title) {
  save_pdf(file, function() {
    set_plot_text()
    ymax <- max(path$subspace_distance, na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) {
      ymax <- 1
    }

    plot(
      path$k,
      path$subspace_distance,
      type = "b",
      pch = 16,
      xlab = "Number of exceedances k",
      ylab = "Projector distance",
      main = title,
      ylim = c(0, ymax * 1.1)
    )
  })
}

plot_loading_contrast <- function(fit, target_contrast, file,
                                  title = "First loading and axis-3 contrast",
                                  target_label = "axis-3 contrast",
                                  target_labels = NULL,
                                  legend_position = "topright") {
  save_pdf(file, function() {
    set_plot_text()

    target_contrast <- as.matrix(target_contrast)
    n_components <- min(ncol(target_contrast), ncol(fit$loadings))
    target_contrast <- target_contrast[, seq_len(n_components), drop = FALSE]
    loading <- fit$loadings[, seq_len(n_components), drop = FALSE]

    for (j in seq_len(n_components)) {
      if (drop(crossprod(loading[, j], target_contrast[, j])) < 0) {
        loading[, j] <- -loading[, j]
      }
    }

    if (is.null(target_labels)) {
      target_labels <- if (n_components == 1L) {
        target_label
      } else {
        paste0("population AGC ", seq_len(n_components))
      }
    }

    ylim <- range(c(loading, target_contrast))
    y_span <- diff(ylim)
    if (!is.finite(y_span) || y_span <= 0) {
      y_span <- max(1, max(abs(ylim)))
    }
    top_pad <- if (identical(legend_position, "top")) 0.38 else 0.15
    ylim <- c(ylim[1L] - 0.12 * y_span, ylim[2L] + top_pad * y_span)

    if (n_components == 1L) {
      centers <- barplot(
        loading[, 1L],
        names.arg = paste0("X", seq_len(nrow(loading))),
        ylim = ylim,
        ylab = "Loading",
        main = title,
        col = "gray75",
        border = "gray35"
      )
      points(centers, target_contrast[, 1L], pch = 16, col = "#D95F02")
      legend(
        legend_position,
        legend = c("estimated first loading", target_labels[1L]),
        fill = c("gray75", NA),
        border = c("gray35", NA),
        pch = c(NA, 16),
        col = c("gray35", "#D95F02"),
        bty = "n",
        horiz = identical(legend_position, "top"),
        cex = plot_text_cex$legend
      )
    } else {
      loading_cols <- c("gray75", "gray45", "#9ECAE1", "#6BAED6")
      target_cols <- c("#D95F02", "#CC6677", "#117733", "#332288")
      loading_cols <- loading_cols[seq_len(n_components)]
      target_cols <- target_cols[seq_len(n_components)]

      centers <- barplot(
        t(loading),
        beside = TRUE,
        names.arg = paste0("X", seq_len(nrow(loading))),
        ylim = ylim,
        ylab = "Loading",
        main = title,
        col = loading_cols,
        border = "gray35"
      )
      for (j in seq_len(n_components)) {
        points(centers[j, ], target_contrast[, j], pch = 16, col = target_cols[j])
      }

      legend(
        legend_position,
        legend = c(
          paste0("est. AGC ", seq_len(n_components)),
          sub("^population", "pop.", target_labels)
        ),
        fill = c(loading_cols, rep(NA, n_components)),
        border = c(rep("gray35", n_components), rep(NA, n_components)),
        pch = c(rep(NA, n_components), rep(16, n_components)),
        col = c(rep("gray35", n_components), target_cols),
        bty = "n",
        horiz = identical(legend_position, "top"),
        x.intersp = 0.55,
        cex = if (identical(legend_position, "top")) {
          0.82 * plot_text_cex$legend
        } else {
          plot_text_cex$legend
        }
      )
    }
  })
}

plot_loading <- function(fit, file, component, loading_intervals = NULL,
                         population_loading = NULL,
                         title = paste0("AGC ", component),
                         ylab = "Loading",
                         bar_width = NULL,
                         plot_width = 5.6,
                         plot_height = 4.9,
                         text_scale = 1,
                         title_scale = text_scale) {
  save_pdf(file, function() {
    left_margin <- if (nzchar(ylab)) 5.2 else 3.5
    set_plot_text(
      mar = c(5.5, left_margin, 3.4, 1.2),
      text_scale = text_scale
    )
    par(cex.main = plot_text_cex$main * title_scale)
    variable_labels <- paste0("X", seq_len(nrow(fit$loadings)))
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

    if (is.null(bar_width)) {
      centers <- barplot(
        loading,
        names.arg = variable_labels,
        ylim = ylim,
        col = "gray75",
        border = "gray35",
        ylab = ylab,
        main = title
      )
    } else {
      bar_gap <- (1 - bar_width) / 4
      center_step <- bar_width + bar_gap
      centers <- seq_along(loading) * center_step
      x_pad <- bar_width / 2 + bar_gap / 2
      plot(
        centers,
        loading,
        type = "n",
        xaxt = "n",
        xlab = "",
        ylab = ylab,
        xlim = range(centers) + c(-x_pad, x_pad),
        ylim = ylim,
        main = title
      )
      axis(1, at = centers, labels = variable_labels)
      rect(
        centers - bar_width / 2,
        pmin(0, loading),
        centers + bar_width / 2,
        pmax(0, loading),
        col = "gray75",
        border = "gray35"
      )
    }
    abline(h = 0, col = "gray70")

    legend_items <- character(0)
    legend_lty <- numeric(0)
    legend_lwd <- numeric(0)
    legend_pch <- numeric(0)
    legend_col <- character(0)

    if (!is.null(current_intervals) &&
        nrow(current_intervals) == length(variable_labels)) {
      segments(
        centers,
        current_intervals$q025,
        centers,
        current_intervals$q975,
        lwd = 1.3,
        col = "black"
      )
      segments(
        centers - 0.055,
        current_intervals$q025,
        centers + 0.055,
        current_intervals$q025,
        lwd = 1.3,
        col = "black"
      )
      segments(
        centers - 0.055,
        current_intervals$q975,
        centers + 0.055,
        current_intervals$q975,
        lwd = 1.3,
        col = "black"
      )
      points(centers, current_intervals$main_loading, pch = 16, cex = 0.55)
      legend_items <- c(legend_items, "bootstrap 95% DI")
      legend_lty <- c(legend_lty, 1)
      legend_lwd <- c(legend_lwd, 1.3)
      legend_pch <- c(legend_pch, NA)
      legend_col <- c(legend_col, "black")
    }

    if (!is.null(population_loading)) {
      points(centers, population_loading, pch = 16, cex = 1.0, col = "#D95F02")
      legend_items <- c(legend_items, "population loading")
      legend_lty <- c(legend_lty, NA)
      legend_lwd <- c(legend_lwd, NA)
      legend_pch <- c(legend_pch, 16)
      legend_col <- c(legend_col, "#D95F02")
    }

    if (length(legend_items) > 0L) {
      scaled_legend <- text_scale > 1.05
      legend_labels <- if (scaled_legend) {
        sub(" loading$", "", legend_items)
      } else {
        legend_items
      }
      legend(
        "top",
        legend = legend_labels,
        lty = legend_lty,
        lwd = legend_lwd,
        pch = legend_pch,
        col = legend_col,
        bty = "n",
        horiz = !scaled_legend,
        x.intersp = 0.60,
        y.intersp = if (scaled_legend) 0.85 else 1.00,
        cex = if (scaled_legend) {
          0.70 * plot_text_cex$legend * text_scale
        } else {
          0.78 * plot_text_cex$legend
        }
      )
    }
  }, width = plot_width, height = plot_height)
}

plot_anchor_sensitivity_distance <- function(anchor_summary, file,
                                             title = "Anchor distance") {
  save_pdf(file, function() {
    set_plot_text(mar = c(4.9, 4.8, 3.4, 1.0), xpd = NA)
    values <- anchor_summary$anchor_distance_to_canonical
    centers <- barplot(
      values,
      names.arg = FALSE,
      las = 2,
      ylim = extendrange(c(0, values), f = 0.10),
      width = 0.65,
      space = 0,
      col = c("gray70", "#2166AC", "#B2182B")[seq_len(nrow(anchor_summary))],
      border = "gray35",
      ylab = "Geodesic distance",
      main = title
    )
    axis(1, at = centers, labels = FALSE)
    text(
      centers,
      par("usr")[3L] - 0.035 * diff(par("usr")[3:4]),
      labels = anchor_summary$anchor,
      srt = 35,
      adj = 1,
      cex = 0.95 * plot_text_cex$axis
    )
  }, width = 3.6, height = 4.8)
}

plot_oracle_alignment <- function(alignment, file,
                                  title = "Oracle loading alignment") {
  save_pdf(file, function() {
    set_plot_text(mar = c(5.2, 8.2, 3.2, 1.0))
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
      main = title
    )
    axis(1, at = seq_len(nrow(z)), labels = rownames(z))
    axis(2, at = seq_len(ncol(z)), labels = colnames(z), las = 1,
         cex.axis = 0.88)
    mtext("Population loading", side = 2, line = 6.4, cex = 1.05)
    for (i in seq_len(nrow(z))) {
      for (j in seq_len(ncol(z))) {
        text(i, j, labels = sprintf("%.2f", z[i, j]), cex = 0.82)
      }
    }
  }, width = 7.0, height = 5.4)
}

geodesic_segment <- function(mu, loading, n = 1000L) {
  loading <- unit_vector(project_to_tangent(matrix(loading, nrow = 1L), mu),
                         "loading")
  t_grid <- seq(-pi / 2 + 1e-4, pi / 2 - 1e-4, length.out = n)
  segment <- cos(t_grid) * matrix(mu, n, length(mu), byrow = TRUE) +
    sin(t_grid) * matrix(loading, n, length(mu), byrow = TRUE)
  segment[rowSums(segment >= -1e-10) == length(mu), , drop = FALSE]
}

plot_sphere_geodesic <- function(fit, file, title, labels = NULL,
                                 target_contrasts = NULL,
                                 target_labels = NULL,
                                 show_estimated_second = FALSE,
                                 anchor_label = "anchor",
                                 reference_anchors = NULL,
                                 reference_anchor_labels = NULL,
                                 title_scale = 1,
                                 max_points = 800L) {
  save_pdf(file, function() {
    set_plot_text(mar = c(4.0, 4.0, 3.4, 1.2), xpd = NA)
    par(cex.main = plot_text_cex$main * title_scale)
    face_on_phi <- 23

    grid <- seq(0, 1, length.out = 60L)
    z <- outer(grid, grid, function(x, y) {
      inside <- x^2 + y^2 <= 1
      out <- sqrt(pmax(0, 1 - x^2 - y^2))
      out[!inside] <- NA_real_
      out
    })

    pmat <- persp(
      grid, grid, z,
      theta = 75,
      phi = face_on_phi,
      r = 4,
      d = 1.5,
      expand = 0.90,
      col = adjustcolor("gray95", alpha.f = 0.55),
      border = "gray85",
      xlab = "X1",
      ylab = "X2",
      zlab = "X3",
      ticktype = "detailed",
      nticks = 1,
      main = title
    )

    g <- fit$g
    if (nrow(g) > max_points) {
      g <- g[seq_len(max_points), , drop = FALSE]
      labels <- if (is.null(labels)) NULL else labels[seq_len(max_points)]
    }

    if (is.null(labels)) {
      point_cols <- adjustcolor("#2C6AA0", alpha.f = 0.45)
    } else {
      label_cols <- c(
        shared_12 = adjustcolor("#1B9E77", alpha.f = 0.45),
        axis_3 = adjustcolor("#D95F02", alpha.f = 0.60)
      )
      point_cols <- label_cols[labels]
    }

    projected_points <- trans3d(g[, 1L], g[, 2L], g[, 3L], pmat)
    points(projected_points, pch = 16, cex = 0.68, col = point_cols)

    estimated1 <- geodesic_segment(fit$mu, fit$loadings[, 1L])
    estimated1_points <- trans3d(estimated1[, 1L], estimated1[, 2L],
                                 estimated1[, 3L], pmat)
    lines(estimated1_points, col = "#542788", lwd = 3)

    if (show_estimated_second && ncol(fit$loadings) >= 2L) {
      estimated2 <- geodesic_segment(fit$mu, fit$loadings[, 2L])
      estimated2_points <- trans3d(estimated2[, 1L], estimated2[, 2L],
                                   estimated2[, 3L], pmat)
      lines(estimated2_points, col = "#1B9E77", lwd = 2.5)
    }

    if (!is.null(target_contrasts)) {
      target_contrasts <- as.matrix(target_contrasts)
      if (is.null(target_labels)) {
        target_labels <- paste0("population AGC ", seq_len(ncol(target_contrasts)))
      }
      target_cols <- c("#D95F02", "#CC6677", "#117733")
      target_lty <- c(2, 4, 3)
      for (j in seq_len(ncol(target_contrasts))) {
        target <- geodesic_segment(fit$mu, target_contrasts[, j])
        target_points <- trans3d(target[, 1L], target[, 2L], target[, 3L], pmat)
        lines(
          target_points,
          col = target_cols[(j - 1L) %% length(target_cols) + 1L],
          lwd = 2.5,
          lty = target_lty[(j - 1L) %% length(target_lty) + 1L]
        )
      }
    }

    anchor <- trans3d(fit$mu[1L], fit$mu[2L], fit$mu[3L], pmat)
    points(anchor, pch = 17, cex = 1.3, col = "black")

    if (!is.null(reference_anchors)) {
      reference_anchors <- as.matrix(reference_anchors)
      if (ncol(reference_anchors) != length(fit$mu) &&
          nrow(reference_anchors) == length(fit$mu)) {
        reference_anchors <- t(reference_anchors)
      }
      if (ncol(reference_anchors) != length(fit$mu)) {
        stop("reference_anchors has incompatible dimensions.", call. = FALSE)
      }
      if (is.null(reference_anchor_labels)) {
        reference_anchor_labels <- paste0("reference anchor ", seq_len(nrow(reference_anchors)))
      }
      reference_anchor_cols <- c("#E7298A", "#7570B3", "#66A61E")
      reference_anchor_pch <- c(4, 8, 3)
      for (j in seq_len(nrow(reference_anchors))) {
        reference_anchor <- trans3d(
          reference_anchors[j, 1L],
          reference_anchors[j, 2L],
          reference_anchors[j, 3L],
          pmat
        )
        points(
          reference_anchor,
          pch = reference_anchor_pch[(j - 1L) %% length(reference_anchor_pch) + 1L],
          cex = 1.35,
          lwd = 2,
          col = reference_anchor_cols[(j - 1L) %% length(reference_anchor_cols) + 1L]
        )
      }
    }

    legend_items <- c("selected directions", "estimated AGC 1", anchor_label)
    legend_cols <- c("#2C6AA0", "#542788", "black")
    legend_lty <- c(NA, 1, NA)
    legend_pch <- c(16, NA, 17)

    if (!is.null(labels)) {
      legend_items <- c("shared_12", "axis_3", "estimated AGC 1", anchor_label)
      legend_cols <- c("#1B9E77", "#D95F02", "#542788", "black")
      legend_lty <- c(NA, NA, 1, NA)
      legend_pch <- c(16, 16, NA, 17)
    }

    if (show_estimated_second && ncol(fit$loadings) >= 2L) {
      legend_items <- c(legend_items, "estimated AGC 2")
      legend_cols <- c(legend_cols, "#1B9E77")
      legend_lty <- c(legend_lty, 1)
      legend_pch <- c(legend_pch, NA)
    }

    if (!is.null(reference_anchors)) {
      reference_anchor_cols <- c("#E7298A", "#7570B3", "#66A61E")
      reference_anchor_pch <- c(4, 8, 3)
      for (j in seq_len(nrow(reference_anchors))) {
        legend_items <- c(legend_items, reference_anchor_labels[j])
        legend_cols <- c(
          legend_cols,
          reference_anchor_cols[(j - 1L) %% length(reference_anchor_cols) + 1L]
        )
        legend_lty <- c(legend_lty, NA)
        legend_pch <- c(
          legend_pch,
          reference_anchor_pch[(j - 1L) %% length(reference_anchor_pch) + 1L]
        )
      }
    }

    if (!is.null(target_contrasts)) {
      target_contrasts <- as.matrix(target_contrasts)
      if (is.null(target_labels)) {
        target_labels <- paste0("population AGC ", seq_len(ncol(target_contrasts)))
      }
      target_cols <- c("#D95F02", "#CC6677", "#117733")
      target_lty <- c(2, 4, 3)
      for (j in seq_len(ncol(target_contrasts))) {
        legend_items <- c(legend_items, target_labels[j])
        legend_cols <- c(legend_cols, target_cols[(j - 1L) %% length(target_cols) + 1L])
        legend_lty <- c(legend_lty, target_lty[(j - 1L) %% length(target_lty) + 1L])
        legend_pch <- c(legend_pch, NA)
      }
    }

    legend(
      "topright",
      legend = legend_items,
      col = legend_cols,
      lty = legend_lty,
      pch = legend_pch,
      lwd = 2,
      bty = "n",
      cex = if (length(legend_items) > 6L) {
        0.78 * plot_text_cex$legend
      } else {
        plot_text_cex$legend
      }
    )
  }, width = 7.0, height = 6.2)
}

plot_sphere_agc1_anchor_overlay <- function(canonical_fit, frechet_fit, file,
                                            title, labels = NULL,
                                            max_points = 800L) {
  save_pdf(file, function() {
    set_plot_text(mar = c(4.0, 4.0, 3.4, 1.2), xpd = NA)
    face_on_phi <- 23
    frechet_red <- "#B2182B"

    grid <- seq(0, 1, length.out = 60L)
    z <- outer(grid, grid, function(x, y) {
      inside <- x^2 + y^2 <= 1
      out <- sqrt(pmax(0, 1 - x^2 - y^2))
      out[!inside] <- NA_real_
      out
    })

    pmat <- persp(
      grid, grid, z,
      theta = 75,
      phi = face_on_phi,
      r = 4,
      d = 1.5,
      expand = 0.90,
      col = adjustcolor("gray95", alpha.f = 0.55),
      border = "gray85",
      xlab = "X1",
      ylab = "X2",
      zlab = "X3",
      ticktype = "detailed",
      nticks = 1,
      main = title
    )

    g <- canonical_fit$g
    if (nrow(g) > max_points) {
      g <- g[seq_len(max_points), , drop = FALSE]
      labels <- if (is.null(labels)) NULL else labels[seq_len(max_points)]
    }

    if (is.null(labels)) {
      point_cols <- adjustcolor("#2C6AA0", alpha.f = 0.45)
      legend_items <- "selected directions"
      legend_cols <- "#2C6AA0"
      legend_lty <- NA
      legend_pch <- 16
    } else {
      label_cols <- c(
        shared_12 = adjustcolor("#1B9E77", alpha.f = 0.45),
        axis_3 = adjustcolor("#D95F02", alpha.f = 0.60)
      )
      point_cols <- label_cols[labels]
      legend_items <- c("shared X1-X2", "axis X3")
      legend_cols <- c("#1B9E77", "#D95F02")
      legend_lty <- c(NA, NA)
      legend_pch <- c(16, 16)
    }

    projected_points <- trans3d(g[, 1L], g[, 2L], g[, 3L], pmat)
    points(projected_points, pch = 16, cex = 0.68, col = point_cols)

    canonical_agc1 <- geodesic_segment(canonical_fit$mu,
                                       canonical_fit$loadings[, 1L])
    canonical_points <- trans3d(canonical_agc1[, 1L], canonical_agc1[, 2L],
                                canonical_agc1[, 3L], pmat)
    lines(canonical_points, col = "black", lwd = 3, lty = 1)

    frechet_agc1 <- geodesic_segment(frechet_fit$mu,
                                     frechet_fit$loadings[, 1L])
    frechet_points <- trans3d(frechet_agc1[, 1L], frechet_agc1[, 2L],
                              frechet_agc1[, 3L], pmat)
    lines(frechet_points, col = frechet_red, lwd = 3, lty = 5)

    canonical_anchor <- trans3d(canonical_fit$mu[1L], canonical_fit$mu[2L],
                                canonical_fit$mu[3L], pmat)
    points(canonical_anchor, pch = 17, cex = 1.45, col = "black")

    frechet_anchor <- trans3d(frechet_fit$mu[1L], frechet_fit$mu[2L],
                              frechet_fit$mu[3L], pmat)
    points(frechet_anchor, pch = 15, cex = 1.35, col = frechet_red)

    legend_items <- c(
      legend_items,
      "canonical AGC 1",
      "canonical anchor",
      "Frechet AGC 1",
      "Frechet anchor"
    )
    legend_cols <- c(legend_cols, "black", "black", frechet_red, frechet_red)
    legend_lty <- c(legend_lty, 1, NA, 5, NA)
    legend_pch <- c(legend_pch, NA, 17, NA, 15)

    legend(
      "topright",
      legend = legend_items,
      col = legend_cols,
      lty = legend_lty,
      pch = legend_pch,
      lwd = 2.2,
      bty = "n",
      cex = 0.82 * plot_text_cex$legend
    )
  }, width = 7.0, height = 5.8)
}

prepare_display_diagnostics <- function(out_dir, main, p, target_contrasts,
                                        target_labels,
                                        bootstrap_seed) {
  anchor_fits <- anchor_fits_for_directions(main$fit$g, p = p)
  anchor_summary <- anchor_distance_summary(anchor_fits)
  alignment <- oracle_alignment(
    main$fit,
    target_contrasts,
    target_labels = target_labels,
    n_components = 2L
  )

  files <- c(
    anchor_sensitivity = file.path(out_dir, "anchor_sensitivity.csv"),
    oracle_alignment = file.path(out_dir, "oracle_alignment.csv")
  )
  write.csv(anchor_summary, files[["anchor_sensitivity"]], row.names = FALSE)
  write.csv(alignment, files[["oracle_alignment"]], row.names = FALSE)

  loading_intervals <- NULL
  if (!skip_bootstrap) {
    loading_intervals <- bootstrap_agca_loadings(
      main$fit$g,
      main$fit,
      loading_components = 2L,
      seed_value = bootstrap_seed
    )
    files <- c(
      files,
      bootstrap_loading_intervals = file.path(out_dir, "bootstrap_loading_intervals.csv")
    )
    write.csv(
      loading_intervals,
      files[["bootstrap_loading_intervals"]],
      row.names = FALSE
    )
  }

  list(
    anchor_fits = anchor_fits,
    anchor_summary = anchor_summary,
    oracle_alignment = alignment,
    loading_intervals = loading_intervals,
    target_contrasts = as.matrix(target_contrasts),
    target_labels = target_labels,
    result_files = files
  )
}

plot_display_diagnostics <- function(main, diagnostics, plot_files,
                                     variation_title, anchor_title,
                                     loading_title_prefix, oracle_title,
                                     loading_ylabs = NULL,
                                     loading_bar_width = NULL,
                                     loading_plot_height = 4.9,
                                     loading_text_scale = 1,
                                     loading_title_scale = loading_text_scale,
                                     variation_text_scale = 1) {
  plot_variation_explained(
    diagnostics$anchor_fits,
    plot_files[["variation_explained"]],
    variation_title,
    text_scale = variation_text_scale
  )
  for (component in seq_len(min(2L, ncol(main$fit$loadings)))) {
    ylab <- if (is.null(loading_ylabs)) {
      "Loading"
    } else {
      loading_ylabs[[component]]
    }
    plot_loading(
      main$fit,
      plot_files[[paste0("loading_agc", component)]],
      component = component,
      loading_intervals = diagnostics$loading_intervals,
      population_loading = diagnostics$target_contrasts[, component],
      title = paste0(loading_title_prefix, ": AGC ", component),
      ylab = ylab,
      bar_width = loading_bar_width,
      plot_height = loading_plot_height,
      text_scale = loading_text_scale,
      title_scale = loading_title_scale
    )
  }
  plot_anchor_sensitivity_distance(
    diagnostics$anchor_summary,
    plot_files[["anchor_sensitivity_distance"]],
    anchor_title
  )
  plot_oracle_alignment(
    diagnostics$oracle_alignment,
    plot_files[["oracle_alignment"]],
    oracle_title
  )
}

prepare_frechet_outputs <- function(out_dir, frechet_fit, target_contrasts,
                                    bootstrap_seed) {
  files <- character(0)
  loading_intervals <- NULL
  if (!skip_bootstrap) {
    loading_intervals <- bootstrap_agca_loadings(
      frechet_fit$g,
      frechet_fit,
      loading_components = 2L,
      seed_value = bootstrap_seed
    )
    files <- c(
      files,
      bootstrap_loading_intervals_frechet = file.path(
        out_dir,
        "bootstrap_loading_intervals_frechet.csv"
      )
    )
    write.csv(
      loading_intervals,
      files[["bootstrap_loading_intervals_frechet"]],
      row.names = FALSE
    )
  }

  list(
    fit = frechet_fit,
    target_contrasts = as.matrix(target_contrasts),
    loading_intervals = loading_intervals,
    result_files = files
  )
}

plot_frechet_outputs <- function(frechet_outputs, plot_files,
                                 loading_title_prefix, sphere_title,
                                 target_labels, labels = NULL,
                                 population_anchor_fit = NULL,
                                 loading_ylabs = c("Loading", ""),
                                 loading_bar_width = 0.25,
                                 loading_plot_width = 6.6,
                                 loading_plot_height = 6.5,
                                 loading_text_scale = 2.0,
                                 loading_title_scale = 2.35) {
  for (component in seq_len(min(2L, ncol(frechet_outputs$fit$loadings)))) {
    ylab <- if (length(loading_ylabs) >= component) {
      loading_ylabs[[component]]
    } else {
      "Loading"
    }
    plot_loading(
      frechet_outputs$fit,
      plot_files[[paste0("loading_agc", component, "_frechet")]],
      component = component,
      loading_intervals = frechet_outputs$loading_intervals,
      population_loading = frechet_outputs$target_contrasts[, component],
      title = paste0(loading_title_prefix, ": AGC ", component),
      ylab = ylab,
      bar_width = loading_bar_width,
      plot_width = loading_plot_width,
      plot_height = loading_plot_height,
      text_scale = loading_text_scale,
      title_scale = loading_title_scale
    )
  }

  plot_sphere_geodesic(
    frechet_outputs$fit,
    plot_files[["sphere_geodesic_frechet"]],
    sphere_title,
    labels = labels,
    target_contrasts = frechet_outputs$target_contrasts,
    target_labels = target_labels,
    show_estimated_second = TRUE,
    anchor_label = "sample Frechet anchor",
    reference_anchors = if (is.null(population_anchor_fit)) {
      NULL
    } else {
      matrix(population_anchor_fit$anchor, nrow = 1L)
    },
    reference_anchor_labels = if (is.null(population_anchor_fit)) {
      NULL
    } else {
      "population Frechet anchor"
    }
  )
}

run_model1 <- function(mu, contrasts, population_anchor_fit = NULL) {
  n <- main_n
  k <- main_k
  k_values <- main_k_values
  p <- 1L
  true_space <- contrasts[, 1L, drop = FALSE]

  observations <- simulate_standard_model1(n)
  main <- fit_at_threshold(observations$x, k = k, mu = mu, p = p,
                           true_space = true_space,
                           target_contrast = contrasts[, 1L])
  path <- threshold_path(
    observations$x,
    k_values = k_values,
    mu = mu,
    p = p,
    true_space = true_space,
    target_contrast = contrasts[, 1L]
  )

  print_fit_summary("Standard Model 1: logistic block with weak residual direction",
                    main$fit, p = p)
  cat("Projector distance to target geodesic span(b1): ",
      round(main$summary$subspace_distance, 4), "\n", sep = "")

  out_dir <- file.path(output_dir, "model1")
  scores <- as.data.frame(main$fit$scores)
  names(scores) <- paste0("score", seq_len(ncol(scores)))

  result_files <- write_common_results(
    out_dir = out_dir,
    fit = main$fit,
    main_summary = main$summary,
    path = path,
    scores = scores
  )
  target_labels <- c("population AGC 1", "population AGC 2")
  diagnostics <- prepare_display_diagnostics(
    out_dir = out_dir,
    main = main,
    p = p,
    target_contrasts = contrasts,
    target_labels = target_labels,
    bootstrap_seed = seed + 210001L
  )
  result_files <- c(result_files, diagnostics$result_files)
  frechet_target_contrasts <- target_contrasts_for_anchor(
    diagnostics$anchor_fits[["Frechet"]]$mu,
    contrasts
  )
  frechet_outputs <- prepare_frechet_outputs(
    out_dir = out_dir,
    frechet_fit = diagnostics$anchor_fits[["Frechet"]],
    target_contrasts = frechet_target_contrasts,
    bootstrap_seed = seed + 215001L
  )
  result_files <- c(result_files, frechet_outputs$result_files)

  plot_files <- c(
    eigenvalues = file.path(out_dir, "eigenvalues.pdf"),
    variation_explained = file.path(out_dir, "anchored_variation_explained.pdf"),
    reconstruction_risk = file.path(out_dir, "reconstruction_risk.pdf"),
    scores = file.path(out_dir, "scores.pdf"),
    loading_agc1 = file.path(out_dir, "loading_agc1.pdf"),
    loading_agc2 = file.path(out_dir, "loading_agc2.pdf"),
    loading_agc1_frechet = file.path(out_dir, "loading_agc1_frechet.pdf"),
    loading_agc2_frechet = file.path(out_dir, "loading_agc2_frechet.pdf"),
    anchor_sensitivity_distance = file.path(out_dir, "anchor_sensitivity_distance.pdf"),
    oracle_alignment = file.path(out_dir, "oracle_alignment.pdf"),
    subspace_distance = file.path(out_dir, "subspace_distance.pdf"),
    sphere_agc1_anchor_overlay = file.path(out_dir, "sphere_agc1_anchor_overlay.pdf"),
    sphere_geodesic = file.path(out_dir, "sphere_geodesic.pdf"),
    sphere_geodesic_frechet = file.path(out_dir, "sphere_geodesic_frechet.pdf")
  )

  plot_eigenvalues(main$fit, plot_files[["eigenvalues"]],
                   "Standard Model 1: AGCA eigenvalues")
  plot_display_diagnostics(
    main,
    diagnostics,
    plot_files,
    variation_title = "Model 1: anchored variation explained",
    anchor_title = "Model 1: anchor distance",
    loading_title_prefix = "Model 1",
    oracle_title = "Model 1: oracle loading alignment",
    loading_ylabs = c("Loading", ""),
    loading_bar_width = 0.25,
    loading_plot_height = 6.5,
    loading_text_scale = 2.0,
    loading_title_scale = 2.35,
    variation_text_scale = 1.25
  )
  plot_reconstruction_risk(main$fit, plot_files[["reconstruction_risk"]],
                           "Standard Model 1: reconstruction risk")
  plot_scores(main$fit, plot_files[["scores"]])
  plot_subspace_path(path, plot_files[["subspace_distance"]],
                     "Estimated geodesic vs target")
  plot_sphere_geodesic(
    main$fit,
    plot_files[["sphere_geodesic"]],
    "Model 1: canonical-anchor AGCs",
    target_contrasts = contrasts,
    target_labels = c("population AGC 1", "population AGC 2"),
    title_scale = 1.45
  )
  plot_sphere_agc1_anchor_overlay(
    main$fit,
    diagnostics$anchor_fits[["Frechet"]],
    plot_files[["sphere_agc1_anchor_overlay"]],
    "Model 1: low-dimensional angular law"
  )
  plot_frechet_outputs(
    frechet_outputs,
    plot_files,
    loading_title_prefix = "Frechet M1",
    sphere_title = "Model 1: Frechet-anchor AGCs",
    target_labels = target_labels,
    population_anchor_fit = population_anchor_fit,
    loading_ylabs = c("Loading", ""),
    loading_bar_width = 0.25,
    loading_plot_width = 6.6,
    loading_plot_height = 6.5,
    loading_text_scale = 2.0,
    loading_title_scale = 2.35
  )

  list(
    main = main,
    path = path,
    observations = observations,
    out_dir = out_dir,
    result_files = result_files,
    plot_files = plot_files
  )
}

run_model1_adaptive_anchor <- function(observations, contrasts,
                                       anchor_type = adaptive_anchor_type,
                                       population_anchor_fit = NULL) {
  k <- main_k
  k_values <- main_k_values
  p <- 1L

  main_threshold <- threshold_directions(observations$x, k = k)
  mu <- adaptive_anchor(main_threshold$g, type = anchor_type)
  target_contrasts <- target_contrasts_for_anchor(mu, contrasts)
  true_space <- target_contrasts[, 1L, drop = FALSE]

  main <- fit_at_threshold(
    observations$x,
    k = k,
    mu = mu,
    p = p,
    true_space = true_space,
    target_contrast = target_contrasts[, 1L]
  )
  path <- threshold_path(
    observations$x,
    k_values = k_values,
    mu = mu,
    p = p,
    true_space = true_space,
    target_contrast = target_contrasts[, 1L]
  )

  title_prefix <- paste0("Standard Model 1: ", anchor_type, " anchor")
  print_fit_summary(title_prefix, main$fit, p = p)
  cat("Adaptive anchor: ",
      paste(round(mu, 4), collapse = ", "), "\n", sep = "")
  cat("Projector distance to projected target geodesic: ",
      round(main$summary$subspace_distance, 4), "\n", sep = "")

  out_dir <- file.path(output_dir, paste0("model1_", anchor_type, "_anchor"))
  scores <- as.data.frame(main$fit$scores)
  names(scores) <- paste0("score", seq_len(ncol(scores)))

  result_files <- write_common_results(
    out_dir = out_dir,
    fit = main$fit,
    main_summary = main$summary,
    path = path,
    scores = scores
  )
  result_files <- c(
    result_files,
    anchor = file.path(out_dir, "anchor.csv")
  )
  write.csv(
    anchor_table(anchor_type, mu, population_anchor_fit = population_anchor_fit),
    result_files[["anchor"]],
    row.names = FALSE
  )
  target_labels <- c("population AGC 1", "population AGC 2")
  diagnostics <- prepare_display_diagnostics(
    out_dir = out_dir,
    main = main,
    p = p,
    target_contrasts = target_contrasts,
    target_labels = target_labels,
    bootstrap_seed = seed + 220001L
  )
  result_files <- c(result_files, diagnostics$result_files)

  plot_files <- c(
    eigenvalues = file.path(out_dir, "eigenvalues.pdf"),
    variation_explained = file.path(out_dir, "anchored_variation_explained.pdf"),
    reconstruction_risk = file.path(out_dir, "reconstruction_risk.pdf"),
    scores = file.path(out_dir, "scores.pdf"),
    loading_agc1 = file.path(out_dir, "loading_agc1.pdf"),
    loading_agc2 = file.path(out_dir, "loading_agc2.pdf"),
    anchor_sensitivity_distance = file.path(out_dir, "anchor_sensitivity_distance.pdf"),
    oracle_alignment = file.path(out_dir, "oracle_alignment.pdf"),
    subspace_distance = file.path(out_dir, "subspace_distance.pdf"),
    sphere_geodesic = file.path(out_dir, "sphere_geodesic.pdf")
  )

  plot_eigenvalues(main$fit, plot_files[["eigenvalues"]],
                   paste0(title_prefix, ": AGCA eigenvalues"))
  plot_display_diagnostics(
    main,
    diagnostics,
    plot_files,
    variation_title = paste0("Model 1: ", anchor_type, "-anchor variation"),
    anchor_title = "Model 1: anchor distance",
    loading_title_prefix = paste0("Model 1: ", anchor_type, " anchor"),
    oracle_title = paste0("Model 1: ", anchor_type, "-anchor oracle alignment")
  )
  plot_reconstruction_risk(main$fit, plot_files[["reconstruction_risk"]],
                           paste0(title_prefix, ": reconstruction risk"))
  plot_scores(main$fit, plot_files[["scores"]])
  plot_subspace_path(path, plot_files[["subspace_distance"]],
                     "Estimated geodesic vs projected target")
  plot_sphere_geodesic(
    main$fit,
    plot_files[["sphere_geodesic"]],
    paste0("Model 1: ", anchor_type, "-anchor AGCs"),
    target_contrasts = target_contrasts,
    target_labels = c("population AGC 1", "population AGC 2"),
    show_estimated_second = TRUE,
    anchor_label = if (identical(anchor_type, "frechet")) {
      "sample Frechet anchor"
    } else {
      paste0(anchor_type, " anchor")
    },
    reference_anchors = if (is.null(population_anchor_fit)) {
      NULL
    } else {
      matrix(population_anchor_fit$anchor, nrow = 1L)
    },
    reference_anchor_labels = if (is.null(population_anchor_fit)) {
      NULL
    } else {
      "population Frechet anchor"
    }
  )

  list(
    main = main,
    path = path,
    anchor = mu,
    out_dir = out_dir,
    result_files = result_files,
    plot_files = plot_files
  )
}

run_model2 <- function(mu, population_anchor_fit = NULL) {
  n <- main_n
  k <- main_k
  k_values <- main_k_values
  p <- 1L
  target <- axis_contrast(3L, axis = 3L, mu = mu)
  target_contrasts <- target_contrasts_for_anchor(
    mu,
    cbind(c3 = target, b1 = contrast_basis_d3()[, "b1"])
  )

  observations <- simulate_standard_model2(n)
  main <- fit_at_threshold(
    observations$x,
    k = k,
    mu = mu,
    p = p,
    target_contrast = target,
    labels = observations$label
  )
  path <- threshold_path(
    observations$x,
    k_values = k_values,
    mu = mu,
    p = p,
    target_contrast = target,
    labels = observations$label
  )

  print_fit_summary(
    "Standard Model 2: logistic pair plus asymptotically independent X3",
    main$fit,
    p = p
  )
  cat("Selected exceedance types:\n")
  print(table(main$selected_labels))
  cat("Alignment of first loading with variable-3 axis contrast: ",
      round(main$summary$loading1_target_alignment, 4), "\n", sep = "")

  out_dir <- file.path(output_dir, "model2")
  scores <- as.data.frame(main$fit$scores)
  names(scores) <- paste0("score", seq_len(ncol(scores)))
  scores$label <- main$selected_labels

  result_files <- write_common_results(
    out_dir = out_dir,
    fit = main$fit,
    main_summary = main$summary,
    path = path,
    scores = scores
  )
  result_files <- c(
    result_files,
    selected_label_counts = file.path(out_dir, "selected_label_counts.csv")
  )
  write.csv(
    as.data.frame(table(main$selected_labels)),
    result_files[["selected_label_counts"]],
    row.names = FALSE
  )
  target_labels <- c("population AGC 1", "population AGC 2")
  diagnostics <- prepare_display_diagnostics(
    out_dir = out_dir,
    main = main,
    p = p,
    target_contrasts = target_contrasts,
    target_labels = target_labels,
    bootstrap_seed = seed + 230001L
  )
  result_files <- c(result_files, diagnostics$result_files)
  frechet_target <- axis_contrast(
    3L,
    axis = 3L,
    mu = diagnostics$anchor_fits[["Frechet"]]$mu
  )
  frechet_target_contrasts <- target_contrasts_for_anchor(
    diagnostics$anchor_fits[["Frechet"]]$mu,
    cbind(c3 = frechet_target, b1 = contrast_basis_d3()[, "b1"])
  )
  frechet_outputs <- prepare_frechet_outputs(
    out_dir = out_dir,
    frechet_fit = diagnostics$anchor_fits[["Frechet"]],
    target_contrasts = frechet_target_contrasts,
    bootstrap_seed = seed + 235001L
  )
  result_files <- c(result_files, frechet_outputs$result_files)

  plot_files <- c(
    eigenvalues = file.path(out_dir, "eigenvalues.pdf"),
    variation_explained = file.path(out_dir, "anchored_variation_explained.pdf"),
    reconstruction_risk = file.path(out_dir, "reconstruction_risk.pdf"),
    scores_by_type = file.path(out_dir, "scores_by_type.pdf"),
    loading_agc1 = file.path(out_dir, "loading_agc1.pdf"),
    loading_agc2 = file.path(out_dir, "loading_agc2.pdf"),
    loading_agc1_frechet = file.path(out_dir, "loading_agc1_frechet.pdf"),
    loading_agc2_frechet = file.path(out_dir, "loading_agc2_frechet.pdf"),
    anchor_sensitivity_distance = file.path(out_dir, "anchor_sensitivity_distance.pdf"),
    oracle_alignment = file.path(out_dir, "oracle_alignment.pdf"),
    sphere_agc1_anchor_overlay = file.path(out_dir, "sphere_agc1_anchor_overlay.pdf"),
    sphere_geodesic = file.path(out_dir, "sphere_geodesic.pdf"),
    sphere_geodesic_frechet = file.path(out_dir, "sphere_geodesic_frechet.pdf")
  )

  plot_eigenvalues(main$fit, plot_files[["eigenvalues"]],
                   "Standard Model 2: AGCA eigenvalues")
  plot_display_diagnostics(
    main,
    diagnostics,
    plot_files,
    variation_title = "Model 2: anchored variation explained",
    anchor_title = "Model 2: anchor distance",
    loading_title_prefix = "Model 2",
    oracle_title = "Model 2: oracle loading alignment",
    loading_ylabs = c("", ""),
    loading_bar_width = 0.25,
    loading_plot_height = 6.5,
    loading_text_scale = 2.0,
    loading_title_scale = 2.35,
    variation_text_scale = 1.25
  )
  plot_reconstruction_risk(main$fit, plot_files[["reconstruction_risk"]],
                           "Standard Model 2: reconstruction risk")
  plot_scores(main$fit, plot_files[["scores_by_type"]],
              labels = main$selected_labels)
  plot_sphere_geodesic(
    main$fit,
    plot_files[["sphere_geodesic"]],
    "Model 2: canonical-anchor AGCs",
    labels = main$selected_labels,
    target_contrasts = target_contrasts,
    target_labels = c("population AGC 1", "population AGC 2"),
    show_estimated_second = TRUE,
    title_scale = 1.45
  )
  plot_sphere_agc1_anchor_overlay(
    main$fit,
    diagnostics$anchor_fits[["Frechet"]],
    plot_files[["sphere_agc1_anchor_overlay"]],
    "Model 2: near-axis angular regime",
    labels = main$selected_labels
  )
  plot_frechet_outputs(
    frechet_outputs,
    plot_files,
    loading_title_prefix = "Frechet M2",
    sphere_title = "Model 2: Frechet-anchor AGCs",
    target_labels = target_labels,
    labels = main$selected_labels,
    population_anchor_fit = population_anchor_fit,
    loading_ylabs = c("", ""),
    loading_bar_width = 0.25,
    loading_plot_width = 6.6,
    loading_plot_height = 6.5,
    loading_text_scale = 2.0,
    loading_title_scale = 2.35
  )

  list(
    main = main,
    path = path,
    observations = observations,
    out_dir = out_dir,
    result_files = result_files,
    plot_files = plot_files
  )
}

run_model2_adaptive_anchor <- function(observations,
                                       anchor_type = adaptive_anchor_type,
                                       population_anchor_fit = NULL) {
  k <- main_k
  k_values <- main_k_values
  p <- 1L

  main_threshold <- threshold_directions(observations$x, k = k)
  mu <- adaptive_anchor(main_threshold$g, type = anchor_type)
  target <- axis_contrast(3L, axis = 3L, mu = mu)
  target_contrasts <- target_contrasts_for_anchor(
    mu,
    cbind(c3 = target, b1 = contrast_basis_d3()[, "b1"])
  )

  main <- fit_at_threshold(
    observations$x,
    k = k,
    mu = mu,
    p = p,
    target_contrast = target_contrasts[, 1L],
    labels = observations$label
  )
  path <- threshold_path(
    observations$x,
    k_values = k_values,
    mu = mu,
    p = p,
    target_contrast = target_contrasts[, 1L],
    labels = observations$label
  )

  title_prefix <- paste0("Standard Model 2: ", anchor_type, " anchor")
  print_fit_summary(title_prefix, main$fit, p = p)
  cat("Adaptive anchor: ",
      paste(round(mu, 4), collapse = ", "), "\n", sep = "")
  cat("Selected exceedance types:\n")
  print(table(main$selected_labels))
  cat("Alignment of first loading with projected variable-3 axis contrast: ",
      round(main$summary$loading1_target_alignment, 4), "\n", sep = "")

  out_dir <- file.path(output_dir, paste0("model2_", anchor_type, "_anchor"))
  scores <- as.data.frame(main$fit$scores)
  names(scores) <- paste0("score", seq_len(ncol(scores)))
  scores$label <- main$selected_labels

  result_files <- write_common_results(
    out_dir = out_dir,
    fit = main$fit,
    main_summary = main$summary,
    path = path,
    scores = scores
  )
  result_files <- c(
    result_files,
    anchor = file.path(out_dir, "anchor.csv"),
    selected_label_counts = file.path(out_dir, "selected_label_counts.csv")
  )
  write.csv(
    anchor_table(anchor_type, mu, population_anchor_fit = population_anchor_fit),
    result_files[["anchor"]],
    row.names = FALSE
  )
  write.csv(
    as.data.frame(table(main$selected_labels)),
    result_files[["selected_label_counts"]],
    row.names = FALSE
  )
  target_labels <- c("population AGC 1", "population AGC 2")
  diagnostics <- prepare_display_diagnostics(
    out_dir = out_dir,
    main = main,
    p = p,
    target_contrasts = target_contrasts,
    target_labels = target_labels,
    bootstrap_seed = seed + 240001L
  )
  result_files <- c(result_files, diagnostics$result_files)

  plot_files <- c(
    eigenvalues = file.path(out_dir, "eigenvalues.pdf"),
    variation_explained = file.path(out_dir, "anchored_variation_explained.pdf"),
    reconstruction_risk = file.path(out_dir, "reconstruction_risk.pdf"),
    scores_by_type = file.path(out_dir, "scores_by_type.pdf"),
    loading_agc1 = file.path(out_dir, "loading_agc1.pdf"),
    loading_agc2 = file.path(out_dir, "loading_agc2.pdf"),
    anchor_sensitivity_distance = file.path(out_dir, "anchor_sensitivity_distance.pdf"),
    oracle_alignment = file.path(out_dir, "oracle_alignment.pdf"),
    sphere_geodesic = file.path(out_dir, "sphere_geodesic.pdf")
  )

  plot_eigenvalues(main$fit, plot_files[["eigenvalues"]],
                   paste0(title_prefix, ": AGCA eigenvalues"))
  plot_display_diagnostics(
    main,
    diagnostics,
    plot_files,
    variation_title = paste0("Model 2: ", anchor_type, "-anchor variation"),
    anchor_title = "Model 2: anchor distance",
    loading_title_prefix = paste0("Model 2: ", anchor_type, " anchor"),
    oracle_title = paste0("Model 2: ", anchor_type, "-anchor oracle alignment")
  )
  plot_reconstruction_risk(main$fit, plot_files[["reconstruction_risk"]],
                           paste0(title_prefix, ": reconstruction risk"))
  plot_scores(main$fit, plot_files[["scores_by_type"]],
              labels = main$selected_labels)
  plot_sphere_geodesic(
    main$fit,
    plot_files[["sphere_geodesic"]],
    paste0("Model 2: ", anchor_type, "-anchor AGCs"),
    labels = main$selected_labels,
    target_contrasts = target_contrasts,
    target_labels = c("population AGC 1", "population AGC 2"),
    show_estimated_second = TRUE,
    anchor_label = if (identical(anchor_type, "frechet")) {
      "sample Frechet anchor"
    } else {
      paste0(anchor_type, " anchor")
    },
    reference_anchors = if (is.null(population_anchor_fit)) {
      NULL
    } else {
      matrix(population_anchor_fit$anchor, nrow = 1L)
    },
    reference_anchor_labels = if (is.null(population_anchor_fit)) {
      NULL
    } else {
      "population Frechet anchor"
    }
  )

  list(
    main = main,
    path = path,
    anchor = mu,
    out_dir = out_dir,
    result_files = result_files,
    plot_files = plot_files
  )
}

fit_summary_for_sensitivity <- function(x, k, mu, target_contrast,
                                        true_space = NULL, labels = NULL) {
  thr <- threshold_directions(x, k = k)
  fit <- agca_fit(thr$g, mu = mu, p = 1L)
  variation <- agca_variation_explained(fit)

  alignment <- abs(drop(crossprod(fit$loadings[, 1L], target_contrast)))
  subspace_err <- if (is.null(true_space)) {
    sqrt(max(0, 1 - alignment^2))
  } else {
    subspace_distance(fit$loadings[, 1L, drop = FALSE], true_space)
  }

  selected_labels <- if (is.null(labels)) NULL else labels[thr$index]
  axis_fraction <- if (is.null(selected_labels)) {
    NA_real_
  } else {
    mean(selected_labels == "axis_3")
  }

  data.frame(
    k = k,
    tail_fraction = k / nrow(x),
    threshold = thr$threshold,
    eig1 = fit$eigenvalues[1L],
    eig2 = fit$eigenvalues[2L],
    eigengap = fit$eigenvalues[1L] - fit$eigenvalues[2L],
    variation_explained_1 = variation[1L],
    residual_risk_1 = agca_mean_residual(fit, p = 1L),
    loading1_target_alignment = alignment,
    subspace_distance = subspace_err,
    axis_fraction = axis_fraction
  )
}

simulate_model1_once <- function(n, k, mu, contrasts) {
  observations <- simulate_standard_model1(n)
  fit_summary_for_sensitivity(
    x = observations$x,
    k = k,
    mu = mu,
    target_contrast = contrasts[, 1L],
    true_space = contrasts[, 1L, drop = FALSE]
  )
}

simulate_model2_once <- function(n, k, mu, target) {
  observations <- simulate_standard_model2(n)
  fit_summary_for_sensitivity(
    x = observations$x,
    k = k,
    mu = mu,
    target_contrast = target,
    labels = observations$label
  )
}

quantile_row <- function(x) {
  stats <- quantile(x, probs = c(0.10, 0.50, 0.90), na.rm = TRUE, names = FALSE)
  data.frame(q10 = stats[1L], median = stats[2L], q90 = stats[3L])
}

summarize_results <- function(results) {
  metrics <- c(
    "eig1",
    "eig2",
    "eigengap",
    "variation_explained_1",
    "residual_risk_1",
    "loading1_target_alignment",
    "subspace_distance",
    "axis_fraction"
  )

  rows <- list()
  idx <- 1L
  for (model in unique(results$model)) {
    for (n in sort(unique(results$n[results$model == model]))) {
      subset <- results[results$model == model & results$n == n, , drop = FALSE]
      for (metric in metrics) {
        qs <- quantile_row(subset[[metric]])
        rows[[idx]] <- data.frame(
          model = model,
          n = n,
          k = unique(subset$k),
          metric = metric,
          qs
        )
        idx <- idx + 1L
      }
    }
  }
  do.call(rbind, rows)
}

summary_for_plot <- function(summary, model, metric) {
  out <- summary[summary$model == model & summary$metric == metric, , drop = FALSE]
  out[order(out$n), , drop = FALSE]
}

plot_metric_band <- function(summary, model, metric, ylab, main,
                             ylim = NULL, h = NULL) {
  s <- summary_for_plot(summary, model, metric)
  if (is.null(ylim)) {
    ylim <- range(c(s$q10, s$q90, h), finite = TRUE)
    pad <- 0.05 * diff(ylim)
    if (!is.finite(pad) || pad == 0) {
      pad <- 0.05 * max(1, abs(ylim[1L]))
    }
    ylim <- ylim + c(-pad, pad)
  }

  plot(
    s$n,
    s$median,
    type = "n",
    log = "x",
    xaxt = "n",
    ylim = ylim,
    xlab = "Sample size n",
    ylab = ylab,
    main = main
  )
  axis(1, at = s$n, labels = s$n)
  if (!is.null(h)) {
    abline(h = h, col = "gray80", lty = 2)
  }
  polygon(
    c(s$n, rev(s$n)),
    c(s$q10, rev(s$q90)),
    col = adjustcolor("#4C78A8", alpha.f = 0.18),
    border = NA
  )
  lines(s$n, s$median, type = "b", pch = 16, col = "#2C6AA0", lwd = 2)
  lines(s$n, s$q10, col = "#2C6AA0", lty = 3)
  lines(s$n, s$q90, col = "#2C6AA0", lty = 3)
  mtext("median with 10-90% band", side = 3, line = 0.2, cex = 0.78, col = "gray35")
}

make_recovery_plot <- function(summary, file) {
  save_pdf(file, function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4.8, 4.9, 3.3, 1.1), cex.axis = 1.05,
        cex.lab = 1.10, cex.main = 1.00)
    plot_metric_band(
      summary, "model1", "loading1_target_alignment",
      ylab = "Alignment with b1",
      main = "Standard Model 1: leading loading",
      ylim = c(0, 1),
      h = 1
    )
    plot_metric_band(
      summary, "model1", "subspace_distance",
      ylab = "Projector distance",
      main = "Standard Model 1: component-space error",
      ylim = c(0, max(summary$q90[summary$model == "model1" &
                                 summary$metric == "subspace_distance"], na.rm = TRUE) * 1.05),
      h = 0
    )
    plot_metric_band(
      summary, "model2", "loading1_target_alignment",
      ylab = "Alignment with c3",
      main = "Standard Model 2: leading loading",
      ylim = c(0, 1),
      h = 1
    )
    plot_metric_band(
      summary, "model2", "axis_fraction",
      ylab = "Selected axis-3 fraction",
      main = "Standard Model 2: selected regimes",
      ylim = c(0, 1)
    )
  }, width = 9.0, height = 7.2)
}

make_variation_plot <- function(summary, file) {
  save_pdf(file, function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4.8, 4.9, 3.3, 1.1), cex.axis = 1.05,
        cex.lab = 1.10, cex.main = 1.00)
    plot_metric_band(
      summary, "model1", "variation_explained_1",
      ylab = "Rank-one variation explained",
      main = "Standard Model 1: explained variation",
      ylim = c(0, 1)
    )
    plot_metric_band(
      summary, "model1", "residual_risk_1",
      ylab = "Rank-one residual risk",
      main = "Standard Model 1: residual risk"
    )
    plot_metric_band(
      summary, "model2", "variation_explained_1",
      ylab = "Rank-one variation explained",
      main = "Standard Model 2: explained variation",
      ylim = c(0, 1)
    )
    plot_metric_band(
      summary, "model2", "residual_risk_1",
      ylab = "Rank-one residual risk",
      main = "Standard Model 2: residual risk"
    )
  }, width = 9.0, height = 7.2)
}

run_sensitivity <- function() {
  mu <- canonical_anchor(3L)
  contrasts <- contrast_basis_d3()
  target_c3 <- axis_contrast(3L, axis = 3L, mu = mu)

  n_values <- c(300L, 600L, 1200L, 2400L, 6000L, 12000L)
  tail_fraction <- 800 / 12000
  k_values <- pmax(5L, as.integer(round(n_values * tail_fraction)))

  sensitivity_dir <- file.path(output_dir, "n_sensitivity")
  dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

  total <- length(n_values) * n_rep * 2L
  rows <- vector("list", total)
  idx <- 1L

  for (j in seq_along(n_values)) {
    n <- n_values[j]
    k <- k_values[j]
    message("Running standard n = ", n, ", k = ", k,
            " (", n_rep, " replicates per model)")
    for (rep in seq_len(n_rep)) {
      rows[[idx]] <- cbind(
        data.frame(model = "model1", n = n, rep = rep),
        simulate_model1_once(n, k, mu, contrasts)
      )
      idx <- idx + 1L

      rows[[idx]] <- cbind(
        data.frame(model = "model2", n = n, rep = rep),
        simulate_model2_once(n, k, mu, target_c3)
      )
      idx <- idx + 1L
    }
  }

  results <- do.call(rbind, rows)
  summary <- summarize_results(results)

  result_file <- file.path(sensitivity_dir, "replicate_results.csv")
  summary_file <- file.path(sensitivity_dir, "summary_by_n.csv")
  recovery_plot <- file.path(sensitivity_dir, "finite_n_recovery.pdf")
  variation_plot <- file.path(sensitivity_dir, "finite_n_variation_risk.pdf")

  write.csv(results, result_file, row.names = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)
  make_recovery_plot(summary, recovery_plot)
  make_variation_plot(summary, variation_plot)

  metadata <- data.frame(
    seed = seed,
    n_rep = n_rep,
    tail_fraction = tail_fraction,
    n = n_values,
    k = k_values,
    logistic_theta_model1 = logistic_theta_model1,
    logistic_theta_model2 = logistic_theta_model2,
    embedding_rho = embedding_rho,
    model1_third_eta = model1_third_eta,
    model1_third_scale = model1_third_scale,
    finite_tau = finite_tau,
    axis_scale = axis_scale,
    adaptive_anchor_type = adaptive_anchor_type
  )
  metadata_file <- file.path(sensitivity_dir, "metadata.csv")
  write.csv(metadata, metadata_file, row.names = FALSE)

  cat("\nStandard finite-sample sensitivity outputs:\n")
  cat("  ", result_file, "\n", sep = "")
  cat("  ", summary_file, "\n", sep = "")
  cat("  ", metadata_file, "\n", sep = "")
  cat("  ", recovery_plot, "\n", sep = "")
  cat("  ", variation_plot, "\n", sep = "")

  invisible(list(
    results = results,
    summary = summary,
    files = c(
      results = result_file,
      summary = summary_file,
      metadata = metadata_file,
      recovery_plot = recovery_plot,
      variation_plot = variation_plot
    )
  ))
}

write_design_metadata <- function(file) {
  metadata <- data.frame(
    seed = seed,
    logistic_theta_model1 = logistic_theta_model1,
    logistic_theta_model2 = logistic_theta_model2,
    embedding_rho = embedding_rho,
    model1_third_eta = model1_third_eta,
    model1_third_scale = model1_third_scale,
    finite_tau = finite_tau,
    axis_scale = axis_scale,
    adaptive_anchor_type = adaptive_anchor_type,
    display_n = main_n,
    display_k = main_k,
    bootstrap_reps = if (skip_bootstrap) 0L else bootstrap_reps,
    threshold_path_k = paste(main_k_values, collapse = ","),
    population_anchor_n = population_anchor_n,
    population_anchor_seed = population_anchor_seed,
    population_anchor_tail_fraction = main_k / main_n
  )
  write.csv(metadata, file, row.names = FALSE)
}

run_all_standard_3d_simulations <- function() {
  mu <- canonical_anchor(3L)
  contrasts <- contrast_basis_d3()

  metadata_file <- file.path(output_dir, "design_metadata.csv")
  write_design_metadata(metadata_file)

  message(
    "Estimating population Frechet anchors with n = ",
    population_anchor_n,
    " and tail fraction = ",
    signif(main_k / main_n, 4)
  )
  population_anchors <- list(
    model1 = estimate_population_frechet_anchor(
      simulate_standard_model1,
      seed_value = population_anchor_seed + 1L
    ),
    model2 = estimate_population_frechet_anchor(
      simulate_standard_model2,
      seed_value = population_anchor_seed + 2L
    )
  )

  model1 <- run_model1(
    mu,
    contrasts,
    population_anchor_fit = population_anchors$model1
  )
  model2 <- run_model2(
    mu,
    population_anchor_fit = population_anchors$model2
  )
  sensitivity <- if (skip_sensitivity) NULL else run_sensitivity()

  cat("\nStandard 3D simulation output directories:\n")
  cat("  ", model1$out_dir, "\n", sep = "")
  cat("  ", model2$out_dir, "\n", sep = "")
  cat("  ", metadata_file, "\n", sep = "")
  if (!is.null(sensitivity)) {
    cat("  ", file.path(output_dir, "n_sensitivity"), "\n", sep = "")
  }

  invisible(list(
    model1 = model1,
    model2 = model2,
    population_anchors = population_anchors,
    sensitivity = sensitivity,
    metadata_file = metadata_file
  ))
}

standard_simulation_3d_results <- run_all_standard_3d_simulations()
