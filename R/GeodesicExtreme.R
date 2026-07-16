# Author: Alberto Quaini

# Anchored geodesic analysis for Euclidean-spherical extreme directions.
# The functions in this file are intentionally dependency-free.

row_norms <- function(x) {
  x <- as.matrix(x)
  sqrt(rowSums(x^2))
}

unit_vector <- function(x, name = "x", eps = sqrt(.Machine$double.eps)) {
  x <- as.numeric(x)
  nrm <- sqrt(sum(x^2))
  if (!is.finite(nrm) || nrm <= eps) {
    stop(name, " must have positive Euclidean norm.", call. = FALSE)
  }
  x / nrm
}

normalize_rows <- function(x, eps = sqrt(.Machine$double.eps)) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  nrms <- row_norms(x)
  if (any(!is.finite(nrms)) || any(nrms <= eps)) {
    stop("All rows must have positive finite Euclidean norm.", call. = FALSE)
  }
  sweep(x, 1L, nrms, "/")
}

canonical_anchor <- function(d) {
  if (length(d) != 1L || d < 2L) {
    stop("d must be an integer at least 2.", call. = FALSE)
  }
  rep(1 / sqrt(d), d)
}

tangent_basis <- function(mu) {
  mu <- unit_vector(mu, "mu")
  d <- length(mu)
  if (d < 2L) {
    stop("The ambient dimension must be at least 2.", call. = FALSE)
  }

  sv <- svd(matrix(mu, nrow = 1L), nu = 0L, nv = d)
  basis <- sv$v[, 2L:d, drop = FALSE]
  basis <- qr.Q(qr(basis), complete = FALSE)
  basis
}

project_to_tangent <- function(x, mu) {
  x <- as.matrix(x)
  mu <- unit_vector(mu, "mu")
  if (ncol(x) != length(mu)) {
    stop("x and mu have incompatible dimensions.", call. = FALSE)
  }

  anchor_coordinate <- drop(x %*% mu)
  x - tcrossprod(anchor_coordinate, mu)
}

anchored_departures <- function(g, mu, normalize = TRUE) {
  g <- as.matrix(g)
  if (normalize) {
    g <- normalize_rows(g)
  }

  mu <- unit_vector(mu, "mu")
  if (ncol(g) != length(mu)) {
    stop("g and mu have incompatible dimensions.", call. = FALSE)
  }

  anchor_coordinate <- drop(g %*% mu)
  if (any(anchor_coordinate <= 0)) {
    stop("All directions must lie in the open hemisphere centered at mu.", call. = FALSE)
  }

  u <- g - tcrossprod(anchor_coordinate, mu)
  list(g = g, mu = mu, a = anchor_coordinate, u = u)
}

gnomonic_inverse <- function(g, mu, normalize = TRUE) {
  dep <- anchored_departures(g, mu, normalize = normalize)
  sweep(dep$u, 1L, dep$a, "/")
}

gnomonic_map <- function(v, mu) {
  v <- as.matrix(v)
  mu <- unit_vector(mu, "mu")
  if (ncol(v) != length(mu)) {
    stop("v and mu have incompatible dimensions.", call. = FALSE)
  }

  normalize_rows(v + matrix(mu, nrow(v), length(mu), byrow = TRUE))
}

sphere_geodesic_distance <- function(x, y, normalize = TRUE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (normalize) {
    x <- normalize_rows(x)
    y <- normalize_rows(y)
  }
  if (ncol(x) != ncol(y)) {
    stop("x and y have incompatible dimensions.", call. = FALSE)
  }
  if (nrow(y) == 1L && nrow(x) > 1L) {
    y <- matrix(y, nrow(x), ncol(y), byrow = TRUE)
  }
  if (nrow(x) != nrow(y)) {
    stop("x and y must have the same number of rows, unless y has one row.", call. = FALSE)
  }

  inner <- rowSums(x * y)
  acos(pmax(-1, pmin(1, inner)))
}

orient_columns <- function(x) {
  x <- as.matrix(x)
  for (j in seq_len(ncol(x))) {
    i <- which.max(abs(x[, j]))
    if (x[i, j] < 0) {
      x[, j] <- -x[, j]
    }
  }
  x
}

validate_rank <- function(p, max_rank) {
  if (length(p) != 1L || is.na(p) || p < 0L || p > max_rank || p != as.integer(p)) {
    stop("p must be an integer between 0 and ", max_rank, ".", call. = FALSE)
  }
  as.integer(p)
}

agca_fit <- function(g, mu = canonical_anchor(ncol(as.matrix(g))), p = NULL,
                     normalize = TRUE) {
  dep <- anchored_departures(g, mu, normalize = normalize)
  g <- dep$g
  mu <- dep$mu
  u <- dep$u
  n <- nrow(g)
  d <- ncol(g)

  basis <- tangent_basis(mu)
  tangent_scores <- u %*% basis
  sigma_tangent <- crossprod(tangent_scores) / n
  eig <- eigen(sigma_tangent, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)

  eigenvalues <- eig$values[ord]
  eigenvalues[abs(eigenvalues) < 100 * .Machine$double.eps] <- 0
  eigenvalues <- pmax(eigenvalues, 0)

  loadings <- basis %*% eig$vectors[, ord, drop = FALSE]
  loadings <- orient_columns(loadings)
  scores <- u %*% loadings

  if (is.null(p)) {
    p <- d - 1L
  }
  p <- validate_rank(p, d - 1L)

  structure(
    list(
      g = g,
      mu = mu,
      anchor_coordinate = dep$a,
      u = u,
      basis = basis,
      sigma = tcrossprod(basis %*% sigma_tangent, basis),
      sigma_tangent = sigma_tangent,
      eigenvalues = eigenvalues,
      loadings = loadings,
      scores = scores,
      p = p
    ),
    class = "agca_fit"
  )
}

agca_variation_explained <- function(fit) {
  total <- sum(fit$eigenvalues)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(fit$eigenvalues)))
  }
  cumsum(fit$eigenvalues) / total
}

agca_reconstruct <- function(fit, p = fit$p) {
  p <- validate_rank(p, length(fit$eigenvalues))
  n <- nrow(fit$g)
  d <- ncol(fit$g)

  if (p == 0L) {
    u_hat <- matrix(0, n, d)
  } else {
    u_hat <- fit$scores[, seq_len(p), drop = FALSE] %*%
      t(fit$loadings[, seq_len(p), drop = FALSE])
  }

  numerator <- tcrossprod(fit$anchor_coordinate, fit$mu) + u_hat
  normalize_rows(numerator)
}

agca_residuals <- function(fit, p = fit$p) {
  p <- validate_rank(p, length(fit$eigenvalues))
  n <- nrow(fit$g)
  d <- ncol(fit$g)

  if (p == 0L) {
    u_hat <- matrix(0, n, d)
  } else {
    u_hat <- fit$scores[, seq_len(p), drop = FALSE] %*%
      t(fit$loadings[, seq_len(p), drop = FALSE])
  }

  rowSums((fit$u - u_hat)^2)
}

agca_mean_residual <- function(fit, p = fit$p) {
  mean(agca_residuals(fit, p = p))
}

agca_rank_summary <- function(fit) {
  max_rank <- length(fit$eigenvalues)
  residual_risk <- vapply(
    0L:max_rank,
    function(p) agca_mean_residual(fit, p = p),
    numeric(1L)
  )

  data.frame(
    rank = 0L:max_rank,
    residual_risk = residual_risk,
    variation_explained = c(0, agca_variation_explained(fit))
  )
}

threshold_directions <- function(x, k = NULL, threshold = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  radius <- row_norms(x)

  if (!is.null(k) && !is.null(threshold)) {
    stop("Specify either k or threshold, not both.", call. = FALSE)
  }
  if (is.null(k) && is.null(threshold)) {
    stop("Specify k or threshold.", call. = FALSE)
  }

  if (!is.null(k)) {
    if (length(k) != 1L || k < 1L || k > nrow(x) || k != as.integer(k)) {
      stop("k must be an integer between 1 and nrow(x).", call. = FALSE)
    }
    index <- order(-radius, seq_along(radius))[seq_len(k)]
    threshold <- min(radius[index])
  } else {
    index <- which(radius > threshold)
  }

  if (length(index) == 0L) {
    stop("No observations exceed the threshold.", call. = FALSE)
  }

  list(
    g = normalize_rows(x[index, , drop = FALSE]),
    radius = radius[index],
    index = index,
    threshold = threshold,
    all_radius = radius
  )
}

orthonormalize_columns <- function(x, tol = 1e-10) {
  x <- as.matrix(x)
  if (ncol(x) == 0L) {
    return(matrix(numeric(0L), nrow(x), 0L))
  }

  qrx <- qr(x, tol = tol)
  rank <- qrx$rank
  if (rank == 0L) {
    return(matrix(numeric(0L), nrow(x), 0L))
  }

  qr.Q(qrx, complete = FALSE)[, seq_len(rank), drop = FALSE]
}

projection_matrix <- function(x) {
  q <- orthonormalize_columns(x)
  if (ncol(q) == 0L) {
    return(matrix(0, nrow(x), nrow(x)))
  }
  q %*% t(q)
}

subspace_distance <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (nrow(x) != nrow(y)) {
    stop("Subspaces must be represented in the same ambient dimension.", call. = FALSE)
  }

  px <- projection_matrix(x)
  py <- projection_matrix(y)
  diff <- (px - py + t(px - py)) / 2
  max(abs(eigen(diff, symmetric = TRUE, only.values = TRUE)$values))
}

principal_angle_cosines <- function(x, y) {
  qx <- orthonormalize_columns(x)
  qy <- orthonormalize_columns(y)
  if (nrow(qx) != nrow(qy)) {
    stop("Subspaces must be represented in the same ambient dimension.", call. = FALSE)
  }
  if (ncol(qx) == 0L || ncol(qy) == 0L) {
    return(numeric(0L))
  }
  pmin(1, pmax(0, svd(t(qx) %*% qy, nu = 0L, nv = 0L)$d))
}
