#' Fit GBLUP by REML
#'
#' Fits the genomic best linear unbiased prediction model
#' `y = X b + g + e`, with `g ~ N(0, G s2u)` and `e ~ N(0, I s2e)`, estimating
#' the variance components by restricted maximum likelihood (REML). The solver
#' uses the spectral / eigendecomposition approach of Endelman (2011) (the same
#' method as `rrBLUP::mixed.solve`): the REML log-likelihood is profiled to a
#' one-dimensional optimisation over the variance ratio `lambda = s2e / s2u`.
#'
#' When `geno` is supplied, the equivalent ridge-regression marker effects are
#' also returned, so the fit can predict breeding values for *new* genotypes via
#' [predict.gblup()].
#'
#' @param y Numeric phenotype vector (length n), no missing values.
#' @param geno Optional marker matrix (n x m, coded 0/1/2). Used to build `K`
#'   and to derive marker effects for out-of-sample prediction.
#' @param K Optional n x n genomic relationship matrix. If `NULL`, built from
#'   `geno`. If you pass `K` directly (no `geno`), out-of-sample prediction is
#'   not available.
#' @param X Optional fixed-effects design matrix (n x p). Defaults to an
#'   intercept.
#' @return An object of class `gblup`: a list with `Vu`, `Ve`, `h2`, `beta`
#'   (fixed effects), `gebv` (length-n GEBVs), `lambda`, `K`, and — when `geno`
#'   was given — `marker_effects`, `marker_means` (2p) and `intercept`.
#' @references Endelman, J. B. (2011) "Ridge regression and other kernels for
#'   genomic selection with R package rrBLUP." The Plant Genome 4, 250-255.
#'   \doi{10.3835/plantgenome2011.08.0024}
#' @examples
#' sim <- simulate_population(n = 120, m = 500, seed = 1)
#' fit <- gblup(sim$pheno, sim$geno)
#' fit$h2
#' cor(fit$gebv, sim$bv)
#' @export
gblup <- function(y, geno = NULL, K = NULL, X = NULL) {
  y <- as.numeric(y)
  n <- length(y)
  if (anyNA(y)) stop("`y` must not contain missing values.", call. = FALSE)

  W <- NULL; marker_means <- NULL; c_denom <- NULL
  if (is.null(K)) {
    if (is.null(geno)) stop("Supply either `geno` or `K`.", call. = FALSE)
    .check_geno(geno)
    if (anyNA(geno)) stop("`geno` has missing values; impute first.",
                          call. = FALSE)
    p <- colMeans(geno) / 2
    marker_means <- 2 * p
    W <- sweep(geno, 2, marker_means, "-")
    c_denom <- 2 * sum(p * (1 - p))
    if (c_denom <= 0) stop("All markers monomorphic; cannot build G.",
                           call. = FALSE)
    K <- tcrossprod(W) / c_denom
  }
  if (nrow(K) != n || ncol(K) != n) {
    stop("`K` must be n x n and conformable with `y`.", call. = FALSE)
  }
  if (is.null(X)) X <- matrix(1, n, 1)
  p_fix <- ncol(X)
  stopifnot(n - p_fix > 0)

  # Spectral REML (Endelman 2011 / EMMA).
  offset <- sqrt(n)
  Hb <- K + offset * diag(n)
  XtXinv <- solve(crossprod(X))
  S <- diag(n) - X %*% tcrossprod(XtXinv, X)
  eig <- eigen(S %*% Hb %*% S, symmetric = TRUE)
  phi <- eig$values[seq_len(n - p_fix)] - offset
  phi[phi < 0] <- 0
  U <- eig$vectors[, seq_len(n - p_fix), drop = FALSE]
  omega_sq <- as.numeric(crossprod(U, y))^2

  neg_reml <- function(log_lambda) {
    lambda <- exp(log_lambda)
    (n - p_fix) * log(sum(omega_sq / (phi + lambda))) + sum(log(phi + lambda))
  }
  opt <- stats::optimize(neg_reml, interval = c(log(1e-9), log(1e9)),
                         tol = 1e-10)
  lambda <- exp(opt$minimum)

  Vu <- sum(omega_sq / (phi + lambda)) / (n - p_fix)
  Ve <- lambda * Vu

  Hinv <- solve(K + lambda * diag(n))
  XtHinv <- crossprod(X, Hinv)
  beta <- solve(XtHinv %*% X, XtHinv %*% y)
  resid_gls <- Hinv %*% (y - X %*% beta)
  gebv <- as.numeric(K %*% resid_gls)
  names(gebv) <- if (!is.null(names(y))) names(y) else rownames(K)

  out <- list(Vu = Vu, Ve = Ve, h2 = Vu / (Vu + Ve),
              beta = as.numeric(beta), gebv = gebv, lambda = lambda, K = K)
  if (!is.null(W)) {
    # Ridge-regression marker effects: gebv = W %*% marker_effects exactly.
    out$marker_effects <- as.numeric(crossprod(W, resid_gls) / c_denom)
    out$marker_means <- marker_means
    out$intercept <- as.numeric(beta[1])
  }
  structure(out, class = "gblup")
}

#' Predict breeding values for new genotypes from a GBLUP fit
#'
#' @param object A `gblup` fit created with `geno` supplied.
#' @param newgeno Marker matrix for new individuals (markers must match the
#'   training markers, in the same order).
#' @param ... Ignored.
#' @return A numeric vector of predicted values (intercept + marker effects),
#'   one per row of `newgeno`.
#' @examples
#' sim <- simulate_population(n = 150, m = 400, seed = 1)
#' fit <- gblup(sim$pheno[1:120], sim$geno[1:120, ])
#' predict(fit, sim$geno[121:150, ])
#' @export
predict.gblup <- function(object, newgeno, ...) {
  if (is.null(object$marker_effects)) {
    stop("This gblup fit has no marker effects; refit with `geno` supplied.",
         call. = FALSE)
  }
  if (ncol(newgeno) != length(object$marker_effects)) {
    stop("`newgeno` must have the same markers as the training data.",
         call. = FALSE)
  }
  Wnew <- sweep(newgeno, 2, object$marker_means, "-")
  as.numeric(Wnew %*% object$marker_effects) + object$intercept
}

#' @export
print.gblup <- function(x, ...) {
  cat("<gblup>\n")
  cat(sprintf("  Vu = %.4g, Ve = %.4g, h2 = %.3f\n", x$Vu, x$Ve, x$h2))
  cat(sprintf("  %d individuals; intercept = %.4g%s\n",
              length(x$gebv), x$beta[1],
              if (!is.null(x$marker_effects)) "; marker effects available" else ""))
  invisible(x)
}
