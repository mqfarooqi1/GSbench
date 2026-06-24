#' Fit GBLUP by REML
#'
#' Fits the genomic best linear unbiased prediction model
#' `y = X b + g + e`, with `g ~ N(0, G s2u)` and `e ~ N(0, I s2e)`, estimating
#' the variance components by restricted maximum likelihood (REML). The solver
#' uses the spectral / eigendecomposition approach of Endelman (2011) (the same
#' method as `rrBLUP::mixed.solve`): the REML log-likelihood is profiled to a
#' one-dimensional optimisation over the variance ratio `lambda = s2e / s2u`.
#'
#' @param y Numeric phenotype vector (length n). Missing values are not allowed;
#'   subset to phenotyped individuals first.
#' @param geno Optional marker matrix (n x m, coded 0/1/2). Used to build `K`
#'   with [Gmatrix()] if `K` is not supplied.
#' @param K Optional n x n genomic relationship matrix. If `NULL`, built from
#'   `geno`.
#' @param X Optional fixed-effects design matrix (n x p). Defaults to an
#'   intercept.
#' @return An object of class `gblup`: a list with `Vu` (genetic variance), `Ve`
#'   (residual variance), `h2` (`Vu/(Vu+Ve)`), `beta` (fixed effects), `gebv`
#'   (genomic estimated breeding values, length n), and `K`.
#' @references Endelman, J. B. (2011) "Ridge regression and other kernels for
#'   genomic selection with R package rrBLUP." The Plant Genome 4, 250-255.
#'   \doi{10.3835/plantgenome2011.08.0024}
#' @examples
#' sim <- simulate_population(n = 120, m = 500, seed = 1)
#' fit <- gblup(sim$pheno, sim$geno)
#' fit$h2
#' cor(fit$gebv, sim$bv)   # accuracy vs the true breeding values
#' @export
gblup <- function(y, geno = NULL, K = NULL, X = NULL) {
  y <- as.numeric(y)
  n <- length(y)
  if (anyNA(y)) stop("`y` must not contain missing values.", call. = FALSE)
  if (is.null(K)) {
    if (is.null(geno)) stop("Supply either `geno` or `K`.", call. = FALSE)
    K <- Gmatrix(geno)
  }
  if (nrow(K) != n || ncol(K) != n) {
    stop("`K` must be n x n and conformable with `y`.", call. = FALSE)
  }
  if (is.null(X)) X <- matrix(1, n, 1)
  p <- ncol(X)
  stopifnot(n - p > 0)

  # Spectral REML (Endelman 2011 / EMMA). Offset keeps the projected matrix
  # positive definite for the eigendecomposition.
  offset <- sqrt(n)
  Hb <- K + offset * diag(n)
  XtXinv <- solve(crossprod(X))
  S <- diag(n) - X %*% tcrossprod(XtXinv, X)
  eig <- eigen(S %*% Hb %*% S, symmetric = TRUE)
  phi <- eig$values[seq_len(n - p)] - offset
  phi[phi < 0] <- 0
  U <- eig$vectors[, seq_len(n - p), drop = FALSE]
  omega_sq <- as.numeric(crossprod(U, y))^2

  neg_reml <- function(log_lambda) {
    lambda <- exp(log_lambda)
    (n - p) * log(sum(omega_sq / (phi + lambda))) + sum(log(phi + lambda))
  }
  opt <- stats::optimize(neg_reml, interval = c(log(1e-9), log(1e9)),
                         tol = 1e-10)
  lambda <- exp(opt$minimum)

  Vu <- sum(omega_sq / (phi + lambda)) / (n - p)
  Ve <- lambda * Vu

  Hinv <- solve(K + lambda * diag(n))
  XtHinv <- crossprod(X, Hinv)                       # p x n
  beta <- solve(XtHinv %*% X, XtHinv %*% y)          # GLS fixed effects
  gebv <- as.numeric(K %*% (Hinv %*% (y - X %*% beta)))
  names(gebv) <- if (!is.null(names(y))) names(y) else rownames(K)

  structure(list(Vu = Vu, Ve = Ve, h2 = Vu / (Vu + Ve),
                 beta = as.numeric(beta), gebv = gebv, lambda = lambda, K = K),
            class = "gblup")
}

#' @export
print.gblup <- function(x, ...) {
  cat("<gblup>\n")
  cat(sprintf("  Vu = %.4g, Ve = %.4g, h2 = %.3f\n", x$Vu, x$Ve, x$h2))
  cat(sprintf("  %d individuals; intercept = %.4g\n",
              length(x$gebv), x$beta[1]))
  invisible(x)
}
