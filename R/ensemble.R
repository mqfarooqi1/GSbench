#' @keywords internal
#' @noRd
.nnls_weights <- function(P, y) {
  ok <- stats::complete.cases(P) & !is.na(y)
  P <- P[ok, , drop = FALSE]; y <- y[ok]
  k <- ncol(P)
  sse <- function(w) sum((y - P %*% w)^2)
  opt <- stats::optim(rep(1 / k, k), sse, method = "L-BFGS-B",
                      lower = rep(0, k))
  w <- pmax(opt$par, 0)
  if (sum(w) > 0) w / sum(w) else rep(1 / k, k)
}

#' Stacked super-learner ensemble of genomic prediction models
#'
#' Combines several base models into one predictor by *stacking*: each base
#' model's out-of-fold cross-validated predictions are used to learn a set of
#' non-negative weights (constrained to sum to one), and the final prediction is
#' that weighted average of the base models refit on all the data. This is the
#' Breiman / van der Laan stacked-regression (super-learner) idea applied to
#' genomic selection; in practice it tends to match or beat the best single
#' model without having to know in advance which that is.
#'
#' @param y Numeric phenotype vector (length n), no missing values.
#' @param geno Marker matrix (n x m, 0/1/2, no missing values).
#' @param base_models Character vector of base model names. Defaults to every
#'   available model except the ensemble itself.
#' @param inner_k Folds for the inner stacking cross-validation. Default 5.
#' @param seed Optional seed for the inner folds (via [withr::with_seed()]).
#' @param ... Passed to [gs_fit()] for the base models.
#' @return An object of class `gs_ensemble` (and `gs_model`): a list with
#'   `base_names`, `weights` (named, summing to 1), the refit `base_fits`, and
#'   the out-of-fold prediction matrix `oof`.
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2007)
#'   "Super Learner." Statistical Applications in Genetics and Molecular Biology
#'   6, Article 25. \doi{10.2202/1544-6115.1309}
#' @examples
#' sim <- simulate_population(n = 100, m = 300, seed = 1)
#' ens <- gs_ensemble(sim$pheno, sim$geno, base_models = "gblup", seed = 1)
#' ens$weights
#' @export
gs_ensemble <- function(y, geno, base_models = NULL, inner_k = 5,
                        seed = NULL, ...) {
  y <- as.numeric(y); n <- length(y)
  if (anyNA(y)) stop("`y` must not contain missing values.", call. = FALSE)
  .check_geno(geno)
  if (is.null(base_models)) base_models <- setdiff(available_models(), "ensemble")
  base_models <- setdiff(base_models, "ensemble")
  base_models <- intersect(base_models, available_models())
  if (length(base_models) == 0) stop("No usable base models.", call. = FALSE)

  folds <- .make_folds(n, k = inner_k, reps = 1, seed = seed)
  P <- matrix(NA_real_, n, length(base_models),
              dimnames = list(NULL, base_models))
  for (test in folds) {
    train <- setdiff(seq_len(n), test)
    for (j in seq_along(base_models)) {
      pred <- tryCatch({
        f <- gs_fit(y[train], geno[train, , drop = FALSE],
                    model = base_models[j], ...)
        predict(f, geno[test, , drop = FALSE])
      }, error = function(e) rep(NA_real_, length(test)))
      P[test, j] <- pred
    }
  }

  weights <- .nnls_weights(P, y)
  names(weights) <- base_models

  base_fits <- lapply(base_models, function(m) gs_fit(y, geno, model = m, ...))
  names(base_fits) <- base_models

  structure(list(model = "ensemble", base_names = base_models,
                 weights = weights, base_fits = base_fits, oof = P),
            class = c("gs_ensemble", "gs_model"))
}

#' @rdname predict.gs_model
#' @export
predict.gs_ensemble <- function(object, newgeno, ...) {
  preds <- vapply(object$base_fits, function(f) predict(f, newgeno),
                  numeric(nrow(newgeno)))
  as.numeric(preds %*% object$weights)
}

#' @export
print.gs_ensemble <- function(x, ...) {
  cat("<gs_ensemble> stacked super-learner\n")
  w <- round(x$weights, 3)
  for (nm in names(w)) cat(sprintf("  %-14s weight %.3f\n", nm, w[nm]))
  invisible(x)
}
