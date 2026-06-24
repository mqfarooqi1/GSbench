#' @keywords internal
#' @noRd
.make_folds <- function(n, k = 5, reps = 1, groups = NULL, seed = NULL) {
  build <- function() {
    if (!is.null(groups)) {
      lapply(unique(groups), function(g) which(groups == g))   # leave-group-out
    } else {
      out <- list()
      for (r in seq_len(reps)) {
        fold_id <- sample(rep_len(seq_len(k), n))
        for (j in seq_len(k)) out[[length(out) + 1L]] <- which(fold_id == j)
      }
      out
    }
  }
  if (!is.null(seed)) withr::with_seed(seed, build()) else build()
}

#' Cross-validate genomic prediction models
#'
#' Runs a breeding-relevant cross-validation and reports predictive ability
#' (the correlation between predictions and observed phenotypes on held-out
#' individuals) for each model. Two schemes are supported:
#' `"kfold"` (random k-fold, i.e. predicting untested lines, repeated `reps`
#' times) and `"leave_group_out"` (hold out one family/environment at a time,
#' via `groups`).
#'
#' Fitting is wrapped so that a model which errors on a fold records `NA` for
#' that fold rather than aborting the whole run.
#'
#' @param y Numeric phenotype vector (length n), no missing values.
#' @param geno Marker matrix (n x m, 0/1/2, no missing values).
#' @param models Models to evaluate; defaults to all available (see
#'   [available_models()]). Unavailable models are dropped with a warning.
#' @param k Number of folds for `"kfold"`. Default 5.
#' @param reps Number of repeats for `"kfold"`. Default 1.
#' @param scheme `"kfold"` or `"leave_group_out"`.
#' @param groups Grouping vector (length n) for `"leave_group_out"`.
#' @param seed Optional seed (applied via [withr::with_seed()]).
#' @param ... Passed to [gs_fit()] (model hyperparameters).
#' @return An object of class `gs_cv`: a list with `accuracy` (a data frame of
#'   `model`, `mean`, `sd`, `n_folds`), `per_fold` (the raw fold results), and
#'   the call settings.
#' @examples
#' sim <- simulate_population(n = 120, m = 400, seed = 1)
#' cv <- gs_cv(sim$pheno, sim$geno, models = "gblup", k = 5, seed = 1)
#' cv
#' @export
gs_cv <- function(y, geno, models = available_models(), k = 5, reps = 1,
                  scheme = c("kfold", "leave_group_out"), groups = NULL,
                  seed = NULL, ...) {
  scheme <- match.arg(scheme)
  y <- as.numeric(y)
  n <- length(y)
  if (anyNA(y)) stop("`y` must not contain missing values.", call. = FALSE)
  .check_geno(geno)
  if (anyNA(geno)) stop("`geno` has missing values; impute first.", call. = FALSE)
  if (nrow(geno) != n) stop("`geno` rows must match length(y).", call. = FALSE)

  avail <- available_models()
  drop <- setdiff(models, avail)
  if (length(drop) > 0) {
    warning("Skipping unavailable models: ", paste(drop, collapse = ", "),
            call. = FALSE)
  }
  models <- intersect(models, avail)
  if (length(models) == 0) stop("No usable models requested.", call. = FALSE)

  if (scheme == "leave_group_out") {
    if (is.null(groups) || length(groups) != n) {
      stop("`groups` (length n) is required for leave_group_out.", call. = FALSE)
    }
    test_sets <- .make_folds(n, groups = groups, seed = seed)
  } else {
    test_sets <- .make_folds(n, k = k, reps = reps, seed = seed)
  }

  rows <- list()
  for (i in seq_along(test_sets)) {
    test <- test_sets[[i]]
    train <- setdiff(seq_len(n), test)
    if (length(train) < 5L || length(test) < 2L) next
    for (mdl in models) {
      acc <- NA_real_; rmse <- NA_real_
      res <- tryCatch({
        fit <- gs_fit(y[train], geno[train, , drop = FALSE], model = mdl, ...)
        pred <- predict(fit, geno[test, , drop = FALSE])
        list(acc = suppressWarnings(stats::cor(pred, y[test])),
             rmse = sqrt(mean((pred - y[test])^2)))
      }, error = function(e) NULL)
      if (!is.null(res)) { acc <- res$acc; rmse <- res$rmse }
      rows[[length(rows) + 1L]] <- data.frame(fold = i, model = mdl,
                                              accuracy = acc, rmse = rmse,
                                              stringsAsFactors = FALSE)
    }
  }
  per_fold <- do.call(rbind, rows)

  agg <- lapply(split(per_fold, per_fold$model), function(d) {
    data.frame(model = d$model[1],
               mean = mean(d$accuracy, na.rm = TRUE),
               sd = stats::sd(d$accuracy, na.rm = TRUE),
               n_folds = sum(!is.na(d$accuracy)),
               stringsAsFactors = FALSE)
  })
  accuracy <- do.call(rbind, agg)
  accuracy <- accuracy[order(-accuracy$mean), , drop = FALSE]
  rownames(accuracy) <- NULL

  structure(list(accuracy = accuracy, per_fold = per_fold, scheme = scheme,
                 k = k, reps = reps, models = models),
            class = "gs_cv")
}

#' @export
print.gs_cv <- function(x, ...) {
  cat(sprintf("<gs_cv: %s>\n", x$scheme))
  if (x$scheme == "kfold") cat(sprintf("  %d-fold x %d rep(s)\n", x$k, x$reps))
  acc <- x$accuracy
  acc$mean <- round(acc$mean, 3); acc$sd <- round(acc$sd, 3)
  print(acc, row.names = FALSE)
  cat("  (accuracy = predictive ability, cor(pred, observed) on held-out data)\n")
  invisible(x)
}
