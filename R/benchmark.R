#' Benchmark all available genomic prediction models
#'
#' Convenience wrapper around [gs_cv()] that evaluates every model available in
#' the session (including the stacked `"ensemble"`) under one cross-validation,
#' so they can be compared on an equal footing. Returns a `gs_cv` object; use
#' [plot()][plot.gs_cv] to draw the comparison.
#'
#' @inheritParams gs_cv
#' @return A `gs_cv` object (see [gs_cv()]).
#' @examples
#' sim <- simulate_population(n = 120, m = 400, seed = 1)
#' bench <- gs_benchmark(sim$pheno, sim$geno, models = "gblup", k = 5, seed = 1)
#' bench$accuracy
#' @export
gs_benchmark <- function(y, geno, models = available_models(), k = 5, reps = 1,
                         scheme = c("kfold", "leave_group_out"), groups = NULL,
                         seed = NULL, ...) {
  scheme <- match.arg(scheme)
  gs_cv(y, geno, models = models, k = k, reps = reps, scheme = scheme,
        groups = groups, seed = seed, ...)
}

#' Summary of a cross-validation / benchmark
#'
#' @param object A `gs_cv` object.
#' @param ... Ignored.
#' @return The accuracy data frame, invisibly.
#' @examples
#' sim <- simulate_population(n = 100, m = 300, seed = 1)
#' summary(gs_cv(sim$pheno, sim$geno, models = "gblup", seed = 1))
#' @export
summary.gs_cv <- function(object, ...) {
  print(object)
  invisible(object$accuracy)
}

#' Plot a model comparison
#'
#' Bar chart of mean predictive ability per model, with +/- 1 SD whiskers across
#' folds.
#'
#' @param x A `gs_cv` object.
#' @param ... Passed to [graphics::barplot()].
#' @return `x`, invisibly.
#' @examples
#' sim <- simulate_population(n = 120, m = 400, seed = 1)
#' \donttest{
#' plot(gs_cv(sim$pheno, sim$geno, models = "gblup", k = 5, seed = 1))
#' }
#' @export
plot.gs_cv <- function(x, ...) {
  acc <- x$accuracy[order(x$accuracy$mean), , drop = FALSE]
  m <- acc$mean
  sds <- ifelse(is.na(acc$sd), 0, acc$sd)
  names(m) <- acc$model
  top <- max(m + sds, na.rm = TRUE)
  bp <- graphics::barplot(m, horiz = TRUE, las = 1, xlim = c(0, top * 1.15),
                          xlab = "predictive ability (cor on held-out data)",
                          main = "GSbench model comparison",
                          col = "#4c8cbf", ...)
  graphics::arrows(m - sds, bp, m + sds, bp, angle = 90, code = 3,
                   length = 0.04, col = "grey30")
  invisible(x)
}
