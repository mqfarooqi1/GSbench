#' @keywords internal
#' @noRd
.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Model needs the '%s' package. Install it with install.packages('%s').",
                 pkg, pkg), call. = FALSE)
  }
}

#' Models available in this session
#'
#' Returns the genomic-prediction models GSbench can currently run. `"gblup"` is
#' always available; the machine-learning models require their (suggested)
#' package to be installed.
#'
#' @return A character vector of usable model names.
#' @examples
#' available_models()
#' @export
available_models <- function() {
  m <- c("gblup")
  if (requireNamespace("glmnet", quietly = TRUE)) m <- c(m, "elastic_net")
  if (requireNamespace("ranger", quietly = TRUE)) m <- c(m, "random_forest")
  if (requireNamespace("xgboost", quietly = TRUE)) m <- c(m, "xgboost")
  c(m, "ensemble")
}

# Internal learner registry: each entry has fit(X, y, ...) and predict(fit, Xnew).
.learners <- list(
  gblup = list(
    fit = function(X, y, ...) gblup(y, geno = X),
    predict = function(fit, Xnew) predict(fit, Xnew)
  ),
  elastic_net = list(
    fit = function(X, y, alpha = 0.5, nfolds = 5, ...) {
      .need("glmnet")
      glmnet::cv.glmnet(as.matrix(X), y, alpha = alpha, nfolds = nfolds)
    },
    predict = function(fit, Xnew) {
      as.numeric(stats::predict(fit, newx = as.matrix(Xnew), s = "lambda.min"))
    }
  ),
  random_forest = list(
    fit = function(X, y, num.trees = 500, ...) {
      .need("ranger")
      ranger::ranger(x = as.matrix(X), y = y, num.trees = num.trees,
                     num.threads = 1)
    },
    predict = function(fit, Xnew) {
      stats::predict(fit, data = as.matrix(Xnew))$predictions
    }
  ),
  xgboost = list(
    fit = function(X, y, nrounds = 150, eta = 0.3, max_depth = 3,
                   subsample = 0.8, ...) {
      .need("xgboost")
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
      xgboost::xgb.train(
        params = list(eta = eta, max_depth = max_depth, subsample = subsample,
                      objective = "reg:squarederror", nthread = 1),
        data = dtrain, nrounds = nrounds, verbose = 0)
    },
    predict = function(fit, Xnew) {
      as.numeric(stats::predict(fit, xgboost::xgb.DMatrix(as.matrix(Xnew))))
    }
  )
)

#' Fit a genomic prediction model
#'
#' One interface for several model families: `"gblup"` (always available),
#' `"elastic_net"` (\pkg{glmnet}), `"random_forest"` (\pkg{ranger}),
#' `"xgboost"` (\pkg{xgboost}), and `"ensemble"` (a stacked super-learner; see
#' [gs_ensemble()]). The returned object has a [predict()][predict.gs_model]
#' method that takes a new marker matrix.
#'
#' @param y Numeric phenotype vector (length n), no missing values.
#' @param geno Marker matrix (n x m, coded 0/1/2, no missing values).
#' @param model Model name; see [available_models()].
#' @param ... Model-specific hyperparameters (e.g. `alpha` for elastic net,
#'   `num.trees` for random forest, `nrounds`/`eta`/`max_depth` for xgboost,
#'   `base_models` for the ensemble).
#' @return An object of class `gs_model` wrapping the fitted model.
#' @examples
#' sim <- simulate_population(n = 120, m = 400, seed = 1)
#' fit <- gs_fit(sim$pheno, sim$geno, model = "gblup")
#' head(predict(fit, sim$geno))
#' @export
gs_fit <- function(y, geno, model = "gblup", ...) {
  y <- as.numeric(y)
  if (anyNA(y)) stop("`y` must not contain missing values.", call. = FALSE)
  .check_geno(geno)
  if (anyNA(geno)) stop("`geno` has missing values; impute first.", call. = FALSE)
  if (nrow(geno) != length(y)) {
    stop("`geno` rows must match length(y).", call. = FALSE)
  }
  if (identical(model, "ensemble")) {
    return(gs_ensemble(y, geno, ...))
  }
  if (!model %in% names(.learners)) {
    stop(sprintf("Unknown model %s. See available_models().", sQuote(model)),
         call. = FALSE)
  }
  fitted <- .learners[[model]]$fit(geno, y, ...)
  structure(list(model = model, fitted = fitted), class = "gs_model")
}

#' Predict from a fitted genomic prediction model
#'
#' @param object A `gs_model` from [gs_fit()].
#' @param newgeno Marker matrix for the individuals to predict (same markers as
#'   training).
#' @param ... Ignored.
#' @return A numeric vector of predictions, one per row of `newgeno`.
#' @examples
#' sim <- simulate_population(n = 120, m = 400, seed = 1)
#' fit <- gs_fit(sim$pheno[1:90], sim$geno[1:90, ], model = "gblup")
#' predict(fit, sim$geno[91:120, ])
#' @export
predict.gs_model <- function(object, newgeno, ...) {
  .learners[[object$model]]$predict(object$fitted, newgeno)
}

#' @export
print.gs_model <- function(x, ...) {
  cat(sprintf("<gs_model: %s>\n", x$model))
  invisible(x)
}
