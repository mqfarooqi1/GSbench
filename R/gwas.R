#' Single-marker genome-wide association scan
#'
#' Tests each marker for additive association with a phenotype, optionally
#' correcting for population structure with genotype principal components. This
#' complements the whole-genome prediction models in GSbench: prediction asks
#' *how well can I predict*, whereas GWAS asks *which markers are associated*.
#'
#' For speed the phenotype and the markers are projected onto the residual space
#' of the covariates (intercept, optional principal components and extra
#' covariates), and the partial regression slope of each marker is tested with a
#' t-test. This is algebraically equivalent to fitting one linear model per
#' marker but avoids the per-marker loop.
#'
#' @param y Numeric phenotype vector (length n). Missing values are dropped.
#' @param geno Marker matrix (n x m), coded 0/1/2, with no missing values.
#' @param n_pc Number of genotype principal components to include as covariates
#'   for structure correction (default 0).
#' @param covariates Optional numeric matrix of additional fixed covariates
#'   (n x k).
#' @param map Optional data frame of per-marker annotation; columns `chrom` and
#'   `pos` (if present) are carried through and enable chromosome-aware Manhattan
#'   plots.
#' @param min_maf Drop markers with minor-allele frequency below this before
#'   testing (default 0).
#' @return An object of class `gs_gwas` (a data frame) with one row per marker:
#'   `marker`, optionally `chrom` and `pos`, `effect`, `se`, `statistic`,
#'   `p_value` and `log10p`. It has `print()` and `plot()` methods (Manhattan
#'   and QQ plots); see [plot.gs_gwas()].
#' @references Yu, J. et al. (2006) "A unified mixed-model method for association
#'   mapping that accounts for multiple levels of relatedness." Nature Genetics
#'   38, 203-208. \doi{10.1038/ng1702}
#' @seealso [gblup()] for whole-genome prediction; [plot.gs_gwas()].
#' @examples
#' sim <- simulate_population(n = 150, m = 500, n_qtl = 10, h2 = 0.6, seed = 1)
#' scan <- gwas(sim$pheno, sim$geno, n_pc = 3)
#' head(scan[order(scan$p_value), ])
#' @importFrom stats prcomp pt ppoints
#' @export
gwas <- function(y, geno, n_pc = 0, covariates = NULL, map = NULL, min_maf = 0) {
  .check_geno(geno)
  y <- as.numeric(y)
  if (length(y) != nrow(geno)) {
    stop("length(y) must equal nrow(geno).", call. = FALSE)
  }
  keep_i <- !is.na(y)
  y <- y[keep_i]
  geno <- geno[keep_i, , drop = FALSE]
  if (anyNA(geno)) {
    stop("`geno` has missing values; impute them first (see impute_markers()).",
         call. = FALSE)
  }

  p <- colMeans(geno) / 2
  if (min_maf > 0) {
    keep_m <- pmin(p, 1 - p) >= min_maf
    geno <- geno[, keep_m, drop = FALSE]
    if (!is.null(map)) map <- map[keep_m, , drop = FALSE]
  }
  n <- nrow(geno)
  m <- ncol(geno)
  if (m == 0) stop("No markers left to test.", call. = FALSE)

  C <- matrix(1, n, 1)
  if (!is.null(covariates)) {
    C <- cbind(C, as.matrix(covariates)[keep_i, , drop = FALSE])
  }
  if (n_pc > 0) {
    pcs <- stats::prcomp(geno, center = TRUE, scale. = FALSE)$x
    C <- cbind(C, pcs[, seq_len(min(n_pc, ncol(pcs))), drop = FALSE])
  }

  Cproj <- C %*% solve(crossprod(C), t(C))
  ry <- as.numeric(y - Cproj %*% y)
  RM <- geno - Cproj %*% geno
  df <- n - ncol(C) - 1
  ss <- colSums(RM^2)
  ss[ss == 0] <- NA_real_
  beta <- colSums(RM * ry) / ss
  resid_var <- colSums((ry - sweep(RM, 2, beta, "*"))^2) / df
  se <- sqrt(resid_var / ss)
  statistic <- beta / se
  p_value <- 2 * stats::pt(-abs(statistic), df)

  marker <- colnames(geno)
  if (is.null(marker)) marker <- paste0("M", seq_len(m))
  out <- data.frame(marker = marker, effect = beta, se = se,
                    statistic = statistic, p_value = p_value,
                    log10p = -log10(p_value),
                    stringsAsFactors = FALSE, row.names = NULL)
  if (!is.null(map)) {
    mp <- as.data.frame(map)
    if ("chrom" %in% names(mp)) {
      out <- cbind(out["marker"],
                   chrom = mp$chrom,
                   pos = if ("pos" %in% names(mp)) mp$pos else NA,
                   out[setdiff(names(out), "marker")])
    }
  }
  class(out) <- c("gs_gwas", "data.frame")
  attr(out, "n") <- n
  out
}

#' @describeIn gwas Print a short summary and the top associations.
#' @param x A `gs_gwas` object.
#' @param top Number of top associations to show (default 6).
#' @param ... Further arguments passed to plotting functions.
#' @export
print.gs_gwas <- function(x, top = 6, ...) {
  thr <- 0.05 / nrow(x)
  cat(sprintf("<gs_gwas: %d markers, n = %d>\n", nrow(x), attr(x, "n")))
  cat(sprintf("  %d markers pass a Bonferroni threshold (p < %.2g)\n",
              sum(x$p_value < thr, na.rm = TRUE), thr))
  cat("  top associations:\n")
  ord <- order(x$p_value)
  idx <- ord[seq_len(min(top, nrow(x)))]
  print(as.data.frame(x)[idx, , drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Plot a GWAS scan (Manhattan or QQ)
#'
#' @param x A `gs_gwas` object from [gwas()].
#' @param type `"manhattan"` (default) or `"qq"`.
#' @param alpha Genome-wide significance level for the Bonferroni line on the
#'   Manhattan plot (default 0.05).
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return The `gs_gwas` object, invisibly.
#' @importFrom graphics abline axis
#' @examples
#' sim <- simulate_population(n = 150, m = 500, n_qtl = 10, h2 = 0.6, seed = 1)
#' scan <- gwas(sim$pheno, sim$geno, n_pc = 3)
#' plot(scan)
#' plot(scan, type = "qq")
#' @export
plot.gs_gwas <- function(x, type = c("manhattan", "qq"), alpha = 0.05, ...) {
  type <- match.arg(type)
  logp <- x$log10p
  if (type == "qq") {
    o <- sort(logp[is.finite(logp)], decreasing = TRUE)
    e <- -log10(stats::ppoints(length(o)))
    plot(e, o, pch = 19, cex = 0.5, col = "#1b9e77",
         xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
         main = "GWAS QQ plot", ...)
    graphics::abline(0, 1, col = "grey60")
    return(invisible(x))
  }
  thr <- -log10(alpha / nrow(x))
  if (!is.null(x$chrom)) {
    ch <- as.character(x$chrom)
    uch <- unique(ch)
    num <- suppressWarnings(as.numeric(uch))
    uch <- uch[order(is.na(num), num, uch)]
    xx <- numeric(nrow(x)); cols <- integer(nrow(x))
    ticks <- numeric(length(uch)); off <- 0
    pos <- if (!is.null(x$pos)) as.numeric(x$pos) else NULL
    for (i in seq_along(uch)) {
      idx <- which(ch == uch[i])
      pp <- if (!is.null(pos)) pos[idx] - min(pos[idx]) else seq_along(idx)
      xx[idx] <- off + pp
      ticks[i] <- off + mean(pp)
      cols[idx] <- i %% 2
      off <- off + max(pp) + 1
    }
    plot(xx, logp, pch = 19, cex = 0.4,
         col = ifelse(cols == 0, "#2c7fb8", "#8fbcd4"),
         xaxt = "n", xlab = "chromosome", ylab = "-log10(p)",
         main = "GWAS Manhattan", ...)
    graphics::axis(1, at = ticks, labels = uch)
  } else {
    plot(seq_len(nrow(x)), logp, pch = 19, cex = 0.4, col = "#2c7fb8",
         xlab = "marker index", ylab = "-log10(p)",
         main = "GWAS Manhattan", ...)
  }
  graphics::abline(h = thr, col = "#d95f02", lty = 2)
  invisible(x)
}
