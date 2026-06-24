#' Validate a marker matrix
#'
#' @param geno Object to check.
#' @keywords internal
#' @noRd
.check_geno <- function(geno) {
  if (!is.matrix(geno) || !is.numeric(geno)) {
    stop("`geno` must be a numeric matrix of allele dosages (0/1/2).",
         call. = FALSE)
  }
  rng <- range(geno, na.rm = TRUE)
  if (is.finite(rng[1]) && (rng[1] < 0 || rng[2] > 2)) {
    stop("`geno` values must be allele dosages in [0, 2].", call. = FALSE)
  }
  invisible(TRUE)
}

#' Impute missing marker genotypes
#'
#' Replaces missing values in each marker (column) with that marker's mean
#' dosage (i.e. twice the estimated allele frequency). Simple and fast; for
#' production work a model-based imputation upstream is preferable.
#'
#' @param geno A numeric marker matrix (individuals x markers), 0/1/2, possibly
#'   with `NA`s.
#' @return The matrix with `NA`s filled by column means. Columns that are
#'   entirely missing are filled with 0.
#' @examples
#' g <- matrix(c(0, 1, NA, 2, 2, 0), nrow = 3)
#' impute_markers(g)
#' @export
impute_markers <- function(geno) {
  .check_geno(geno)
  mu <- colMeans(geno, na.rm = TRUE)
  mu[is.nan(mu)] <- 0
  na <- which(is.na(geno), arr.ind = TRUE)
  if (nrow(na) > 0) geno[na] <- mu[na[, "col"]]
  geno
}

#' Quality-control filter for marker data
#'
#' Drops markers failing a call-rate or minor-allele-frequency threshold, and
#' (optionally) monomorphic markers, then imputes any remaining missing values.
#'
#' @param geno A numeric marker matrix (individuals x markers), coded 0/1/2.
#' @param maf Minimum minor allele frequency to keep a marker. Default 0.05.
#' @param max_missing Maximum fraction of missing calls to keep a marker.
#'   Default 0.1.
#' @param impute Whether to mean-impute remaining missing values. Default `TRUE`.
#' @return A list with `geno` (the filtered, optionally imputed matrix) and
#'   `removed` (a named integer vector counting markers dropped by each rule).
#' @examples
#' sim <- simulate_population(n = 50, m = 200, seed = 1)
#' qc <- qc_markers(sim$geno)
#' dim(qc$geno)
#' qc$removed
#' @export
qc_markers <- function(geno, maf = 0.05, max_missing = 0.1, impute = TRUE) {
  .check_geno(geno)
  n <- nrow(geno)

  call_rate <- colMeans(!is.na(geno))
  drop_missing <- call_rate < (1 - max_missing)

  p <- colMeans(geno, na.rm = TRUE) / 2
  p[is.nan(p)] <- 0
  maf_vec <- pmin(p, 1 - p)
  drop_maf <- maf_vec < maf

  drop_mono <- maf_vec == 0

  keep <- !(drop_missing | drop_maf | drop_mono)
  removed <- c(
    call_rate = sum(drop_missing),
    maf = sum(drop_maf & !drop_missing),
    monomorphic = sum(drop_mono & !drop_missing & !drop_maf)
  )
  out <- geno[, keep, drop = FALSE]
  if (impute) out <- impute_markers(out)
  list(geno = out, removed = removed)
}
