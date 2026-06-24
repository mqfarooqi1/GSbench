#' Genomic relationship matrix (VanRaden)
#'
#' Builds the additive genomic relationship matrix from allele dosages using
#' VanRaden's first method: with markers centred by twice the allele frequency,
#' `G = W W' / (2 * sum(p (1 - p)))`. This is the standard additive GRM used in
#' GBLUP.
#'
#' @param geno A numeric marker matrix (individuals x markers), coded 0/1/2,
#'   with no missing values (run [qc_markers()] or [impute_markers()] first).
#' @param min_maf Markers with minor allele frequency below this are dropped
#'   before building the matrix. Default 0 (keep all supplied markers).
#' @return An `n` x `n` symmetric genomic relationship matrix with the row/column
#'   names taken from `rownames(geno)`.
#' @references VanRaden, P. M. (2008) "Efficient methods to compute genomic
#'   predictions." Journal of Dairy Science 91, 4414-4423.
#'   \doi{10.3168/jds.2007-0980}
#' @examples
#' sim <- simulate_population(n = 80, m = 400, seed = 1)
#' G <- Gmatrix(sim$geno)
#' dim(G)
#' @export
Gmatrix <- function(geno, min_maf = 0) {
  .check_geno(geno)
  if (anyNA(geno)) {
    stop("`geno` has missing values; impute them first (see impute_markers()).",
         call. = FALSE)
  }
  p <- colMeans(geno) / 2
  if (min_maf > 0) {
    keep <- pmin(p, 1 - p) >= min_maf
    geno <- geno[, keep, drop = FALSE]
    p <- p[keep]
  }
  if (ncol(geno) == 0) stop("No markers left to build G.", call. = FALSE)

  W <- sweep(geno, 2, 2 * p, "-")          # centre by 2p
  denom <- 2 * sum(p * (1 - p))
  if (denom <= 0) stop("All markers are monomorphic; cannot build G.",
                       call. = FALSE)
  G <- tcrossprod(W) / denom
  dimnames(G) <- list(rownames(geno), rownames(geno))
  G
}
