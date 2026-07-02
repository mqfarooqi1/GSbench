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

#' Dominance genomic relationship matrix
#'
#' Builds the dominance genomic relationship matrix from allele dosages using the
#' genotypic parameterisation of Vitezica et al. (2013): homozygotes and
#' heterozygotes are coded so the matrix captures dominance (heterozygote versus
#' homozygote) deviations, complementing the additive relationships from
#' [Gmatrix()]. Supplying it to [gblup()] as `K` fits a dominance GBLUP model,
#' and the additive and dominance GEBVs can be compared.
#'
#' @inheritParams Gmatrix
#' @return An `n` x `n` symmetric dominance relationship matrix, with row/column
#'   names taken from `rownames(geno)`.
#' @references Vitezica, Z. G., Varona, L. & Legarra, A. (2013) "On the additive
#'   and dominant variance and covariance of individuals within the genomic
#'   selection scope." Genetics 195, 1223-1230. \doi{10.1534/genetics.113.155176}
#' @seealso [Gmatrix()], [gblup()]
#' @examples
#' sim <- simulate_population(n = 80, m = 400, seed = 1)
#' D <- Dmatrix(sim$geno)
#' dim(D)
#' # dominance GBLUP:
#' fit <- gblup(sim$pheno, K = D)
#' @export
Dmatrix <- function(geno, min_maf = 0) {
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
  if (ncol(geno) == 0) stop("No markers left to build D.", call. = FALSE)
  q <- 1 - p
  n <- nrow(geno)
  # genotypic dominance coding (Vitezica et al. 2013):
  #   AA (dose 2) -> -2 q^2 ; Aa (dose 1) -> 2 p q ; aa (dose 0) -> -2 p^2
  AA <- matrix(rep(-2 * q^2,  each = n), n)
  Aa <- matrix(rep( 2 * p * q, each = n), n)
  aa <- matrix(rep(-2 * p^2,  each = n), n)
  H <- AA * (geno == 2) + Aa * (geno == 1) + aa * (geno == 0)
  denom <- sum((2 * p * q)^2)
  if (denom <= 0) stop("All markers are monomorphic; cannot build D.",
                       call. = FALSE)
  D <- tcrossprod(H) / denom
  dimnames(D) <- list(rownames(geno), rownames(geno))
  D
}
