#' Simulate a genomic prediction dataset
#'
#' Generates a small marker matrix and a phenotype with a known narrow-sense
#' heritability, for examples, tests and demonstrations. Markers are coded as
#' allele dosages (0, 1, 2). A random subset of markers are QTL with additive
#' effects; the phenotype is the resulting breeding value plus normal noise
#' scaled to the target heritability.
#'
#' @param n Number of individuals. Default 200.
#' @param m Number of markers. Default 1000.
#' @param n_qtl Number of markers with a non-zero (QTL) effect. Default 50.
#' @param h2 Target narrow-sense heritability in (0, 1]. Default 0.5.
#' @param seed Optional integer seed for reproducibility. Applied with
#'   [withr::with_seed()], which scopes the seed locally and leaves the caller's
#'   random-number state unchanged.
#' @return A list with `geno` (an `n` x `m` 0/1/2 matrix, row/column names set),
#'   `pheno` (length-`n` numeric), `bv` (true breeding values), `qtl` (indices of
#'   the causal markers) and `h2` (the realised heritability).
#' @examples
#' sim <- simulate_population(n = 100, m = 300, seed = 1)
#' dim(sim$geno)
#' length(sim$pheno)
#' @export
simulate_population <- function(n = 200, m = 1000, n_qtl = 50, h2 = 0.5,
                                seed = NULL) {
  if (!is.null(seed)) {
    return(withr::with_seed(
      seed, simulate_population(n = n, m = m, n_qtl = n_qtl, h2 = h2)))
  }
  stopifnot(n >= 2, m >= 2, n_qtl >= 1, n_qtl <= m, h2 > 0, h2 <= 1)

  freq <- stats::runif(m, 0.05, 0.95)
  geno <- vapply(freq, function(p) stats::rbinom(n, 2, p), numeric(n))
  dimnames(geno) <- list(paste0("ind", seq_len(n)), paste0("snp", seq_len(m)))

  qtl <- sort(sample.int(m, n_qtl))
  effects <- stats::rnorm(n_qtl)
  bv <- as.numeric(scale(geno[, qtl, drop = FALSE]) %*% effects)
  bv[is.na(bv)] <- 0
  var_g <- stats::var(bv)
  var_e <- if (var_g > 0) var_g * (1 - h2) / h2 else 1
  pheno <- bv + stats::rnorm(n, 0, sqrt(var_e))
  names(pheno) <- rownames(geno)
  names(bv) <- rownames(geno)

  list(geno = geno, pheno = pheno, bv = bv, qtl = qtl,
       h2 = var_g / (var_g + var_e))
}
