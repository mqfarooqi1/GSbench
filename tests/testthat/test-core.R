test_that("simulate_population returns sane dimensions and heritability", {
  sim <- simulate_population(n = 120, m = 400, n_qtl = 40, h2 = 0.5, seed = 1)
  expect_equal(dim(sim$geno), c(120, 400))
  expect_length(sim$pheno, 120)
  expect_true(all(sim$geno %in% c(0, 1, 2)))
  expect_true(sim$h2 > 0.2 && sim$h2 < 0.8)        # realised h2 near target
})

test_that("simulate_population is reproducible and leaves RNG untouched", {
  a <- simulate_population(n = 50, m = 100, seed = 7)
  b <- simulate_population(n = 50, m = 100, seed = 7)
  expect_identical(a$geno, b$geno)
  expect_identical(a$pheno, b$pheno)

  set.seed(123); before <- runif(1)
  set.seed(123); simulate_population(n = 30, m = 60, seed = 99); after <- runif(1)
  expect_equal(before, after)                       # caller's stream unaffected
})

test_that("qc_markers filters and imputes", {
  sim <- simulate_population(n = 60, m = 200, seed = 2)
  g <- sim$geno
  g[, 1] <- 0                                       # monomorphic -> dropped
  g[1:30, 2] <- NA                                  # high missingness -> dropped
  g[1, 3] <- NA                                     # a stray NA -> imputed
  qc <- qc_markers(g, maf = 0.05, max_missing = 0.1)
  expect_false(anyNA(qc$geno))
  expect_lt(ncol(qc$geno), ncol(g))
  expect_true(qc$removed["monomorphic"] >= 1 || qc$removed["maf"] >= 1)
})

test_that("Gmatrix is symmetric, n x n, with the right scale", {
  sim <- simulate_population(n = 80, m = 500, seed = 3)
  G <- Gmatrix(sim$geno)
  expect_equal(dim(G), c(80, 80))
  expect_equal(G, t(G))
  expect_equal(mean(diag(G)), 1, tolerance = 0.15)  # VanRaden diag averages ~1
})

test_that("gblup runs and recovers signal", {
  sim <- simulate_population(n = 150, m = 800, h2 = 0.6, seed = 4)
  fit <- gblup(sim$pheno, sim$geno)
  expect_s3_class(fit, "gblup")
  expect_true(fit$h2 > 0 && fit$h2 < 1)
  expect_length(fit$gebv, 150)
  expect_gt(cor(fit$gebv, sim$bv), 0.4)             # GEBVs track true BVs
})

test_that("gblup matches rrBLUP::mixed.solve", {
  skip_if_not_installed("rrBLUP")
  sim <- simulate_population(n = 150, m = 800, seed = 5)
  G <- Gmatrix(sim$geno)
  fit <- gblup(sim$pheno, K = G)
  ref <- rrBLUP::mixed.solve(sim$pheno, K = G)
  expect_equal(fit$Vu, ref$Vu, tolerance = 1e-3)
  expect_equal(fit$Ve, ref$Ve, tolerance = 1e-3)
  expect_gt(cor(fit$gebv, ref$u), 0.9999)
})

test_that("gblup validates its inputs", {
  sim <- simulate_population(n = 40, m = 100, seed = 6)
  expect_error(gblup(sim$pheno), "Supply either")
  expect_error(gblup(c(sim$pheno, NA), sim$geno), "missing")
})
