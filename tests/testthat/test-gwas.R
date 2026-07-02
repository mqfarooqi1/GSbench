test_that("Dmatrix builds a valid dominance GRM usable in gblup", {
  sim <- simulate_population(n = 80, m = 300, seed = 1)
  D <- Dmatrix(sim$geno)
  expect_equal(dim(D), c(80, 80))
  expect_true(isSymmetric(unname(D), tol = 1e-8))
  fit <- gblup(sim$pheno, K = D)
  expect_true(fit$h2 >= 0 && fit$h2 <= 1)
})

test_that("gwas returns a valid scan and detects signal", {
  sim <- simulate_population(n = 200, m = 500, n_qtl = 10, h2 = 0.6, seed = 2)
  scan <- gwas(sim$pheno, sim$geno, n_pc = 3)
  expect_s3_class(scan, "gs_gwas")
  expect_equal(nrow(scan), 500)
  expect_true(all(scan$p_value >= 0 & scan$p_value <= 1, na.rm = TRUE))
  causal <- rep(FALSE, 500); causal[sim$qtl] <- TRUE
  expect_gt(mean(scan$log10p[causal]), mean(scan$log10p[!causal]))
  expect_gte(length(intersect(order(scan$p_value)[1:5], sim$qtl)), 1)
})

test_that("gwas map enables chrom columns and print/plot work", {
  sim <- simulate_population(n = 120, m = 200, n_qtl = 8, h2 = 0.5, seed = 3)
  map <- data.frame(chrom = rep(1:2, each = 100), pos = rep(seq_len(100), 2))
  scan <- gwas(sim$pheno, sim$geno, n_pc = 2, map = map)
  expect_true(all(c("chrom", "pos") %in% names(scan)))
  expect_output(print(scan), "gs_gwas")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_invisible(plot(scan))
  expect_invisible(plot(scan, type = "qq"))
})

test_that("gwas drops missing phenotypes", {
  sim <- simulate_population(n = 100, m = 150, seed = 4)
  y <- sim$pheno; y[1:10] <- NA
  scan <- gwas(y, sim$geno)
  expect_equal(nrow(scan), 150)
  expect_equal(attr(scan, "n"), 90)
})
