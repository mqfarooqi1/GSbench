test_that("available_models always offers gblup and ensemble", {
  am <- available_models()
  expect_true("gblup" %in% am)
  expect_true("ensemble" %in% am)
})

test_that("GBLUP marker-effect prediction reproduces the GEBVs on training", {
  sim <- simulate_population(n = 120, m = 400, seed = 1)
  fit <- gblup(sim$pheno, sim$geno)
  pred_train <- predict(fit, sim$geno)
  expect_equal(pred_train - fit$intercept, unname(fit$gebv), tolerance = 1e-6)
})

test_that("GBLUP predicts held-out individuals with signal", {
  sim <- simulate_population(n = 200, m = 600, h2 = 0.6, seed = 2)
  tr <- 1:150; te <- 151:200
  fit <- gs_fit(sim$pheno[tr], sim$geno[tr, ], model = "gblup")
  pred <- predict(fit, sim$geno[te, ])
  expect_length(pred, length(te))
  expect_gt(cor(pred, sim$bv[te]), 0.3)
})

test_that("each available ML learner fits and predicts", {
  sim <- simulate_population(n = 120, m = 300, seed = 3)
  for (mdl in setdiff(available_models(), c("gblup", "ensemble"))) {
    fit <- gs_fit(sim$pheno, sim$geno, model = mdl)
    pred <- predict(fit, sim$geno)
    expect_length(pred, nrow(sim$geno))
    expect_true(is.numeric(pred) && !anyNA(pred))
  }
})

test_that("gs_cv returns an accuracy table with signal", {
  sim <- simulate_population(n = 150, m = 400, h2 = 0.5, seed = 4)
  cv <- gs_cv(sim$pheno, sim$geno, models = "gblup", k = 5, seed = 1)
  expect_s3_class(cv, "gs_cv")
  expect_true("gblup" %in% cv$accuracy$model)
  expect_gt(cv$accuracy$mean[cv$accuracy$model == "gblup"], 0.25)
})

test_that("leave_group_out scheme works", {
  sim <- simulate_population(n = 120, m = 300, seed = 5)
  grp <- rep(1:4, length.out = 120)
  cv <- gs_cv(sim$pheno, sim$geno, models = "gblup",
              scheme = "leave_group_out", groups = grp)
  expect_equal(cv$accuracy$n_folds[1], 4)
})

test_that("stacked ensemble produces valid weights and predictions", {
  sim <- simulate_population(n = 120, m = 300, seed = 6)
  ens <- gs_ensemble(sim$pheno, sim$geno,
                     base_models = c("gblup", "elastic_net"), seed = 1)
  expect_s3_class(ens, "gs_ensemble")
  expect_equal(sum(ens$weights), 1, tolerance = 1e-6)
  expect_true(all(ens$weights >= 0))
  pred <- predict(ens, sim$geno)
  expect_length(pred, 120)
})

test_that("gs_benchmark runs across models", {
  sim <- simulate_population(n = 120, m = 300, seed = 7)
  bench <- gs_benchmark(sim$pheno, sim$geno,
                        models = c("gblup", "ensemble"), k = 3, seed = 1)
  expect_s3_class(bench, "gs_cv")
  expect_setequal(bench$accuracy$model, c("gblup", "ensemble"))
})
