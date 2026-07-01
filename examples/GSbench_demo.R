###############################################################################
## GSbench — comprehensive demonstration script
##
## Exercises every exported function of GSbench on a small DUMMY SNP dataset:
##   simulate_population, qc_markers, impute_markers, Gmatrix, gblup,
##   predict(gblup), available_models, gs_fit, predict(gs_model), gs_cv,
##   summary(gs_cv), plot(gs_cv), gs_ensemble, predict(gs_ensemble),
##   gs_benchmark.
##
## How to use in RStudio: open this file and click "Source" (Ctrl+Shift+S),
## or step through it section by section (Ctrl+Enter).
###############################################################################

## ---- 0. Setup --------------------------------------------------------------
# install.packages("GSbench")                          # core (from CRAN)
# install.packages(c("glmnet", "ranger", "xgboost"))   # for the ML models
library(GSbench)
set.seed(1)

## ---- 1. Create a dummy SNP dataset -----------------------------------------
## Convention: rows = individuals, columns = SNP markers coded as allele
## dosage 0 / 1 / 2 (copies of one allele).
n_ind    <- 150
n_marker <- 400

freqs <- runif(n_marker, 0.05, 0.95)                   # allele frequencies
geno  <- sapply(freqs, function(p) rbinom(n_ind, 2, p))
colnames(geno) <- sprintf("SNP%03d", seq_len(n_marker))
rownames(geno) <- sprintf("ID%03d",  seq_len(n_ind))

## a quantitative trait controlled by 30 causal markers + noise (h2 ~= 0.5)
qtl     <- sample(n_marker, 30)
gv      <- scale(geno[, qtl]) %*% rnorm(30)
h2_true <- 0.5
ve      <- as.numeric(var(gv)) * (1 - h2_true) / h2_true
pheno   <- as.numeric(gv + rnorm(n_ind, sd = sqrt(ve)))

## sprinkle in missing genotypes to demonstrate QC / imputation
geno[cbind(sample(n_ind, 60, TRUE), sample(n_marker, 60, TRUE))] <- NA
cat(sprintf("Dummy data: %d individuals x %d markers; %d missing genotypes\n",
            n_ind, n_marker, sum(is.na(geno))))

## ---- 2. Quality control & imputation ---------------------------------------
## qc_markers filters by minor-allele frequency and missingness (and imputes).
geno_qc <- qc_markers(geno, maf = 0.05, max_missing = 0.10, impute = TRUE)
X <- if (is.list(geno_qc)) geno_qc$geno else geno_qc   # cleaned marker matrix
cat(sprintf("After QC: %d markers retained; missing = %d\n",
            ncol(X), sum(is.na(X))))

## impute_markers imputes missing values without filtering
geno_imp <- impute_markers(geno)
cat("After imputation only, missing =", sum(is.na(geno_imp)), "\n")

## ---- 3. Genomic relationship matrix (VanRaden 2008) ------------------------
G <- Gmatrix(X, min_maf = 0.05)
cat(sprintf("GRM: %d x %d, mean diagonal = %.3f\n",
            nrow(G), ncol(G), mean(diag(G))))
image(G, axes = FALSE, main = "Genomic relationship matrix")

## ---- 4. GBLUP by REML ------------------------------------------------------
## hold out 30 individuals to show out-of-sample prediction
test_idx  <- sample(n_ind, 30)
train_idx <- setdiff(seq_len(n_ind), test_idx)

fit <- gblup(pheno[train_idx], geno = X[train_idx, ])
print(fit)                                             # print.gblup
cat(sprintf("Estimated h2 = %.3f  (Vu = %.3f, Ve = %.3f)\n",
            fit$h2, fit$Vu, fit$Ve))
head(fit$gebv)                                         # training GEBVs

pred_gblup <- predict(fit, X[test_idx, ])              # predict.gblup
cat("GBLUP test accuracy (cor):",
    round(cor(pred_gblup, pheno[test_idx]), 3), "\n")

## ---- 5. Which models are available? ----------------------------------------
print(available_models())     # gblup, elastic_net, random_forest, xgboost, ensemble

## ---- 6. Fit a single model via the common interface ------------------------
rf <- gs_fit(pheno[train_idx], X[train_idx, ], model = "random_forest")
pred_rf <- predict(rf, X[test_idx, ])                  # predict.gs_model
cat("Random-forest test accuracy:",
    round(cor(pred_rf, pheno[test_idx]), 3), "\n")

## ---- 7. Cross-validation across models -------------------------------------
cv <- gs_cv(pheno, X,
            models = c("gblup", "elastic_net", "random_forest", "xgboost"),
            k = 5, reps = 1, scheme = "kfold", seed = 1)
print(summary(cv))            # summary.gs_cv: mean +/- sd predictive ability
plot(cv)                      # plot.gs_cv: barplot with SD whiskers

## ---- 8. Leave-group-out CV (families / environments) -----------------------
groups <- factor(rep(1:5, length.out = n_ind))         # 5 dummy groups
cv_lgo <- gs_cv(pheno, X, models = c("gblup", "elastic_net"),
                scheme = "leave_group_out", groups = groups, seed = 1)
print(summary(cv_lgo))

## ---- 9. Stacked ensemble ---------------------------------------------------
ens <- gs_ensemble(pheno[train_idx], X[train_idx, ], seed = 1)
pred_ens <- predict(ens, X[test_idx, ])                # predict.gs_ensemble
cat("Ensemble test accuracy:",
    round(cor(pred_ens, pheno[test_idx]), 3), "\n")

## ---- 10. Benchmark all models ----------------------------------------------
bench <- gs_benchmark(pheno, X, k = 5, seed = 1)
print(summary(bench))
plot(bench)

## ---- 11. Built-in simulator ------------------------------------------------
sim <- simulate_population(n = 120, m = 500, n_qtl = 50, h2 = 0.5, seed = 1)
str(sim, max.level = 1)
sim_fit <- gblup(sim$pheno, sim$geno)
cat("Simulated-data GBLUP h2 =", round(sim_fit$h2, 3), "\n")

cat("\n==== GSbench demonstration complete ====\n")
