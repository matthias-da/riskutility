# Builds tests/testthat/fixtures/rapid_ot_fixture.rds: reference values produced
# by Oscar Thees' RAPID package (repository RAPID_Dev, branch numeric-band,
# commit f436359). The RAPID tests of riskutility compare against these values.
#
# Run once from the package root:
#   Rscript tests/testthat/fixtures/make-rapid-ot-fixture.R
#
# Needs git, pkgload, ranger, rpart, xgboost and a clone of RAPID_Dev that
# contains the commit (default ~/workspace26/RAPID_Dev, or the environment
# variable RAPID_DEV). The clone is only read through `git archive`; its
# checkout is not touched.
#
# Not part of the test suite: testthat does not run files in fixtures/.

stopifnot(file.exists("DESCRIPTION"), file.exists("tests/testthat/helper-rapid.R"))
source("tests/testthat/helper-rapid.R")

repo <- path.expand(Sys.getenv("RAPID_DEV", "~/workspace26/RAPID_Dev"))
commit <- "f436359"
src <- file.path(tempdir(), "rapid_ot")
dir.create(src)
status <- system(sprintf("git -C %s archive %s R DESCRIPTION NAMESPACE | tar -x -C %s",
                         shQuote(repo), commit, shQuote(src)))
stopifnot(status == 0)
ot <- pkgload::load_all(src, export_all = TRUE, quiet = TRUE)$env

qi <- c("x1", "x2")
models <- c("rf", "cart", "lm", "gbm")

# ---- Data -----------------------------------------------------------------
num_o <- sim_numeric(n = 150, seed = 1, signal = 0.8)
num_s <- noisy_copy(num_o, sd = 2)                 # continues the RNG stream of seed 1
cat_o <- to_classes(num_o)
cat_s <- to_classes(num_s, ref = num_o)

o_sub <- num_o[, c("y", qi)]
s_sub <- num_s[, c("y", qi)]
y <- o_sub$y

# ---- Band method per attacker model ----------------------------------------
D <- point <- band <- rapid <- chance <- list()
for (m in models) {
  set.seed(1)
  fit <- ot$fit_model(m, y ~ x1 + x2, s_sub, o_sub, "y",
                      ranger.control = if (m == "rf") list(quantreg = TRUE) else list())
  point[[m]] <- as.numeric(ot$predict_sensitive(m, fit, o_sub, "y"))
  D[[m]] <- ot$predict_distribution(m, fit, o_sub, synthetic_data = s_sub)
  band[[m]] <- ot$band_risk(y, D[[m]], "symmetric", 0.1, 0.01)

  res <- ot$rapid(num_o, num_s, qi, "y", model_type = m, seed = 1)
  rapid[[m]] <- res$risk$confidence_rate
  chance[[m]] <- res$risk$chance_rate
}
# lm and cart are deterministic: the package result must equal the scored grid
for (m in c("lm", "cart")) {
  stopifnot(isTRUE(all.equal(rapid[[m]], mean(band[[m]]$r_hat > 0.3))))
}

# ---- Band bounds -----------------------------------------------------------
y_bounds <- c(-120, -3, -0.004, 0, 0.004, 3, 50, 50000)
bounds <- do.call(rbind, lapply(
  c("symmetric", "stabilised_relative", "absolute"),
  function(metric) {
    epsilon <- if (metric == "absolute") 2 else 0.1
    b <- ot$band_bounds(y_bounds, metric, epsilon, 0.01)
    data.frame(y = y_bounds, metric = metric, epsilon = epsilon, delta = 0.01,
               lo = b$lo, hi = b$hi, stringsAsFactors = FALSE)
  }
))

# ---- Categorical target (rf; no gbm, see D12) --------------------------------
co <- cat_o[, c("y", qi)]
cs <- cat_s[, c("y", qi)]
set.seed(1)
fit_cat <- ot$fit_model("rf", y ~ x1 + x2, cs, co, "y")
P <- as.matrix(ot$predict_sensitive("rf", fit_cat, co, "y"))
ev_cat <- ot$evaluate_categorical(co$y, P, 0.3, co, "RCS_marginal", "y",
                                  return_all_records = TRUE)

# ---- Point rule with common-value marking (lm, deterministic) ----------------
pr <- ot$rapid(num_o, num_s, qi, "y", model_type = "lm", num_eval_method = "point",
               num_common_flags = "mark", return_all_records = TRUE, seed = 1)

fixture <- list(
  num = list(o = num_o, s = num_s),
  cat = list(o = cat_o, s = cat_s, P = P, gain = ev_cat$normalized_gain,
             rapid = ev_cat$confidence_rate),
  bounds = bounds,
  D = D, point = point, band = band, rapid = rapid, chance = chance,
  point_rule = list(
    at_risk = pr$risk$rows_risk_df$at_risk,
    marginal_share = pr$risk$rows_risk_df$marginal_share,
    rare = pr$risk$rows_risk_df$rare,
    rare_hit_rate = pr$risk$rare_hit_rate,
    rare_confidence_rate = pr$risk$rare_confidence_rate,
    rapid = pr$risk$confidence_rate
  ),
  meta = list(
    source = "RAPID_Dev (Oscar Thees), branch numeric-band", commit = commit,
    settings = list(metric = "symmetric", epsilon = 0.1, delta = 0.01, tau = 0.3,
                    n_quantiles = 100, n_perm = 10, seed = 1),
    versions = vapply(c("ranger", "rpart", "xgboost"),
                      function(p) as.character(utils::packageVersion(p)), character(1)),
    r_version = R.version.string
  )
)

out <- "tests/testthat/fixtures/rapid_ot_fixture.rds"
saveRDS(fixture, out, compress = "xz")
message("Wrote ", out, " (", round(file.size(out) / 1024), " KB)")
