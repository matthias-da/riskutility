## common.R -- shared helpers for the EPJ DS revision experiments (A, B, C)
suppressPackageStartupMessages({library(riskutility); library(synthpop)})

here <- function(...) file.path("paper/rapid", ...)

load_adult <- function() {
  f <- here("data", "adult.csv")
  if (!file.exists(f)) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    cn <- c("age","workclass","fnlwgt","education","education.num","marital.status",
            "occupation","relationship","race","sex","capital.gain","capital.loss",
            "hours.per.week","native.country","income")
    a <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data",
                  header = FALSE, col.names = cn, strip.white = TRUE, na.strings = "?")
    a <- na.omit(a)
    a <- a[, !names(a) %in% c("fnlwgt", "education")]
    write.csv(a, f, row.names = FALSE)
  }
  a <- read.csv(f, stringsAsFactors = TRUE)
  a
}

## RAPID as defined in the paper: normalized gain vs marginal baseline, tau = 0.3
rapid_paper <- function(orig, syn, qi, target, tau = 0.3, seed = 1) {
  r <- rapid(X = orig, Y = syn, key_vars = qi, target_var = target, model_type = "rf",
             cat_eval_method = "RCS_marginal", cat_tau = tau, seed = seed)
  ci <- confint(r, method = "wilson")
  g <- r$records$normalized_gain
  data.frame(rapid = r$rapid, lo = ci[1], hi = ci[2], acc = r$model_metrics$accuracy,
             rapid_t50 = mean(g > 0.5, na.rm = TRUE), rapid_t70 = mean(g > 0.7, na.rm = TRUE))
}

## Confidence-interval overlap (Karr et al. 2006), averaged over coefficients
ci_overlap <- function(orig, syn, formula) {
  fam <- binomial()
  fo <- glm(formula, data = orig, family = fam); fs <- glm(formula, data = syn, family = fam)
  co <- intersect(names(coef(fo))[!is.na(coef(fo))], names(coef(fs))[!is.na(coef(fs))])
  co <- setdiff(co, "(Intercept)")
  ci_o <- suppressMessages(confint.default(fo))[co, , drop = FALSE]
  ci_s <- suppressMessages(confint.default(fs))[co, , drop = FALSE]
  L <- pmax(ci_o[, 1], ci_s[, 1]); U <- pmin(ci_o[, 2], ci_s[, 2])
  ov <- 0.5 * ((U - L) / (ci_o[, 2] - ci_o[, 1]) + (U - L) / (ci_s[, 2] - ci_s[, 1]))
  mean(pmax(ov, 0))
}
