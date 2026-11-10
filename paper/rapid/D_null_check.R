## Diagnostic (statistician review): is RAPID in low-accuracy cells driven by overconfidence?
## Permutation null (fixed model, Algorithm 2) and RAPID split by true class.
source("paper/rapid/common.R")
adult <- load_adult()
train <- read.csv(here("data", "A", "train.csv"), stringsAsFactors = TRUE)
hold  <- read.csv(here("data", "A", "holdout.csv"), stringsAsFactors = TRUE)
align <- function(s, ref) { for (v in names(ref)) if (is.factor(ref[[v]]))
  s[[v]] <- factor(as.character(s[[v]]), levels = levels(ref[[v]])); s[names(ref)] }
diag1 <- function(label, orig, syn, qi, target, B = 199) {
  r <- rapid(X = orig, Y = syn, key_vars = qi, target_var = target, model_type = "rf",
             cat_eval_method = "RCS_marginal", cat_tau = 0.3, seed = 1)
  g <- r$records$true_prob; y <- orig[[target]]
  p_true <- g  # prob of true class
  set.seed(1)
  ## fixed-model null: permute labels, recompute g_i from stored class-prob matrix is not
  ## exposed, so use rapid_test (Algorithm 2) for the null mean
  t <- rapid_test(original_data = orig, synthetic_data = syn, quasi_identifiers = qi,
                  sensitive_attribute = target, model_type = "rf",
                  cat_eval_method = "RCS_marginal", n_permutations = B, seed = 1)
  byc <- tapply(r$records$at_risk, y, mean)
  data.frame(cell = label, rapid = r$rapid, null_mean = mean(t$null_distribution),
             null_q95 = quantile(t$null_distribution, .95), p = t$p_value,
             by_class = paste(sprintf("%s=%.2f", names(byc), byc), collapse = "; "),
             acc = r$model_metrics$accuracy)
}
syn_full <- syn(adult, method = "cart", m = 1, seed = 2026, print.flag = FALSE)$syn
out <- rbind(
  diag1("B income demographic", adult, syn_full, c("age","sex","race","native.country"), "income"),
  diag1("B income full", adult, syn_full, setdiff(names(adult), "income"), "income"),
  diag1("A MST eps=10 holdout", hold, align(read.csv(here("data","A","syn_mst10_1.csv"), stringsAsFactors=TRUE), train),
        setdiff(names(train), "income"), "income"),
  diag1("A MST eps=0.1 holdout", hold, align(read.csv(here("data","A","syn_mst0.1_1.csv"), stringsAsFactors=TRUE), train),
        setdiff(names(train), "income"), "income"),
  diag1("A CART holdout", hold, align(read.csv(here("data","A","syn_cart_1.csv"), stringsAsFactors=TRUE), train),
        setdiff(names(train), "income"), "income"))
print(out, digits = 3)
saveRDS(out, here("results", "D_null_check.rds"))
