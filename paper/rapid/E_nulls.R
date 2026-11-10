## Permutation-null means (fixed model, Algorithm 2, B = 199) for every cell of Tables 3 and 4,
## evaluated on release 1 of each scheme (added after the statistician review: RAPID's level
## includes a confidence component that the null removes).
source("paper/rapid/common.R")
adult <- load_adult()
train <- read.csv(here("data", "A", "train.csv"), stringsAsFactors = TRUE)
hold  <- read.csv(here("data", "A", "holdout.csv"), stringsAsFactors = TRUE)
align <- function(s) { for (v in names(train)) if (is.factor(train[[v]]))
  s[[v]] <- factor(as.character(s[[v]]), levels = levels(train[[v]])); s[names(train)] }
nullm <- function(orig, syn, qi, target) {
  t <- rapid_test(original_data = orig, synthetic_data = syn, quasi_identifiers = qi,
                  sensitive_attribute = target, model_type = "rf", cat_tau = 0.3,
                  cat_eval_method = "RCS_marginal", n_permutations = 199, seed = 1)
  c(null = mean(t$null_distribution), p = t$p_value)
}
syn1 <- syn(adult, method = "cart", m = 1, seed = 2026, print.flag = FALSE)$syn  # = B release 1
q_demo <- c("age", "sex", "race", "native.country")
q_ext  <- c(q_demo, "education.num", "workclass", "occupation", "hours.per.week")
B <- list()
for (s in c("income", "marital.status")) {
  sc <- list(demographic = q_demo, extended = q_ext, full = setdiff(names(adult), s))
  for (nm in names(sc)) { cat(s, nm, "\n")
    B[[length(B) + 1]] <- data.frame(target = s, scenario = nm, t(nullm(adult, syn1, setdiff(sc[[nm]], s), s))) }
}
B <- do.call(rbind, B); print(B, digits = 3)
qi <- setdiff(names(train), "income"); A <- list()
for (rel in c("original_1", "cart_1", "mst10_1", "mst1_1", "mst0.1_1")) { cat(rel, "\n")
  s <- align(read.csv(here("data", "A", paste0("syn_", rel, ".csv")), stringsAsFactors = TRUE))
  A[[length(A) + 1]] <- data.frame(release = rel, t(nullm(hold, s, qi, "income"))) }
A <- do.call(rbind, A); print(A, digits = 3)
saveRDS(list(B = B, A = A), here("results", "E_nulls.rds"))
