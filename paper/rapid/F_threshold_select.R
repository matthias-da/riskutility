## tau*_alpha (Algorithm 3) for the six scenario cells of Table 3, release 1 (B = 199).
## Reported in Section 5.1 / 7.3: tau* reaches the grid maximum 0.95 in all cells.
source("paper/rapid/common.R")
adult <- load_adult()
syn1 <- syn(adult, method = "cart", m = 1, seed = 2026, print.flag = FALSE)$syn
q_demo <- c("age", "sex", "race", "native.country")
q_ext  <- c(q_demo, "education.num", "workclass", "occupation", "hours.per.week")
out <- list()
for (s in c("income", "marital.status")) {
  sc <- list(demographic = q_demo, extended = q_ext, full = setdiff(names(adult), s))
  for (nm in names(sc)) {
    th <- rapid_threshold_select(original_data = adult, synthetic_data = syn1,
                                 quasi_identifiers = setdiff(sc[[nm]], s), sensitive_attribute = s,
                                 model_type = "rf", cat_eval_method = "RCS_marginal",
                                 n_permutations = 199, alpha = 0.05, seed = 1)
    out[[length(out) + 1]] <- data.frame(target = s, scenario = nm, tau_star = th$threshold_star)
  }
}
out <- do.call(rbind, out); print(out)
saveRDS(out, here("results", "F_threshold_select.rds"))
