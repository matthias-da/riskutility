## Experiment C: simulation 1 extended by attacker knowledge and a utility metric (R1.9, R1.10)
## DGP and settings as in Section 6 / Appendix (n = 1000, CART, RF attacker, tau = 0.3).
source("paper/rapid/common.R")
source("paper/rapid/simulate_microdata.R")

kappas <- c(0, 2, 5, 10, 20, 50, 100); n_sim <- 10
qi_all  <- c("gender", "age", "education", "income", "health_score")
qi_demo <- c("gender", "age")
fo_util <- I(disease_status != "healthy") ~ gender + age + education + log(income) + health_score

res <- list()
for (k in kappas) for (sim in 1:n_sim) {
  cat("kappa", k, "sim", sim, "\n")
  o <- simulate_microdata(n = 1000, dep = k, seed = 2025 + sim)
  s <- syn(o, method = "cart", print.flag = FALSE, seed = 2025 + sim + 1000)$syn
  a <- rapid_paper(o, s, qi_all,  "disease_status", seed = sim)
  d <- rapid_paper(o, s, qi_demo, "disease_status", seed = sim)
  u <- suppressWarnings(ci_overlap(o, s, fo_util))
  res[[length(res) + 1]] <- data.frame(kappa = k, sim = sim, rapid_all = a$rapid,
                                       rapid_demo = d$rapid, acc_all = a$acc, ci_overlap = u)
}
res <- do.call(rbind, res)
saveRDS(res, here("results", "C_simulation.rds"))
summ <- aggregate(. ~ kappa, data = res[, -2], FUN = function(x) c(m = mean(x), sd = sd(x)))
print(summ, digits = 3)
