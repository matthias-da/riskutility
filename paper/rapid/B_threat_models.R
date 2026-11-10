## Experiment B: RAPID under alternative threat models on UCI Adult (R1.3, R1.8, ED4, R2.5)
## 3 attacker-knowledge scenarios x 2 sensitive attributes, CART releases as in Section 5.
source("paper/rapid/common.R")
adult <- load_adult()

M <- 5
syn_list <- lapply(1:M, function(i)
  syn(adult, method = "cart", m = 1, seed = 2025 + i, print.flag = FALSE)$syn)

q_demo <- c("age", "sex", "race", "native.country")
q_ext  <- c(q_demo, "education.num", "workclass", "occupation", "hours.per.week")
scen <- list(demographic = q_demo, extended = q_ext)
targets <- c("income", "marital.status")

res <- list(); thr <- list()
for (s in targets) {
  sc <- c(scen, list(full = setdiff(names(adult), s)))
  for (nm in names(sc)) {
    qi <- setdiff(sc[[nm]], s)
    for (i in 1:M) {
      cat(s, nm, i, "\n")
      res[[length(res) + 1]] <- cbind(target = s, scenario = nm, rep = i, n_qi = length(qi),
                                      rapid_paper(adult, syn_list[[i]], qi, s, seed = i))
    }
  }
}
## Note: a first run with rapid_threshold_select (B = 199) gave tau*_alpha = 0.95 (grid maximum)
## in all six cells: at n ~ 30k the permutation calibration is non-binding.
res <- do.call(rbind, res)
saveRDS(res, here("results", "B_threat_models.rds"))

summ <- aggregate(cbind(rapid, lo, hi, rapid_t50, rapid_t70, acc) ~ target + scenario + n_qi,
                  data = res, FUN = mean)
summ$baseline_acc <- sapply(summ$target, function(s) max(prop.table(table(adult[[s]]))))
print(summ[order(summ$target, summ$n_qi), ], digits = 3)
write.csv(summ, here("results", "B_summary.csv"), row.names = FALSE)
