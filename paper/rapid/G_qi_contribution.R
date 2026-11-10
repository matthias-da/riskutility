## Experiment G: risk contribution of each quasi-identifier on UCI Adult (R1.9, R1.11, ED4).
## Drop-one from the full scenario, attacker refitted without QI j:
##   Delta_j = RAPID(Q) - RAPID(Q \ {j}), and the same for the excess over the fixed-model null.
## Plus one grouped drop {relationship, marital.status} (substitutes) for the income target.
## CART releases 1..M as in Experiment B / Table 3 (seeds 2025 + i).
## Usage: Rscript paper/rapid/G_qi_contribution.R [M]   (M = 1 was the pilot)
source("paper/rapid/common.R")
adult <- load_adult()
args <- commandArgs(trailingOnly = TRUE)
M <- if (length(args)) as.integer(args[1]) else 5

nullm <- function(orig, syn, qi, target, seed) {
  t <- rapid_test(original_data = orig, synthetic_data = syn, quasi_identifiers = qi,
                  sensitive_attribute = target, model_type = "rf", cat_tau = 0.3,
                  cat_eval_method = "RCS_marginal", n_permutations = 199, seed = seed)
  mean(t$null_distribution)
}

res <- list()
for (i in 1:M) {
  syn_i <- syn(adult, method = "cart", m = 1, seed = 2025 + i, print.flag = FALSE)$syn
  one <- function(qi, target, label) {
    cat("release", i, target, label, "\n")
    r <- rapid_paper(adult, syn_i, qi, target, seed = i)
    data.frame(rep = i, target = target, dropped = label, n_qi = length(qi), rapid = r$rapid,
               acc = r$acc, null = nullm(adult, syn_i, qi, target, seed = i))
  }
  for (s in c("income", "marital.status")) {
    q_full <- setdiff(names(adult), s)
    res[[length(res) + 1]] <- one(q_full, s, "(none)")
    for (j in q_full) res[[length(res) + 1]] <- one(setdiff(q_full, j), s, j)
    if (s == "income")
      res[[length(res) + 1]] <- one(setdiff(q_full, c("relationship", "marital.status")), s,
                                    "relationship+marital.status")
  }
}
res <- do.call(rbind, res)
res$excess <- res$rapid - res$null
## paired differences within release and target
key <- paste(res$rep, res$target)
ref <- res[res$dropped == "(none)", ]; rk <- paste(ref$rep, ref$target)
res$d_rapid  <- ref$rapid[match(key, rk)]  - res$rapid
res$d_excess <- ref$excess[match(key, rk)] - res$excess
res$d_acc    <- ref$acc[match(key, rk)]    - res$acc

out <- if (M == 1) "G_qi_contribution_pilot" else "G_qi_contribution"
saveRDS(res, here("results", paste0(out, ".rds")))
summ <- do.call(rbind, lapply(split(res, list(res$target, res$dropped), drop = TRUE), function(d)
  data.frame(target = d$target[1], dropped = d$dropped[1],
             d_rapid = mean(d$d_rapid), d_excess = mean(d$d_excess),
             d_excess_min = min(d$d_excess), d_excess_max = max(d$d_excess),
             d_acc = mean(d$d_acc))))
summ <- summ[order(summ$target, -summ$d_excess), ]
print(summ, digits = 3, row.names = FALSE)
write.csv(summ, here("results", paste0(out, "_summary.csv")), row.names = FALSE)
