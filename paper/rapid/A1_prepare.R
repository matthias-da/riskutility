## Experiment A, step 1: train/holdout split of Adult; synthpop CART releases on train (parametric dropped: polyreg too slow with 41 countries)
source("paper/rapid/common.R")
adult <- load_adult()
set.seed(2026)
idx <- sample(nrow(adult), floor(nrow(adult) / 2))
dir.create(here("data", "A"), showWarnings = FALSE)
write.csv(adult[idx, ],  here("data", "A", "train.csv"),   row.names = FALSE)
write.csv(adult[-idx, ], here("data", "A", "holdout.csv"), row.names = FALSE)
train <- adult[idx, ]
for (m in "cart") for (i in 1:5) {
  s <- syn(train, method = m, m = 1, seed = 3000 + i, print.flag = FALSE)$syn
  write.csv(s, here("data", "A", sprintf("syn_%s_%d.csv", m, i)), row.names = FALSE)
}
## Reference anchor "no protection": the training records themselves released as-is
file.copy(here("data", "A", "train.csv"), here("data", "A", "syn_original_1.csv"), overwrite = TRUE)
