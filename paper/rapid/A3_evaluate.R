## Experiment A, step 3: RAPID on training vs holdout records, utility, merge with Anonymeter
source("paper/rapid/common.R")
train <- read.csv(here("data", "A", "train.csv"), stringsAsFactors = TRUE)
hold  <- read.csv(here("data", "A", "holdout.csv"), stringsAsFactors = TRUE)
qi <- setdiff(names(train), "income")
align <- function(s) { for (v in names(train)) if (is.factor(train[[v]]))
  s[[v]] <- factor(as.character(s[[v]]), levels = levels(train[[v]])); s[names(train)] }
fo <- income ~ age + education.num + sex + hours.per.week

files <- list.files(here("data", "A"), pattern = "^syn_.*\\.csv$", full.names = TRUE)
res <- do.call(rbind, lapply(files, function(f) {
  rel <- sub("^syn_(.*)\\.csv$", "\\1", basename(f)); cat(rel, "\n")
  s <- align(read.csv(f, stringsAsFactors = TRUE))
  a <- rapid_paper(train, s, qi, "income"); h <- rapid_paper(hold, s, qi, "income")
  data.frame(release = rel, scheme = sub("_[0-9]+$", "", rel),
             rapid_train = a$rapid, rapid_hold = h$rapid, delta = a$rapid - h$rapid,
             tstr_acc = h$acc, ci_overlap = ci_overlap(train, s, fo))
}))
an <- read.csv(here("results", "A_anonymeter.csv"))
res <- merge(res, an, by = "release")
saveRDS(res, here("results", "A_schemes.rds"))
summ <- aggregate(cbind(rapid_train, rapid_hold, delta, anon_risk, tstr_acc, ci_overlap) ~ scheme,
                  data = res, FUN = mean)
sds  <- aggregate(cbind(rapid_train, rapid_hold, delta, anon_risk) ~ scheme, data = res, FUN = sd)
print(summ, digits = 3); print(sds, digits = 2)
write.csv(summ, here("results", "A_summary.csv"), row.names = FALSE)
