## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  prompt = TRUE,
  comment = NA,
  collapse = FALSE,
  fig.width = 7,
  fig.height = 5,
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)
options(prompt = "R> ", continue = "+  ", width = 76, digits = 4,
        knitr.table.format = "latex")
library("riskutility")


## ----intro-hook, eval=FALSE---------------------------------------------------
# library(riskutility)
# pair <- synth_pair(original, synthetic,
#                    key_vars = c("age", "sex", "region"),
#                    target_var = "income")
# report <- disclosure_report(pair)
# print(report)


## ----table1-scope, echo=FALSE-------------------------------------------------
scope <- data.frame(
  Category = c("Privacy models", "Attribution (CAP)",
                "ML-based (RAPID)", "Distance-based", "Record linkage",
                "Membership inference", "Utility measures", "Frameworks"),
  Functions = c(8, 4, 4, 6, "1 (8 methods)", 6, 15, 3),
  Paradigm = c("Equivalence class", "Matching", "Prediction",
                "Nearest neighbor", "Linkage", "Attack simulation",
                "Various", "Composite"),
  `Applicable to` = c(rep("Both", 8)),
  check.names = FALSE
)
knitr::kable(scope, caption = "Package scope at a glance.
             Both = applicable to traditionally anonymized and synthetic data.")


## ----table2-threats, echo=FALSE-----------------------------------------------
threats <- data.frame(
  Threat = c("Identity", "Attribute", "Membership", "Memorization"),
  Definition = c(
    "Attacker links record to individual",
    "Attacker learns sensitive value via linkage",
    "Attacker determines if individual is in dataset",
    "Generator reproduces training records"
  ),
  `Key measures` = c(
    "recordLinkage, kanonymity, individual_risk",
    "dcap, tcap, weap, disco, rapid",
    "mia_classifier, domias, nnaa, singling_out",
    "ims, dcr, nndr"
  ),
  check.names = FALSE
)
knitr::kable(threats, caption = "Disclosure threat taxonomy.")


## ----table3-comparison, echo=FALSE--------------------------------------------
comp <- data.frame(
  Measure = c("k-Anonymity", "l-Diversity", "t-Closeness",
              "DCAP/TCAP", "RAPID", "DCR/NNDR", "Record linkage",
              "MIA / GDPR criteria", "pMSE / SPECKS",
              "R-U map"),
  sdcMicro = c("freq()", "ldiversity()", "--",
               "--", "--", "dRisk()", "--",
               "--", "--", "--"),
  synthpop = c("--", "--", "--",
               "disclosure()", "--", "--", "--",
               "--", "utility.gen()", "--"),
  `SDMetrics (Py)` = c("--", "--", "--",
                        "--", "--", "Yes", "--",
                        "--", "--", "--"),
  `Anonymeter (Py)` = c("--", "--", "--",
                         "--", "--", "--", "--",
                         "Yes", "--", "--"),
  riskutility = c("kanonymity()", "ldiversity()", "tcloseness()",
                   "dcap(), tcap()", "rapid()", "dcr(), nndr()",
                   "recordLinkage()",
                   "mia_classifier(), singling_out()", "propscore(), specks()",
                   "rumap()"),
  check.names = FALSE
)
knitr::kable(comp,
             caption = "Risk/utility measures across R and Python packages.")


## ----synthpop-comparison, warning=FALSE---------------------------------------
cmp <- tryCatch({
  if (!requireNamespace("synthpop", quietly = TRUE))
    stop("synthpop not available")
  keys_cmp <- c("age", "sex", "region")
  syn_sp <- synthpop::syn(original, seed = 2024, print.flag = FALSE)
  pair_sp <- synth_pair(original, syn_sp$syn,
                        key_vars = keys_cmp, target_var = "education")
  cap_ru  <- 100 * dcap(pair_sp)$cap
  disc_sp <- synthpop::disclosure(syn_sp, original, keys = keys_cmp,
                                  target = "education", print.flag = FALSE)
  cap_sp  <- as.numeric(disc_sp$allCAPs[["CAPd"]])
  data.frame(Tool = c("riskutility::dcap()", "synthpop::disclosure()"),
             `Mean CAP (%)` = round(c(cap_ru, cap_sp), 2), check.names = FALSE)
}, error = function(e)
  data.frame(Note = paste("Comparison skipped:", conditionMessage(e))))
knitr::kable(cmp, row.names = FALSE,
             caption = "Mean CAP from riskutility and synthpop on identical data.")


## ----synth-pair-demo, eval=FALSE----------------------------------------------
# pair <- synth_pair(original, synthetic,
#                    key_vars = c("age", "gender", "region"),
#                    target_var = "income",
#                    holdout = holdout_data)


## ----synth-pair-dispatch, eval=FALSE------------------------------------------
# # All functions accept synth_pair --- no parameter repetition:
# dcap(pair)                    # Attribution risk
# rapid(pair, model_type = "rf") # ML-based risk
# propscore(pair)                # Propensity score utility
# disclosure_report(pair)        # Full risk report
# rumap(pair)                    # Risk-Utility map


## ----s3-pattern, eval=FALSE---------------------------------------------------
# # 1. Two equivalent calling conventions:
# result <- dcap(pair)                                      # synth_pair method
# result <- dcap(X, Y, key_vars = ..., target_var = ...)    # default method
# 
# # 2. Inspection:
# print(result)                # One-screen summary with key statistic
# s <- summary(result)         # Detailed statistics (returns summary.dcap)
# print(s)                     # Formatted multi-line output
# 
# # 3. Visualization:
# plot(result, which = 1)      # Plot type 1
# plot(result, which = 1:2)    # Multiple plot types


## ----integration, eval=FALSE--------------------------------------------------
# # From synthpop: pass synds object + original data
# pair <- from_synthpop(synds_object, original_data,
#                       key_vars = c("age", "sex"),
#                       target_var = "income")
# 
# # From simPop: original data extracted automatically from simPopObj
# pair <- from_simPop(simPopObj,
#                     key_vars = c("age", "sex"),
#                     target_var = "income")
# 
# # From sdcMicro: variable roles extracted from sdcMicroObj
# pair <- from_sdcMicro(sdcMicroObj)


## ----running-data-------------------------------------------------------------
set.seed(42)
n <- 500
original <- data.frame(
  age = sample(18:85, n, replace = TRUE),
  sex = factor(sample(c("M", "F"), n, replace = TRUE)),
  education = factor(sample(c("Primary", "Secondary", "Tertiary"), n,
                            replace = TRUE, prob = c(0.3, 0.5, 0.2))),
  region = factor(sample(paste0("R", 1:5), n, replace = TRUE)),
  income = round(rlnorm(n, log(40000), 0.5))
)

# Synthetic: independent draws (low risk expected)
synthetic <- data.frame(
  age = sample(18:85, n, replace = TRUE),
  sex = factor(sample(c("M", "F"), n, replace = TRUE)),
  education = factor(sample(c("Primary", "Secondary", "Tertiary"), n,
                            replace = TRUE, prob = c(0.3, 0.5, 0.2))),
  region = factor(sample(paste0("R", 1:5), n, replace = TRUE)),
  income = round(rlnorm(n, log(40000), 0.5))
)

key_vars <- c("age", "sex", "education", "region")
target_var <- "income"

pair <- synth_pair(original, synthetic,
                   key_vars = key_vars, target_var = target_var)

# Train/holdout split for distance-based metrics
set.seed(123)
train_idx <- sample(n, size = floor(0.7 * n))
train_data <- original[train_idx, ]
holdout_data <- original[-train_idx, ]


## ----privacy-models-----------------------------------------------------------
# k-Anonymity: minimum equivalence class size
k_res <- riskutility::kanonymity(synthetic, key_vars = key_vars)
k_res

# l-Diversity: sensitive attribute diversity per EC
l_res <- riskutility::ldiversity(synthetic, key_vars = key_vars,
                                 sensitive_var = target_var)
print(l_res)

# t-Closeness: EMD between EC and overall distribution
t_res <- riskutility::tcloseness(synthetic, key_vars = key_vars,
                                 sensitive_var = target_var)
t_res


## ----table5-privacy, echo=FALSE-----------------------------------------------
privacy <- data.frame(
  Function = c("kanonymity()", "ldiversity()", "tcloseness()",
               "suda()", "individual_risk()", "population_uniqueness()",
               "attacker_risk()", "epsilon_identifiability()"),
  Input = rep("Single dataset", 8),
  `Key output` = c("Min EC size", "Min distinct values per EC",
                    "Max EMD across ECs", "SUDA scores",
                    "Per-record frequency risk", "Estimated pop. uniques",
                    "Scenario-based risk", "Identifiability fraction"),
  `Threats` = c("Identity", "Attribute", "Attribute", "Identity",
                "Identity", "Identity", "Identity", "Identity"),
  check.names = FALSE
)
knitr::kable(privacy,
             caption = "Privacy models overview.")


## ----cap-demo, fig.cap="Per-record TCAP attribution-risk scores for the running example.\\label{fig:tcap}"----
# TCAP: per-record risk (most informative member of CAP family)
tcap_res <- tcap(pair)
summary(tcap_res)
plot(tcap_res)


## ----cap-table, echo=FALSE----------------------------------------------------
cap <- data.frame(
  Metric = c("DCAP", "TCAP", "WEAP", "DiSCO"),
  `Requires original?` = c("Yes", "Yes", "No", "Yes"),
  `Per-record?` = c("No", "Yes", "Yes", "Yes"),
  `Measures` = c("Mean attribution probability",
                 "Individual attribution risk",
                 "Within-EC homogeneity",
                 "Correct + confident attribution"),
  `Low risk` = c("ratio < 1.5", "< 0.1 per record",
                 "< 0.1", "< 5%"),
  check.names = FALSE
)
knitr::kable(cap,
             caption = "CAP family comparison with interpretation thresholds.")


## ----rapid-demo, warning=FALSE, fig.cap="RAPID risk histogram (left) and threshold-sensitivity diagnostic (right).\\label{fig:rapid}"----
rapid_res <- rapid(pair, model_type = "lm")
summary(rapid_res)
plot(rapid_res, which = c(1, 3))


## ----rapid-models, echo=FALSE-------------------------------------------------
models <- data.frame(
  Model = c("lm", "rf", "cart", "gbm", "logit"),
  Package = c("stats", "ranger", "rpart", "xgboost", "stats"),
  Numeric = c("Yes", "Yes", "Yes", "Yes", "No"),
  Categorical = c("No", "Yes", "Yes", "Yes", "Yes"),
  Interactions = c("Manual", "Automatic", "Automatic", "Automatic", "Manual"),
  check.names = FALSE
)
knitr::kable(models, caption = "RAPID model backends.")


## ----distance-demo, warning=FALSE, fig.cap="Distance to closest record (DCR): synthetic-to-training versus synthetic-to-holdout distances.\\label{fig:dcr}"----
dcr_res <- dcr(pair, holdout_fraction = 0.2)
summary(dcr_res)
plot(dcr_res, which = 1)


## ----distance-table, echo=FALSE-----------------------------------------------
dist <- data.frame(
  Metric = c("DCR", "NNDR", "IMS", "RF proximity", "dRisk", "Hitting rate",
             "Epsilon ID", "Delta-presence"),
  Holdout = c("Yes", "Yes", "No", "Yes", "No", "No", "No", "No"),
  Detects = c("Memorization", "Memorization", "Exact copies",
              "Memorization (non-linear)", "Close records", "Close records",
              "Identifiability", "Membership bounds"),
  `Low risk` = c("share < 0.55", "share < 0.55", "< 0.01",
                 "ratio near 1", "< 0.05", "< 0.05",
                 "< 0.01", "> 0.5"),
  check.names = FALSE
)
knitr::kable(dist,
             caption = "Distance-based and proximity risk measures.")


## ----recordlinkage-demo-------------------------------------------------------
rl_res <- recordLinkage(pair, method = "deterministic")
print(rl_res)


## ----recordlinkage-table, echo=FALSE------------------------------------------
rl <- data.frame(
  Method = c("Deterministic", "Probabilistic", "PRAM", "Predictive",
             "RF", "RBRL", "Mahalanobis", "Embedding"),
  Distance = c("Gower", "Fellegi-Sunter", "Transition prob.", "Propensity",
               "RF proximity", "Rank-based", "Mahalanobis", "Autoencoder"),
  `Mixed types` = c("Yes", "Yes", "Categorical", "Yes",
                     "Yes", "Yes", "Numeric", "Yes"),
  Matching = c("All 3", "All 3", "All 3", "All 3",
               "All 3", "Independent", "All 3", "All 3"),
  check.names = FALSE
)
knitr::kable(rl,
             caption = "Record linkage methods. All 3 = independent, bijective, OT.")


## ----nnaa-demo----------------------------------------------------------------
nnaa_res <- nnaa(train_data, synthetic, holdout = holdout_data,
                 method = "gower", seed = 42)
print(nnaa_res)


## ----membership-demo----------------------------------------------------------
so_res <- singling_out(original, synthetic,
                       n_attacks = 500, n_cols = 3,
                       mode = "multivariate", seed = 42)
print(so_res)

link_res <- linkability(original, synthetic,
                        n_attacks = 500, n_neighbors = 1, seed = 42)
print(link_res)


## ----membership-table, echo=FALSE---------------------------------------------
mia <- data.frame(
  Metric = c("MIA classifier", "DOMIAS", "NNAA",
             "Singling out", "Linkability", "delta-Presence"),
  `Attack type` = c("Shadow model", "Density overfitting",
                     "Nearest neighbor", "Predicate-based",
                     "Record linkage", "Membership bounds"),
  Holdout = c("Yes", "Yes", "Yes", "Yes", "Yes", "No"),
  `GDPR criterion` = c("--", "--", "--",
                        "Art. 29 WP", "Art. 29 WP", "--"),
  `Low risk` = c("< 0.55", "< 0.6", "< 0.05",
                 "< 0.1", "< 0.1", "> 0.5"),
  check.names = FALSE
)
knitr::kable(mia,
             caption = "Membership inference and GDPR measures.")


## ----rosetta------------------------------------------------------------------
# Near-copy: original + small noise (high risk expected)
set.seed(99)
near_copy <- original
near_copy$age <- near_copy$age + sample(-1:1, n, replace = TRUE)
near_copy$income <- near_copy$income + round(rnorm(n, 0, 500))
pair_risky <- synth_pair(original, near_copy,
                         key_vars = key_vars, target_var = target_var)

# Compare key metrics across the two datasets
comparison <- data.frame(
  Metric = c("DCAP", "RAPID (lm)", "IMS"),
  Safe = c(
    dcap(pair)$dcap,
    rapid(pair, model_type = "lm", verbose = FALSE)$rapid,
    ims(pair)$ims
  ),
  Risky = c(
    dcap(pair_risky)$dcap,
    rapid(pair_risky, model_type = "lm", verbose = FALSE)$rapid,
    ims(pair_risky)$ims
  )
)
comparison$Safe <- round(comparison$Safe, 4)
comparison$Risky <- round(comparison$Risky, 4)
knitr::kable(comparison,
             caption = "Cross-family comparison: safe vs. risky synthetic data.")


## ----utility-quick, warning=FALSE---------------------------------------------
prop_res <- propscore(pair)
summary(prop_res)


## ----utility-univariate-------------------------------------------------------
# Hellinger distance for categorical variables
h_res <- hellinger(original, synthetic, vars = c("sex", "education"))
print(h_res)

# CI proximity: confidence interval overlap for means
cip_res <- ci_proximity(original, synthetic, vars = c("age", "income"))
print(cip_res)


## ----utility-structural-------------------------------------------------------
e_res <- energy_distance(original[, c("age", "income")],
                         synthetic[, c("age", "income")],
                         seed = 42)
print(e_res)


## ----mmd-demo-----------------------------------------------------------------
mmd_res <- mmd(original[, c("age", "income")],
               synthetic[, c("age", "income")],
               kernel = "gaussian", method = "rff",
               n_features = 500, seed = 42)
print(mmd_res)


## ----fidelity-demo------------------------------------------------------------
cop_res <- copula_fidelity(original, synthetic, vars = c("age", "income"))
print(cop_res)

ctf_res <- contingency_fidelity(original, synthetic,
                                vars = c("sex", "education", "region"))
print(ctf_res)


## ----tstr-demo, warning=FALSE, eval=requireNamespace("ranger", quietly=TRUE)----
set.seed(42)
tstr_res <- tstr(pair, target_var = "income", model = "rf",
                 test_fraction = 0.3, seed = 42)
print(tstr_res)


## ----regression-demo, fig.cap="Regression-coefficient fidelity: original versus synthetic estimates with confidence intervals.\\label{fig:reg}"----
reg_res <- regression_fidelity(original, synthetic,
                               formula = income ~ age + sex + education)
summary(reg_res)
plot(reg_res, which = 1)


## ----tail-demo----------------------------------------------------------------
tail_res <- tail_fidelity(original, synthetic, vars = c("age", "income"),
                          percentile = 95, tails = "both")
print(tail_res)


## ----subgroup-demo------------------------------------------------------------
su_res <- subgroup_utility(original, synthetic, group_var = "region",
                           utility_fun = energy_distance,
                           threshold = 0.5, seed = 42)
print(su_res)


## ----table7-utility, echo=FALSE-----------------------------------------------
util <- data.frame(
  `Use case` = c(rep("Quick assessment", 2),
                 rep("Univariate", 3),
                 rep("Multivariate", 4),
                 rep("Predictive", 3),
                 "Subgroup"),
  Function = c("propscore()", "specks()",
               "compare_wasserstein()", "hellinger()", "ci_proximity()",
               "energy_distance()", "mmd()",
               "copula_fidelity()", "contingency_fidelity()",
               "tstr()", "regression_fidelity()",
               "compare_feature_importance()",
               "subgroup_utility()"),
  `Data type` = c("Mixed", "Mixed",
                   "Numeric", "Categorical", "Numeric",
                   "Numeric", "Numeric",
                   "Numeric", "Categorical",
                   "Mixed", "Mixed", "Mixed",
                   "Mixed"),
  Interpretation = c("< 0.1: good", "< 0.05: good",
                     "Lower = better", "< 0.1: good", "> 0.8: good",
                     "Lower = better", "Lower = better",
                     "< 0.1: good", "< 0.05: good",
                     "ratio near 1: good", "overlap > 0.8: good",
                     "High corr: good",
                     "min > 0.5: good"),
  check.names = FALSE
)
knitr::kable(util,
             caption = "Utility measures by use case.")


## ----case-data----------------------------------------------------------------
set.seed(123)
N <- 1000
edu_levels <- c("Primary", "Secondary", "Tertiary")
age_groups <- c("20-29", "30-39", "40-49", "50-59", "60-69")
orig <- data.frame(
  age_group = factor(sample(age_groups, N, replace = TRUE)),
  sex = factor(sample(c("M", "F"), N, replace = TRUE)),
  education = factor(sample(edu_levels, N, replace = TRUE,
                            prob = c(0.25, 0.50, 0.25))),
  region = factor(sample(paste0("R", 1:4), N, replace = TRUE))
)
edu_effect <- c(Primary = 0, Secondary = 0.3, Tertiary = 0.7)
age_effect <- c("20-29" = 0, "30-39" = 0.15, "40-49" = 0.3,
                "50-59" = 0.4, "60-69" = 0.35)
orig$income <- round(exp(
  10 + age_effect[as.character(orig$age_group)] +
    edu_effect[as.character(orig$education)] + rnorm(N, 0, 0.4)
))

qi <- c("age_group", "sex", "education", "region")
sens <- "income"


## ----case-synthesis-----------------------------------------------------------
set.seed(456)

# Method A: Independent marginals (safest, but destroys correlations)
synA <- data.frame(
  age_group = factor(sample(age_groups, N, replace = TRUE)),
  sex = factor(sample(c("M", "F"), N, replace = TRUE)),
  education = factor(sample(edu_levels, N, replace = TRUE,
                            prob = c(0.25, 0.50, 0.25))),
  region = factor(sample(paste0("R", 1:4), N, replace = TRUE)),
  income = sample(orig$income, N, replace = TRUE)
)

# Method B: Category-preserving bootstrap with income noise
idx_B <- sample(N, N, replace = TRUE)
synB <- orig[idx_B, ]
rownames(synB) <- NULL
synB$income <- round(synB$income * exp(rnorm(N, 0, 0.15)))
swap_idx <- sample(N, round(0.2 * N))
synB$age_group[swap_idx] <- factor(sample(age_groups,
                                          length(swap_idx), replace = TRUE))

# Method C: Near-copy with minimal perturbation (risky)
synC <- orig
synC$income <- round(synC$income * exp(rnorm(N, 0, 0.03)))


## ----case-report, warning=FALSE-----------------------------------------------
pair_A <- synth_pair(orig, synA, key_vars = qi, target_var = sens)
pair_B <- synth_pair(orig, synB, key_vars = qi, target_var = sens)
pair_C <- synth_pair(orig, synC, key_vars = qi, target_var = sens)

rep_A <- disclosure_report(pair_A, compute = c("attribution", "privacy"),
                           seed = 42, verbose = FALSE)
rep_B <- disclosure_report(pair_B, compute = c("attribution", "privacy"),
                           seed = 42, verbose = FALSE)
rep_C <- disclosure_report(pair_C, compute = c("attribution", "privacy"),
                           seed = 42, verbose = FALSE)

verdicts <- data.frame(
  Method = c("A: Independent", "B: Bootstrap+noise", "C: Near-copy"),
  Overall = c(rep_A$overall_risk, rep_B$overall_risk, rep_C$overall_risk),
  Pass = c(rep_A$n_pass, rep_B$n_pass, rep_C$n_pass),
  Warn = c(rep_A$n_warn, rep_B$n_warn, rep_C$n_warn)
)
knitr::kable(verdicts, caption = "Quick risk screening across three methods.")


## ----case-rumap, warning=FALSE------------------------------------------------
set.seed(42)
ru <- rumap(orig,
            list("A: Independent" = synA,
                 "B: Bootstrap+noise" = synB,
                 "C: Near-copy" = synC),
            risk_measures = c("dcap", "tcap", "ims"),
            utility_measures = c("pmse", "wasserstein"),
            key_vars = qi, target_var = sens,
            seed = 42)
print(ru)


## ----case-rumap-scatter, fig.width=8, fig.height=6, fig.cap="Risk-Utility map of the three synthesizers; the Pareto-optimal front is highlighted.\\label{fig:rumap-scatter}"----
plot(ru, which = 1)  # R-U scatterplot with Pareto front


## ----case-rumap-heatmap, fig.width=8, fig.height=5, fig.cap="Heatmap of normalized risk and utility measures across the three synthesizers.\\label{fig:rumap-heatmap}"----
plot(ru, which = 2)  # Heatmap of individual measures


## ----scalability-benchmark, echo=FALSE----------------------------------------
bench_fun <- function(f, p) {
  round(median(replicate(3, system.time(f(p))[["elapsed"]])), 3)
}
mk_pair <- function(n) {
  set.seed(1)
  d1 <- data.frame(
    age = sample(18:85, n, TRUE),
    sex = factor(sample(c("M", "F"), n, TRUE)),
    education = factor(sample(c("Primary", "Secondary", "Tertiary"), n, TRUE)),
    region = factor(sample(paste0("R", 1:5), n, TRUE)),
    income = round(rlnorm(n, log(40000), 0.5)))
  d2 <- d1
  d2$age <- sample(18:85, n, TRUE)
  d2$income <- round(rlnorm(n, log(40000), 0.5))
  synth_pair(d1, d2, key_vars = c("age", "sex", "education", "region"),
             target_var = "income")
}
fns <- list("dcap()" = function(p) dcap(p),
            "tcap()" = function(p) tcap(p),
            "kanonymity()" = function(p) kanonymity(p),
            "ims()" = function(p) ims(p))
sizes <- c(1000, 5000, 20000)
timing <- sapply(sizes, function(n) {
  p <- mk_pair(n)
  vapply(fns, bench_fun, numeric(1), p = p)
})
timing <- data.frame(Metric = names(fns), timing, check.names = FALSE)
colnames(timing)[-1] <- paste0("n = ", format(sizes, big.mark = ","), " (s)")
knitr::kable(timing, row.names = FALSE,
  caption = paste0("Median elapsed time over three runs, measured on ",
                   R.version$platform, " (R ", getRversion(),
                   "). The attribution, frequency, and exact-match metrics scale near-linearly."))


## ----session-info-------------------------------------------------------------
sessionInfo()

