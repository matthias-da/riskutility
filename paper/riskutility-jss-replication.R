## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  prompt = TRUE,
  comment = NA,
  collapse = FALSE,
  fig.width = 4.9,
  fig.height = 3.675,
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)
options(prompt = "R> ", continue = "+  ", width = 76, digits = 4,
        knitr.table.format = "latex")
library("riskutility")

## Helper: format a regression coefficient for inline reporting (no sci. notation)
fmt_coef <- function(x, term, which = c("orig", "synth")) {
  which <- match.arg(which)
  i <- grep(term, x$coefficients$term)[1]
  v <- x$coefficients[[paste0("estimate_", which)]][i]
  formatC(round(v), format = "d", big.mark = ",")
}

## Helper: render a Graphviz DOT specification to PDF for inclusion as a figure
## (requires DiagrammeR, DiagrammeRsvg, and rsvg).
render_dot <- function(dot) {
  f <- tempfile(fileext = ".pdf")
  svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
  rsvg::rsvg_pdf(charToRaw(svg), f)
  knitr::include_graphics(f)
}


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
  check.names = FALSE
)
knitr::kable(scope, caption = "Package scope at a glance. All families apply to both traditionally anonymized and synthetic data.\\label{tab:scope}")


## ----fig-decision, echo=FALSE, out.width="95%", fig.cap="Decision guide: choosing a riskutility function by assessment goal and disclosure threat.\\label{fig:decision}"----
render_dot('
digraph decision {
  graph [rankdir=LR, fontsize=11];
  node [shape=box, style="rounded,filled", fillcolor="#EAF2F8", fontname="Helvetica", fontsize=9];
  edge [fontname="Helvetica", fontsize=8];
  Q  [label="What do you want\nto assess?", fillcolor="#D4E6F1"];
  R  [label="Disclosure risk"];
  U  [label="Data utility"];
  M  [label="Compare multiple\nsynthesizers?", fillcolor="#D4E6F1"];
  Q -> R; Q -> U; Q -> M;
  R -> Tid [label="identity"];
  R -> Tat [label="attribute"];
  R -> Tme [label="membership"];
  R -> Tmo [label="memorization"];
  Tid [label="recordLinkage()\nkanonymity()\nindividual_risk()", fillcolor="#FDEDEC"];
  Tat [label="dcap(), tcap()\nweap(), disco(), rapid()", fillcolor="#FDEDEC"];
  Tme [label="domias(), nnaa()\nsingling_out(), linkability()", fillcolor="#FDEDEC"];
  Tmo [label="dcr(), nndr(), ims()", fillcolor="#FDEDEC"];
  U -> Ug [label="global"];
  U -> Ud [label="distributional"];
  U -> Us [label="structural"];
  U -> Up [label="predictive"];
  Ug [label="propscore()\nspecks()", fillcolor="#E8F8F5"];
  Ud [label="hellinger()\nenergy_distance(), mmd()", fillcolor="#E8F8F5"];
  Us [label="copula_fidelity()\ncontingency_fidelity()", fillcolor="#E8F8F5"];
  Up [label="tstr()\nregression_fidelity()", fillcolor="#E8F8F5"];
  M -> RU [label="yes"];
  RU [label="rumap()", fillcolor="#FCF3CF"];
}
')


## ----table2-threats, echo=FALSE-----------------------------------------------
threats <- data.frame(
  Threat = c("Identity", "Attribute", "Membership", "Memorization"),
  Definition = c(
    "Link a record to an individual",
    "Learn a sensitive value via linkage",
    "Determine whether an individual is in the data",
    "Generator reproduces training records"
  ),
  `Key measures` = c(
    "recordLinkage, kanonymity",
    "dcap, tcap, weap, disco, rapid",
    "mia_classifier, domias, nnaa",
    "ims, dcr, nndr"
  ),
  check.names = FALSE
)
knitr::kable(threats, caption = "Disclosure threat taxonomy.\\label{tab:threats}")


## ----table3-comparison, echo=FALSE--------------------------------------------
comp <- data.frame(
  Measure = c("k-Anonymity", "l-Diversity", "t-Closeness",
              "CAP family", "RAPID", "DCR / NNDR", "Record linkage",
              "WP29 criteria / MIA", "pMSE / SPECKS",
              "R-U map"),
  sdcMicro = c("kAnon()", "ldiversity()", "--",
               "--", "--", "dRisk()*", "recordLinkage()",
               "--", "dUtility()", "--"),
  synthpop = c("--", "--", "--",
               "disclosure()", "--", "--", "--",
               "--", "utility.gen()", "--"),
  Python = c("--", "--", "--",
             "SDMetrics", "--", "SDMetrics", "--",
             "Anonymeter", "--", "--"),
  riskutility = c("kanonymity()", "ldiversity()", "tcloseness()",
                  "dcap(), tcap()", "rapid()",
                  "dcr(), nndr()", "recordLinkage()",
                  "singling_out()",
                  "propscore(), specks()",
                  "rumap()"),
  check.names = FALSE
)
knitr::kable(comp,
             caption = "Risk and utility measures across the main R and Python tools. *sdcMicro's dRisk() is an interval/RMD risk for paired perturbed records, not a train-versus-holdout comparison.\\label{tab:software}")


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


## ----fig-integration, echo=FALSE, out.width="90%", fig.cap="Ecosystem integration: converters map synthpop, simPop, and sdcMicro objects (and plain data frames) into the paired-data container consumed by all measures.\\label{fig:integration}"----
render_dot('
digraph integration {
  graph [rankdir=LR, fontsize=10];
  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=9];
  edge [fontname="Helvetica", fontsize=8];
  synthpop [label="synthpop\n(synds)", fillcolor="#FDEBD0"];
  simPop   [label="simPop\n(simPopObj)", fillcolor="#FDEBD0"];
  sdcMicro [label="sdcMicro\n(sdcMicroObj)", fillcolor="#FDEBD0"];
  raw      [label="data frames\n(original, synthetic)", fillcolor="#FDEBD0"];
  sp [label="synth_pair", fillcolor="#D4E6F1"];
  ru [label="riskutility\nrisk + utility measures\ndisclosure_report(), rumap()", fillcolor="#D5F5E3"];
  synthpop -> sp [label="from_synthpop()"];
  simPop   -> sp [label="from_simPop()"];
  sdcMicro -> sp [label="from_sdcMicro()"];
  raw      -> sp [label="synth_pair()"];
  sp -> ru;
}
')


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
edu_levels <- c("Primary", "Secondary", "Tertiary")
regions <- paste0("R", 1:5)

age <- sample(25:75, n, replace = TRUE)
sex <- factor(sample(c("M", "F"), n, replace = TRUE))

# Education depends on age: younger cohorts are better educated.
education <- factor(ifelse(runif(n) < 0.5 * plogis(-0.05 * (age - 45)),
                           "Tertiary",
                           sample(edu_levels[1:2], n, replace = TRUE,
                                  prob = c(0.4, 0.6))),
                    levels = edu_levels)

# Region depends on education: R1 and R2 are urban and attract graduates.
region <- factor(ifelse(education == "Tertiary",
                        sample(regions, n, TRUE,
                               prob = c(.35, .30, .15, .10, .10)),
                        sample(regions, n, TRUE,
                               prob = c(.10, .15, .25, .25, .25))),
                 levels = regions)

# Income depends on age, education and sex.
edu_effect <- c(Primary = 0, Secondary = 0.35, Tertiary = 0.75)
income <- round(exp(log(20000) + 0.020 * age +
                      edu_effect[as.character(education)] -
                      0.15 * (sex == "F") + rnorm(n, 0, 0.18)))

original <- data.frame(age, sex, education, region, income)

# Synthetic: every variable is resampled independently from its own marginal.
# Marginal distributions are preserved; the joint structure is destroyed.
synthetic <- as.data.frame(lapply(original, function(v)
  sample(v, n, replace = TRUE)))

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
  Input = c(rep("Single dataset", 7), "Both datasets"),
  `Key output` = c("Min EC size", "Min distinct values per EC",
                    "Max EMD across ECs", "SUDA scores",
                    "Per-record frequency risk", "Estimated pop. uniques",
                    "Scenario-based risk", "Identifiability fraction"),
  `Threats` = c("Identity", "Attribute", "Attribute", "Identity",
                "Identity", "Identity", "Identity", "Identity"),
  check.names = FALSE
)
knitr::kable(privacy,
             caption = "Privacy models overview.\\label{tab:privacy}")


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
  `Low risk` = c("ratio <= 1.5", "ratio <= 1.5",
                 "<= 5% disclosive", "ratio <= 1.5"),
  check.names = FALSE
)
knitr::kable(cap,
             caption = "CAP family comparison with interpretation thresholds.\\label{tab:cap}")


## ----synthpop-comparison, warning=FALSE, eval=requireNamespace("synthpop", quietly = TRUE)----
keys_cmp <- c("sex", "region")
syn_sp   <- synthpop::syn(original, seed = 2024, print.flag = FALSE)
pair_sp  <- synth_pair(original, syn_sp$syn,
                       key_vars = keys_cmp, target_var = "education")
disc_sp  <- synthpop::disclosure(syn_sp, original, keys = keys_cmp,
                                 target = "education", print.flag = FALSE)

# from_synthpop() builds the paired container directly from a synds object,
# which is Contribution 5 of Section 7.1 in one line:
pair_conv <- from_synthpop(syn_sp, original,
                           key_vars = keys_cmp, target_var = "education")
stopifnot(identical(dcap(pair_conv)$cap, dcap(pair_sp)$cap))

cmp <- data.frame(
  Quantity = c("riskutility: dcap()$cap", "synthpop: DCAP", "synthpop: CAPd"),
  Value = round(c(100 * dcap(pair_sp)$cap,
                  as.numeric(disc_sp$allCAPs[["DCAP"]]),
                  as.numeric(disc_sp$allCAPs[["CAPd"]])), 4)
)


## ----synthpop-comparison-tab, echo=FALSE--------------------------------------
knitr::kable(cmp, row.names = FALSE, col.names = c("Quantity", "Value (%)"),
             caption = "Mean CAP from riskutility and synthpop on identical data. The first two rows are the same estimand and agree to all printed digits.\\label{tab:synthpop}")


## ----rapid-demo, warning=FALSE, fig.cap="RAPID threshold-sensitivity diagnostic for the running example.\\label{fig:rapid}"----
rapid_res <- rapid(pair, model_type = "lm")
summary(rapid_res)

# The raw score is not interpretable on its own: a fixed error band captures a
# non-zero share of records by chance. rapid_test() supplies the permutation
# null against which the score must be read.
rapid_null <- rapid_test(original, synthetic,
                         quasi_identifiers = key_vars,
                         sensitive_attribute = target_var,
                         model_type = "lm", n_permutations = 199, seed = 42)
print(rapid_null)

plot(rapid_res, which = 3)


## ----rapid-models, echo=FALSE-------------------------------------------------
models <- data.frame(
  Model = c("lm", "rf", "cart", "gbm", "logit"),
  Package = c("stats", "ranger", "rpart", "xgboost", "stats"),
  Numeric = c("Yes", "Yes", "Yes", "Yes", "No"),
  Categorical = c("No", "Yes", "Yes", "Yes", "Binary only"),
  Interactions = c("Manual", "Automatic", "Automatic", "Automatic", "Manual"),
  check.names = FALSE
)
knitr::kable(models, caption = "RAPID model backends.")


## ----distance-demo, warning=FALSE, fig.cap="Distance to closest record (DCR): synthetic-to-training versus synthetic-to-holdout distances.\\label{fig:dcr}"----
dcr_res <- dcr(pair, holdout_fraction = 0.5)
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
  `Low risk` = c("share < 0.55", "ratio >= 0.8", "< 0.01",
                 "ratio near 1", "< 0.05", "< 0.05",
                 "< 0.01", "delta_max well below 1"),
  check.names = FALSE
)
knitr::kable(dist,
             caption = "Distance-based and proximity risk measures.\\label{tab:distance}")


## ----recordlinkage-demo-------------------------------------------------------
# Namespace-qualified: sdcMicro (>= 5.8.2) exports an unrelated
# recordLinkage() with a different signature, which masks this one whenever
# sdcMicro is attached afterwards.
rl_res <- riskutility::recordLinkage(pair, method = "deterministic")
print(rl_res)


## ----recordlinkage-table, echo=FALSE------------------------------------------
rl <- data.frame(
  Method = c("Deterministic", "Probabilistic", "PRAM", "Predictive",
             "RF", "RBRL", "Mahalanobis", "Embedding"),
  Distance = c("Gower", "Fellegi-Sunter", "Transition prob.", "Propensity",
               "RF proximity", "Rank-based", "Mahalanobis", "Autoencoder"),
  `Mixed types` = c("Yes", "Yes", "Categorical", "Yes",
                     "Yes", "Yes", "Numeric", "Yes"),
  check.names = FALSE
)
knitr::kable(rl,
             caption = "Record linkage methods. All eight support all three matching modes (independent, bijective, optimal transport).\\label{tab:rl}")


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
             caption = "Membership inference and GDPR measures.\\label{tab:membership}")


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
comparison$Range <- c("(-1, 1), differential", "[0, 1]", "[0, 1]")
comparison$Safe <- sprintf("%.3f", comparison$Safe)
comparison$Risky <- sprintf("%.3f", comparison$Risky)
comparison <- comparison[, c("Metric", "Range", "Safe", "Risky")]
knitr::kable(comparison,
             caption = "Cross-family comparison: safe versus risky synthetic data.")


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
reg_res <- regression_fidelity(pair,
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
             caption = "Utility measures by use case.\\label{tab:utility}")


## ----fig-workflow, echo=FALSE, out.width="95%", fig.cap="The iterative risk-utility evaluation workflow.\\label{fig:workflow}"----
render_dot('
digraph workflow {
  graph [rankdir=LR, fontsize=10];
  node [shape=box, style="rounded,filled", fillcolor="#EAF2F8", fontname="Helvetica", fontsize=9];
  edge [fontname="Helvetica", fontsize=8];
  gen  [label="Generate /\nanonymize data"];
  risk [label="Assess risk\n(disclosure_report())"];
  util [label="Assess utility"];
  dec  [label="Acceptable\ntrade-off?", shape=diamond, fillcolor="#D4E6F1"];
  rel  [label="Release", fillcolor="#D5F5E3"];
  ref  [label="Refine synthesis /\nprotection parameters", fillcolor="#FCF3CF"];
  cmp  [label="Compare candidates\n(rumap())", fillcolor="#FCF3CF"];
  gen -> risk -> util -> dec;
  dec -> rel [label="yes"];
  dec -> ref [label="no"];
  ref -> cmp -> gen;
}
')


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

# The report's own print method is the practitioner-facing artefact:
print(rep_A)

verdicts <- data.frame(
  Method = c("A: Independent", "B: Bootstrap+noise", "C: Near-copy"),
  Overall = c(rep_A$overall_risk, rep_B$overall_risk, rep_C$overall_risk),
  Pass = c(rep_A$n_pass, rep_B$n_pass, rep_C$n_pass),
  Warn = c(rep_A$n_warn, rep_B$n_warn, rep_C$n_warn)
)


## ----case-report-tab, echo=FALSE----------------------------------------------
knitr::kable(verdicts,
             caption = "Quick risk screening across the three methods.\\label{tab:verdicts}")


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


## ----case-rumap-scatter, fig.height=4.2, fig.cap="Risk-Utility map of the three synthesizers; the Pareto-optimal front is highlighted.\\label{fig:rumap-scatter}"----
plot(ru, which = 1)  # R-U scatterplot with Pareto front


## ----case-rumap-heatmap, fig.height=3.6, fig.cap="Heatmap of normalized risk and utility measures across the three synthesizers.\\label{fig:rumap-heatmap}"----
plot(ru, which = 2)  # Heatmap of individual measures


## ----case-stability, warning=FALSE--------------------------------------------
make_case <- function(seed_data, seed_syn, N = 1000) {
  set.seed(seed_data)
  d <- data.frame(
    age_group = factor(sample(age_groups, N, replace = TRUE)),
    sex       = factor(sample(c("M", "F"), N, replace = TRUE)),
    education = factor(sample(edu_levels, N, replace = TRUE,
                              prob = c(0.25, 0.50, 0.25))),
    region    = factor(sample(paste0("R", 1:4), N, replace = TRUE)))
  d$income <- round(exp(10 + age_effect[as.character(d$age_group)] +
                          edu_effect[as.character(d$education)] +
                          rnorm(N, 0, 0.4)))
  set.seed(seed_syn)
  a <- as.data.frame(lapply(d, function(v) sample(v, N, replace = TRUE)))
  b <- d[sample(N, N, replace = TRUE), ]; rownames(b) <- NULL
  b$income <- round(b$income * exp(rnorm(N, 0, 0.15)))
  cc <- d; cc$income <- round(cc$income * exp(rnorm(N, 0, 0.03)))
  list(orig = d, A = a, B = b, C = cc)
}

seed_pairs <- list(c(123, 456), c(1, 2), c(11, 22), c(7, 8), c(99, 100),
                   c(2024, 2025), c(31, 41), c(5, 15), c(77, 88), c(300, 400))
front <- t(vapply(seed_pairs, function(sp) {
  cs <- make_case(sp[1], sp[2])
  set.seed(42)
  r <- rumap(cs$orig, list(A = cs$A, B = cs$B, C = cs$C),
             risk_measures = c("dcap", "tcap", "ims"),
             utility_measures = c("pmse", "wasserstein"),
             key_vars = qi, target_var = sens, seed = 42)
  r$pareto
}, logical(3)))
colnames(front) <- c("A", "B", "C")
round(100 * colMeans(front))


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
                   "). The attribution, frequency, and exact-match metrics scale near-linearly.\\label{tab:benchmark}"))


## ----session-info, echo=FALSE, results="asis"---------------------------------
si <- sessionInfo()
cat(sprintf(paste0("All results were obtained using %s on %s, with the ",
                   "\\pkg{riskutility} package version %s, ",
                   "\\pkg{data.table} %s, \\pkg{ggplot2} %s, ",
                   "\\pkg{ranger} %s and \\pkg{synthpop} %s. ",
                   "\\proglang{R} and all packages used are available from ",
                   "CRAN at https://CRAN.R-project.org/.\n"),
            si$R.version$version.string, si$running,
            packageVersion("riskutility"), packageVersion("data.table"),
            packageVersion("ggplot2"), packageVersion("ranger"),
            packageVersion("synthpop")))

