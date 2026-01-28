## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----install, eval=FALSE------------------------------------------------------
# # Install from CRAN (when available)
# install.packages("riskutility")
# 
# # Or install development version from GitHub
# # devtools::install_github("mtempl/riskutility")

## ----load---------------------------------------------------------------------
library(riskutility)
library(ggplot2)

## ----example-data-------------------------------------------------------------
set.seed(42)
n <- 500

# Original data
original <- data.frame(
  age = as.integer(sample(18:80, n, replace = TRUE)),
  gender = factor(sample(c("Male", "Female"), n, replace = TRUE)),
  education = factor(sample(c("Primary", "Secondary", "Tertiary"), n,
                           replace = TRUE, prob = c(0.2, 0.5, 0.3))),
  region = factor(sample(paste0("Region_", 1:5), n, replace = TRUE)),
  income = round(rgamma(n, shape = 2, scale = 15000), -2),
  health_status = factor(sample(c("Good", "Fair", "Poor"), n,
                                replace = TRUE, prob = c(0.6, 0.3, 0.1)))
)

# Synthetic data (with some noise and distributional shifts)
synthetic <- data.frame(
  age = as.integer(pmin(80, pmax(18, original$age + sample(-3:3, n, replace = TRUE)))),
  gender = factor(sample(c("Male", "Female"), n, replace = TRUE)),
  education = factor(sample(c("Primary", "Secondary", "Tertiary"), n,
                           replace = TRUE, prob = c(0.22, 0.48, 0.30))),
  region = factor(sample(paste0("Region_", 1:5), n, replace = TRUE)),
  income = round(rgamma(n, shape = 2.1, scale = 14500), -2),
  health_status = factor(sample(c("Good", "Fair", "Poor"), n,
                                replace = TRUE, prob = c(0.58, 0.32, 0.10)))
)

head(original)

## ----dcap---------------------------------------------------------------------
# Define quasi-identifiers (key variables) and sensitive attribute (target)
key_vars <- c("age", "gender", "education", "region")
target_var <- "health_status"

# Compute DCAP
dcap_result <- dcap(
  X = original,
  Y = synthetic,
  key_vars = key_vars,
  target_var = target_var
)

print(dcap_result)

## ----dcap-summary-------------------------------------------------------------
summary(dcap_result)

## ----tcap---------------------------------------------------------------------
tcap_result <- tcap(
  X = original,
  Y = synthetic,
  key_vars = key_vars,
  target_var = target_var
)

print(tcap_result)

## ----tcap-plot, fig.cap="Distribution of per-record TCAP scores"--------------
plot(tcap_result)

## ----weap---------------------------------------------------------------------
# WEAP only uses synthetic data (Y) - it measures disclosure risk
# within the synthetic dataset itself
weap_result <- weap(
  Y = synthetic,
  key_vars = key_vars,
  target_var = target_var
)

print(weap_result)
summary(weap_result)

## ----disco--------------------------------------------------------------------
disco_result <- disco(
  X = original,
  Y = synthetic,
  key_vars = key_vars,
  target_var = target_var
)

print(disco_result)

## ----rapid, eval=requireNamespace("ranger", quietly = TRUE)-------------------
rapid_result <- rapid(
  X = original,
  Y = synthetic,
  key_vars = key_vars,
  target_var = target_var,
  model_type = "rf",
  cat_tau = 0.3,
  cat_eval_method = "RCS_marginal"
)

print(rapid_result)
summary(rapid_result)

## ----dcr----------------------------------------------------------------------
# Split original data into training and holdout
set.seed(123)
n_orig <- nrow(original)
train_idx <- sample(n_orig, size = floor(0.7 * n_orig))

train_data <- original[train_idx, ]
holdout_data <- original[-train_idx, ]

# For this example, synthetic was generated from full original
# In practice, synthetic should be generated only from train_data

dcr_result <- dcr(
  X = train_data,
  Y = synthetic,
  holdout = holdout_data,
  progress = FALSE
)

print(dcr_result)

## ----dcr-plot, fig.cap="DCR distribution comparison"--------------------------
plot(dcr_result)

## ----nndr---------------------------------------------------------------------
nndr_result <- nndr(
  X = train_data,
  Y = synthetic,
  holdout = holdout_data,
  progress = FALSE
)

print(nndr_result)
summary(nndr_result)

## ----ims----------------------------------------------------------------------
ims_result <- ims(
  X = train_data,
  Y = synthetic
)

print(ims_result)

## ----propscore----------------------------------------------------------------
ps_result <- propscore(
  X = original,
  Y = synthetic,
  form = ~ age + income + education + gender
)

print(ps_result)

## ----propscore-summary--------------------------------------------------------
summary(ps_result)

## ----compare-dist-------------------------------------------------------------
dist_result <- compare_distributions_cont(
  X = original,
  Y = synthetic,
  variables = c("age", "income")
)

print(dist_result)

## ----compare-dist-plot, fig.cap="Distribution comparison for continuous variables"----
plot(dist_result)

## ----compare-cat--------------------------------------------------------------
chisq_result <- compare_chisq_gof(
  X = original,
  Y = synthetic,
  cat_vars = "gender"
)

print(chisq_result)

## ----summary-stats------------------------------------------------------------
stats_result <- compare_multivariate_summary_statistics(
  X = original,
  Y = synthetic,
  cont_vars = c("age", "income"),
  cat_vars = c("gender", "education")
)

print(stats_result$continuous)

## ----disclosure-report--------------------------------------------------------
report <- disclosure_report(
  X = original,
  Y = synthetic,
  key_vars = key_vars,
  target_var = target_var
)

print(report)

## ----report-summary-----------------------------------------------------------
summary(report)

## ----report-plot, fig.cap="Comprehensive disclosure risk visualization"-------
plot(report)

## ----synthpop-integration, eval=FALSE-----------------------------------------
# library(synthpop)
# 
# # Generate synthetic data with synthpop
# synth_obj <- syn(original, seed = 123)
# 
# # Extract synthetic data for riskutility
# synthetic_sp <- synth_obj$syn
# 
# # Run disclosure assessment
# dcap_sp <- dcap(
#   X = original,
#   Y = synthetic_sp,
#   key_vars = key_vars,
#   target_var = target_var
# )
# 
# # Or use the from_synthpop helper (if available)
# # pair <- from_synthpop(synth_obj, original)

## ----simPop-integration, eval=FALSE-------------------------------------------
# library(simPop)
# 
# # Assuming you have a simPopObj object
# # synthetic_data <- simPop::pop(simPopObj)
# 
# # Use with riskutility
# # dcap_result <- dcap(X = original, Y = synthetic_data, ...)

## ----metric-comparison--------------------------------------------------------
# Extract scalar values from results
dcap_val <- dcap_result$dcap[1]
weap_val <- weap_result$weap_mean[1]
disco_val <- disco_result$pct_disco[1]
ims_val <- ims_result$ims[1]
ps_score_val <- ps_result$ps_score[1]

# Include RAPID if available
if (exists("rapid_result")) {
  rapid_val <- rapid_result$rapid[1]
  metrics <- c("DCAP", "Mean WEAP", "DiSCO %", "RAPID", "IMS", "PS Score")
  values <- round(c(dcap_val, weap_val, disco_val, rapid_val, ims_val, ps_score_val), 4)
  interp <- c(
    ifelse(dcap_val < 0.05, "Low risk",
           ifelse(dcap_val < 0.15, "Moderate", "Elevated")),
    ifelse(weap_val < 0.5, "Low risk",
           ifelse(weap_val < 0.8, "Moderate", "Elevated")),
    ifelse(disco_val < 1, "Low risk",
           ifelse(disco_val < 5, "Moderate", "Elevated")),
    ifelse(rapid_val < 0.05, "Low risk",
           ifelse(rapid_val < 0.15, "Moderate", "Elevated")),
    ifelse(ims_val < 0.01, "No memorization", "Potential memorization"),
    ifelse(ps_score_val < 0.1, "Good utility",
           ifelse(ps_score_val < 0.3, "Acceptable", "Poor utility"))
  )
} else {
  metrics <- c("DCAP", "Mean WEAP", "DiSCO %", "IMS", "PS Score")
  values <- round(c(dcap_val, weap_val, disco_val, ims_val, ps_score_val), 4)
  interp <- c(
    ifelse(dcap_val < 0.05, "Low risk",
           ifelse(dcap_val < 0.15, "Moderate", "Elevated")),
    ifelse(weap_val < 0.5, "Low risk",
           ifelse(weap_val < 0.8, "Moderate", "Elevated")),
    ifelse(disco_val < 1, "Low risk",
           ifelse(disco_val < 5, "Moderate", "Elevated")),
    ifelse(ims_val < 0.01, "No memorization", "Potential memorization"),
    ifelse(ps_score_val < 0.1, "Good utility",
           ifelse(ps_score_val < 0.3, "Acceptable", "Poor utility"))
  )
}

comparison <- data.frame(
  Metric = metrics,
  Value = values,
  Interpretation = interp
)

knitr::kable(comparison, caption = "Summary of disclosure risk and utility metrics")

## ----computational, eval=FALSE------------------------------------------------
# # Enable progress bars for long computations
# dcr_result <- dcr(
#   X = train_data,
#   Y = synthetic,
#   holdout = holdout_data,
#   progress = TRUE  # Show progress bar
# )

## ----session-info-------------------------------------------------------------
sessionInfo()

