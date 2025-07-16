library(riskutility)
library(sdcMicro)

simulate_data <- function(n = 1000,
                          n_quasi = 3,
                          cat_levels_quasi = rep(3, 3),
                          corr_quasi = 0.3,
                          sens_cat_levels = 3,
                          sens_cat_corr = 0.2,
                          sens_cont_corr = 0.5,
                          seed = 42) {

  set.seed(seed)

  # Helper function: Generate correlated ordinal data
  generate_correlated_quasi <- function(n, n_quasi, cat_levels, rho) {
    library(MASS)
    Sigma <- matrix(rho, n_quasi, n_quasi)
    diag(Sigma) <- 1
    z <- mvrnorm(n, mu = rep(0, n_quasi), Sigma = Sigma)

    quasi_list <- vector("list", n_quasi)
    names(quasi_list) <- paste0("QI_", seq_len(n_quasi))

    for (j in seq_len(n_quasi)) {
      q <- cut(z[, j],
               breaks = quantile(z[, j], probs = seq(0, 1, length.out = cat_levels[j] + 1)),
               include.lowest = TRUE, labels = FALSE)
      quasi_list[[j]] <- as.factor(q)
    }

    quasi_df <- as.data.frame(quasi_list)
    return(quasi_df)
  }

  # Generate quasi-identifiers
  if (length(cat_levels_quasi) != n_quasi) {
    cat_levels_quasi <- rep(cat_levels_quasi[1], n_quasi)
  }
  quasi <- generate_correlated_quasi(n, n_quasi, cat_levels_quasi, corr_quasi)

  # Generate numeric latent variable for correlation with sensitive vars
  latent <- rowMeans(model.matrix(~ . -1, data = quasi))  # numeric representation

  # Create sensitive categorical variable with correlation to latent
  z_cat <- latent + rnorm(n, sd = sqrt((1 - sens_cat_corr^2) / sens_cat_corr^2))
  sens_cat <- cut(z_cat,
                  breaks = quantile(z_cat, probs = seq(0, 1, length.out = sens_cat_levels + 1)),
                  include.lowest = TRUE, labels = paste0("S", seq_len(sens_cat_levels)))
  sens_cat <- factor(sens_cat)

  # Create sensitive continuous variable with correlation to latent
  z_cont <- scale(latent)
  noise <- rnorm(n, mean = 0, sd = sqrt(1 - sens_cont_corr^2))
  sens_cont <- scale(sens_cont_corr * z_cont + noise)

  # Final dataset
  df <- cbind(quasi,
              sensitive_cat = sens_cat,
              sensitive_cont = as.numeric(sens_cont))

  return(as.data.frame(df))
}

df1 <- simulate_data(n = 1000)


anonymize_data <- function(data,
                           method = c("kAnon", "PRAM", "swap", "recoding", "synthpop", "simPop"),
                           qid = 1:(ncol(data) - 2),
                           k = 3,
                           pram_vars = NULL,
                           swap_rate = 0.07,
                           recode_vars = NULL,
                           synthpop_args = list(),
                           simPop_args = list()) {
  stopifnot(requireNamespace("sdcMicro", quietly = TRUE))

  method <- match.arg(method)

  if (method == "kAnon") {
    sdc <- createSdcObj(dat = data, keyVars = qid)
    sdc <- localSuppression(sdc, k = k)
    return(extractManipData(sdc))

  } else if (method == "PRAM") {
    if(is.null(pram_vars)) pram_vars <- colnames(data[, 1:(ncol(data)-1)])
    sdc <- createSdcObj(dat = data, keyVars = pram_vars)
    sdc <- pram(sdc, variables = pram_vars)
    return(extractManipData(sdc))

  } else if (method == "swap") {
   return(rankSwap(data, variables = colnames(data[, 1:(ncol(data) - 1)])))

  } else if (method == "recoding") {
    stopifnot(!is.null(recode_vars))
    sdc <- createSdcObj(dat = data, keyVars = recode_vars)
    for (v in recode_vars) {
      sdc <- globalRecode(sdc, column = v, breaks = 4)
    }
    return(extractManipData(sdc))

  } else if (method == "synthpop") {
    stopifnot(requireNamespace("synthpop", quietly = TRUE))
    syn_args <- modifyList(list(data = data), synthpop_args)
    syn_result <- do.call(synthpop::syn, syn_args)
    return(syn_result$syn)
  }
}


# Generate simulated data first
sim_data <- simulate_data(n = 200, n_quasi = 3, cat_levels_quasi = 3,
                          corr_quasi = 0.3,
                                    sens_cat_levels = 4, sens_cat_corr = 0.3,
                          sens_cont_corr = 0.3)

# Anonymize using k-anonymity
anon_data_k <- anonymize_data(sim_data, method = "kAnon", k = 3)

# Anonymize using PRAM
anon_data_pram <- anonymize_data(sim_data, method = "PRAM", k = 3)

# recordSwap
anon_data_swap <- anonymize_data(sim_data, method = "swap", k = 3)

# Anonymize using synthetic data with synthpop
anon_data_syn_cart <- anonymize_data(sim_data, method = "synthpop")

# Anonymize using synthetic data with synthpop
anon_data_syn_rf <- anonymize_data(sim_data, method = "synthpop", synthpop_args = list(method = "rf"))

# Anonymize using synthetic data with synthpop
anon_data_syn_ranger <- anonymize_data(sim_data, method = "synthpop", synthpop_args = list(method = "ranger"))

# Anonymize using synthetic data with synthpop
anon_data_syn_bag <- anonymize_data(sim_data, method = "synthpop", synthpop_args = list(method = "bag"))

## Calculate risk

source("~/workspace25/SDG_Comparison/code/auxiliary_functions/risk_inferential.R", echo=TRUE)

risk_inferential(sim_data,
                 synth = anon_data_k,
                 model = "rf",
                 known_vars = c("QI_1", "QI_2", "QI_3"),
                 sensitive_var = c("sensitive_cat"))


