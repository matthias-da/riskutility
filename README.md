# riskutility

Risk and utility measurement for anonymized and synthetic data.

## Installation

```r
# Install from GitHub
devtools::install_github("username/riskutility")
```

## Overview

The `riskutility` package provides comprehensive methods to measure disclosure risk and data utility for anonymized and synthetic data. It is designed to work standalone or integrated with the `simPop` package's S4 class structure.

## Disclosure Risk Metrics

### Attribution-Based Measures (CAP Family)

These metrics measure the probability that an adversary can correctly infer sensitive attributes from quasi-identifiers.

| Function | Description |
|----------|-------------|
| `dcap()` | Differential Correct Attribution Probability - overall disclosure risk measure |
| `tcap()` | Targeted CAP - per-record attribution probability with risk categories |
| `weap()` | Within Equivalence Class Attribution Probability - identifies risky synthetic records |
| `disco()` | Disclosive in Synthetic Correct Original - counts records leaking original information |

```r
# Example: Compute DCAP
result <- dcap(original_data, synthetic_data,
               key_vars = c("age", "gender", "region"),
               target_var = "income")
print(result)
summary(result)
plot(result)
```

### Distance-Based Measures (Holdout Method)

These metrics detect memorization by comparing distances to training vs. holdout data.

| Function | Description |
|----------|-------------|
| `dcr()` | Distance to Closest Record - detects if synthetic is too close to training |
| `nndr()` | Nearest Neighbor Distance Ratio - detects suspicious proximity patterns |
| `ims()` | Identical Match Share - percentage of exact copies |
| `repu()` | Replicated Uniques - copies of unique (singleton) training records |

```r
# Example: Compute DCR with holdout
result <- dcr(training_data, synthetic_data,
              holdout = holdout_data)
print(result)
# Or with automatic holdout split:
result <- dcr(original_data, synthetic_data,
              holdout_fraction = 0.5, seed = 42)
```

### Information-Theoretic Measures

| Function | Description |
|----------|-------------|
| `mutualInformation()` | Mutual information between variables |
| `max_info_leakage()` | Maximum information leakage measure |
| `information_surprisal()` | Information surprisal metric |
| `positive_information_disclosure()` | Positive information disclosure |
| `privacy_score()` | Overall privacy score |
| `systemAnonymityLevel()` | System anonymity level |

### Entropy Measures

| Function | Description |
|----------|-------------|
| `KLDiv()` / `KLDiv_bayes()` | Kullback-Leibler divergence |
| `JSDiv()` / `JSDiv_bayes()` | Jensen-Shannon divergence |
| `MaxEntropy()` | Maximum entropy |
| `MinEntropy()` | Minimum entropy |
| `RenyiEntropy()` | Renyi entropy |
| `ConditionalEntropy()` | Conditional entropy |
| `CumulativeEntropy()` | Cumulative entropy |
| `NormalizedEntropy()` | Normalized entropy |

## Utility Metrics

### Distribution Comparison

| Function | Description |
|----------|-------------|
| `compare_distributions_cont()` | Compare continuous distributions |
| `compare_histograms()` | Visual histogram comparison |
| `compare_boxplots()` | Boxplot comparison |
| `compare_ks_test()` | Kolmogorov-Smirnov test |
| `compare_wasserstein()` | Wasserstein distance |
| `compare_chisq_gof()` | Chi-squared goodness of fit |
| `compare_means_frequencies()` | Compare means and frequencies |

### Multivariate Comparison

| Function | Description |
|----------|-------------|
| `compare_multivariate_distribution()` | Multivariate distribution comparison |
| `compare_multivariate_summary_statistics()` | Summary statistics comparison |
| `compare_correlation_matrices()` | Correlation matrix comparison |
| `compare_pca()` | PCA-based comparison |
| `compare_embedding()` | Embedding-based comparison (t-SNE, UMAP) |

### Model-Based Utility

| Function | Description |
|----------|-------------|
| `propscore()` | Propensity score utility measure |
| `compare_model_performance()` | Compare predictive model performance |
| `compare_feature_importance()` | Compare feature importance |

### Other Utility Functions

| Function | Description |
|----------|-------------|
| `compare_missing_values()` | Compare missing value patterns |
| `compare_outliers()` | Compare outlier patterns |
| `ci_overlap()` | Confidence interval overlap |
| `gower()` | Gower distance computation |
| `mqs()` | Multivariate quality score |
| `densitydiff_1d_num()` | 1D density difference (numerical) |
| `densitydiff_kl_num()` | KL-based density difference |
| `densitydiff_pca()` | PCA-based density difference |

### Evaluation Statistics

| Function | Description |
|----------|-------------|
| `mae()` | Mean Absolute Error |
| `mse()` | Mean Squared Error |
| `rmse()` | Root Mean Squared Error |
| `mape()` | Mean Absolute Percentage Error |
| `ait()` | Average Information Transfer |

## Quick Start Example

```r
library(riskutility)

# Load or create your data
data(eusilc13puf, package = "simPop")
original <- eusilc13puf[1:500, ]

# Create synthetic data (example: shuffle sensitive variable)
synthetic <- original
synthetic$pb220a <- sample(synthetic$pb220a)

# Define quasi-identifiers and target
key_vars <- c("age", "rb090", "db040")
target_var <- "pb220a"

# === Disclosure Risk Assessment ===

# 1. Attribution probability
dcap_result <- dcap(original, synthetic, key_vars, target_var)
print(dcap_result)

# 2. Distance-based privacy check
dcr_result <- dcr(original, synthetic, holdout_fraction = 0.5, seed = 123)
print(dcr_result)

# 3. Check for exact copies
ims_result <- ims(original, synthetic)
print(ims_result)

# === Utility Assessment ===

# 4. Propensity score
ps_result <- propscore(original, synthetic, na = "remove")
print(ps_result)

# 5. Distribution comparison
compare_histograms(original, synthetic, var = "age")
```

## Interpretation Guide

### DCAP/TCAP
- **DCAP close to baseline**: Low risk - synthetic doesn't leak target info
- **DCAP >> baseline**: Elevated risk - better-than-random inference possible
- **Risk ratio > 1.5**: Generally considered elevated risk

### DCR (Distance to Closest Record)
- **DCR ratio ~ 1.0**: Good privacy - synthetic equally distant from train/holdout
- **DCR ratio < 0.9**: Privacy concern - synthetic too close to training
- **DCR share ~ 50%**: Ideal - half closer to train, half to holdout

### NNDR (Nearest Neighbor Distance Ratio)
- **NNDR ~ 1**: Good - nearest neighbors similarly distant
- **NNDR ~ 0**: Suspicious - one neighbor much closer (potential copy)
- **Many NNDR < 0.1**: Privacy concern

### IMS (Identical Match Share)
- **IMS = 0%**: Ideal - no exact copies
- **IMS < 1%**: Acceptable - likely coincidental
- **IMS > 5%**: Serious concern - significant copying

## References

- Taub, J., Elliot, M., Pampaka, M., & Smith, D. (2018). Differential Correct Attribution Probability for Synthetic Data: An Exploration. *Privacy in Statistical Databases*, 122-137.
- Templ, M. (2017). Statistical Disclosure Control for Microdata: Methods and Applications in R. *Springer International Publishing*.
- MOSTLY AI (2024). Synthetic Data Quality Metrics. https://docs.mostly.ai/

## Development Status

### Ready
- All disclosure risk metrics (DCAP, TCAP, WEAP, DiSCO, DCR, NNDR, IMS)
- Propensity score utility measure
- Distribution comparison functions
- Entropy and divergence measures

### To Be Revised
- `densitydiff_1d_num` - consider merging with categorical version
- `densitydiff_pca` - plot method needed
- `gower_density` - under development

### Planned Features
- SPECKS (KS-based synthetic data measure)
- Inferential disclosure measures
- Record linkage based identity disclosure
- Longitudinal data metrics
- Comprehensive reporting function

## License

GPL-3

## Author

Matthias Templ
