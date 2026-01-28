# riskutility 0.1.0

* Initial CRAN submission

## Disclosure Risk Metrics

* `dcap()`: Differential Correct Attribution Probability
* `tcap()`: Targeted Correct Attribution Probability with risk categories
* `weap()`: Within Equivalence Class Attribution Probability
* `disco()`: Disclosive in Synthetic Correct Original
* `dcr()`: Distance to Closest Record
* `nndr()`: Nearest Neighbor Distance Ratio
* `ims()`: Identical Match Share
* `repu()`: Replicated Uniques
* `disclosure_report()`: Comprehensive disclosure risk report

## Utility Metrics

* `propscore()`: Propensity score utility measure
* Distribution comparison functions: `compare_distributions_cont()`,
  `compare_histograms()`, `compare_boxplots()`, etc.
* Statistical tests: `compare_ks_test()`, `compare_chisq_gof()`,
  `compare_wasserstein()`
* Multivariate comparisons: `compare_correlation_matrices()`,
  `compare_pca()`, `compare_embedding()`

## Information-Theoretic Measures

* Entropy functions: `KLDiv()`, `JSDiv()`, `MaxEntropy()`, `MinEntropy()`, etc.
* `mutualInformation()`: Mutual information between variables
* `privacy_score()`: Overall privacy composite score

## S3 Methods

* All major functions return S3 objects with `print()`, `summary()`, and
  `plot()` methods for consistent user experience
