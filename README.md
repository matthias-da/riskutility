---
editor_options: 
  markdown: 
    wrap: 72
---

# riskutility

## Functions almost ready:

`compare` ... maybe rename to `compareCDF` and split to several
functions, such as `compareDensity`, `compareTable`, ..., i.e. must noch
solve everything. `ci_overlap` ... almost ready\
`divergence_measures` ... Add more documentation from
`densitydiff_kl_num`.\
`gower` ... ready\
`mqs` ... ready\
`plot.propscore` ... ready\
`print.propscore` ... ready\
`propscore` ... ready\
`summary.propscore` ... ready

## Functions to be revised:

`densitydiff_1d_num` ... probably another name for the function, and
maybe to join together with `densitydiff_1d_cat`. New function arguments
`conditional` instead of `strata_x` and `strata_y`. Also for more than
one variable, such as function `compare`. `densitydiff_pca` ... TBD:
plot method\
`plot.denratio` ... adapt similar to `plot.compare` `densitydiff_kl_num`
... Not ready, TBD\
`gower_density` ... TBD\
`plot.denpca` ... TBD

## Missing features:

Risk measurement for synthetic data 

-   Differential Correct Attribution Probability (DCAP)
-   Within Equivalence Class Attribution Probability (WEAP) and Targeted
    Correct Attribution Probability (TCAP)
-   Inferencial disclosure
-   Holdout method
-   Mostly AI methods
-   *I do not believe that we need this: Risks from membership attacks /
    identity disclosure -assume that the attacker already knows the true
    values for some target records and uses this information to learn
    whether these units are included in the original data. (sometimes
    the fact that someone is contained in a database already reveals
    sensitive information, if the database only contains a specific
    subgroup of the population such as the Survey of Prison Inmates.*
- Kolmogorov-Smirnov distance between the two distributions 
(they call this measure SPECKS for Synthetic data generation)
-   and many more

Risk measures for traditionally anonymized data

- Identity and attribute disclosure risk measures based on record linkage
- Better formula for individual risk as in sdcMicro
- Bayesian approaches?

Risk measures for longitudinal data

- ...


Utility

-   chi\^2 statitics for comparison of tables (maybe directly included
    in mosaic-plots)
-   correlation statistics similar to synthpop.
-   PCA-biplot comparison. There are versions for categorical variables as well.
-   Utility for longitudinal data (random mixed effects models, ...)

