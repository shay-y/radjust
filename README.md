
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/radjust)](https://cran.r-project.org/package=radjust)
[![Build
Status](https://travis-ci.org/shay-y/radjust.svg?branch=master)](https://travis-ci.org/shay-y/radjust)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# radjust: Replicability Adjusted p-values for Two Independent Studies with Multiple Endpoints

Given p-values from two independent studies with multiple endpoints
(features), the functions in the package return the adjusted p-values
for false discovery rate control on replicability claims.

In replicability analysis we seek to reject the null hypothesis of no
replicability in favor of the alternative hypothesis of replicability:
that the finding replicated across the studies.  
We do so by testing for signal in both studies. This is in contrast to a
typical meta-analysis, where the test can also reject when only a single
study has signal.

The procedures implemented in the functions compute the adjusted
p-values for FDR control on replicability claims. By declaring as
replicability discoveries the features with adjusted p-values (termed
r-values) below the desired nominal level (e.g., 0.05), the FDR on
replicability claims is controlled at the nominal level. See Bogomolov
and Heller (2013), Heller, Bogomolov and Benjamini (2014), and Bogomolov
and Heller (2018) for details.

The function `radjust_sym` should be used for replicability analysis of
two independent studies, each examining multiple features. The features
for replicability are first selected in each study separately based on
the results of that study.

The function `radjust_pf` should be used for replicability analysis of a
primary study and an independent follow-up study, where the features in
the follow-up study are selected from the primary study.

## Examples

Using `radjust_sym`:

``` r
library(radjust)

## transform the example two-sided p-values to one-sided in the same direction (left):
## (we use the direction of the test statistic to do so and assume that it is continuous)
pv1 <- ifelse(mice$dir_is_left1, mice$twosided_pv1/2, 1-mice$twosided_pv1/2)
pv2 <- ifelse(mice$dir_is_left2, mice$twosided_pv2/2, 1-mice$twosided_pv2/2)

radjust_sym(pv1, pv2, input_type = "all", directional_rep_claim = TRUE, variant = "adaptive", alpha=0.05)
```

    > Note: pv1 and pv2 have the same length and don't have names
    >   -> matching features by location.

    > Note: Directional replicability claim option is set to TRUE.
    >   Make sure you have entered the *left* sided p-values.

    > 
    >   Replicability Analysis
    > 
    > Call:
    > radjust_sym(pv1 = pv1, pv2 = pv2, input_type = "all", directional_rep_claim = TRUE, 
    >     variant = "adaptive", alpha = 0.05)
    > 
    > Selection (adaptive):
    > 20 features selected in study 1.
    > 19 features selected in study 2.
    > 12 features selected in both studies.
    > 
    > Estimates for fraction of nulls among the selected in the other study:
    > 0.4432133 in study 1.
    > 0.4736842 in study 2.
    > 
    > Features selected in both studies:
    >  name    p_value1    p_value2     r_value Direction Significant
    >     2 1.18873e-03 1.61210e-06 0.004004153      Left           *
    >     9 6.11236e-03 3.16097e-08 0.012868127      Left           *
    >    14 4.34268e-05 4.77527e-03 0.012868127      Left           *
    >    16 5.88782e-03 1.96218e-04 0.012868127      Left           *
    >    17 1.75750e-02 3.26740e-04 0.026909119     Right           *
    >    20 1.57223e-02 6.52192e-05 0.026479584      Left           *
    >    21 2.64690e-06 2.34075e-02 0.036959205      Left           *
    >    23 3.32734e-09 5.37832e-05 0.000509525      Left           *
    >    24 6.65468e-09 7.59238e-03 0.015983952      Left           *
    >    25 3.32734e-09 1.37186e-05 0.000259932      Left           *
    >    26 6.65468e-09 3.15068e-04 0.001492426      Left           *
    >    27 6.65468e-09 9.48060e-05 0.000598774      Left           *
    > 
    > 12 features are discovered in the directional replicability analysis (alpha = 0.05).

Primary and follow-up studies (`radjust_pf`):

``` r
rv  <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547)
head(rv)
```

    > [1] 6.419025e-30 2.027395e-28 5.719923e-19 6.380892e-17 6.380892e-17
    > [6] 2.711667e-16

## Installation

You can install radjust from github with:

``` r
# install.packages("devtools")
devtools::install_github("shay-y/radjust")
```

## How to cite

Use the `citation()` R function:

``` r
citation("radjust")
```

    > 
    > To cite radjust in publications, please use:
    > 
    >   Shay Yaacoby, Marina Bogomolov and Ruth Heller (2018). radjust:
    >   Replicability Adjusted p-values for Two Independent Studies with
    >   Multiple Endpoints. R package version 0.1.0.
    > 
    > To cite radjust_sym(), add:
    > 
    >   Bogomolov, M. and Heller, R. (2018). Assessing replicability of
    >   findings across two studies of multiple features. Biometrika.
    > 
    > To cite radjust_pf(), add:
    > 
    >   Bogomolov, M. and Heller, R. (2013). Discovering findings that
    >   replicate from a primary study of high dimension to a follow-up
    >   study. Journal of the American Statistical Association, Vol.
    >   108, No. 504, Pp. 1480-1492.
    > 
    >   Heller, R., Bogomolov, M., & Benjamini, Y. (2014). Deciding
    >   whether follow-up studies have replicated findings in a
    >   preliminary large-scale omics study. Proceedings of the National
    >   Academy of Sciences of the United States of America, Vol. 111,
    >   No. 46, Pp. 16262–16267.
    > 
    > To see these entries in BibTeX format, use 'print(<citation>,
    > bibtex=TRUE)', 'toBibtex(.)', or set
    > 'options(citation.bibtex.max=999)'.
