
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/radjust)](https://cran.r-project.org/package=radjust)

<!-- README.md is generated from README.Rmd. Please edit that file -->

### radjust: Replicability Adjusted p-values for Two Independent Studies with Multiple Endpoints

This package provides the adjusted p-values for the null hypothesis of
no replicability across studies for two study designs: a primary and
follow-up study, where the features in the follow-up study are selected
from the primary study; two independent studies, where the features for
replicability are first selected in each study separately. The latter
design is the one encountered in typical meta-analysis of two studies,
but the inference is for replicability rather than for identifying the
features that are nonnull in at least one study.

#### Installation

You can install radjust from github with:

``` r
# install.packages("devtools")
devtools::install_github("shay-y/radjust")
```

#### Examples

Two independent studies:

``` r
library(radjust)

## transform the example two-sided p-values to one-sided in the same direction (left):
## (we use the direction of the test statistic to do so and assume that it is continuous)
pv1 <- ifelse(mice$dir_is_left1, mice$twosided_pv1/2, 1-mice$twosided_pv1/2)
pv2 <- ifelse(mice$dir_is_left2, mice$twosided_pv2/2, 1-mice$twosided_pv2/2)

radjust_sym(pv1, pv2, input_type = "all", directional_rep_claim = TRUE, variant = "adaptive", alpha=0.025)
```

    #> Note:    pv1 and pv2 have the same length and don't have names
    #>  -> matching features by location.

    #> Note:    Directional replicability claim option is set to TRUE,
    #>  make sure you have entered the *left* sided p-values.

    #> 
    #>  Replicability Analysis
    #> 
    #> Call:
    #> radjust_sym(pv1 = pv1, pv2 = pv2, input_type = "all", directional_rep_claim = TRUE, 
    #>     variant = "adaptive", alpha = 0.025)
    #> 
    #> Selection (adaptive):
    #> 20 features selected in study 1.
    #> 19 features selected in study 2.
    #> 12 features selected in both studies.
    #> 
    #> Estimates for fraction of nulls among the selected:
    #> 0.4318489 in study 1.
    #> 0.4615385 in study 2.
    #> 
    #> Features selected in both studies:
    #>  name    p_value1    p_value2     r_value Direction Significant
    #>     2 1.18873e-03 1.61210e-06 0.003901482      Left           *
    #>     9 6.11236e-03 3.16097e-08 0.012538175      Left           *
    #>    14 4.34268e-05 4.77527e-03 0.012538175      Left           *
    #>    16 5.88782e-03 1.96218e-04 0.012538175      Left           *
    #>    17 1.75750e-02 3.26740e-04 0.026219142     Right            
    #>    20 1.57223e-02 6.52192e-05 0.025800620      Left            
    #>    21 2.64690e-06 2.34075e-02 0.036011533      Left            
    #>    23 3.32734e-09 5.37832e-05 0.000496460      Left           *
    #>    24 6.65468e-09 7.59238e-03 0.015574107      Left           *
    #>    25 3.32734e-09 1.37186e-05 0.000253267      Left           *
    #>    26 6.65468e-09 3.15068e-04 0.001454158      Left           *
    #>    27 6.65468e-09 9.48060e-05 0.000583421      Left           *
    #> 
    #> 12 features are significant for directional replicability claims (alpha = 0.025).

Primary and follow-up studies:

``` r
rv  <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547)
head(rv)
```

    #> [1] 6.419025e-30 2.027395e-28 5.719923e-19 6.380892e-17 6.380892e-17
    #> [6] 2.711667e-16
