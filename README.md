
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/radjust)](https://cran.r-project.org/package=radjust)

<!-- README.md is generated from README.Rmd. Please edit that file -->

### radjust: Replicability Adjusted p-values for Two Independent Studies with Multiple Endpoints

Given p-values from two independent studies with multiple endpoints
(features), the functions in the package return the adjusted p-values
for false discovery rate control on replicability claims.

In replicability analysis we seek to reject the null hypothesis of no
replicability in favor of the alternative hypothesis of replicability:
that the findings were replicated across the studies.  
We do so by testing for rejection in **both** studies. This is in
contrast to a typical meta-analysis, where the test can also reject when
only a single finding is extreme enough.

The procedures implemented in the functions guarantee false discovery
rate control (when testing multiple endpoints in each study) by
comparing the adjusted p-values to the rate threshold (typically alpha =
0.05).

The function `radjust_sym` fits to a design of two studies, where the
features for replicability are first selected in each study separately.

The function `radjust_pf` fits to a design of primary and follow-up
studies, where the features in the follow-up study are selected from the
primary study.

#### Examples

Using `radjust_sym`:

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

Primary and follow-up studies (`radjust_pf`):

``` r
rv  <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547)
head(rv)
```

    #> [1] 6.419025e-30 2.027395e-28 5.719923e-19 6.380892e-17 6.380892e-17
    #> [6] 2.711667e-16

#### Installation

You can install radjust from github with:

``` r
# install.packages("devtools")
devtools::install_github("shay-y/radjust")
```

#### How to cite

Use the `citation()` R function:

``` r
citation("radjust")
```

    #> 
    #> To cite radjust in publications, please use:
    #> 
    #>   Shay Yaacoby, Marina Bogomolov and Ruth Heller (2018). radjust:
    #>   Replicability Adjusted p-values for Two Independent Studies with
    #>   Multiple Endpoints. R package version 0.1.0.
    #> 
    #> To cite radjust_sym(), add:
    #> 
    #>   Bogomolov, M. and Heller, R. (2018). Assessing replicability of
    #>   findings across two studies of multiple features. Biometrika.
    #> 
    #> To cite radjust_pf(), add:
    #> 
    #>   Bogomolov, M. and Heller, R. (2013). Discovering findings that
    #>   replicate from a primary study of high dimension to a follow-up
    #>   study. Journal of the American Statistical Association, Vol.
    #>   108, No. 504, Pp. 1480-1492.
    #> 
    #> To see these entries in BibTeX format, use 'print(<citation>,
    #> bibtex=TRUE)', 'toBibtex(.)', or set
    #> 'options(citation.bibtex.max=999)'.
