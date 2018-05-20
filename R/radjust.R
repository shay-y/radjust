#' @title radjust: Replicability Adjusted p-values for Two Independent Studies with Multiple Endpoints
#'
#' @description This package provides the adjusted p-values for the null hypothesis of no replicability across studies for two study designs: a primary and follow-up study, where the features in the follow-up study are selected from the primary study; two independent studies, where the features for replicability are first selected in each study separately. The latter design is the one encountered in typical meta-analysis of two studies, but the inference is for replicability rather than for identifying the features that are nonnull in at least one study.
#'
#' @docType package
#' @name radjust
NULL

#' @title p-values of 29 Behavioural Measures in Two Studies
#'
#' @description In different laboratories, the comparison of behaviours of the same two strains of mice may lead to opposite conclusions that are both statistically significant.
#' An explanation may be the different laboratory environment, i.e., personnel, equipment, or measurement techniques,
#' affecting differently the study strains.
#' This data set provides the p-values for testing the association of mice strain with 29 behavioural measures from five commonly used behavioural tests  in
#' two laboratories: the laboratory of H. Wurbel at the University of Giessen, and the laboratory of P. Gass at the  Central Institute of Mental Health, Mannheim.
#' The data table contains two-sided p-values. To transform all the two-sided p-values to one sided in the same direction, see the example in \code{\link{radjust_sym}}.
#'
#' @format A data frame with 29 rows and 5 columns:
#' \tabular{lll}{
#' \code{feature_name} \tab char.   \tab The name of the measure and the test, concatenated. \cr
#' \code{twosided_pv1} \tab numeric \tab the \emph{two-sided} p-value from study 1. \cr
#' \code{twosided_pv2} \tab numeric \tab the \emph{two-sided} p-value from study 2. \cr
#' \code{dir_is_left1} \tab logical \tab whether the direction of the test statistic from study 1 is \emph{left}. \cr
#' \code{dir_is_left2} \tab logical \tab  whether the direction of the test statistic from study 2 is \emph{right}. \cr
#' }
#'
#' @source Richter, S. Helene, et al. "Effect of population heterogenization on the reproducibility of mouse behavior:
#'  a multi-laboratory study." PLoS One 6.1 (2011): e16461.
#'
#' @seealso \code{\link{radjust_sym}}
#'
#' @references  Bogomolov, M. and Heller, R. (2018). Assessing replicability of findings across two studies of multiple
#' features. Biometrika.
#'
#'
"mice"

#' @title p-values of 126 SNPs followed from a primary study to a follow-up study for testing their association
#' with Crohn's disease.
#'
#' @description To discover the associations between SNPs and Crohn's disease, 635547 SNPs were examined in a primary study. For follow-up, 126 SNPs were measured in an independent study.
#' The criteria for follow-up were as follows: the two smallest p-values in each distinct region with primary study p-values below 0.00005.
#'
#' @format A data frame with 126 rows and 3 columns:
#' \tabular{lll}{
#' \code{index} \tab integer \tab just the row number. \cr
#' \code{pv1}   \tab numeric \tab p-value from study 1. \cr
#' \code{pv2}   \tab numeric \tab p-value from study 2.
#' }
#' @source Barrett, Jeffrey C., et al. "Genome-wide association defines more than 30 distinct susceptibility loci for
#'  Crohn's disease." Nature genetics 40.8 (2008): 955.
#'
#' @seealso the example in \code{\link{radjust_pf}} uses this data.
#'
#' @references Bogomolov, M. and Heller, R. (2013). Discovering findings that replicate from a primary study of high dimension to a follow-up study.
#' Journal of the American Statistical Association, Vol. 108, No. 504, Pp. 1480-1492.
#'
"crohn"
