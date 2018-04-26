#' @title radjust: Replicability Adjusted p-values for Two Independent Studies with Multiple Endpoints
#'
#' @docType package
#' @name radjust
NULL

#' @title p-values of 29 Behavioural Measures in Two Studies
#'
#' @description EDIT ... note that the table contains two-sided p-values. To transform all
#'  the two-sided p-values to one sided in the same direction, see the example in \code{\link{radjust_sym}}.
#'
#' @format A data frame with 29 rows and 5 columns:
#' \tabular{lll}{
#' \code{feature_name} \tab char.   \tab The name of the measure and the test, concatenated. \cr
#' \code{twosided_pv1} \tab numeric \tab the \emph{two-sided} p-value from study 1. \cr
#' \code{twosided_pv2} \tab numeric \tab EDIT \cr
#' \code{dir_is_left1} \tab logical \tab whether the direction of the test statistic from study 1 is \emph{left}. \cr
#' \code{dir_is_left2} \tab logical \tab EDIT
#' }
#'
#' @source Richter, S. Helene, et al. "Effect of population heterogenization on the reproducibility of mouse behavior:
#'  a multi-laboratory study." PLoS One 6.1 (2011): e16461.
#'
#' @seealso \code{\link{radjust_sym}}
#' @references  Bogomolov, M. and Heller, R. (2018). Assessing replicability of findings across two studies of multiple
#' features. Biometrika.
#'
#'
"mice"

#' @title EDIT
#'
#' @description EDIT ...
#'
#' @format A data frame with 126 rows and 3 columns:
#' \tabular{lll}{
#' \code{index} \tab integer \tab just the row number. \cr
#' \code{pv1}   \tab numeric \tab p-value from study 1. \cr
#' \code{pv2}   \tab numeric \tab p-value from study 2.
#' }
#' @source Barrett, Jeffrey C., et al. "Genome-wide association defines more than 30 distinct susceptibility loci for
#'  Crohn's disease." Nature genetics 40.8 (2008): 955.
#' @seealso the example in \code{\link{radjust_pf}} uses this data.
#' @references Bogomolov, M. and Heller, R. (2013). Discovering findings that replicate from a primary study of high dimension to a follow-up study.
#' Journal of the American Statistical Association, Vol. 108, No. 504, Pp. 1480-1492.
#'
"crohn"
