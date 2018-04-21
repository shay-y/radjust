#' @title Adjust p-values for Replicability across Two Independent Studies with Multiple Endpoints
#'
#' @description  Given two vectors of p-values from two independent studies, returns the adjusted p-values for
#'  false discovery rate control on replicability claims.
#'
#' @param pv1,pv2 numeric vectors of p-values. Can be the selected features from each study (default input type) from the two studies, corresponding to the selected features
#'  (the default) or all the features from each study. Can be either of the same length (so the same location in each vector corresponds to the same feature) or with names (so the same name in each vector correspond to the same feature).
#' @param w1 fraction between zero and one, of the relative weight for the p-values from study 1. Default value is 0.5
#'  (see Details for other values).
#' @param input_type whether \code{pv1} and \code{pv2} contain all the p-values from each study or only selected ones(the default).
#' @param general_dependency TRUE or FALSE, indicating whether to correct for general dependency.
#' The recommended default value is FALSE (see Details).
#' @param directional_rep_claim EDIT
#' @param variant A character string specifying the chosen variant for a potential increase in the number of discoveries.
#' Must be one of \code{"non-adaptive-with-alpha-selection"} (default), \code{"adaptive"}, or \code{"non-adaptive"} (see Details).
#' @param alpha the threshold for the selection method.
#'
#' @return The function returns a list with the following elements:
#'
#'   \tabular{ll}{
#'     \code{call}  \tab  the function call. \cr
#'     \code{inputs}  \tab  a list with the function's input parameters (except \code{pv1} and \code{pv2}). \cr
#'     \code{results_table} \tab  a data frame with the features selected in both studies and their r-values. The data frame includes 6 columns: \cr
#'     \code{selected1}  \tab the features selected in study 1 (when the variant is either \code{"adaptive"} or \code{"non-adaptive-with-alpha-selection"}). \cr
#'     \code{selected2}  \tab the features selected in study 2, same as above. \cr
#'     \code{n_selected1} \tab  the number of features in study 1. \cr
#'     \code{n_selected2} \tab  the number of features in study 2. \cr
#'     \code{pi1}  \tab  the estimate of the true-nulls fraction in the study1.\cr
#'     \code{pi2}  \tab  the estimate of the true-nulls fraction in the study2.
#'   }
#'
#' The third element in the list, \code{results_table}, includes the following columns:#'
#'
#'   \tabular{lll}{
#'     \code{name}         \tab char.   \tab the name of the feature as extracted from the named vectors, or the location, if the input vectors are not named. \cr
#'     \code{p.value.1}    \tab numeric \tab the one-sided p-value from study 1 as inputed. In the case of \code{directional_rep_claim==TRUE} the 'real' one-sided p-values is presented, i.e: \code{pmin(pv1,1-pv1}. \cr
#'     \code{p.value.2}    \tab numeric \tab same as in \code{p.value.1}. \cr
#'     \code{r.value}      \tab numeric \tab the replicability adjusted p-value (= r-value). \cr
#'     \code{Direction}    \tab char.   \tab the direction of the replicability claim, when \code{directional_rep_claim==TRUE}. \cr
#'     \code{Significant}  \tab char.   \tab a notion whether the replicability claim is significant, when \code{variant!="non-adaptive"}.
#'   }
#'
#' @details For FDR control at level \eqn{\alpha} on replicability claims, declare all features with \strong{\eqn{r}-value} at most \eqn{\alpha}
#' as replicated. In addition, the discoveries from study 1 among the replicability claims have an FDR control guarantee at level \eqn{w_{1}\alpha}{w1 * \alpha}.
#'  Similarly, the discoveries from study 2 among the replicability claims have an FDR control guarantee at level \eqn{(1-w_{1})\alpha}{(1-w1) * \alpha}. Setting
#'  a value of \eqn{w_{1}}{w1} different than half is appropriate if stricter FDR control is desired for one of the studies.
#'  For example, if study two has a much larger sample size than study one (and both studies examine the same problem), then
#'  setting \eqn{w_{1} > 0.5}{w1 > 0.5} will provide a stricter FDR control for the larger study and greater power for the replicability analysis,
#'   see Bogomolov and Heller (2018) for details.
#'
#'   The theoretical FDR control guarantees assume independence within each vector of p-values. However, empirical
#'   investigations suggest that the method is robust to deviations from independence. In practice, we recommend using it whenever the
#'   Benjamini-Hochberg procedure is appropriate for use with single studies, as this procedure can be viewed as a two-dimensional
#'  Benjamini-Hochberg procedure which enjoys similar robustness properties. For general dependence, we provide the option to apply
#'  a more conservative procedure with theoretical FDR control guarantee for any type of dependence for the non-adaptive procedure,
#'  by setting  \code{general_dependency} to TRUE.
#'
#' If \code{variant} is \code{"non-adaptive-with-alpha-selection"}, then for a user specified \code{alpha} (default 0.05) only p-values from
#' study one below \eqn{w_{1}\alpha}{w1 * \alpha} and from study
#' two below \eqn{(1-w_{1})\alpha}{(1-w1) * \alpha} are considered for replicability analysis. This additional step prevents
#' including in the selected sets features that cannot be discovered as replicability claims at the nominal FDR level
#' \eqn{\alpha}, thus reducing the multiple adjustment necessary for replicability analysis.  If \code{variant} is \code{"adaptive"}, then for a user specified \code{alpha}
#' the adaptive replicability analysis procedure is applied on the dataset, see Bogomolov and Heller (2018) for details.
#'
#' @references Bogomolov, M. and Heller, R. (2018). Assessing replicability of findings across two studies of multiple
#' features. Biometrika.
#'
#' @seealso \code{\link{radjust_pf}} for replicability analysis in primary and follow-up studies.
#'
#' @export
#'
#' @examples
#'
#' data(mice)
#' ## transform the two-sided p-values to one-sided in the same direction (left):
#' ## (we use the direction of the test statistic to do so and assume that it is continuous)
#'
#' pv1 <- ifelse(mice$dir_is_left1, mice$twosided_pv1/2, 1-mice$twosided_pv1/2)
#' pv2 <- ifelse(mice$dir_is_left2, mice$twosided_pv2/2, 1-mice$twosided_pv2/2)
#'
#' ## run the examples as in the article:
#'
#' mice_rv_adaptive <- radjust_sym(pv1, pv2, input_type = "all", directional_rep_claim = TRUE,
#'                                 variant = "adaptive", alpha=0.025)
#' print(mice_rv_adaptive)
#'
#' mice_rv_non_adpt_sel <- radjust_sym(pv1, pv2, input_type = "all", directional_rep_claim = TRUE,
#'                                     variant = "non-adaptive-with-alpha-selection", alpha=0.025)
#' print(mice_rv_non_adpt_sel)
#'
#' mice_rv_non_adpt <- radjust_sym(pv1, pv2, input_type = "selected", directional_rep_claim = TRUE,
#'                                 variant = "non-adaptive")
#' str(mice_rv_non_adpt)
#'

radjust_sym <- function(pv1,
                         pv2,
                         w1 = 0.5,
                         input_type = c("selected_features", "all_features"),
                         general_dependency = FALSE,
                         directional_rep_claim = FALSE,
                         variant = c("non-adaptive-with-alpha-selection", "adaptive", "non-adaptive"),
                         alpha = if (variant == "non-adaptive") NULL else 0.05
)
{

  ## ---- validate inputs ----
  stopifnot(is.vector(pv1, mode = "numeric"), pv1 <= 1, 0 < pv1,
            is.vector(pv2, mode = "numeric"), pv2 <= 1, 0 < pv2,
            is.vector(w1, mode = "numeric"),  w1 < 1, 0 < w1, length(w1) == 1)
  stopifnot(is.logical(general_dependency),
            is.logical(directional_rep_claim))

  variant <- match.arg(variant)
  if (variant == "non-adaptive" && !is.null(alpha))
    warning("variant == 'non-adaptive', alpha value is ignored")
  if (variant != "non-adaptive" && is.null(alpha))
    stop("variant != 'non-adaptive' but alpha value is not specified")
  input_type <- match.arg(input_type)

  ## ---- generate messages regarding the matching method (i.e by names or positions)----
  if (has_some_names(pv1) || has_some_names(pv2))
  {
    if (has_all_names(pv1) && has_all_names(pv2))
    {
      if (length(pv1) != length(pv2))
      message("Note:\tpv1 and pv2 differ in length and have names
\t-> matching features by names.")
      else
      message("Note:\tpv1 and pv2 have the same length and have names
\t-> matching features by names.")
    } else
        stop("Some of the names are missing. Check with names(pv1) or names(pv2)")
  } else {
    if (length(pv1) != length(pv2))
      stop("\tpv1 and pv2 differ in length and don't have names
\tCannot match features.")
    else
      message(
        "Note:\tpv1 and pv2 have the same length and don't have names
\t-> matching features by location."
      )
  }

  ## ---- do the matching (by names or positions) ----
  if (has_all_names(pv1) && has_all_names(pv2))
  {
    # The following will produce NAs where there is no match, so to get lengths and sums in the next sections we use sum( , na.rm = TRUE).
    pv_union_df <- merge(pv1, pv2 , by.x = names(pv1), by.y = names(pv2))
    pv1 <- pv_union_df$pv1
    names(pv1) <- pv_union_df$names.pv1
    pv2 <- pv_union_df$pv2
    names(pv2) <- pv_union_df$names.pv2
  } else
    names(pv1) <- names(pv2) <- (1:length(pv1)) # If there are no names, then we assign indices as names.
  # From now on pv1, pv2 have the same length.

  ## ---- handle directional claims case----
  if (directional_rep_claim)
  {
    message(
      "Note:\tDirectional replicability claim option is set to TRUE,
\tmake sure you have entered the *left* sided p-values."
  # TODO: should we allow any side, given it is the same one? this while change line 119 also...
    )

    # since the user enters left p-values only, we need the 'real' p-values for some of the calculations.
    pv1real  <- ifelse(pv1 <= 0.5, pv1, 1 - pv1)
    pv2real  <- ifelse(pv2 <= 0.5, pv2, 1 - pv2)

    rep_claim_direction <- ifelse(pv1 <= 0.5 & pv2 <= 0.5, "Left",ifelse(pv1 > 0.5 & pv2 > 0.5, "Right","Opposite"))

    pv1tag <- ifelse(pv2 <= 0.5, pv1, 1 - pv1) # as in the article
    pv2tag <- ifelse(pv1 <= 0.5, pv2, 1 - pv2)

  } else {
    pv1real  <- pv1 # we allow pv > 0.5 in the non-directional claims case
    pv2real  <- pv2
    pv1tag <- pv1
    pv2tag <- pv2
  }

  ## ---- generate the logical s1, s2 and s12 vectors - the selections in each study and in both ----
  # Note that NAs in p1 and p2 are kept NAs in the selections vectors
  if (variant == "non-adaptive")
  {
    if (input_type == "all_features")
      stop("The variant cannot be 'non-adaptive' when input_type equals to 'all_features'")
    s1  <- ifelse(is.na(pv1), NA, TRUE) # "select" all except NAs
    s2  <- ifelse(is.na(pv2), NA, TRUE)
    R1 <- sum(s1, na.rm = TRUE)
    R2 <- sum(s2, na.rm = TRUE)
    pi1 <- 1
    pi2 <- 1
  } else
   if (variant == "non-adaptive-with-alpha-selection")
    {
      s1 <- (pv1real <= w1 * alpha)
      s2 <- (pv2real <= (1 - w1) * alpha)
      R1 <- sum(s1, na.rm = TRUE)
      R2 <- sum(s2, na.rm = TRUE)
      pi1 <- 1
      pi2 <- 1
    } else
      if (variant == "adaptive")
      {
        s1 <- (pv1real <= alpha)
        s2 <- (pv2real <= alpha)

        R1 <- sum(s1, na.rm = TRUE)
        R2 <- sum(s2, na.rm = TRUE)

        # the fraction of true nulls in study 1 among the selected in 2
        pi1 <- (1 + sum(pv1tag[s2] > alpha, na.rm = TRUE)) / (R2 * (1 - alpha))
        # the fraction of true nulls in study 2 among the selected in 1
        pi2 <- (1 + sum(pv2tag[s1] > alpha, na.rm = TRUE)) / (R1 * (1 - alpha))
      }

  if (directional_rep_claim)
    s12 <- (s1 & s2 & rep_claim_direction != "Opposite")
  else
    s12 <- (s1 & s2)

   ## ---- estimate pi1 and pi2 (when adaptive method is selected) ----

  if (general_dependency)
  {
    d1 <- sum(1 / (1:R2))
    d2 <- sum(1 / (1:R1))
  }
  else
  {
    d1 <- 1
    d2 <- 1
  }

  ## ---- the adjustment procedure ----
  z <- pmax(pi1 * d1 * pv1tag * R2 / w1, pi2 * d2 * pv2tag * R1 / (1 - w1))[s12] #TODO: verify that NAs are dropped here
  oz <- order(z, decreasing = TRUE)
  ozr <- order(oz)
  rv <- cummin((z / rank(z, ties.method = "max"))[oz])[ozr]
  rv <- pmin(rv, 1)
  names(rv) <- names(pv1)[s12]

  ## ---- the table of the intersection ----
  tbl_s12 <- data.frame(
      `name` = names(pv1)[s12],
      `p_value1` = pv1tag[s12],
      `p_value2` = pv2tag[s12],
      `r_value`  = rv
    )

  if (directional_rep_claim) tbl_s12$`Direction` <- rep_claim_direction[s12]

  if (!(variant == "non-adaptive")) tbl_s12$`Significant` <- ifelse(rv <= alpha,"*","")
  # TODO: consider to replace to TRUE\FALSE (and use symnum() in th print method).

  ## ---- generate the output ----

  cl <- match.call()
  inputs <- list(
    w1 = w1,
    input_type = input_type,
    general_dependency = general_dependency,
    directional_rep_claim = directional_rep_claim,
    variant = variant,
    alpha = alpha)

  output <- list(
    call = cl,
    inputs = inputs,
    results_table = tbl_s12,
    selected1 = pv1real[s1],
    selected2 = pv2real[s2],
    n_selected1 = R1,
    n_selected2 = R2,
    pi1 = if (variant == "adaptive") pi1 else NULL,
    pi2 = if (variant == "adaptive") pi2 else NULL
  )

  class(output) <- "radjust"
  return(output)
  ## ---- END ----
}

#' @exportClass radjust

# TODO: document the print method or add to radjust_sym documentation.
#' @export
print.radjust <- function (x, digits_df = max(3L, getOption("digits") - 1L), ...)
{
  if (!inherits(x, "radjust"))
    warning("calling print.radjust(<fake-radjust-object>) ...")
  cat("\n")
  cat("\tReplicability Analysis\n")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Selection (", if (x$inputs$variant == "non-adaptive") "as-inputed"  else x$inputs$variant,"):\n", sep = "")
  cat(x$n_selected1," features selected in study 1.\n", sep = "")
  cat(x$n_selected2," features selected in study 2.\n", sep = "")
  cat(nrow(x$results_table)," features selected in both studies.\n", sep = "")
  if (x$inputs$variant == "adaptive")
  {
    cat("\nEstimates for fraction of nulls among the selected:\n", sep = "")
    cat(x$pi1," in study 1.\n", sep = "")
    cat(x$pi2," in study 2.\n", sep = "")
  }
  cat("\nFeatures selected in both studies:\n", sep = "")
  print.data.frame(x$results_table, row.names = F, digits = digits_df,...)
  if (x$inputs$variant != "non-adaptive")
    cat("\n", nrow(x$results_table)," features are significant for",
        if (x$inputs$directional_rep_claim) " directional" else "",
        " replicability claims (alpha = ",x$inputs$alpha,").\n", sep = "")
  invisible(x)
}

#' @title Check Whether a Vector Has All\\Some Names
#' @description NAs in the values are allowed.
#' @examples
#' radjust:::has_some_names(c(a = 1, b = 2, 3))
#' radjust:::has_some_names(c(a = 1, b = 2, NA))
#' radjust:::has_all_names(c(a = 1, b = 2, 3))
#' radjust:::has_all_names(c(a = 1, b = 2, NA))
#' @keywords internal
has_some_names <- function(x) {!is.null(names(x))}

#' @rdname has_some_names
#' @keywords internal
has_all_names <- function(x) {!is.null(names(x)) & all(names(x)!="")}

