#' @title Adjust p-values for Replicability across Two Independent Studies with Multiple Endpoints
#'
#' @description  The function receives two vectors of p-values from two independent studies and adjust
#'  them to control the False Discovery Rate of replicability claims.
#'
#' @param pv1,pv2 numeric vectors of p-values from the two studies, corresponding to the selected features (the default) or all the features; can be either of the same length or named.
#' @param w1 the relative weight for the p-values from study 1.
#' @param input_type whether \code{pv1} and \code{pv2} contain all the p-values from the studies or only the selected (the default).
#' @param general_dependency whether the correct for general dependency.
#' @param selection_method ...
#' @param alpha the threshold for the selection method ...
#'
#' @return a list with the following elements:
#'   \tabular{ll}{
#'     \code{call}  \tab  the function call \cr
#'     \code{inputs}  \tab  a list with the function's input parameters (except \code{pv1} and \code{pv2}) \cr
#'     \code{results_table} \tab  a data frame with the features selected in both studies and their r-values. The data frame includes 6 columns: \cr
#'     \code{selected1, selected2}  \tab the named or indexed features selected in study 1, when \code{selection_method!="none"} \cr
#'     \code{n_selected1, n_selected2} \tab  the named or indexed features selected in study 1, when \code{selection_method!="none"} \cr
#'     \code{pi1, pi2}  \tab  the estimate of the true-nulls fraction in the study1\\study2
#'   }
#'
#'   \tabular{lll}{
#'     \code{name}         \tab char.   \tab the name of the feature as extracted from the named vectors, or the location, if the input vectors are not named. \cr
#'     \code{p.value.1}    \tab numeric \tab the one-sided p-value from study 1 as inputed. In the case of \code{directional_rep_claim==TRUE} the 'real' one-sided p-values is presented, i.e: \code{pmin(pv1,1-pv1}. \cr
#'     \code{p.value.2}    \tab numeric \tab same as in \code{p.value.1}. \cr
#'     \code{r.value}      \tab numeric \tab the replicability adjusted p-value (= r-value). \cr
#'     \code{Direction}    \tab char.   \tab the direction of the replicability claim, when \code{directional_rep_claim==TRUE}. \cr
#'     \code{Significant}  \tab char.   \tab a notion whether the replicability claim is significant, when \code{selection_method!="none"}.
#'   }
#'
#' @details   EDIT extended details about the function.
#'
#' @references Add here the reference to the article.
#'
#' @seealso \code{\link{radjust_pf}} for replicability analysis in primary and follow-up studies.
#'
#' @export
#'
#' @examples
#'
#' data(mice)
#' ## transform the two-sided p-values to one-sided in the same direction (left);
#' ## we use the direction of the test statistic to do so and assume that it is continuous
#' pv1 <- ifelse(mice$dir_is_left1, mice$twosided_pv1/2, 1-mice$twosided_pv1/2)
#' pv2 <- ifelse(mice$dir_is_left2, mice$twosided_pv2/2, 1-mice$twosided_pv2/2)
#' ## run the examples as in the article
#'
#' mice_rv_adaptive <- radjust_sym(pv1, pv2, input_type = "all", directional_rep_claim = T, selection_method = "alpha-adaptive")
#' mice_rv_adjusted <- radjust_sym(pv1, pv2, input_type = "all", directional_rep_claim = T, selection_method = "alpha-adjusted")
#'
#' mice_rv_no_sel <- radjust_sym(pv1, pv2, input_type = "selected", directional_rep_claim = T, selection_method = "none")
#' mice_rv_no_sel
#'
#' print(mice_rv_adaptive)
radjust_sym <- function(pv1,
                         pv2,
                         w1 = 0.5,
                         input_type = c("selected_features", "all_features"),
                         general_dependency = FALSE,
                         directional_rep_claim = FALSE,
                         selection_method = c("alpha-adjusted", "alpha-adaptive", "none"),
                         alpha = if (selection_method == "none") NULL else 0.05
)
{



  ## ---- validate inputs ----
  stopifnot(is.vector(pv1, mode = "numeric"), pv1 <= 1, 0 < pv1,
            is.vector(pv2, mode = "numeric"), pv2 <= 1, 0 < pv2,
            is.vector(w1, mode = "numeric"),  w1 < 1, 0 < w1, length(w1) == 1)
  stopifnot(is.logical(general_dependency),
            is.logical(directional_rep_claim))

  selection_method <- match.arg(selection_method)
  if (selection_method == "none" && !is.null(alpha))
    warning("selection_method == 'none' so alpha value is ignored")
  if (selection_method != "none" && is.null(alpha))
    stop("selection_method != 'none' but alpha value is not specified")
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
    # The following will produce NAs where there is no match, so to get lengths and sums we use sum( , na.rm = T).
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
      "Note:\tDirectional replicability claim option is TRUE,
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
  # Note that NAs are kept NAs
  if (selection_method == "none")
  {
    if (input_type == "all_features")
      stop("The selection_method cannot be 'none' when input_type equals to 'all_features'")
    s1  <- rep(TRUE,length(pv1)) # "select" all
    s2  <- rep(TRUE,length(pv2))
    R1 <- sum(s1, na.rm = T)
    R2 <- sum(s2, na.rm = T)
    pi1 <- 1
    pi2 <- 1
  } else
   if (selection_method == "alpha-adjusted")
    {
      s1 <- (pv1real <= w1 * alpha)
      s2 <- (pv2real <= (1 - w1) * alpha)
      R1 <- sum(s1, na.rm = T)
      R2 <- sum(s2, na.rm = T)
      pi1 <- 1
      pi2 <- 1
    } else
      if (selection_method == "alpha-adaptive")
      {
        s1 <- (pv1real <= alpha)
        s2 <- (pv2real <= alpha)

        R1 <- sum(s1, na.rm = T)
        R2 <- sum(s2, na.rm = T)

        # the fraction of true nulls in study 1 among the selected in 2
        pi1 <- (1 + sum(pv1tag[s2] > alpha, na.rm = T)) / (R2 * (1 - alpha))
        # the fraction of true nulls in study 2 among the selected in 1
        pi2 <- (1 + sum(pv2tag[s1] > alpha, na.rm = T)) / (R1 * (1 - alpha))
      }

  if (directional_rep_claim)
    s12 <- (s1 & s2 & rep_claim_direction != "Opposite")
    # TODO: Should we not present "Opposite" direction in the intersection table?
  else
    s12 <- (s1 & s2)

   ## ---- estimate pi1 and pi2 (when alpha-adaptive method is selected) ----

  if (general_dependency)
  {
    d1 <- sum(1 / (1:R2), na.rm = T) # <-- verify this
    d2 <- sum(1 / (1:R1), na.rm = T) # <-- verify this
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

  if (!(selection_method == "none")) tbl_s12$`Significant` <- ifelse(rv <= alpha,"*","")
  # TODO: consider to replace with TRUE\FALSE (and use symnum() in th print method).

  ## ---- generate the output ----

  cl <- match.call()
  inputs <- list(
    w1 = w1,
    input_type = input_type,
    general_dependency = general_dependency,
    directional_rep_claim = directional_rep_claim,
    selection_method = selection_method,
    alpha = alpha)

  output <- list(
    call = cl,
    inputs = inputs,
    results_table = tbl_s12,
    selected1 = pv1real[s1],
    selected2 = pv2real[s2],
    n_selected1 = R1,
    n_selected2 = R2,
    pi1 = pi1,
    pi2 = pi2
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
  cat("Selection (", if (x$inputs$selection_method == "none") "as-inputed"  else x$inputs$selection_method,"):\n", sep = "")
  cat(x$n_selected1," features selected in study 1.\n", sep = "")
  cat(x$n_selected2," features selected in study 2.\n", sep = "")
  cat(nrow(x$results_table)," features selected in both studies.\n", sep = "")
  if (x$inputs$selection_method == "alpha-adaptive")
  {
    cat("\nEstimates for fraction of nulls among the selected:\n", sep = "")
    cat(x$pi1," in study 1.\n", sep = "")
    cat(x$pi2," in study 2.\n", sep = "")
    cat("\nFeatures selected in both studies:\n", sep = "")
  }
  print.data.frame(x$results_table, row.names = F, digits = digits_df,...)
  if (x$inputs$selection_method != "none")
    cat("\n", nrow(x$results_table)," features are significant for",
        if (x$inputs$directional_rep_claim) " directional" else "",
        " replicability claims (alpha = ",x$inputs$alpha,").\n", sep = "")
  invisible(x)
}

#' @title Check Whether a Vector Has All\\Some Names
#' @description NAs in the values are allowed.
#' @examples
#' has_some_names(c(a = 1, b = 2, 3))
#' has_some_names(c(a = 1, b = 2, NA))
#' has_all_names(c(a = 1, b = 2, 3))
#' has_all_names(c(a = 1, b = 2, NA))
#' @keywords internal
has_some_names <- function(x) {!is.null(names(x))}

#' @rdname has_some_names
#' @keywords internal
has_all_names <- function(x) {!is.null(names(x)) & all(names(x)!="")}

