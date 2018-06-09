#' @title Adjust p-values for Replicability across Two Independent, Primary and Follow-up, Studies with Multiple Endpoints
#'
#' @description Given two vectors of p-values from the primary and follow-up studies, returns the adjusted p-values for false
#' discovery rate control on replicability claims. The p-value vectors are only for features selected for follow-up.
#'
#' @param pv1 numeric vector of p-values from the primary study  which
#'  correspond to the p-values from the follow-up study (\code{pv2}).
#' @param pv2 numeric vector of p-values from the follow-up study.
#' @param m     the number of features examined in the primary study (> \code{length(pv1)}).
#' @param c2    the relative boost to the p-values from the \strong{follow-up} study.
#'  \code{c2 = 0.5} (the default) is recommended. It was observed in simulations to yield
#'  similar power to procedure with the optimal value (which is unknown for real data).
#' @param l00   a lower bound of the fraction of features (out of m) with true null hypotheses in both studies.
#'   For example, for GWAS on the whole genome, the choice of 0.8 is conservative in typical applications.
#' @param variant
#'  \describe{
#'     \item{none}{the default.}
#'     \item{general_dependency}{use \eqn{m^*=m\sum_{i=1}^{m}\frac{1}{i}}{m*=m*sum(1/i)} instead of \code{m}.}
#'     \item{use_threshold}{c1 is computed given the \code{threshold} value.}
#'  }
#'  Both variants guarantee that the procedure that declares as replicated all features with r-values below \code{alpha},
#'  controls the FDR at level \code{alpha}, for any type of dependency of the p-values in the primary study.
#' @param threshold  the selection threshold for p-values from the primary study; must be supplied when
#'  variant 'use_threshold' is selected, otherwise ignored.
#' @param alpha The FDR level to control.
#'
#' @return vector of length of \code{pv1} and \code{pv2}, containing the r-values.
#'
#' @details When many hypotheses are simultaneously examined in a primary study, and then a subset of hypotheses
#'  are forwarded for follow-up in an independent study, it is of interest to know which findings are replicated across studies.
#'  As a measure of replicability of significance, we compute the r-value, i.e. the FDR adjusted replicability p-value, for each hypothesis followed-up.
#'  This measure is different than the FDR adjusted p-value in a typical meta-analysis, where a p-value close to zero in one of
#'  the studies is enough to declare the finding as highly significant. The FDR r-value for a feature is the smallest FDR level
#'  at which we can say that the finding is among the replicated ones.
#'
#' @note The function is also available as a web applet:  \url{http://www.math.tau.ac.il/~ruheller/App.html}
#'
#' @references Bogomolov, M. and Heller, R. (2013). Discovering findings that replicate from a primary study of high dimension to a follow-up study.
#' Journal of the American Statistical Association, Vol. 108, No. 504, Pp. 1480-1492.
#'
#' @examples
#'  data(crohn)
#'  rv  <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547, l00 = 0.8)
#'  rv2 <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547, l00 = 0.8,
#'                    variant="use_threshold",threshold = 1e-5)
#'
#' @seealso \code{\link{radjust_sym}} for replicability analysis in two symmetric studies.
#' @importFrom stats uniroot
#' @export

radjust_pf <- function (pv1, pv2, m, c2 = 0.5, l00= 0, variant = c("none","general_dependency","use_threshold"), threshold = NULL, alpha = 0.05)
{
  variant <- match.arg(variant)

  if (variant != "use_threshold" & !is.null(threshold))
    warning("threshold value is ignored")
  if (variant == "use_threshold" & !is.null(threshold))
  {
    if (threshold <= (1-c2)/(1-l00*(1-c2*alpha)) * alpha/m )
    {
      warning("since t < c(q)q/m, no modification to the original r-value computation was necessary (see section Derivation & Properties in the article)")
      variant <- "none"
    }
    if (threshold >= (1-c2)/(1-l00*(1-c2*alpha))*alpha/(1+sum(1/(1:(m-1)))))
      warning("for the selected threshold t, the 'use_threshold' variant won't lead\nto more discoveries than the 'general_dependency' variant.")
  }
  if (variant == "general_dependency") m <- m*sum(1L/(1L:m))

  k <- length(pv1)

  # ---- input validations: ----
  if (k==0 | length(pv2)==0)
    stop("p-value vectors cannot have length zero")
  if (k!=length(pv2))
    stop("p-value vectors must be of equal lengths")
  if ((1<=c2)|(c2<=0))
    stop("c2 value should be in the interval (0,1)")
  if ((1<=l00)|(l00<0))
    stop("l00 value should be in the interval [0,1)")
  if (m<k)
    stop("Number of features in the primary stage (m) must be at least as large as the number of features followed-up.")
  if (any(is.na(c(pv1,pv2))))
    stop("NA's are not allowed as p-values")
  if (any(c(pv1,pv2)>1) | any(c(pv1,pv2)<=0))
    stop("p-values must be in the interval (0,1]")
  if (variant == "use_threshold" & is.null(threshold))
    stop("specify threshold value")

  # ---- function definition: compute r-value of given value of x: ----
  radjust_pf_x <- function (x) {
    c1 <- switch(variant,
                 none       = (1-c2)/(1-l00*(1-c2*x)),
                 general_dependency = (1-c2)/(1-l00*(1-c2*x)),
                 use_threshold = {
                   if (threshold <= (1-c2)/(1-l00*(1-c2*x)) * x/m)
                     (1-c2)/(1-l00*(1-c2*x))
                   else
                   {
                     lower <- 1e-6 ; upper <- (1-c2)/(1-l00*(1-c2*x))
                     f <- function(a)
                     {
                       if (threshold*m/(a*x) < 10)
                       {
                         a*(1+sum(1/(1:(ceiling(threshold*m/(a*x)-1))))) - (1-c2)/(1-l00*(1-c2*x))
                       }
                       else
                       {
                         a*(1-digamma(1)+1/(2*ceiling(threshold*m/(a*x)-1))+log(ceiling(threshold*m/(a*x)-1))) - (1-c2)/(1-l00*(1-c2*x))
                       }
                     }

                     c1_sol1 <- uniroot(f,c(lower,upper))
                     next_step <- threshold*m/(ceiling(threshold*m/(c1_sol1$root*x))*x)
                     if (f(next_step)<=0)
                     {
                       c1_sol2 <- uniroot(f,c(next_step,upper))
                       c1_sol2$root
                     }
                     else
                     {
                       c1_sol1$root
                     }
                   }
                 })

    E   <- pmax(m*pv1/c1, k*pv2/c2)
    oe  <- order(E, decreasing =TRUE)
    oer <- order(oe)
    r   <- cummin((E/rank(E, ties.method= "max"))[oe])[oer]
    r   <- pmin(r,1)
    return(r)
  }

  rv <- rep(NA,length(pv2))
  tol = min(c(pv1,pv2)[c(pv1,pv2)!=0], 0.0001)
  for (i in 1:length(pv2))
  {
    aux <- function(x) return(radjust_pf_x(x)[i]-x)
    if (aux(1)>=0) rv[i] <- 1
    else
    {
      sol <- uniroot(aux,c(tol,1),tol = tol)
      rv[i] <- sol$root
    }
  }
  return(rv)
}
