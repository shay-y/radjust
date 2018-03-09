#' @title Adjust p-values for Replicability across (Independent) Primary and Follow-up Studies with Multiple Endpoints
#'
#' @description The function computes r-values given two vectors of p-values from primary and
#'  follow-up studies. The r-values assess the False Discovery Rate (FDR) of repilcability
#'  claims across the primary and follow-up studies.

#' @param pv1,pv2 numeric vectors, of the same length, of the p-values from the from the primary study (\code{pv1}) with
#'  the corresponding p-values from the follow-up study (\code{pv2}).
#' @param m     the number of features examined in the primary study (can be bigger than \code{length(pv1)}).
#' @param c2    the relative boost to the p-values from the \strong{primary} study.
#'  \code{c2 = 0.5} (the default) is recommended. It was observed in simulations to yield
#'  similar power to procedure with the optimal value (which is unknown for real data).
#' @param l00   a lower bound of the fraction of features (out of m) with true null hypotheses in both studies.
#'   For example, for GWAS on the whole genome, the choice of 0.8 is conservative
#'   in typical applications.
#' @param variation
#'  \describe{
#'     \item{none}{the default.}
#'     \item{use_m_star}{use \eqn{m^*=m\sum_{i=1}^{m}\frac{1}{i}}{m*=m*sum(1/i)} modification of \code{m}.}
#'     \item{use_threshold}{c1 is computed given the threshold \code{tt}.}
#'  }
#'  Both variations guarantee that the procedure that decleares all r-values below \code{alpha} as replicability claims,
#'  controls the FDR at level \code{alpha}, for any type of dependency of the p-values in the primary study.
#' @param tt  the selection rule threshold for p-values from the primary study; must be supplied when
#'  variation 'use_threshold' is selected, otherwise ignored.
#'
#' @return vector of length of \code{pv2} and \code{pv2}, containing the r-values.
#'
#' @details EDIT  extended details about the function.
#'  
#' @note The function is also available as a web applet:  \url{http://www.math.tau.ac.il/~ruheller/App.html}
#'
#'
#' @examples
#'  data(crohn)
#'  rv  <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547, l00 = 0.8)
#'  rv2 <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547, l00 = 0.8, variation="use_threshold",tt = 1e-5)
#'
#' @seealso \code{\link{radjust_sym}} for replicability analysis in two symmetric design (EDIT)
#' @export

radjust_pf <- function (pv1, pv2, m, c2 = 0.5, l00= 0, variation = c("none","use_m_star","use_threshold"), tt = NULL, alpha = 0.05)
{
  variation <- match.arg(variation)

  if (variation != "use_threshold" & !is.null(tt))
    warning("threshold tt is ignored")
  if (variation == "use_threshold" & !is.null(tt))
  {
    if (tt <= (1-c2)/(1-l00*(1-c2*alpha)) * alpha/m )
    {
      warning("since t < c(q)q/m, no modification to the original r-value computation was necessary (see section Derivation & Properties in the article)")
      variation <- "none"
    }
    if (tt >= (1-c2)/(1-l00*(1-c2*alpha))*alpha/(1+sum(1/(1:(m-1)))))
      warning("for the selected threshold t, the 'use_threshold' variation won't lead\nto more discoveries than the 'use_m_star' variation.")
  }
  if (variation == "use_m_star") m <- m*sum(1L/(1L:m))
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
  if (variation == "use_threshold" & is.null(tt))
    stop("specify threshold tt")

  # ---- function definition: compute r-value of given value of x: ----
  radjust_pf_x <- function (x) {
    c1 <- switch(variation,
                 none       = (1-c2)/(1-l00*(1-c2*x)),
                 use_m_star = (1-c2)/(1-l00*(1-c2*x)),
                 use_threshold = {
                   if (tt <= (1-c2)/(1-l00*(1-c2*x)) * x/m)
                     (1-c2)/(1-l00*(1-c2*x))
                   else
                   {
                     lower <- 1e-6 ; upper <- (1-c2)/(1-l00*(1-c2*x))
                     f <- function(a)
                     {
                       if (tt*m/(a*x) < 10)
                       {
                         a*(1+sum(1/(1:(ceiling(tt*m/(a*x)-1))))) - (1-c2)/(1-l00*(1-c2*x))
                       }
                       else
                       {
                         a*(1-digamma(1)+1/(2*ceiling(tt*m/(a*x)-1))+log(ceiling(tt*m/(a*x)-1))) - (1-c2)/(1-l00*(1-c2*x))
                       }
                     }

                     c1_sol1 <- uniroot(f,c(lower,upper))
                     next_step <- tt*m/(ceiling(tt*m/(c1_sol1$root*x))*x)
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
