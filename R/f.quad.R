#' Fitting quadratic model using mutiple initial vaues
#'
#' \code{f.quad} fits quadratic model using multiple initial values. The multiple initial values are randomly sampled in  a "cube" of parameter space. More precisely, quadratic model assumes
#'     y =  a+b*x+c*x^2,
#'
#' @param d data drame with two columns (\code{x} and \code{y})
#' @param start initial estimate for non-linear least square (default value: \code{ list(a = 1, b = 1, c = 1)})
#' @param plus_minus radius of random initial values (default: \code{100})
#' @param n.start total number of initial points considered (deafult: \code{1000})
#' @param msg logical flag whehter printing progress
#'
#' @return \code{f.quad} returns a list of two components (if converged): \code{nls.summary}: summary of the fitted model; \code{nls.model}: nls object
#' @import stats nls.multstart simpleboot
#' @examples
#'
#' set.seed(1)
#' x <- rep(1:300, each=2)
#' a <- 8; b <- 0.05; c <- -1e-3
#' y <- a + b*x + c*x^2 + rnorm(length(x), sd=0.1)
#' d <- cbind(x,y)
#'
#' # a converged example:
#' ans <- f.quad(d, start=list(a = 7, b = 0.02, c = 1e-5),
#'     plus_minus=10, n.start=10, msg=FALSE)
#'
#' summary(ans$nls.model)
#'
#'
#'
#' @export
#'
f.quad <- function(d, start = list(a=1, b =1, c=1), plus_minus = 1, n.start=10, msg=FALSE) {

  if ("matrix" %in% class(d)) {d <- data.frame(d)}

  names(d) = c("x", "y")
  mq = lm(y ~ x + I(x^2), data=d)
  c3 = coef(mq)[1]
  c4 = coef(mq)[2]
  c5 = coef(mq)[3]

  a = start[[1]]
  b = start[[2]]
  c = start[[3]]
  a11 = ifelse(a == 1, c3, a)
  b11 = ifelse(b == 1, c4, b)
  c11 = ifelse(c == 1, c5, c)


  m <-  nls.multstart::nls_multstart(y ~ a+b*x+c*x^2,
                      start_lower=list(a = a11 - plus_minus, b=b11 - plus_minus, c=c11 - plus_minus),
                      start_upper=list(a = a11 + plus_minus, b=b11 + plus_minus, c=c11 + plus_minus),
                      iter=n.start, data = d, supp_errors ="Y")
  if(!is.null(m)) {
    c = coef(m)
    a = c[1]
    b = c[2]
    c = c[3]
    pmm = a - (b^2)/(4 * c)
    pcc = -0.5 * b/c
    nls.summary = c(a, b, c, summary(m)[[11]][1,4], summary(m)[[11]][2,4],
                    summary(m)[[11]][3,4], AIC(m), BIC(m),
                    pmm, pcc)
    # nls.summary = round(nls.summary, 8)
    nls.summary = as.data.frame(nls.summary)
    rownames(nls.summary) = c("coefficient a", "coefficient b", "coefficient c",
                              "p-value t.test for a", "p-value t.test for b", "p-value t.test for c",
                              "AIC", "BIC",
                              "maximum or minimum value for y", "critical point in x")
    res <- list(nls.summary=nls.summary, nls.model=m)
  } else{
    res <- NULL
    if(msg) cat("Not converged! Should try different start!\n")
  }
  res
}

