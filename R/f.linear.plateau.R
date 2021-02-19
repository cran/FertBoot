#'
#' Fitting linear plateau model using multiple initial values
#'
#' \code{f.linear.plateau} fits linear plateau model using multiple initial values. The multiple initial values are randomly sampled in  a "cube" of parameter space.
#'     More precisely, linear plateau model assumes y ~ a + b * (x - c) * (x <= c).
#'
#' @param d data frame with two columns (\code{x} and \code{y})
#' @param start initial estimate for non-linear least square (default value: \code{ list(a = 1, b = 1, c = 1)})
#' @param plus_minus radius of random initial values (default: \code{100})
#' @param n.start total number of initial points considered (default: \code{1000})
#' @param msg logical flag whether printing progress
#'
#' @return \code{f.linear.plateau} returns a list of two components (if converged): \code{nls.summary}: summary of the fitted model; \code{nls.model}: nls object
#' @import stats nls.multstart simpleboot
#' @examples
#'
#' set.seed(4)
#' x <- rep(1:300, each=4)
#' a <- 8; b <- 0.05; c <- 100
#' y <- a + b * (x - c) * (x <= c) +
#'     rnorm(length(x), sd=0.1)
#' d <- cbind(x,y)
#'
#' # a converged example:
#' ans <- f.linear.plateau(d, start=list(a = 7, b = 0.1, c = 150),
#'     plus_minus=10, n.start=10, msg=FALSE)
#'
#' summary(ans$nls.model)
#'
#'
#'
#' @export
#'
f.linear.plateau <- function(d, start = list(a = 1, b = 1, c = 1),
                           plus_minus = 1e2, n.start=1000, msg=FALSE) {

  if ("matrix" %in% class(d)) {d <- data.frame(d)}

  names(d) = c("x", "y")

  ml = lm(d[, 2] ~ d[, 1])
  mq = lm(d[, 2] ~ d[, 1] + I(d[, 1]^2))
  c1 = coef(ml)[[1]]
  c2 = coef(ml)[[2]]
  c3 = coef(mq)[[1]]
  c4 = coef(mq)[[2]]
  c5 = coef(mq)[[3]]
  pc = -0.5 * c4/c5
  ff = function(x) {
    c3 + c4 * x + c5 * x^2
  }
  pp = ff(pc)
  a = start[[1]]
  b = start[[2]]
  c = start[[3]]
  a11 = ifelse(a == 1, pp, a)
  b11 = ifelse(b == 1, c2, b)
  c11 = ifelse(c == 1, pc, c)

  m <-  nls.multstart::nls_multstart(y ~ a + b * (x - c) * (x <= c),
                                     start_lower=list(a = a11 - plus_minus, b=b11 - plus_minus, c=c11 - plus_minus),
                                     start_upper=list(a = a11 + plus_minus, b=b11 + plus_minus, c=c11 + plus_minus),
                                     iter=n.start, data = d, supp_errors ="Y")
  if(!is.null(m)) {

    c = coef(m)
    a = c[1]
    b = c[2]
    c = c[3]
    nls.summary = c(a, b, c, summary(m)[11][[1]][10], summary(m)[11][[1]][11],
          summary(m)[11][[1]][12], AIC(m), BIC(m),
          a, c)

    nls.summary = as.data.frame(nls.summary)
    rownames(nls.summary) = c("coefficient a", "coefficient b", "coefficient c",
                              "p-value t.test for a", "p-value t.test for b", "p-value t.test for c",
                              "AIC", "BIC",
                              "maximum or minimum value for y", "critical point in x")
    res <- list(nls.summary=nls.summary, nls.model=m)
  } else{
    res <- NULL
    if (msg) cat("Not converged! may try larger n.start!\n")
  }
  res
}
