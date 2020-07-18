#' Quadratic plateau model estimation by bootstrapping residuals
#'
#' \code{boot.resid.quad.plateau} is the core function that implement bootstrapping residuals on quadratic plateau model, which assumes
#'     y = (a + b * x + c *x^2) * (x <= -0.5*b/c) + (a + -b^2/(4 * c)) * (x > -0.5 * b/c). Note that this functions may take minutes up to days. Parallel computing could be considered when necessary. We suggest users start with a smaller \code{B} and moderate \code{n.start} to see if the bootstrap models can convergence.  In general, increasing \code{n.start} and \code{plus_minus} may help convergence. For rigorous statistical inference, \code{B} should be reach order of thousand.
#'
#'
#' @param mod a full model list, probably from \code{f.quad.plateau()}
#' @param data data drame with two columns (\code{x} and \code{y})
#' @param x.range vector of data.frame with one column for range of N rate of interested for prediction interval
#' @param B bootstrap sample size
#' @param plus_minus radius of random initial values (default: \code{100})
#' @param n.start total number of initial points considered (deafult: \code{1000})
#' @param print.progress logical flag whehter printing progress
#'
#' @return \code{boot.resid.quad.plateau} returns a list of two elements:
#' \code{result}: matrix with B rows and columns containing bootstrap sample for parameter (\code{a,b,c}), optimal N and yield (\code{max_x, max_y}), log-likelihood (\code{logLik}) and N values of interest;
#' \code{x.range}: range of x considered for prediction interval (same as \code{x.range} in vector form)
#'
#' @import stats nls.multstart simpleboot
#'
#' @examples
#'
#'\donttest{
#' set.seed(1)
#' x <- rep(1:300, each=5)
#' a <- 8; b <- 0.05; c <- -1e-4
#' y <- (a + b * x + c *x^2) * (x <= -0.5*b/c) + (a + -b^2/(4 * c)) * (x > -0.5 * b/c) +
#'     rnorm(length(x), sd=0.1)
#' d <- cbind(x,y)
#'
#' ans <- f.quad.plateau(d, start=list(a = 7, b = 0.02, c = 1e-5),
#'     plus_minus=10, n.start=10, msg=FALSE)
#'
#'
#' boot.resid.quad.plateau(ans, d, x.range=seq(0,280,by=40),
#'     B=1e2-1, plus_minus = 1e2, n.start=1000, print.progress=TRUE)
#'
#' }
#'
#'
#'
#' @export
#'
boot.resid.quad.plateau <- function(mod, data, x.range=data.frame(x=seq(0,280,by=40)), B=1e2-1, plus_minus = 1e2, n.start=5000, print.progress=TRUE) {

  if (class(x.range) == "numeric")  x.range <- data.frame(x=x.range)

  # a, b, c, max x value, max y value, logLik
  result <- data.frame(matrix(NA, nrow= B + 1, ncol=6 + NROW(x.range)))
  names(result) <- c("a", "b", "c", "max_x", "max_y", "logLik", paste0("x_", x.range[,1]))

  data <- data.frame(data)

  data.tmp <- data ; names(data.tmp) <- c("x", "y")


  fit.value <- fitted(mod$nls.model)
  res.value <- residuals(mod$nls.model)
  start.value <- list(a=coef(mod$nls.model)[1], b=coef(mod$nls.model)[2], c=coef(mod$nls.model)[3])

  result[1, ] <- c(mod$nls.summary[c(1:3,10,9),], logLik(mod$nls.model),
                   predict(mod$nls.model, newdata=x.range))
  i <- 1
  while (i <= B) {
    if (print.progress) cat("Bootstrap residuals:", i,"out of",B,"\n")

    # Bootstrap residual
    data.tmp$y <- as.numeric(fit.value + sample(res.value, replace=TRUE))

    m.tmp <- f.quad.plateau(d=data.tmp, start=start.value, plus_minus=plus_minus, n.start=n.start)

    if(!is.null(m.tmp$nls.model)) {
      result[i+1, ] <- c(m.tmp$nls.summary[c(1:3,10,9),], logLik(m.tmp$nls.model),
                         predict(m.tmp$nls.model, newdata=x.range))
      i <- i + 1
    }else{ if (print.progress) cat("Not converged! Retry...\n")}
  }
  list(result=result, x.range = x.range[,1])
}

