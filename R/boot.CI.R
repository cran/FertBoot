#'
#' Bootstrap confidence intervals of mean
#' @param x a vector of observation
#' @param alpha significance level (default: 0.05)
#' @param CI.type type of CI required (default: "all")
#'
#' @return \code{boot.CI} return list of confidence intervals of mean (\code{CI.percent}: percentile, \code{CI.BC}: bias-corrected and \code{CI.BCa}: bias-corrected and accelerated).
#' @import stats nls.multstart simpleboot
#' @examples
#'
#' set.seed(12)
#' boot.CI(rnorm(1000, mean=0, sd=1), alpha=0.05, CI.type="per") # example of wrong input for type
#' boot.CI(rnorm(1000, mean=0, sd=1), alpha=0.05, CI.type="all") # require all type
#'
#' @export
#'
boot.CI <- function(x, alpha=0.05, CI.type="all") {
  B <- length(x)

  result <- NULL

  if ("percentile" %in% CI.type || "all" %in% CI.type) {
  # Percentile CI
  CI.percent <- quantile(x, c(alpha/2, 1 - alpha/2))

  result <- c(result, list(CI.percent=CI.percent))
  }


  if ("BC" %in% CI.type || "all" %in% CI.type) {
  # Bias-corrected CI
  z0 <- qnorm(sum(x < mean(x))/B)
  CI.BC <- quantile(x, c(pnorm(qnorm(alpha/2) + 2*z0), pnorm(qnorm(1-alpha/2) + 2*z0)))

  result <- c(result, list(CI.BC=CI.BC))
  }

  if ("BCa" %in% CI.type || "all" %in% CI.type) {
  # Accelerated bias-corrected CI
  theta <- rep(NA, B)
  for (i in 1:B) {theta[i] <- mean(x[-i])}
  theta.mean <- mean(theta)

  a <- sum((theta.mean - theta)^3)/(6*sum((theta.mean - theta)^2)^1.5)

  alpha1 <- pnorm(z0 + (z0 + qnorm(alpha/2))/(1-a*(z0 + qnorm(alpha/2))))
  alpha2 <- pnorm(z0 + (z0 + qnorm(1 - alpha/2))/(1-a*(z0 + qnorm(1 - alpha/2))))
  CI.BCa <- quantile(x, c(alpha1, alpha2))

  result <- c(result, list(CI.BCa=CI.BCa))
  }

  if (is.null(result)) {
    cat("Auto-corrected with percentile CI!\n")
    result <- boot.CI(x, alpha=0.05, CI.type="percentile")
  }

  result
}
