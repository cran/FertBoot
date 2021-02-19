#'
#' Two sample bootstrap test for comparing different in \code{sample1} and \code{sample2}, not necessary with same sample size
#'
#' @param sample1 first sample
#' @param sample2 second sample
#' @param fun statistic (univariate) to be compared (default: \code{mean})
#' @param R number of resamples (default: \code{1000})
#' @return \code{compare.two.sample} return a list with two components, namely,
#' \code{p.value}: two tailed p-value for the bootstrap test
#' \code{object}: a "\code{simpleboot}" object allowing further analysis using other R packages, such as \code{boot})
#'
#' @import stats nls.multstart simpleboot
#' @examples
#'
#' set.seed(1203)
#' # compare median of two expontential r.v.
#' compare.two.sample(rexp(100, rate=1), rexp(100, rate=2), fun=median, R=1e3)$p.value
#'
#' f.Q1 <- function(x) quantile(x, probs=0.25)
#' compare.two.sample(rnorm(100, mean=0), rnorm(200, mean=0.5), fun=f.Q1, R=1e3)$p.value
#'
#' @export
#'
compare.two.sample <- function(sample1, sample2, fun=mean, R=1000) {

  b <- simpleboot::two.boot(sample1, sample2, fun, R = R)

  p.value <- max(min(mean(b$t >=0), mean(b$t <= 0)),1/R)*2

  list(object=b, p.value=p.value)
}
