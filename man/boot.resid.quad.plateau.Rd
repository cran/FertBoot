% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot.resid.quad.plateau.R
\name{boot.resid.quad.plateau}
\alias{boot.resid.quad.plateau}
\title{Quadratic plateau model estimation by bootstrapping residuals}
\usage{
boot.resid.quad.plateau(
  mod,
  data,
  x.range = data.frame(x = seq(0, 280, by = 40)),
  B = 100 - 1,
  plus_minus = 100,
  n.start = 5000,
  print.progress = TRUE
)
}
\arguments{
\item{mod}{a full model list, probably from \code{f.quad.plateau()}}

\item{data}{data frame with two columns (\code{x} and \code{y})}

\item{x.range}{vector of data.frame with one column for range of N rate of interested for prediction interval}

\item{B}{bootstrap sample size}

\item{plus_minus}{radius of random initial values (default: \code{100})}

\item{n.start}{total number of initial points considered (default: \code{1000})}

\item{print.progress}{logical flag whether printing progress}
}
\value{
\code{boot.resid.quad.plateau} returns a list of two elements:
\code{result}: matrix with B rows and columns containing bootstrap sample for parameter (\code{a,b,c}), optimal N and yield (\code{max_x, max_y}), log-likelihood (\code{logLik}) and N values of interest;
\code{x.range}: range of x considered for prediction interval (same as \code{x.range} in vector form)
}
\description{
\code{boot.resid.quad.plateau} is the core function to implement bootstrapping residuals on quadratic plateau models, which assumes
    y = (a + b * x + c *x^2) * (x <= -0.5*b/c) + (a + -b^2/(4 * c)) * (x > -0.5 * b/c). Note that this function may take minutes up to days. Parallel computing may be necessary.
    We suggest users start with a smaller \code{B} and moderate n.start to see if the bootstrap models can converge.
    In general, increasing n.start and plus_minus may help with ease of convergence. For rigorous statistical inference, B should be on the order of a thousand.
}
\examples{

\donttest{
set.seed(1)
x <- rep(1:300, each=5)
a <- 8; b <- 0.05; c <- -1e-4
y <- (a + b * x + c *x^2) * (x <= -0.5*b/c) + (a + -b^2/(4 * c)) * (x > -0.5 * b/c) +
    rnorm(length(x), sd=0.1)
d <- cbind(x,y)

ans <- f.quad.plateau(d, start=list(a = 7, b = 0.02, c = 1e-5),
    plus_minus=10, n.start=10, msg=FALSE)


boot.resid.quad.plateau(ans, d, x.range=seq(0,280,by=40),
    B=1e1-1, plus_minus = 1e2, n.start=1000, print.progress=TRUE) # use larger B for inference

}



}
