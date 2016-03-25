\name{multiESS}
\alias{multiESS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effective Sample Size of a multivariate Markov chain.
}
\description{
Calculate the effective sample size of the Markov chain, using the multivariate dependence structure of the process.
}
\usage{
multiESS(x, covmat = NULL, g = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{an n by p matrix that represents the Markov chain output.}
\item{covmat}{optional matrix estimate obtained using \code{mcse.multi}.}
\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}
}
\details{
Effective sample size is the size of an iid sample with the same variance
  as the current sample. ESS is given by
  \deqn{\mbox{ESS}= n \frac{|\Lambda|^{1/p}}{ |\Sigma|^{1/p},}}{ESS = n |\Lambda|^{1/p}/ |\Sigma|^{1/p},}
  where \eqn{\Lambda}{\Lambda} is the sample covariance matrix for g and \eqn{\Sigma}{\Sigma} is an estimate of the Monte Carlo standard error for g.
}
\value{
  The function returns the estimated effective sample size.
}


\examples{
library(mAr)
p <- 3
n <- 1e3
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
set.seed(100)
foo <- matrix(rnorm(p^2), nrow = p)
foo <- foo \%*\% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)
  
out <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))

multiESS(out)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


\seealso{
  \code{\link{minESS}}, which calculates the minimum effective samples required for the problem.

  \code{\link{ess}} which calculates univariate effective sample size using a Markov chain and a function g.
}
