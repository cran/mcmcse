\name{qqTest}
\alias{qqTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QQplot for Markov chains 
}
\description{
QQplot for Markov chains using a consistent estimate of the Markov Chain CLT covariance matrix.
}
\usage{
qqTest(x, covmat, g = NULL)
}
\arguments{
\item{x}{an n by p matrix that represents the Markov chain output}
\item{covmat}{a p by p matrix estimate obtained through the \code{mcse.multi} function}
\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}
}
\examples{
library(mAr)
p <- 35
n <- 1e4
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
set.seed(100)
foo <- matrix(rnorm(p^2), nrow = p)
foo <- foo \%*\% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)
  
out <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))

mcse.bm <- mcse.multi(x = out)
qqTest(out, covmat = mcse.bm$cov)
}
