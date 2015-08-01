\name{mcse.multi}
\alias{mcse.multi}
\docType{package}
\title{
Multivariate Monte Carlo Standard Errors Estimation.
}
\description{
Function returns the estimate of the covariance matrix in the Markov Chain CLT using batch means or spectral variance methods (with different lag windows). The function also returns the volume of the resulting ellipsoidal confidence regions. 
}
\usage{
mcse.multi(x, method = "bm", size = "sqroot", g = NULL, level = 0.95, large = FALSE)	
}
\arguments{
\item{x}{an n by p matrix that of the Markov chain output}

\item{method}{any of \code{``bm'', ``bartlett'', ``tukey''}. \code{``bm''} represents batch means estimator,  \code{``bartlett''} and \code{``tukey''} represents the modified-Bartlett window and the Tukey-Hanning windows for the spectral variance estimators. }

\item{size}{can take character values of \code{``sqroot''} and \code{``cuberoot''} or any numeric value between 1 and n. Size represents the batch size in bm and the truncation point in bartlett and tukey. \code{sqroot} means size is floor(n^(1/2) and cuberoot means size is floor(n^(1/3)).}

\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}

\item{level}{confidence level of the confidence ellipsoid.}

\item{large}{if TRUE, returns the volume of the large sample confidence region using a chi square critical value.}
}
\value{
A list is returned with the following components,
\item{cov}{a p by p covariance matrix estimate.}
\item{vol}{volume of the confidence ellipsoid to the pth root.}
\item{est}{estimate of g(x).}
\item{nsim}{number of rows of the input \code{x}.}
\item{method}{method used to calculate matrix \code{cov}.}
\item{large}{logical of whether a large sample confidence region volume was calculated.}
\item{size}{value of size used to calculate \code{cov}.}
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

mcse.bm <- mcse.multi(x = out)
mcse.tuk <- mcse.multi(x = out, method = "tukey")

# If we are only estimating the mean of the first component, 
# and the second moment of the second component

g <- function(x) return(c(x[1], x[2]^2))
mcse <- mcse.multi(x = out, g = g)
}
