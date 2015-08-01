\name{confRegion}
\alias{confRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence regions (ellipses) from Markov chain
}
\description{
Function constructs confidence regions (ellipses) from the Markov chain output for the features of interest. Function uses the ellipse package.
}
\usage{
confRegion(mcse.obj, which = c(1,2), level = .95)
}

\arguments{
\item{mcse.obj}{the list returned by the \code{mcse.multi} command}
% \item{x}{an n by p matrix that represents the Markov chain output}
\item{which}{integer vector of length 2 indicating the component for which to make the confidence ellipse. Chooses the first twp by default.}
% \item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}
% \item{covmat}{optional matrix estimate obtained using \code{mcse.multi}}
% \item{method}{any of \code{``bm'', ``bartlett'', ``tukey''}. bm represents batch means estimator,  bartlett and tukey represents the modified-Bartlett window and the Tukey-Hanning windows for the spectral variance estimators. }
\item{level}{confidence level for the ellipse}
% \item{...}{arguments to pass to the \code{plot} function}
}
\details{
Plots the resulting confidence region
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
mcerror <- mcse.multi(out)
confRegion(mcerror)
}

