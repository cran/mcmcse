\name{mcmcse-package}
\alias{mcmcse-package}
\alias{mcmcse}
\docType{package}
\title{
Monte Carlo standard errors for MCMC
}
\description{
Tools for calculating Monte Carlo standard errors (MCSE) in Markov chain Monte Carlo (MCMC) settings.  This package considers standard errors for expectation and quantile estimators.   
}
\details{
\tabular{ll}{
Package: \tab mcse\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-05-14\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
James M. Flegal <jflegal@ucr.edu>

Maintainer: James M. Flegal <jflegal@ucr.edu>
}
\references{
Flegal, J. M. (2012).  Applicability of subsampling bootstrap methods in {M}arkov chain {M}onte {C}arlo.  In Wozniakowski, H. and Plaskota, L., editors, \emph{Monte Carlo and Quasi-Monte Carlo Methods 2010} (to appear). Springer-Verlag.

Flegal, J. M. and Jones, G. L. (2010). Batch means and spectral variance estimators in {M}arkov chain {M}onte {C}arlo.  \emph{The Annals of Statistics}, 38:1034--1070.

Flegal, J. M. and Jones, G. L. (2011).  Implementing {M}arkov chain {M}onte {C}arlo: Estimating with confidence.  In Brooks, S., Gelman, A., Jones, G., and Meng, X., editors, \emph{Handbook of {M}arkov Chain {M}onte {C}arlo}, pages 175--197. Chapman & Hall/CRC Press.

Flegal, J. M., Jones, G. L., and Neath, R. (2012).  Quantile Estimation via {M}arkov chain {M}onte {C}arlo.  \emph{University of California, Riverside, Technical Report}.

Jones, G. L., Haran, M., Caffo, B. S., and Neath, R. (2006).  Fixed-width output analysis for {M}arkov chain {M}onte {C}arlo.  \emph{Journal of the American Statistical Association}, 101:1537--1547.

}
\examples{
# Creates 10000 iterations in an AR(1) Markov chain with rho = .9
p <- 10000
tau <- 1
mc <- double(p)
mc[1] <- 2
for(i in 1:(p-1)){
mc[(i+1)] <- .9 * mc[i] + rnorm(1, 0, tau)
}

# Estimates the mean, .1 quantile and .9 quantile with MCSE using BM
mcse(mc)
mcse.q(mc, .1)
mcse.q(mc, .9)

# Estimates the mean, .1 quantile and .9 quantile with MCSE using OBM
mcse(mc, meth="OBM")
mcse.q(mc, .1, meth="OBM")
mcse.q(mc, .9, meth="OBM")

# Estimates E(x^2) with MCSE using spectral methods
g.fun <- function(x){x^2}
mcse(mc, g=g.fun, meth="TukeyHanning")

}
