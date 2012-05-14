\name{mcse}
\alias{mcse}
\alias{mcse.q}
\title{
Monte Carlo standard errors in MCMC
}
\description{
\code{mcse} and \code{mcse.q} calculate standard errors in Markov chain Monte Carlo settings for expectations and quantiles, respectively.
}
\usage{
mcse(vals, bs = "sqroot", g = NULL, meth = "BM", warn = FALSE)

mcse.q(vals, qval, bs = "sqroot", g = NULL, 
	     meth = "BM", warn = FALSE)
}
\arguments{
  \item{vals}{A vector of values from a Markov chain of length \eqn{N}.}
  \item{qval}{The quantile of interest when using \code{mcse.q}.}
  \item{bs}{The batch size used for the calculation (equal to the truncation point or subsampling size).  The default value of \code{sqroot} uses \eqn{\sqrt{N}}.}
  \item{g}{A function \code{g} for \code{mcse} returns an estimate of \eqn{E(g(x))} with an estimate of the MCSE and \code{mcse.q} returns an estimate of the qth quantile of the univariate distribution function of \eqn{g(x)} with an estimate of the MCSE.  The default value of \code{g = NULL} will use the original values from \code{vals}.}
  \item{meth}{Specifies the method used in calculating the MCSE.  \code{mcse} contains options for \code{"BM"} (batch means), \code{"OBM"} (overlapping batch means), \code{"TukeyHanning"} (spectral variance methods with a Tukey-Hanning window), and \code{"Bartlett"} (spectral variance methods with a Bartlett window).  \code{mcse.q} contains options for \code{"BM"} (batch means), \code{"OBM"} (overlapping batch means), and \code{"Subsampling"} (subsampling bootstrap methods).}
  \item{warn}{Provides a warning if the chain length is less than 1000.}
}
\value{
\code{mcse} and \code{mcse.q} return a list with the requested estimate (expectation or quantile) and corresponding MCSE.
\item{\code{$est}}{Expectation or quantile estimate of interest.}
\item{\code{$se}}{Monte Carlo standard error.}
}

\references{
Flegal, J. M. (2012).  Applicability of subsampling bootstrap methods in {M}arkov chain {M}onte {C}arlo.  In Wozniakowski, H. and Plaskota, L., editors, \emph{Monte Carlo and Quasi-Monte Carlo Methods 2010} (to appear). Springer-Verlag.

Flegal, J. M. and Jones, G. L. (2010). Batch means and spectral variance estimators in {M}arkov chain {M}onte {C}arlo.  \emph{The Annals of Statistics}, 38:1034--1070.

Flegal, J. M. and Jones, G. L. (2011).  Implementing {M}arkov chain {M}onte {C}arlo: Estimating with confidence.  In Brooks, S., Gelman, A., Jones, G., and Meng, X., editors, \emph{Handbook of {M}arkov Chain {M}onte {C}arlo}, pages 175--197. Chapman & Hall/CRC Press.

Flegal, J. M., Jones, G. L., and Neath, R. (2012).  Quantile Estimation via {M}arkov chain {M}onte {C}arlo.  \emph{University of California, Riverside, Technical Report}.

Jones, G. L., Haran, M., Caffo, B. S., and Neath, R. (2006).  Fixed-width output analysis for {M}arkov chain {M}onte {C}arlo.  \emph{Journal of the American Statistical Association}, 101:1537--1547.
}
\author{
James M. Flegal <jflegal@ucr.edu>
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
