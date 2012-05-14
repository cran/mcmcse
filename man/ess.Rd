\name{ess}
\alias{ess}
\title{Estimate effective sample size (ESS) as described in Kass et al. (1998) and Robert and Casella (2004; p. 500).}
\usage{
  ess(x, imse = TRUE, verbose = FALSE)
}
\arguments{
  \item{x}{a vector of values from a Markov chain.}

  \item{imse}{logical. If \code{TRUE}, use an approach that
  is analogous to Geyer's initial monotone positive
  sequence estimator (IMSE), where correlations beyond a
  certain lag are removed to reduce noise.}

  \item{verbose}{logical. If \code{TRUE} and \code{imse =
  TRUE}, inform about the lag at which truncation occurs,
  and warn if the lag is probably too small.}
}
\value{
  The function returns the estimated effective sample size.
}
\description{
  Estimate effective sample size (ESS) as described in Kass
  et al. (1998) and Robert and Casella (2004; p. 500).
}
\details{
  ESS is the size of an iid sample with the same variance
  as the current sample. ESS is given by
  \deqn{\mbox{ESS}=T/\eta,}{ESS = T / \eta,} where
  \deqn{\eta=1+2\sum \mbox{all lag autocorrelations}.}{\eta
  = 1 + 2 \sum all lag autocorrelations.}
}
\references{
  Kass, R. E., Carlin, B. P., Gelman, A., and Neal, R.
  (1998) Markov chain Monte Carlo in practice: A roundtable
  discussion. \emph{The American Statistician}, \bold{52},
  93--100.

  Robert, C. P. and Casella, G. (2004) \emph{Monte Carlo
  Statistical Methods}. New York: Springer.

  Geyer, C. J. (1992) Practical Markov chain Monte Carlo.
  \emph{Statistical Science}, \bold{7}, 473--483.
}

