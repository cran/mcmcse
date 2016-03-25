\name{minESS}
\alias{minESS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum effective sample size.
}
\description{
Calculate the minimum effective samples required for appropriate accuracy.
}
\usage{
minESS(p, alpha = .05, eps = .05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{p}{dimension of the estimation problem.}
\item{alpha}{confidence level}
\item{eps}{tolerance level}
}
\details{
The minimum effective samples required when estimating a vector of length p, with 100(1-\eqn{\alpha})\% confidence and tolerance of \eqn{\epsilon} is

\deqn{\mbox{mESS} \geq \frac{2^{2/p} \pi}{(p \Gamma(p/2))^{2/p}} \frac{\chi^2_{1-\alpha, p}}{\epsilon^2} }

}
\value{
  The function returns the estimated effective sample size.
}


\examples{
minESS(p = 5)
}

\references{
Gong, L., and Flegal, J. M. A practical sequential stopping rule for high-dimensional Markov chain Monte Carlo. Journal of Computational and Graphical Statistics (to appear).

Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte Carlo, \emph{arXiv preprint arXiv:1512.07713 (2015)}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\seealso{
  \code{\link{multiESS}}, which calculates multivariate effective sample size using a Markov chain and a function g.

  \code{\link{ess}} which calculates univariate effective sample size using a Markov chain and a function g.
}