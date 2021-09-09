## ----var----------------------------------------------------------------------
library(mcmcse)
mu = c(2, 50)
sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)

# Monte Carlo sample size is N
N <- 5e3
set.seed(100)
chain <- BVN_Gibbs(n = N, mu = mu, sigma = sigma)

## ----foo, echo = FALSE--------------------------------------------------------
  colnames(chain) <- c("Y1", "Y2")

## ----output-------------------------------------------------------------------
#Rows has observations (samples) and each comlumn is a component. 
head(chain)

## ----means--------------------------------------------------------------------
colMeans(chain)

## ----g------------------------------------------------------------------------
g <- function(x)
{
  return(sum(x^2))
}

## ----g_est--------------------------------------------------------------------
# Apply the function g to each row
gofy <- apply(chain, 1, g)

# Monte Carlo estimate
mean(gofy)

## ----mcse---------------------------------------------------------------------
# Batch means estimator
mcerror_bm <- mcse.multi(x = chain, method =  "bm", r = 1,
                         size = NULL, g = NULL, adjust = TRUE, 
                         blather = TRUE)

# Overlapping batch means estimator
mcerror_obm <- mcse.multi(x = chain, method =  "obm", r = 1,
                         size = NULL, g = NULL, adjust = TRUE, 
                         blather = TRUE)

# Spectral variance estimator with Bartlett window
mcerror_bart <- mcse.multi(x = chain, method =  "bartlett", r = 1,
                           size = NULL, g = NULL, adjust = TRUE, 
                           blather = TRUE)

# Spectral variance estimator with Tukey window
mcerror_tuk <- mcse.multi(x = chain, method =  "tukey", r = 1,
                          size = NULL, g = NULL, adjust = TRUE, 
                          blather = TRUE)

# Initial sequence estimator, unadjusted
mcerror_is <- mcse.initseq(x = chain, g = NULL, 
                           adjust = FALSE, blather = TRUE)

# Initial sequence estimator, adjusted
mcerror_isadj <- mcse.initseq(x = chain, g = NULL, 
                              adjust = TRUE, blather = TRUE)

## ----uni----------------------------------------------------------------------
mcse(x = chain[,1], method = "bm", g = NULL)
mcse.mat(x = chain, method = "bm", g = NULL)

## ----sigma_g------------------------------------------------------------------
g
mcerror_g_bm <- mcse.multi(x = chain, g = g, blather = TRUE)
mcerror_g_is <- mcse.initseq(x = chain, g = g, blather = TRUE)

mcerror_g_bm$cov

# Initial Sequence error is larger than batch means, as expected.
mcerror_g_is$cov

# Returned value is asymptotic variance. 
# So we calculate the standard error here.
sqrt(mcerror_g_bm$cov/N) 
sqrt(mcerror_g_is$cov/N)

## ----confRegion, out.height = '8cm'-------------------------------------------
plot(confRegion(mcerror_bm, which = c(1,2), level = .90), type = 'l', asp = 1)
lines(confRegion(mcerror_bart, which = c(1,2), level = .90), col = "red")

## ----comp_region, out.height = '8cm'------------------------------------------
plot(confRegion(mcerror_is, which = c(1,2), level = .95), type = 'l', asp = 1)
lines(confRegion(mcerror_is, which = c(1,2), level = .90), col = "red")

## ----minESS-------------------------------------------------------------------
# For mu
minESS(p = 2, alpha = .05, eps = .05)

#For mu_g
minESS(p = 1, alpha = .05, eps = .05)

## ----eps----------------------------------------------------------------------
# For mu
minESS(p = 2, alpha = .05, ess = 1000)

#For mu_g
minESS(p = 1, alpha = .05, ess = 1000)

## ----multiess-----------------------------------------------------------------
multiESS(chain)

# Using spectral variance estimators
multiESS(chain, covmat = mcerror_bart$cov)

# Using initial sequence estimators
# Since this is a conservative estimator, ess will be smaller
multiESS(chain, covmat = mcerror_is$cov)

## ----moresamples--------------------------------------------------------------
set.seed(100)
chain <- BVN_Gibbs(1e4, mu, sigma)

# larger than 7529
multiESS(chain)

# larger than 7529
multiESS(chain, covmat = mcerror_bart$cov)

# larger than 7529
multiESS(chain, covmat = mcerror_is$cov)

## ----estvssamp, out.width = '8cm'---------------------------------------------
estvssamp(chain[,1])

