## ----var----------------------------------------------------------------------
library(mAr)
p <- 3
A <- diag(c(.1, .5, .8))
C <- diag(rep(2, 3))

# Monte Carlo sample size is N
N <- 1e5
set.seed(100)
chain <- mAr.sim(w = rep(2,p), A = A, C = C, N = N)

## ----foo, echo = FALSE--------------------------------------------------------
colnames(chain) <- c("Y1", "Y2", "Y3")

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
library(mcmcse)

# Batch means estimator
mcerror_bm <- mcse.multi(x = chain, method =  "bm", r = 1,
	size = NULL, g = NULL, adjust = TRUE, blather = TRUE)

# Overlapping batch means estimator
mcerror_bm <- mcse.multi(x = chain, method =  "obm", r = 1,
	size = NULL, g = NULL, adjust = TRUE, blather = TRUE)

# Spectral variance estimator with Bartlett window
mcerror_bart <- mcse.multi(x = chain, method =  "bartlett", r = 1,
	size = NULL, g = NULL, adjust = TRUE, blather = TRUE)

# Spectral variance estimator with Tukey window
mcerror_tuk <- mcse.multi(x = chain, method =  "tukey", r = 1,
	size = NULL, g = NULL, adjust = TRUE, blather = TRUE)

# Initial sequence estimator, unadjusted
mcerror_is <- mcse.initseq(x = chain, g = NULL, 
                           adjust = FALSE, blather = TRUE)

# Initial sequence estimator, adjusted
mcerror_isadj <- mcse.initseq(x = chain, g = NULL, 
                              adjust = TRUE, blather = TRUE)

## ----outputvalue--------------------------------------------------------------
mcerror_bm$cov

mcerror_bart$cov

mcerror_tuk$cov

mcerror_is$cov

mcerror_isadj$cov.adj

rbind(mcerror_bm$est, mcerror_bart$est, mcerror_tuk$est,
      mcerror_is$est, mcerror_isadj$est)

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
minESS(p = 3, alpha = .05, eps = .05)

#For mu_g
minESS(p = 1, alpha = .05, eps = .05)

## ----eps----------------------------------------------------------------------
# For mu
minESS(p = 3, alpha = .05, ess = 1000)

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
chain <- mAr.sim(w = rep(2,p), A = A, C = C, N = 28000)

# larger than 8123
multiESS(chain)

# larger than 8123
multiESS(chain, covmat = mcerror_bart$cov)

# larger than 8123
multiESS(chain, covmat = mcerror_is$cov)

## ----estvssamp, out.width = '8cm'---------------------------------------------
estvssamp(chain[,1])

## ----qqbig--------------------------------------------------------------------
p <- 50
A <- diag(seq(.1, .9, length = p))
C <- diag(rep(2, p))

set.seed(100)
chain <- mAr.sim(w = rep(2,p), A = A, C = C, N = 10000)

## ----qq, out.width = '8cm'----------------------------------------------------
mcerror_bm <- mcse.multi(chain, method = "bm", blather = TRUE)
mcerror_isadj <- mcse.initseq(chain, adjust = TRUE, blather = TRUE)
qqTest(mcerror_bm)
qqTest(mcerror_isadj)

