# library(Rcpp)
# library(RcppArmadillo)
# sourceCpp("mobmc.cpp")
# sourceCpp("mbmc.cpp")
# sourceCpp("msvec.cpp")

#####################################################################
### Function adjusts a non-positive definite estimator to be PD
#####################################################################

adjust_matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 9/10)
{
  mat.adj <- mat
  adj <- epsilon*N^(-b)
  vars <- diag(mat)
  corr <- cov2cor(mat)
  eig <- eigen(corr)
  adj.eigs <- pmax(eig$values, adj)
  mat.adj <- diag(vars^(1/2))%*% eig$vectors %*% diag(adj.eigs) %*% t(eig$vectors) %*% diag(vars^(1/2))
  return(mat.adj)
}

#####################################################################
### Function approximates each column of the output by an AR(p)
### and estimates Gamma and Sigma to calculate batch size
#####################################################################

arp_approx <- function(x)
{
  
  # Fitting a univariate AR(m) model
  ar.fit <- ar(x, aic = TRUE)

  # estimated autocovariances
  gammas <- as.numeric(acf(x, type = "covariance", lag.max = ar.fit$order, plot = FALSE)$acf)
  spec <- ar.fit$var.pred/(1-sum(ar.fit$ar))^2  #asym variance
  
  if(ar.fit$order != 0)
  {
  foo <- 0
  for(i in 1:ar.fit$order)
  {
    for(k in 1:i)
    {
      foo <- foo + ar.fit$ar[i]*k*gammas[abs(k-i)+1]
    }
  }
  Gamma <- 2*(foo + (spec - gammas[1])/2 *sum(1:ar.fit$order * ar.fit$ar)  )/(1-sum(ar.fit$ar))
  } else{
    Gamma <- 0
  }
  rtn <- cbind(Gamma, spec)
  colnames(rtn) <- c("Gamma", "Sigma")
  return(rtn)
}

#####################################################################
#####################################################################
### Functions below will be exported
#####################################################################
#####################################################################


#####################################################################
### Functions estimates the "optimal" batch size using the parametric
### method of Liu et.al
#####################################################################
batchSize <- function(x, method = "bm", g = NULL)
{

  chain <- as.matrix(x)
  if(!is.matrix(chain) && !is.data.frame(chain))
    stop("'x' must be a matrix or data frame.")

  if (is.function(g)) 
  {
    chain <- apply(chain, 1, g)

    if(is.vector(chain))
    {
      chain <- as.matrix(chain)
    }else
    {
      chain <- t(chain)
    }
  }

  n <- dim(chain)[1]
  ar_fit <- apply(chain, 2, arp_approx)^2
  coeff <- ( sum(ar_fit[1,])/sum(ar_fit[2,]) )^(1/3)

  b.const <- (3/2*n)*(method == "obm" || method == "bartlett" || method == "tukey") + (n)*(method == "bm")
  b <- b.const^(1/3) * coeff
  if(b == 0) b <- 1

  b <- floor(b)
  return(b)
}

#####################################################################
### Main function. Estimates the covariance matrix
### Recommend blather = FALSE for users and TRUE for developers
#####################################################################
mcse.multi <- function(x, method = "bm", r = 3, size = NULL, g = NULL, adjust = TRUE, blather = FALSE)
{ 
  
  # at some point the method used may be different
  # from the method asked. Blather will output this
  method.used <- method
  if(method == "lug")   # not releaved to the public. Inside option for developers
  {
    method = "bm"
    r <- 3
  }

  if(!is.numeric(r)) stop("r should be numeric")
  if(method != "bm" &&  method != "obm" && method != "bartlett" && method != "tukey")
  {
    stop("No such method available")
  }
   
  if(r > 5) warning("We recommend using r <=5")
  if(r < 0) stop("r cannot be negative")
  # making matrix compatible and applying g
  chain <- as.matrix(x)
  if(!is.matrix(chain) && !is.data.frame(chain))
    stop("'x' must be a matrix or data frame.")

  if (is.function(g)) 
  {
    chain <- apply(chain, 1, g)

    if(is.vector(chain))
    {
      chain <- as.matrix(chain)
    }else
    {
      chain <- t(chain)
    }
  }

  
  ## Setting dimensions on the mcmc output. 
  n = dim(chain)[1]
  p = dim(chain)[2]

  # Setting batch size based on chosen option
  if(is.null(size))
  {
    b <- batchSize(x = x, method = method, g = g)  # optimal
  }
  else if(size == "sqroot")
  {
    b <- floor(sqrt(n))
  } 
  else if(size == "cuberoot") {
    b <- floor(n^(1/3))
  }
  else {
    if (!is.numeric(size) || size < 1 || size >= n || floor(n/size) <=1) 
        stop("size is either too large, too small, or not a number")

    b <- floor(size)
  }
  a <- floor(n/b)
  if(b == 1 && r != 1)
  {
    r = 1
    message <- "r was set to 1 since b = 1. "
  }
 ########################## 


  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)  # this is based on the full output. not n = a*b
  
  ## Setting matrix sizes to avoid dynamic memory 
  sig.mat = matrix(0, nrow = p, ncol = p)

  if(floor(b/r) < 1) stop("Either decrease r or increase n")

  message <- ""   # will store some info for blather

  ## Batch Means
  if(method == "bm")
  {
    bm.mat <- mbmC(chain, b)
    sig.mat <- bm.mat
    method.used <- "Batch Means"
    if(r > 1)
    {
      sig.mat <- 2*bm.mat - mbmC(chain, floor(b/r))
      method.used <- paste("Lugsail Batch Means with r = ", r)
      if(prod(diag(sig.mat) > 0) == 0)  # If diagonals are negative, cannot use larger values of r
      {
        sig.mat <- bm.mat
        method.used <- "Batch Means"
        message <- paste(message, paste("Diagonals were negative with r = ", r,". r = 1 was used.", sep = ""), sep = "")
      }
    }
  }

  # Overlapping Batch means
  if(method == "obm")
  {
    obm.mat <- mobmC(chain, b)
    sig.mat <- obm.mat
    method.used <- "Overlapping Batch Means"
    if(r > 1)
    {
      sig.mat <- 2*obm.mat - mobmC(chain, floor(b/r))
      method.used <- paste("Lugsail Overlapping Batch Means with r = ", r)
      if(prod(diag(sig.mat) > 0) == 0)  # If diagonals are negative, cannot use larger values of r
      {
        sig.mat <- obm.mat
        method.used <- "Overlapping Batch Means"
        message <- paste(message, paste("Diagonals were negative with r = ", r,". r = 1 was used.", sep = ""), sep = "")
      }
    }
  }

  # SVE with Bartlett window
  if(method == "bartlett")
  {
   chain <- scale(chain, center = mu.hat, scale = FALSE)
   bar.mat <- msveC(chain, b, "bartlett")
   sig.mat <- bar.mat
   method.used <- "Bartlett Spectral Variance"
   if(r > 1)
   {
    sig.mat <- 2*bar.mat - msveC(chain, floor(b/r), "bartlett")
    method.used <- paste("Lugsail Bartlett Spectral Variance with r = ", r)
    if(prod(diag(sig.mat) > 0) == 0)  # If diagonals are negative, cannot use larger values of r
    {
      sig.mat <- bar.mat
      method.used <- "Bartlett Spectral Variance"
      message <- paste(message, paste("Diagonals were negative with r = ", r,". r = 1 was used.", sep = ""), sep = "")
    }
   }
  }

  ## SVE with tukey window
  if(method == "tukey")
  {
   chain <- scale(chain, center = mu.hat, scale = FALSE)
   tuk.mat <- msveC(chain, b, "tukey")
   sig.mat <- tuk.mat
   method.used <- "Tukey Spectral Variance"
   if(r > 1)
   {
    sig.mat <- 2*tuk.mat - msveC(chain, floor(b/r), "tukey")
    method.used <- paste("Lugsail Tukey Spectral Variance with r = ", r)
    if(prod(diag(sig.mat) > 0) == 0)  # If diagonals are negative, cannot use larger values of r
    {
      sig.mat <- tuk.mat
      method.used <- "Tukey Spectral Variance"
      message <- paste(message, paste("Diagonals were negative with r = ", r,". r = 1 was used.", sep = ""), sep = "")
    }
   }
  }

  adjust.used <- FALSE  #whether an adjustment was made. None yet
  
  if(adjust) # if adjust is FALSE, may output non PD estimator
  {
    sig.eigen <- eigen(sig.mat, only.values = TRUE)$values
    if(min(sig.eigen) <= 0)  #needs an adjustment. No need to adjust is not needed
    {
      adjust.used <- TRUE
      sig.mat <- adjust_matrix(sig.mat, N = n) 
    }    
  } 

  if(blather)
  {
    return(list("cov" = sig.mat, "est" = mu.hat, "nsim" = n, 
        "method" = method.used, "size" = b, "Adjustment-Used" = adjust.used, "message" = message))
  } else {
    return(list("cov" = sig.mat, "est" = mu.hat))
  }

}






