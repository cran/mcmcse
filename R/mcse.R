mcse <-
function(vals, bs="sqroot", g=NULL, meth="BM", warn=FALSE)
  {
  	g.check <- is.function(g)
  	if(g.check==FALSE) {
  		g <- function(x) return(x)  # default: identity function
    }
    N <- length(vals)
    if (N<1000)
      {
        if (warn) # if warning
          cat("WARNING: too few samples (less than 1000)\n")
        if (N<10)
          return(NA)
      }

    if (bs=="sqroot") 
      {
        b <- floor(sqrt(N)) # batch size
        a <- floor(N/b) # number of batches
      }
    else # batch size provided
      {
        stopifnot(is.numeric(bs))  
        b <- floor(bs) # batch size
        if (b > 1) # batch size valid
          a <- floor(N/b) # number of batches
        else
          stop("batch size invalid (bs=",bs,")")
      }
    
    if (meth=="BM")
    {
      Ys <- sapply(1:a,function(k) return(mean(g(vals[((k-1)*b+1):(k*b)]))) )
      muhat <- mean(Ys) 
      sigmahatsq <- b*sum((Ys-muhat)^2)/(a-1)

      se <- sqrt(sigmahatsq/N)
      return(list(est=muhat,se=se))
    }
    
    if (meth=="OBM")
    {
    	a <- N - b + 1
    	Ys <- sapply(1:a, function(k) return(mean(g(vals[k:(k+b-1)]))) )
    	muhat <- mean(g(vals))
    	sigmahatsq <- N*b*sum((Ys-muhat)^2)/(a-1)/a

    	se <- sqrt(sigmahatsq/N)
        return(list(est=muhat,se=se))
    }
    	
    if (meth=="TukeyHanning")
    {
	alpha <- seq(1, b, 1)
	alpha <- (1 + cos(pi*alpha/b))/2 * (1 - alpha/N)
    	muhat <- mean(g(vals))
    	R <- sapply(0:b, function(j) return(mean((g(vals[1:(N-j)]) - muhat) * ( g(vals[(j+1):N]) - muhat)) ) )
    	sigmahatsq <- (R[1] + 2* sum(alpha * R[-1]))
    	se <- sqrt(sigmahatsq / N)
        return(list(est=muhat,se=se))
    }
    
    if (meth=="Bartlett")
    {
	alpha <- seq(1, b, 1)
	alpha <- (1 - abs(alpha)/b) * (1 - alpha/N)
    	muhat <- mean(g(vals))
    	R <- sapply(0:b, function(j) return(mean((g(vals[1:(N-j)]) - muhat) * ( g(vals[(j+1):N]) - muhat)) ) )
    	sigmahatsq <- (R[1] + 2* sum(alpha * R[-1]))
    	se <- sqrt(sigmahatsq / N)
        return(list(est=muhat,se=se))
    }
    
    else # method not valid
    {
    	stop("method specified invalid (meth=" , meth, ")")
    }  	
  }

