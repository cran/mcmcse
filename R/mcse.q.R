mcse.q <-
function(vals, qval, bs="sqroot", g=NULL, meth="BM", warn=FALSE)
  {
  	g.check <- is.function(g)
  	if(g.check==FALSE) {
  		g <- function(x) return(x)  # default: identity function
    }
    counting <- function(var.vector, var.number) {
    	return(length(var.vector[var.vector<=var.number]))
    }
    N <- length(vals)
    if (N<1000)
      {
        if (warn) # if warning
          cat("WARNING: too few samples (less than 1000)\n")
        if (N<10)
          return(NA)
      }

    if (qval <= 0 || qval >= 1)
      {
        stop("quantile invalid (qval=",qval,")")
      }
    quant <- function(input){ quantile(input, prob=qval, type=1, names=FALSE) }

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
    	hat.xi <- quant(g(vals))
    	Ys <- sapply(1:a, function(k) return(counting(g(vals[((k-1)*b+1):(k*b)]), hat.xi))) / b
    	muhat <- mean(Ys)
    	sigmahatsq <- b*sum((Ys-muhat)^2)/(a-1)
    	f.hat.junk <- density(g(vals), from=hat.xi, to=hat.xi, n=1)
    	f.hat <- f.hat.junk$y
    	se <- sqrt(sigmahatsq/N) / f.hat
    	return(list(est=hat.xi,se=se))
    }
    
    if (meth=="OBM")
    {
    	hat.xi <- quant(g(vals))
    	a <- N - b + 1
    	Ys <- sapply(1:a, function(k) return(counting(g(vals[k:(k+b-1)]), hat.xi))) / b
    	muhat <- mean(Ys)
    	sigmahatsq <- N*b*sum((Ys-muhat)^2)/(a-1)/a
    	f.hat.junk <- density(g(vals), from=hat.xi, to=hat.xi, n=1)
    	f.hat <- f.hat.junk$y
    	se <- sqrt(sigmahatsq/N) / f.hat
    	return(list(est=hat.xi,se=se))
    }
    	
    if (meth=="Subsampling")
    {
    	hat.xi <- quant(g(vals))
    	a <- N - b + 1
    	Ys <- sapply(1:a, function(k) return(quant(g(vals[k:(k+b-1)]))) )
    	muhat <- mean(Ys)
    	sigmahatsq <- N*b*sum((Ys-muhat)^2)/(a-1)/a
    	se <- sqrt(sigmahatsq / N)
        return(list(est=hat.xi,se=se))
    }
    
    else # method not valid
    {
    	stop("method specified invalid (meth=" , meth, ")")
    }  	
  }

