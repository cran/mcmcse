multiESS <- function(x, covmat = NULL, g = NULL)
{

	chain <- as.matrix(x)
	if(!is.matrix(x) && !is.data.frame(x))
	  stop("'x' must be a matrix or data frame.")

	if (is.function(g)) 
	{
	  chain <- apply(x, 1, g)

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

	if(is.matrix(covmat))
	{
		var_mat <- cov(chain)
		ess <- n*(det(var_mat)/det(covmat))^(1/p)
	} else
	{
		covmat <- mcse.multi(chain)$cov
		var_mat <- cov(chain)
		ess <- n*(det(var_mat)/det(covmat))^(1/p)
	}
return(ess)

}