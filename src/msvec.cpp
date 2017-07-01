#include <RcppArmadillo.h>
#include <math.h> 

using namespace arma;
using namespace Rcpp;

// returns the lag window value for the corresponding window
double lag(int s, double b, String method)
{
  if(method == "bartlett")
  {
    return(1 - s/b);
  }
  else
  {
    return((1 + cos(PI * s/b))/2 );
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat msvec(const arma::mat& chain, double b, String method = "bartlett")
{
  int n = chain.n_rows;
  int p = chain.n_cols;

  // t_chain improves performance a little
  mat t_chain = trans(chain);
  mat out(p,p);
  mat dummy(p,p);
  dummy.zeros();
  out.zeros();

  for(int s = 1; s < b; ++s)
  {
    dummy = trans(chain.rows(0, n-s-1))*chain.rows(s, n-1);
    out +=  lag(s, b, method)*(dummy + trans(dummy));
  }

  out += trans(chain)*chain;
  return(out/n);
}
