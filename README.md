# mcmcse
An R package for computing Monte Carlo standard
errors (MCSE) in Markov chain Monte Carlo (MCMC) settings. MCSE
computation for expectation and quantile estimators is
supported as well as multivariate estimations. The package also provides 
functions for computing effective sample size and for plotting
Monte Carlo estimates versus sample size.


# Installation
This R package is on CRAN, and its preferred URL being https://CRAN.R-project.org/package=mcmcse.

To download this development repo,  through the the `devtools` package:

```{r}
# install.packages("devtools")
library(devtools)
devtools::install_github("dvats/mcmcse")
```
# Citation
Please run `citation("mcmcse")` after loading the package for citation details.
