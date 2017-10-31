# BayesianPGMM

This is an R package to estimates a Bayesian piecewise growth mixture model with linear segments, for a given number of latent classes and a latent number of possible change points in each class. See [1] for methodological details.
This package requires Just Another Gibbs Sampler (JAGS) to be installed on your computer (http://mcmc-jags.sourceforge.net/), and depends on the packages `rjags`  and `label.switching`. 

The `BayesianPGMM` package can then be installed, directly from GitHub, using the devtools library:

```
install.packages(devtools)
library(devtools)
install_github("lockEF/BayesianPGMM")
``` 


[1] Lock, E.F., Kohli, N., & Bose, M. (2017). Detecting multiple random changepoints in Bayesian piecewise growth mixture models. Psychometrika, to appear. https://arxiv.org/abs/1710.10626 .