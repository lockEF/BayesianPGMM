# BayesianPGMM

This is an R package to estimates a Bayesian piecewise growth mixture model with linear segments, for a given number of latent classes and a latent number of possible change points in each class. See [1] for methodological details.
This package requires Just Another Gibbs Sampler (JAGS) to be installed on your computer (http://mcmc-jags.sourceforge.net/), and depends on the packages `rjags`  and `label.switching`. 

The `BayesianPGMM` package can then be installed, directly from GitHub, using the devtools library:

```
install.packages("devtools")
library(devtools)
install_github("lockEF/BayesianPGMM")
``` 
The additional R function BayesBPLMEM.R estimates a bivariate piecewise linear mixed-effects model.  This function was developed by Yadira Peralta Torres and is described in [2]. SimDataBPLMEMExample.R loads the required libraries and includes example code to run the BayesPLMEM function on the simulated data in the SimDataBPLMEM.rda file.  

[1] Lock, E.F., Kohli, N., & Bose, M. (2018). Detecting multiple random changepoints in Bayesian piecewise growth mixture models. <em>Psychometrika</em>, 83 (3): 733-750. Preprint: https://arxiv.org/abs/1710.10626 

Peralta Torres, Y. E. (2018). Bayesian modeling of associations in bivariate mixed-effects models for segmented growth curves (Doctoral dissertation). Retrieved from the University of Minnesota Digital Conservancy, http://hdl.handle.net/11299/201670 .
