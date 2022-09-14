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

The additional R function BayesPCREM.R estimates a piecewise linear mixed-effects model with crossed random effects, developed by Corissa Rohloff [3]. SimDataBPCREM_Example.R loads the required libraries and includes example code to run the function on the data in the SimDataBPCREM.Rdata file.  

[1] Lock, E.F., Kohli, N., & Bose, M. (2018). Detecting multiple random changepoints in Bayesian piecewise growth mixture models. <em>Psychometrika</em>, 83 (3): 733-750. Preprint: https://arxiv.org/abs/1710.10626 

[2] Peralta Y., Kohli N., Lock E.F., & Davison M.L. Bayesian modeling of associations in bivariate piecewise linear mixed-effects models. <em>Psychological Methods</em>, 2020. https://psycnet.apa.org/doiLanding?doi=10.1037%2Fmet0000358

[3] Rohloff C.T., Kohli N., & Lock E.F. Identifiability and Estimability of Bayesian Linear and Nonlinear Crossed Random Effects Models. In Preparation. 

