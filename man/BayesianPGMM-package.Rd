\name{BayesianPGMM-package}
\alias{BayesianPGMM}
\docType{package}
\title{
Bayesian Piecewise Growth Mixture Modeling
}
\description{
Estimates a Bayesian piecewise growth mixture model with linear segments, for a given number of latent classes and a latent number of possible change points  in each class. See [1] for the methodological details, and see ?BayesPGMM for details on the input and output. 
}
\details{
\tabular{ll}{
Package: \tab BayesianPGMM-package\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2020-8-12\cr
License: \tab GPL-3\cr
}

}
\author{
Eric F. Lock

Maintainer: Eric F. Lock
}

\references{
[1] EF Lock, N Kohli & M Bose (2018). Detecting multiple random changepoints in Bayesian piecewise growth mixture models. Psychometrika, 83 (3): 733-750.  
}

\examples{
data(SimData)
plotPGMM(X,Y) ##Plot the data
Results <- BayesPGMM(X,Y) ##Fit PGMM (can take about 5 minutes)
plotPGMM(X,Y,Results) ##Plot results
}


\keyword{ package }
