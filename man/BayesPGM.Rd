\name{BayesPGM}
\alias{BayesPGM}
\title{
Fit Bayesian Piecewise Growth Model
}
\description{
Estimates a Bayesian piecewise growth model with linear segments, with a given latent number of possible change points. See [1] for methodological details.  To fit a mixture model with two or more classes, see BayesPGMM. 
}
\usage{
BayesPGM(X,Y,max_cp=2,cp_prior='binomial',binom_prob=0.5,
		scale_prior='uniform', iters_BurnIn=20000,
		iters_sampling=30000,Thin=15,SaveChains=FALSE)
}
\arguments{
  \item{X}{An array of dimensions N x M,  where the ij'th value gives the j'th time point for subject i.
}
  \item{Y}{An array of dimension N x M, where the ij'th value gives the outcome for subject i at the j'th time point.  Missing measurements should have the value NA.  
}
  \item{max_cp}{Maximum number of possible changepoints. 
}
 \item{cp_prior}{Prior for the number of changepoints, with options 'binomial' or 'uniform'.  Default is 'binomial'. 
}
\item{binom_prob}{Probability for binomial prior, if specified.  Default is 0.5.   
}
\item{scale_prior}{Prior for the scale parameter for the hierarchical random effects.  The default is 'uniform', a scaled uniform prior; the option 'hc' is a scaled half-cauchy prior. See [1] for more details.}
  \item{iters_BurnIn}{(optional) Number of Gibbs sampling iterations to run for the burn in.  Default is 20000.
}
  \item{iters_sampling}{(optional) Number of Gibbs sampling iterations to run for posterior sampling.  Default is 30000.
}
  \item{Thin}{(optional) Will save every k'th posterior sample, where k=Thin. Default is 15.  
}
 \item{SaveChains}{(optional) If TRUE, raw MCMC samples from rjags will be returned. 
}
}
\value{
%%  ~Describe the value returned 
\item{Gelman.msrf}{Gelman multivariate scale reduction factor to assess convergence, with the potential scale reduction factor (psrf) for each parameter}

\item{y.mean}{NxM array giving the fitted value at each time point for each subject.
}
\item{error.sd}{The residual standard deviation.
}
\item{error.sd.CI}{95\% CI for the residual standard deviation.
}
\item{DIC}{Deviance information criterion for fitted model.
}
\item{Save}{If SaveChains=TRUE in input, raw MCMC samples from rjags
}
\item{K_prob}{Vector of posterior probabilities for the number of changepoints [0,1,...,K]
}
\item{K}{List where K[[k+1]] gives results for the model with k changepoints.  
}
\item{K[[k+1]]$b.mean}{Vector giving the coefficient means [intercept,slope,change at cp1,change at cp2,...] for the k-changepoint model
}
\item{K[[k+1]]$b.sd}{Vector giving the coefficient standard deviations for the k-changepoint model
}
\item{K[[k+1]]$b.mean.CI}{95\% CI for the coefficient means for the k-changepoint model
}
\item{K[[k+1]]$b.sd.CI}{95\% CI for the coefficient standard deviations for the k-changepoint model
}
\item{K[[k+1]]$cp.mean}{Vector giving mean location of each changepoint in the k-changepoint model
}
\item{K[[k+1]]$cp.sd}{Vector giving standard deviation of each changepoint in the k-changepoint model
}
}

\details{
	For more information on the priors implemented in this package, see [1].  For more control over the priors and aspects of posterior computation, the source functions can be modified directly (enter BayesPGMM in the R console to view source code).
}

\author{
XXXX
}

\references{
[1] XXXX authors XXXX (2017). Detecting multiple random changepoints in Bayesian piecewise growth mixture models. Preprint.
}

\examples{
data(SimData)  ##load simple simulated dataset
X=X[1:10,];Y=Y[1:10,]  ###select subjects from first class only
plotPGM(X,Y) ##Plot the data
Results <- BayesPGM(X,Y) ##Fit PGM (can take about 5 minutes)
plotPGM(X,Y,Results) ##Plot results
}


