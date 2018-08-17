\name{plotPGM}
\alias{plotPGM}
\title{
Plot a Bayesian Piecewise Growth Model
}
\description{
Generates a "spaghetti plot" for longitudinal data for multiple subjects.  If the results from the BayesBGMMfunction are supplied, each trajectory is colored gray, and the trajectory defined by the mean parameters is shown in bold.\\ 
}
\usage{
plotPGMM <- (X,Y,Results=NA,xlab='X',ylab='Y',...)
}
\arguments{
  \item{X}{
  	An array of dimensions N x M,  where the ij'th value gives the j'th time point for subject i.
}
  \item{Y}{
  	An array of dimension N x M, where the ij'th value gives the outcome for subject i at the j'th time point.
}
  \item{Results}{
  	The output of an application of BayesPGMM() to the data X,Y.  If Results=NA, only a spaghetti plot of the data will be generated.  
}
  \item{xlab}{
  	X-axis label for the generated plot. 
}
 \item{ylab}{
  	Y-axis label for the generated plot. 
}
  \item{...}{
  	(optional) Other parameters to pass to the plot() function.
}
}


\author{
Eric F. Lock
}


\examples{
	data(SimData)
X=X[1:10,]; Y=Y[1:10,] ###Take 10 samples belonging to single class
plotPGM(X,Y) ##Plot the data
Results <- BayesPGM(X,Y) ##Fit PGM (can take about 5 minutes)
plotPGM(X,Y,Results) ##Plot results
}


