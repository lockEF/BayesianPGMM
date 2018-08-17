\name{plotPGMM}
\alias{plotPGMM}
\title{
Plot a Bayesian Piecewise Growth Mixture Model
}
\description{
Generates a "spaghetti plot" for longitudinal data for multiple subjects.  If the results from the BayesBGMM function are supplied, each trajectory is colored by the estimated class membership of its subject, and the trajectory defined by the mean parameters for each class is shown in bold.\\ 
}
\usage{
plotPGMM <- (X,Y,Results=NA,xlab='X',ylab='Y',Colors=NULL,MeanColors=NULL,...)
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
 \item{Colors}{
  	Colors for each class.  Default is blue for class 1, red for class 2, green for class 3, gold for class 4, gray for class 5.   
}
 \item{MeanColors}{
  	Colors for the trajectory defined by the mean parameters for each class.   
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
plotPGMM(X,Y) ##Plot the data
Results <- BayesPGMM(X,Y) ##Fit PGMM (can take about 5 minutes)
plotPGMM(X,Y,Results) ##Plot results
}


