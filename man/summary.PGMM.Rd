\name{summary.PGMM}
\alias{summary.PGMM}
\title{
Summarize the results of piecewise growth mixture model (PGMM)
}
\description{
Provides a summary of a PGMM model, as returned by BayesPGMM(...). 
}
\usage{

\method{summary}{PGMM}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "PGMM" (returned by BayesPGMM(...)).
}
  \item{\dots}{
Additional arguments. 
}
}

\value{
Prints a table with estimates for key parameters in each class identified by BayesPGMM(...).  
}


\author{
Eric F. Lock
}


\examples{
data(SimData)  ##load simple simulated dataset
plotPGMM(X,Y) ##Plot the data
Results <- BayesPGMM(X,Y) ##Fit PGMM (can take about 5 minutes)
summary(Results) ##summarize results
}


