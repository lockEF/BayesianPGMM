\name{summary.PGM}
\alias{summary.PGM}
\title{
  Summarize the results of piecewise growth mixture model (PGM)
}
\description{
  Provides a summary of a PGM model, as returned by BayesPGM(...). 
}
\usage{
  
  \method{summary}{PGM}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object of class "PGM" (returned by BayesPGM(...)).
  }
  \item{\dots}{
    Additional arguments. 
  }
}

\value{
  Prints a table with estimates for key parameters in each class identified by BayesPGM(...).  
}


\author{
  Eric F. Lock
}


\examples{
  data(SimData)  ##load simple simulated dataset
  X=X[1:10,];Y=Y[1:10,]  ###select subjects from first class only
  Results <- BayesPGM(X,Y) ##Fit PGM
  summary(Results) ##summarize results
}