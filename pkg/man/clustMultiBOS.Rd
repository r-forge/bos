\name{clustMultiBOS}
\alias{clustMultiBOS}
\title{
Function to cluster multivariate ordinal data
}
\description{
This function performs the cluster analysis of multivariate ordinal data.
}
\usage{
clustMultiBOS (x,m,k,ntrials=1)
}

\arguments{
  \item{x}{
 ordinal data set [n,d]
}
  \item{m}{
nb of modalities [d]
}
  \item{k}{
number of components
}
  \item{ntrials}{
ntrials nb of random trials for EM OPTIONAL (default is 1)
}

}
\value{ list of the model parameters, the partition and the information criteria. 
}


\examples{
\dontrun{
  data(AERS)

  res <- clustMultiBOS(AERS,c(4,4,4,4),2)

  print(res$bic)
}
}

