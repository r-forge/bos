\name{BOS-package}
\alias{BOS-package}
\alias{BOS}
\docType{package}
\title{Model-Based Clustering of Multivariate Ordinal Data}
\description{
We design the first univariate probability distribution for ordinal data which strictly respects the ordinal nature of data. More precisely, it relies only on order comparisons between modalities. Contrariwise, most competitors either forget the order information or add a nonexistent distance information. The proposed distribution is obtained by modeling the data generating process which is assumed, from optimality arguments, to be a stochastic binary search algorithm in a sorted table. The resulting distribution is natively governed by two meaningful parameters (position and precision) and has very appealing properties: decrease around the mode, shape tuning from uniformity to a Dirac, identifiability. Moreover, it is easily estimated by an EM algorithm since the path in the stochastic binary search algorithm is missing. Using then the classical latent class assumption, the previous univariate ordinal model is straightforwardly extended to model-based clustering for multivariate ordinal data.}
\details{
\tabular{ll}{
Package: \tab BOS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-02-07\cr
License: \tab GPL (>= 2)\cr
}
The main function is clustMultiBOS.
}
\author{
Christophe Biernacki and Julien Jacques
Maintainer: Julien Jacques <julien.jacques@univ-lyon2.fr>
}
\references{
 C. Biernacki and J. Jacques (2014), Model-based clustering of multivariate ordinal data relying on a stochastic binary search algorithm, to appear in Statistics and Computing.
}
\examples{
\dontrun{
  data(AERS)

  res <- clustMultiBOS(AERS,c(4,4,4,4),2)

  print(res$bic)
}
}

