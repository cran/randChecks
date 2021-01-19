\name{getMD}
\alias{getMD}
\title{Mahalanobis Distance}
\description{
\code{getMD} computes the Mahalanobis distance of the covariate means between a treatment and control group.
}
\usage{
getMD(X, indicator, covX.inv = NULL)
}
\arguments{
  \item{X}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates).}
  \item{indicator}{A vector of 1s and 0s (e.g., denoting treatment and control).}
  \item{covX.inv}{The inverse of X's covariance matrix. Almost always this should be set to NULL, and \code{getMD} will compute the inverse of the covariance matrix automatically.}
}
\value{
The Mahalanobis distance of the covariate means between a treatment and control group.
}

\examples{
	#This loads the classic Lalonde (1986) dataset,
	#as well as two matched datasets:
	#one from 1:1 propensity score matching,
	#and one from cardinality matching, where
	#the standardized covariate mean differences are all below 0.1.
	data("lalondeMatches")

	#obtain the covariates for these datasets
	X.lalonde = subset(lalonde, select = -c(treat))
	X.matched.ps = subset(lalonde.matched.ps, select = -c(treat,subclass))
	X.matched.card = subset(lalonde.matched.card, select = -c(treat,subclass))
	#the treatment indicators are
	indicator.lalonde = lalonde$treat
	indicator.matched.ps = lalonde.matched.ps$treat
	indicator.matched.card = lalonde.matched.card$treat
	
	#the Mahalanobis distance for each dataset is:
	getMD(X = X.lalonde, indicator = indicator.lalonde)
	getMD(X = X.matched.ps, indicator = indicator.matched.ps)
	getMD(X = X.matched.card, indicator = indicator.matched.card)
}

\seealso{
	See also \code{lalondeMatches} for details about the Lalonde and matched datasets.
}

\references{
Mahalanobis, P. C. (1936). On the generalized distance in statistics. National Institute of Science of India, 1936.
}
\author{Zach Branson
}


\keyword{Mahalanobis distance}
