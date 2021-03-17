\name{getMD}
\alias{getMD}
\title{Mahalanobis Distance}
\description{
\code{getMD} computes the Mahalanobis distance of the covariate means between a treatment and control group.
}
\usage{
getMD(X.matched, indicator.matched,
	covX.inv = NULL,
	X.full = NULL)
}
\arguments{
  \item{X.matched}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates) for the matched dataset.}
  \item{indicator.matched}{A vector of 1s and 0s (e.g., denoting treatment and control) for the matched dataset.}
  \item{covX.inv}{The inverse of X's covariance matrix. Almost always this should be set to NULL, and \code{getMD} will compute the inverse of the covariance matrix automatically.}
  \item{X.full}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates) for the full, unmatched dataset if available.}
}
\value{
The Mahalanobis distance of the covariate means between a treatment and control group.
}

\details{
	The argument X.full (i.e., the covariate matrix for the full, unmatched dataset) is only used to correctly define the Mahalanobis distance after matching. Technically, the Mahalanobis distance should be standardized by the covariance matrix within the full, unmatched dataset, instead of within the matched dataset. If X.full is unspecified, the covariance matrix within the matched dataset is used instead. This distinction rarely leads to large differences in the resulting distance, and so researchers should feel comfortable only specifying X.matched and indicator.matched if only a matched dataset is available. Furthermore, if one wants to compute the Mahalanobis distance for a full, unmatched dataset, then they should only specify X.matched and indicator.matched.
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
	getMD(X.matched = X.lalonde, indicator.matched = indicator.lalonde)
	getMD(X.matched = X.matched.ps, indicator.matched = indicator.matched.ps,
		X.full = X.lalonde)
	getMD(X.matched = X.matched.card, indicator.matched = indicator.matched.card,
		X.full = X.lalonde)
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
