\name{getCovMeanDiffs}
\alias{getCovMeanDiffs}
\title{Covariate Mean Differences}
\description{
\code{getCovMeanDiffs} computes the covariate mean differences between a treatment and control group.
}
\usage{
getCovMeanDiffs(X, indicator)
}
\arguments{
  \item{X}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates).}
  \item{indicator}{A vector of 1s and 0s (e.g., denoting treatment and control).}
}
\value{
The covariate mean differences between a treatment and control group, defined as treatment minus control.
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
	
	#the covariate mean differences are:
	getCovMeanDiffs(X = X.lalonde, indicator = indicator.lalonde)
	getCovMeanDiffs(X = X.matched.ps, indicator = indicator.matched.ps)
	getCovMeanDiffs(X = X.matched.card, indicator = indicator.matched.card)
}

\seealso{
	See also \code{lalondeMatches} for details about the Lalonde and matched datasets.
}

\author{Zach Branson
}


\keyword{Covariate mean difference.}
