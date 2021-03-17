\name{getStandardizedCovMeanDiffs}
\alias{getStandardizedCovMeanDiffs}
\title{Standardized Covariate Mean Differences}
\description{
\code{getStandardizedCovMeanDiffs} computes the standardized covariate mean differences between a treatment and control group, defined as treatment minus control. The standardized covariate mean differences are defined as the covariate mean differences divided by the square-root of the pooled variance between groups.
}
\usage{
getStandardizedCovMeanDiffs(X.matched, indicator.matched,
	X.full = NULL, indicator.full = NULL)
}
\arguments{
  \item{X.matched}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates) for the matched dataset.}
  \item{indicator.matched}{A vector of 1s and 0s (e.g., denoting treatment and control) for the matched dataset.}
  \item{X.full}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates) for the full, unmatched dataset if available.}
  \item{indicator.full}{A vector of 1s and 0s (e.g., denoting treatment and control) for the full, unmatched dataset if available.}
}
\value{
The standardized covariate mean differences between a treatment and control group, defined as treatment minus control.
}

\details{
	The arguments X.full and indicator.full (i.e., the covariate matrix and indicator for the full, unmatched dataset) are only used to correctly define the standardized covariate mean differences. Technically, the covariate mean differences should be standardized by the pooled variance within the full, unmatched dataset, instead of within the matched dataset. If X.full and indicator.full are unspecified, the pooled variance within the matched dataset is used for standardization instead. This distinction rarely leads to large differences in the resulting standardized covariate mean differences, and so researchers should feel comfortable only specifying X.matched and indicator.matched if only a matched dataset is available. Furthermore, if one wants to compute the standardized mean differences for a full, unmatched dataset, then they should only specify X.matched and indicator.matched.
}

\author{Zach Branson
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
	
	#the standardized covariate mean differences
	#for these three datasets are:
	getStandardizedCovMeanDiffs(
		X.matched = X.lalonde,
		indicator.matched = indicator.lalonde)
	getStandardizedCovMeanDiffs(
		X.matched = X.matched.ps,
		indicator.matched = indicator.matched.ps,
		X.full = X.lalonde,
		indicator.full = indicator.lalonde)
	getStandardizedCovMeanDiffs(
		X.matched = X.matched.card,
		indicator.matched = indicator.matched.card,
		X.full = X.lalonde,
		indicator.full = indicator.lalonde)
}

\seealso{
	See also \code{lalondeMatches} for details about the Lalonde and matched datasets.
}

\keyword{Covariate mean difference.}
