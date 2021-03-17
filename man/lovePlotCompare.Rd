\name{lovePlotCompare}
\alias{lovePlotCompare}
\title{Love Plot of Standardized Covariate Mean Differences (along with Permutation Quantiles)}
\description{
\code{lovePlotCompare} produces a Love plot displaying the standardized covariate mean differences (produced by \code{getStandardizedCovMeanDiffs()}) for two different datasets. The dataset with smaller covariate mean differences is deemed the "more balanced" dataset; this is particularly useful when comparing a full dataset to a matched dataset.
}
\usage{
lovePlotCompare(X1, indicator1, X2, indicator2, dataNames = c("Dataset1", "Dataset2"),
	X.full = NULL, indicator.full = NULL)
}
\arguments{
  \item{X1}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates) for one dataset.}
  \item{indicator1}{A vector of 1s and 0s (e.g., denoting treatment and control) for one dataset.}
  \item{X2}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates) for another dataset.}
  \item{indicator2}{A vector of 1s and 0s (e.g., denoting treatment and control) for another dataset.}
  \item{dataNames}{A two-length vector denoting the names of the datasets (used in the legend of the plot).}
  \item{X.full}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates) for the full, unmatched dataset if available.}
  \item{indicator.full}{A vector of 1s and 0s (e.g., denoting treatment and control) for the full, unmatched dataset if available.}
}
\value{
A Love plot displaying the standardized covariate mean differences for two datasets.
}

\details{
	Note that the covariate matrices X1 and X2 have to have the same number of columns and should correspond to the same covariates. However, they do not have to have the same number of rows (i.e., the same number of subjects/units).

	Furthermore, the arguments X.full and indicator.full (i.e., the covariate matrix and indicator for the full, unmatched dataset) are only used to correctly define the standardized covariate mean differences. Technically, the covariate mean differences should be standardized by the pooled variance within the full, unmatched dataset, instead of within the matched dataset. If X.full and indicator.full are unspecified, the pooled variance within the matched dataset is used for standardization instead. This distinction rarely leads to large differences in the resulting standardized covariate mean differences, and so researchers should feel comfortable only specifying X1, X2, indicator1, and indicator2 if a full, unmatched dataset is not available.
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

	#The following code will display a classic Love plot
	#(with a dot for each standardized covariate mean difference),
	#where there are differently-colored dots for each dataset.
	
	#full lalonde dataset vs ps matched dataset
	lovePlotCompare(X1 = X.lalonde, indicator1 = indicator.lalonde,
	  X2 = X.matched.ps, indicator2 = indicator.matched.ps,
	  X.full = X.lalonde, indicator.full = indicator.lalonde,
	  dataNames = c("unmatched", "ps matched"))
	  
	#ps vs card
	lovePlotCompare(X1 = X.matched.ps, indicator1 = indicator.matched.ps,
	  X2 = X.matched.card, indicator2 = indicator.matched.card,
	  X.full = X.lalonde, indicator.full = indicator.lalonde,
	  dataNames = c("ps matched", "card matched"))
}

\author{Zach Branson
}

\keyword{Love plot}
\keyword{covariate mean difference}
