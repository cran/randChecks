\name{lovePlot}
\alias{lovePlot}
\title{Love Plot of Standardized Covariate Mean Differences (along with Permutation Quantiles)}
\description{
\code{lovePlot} produces a Love plot displaying the standardized covariate mean differences (produced by \code{getStandardizedCovMeanDiffs()}). This function can also produce permutation quantiles for different assignment mechanisms - if a standardized covariate mean difference is outside these quantiles, then that is evidence that the assignment mechanism does not hold. This function supports the following assignment mechanisms:
\itemize{
	\item Complete randomization ("complete"): Corresponds to random permutations of the indicator across units.
	\item Block randomization ("blocked"): Corresponds to random permutations of the indicator within blocks of units.
	\item Constrained-differences randomization ("constrained diffs"): Corresponds to random permutations of the indicator across units, conditional on the standardized covariate mean differences being below some threshold.
	\item Constrained-Mahalanobis randomization ("constrained md"): Corresponds to random permutations of the indicator across units, conditional on the Mahalanobis being below some threshold.
	\item Blocked Constrained-differences randomization ("blocked constrained diffs"): Corresponds to random permutations of the indicator within blocks of units, conditional on the standardized covariate mean differences being below some threshold.
	\item Blocked Constrained-Mahalanobis randomization ("blocked constrained md"): Corresponds to random permutations of the indicator within blocks of units, conditional on the Mahalanobis being below some threshold.
}

}
\usage{
lovePlot(X, indicator,
  permQuantiles = FALSE,
  assignment = "complete",
  subclass = NULL, threshold = NULL,
  alpha = 0.15, perms = 1000)
}
\arguments{
  \item{X}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates).}
  \item{indicator}{A vector of 1s and 0s (e.g., denoting treatment and control).}
  \item{permQuantiles}{Display permutation quantiles? TRUE or FALSE.}
  \item{assignment}{An assignment mechanism that the user wants to visualize. The possible choices are "complete", "blocked", "constrained diffs", "constrained md", "blocked constrained diffs", and "blocked constrained md". See Description for more details on these assignment mechanisms.}
  \item{subclass}{A vector denoting the subclass/block for each subject/unit. This must be specified only if one of the blocked assignment mechanisms are used.}
  \item{threshold}{The threshold used within the constrained assignment mechanisms; thus, this must be specified only if one of the constrained assignment mechanisms are used. This can be a single number or a vector of numbers (e.g., if one wants to use a different threshold for each covariate when testing constrained-differences randomization).}
  \item{alpha}{The alpha-level of the permutation quantiles, where the lower quantile is the alpha/2 quantile and the upper quantile is the 1-alpha/2 quantile. For example, if alpha = 0.15 (the default), then \code{lovePlot()} will display the 7.5-percent and 92.5-percent quantiles of the standardized covariate mean differences.}
  \item{perms}{The number of permutations used to compute the permutation quantiles. A larger number requires more computation time but results in a more consistent p-value.}
}
\value{
A Love plot displaying the standardized covariate mean differences. Can also produce permutation quantiles for different assignment mechanisms.
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

	#the subclass for the matched datasets are
	subclass.matched.ps = lalonde.matched.ps$subclass
	subclass.matched.card = lalonde.matched.card$subclass

	#The following code will display a classic Love plot
	#(with a dot for each standardized covariate mean difference):
	lovePlot(X = X.lalonde, indicator = indicator.lalonde)
	lovePlot(X = X.matched.ps, indicator = indicator.matched.ps)
	lovePlot(X = X.matched.card, indicator = indicator.matched.card)
	
	#The following lines of code create Love plots assessing
	#whether indicator.data follows different assignment mechanisms by
	#plotting the permutation quantiles
	
	#Note that the following examples only use 100 permutations
	#to approximate the randomization distribution.
	#In practice, we recommend setting perms = 1000 or more;
	#in these examples we use perms = 50 to save computation time.
	
	#Assessing complete randomization for the full dataset
	#Here we conclude complete randomization doesn't hold
	#because the standardized covariate mean differences
	#are almost all outside the quantiles.
	lovePlot(X = X.lalonde, indicator = indicator.lalonde,
		permQuantiles = TRUE,
		perms = 50)

	#assessing block (paired) randomization for
	#the 1:1 propensity score matched dataset
	#Many of the standardized covariate mean differences
	#are within the permutation quantiles,
	#but the race covariates (hispan and black)
	#are outside these quantiles.
	lovePlot(X = X.matched.ps, indicator = indicator.matched.ps,
  		permQuantiles = TRUE,
  		perms = 50,
  		assignment = "blocked", subclass = subclass.matched.ps)

	#assessing block (paired) randomization for
	#the cardinality matched dataset
	#All of the standardized covariate mean differences
	#are within the permutation quantiles
	lovePlot(X = X.matched.card, indicator = indicator.matched.card,
  		permQuantiles = TRUE,
  		perms = 50,
  		assignment = "blocked", subclass = subclass.matched.card)

	#assessing constrained randomization,
	#where the Mahalanobis distance is constrained.
	#Note that the Mahalanobis distance is approximately
	#a chi^2_K distribution, where K is the number of covariates.
	#In the Lalonde data, K = 8.
	#Thus, the threshold can be chosen as the quantile of the chi^2_8 distribution.
	#This threshold constrains the Mahalanobis distance to be below the 25-percent quantile:
	a = qchisq(p = 0.25, df = 8)
	#Then, the corresponding Love plot and permutation quantiles are:
	lovePlot(X = X.matched.card, indicator = indicator.matched.card,
  		permQuantiles = TRUE,
  		perms = 50,
  		assignment = "constrained md",
  		threshold = a)
}

\author{Zach Branson
}

\keyword{Love plot}
\keyword{randomization test}
\keyword{complete randomization}
\keyword{block randomization}
\keyword{constrained randomization}
\keyword{Mahalanobis distance}
\keyword{covariate mean difference}
