\name{asIfRandPlot}
\alias{asIfRandPlot}
\title{Graphical Diagnostic of As-If Randomization for Different Assignment Mechanisms}
\description{
\code{asIfRandPlot} produces a plot showing the distribution of the Mahalanobis distance for different assignment mechanisms, along with the observed Mahalanobis distance. If the observed Mahalanobis distance is well within the range of a particular distribution, then that suggests that a particular assignment mechanism holds. This function supports the following assignment mechanisms:
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
asIfRandPlot(X, indicator,
  assignment = c("complete"),
  subclass = NULL, threshold = NULL,
  perms = 1000)
}
\arguments{
  \item{X}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates).}
  \item{indicator}{A vector of 1s and 0s (e.g., denoting treatment and control).}
  \item{assignment}{A vector of assignment mechanisms that the user wants to visualize; the user can test one assignment mechanism or multiple. The possible choices are "complete", "blocked", "constrained diffs", "constrained md", "blocked constrained diffs", and "blocked constrained md". See Description for more details on these assignment mechanisms.}
  \item{subclass}{A vector denoting the subclass/block for each subject/unit. This must be specified only if one of the blocked assignment mechanisms are used.}
  \item{threshold}{The threshold used within the constrained assignment mechanisms; thus, this must be specified only if one of the constrained assignment mechanisms are used. This can be a single number or a vector of numbers (e.g., if one wants to use a different threshold for each covariate when testing constrained-differences randomization).}
  \item{perms}{The number of permutations used within the randomization test. A larger number requires more computation time but results in a more consistent p-value.}
}
\value{
A plot showing the distribution of the Mahalanobis distance for different assignment mechanisms, along with the observed Mahalanobis distance. Also returns a p-value for each assignment mechanism - this is simply the area of the distribution more extreme than the observed Mahalanobis distance. This is the same as \code{asIfRandTest()} using the Mahalanobis distance as a test statistic.
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
	

	#The following lines of code create diagnostic plots assessing
	#whether the treatment follows different assignment mechanisms.
	
	#Note that the following examples only use 100 permutations
	#to approximate the randomization distribution.
	#In practice, we recommend setting perms = 1000 or more;
	#in these examples we use perms = 50 to save computation time.
	
	#Assessing complete randomization for the full dataset
	#Here, complete randomization clearly does not hold,
	#because the observed Mahalanobis distance is far outside
	#the complete randomization distribution.
	asIfRandPlot(X = X.lalonde, indicator = indicator.lalonde, perms = 50)

	#Assessing complete and block (paired) randomization for
	#the propensity score matched dataset
	#Again, complete and block randomization appear to not hold
	#because the observed Mahalanobis distance is far outside
	#the randomization distributions.
	asIfRandTest(X = X.matched.ps, indicator = indicator.matched.ps,
  		assignment = c("complete", "blocked"),
  		subclass = lalonde.matched.ps$subclass,
  		perms = 50)
	
	#Assessing three assignment mechanisms for the
	#cardinality matched dataset:
	# 1) complete randomization
	# 2) blocked (paired) randomization
	# 3) constrained-MD randomization
	#Note that the Mahalanobis distance is approximately a chi^2_K distribution,
	#where K is the number of covariates. In the Lalonde data, K = 8.
	#Thus, the threshold can be chosen as the quantile of the chi^2_8 distribution.
	#This threshold constrains the Mahalanobis distance to be below the 25-percent quantile:
	a = qchisq(p = 0.25, df = 8)
	#Then, we can assess these three assignment mechanisms with the plot below.
	#Here, these assignment mechanisms seem plausible,
	#because the observed Mahalanobis distance is well
	#within the randomization distributions.
	asIfRandPlot(X = X.matched.card, indicator = indicator.matched.card,
  		assignment = c("complete", "blocked", "constrained md"),
  		subclass = lalonde.matched.card$subclass,
  		threshold = a,
  		perms = 50)
}

\author{Zach Branson
}


\keyword{randomization test}
\keyword{complete randomization}
\keyword{block randomization}
\keyword{constrained randomization}
\keyword{Mahalanobis distance}
