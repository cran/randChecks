\name{asIfRandTest}
\alias{asIfRandTest}
\title{As-If Randomization Test for Different Assignment Mechanisms: Global Tests and Covariate-by-Covariate Tests}
\description{
\code{asIfRandTest} computes p-values testing whether an indicator follows a given assignment mechanism, based on observed covariates. This function supports the following assignment mechanisms:
\itemize{
	\item Complete randomization ("complete"): Corresponds to random permutations of the indicator across units.
	\item Block randomization ("blocked"): Corresponds to random permutations of the indicator within blocks of units.
	\item Constrained-differences randomization ("constrained diffs"): Corresponds to random permutations of the indicator across units, conditional on the standardized covariate mean differences being below some threshold.
	\item Constrained-Mahalanobis randomization ("constrained md"): Corresponds to random permutations of the indicator across units, conditional on the Mahalanobis being below some threshold.
	\item Blocked Constrained-differences randomization ("blocked constrained diffs"): Corresponds to random permutations of the indicator within blocks of units, conditional on the standardized covariate mean differences being below some threshold.
	\item Blocked Constrained-Mahalanobis randomization ("blocked constrained md"): Corresponds to random permutations of the indicator within blocks of units, conditional on the Mahalanobis being below some threshold.
}
The null hypothesis is that the assignment mechanism holds. A large p-value does not prove that the assumption holds, but a small p-value implies that the assumption doesn't hold. These p-values are exact, in the sense that they only rely on permutations within the data and not asymptotic approximations.

In addition to specifying different assignment mechanisms, the user can specify two different test statistics:
\itemize{
	\item The Mahalanobis distance ("mahalanobis"). This acts as a global test statistic, and thus only one p-value is computed.
	\item The standardized covariate mean differences ("diffs"). This acts as a covariate-by-covariate test statistic, and thus a p-value for each covariate is computed.
}
}
\usage{
asIfRandTest(X, indicator,
  assignment = c("complete"),
  statistic = "mahalanobis",
  subclass = NULL, threshold = NULL,
  perms = 1000)
}
\arguments{
  \item{X}{A covariate matrix (rows correspond to subjects/units; columns correspond to covariates).}
  \item{indicator}{A vector of 1s and 0s (e.g., denoting treatment and control).}
  \item{assignment}{A vector of assignment mechanisms that the user wants to test; the user can test one assignment mechanism or multiple. The possible choices are "complete", "blocked", "constrained diffs", "constrained md", "blocked constrained diffs", and "blocked constrained md". See Description for more details on these assignment mechanisms.}
  \item{statistic}{The test statistic used in the randomization test. The choices are either "mahalanobis" (the Mahalanobis distance) or "diffs" (the standardized covariate mean differences). The former runs a global test and provides one p-value; the latter runs covariate-by-covariate tests and provides a p-value for each covariate.}
  \item{subclass}{A vector denoting the subclass/block for each subject/unit. This must be specified only if one of the blocked assignment mechanisms are used.}
  \item{threshold}{The threshold used within the constrained assignment mechanisms; thus, this must be specified only if one of the constrained assignment mechanisms are used. This can be a single number or a vector of numbers (e.g., if one wants to use a different threshold for each covariate when testing constrained-differences randomization).}
  \item{perms}{The number of permutations used within the randomization test. A larger number requires more computation time but results in a more consistent p-value.}
}
\value{
p-values assessing as-if randomization of an indicator for different assignment mechanisms. If the Mahalanobis distance is used as a test statistic, then a vector of p-values is reported is reported (one for each assignment mechanism). If the standardized covariate mean differences are used as a test statistic, then a table of p-values is reported, where the rows correspond to assignment mechanisms and the columns correspond to covariates.
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
	
	#Note that the following examples only use 100 permutations
	#to approximate the randomization distribution.
	#In practice, we recommend setting perms = 1000 or more;
	#in these examples we use perms = 50 to save computation time.

	#testing complete randomization for the full dataset
	#using the Mahalanobis distance.
	#We reject complete randomization in this test.
	# asIfRandTest(X = X.lalonde, indicator = indicator.lalonde)
	#testing complete randomization for the full dataset
	#using standardized covariate mean differences.
	#We reject complete randomization for most covariates:
	asIfRandTest(X = X.lalonde, indicator = indicator.lalonde,
  		statistic = "diffs",
  		perms = 50)

	#testing complete randomization and block (paired) randomization
	#for the propensity score matched dataset
	#using the Mahalanobis distance.
	#We reject both assignment mechanisms in this test.
	asIfRandTest(X = X.matched.ps, indicator = indicator.matched.ps,
  		assignment = c("complete", "blocked"),
		subclass = lalonde.matched.ps$subclass,
		perms = 50)
	#testing complete randomization and block (paired) randomization
	#for the propensity score matched dataset
	#using the standardized covariate mean differences.
	#We reject these assignment mechanisms for
	#the race covariates (hispan and black):
	asIfRandTest(X = X.matched.ps, indicator = indicator.matched.ps,
  		assignment = c("complete", "blocked"),
  		subclass = lalonde.matched.ps$subclass,
		statistic = "diffs",
		perms = 50)

	#testing three assignment mechanisms for
	#the cardinality matched dataset:
	# 1) complete randomization
	# 2) blocked (paired) randomization
	# 3) constrained-MD randomization
	#Note that the Mahalanobis distance is approximately a chi^2_K distribution,
	#where K is the number of covariates. In the Lalonde data, K = 8.
	#Thus, the threshold can be chosen as the quantile of the chi^2_8 distribution.
	#This threshold constrains the Mahalanobis distance to be below the 25-percent quantile:
	a = qchisq(p = 0.25, df = 8)
	#First we'll run the test using the Mahalanobis distance.
	#We fail to reject for the first two assignment mechanisms,
	#but reject the third.
	asIfRandTest(X = X.matched.card, indicator = indicator.matched.card,
  		assignment = c("complete", "blocked", "constrained md"),
  		subclass = lalonde.matched.card$subclass,
  		threshold = a,
  		perms = 50)
  	#Now we'll run the test using the standardized covariate mean differences.
  	#Interestingly, you fail to reject for all three assignment mechanisms
  	#for all covariates:
  	asIfRandTest(X = X.matched.card, indicator = indicator.matched.card,
  		assignment = c("complete", "blocked", "constrained md"),
  		subclass = lalonde.matched.card$subclass,
  		threshold = a,
  		statistic = "diffs",
  		perms = 50)
}

\author{Zach Branson
}


\keyword{randomization test}
\keyword{complete randomization}
\keyword{block randomization}
\keyword{constrained randomization}
\keyword{Mahalanobis distance}
