\name{randChecks-internal}
\title{Internal \code{randChecks} functions}

\alias{getBlockPerm}
\alias{permuteData.constrainedStandardizedCovMeanDiffs}
\alias{permuteData.blocked.constrainedStandardizedCovMeanDiffs}
\alias{permuteData.constrainedMD}
\alias{permuteData.blocked.constrainedMD}
\alias{getCompletePerms.balance}
\alias{getBlockPerms.balance}
\alias{getConstrainedDiffsPerms.balance}
\alias{getConstrainedDiffsBlockedPerms.balance}
\alias{getConstrainedMDPerms.balance}
\alias{getConstrainedMDBlockedPerms.balance}
\alias{getCompletePerms.md}
\alias{getBlockPerms.md}
\alias{getConstrainedDiffsPerms.md}
\alias{getConstrainedMDPerms.md}
\alias{getConstrainedDiffsBlockedPerms.md}
\alias{getConstrainedMDBlockedPerms.md}
\alias{getAssignmentColor}

\description{Internal \code{randChecks} functions}
\usage{
  getBlockPerm(subclassIndicatorTable)
  permuteData.constrainedStandardizedCovMeanDiffs(X, indicator, threshold)
  permuteData.blocked.constrainedStandardizedCovMeanDiffs(X,indicator,subclass,threshold)
  permuteData.constrainedMD(X, indicator, threshold)
  permuteData.blocked.constrainedMD(X, indicator, subclass, threshold)
  getCompletePerms.balance(X, indicator, perms = 1000)
  getBlockPerms.balance(X, indicator, subclass, perms = 1000)
  getConstrainedDiffsPerms.balance(X, indicator, threshold, perms = 1000)
  getConstrainedDiffsBlockedPerms.balance(X, indicator, subclass, threshold, perms = 1000)
  getConstrainedMDPerms.balance(X, indicator, threshold, perms = 1000)
  getConstrainedMDBlockedPerms.balance(X, indicator, subclass, threshold, perms = 1000)
  getCompletePerms.md(X, indicator, perms = 1000)
  getBlockPerms.md(X, indicator, subclass, perms = 1000)
  getConstrainedDiffsPerms.md(X, indicator, threshold, perms = 1000)
  getConstrainedMDPerms.md(X, indicator, threshold, perms = 1000)
  getConstrainedDiffsBlockedPerms.md(X, indicator, subclass, threshold, perms = 1000)
  getConstrainedMDBlockedPerms.md(X, indicator, subclass, threshold, perms = 1000)
  getAssignmentColor(assignment)
}
\author{Zach Branson}
\details{
  	\code{getBlockPerm} is a helper function to permute an indicator within blocks (defined by a subclass). 
  	
  	\code{permuteData.constrainedStandardizedCovMeanDiffs} permutes an indicator given a constraint on standardized covariate mean differences.
  	
	\code{permuteData.blocked.constrainedStandardizedCovMeanDiffs} permutes an indicator within blocks given a constraint on standardized covariate mean differences.
	
	\code{permuteData.constrainedMD} permutes an indicator given a constraint on the Mahalanobis distance.
	
	\code{permuteData.blocked.constrainedMD} permutes an indicator within blocks given a constraint on the Mahalanobis distance.
	
	\code{getCompletePerms.balance} returns standardized covariate mean differences across many permutations.
	
	\code{getBlockPerms.balance} returns standardized covariate mean differences across many permutations within blocks.
	
	\code{getConstrainedDiffsPerms.balance} returns standardized covariate mean differences across many permutations, conditional on a constraint on the standardized covariate mean differences.
	
	\code{getConstrainedDiffsBlockedPerms.balance} returns standardized covariate mean differences across many permutations within blocks, conditional on a constraint on the standardized covariate mean differences.
	
	\code{getConstrainedMDPerms.balance} returns standardized covariate mean differences across many permutations, conditional on a constraint on the Mahalanobis distance.
	
	\code{getConstrainedMDBlockedPerms.balance} returns standardized covariate mean differences across many permutations within blcoks, conditional on a constraint on the Mahalanobis distance.
	
	\code{getCompletePerms.md} returns the Mahalanobis distance across many permutations.
	
	\code{getBlockPerms.md} returns the Mahalanobis distance across many permutations within blocks.
	
	\code{getConstrainedDiffsPerms.md} returns the Mahalanobis distance across many permutations, conditional on a constraint on the standardized covariate mean differences.
	
	\code{getConstrainedMDPerms.md} returns the Mahalanobis distance across many permutations, conditional on a constraint on the Mahalanobis distance.
	
	\code{getConstrainedDiffsBlockedPerms.md} returns the Mahalanobis distance across many permutations within blocks, conditional on a constraint on the standardized covariate mean differences.
	
	\code{getConstrainedMDBlockedPerms.md} returns the Mahalanobis distance across many permutations within blocks, conditional on a constraint on the Mahalanobis distance.
	
	\code{getAssignmentColor} is a helper function that pairs different assignment mechanisms with different colors to be used in plots within \code{randChecks}.
}
\keyword{internal}