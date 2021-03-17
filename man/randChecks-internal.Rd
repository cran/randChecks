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
  permuteData.constrainedStandardizedCovMeanDiffs(X.matched, indicator.matched,
  threshold, X.full = NULL, indicator.full = NULL)
  permuteData.blocked.constrainedStandardizedCovMeanDiffs(X.matched, indicator.matched,
  subclass, threshold, X.full = NULL, indicator.full = NULL)
  permuteData.constrainedMD(X.matched, indicator.matched, threshold, X.full = NULL)
  permuteData.blocked.constrainedMD(X.matched, indicator.matched,
  subclass, threshold, X.full = NULL)
  getCompletePerms.balance(X.matched, indicator.matched,
  perms = 1000, X.full = NULL, indicator.full = NULL)
  getBlockPerms.balance(X.matched, indicator.matched,
  subclass, perms = 1000, X.full = NULL, indicator.full = NULL)
  getConstrainedDiffsPerms.balance(X.matched, indicator.matched,
  threshold, perms = 1000, X.full = NULL, indicator.full = NULL)
  getConstrainedDiffsBlockedPerms.balance(X.matched, indicator.matched,
  subclass, threshold, perms = 1000, X.full = NULL, indicator.full = NULL)
  getConstrainedMDPerms.balance(X.matched, indicator.matched,
  threshold, perms = 1000, X.full = NULL, indicator.full = NULL)
  getConstrainedMDBlockedPerms.balance(X.matched, indicator.matched,
  subclass, threshold, perms = 1000, X.full = NULL, indicator.full = NULL)
  getCompletePerms.md(X.matched, indicator.matched,
  perms = 1000, X.full = NULL)
  getBlockPerms.md(X.matched, indicator.matched,
  subclass, perms = 1000, X.full = NULL)
  getConstrainedDiffsPerms.md(X.matched, indicator.matched,
  threshold, perms = 1000, X.full = NULL, indicator.full = NULL)
  getConstrainedMDPerms.md(X.matched, indicator.matched,
  threshold, perms = 1000, X.full = NULL)
  getConstrainedDiffsBlockedPerms.md(X.matched, indicator.matched,
  subclass, threshold, perms = 1000,
  	X.full = NULL, indicator.full = NULL)
  getConstrainedMDBlockedPerms.md(X.matched, indicator.matched,
  subclass, threshold, perms = 1000, X.full = NULL)
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