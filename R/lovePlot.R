#create a Love plot of the standardized covariate mean differences
#across an indicator (usually a treatment)
#Can also display the permutation quantiles for these quantities.
#Currently, this function supports the following assignment mechanisms:
#-complete randomization ("complete")
#-block randomization ("blocked")
#-constrained MD randomization ("constrained md")
#-constrained standardized diffs randomization ("constrained diffs")
#-blocked constrained MD randomization ("blocked constrained md")
#-blocked constrained standardized diffs randomization ("blocked constrained diffs")
lovePlot = function(X.matched, indicator.matched,
  permQuantiles = FALSE,
  assignment = "complete",
  subclass = NULL, threshold = NULL,
  alpha = 0.15, perms = 1000,
  X.full = NULL, indicator.full = NULL){
  #must provide a covariate matrix
  if(is.null(X.matched)){
    stop("Error: Must provide a covariate matrix X.matched.")
  }
  #must provide an indicator
  if(is.null(indicator.matched)){
    stop("Error: Must provide an indicator (indicator.matched) of 1s and 0s.")
  }
  #assignments must be one of the following:
  #-complete
  #-blocked
  #-constrained md
  #-constrained diffs
  #-blocked constrained md
  #-blocked constrained diffs
  if(assignment != "complete" & assignment != "blocked" &
    assignment != "constrained md" & assignment != "constrained diffs" &
    assignment != "blocked constrained md" &
    assignment != "blocked constrained diffs"){
    stop('Error: Available assignemnts are "complete", "blocked", "constrained md", "constrained diffs", "blocked constrained md", "blocked constrained diffs"')
  }
  #If you pick "blocked", you need to provide a
  #subclass vector that specifies the blocks
  if((assignment == "blocked" |
  assignment == "blocked constrained md" |
  assignment == "blocked constrained diffs") & is.null(subclass)){
    stop("Error: For blocked assignments, must provide a subclass vector.")
  }
  #If you pick "constrained diffs" or "constrained md",
  #you need to provide a threshold that specifies the balance constraint.
  if((assignment == "constrained diffs" | assignment == "constrained md" |
    assignment == "blocked constrained diffs" |
    assignment == "blocked constrained md") & is.null(threshold)){
    stop("Error: For constrained assignments, must provide a threshold.")
  }

  #compute standardized covariate mean differences
  covMeanDiffs = getStandardizedCovMeanDiffs(
    X.matched = X.matched, indicator.matched = indicator.matched,
    X.full = X.full, indicator.full = indicator.full)

  #number of covariates
  K = ncol(X.matched)

  #get range of covMeanDiffs for plot limits
  plot.min = min( covMeanDiffs )
  plot.max = max( covMeanDiffs )

  #first, compute the permutation quantiles (if desired)
  if(permQuantiles == TRUE & assignment == "complete"){
    permutations.covMeanDiffs = getCompletePerms.balance(
      X.matched = X.matched, indicator.matched = indicator.matched,
      perms = perms,
      X.full = X.full, indicator.full = indicator.full)
  }
  if(permQuantiles == TRUE & assignment == "blocked"){
    permutations.covMeanDiffs = getBlockPerms.balance(
      X.matched = X.matched, indicator.matched = indicator.matched,
      subclass = subclass, perms = perms,
      X.full = X.full, indicator.full = indicator.full)
  }
  if(permQuantiles == TRUE & assignment == "constrained diffs"){
    permutations.covMeanDiffs = getConstrainedDiffsPerms.balance(
      X.matched = X.matched, indicator.matched = indicator.matched,
      threshold = threshold, perms = perms,
      X.full = X.full, indicator.full = indicator.full)
  }
  if(permQuantiles == TRUE & assignment == "blocked constrained diffs"){
    permutations.covMeanDiffs = getConstrainedDiffsBlockedPerms.balance(
      X.matched = X.matched, indicator.matched = indicator.matched,
      threshold = threshold, subclass = subclass, perms = perms,
      X.full = X.full, indicator.full = indicator.full)
  }
  if(permQuantiles == TRUE & assignment == "constrained md"){
    permutations.covMeanDiffs = getConstrainedMDPerms.balance(
      X.matched = X.matched, indicator.matched = indicator.matched,
      threshold = threshold, perms = perms,
      X.full = X.full, indicator.full = indicator.full)
  }
  if(permQuantiles == TRUE & assignment == "blocked constrained md"){
    permutations.covMeanDiffs = getConstrainedMDBlockedPerms.balance(
      X.matched = X.matched, indicator.matched = indicator.matched,
      threshold = threshold, subclass = subclass, perms = perms,
      X.full = X.full, indicator.full = indicator.full)
  }
  
  if(permQuantiles == TRUE){
    #the quantiles are
    permutations.covMeanDiffs.lowerQuantile = apply(permutations.covMeanDiffs, MARGIN = 2, FUN = stats::quantile, probs = alpha/2)
    permutations.covMeanDiffs.upperQuantile = apply(permutations.covMeanDiffs, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha/2)

    #the permutation quantiles may be larger than the observed covariate mean differences,
    #thus we need to change the plot limits
    plot.min = min( covMeanDiffs, permutations.covMeanDiffs.lowerQuantile )
    plot.max = max( covMeanDiffs, permutations.covMeanDiffs.upperQuantile )
  }

  #make love plot
  graphics::plot(covMeanDiffs, 1:K,
    xlim = c(plot.min, plot.max), xlab = "Standardized Covariate Mean Differences", ylab = "", main = "", yaxt = "n",
    pch = 16)
  graphics::abline(v = 0, col = "gray")
  graphics::axis(side=2, at=1:K, labels = names(covMeanDiffs), las = 1, cex.axis = 0.5)
  #add the quantile lines to the plot (if desired)
  if(permQuantiles == TRUE){
    #color by the assignment mechanism
    assignmentColor = getAssignmentColor(assignment)
    graphics::lines(permutations.covMeanDiffs.lowerQuantile, 1:K, lty = 2, col = assignmentColor)
    graphics::lines(permutations.covMeanDiffs.upperQuantile, 1:K, lty = 2, col = assignmentColor)
  }
  


}