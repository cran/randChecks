#Return an exact p-value testing whether an indicator
#follows a given assignment mechanism.
#Currently, this function supports the following assignment mechanisms:
#-complete randomization ("complete")
#-block randomization ("blocked")
#-constrained MD randomization ("constrained md")
#-constrained standardized diffs randomization ("constrained diffs")
#-block constrained MD randomization ("blocked constrained md")
#-block constrained standardized diffs randomization ("blocked constrained diffs")
#Currently, the function supports the following test statistics:
#-Mahalanobis distance (one p-value as a global test)
#-standardized covariate mean differences (covariate-by-covariate p-values)
asIfRandTest = function(X,
  indicator,
  assignment = c("complete"),
  statistic = "mahalanobis",
  subclass = NULL, threshold = NULL,
  perms = 1000){
  #must provide a covariate matrix
  if(is.null(X)){
    stop("Error: Must provide a covariate matrix X.")
  }
  #must provide an indicator
  if(is.null(indicator)){
    stop("Error: Must provide an indicator of 1s and 0s.")
  }
  #must provide an assignment mechanism
  if(is.null(assignment)){
    stop("Error: Must provide at least one assignment mechanism. 
      Mechanisms available: complete, blocked, constrained md, constrained diffs")
  }

  #If you pick "blocked", you need to provide a
  #subclass vector that specifies the blocks
  if(("blocked" %in% assignment |
  "blocked constrained md" %in% assignment |
  "blocked constrained diffs" %in% assignment) & is.null(subclass)){
    stop("Error: For blocked assignments, must provide a subclass vector.")
  }
  #If you pick "constrained diffs" or "constrained md",
  #you need to provide a threshold that specifies the balance constraint.
  if(("constrained diffs" %in% assignment | "constrained md" %in% assignment |
    "blocked constrained diffs" %in% assignment |
    "blocked constrained md" %in% assignment) & is.null(threshold)){
    stop("Error: For constrained assignments, must provide a threshold.")
  }
  #assignment must be one of the following:
  #-complete
  #-blocked
  #-constrained md
  #-constrained diffs
  #-blocked constrained md
  #-blocked constrained diffs
  if(sum( assignment %in% c("complete", "blocked", "constrained md", "constrained diffs", "blocked constrained md", "blocked constrained diffs") ) != length(assignment) ){
    stop('Error: Available assignments are "complete", "blocked", "constrained md", "constrained diffs", "blocked constrained md", "blocked constrained diffs"')
  }

  #the user must choose either the mahalanobis distance
  #or standardized covariate mean differences
  #as the test statistic
  if(statistic != "mahalanobis" & statistic != "diffs"){
    stop('Error: Must set statistic = "mahalanobis" (global test) or statistic = "diffs" (covariate-by-covariate test)')
  }

  #if the user chose the mahalanobis distance:
  if(statistic == "mahalanobis"){
    #vector of MDs across different assignment mechanisms
    compPerms = vector()
    blockPerms = vector()
    constrainedDiffsPerms = vector()
    constrainedMDPerms = vector()
    blockedConstrainedDiffsPerms = vector()
    blockedConstrainedMDPerms = vector()
    #define these vectors, depending on which
    #assignment mechanisms the user selects.
    if("complete" %in% assignment){
      compPerms = getCompletePerms.md(X = X, indicator = indicator, perms = perms)
    }
    if("blocked" %in% assignment){
      blockPerms = getBlockPerms.md(X = X, indicator = indicator, subclass = subclass, perms = perms)
    }
    if("constrained diffs" %in% assignment){
      constrainedDiffsPerms = getConstrainedDiffsPerms.md(X = X, indicator = indicator, threshold = threshold, perms = perms)
    }
    if("constrained md" %in% assignment){
      constrainedMDPerms = getConstrainedMDPerms.md(X = X, indicator = indicator, threshold = threshold, perms = perms)
    }
    if("blocked constrained diffs" %in% assignment){
      blockedConstrainedDiffsPerms = getConstrainedDiffsBlockedPerms.md(X = X, indicator = indicator, subclass = subclass, threshold = threshold, perms = perms)
    }
    if("blocked constrained md" %in% assignment){
      blockedConstrainedMDPerms = getConstrainedMDBlockedPerms.md(X = X, indicator = indicator, subclass = subclass, threshold = threshold, perms = perms)
    }
    #the observed MD is
    md.obs = as.numeric(getMD(X, indicator = indicator))
    #compute p-values for different assignment mechanisms
    assignment.pvalues = vector(length = length(assignment))
    for(i in 1:length(assignment.pvalues)){
      if(assignment[i] == "complete"){
        assignment.pvalues[i] = mean( compPerms > md.obs )
      }
      if(assignment[i] == "blocked"){
        assignment.pvalues[i] = mean( blockPerms > md.obs )
      }
      if(assignment[i] == "constrained diffs"){
        assignment.pvalues[i] = mean( constrainedDiffsPerms > md.obs )
      }
      if(assignment[i] == "constrained md"){
        assignment.pvalues[i] = mean( constrainedMDPerms > md.obs )
      }
      if(assignment[i] == "blocked constrained diffs"){
        assignment.pvalues[i] = mean( blockedConstrainedDiffsPerms > md.obs )
      }
      if(assignment[i] == "blocked constrained md"){
        assignment.pvalues[i] = mean( blockedConstrainedMDPerms > md.obs )
      }
    }
    #then, the table of p-values is:
    assignment.pvalues = as.matrix(assignment.pvalues)
    rownames(assignment.pvalues) = assignment
    colnames(assignment.pvalues) = "pvalues"
    return(assignment.pvalues)
  }
  #if the user chose the standardized covariate mean differences:
  if(statistic == "diffs"){
    #the number of covariates is
    K = ncol(X)

    #matrix of standardized covariate mean diffs
    #across different assignment mechanisms
    #The rows correspond to permutations,
    #the columns correspond to covariates.
    compPerms = matrix(nrow = perms, ncol = K)
    blockPerms = matrix(nrow = perms, ncol = K)
    constrainedDiffsPerms = matrix(nrow = perms, ncol = K)
    constrainedMDPerms = matrix(nrow = perms, ncol = K)
    blockedConstrainedDiffsPerms = matrix(nrow = perms, ncol = K)
    blockedConstrainedMDPerms = matrix(nrow = perms, ncol = K)
    #define these matrices, depending on which
    #assignment mechanisms the user selects.
    if("complete" %in% assignment){
      compPerms = abs(getCompletePerms.balance(X = X, indicator = indicator, perms = perms))
    }
    if("blocked" %in% assignment){
      blockPerms = abs(getBlockPerms.balance(X = X, indicator = indicator, subclass = subclass, perms = perms))
    }
    if("constrained diffs" %in% assignment){
      constrainedDiffsPerms = abs(getConstrainedDiffsPerms.balance(X = X, indicator = indicator, threshold = threshold, perms = perms))
    }
    if("constrained md" %in% assignment){
      constrainedMDPerms = abs(getConstrainedMDPerms.balance(X = X, indicator = indicator, threshold = threshold, perms = perms))
    }
    if("blocked constrained diffs" %in% assignment){
      blockedConstrainedDiffsPerms = abs(getConstrainedDiffsBlockedPerms.balance(X = X, indicator = indicator, subclass = subclass, threshold = threshold, perms = perms))
    }
    if("blocked constrained md" %in% assignment){
      blockedConstrainedMDPerms = abs(getConstrainedMDBlockedPerms.balance(X = X, indicator = indicator, subclass = subclass, threshold = threshold, perms = perms))
    }
    #the observed standardized covariate mean differences are
    diffs.obs = abs(as.numeric(getStandardizedCovMeanDiffs(X, indicator = indicator)))
    #compute p-values for different assignment mechanisms
    #Each row corresponds to an assignment mechanism
    #Each column corresponds to a covariate.
    assignment.pvalues = matrix(nrow = length(assignment), ncol = K)
    for(i in 1:nrow(assignment.pvalues)){
      if(assignment[i] == "complete"){
        assignment.pvalues[i,] = colMeans(t(apply(compPerms, MARGIN = 1, FUN = function(x) x > diffs.obs)))
      }
      if(assignment[i] == "blocked"){
        assignment.pvalues[i,] = colMeans(t(apply(blockPerms, MARGIN = 1, FUN = function(x) x > diffs.obs)))
      }
      if(assignment[i] == "constrained diffs"){
        assignment.pvalues[i,] = colMeans(t(apply(constrainedDiffsPerms, MARGIN = 1, FUN = function(x) x > diffs.obs)))
      }
      if(assignment[i] == "constrained md"){
        assignment.pvalues[i,] = colMeans(t(apply(constrainedMDPerms, MARGIN = 1, FUN = function(x) x > diffs.obs)))
      }
      if(assignment[i] == "blocked constrained diffs"){
        assignment.pvalues[i,] = colMeans(t(apply(blockedConstrainedDiffsPerms, MARGIN = 1, FUN = function(x) x > diffs.obs)))
      }
      if(assignment[i] == "blocked constrained md"){
        assignment.pvalues[i,] = colMeans(t(apply(blockedConstrainedMDPerms, MARGIN = 1, FUN = function(x) x > diffs.obs)))
      }
    }
    #then, the table of p-values is:
    rownames(assignment.pvalues) = assignment
    colnames(assignment.pvalues) = colnames(X)
    return(assignment.pvalues)
  }
}