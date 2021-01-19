#Internal function for permuting an indicator within a subclass.
#This function takes a table(subclass, indicator) object.
getBlockPerm = function(subclassIndicatorTable){
  #number of subclasses
  Ns = nrow(subclassIndicatorTable)
  #total number of units
  N = sum(subclassIndicatorTable)
  #for each subclass, create a permuted indicator,
  #according to the number of 1s and 0s in that subclass.
  permutedIndicator = vector()
  for(s in 1:Ns){
    permutedIndicator.s = sample(c(rep(0, subclassIndicatorTable[s,"0"]), rep(1, subclassIndicatorTable[s,"1"])))
    permutedIndicator = append(permutedIndicator, permutedIndicator.s)
  }
  return(permutedIndicator)
}

#Internal function that permutes an indicator until
#all absolute standardized covariate mean differences are below some threshold.
permuteData.constrainedStandardizedCovMeanDiffs = function(X, indicator, threshold){
  #First, permute the indicator
  indicator = sample(indicator)

  #Now check if the absolute standardized covariate mean diffs are less than the threshold
  standardizedCovMeanDiff.obs = getStandardizedCovMeanDiffs(X = X, indicator = indicator)
  #while there are absolute standardized covariate mean differences bigger than the threshold...
  while( sum( abs(standardizedCovMeanDiff.obs) > threshold) > 0 ){
    indicator = sample(indicator)
    standardizedCovMeanDiff.obs = getStandardizedCovMeanDiffs(X = X, indicator = indicator)
  }
  return(indicator)
}
#Internal function that permutes an indicator within blocks until
#all absolute standardized covariate mean differences are below some threshold.
permuteData.blocked.constrainedStandardizedCovMeanDiffs = function(X, indicator, subclass, threshold){
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  indicator = data$indicator
  subclass = data$subclass

  #To efficiently get block permutations,
  #we just need a table of the indicator and subclass
  subclassIndicatorTable = table(subclass, indicator)

  #Then, the block permutation is
  indicator = getBlockPerm(subclassIndicatorTable)

  #Now check if the absolute standardized covariate mean diffs are less than the threshold
  standardizedCovMeanDiff.obs = getStandardizedCovMeanDiffs(X = X, indicator = indicator)
  #while there are absolute standardized covariate mean differences bigger than the threshold...
  while( sum( abs(standardizedCovMeanDiff.obs) > threshold) > 0 ){
    indicator = getBlockPerm(subclassIndicatorTable)
    standardizedCovMeanDiff.obs = getStandardizedCovMeanDiffs(X = X, indicator = indicator)
  }
  return(indicator)
}

#Internal function that permutes an indicator until
#the Mahalanobis distance is below some threshold.
permuteData.constrainedMD = function(X, indicator, threshold){
  #First, permute the indicator
  indicator = sample(indicator)

  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #Now check if the MD is less than the threshold
  md.obs = getMD(X = X, indicator = indicator, covX.inv = covX.inv)

  #while the MD is bigger than the threshold...
  while( md.obs > threshold ){
    indicator = sample(indicator)
    md.obs = getMD(X = X, indicator = indicator, covX.inv = covX.inv)
  }
  return(indicator)
}

#Internal function that permutes an indicator within blocks until
#the Mahalanobis distance is below some threshold.
permuteData.blocked.constrainedMD = function(X, indicator, subclass, threshold){
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  indicator = data$indicator
  subclass = data$subclass

  #To efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))


  #To efficiently get block permutations,
  #we just need a table of the indicator and subclass
  subclassIndicatorTable = table(subclass, indicator)

  #Then, the block permutation is
  indicator = getBlockPerm(subclassIndicatorTable)

  #Now check if the MD is less than the threshold
  md.obs = getMD(X = X, indicator = indicator, covX.inv = covX.inv)

  #while the MD is bigger than the threshold...
  while( md.obs > threshold ){
    indicator = getBlockPerm(subclassIndicatorTable)
    md.obs = getMD(X = X, indicator = indicator, covX.inv = covX.inv)
  }
  return(indicator)
}

#Internal function that returns the standardized covariate mean differences
#across many permutations of an indicator.
getCompletePerms.balance = function(X, indicator, perms = 1000){
  #number of covariates
  K = ncol(X)

  #observed standardized covariate mean differences
  standardizedCovMeanDiff.obs = getStandardizedCovMeanDiffs(X, indicator)

  #permutations for the randomization test
  indicator.permutations = replicate(perms, sample(indicator), simplify = FALSE)
  #compute the vector of covariate mean differences for each permutation
  permutations.standardizedCovMeanDiffs = matrix(nrow = perms, ncol = K)
  for(i in 1:perms){
    permutations.standardizedCovMeanDiffs[i,] = getStandardizedCovMeanDiffs(X, indicator.permutations[[i]])
  }
  return(permutations.standardizedCovMeanDiffs)
}
#Internal function that returns the standardized covariate mean differences
#across many block permutations of an indicator within a subclass.
getBlockPerms.balance = function(X, indicator, subclass, perms = 1000){
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  indicator = data$indicator
  subclass = data$subclass

  #To efficiently get block permutations,
  #we just need a table of the indicator and subclass
  subclassIndicatorTable = table(subclass, indicator)

  #Then, the set of block permutations is
  indicator.permutations = t(replicate(perms,
    getBlockPerm(subclassIndicatorTable), simplify = TRUE))
  #compute the vector of covariate mean differences for each permutation
  permutations.covMeanDiffs = matrix(nrow = perms, ncol = ncol(X))
  for(i in 1:perms){
    permutations.covMeanDiffs[i,] = getStandardizedCovMeanDiffs(X, indicator = indicator.permutations[i,])
  }
  return(permutations.covMeanDiffs)
}
#Internal function that returns the standardized covariate mean differences
#across many constrained permutations of an indicator.
#(Here, we're constraining the standardized covariate mean differences.)
getConstrainedDiffsPerms.balance = function(X, indicator, threshold, perms = 1000){
  #number of covariates
  K = ncol(X)

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.constrainedStandardizedCovMeanDiffs(X = X, indicator = indicator, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.standardizedCovMeanDiffs = matrix(nrow = perms, ncol = K)
  for(i in 1:perms){
    permutations.standardizedCovMeanDiffs[i,] = getStandardizedCovMeanDiffs(X, indicator.permutations[[i]])
  }
  return(permutations.standardizedCovMeanDiffs)
}
#Internal function that returns the standardized covariate mean differences
#across many constrained blocked permutations of an indicator.
#(Here, we're constraining the standardized covariate mean differences.)
getConstrainedDiffsBlockedPerms.balance = function(X, indicator, subclass, threshold, perms = 1000){
  #number of covariates
  K = ncol(X)

  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  indicator = data$indicator
  subclass = data$subclass

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.blocked.constrainedStandardizedCovMeanDiffs(X = X, indicator = indicator, subclass = subclass, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.standardizedCovMeanDiffs = matrix(nrow = perms, ncol = K)
  for(i in 1:perms){
    permutations.standardizedCovMeanDiffs[i,] = getStandardizedCovMeanDiffs(X, indicator.permutations[[i]])
  }
  return(permutations.standardizedCovMeanDiffs)
}
#Internal function that returns the standardized covariate mean differences
#across many constrained permutations of an indicator.
#(Here, we're constraining the Mahalanobis distance.)
getConstrainedMDPerms.balance = function(X, indicator, threshold, perms = 1000){
  #number of covariates
  K = ncol(X)

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.constrainedMD(X = X, indicator = indicator, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.standardizedCovMeanDiffs = matrix(nrow = perms, ncol = K)
  for(i in 1:perms){
    permutations.standardizedCovMeanDiffs[i,] = getStandardizedCovMeanDiffs(X, indicator.permutations[[i]])
  }
  return(permutations.standardizedCovMeanDiffs)
}
getConstrainedMDBlockedPerms.balance = function(X, indicator, subclass, threshold, perms = 1000){
  #number of covariates
  K = ncol(X)

  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  indicator = data$indicator
  subclass = data$subclass

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.blocked.constrainedMD(X = X, indicator = indicator, subclass = subclass, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.standardizedCovMeanDiffs = matrix(nrow = perms, ncol = K)
  for(i in 1:perms){
    permutations.standardizedCovMeanDiffs[i,] = getStandardizedCovMeanDiffs(X, indicator.permutations[[i]])
  }
  return(permutations.standardizedCovMeanDiffs)
}

#Internal function that returns the Mahalanobis distance (MD)
#across many permutations of an indicator.
getCompletePerms.md = function(X, indicator, perms = 1000){
  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #permutations for the randomization test
  indicator.permutations = replicate(perms, sample(indicator), simplify = FALSE)
  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[[i]], covX.inv)
  }
  return(permutations.md)
}

#Internal function that returns the Mahalanobis distance (MD)
#across many block permutations of an indicator within a subclass.
getBlockPerms.md = function(X, indicator, subclass, perms = 1000){
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  covX.inv = solve(as.matrix(stats::cov(X)))
  indicator = data$indicator
  subclass = data$subclass

  #To efficiently get block permutations,
  #we just need a table of the indicator and subclass
  subclassIndicatorTable = table(subclass, indicator)

  #Then, the set of block permutations is
  indicator.permutations = t(replicate(perms,
    getBlockPerm(subclassIndicatorTable), simplify = TRUE))
  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[i,], covX.inv)
  }
  return(permutations.md)
}

#Internal function that returns the Mahalanobis distance
#across many constrained permutations of an indicator.
#(Here, we're constraining the standardized covariate mean differences.)
getConstrainedDiffsPerms.md = function(X, indicator, threshold, perms = 1000){
  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.constrainedStandardizedCovMeanDiffs(X = X, indicator = indicator, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[[i]], covX.inv)
  }
  return(permutations.md)
}
#Internal function that returns the Mahalanobis distance
#across many constrained permutations of an indicator.
#(Here, we're constraining the Mahalanobis distance.)
getConstrainedMDPerms.md = function(X, indicator, threshold, perms = 1000){
  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.constrainedMD(X = X, indicator = indicator, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[[i]], covX.inv)
  }
  return(permutations.md)
}

#Internal function that returns the Mahalanobis distance
#across many constrained blocked permutations of an indicator.
#(Here, we're constraining the standardized covariate mean differences.)
getConstrainedDiffsBlockedPerms.md = function(X, indicator, subclass, threshold, perms = 1000){
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  indicator = data$indicator
  subclass = data$subclass

  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.blocked.constrainedStandardizedCovMeanDiffs(X = X, indicator = indicator, subclass = subclass, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[[i]], covX.inv)
  }
  return(permutations.md)
}

#Internal function that returns the Mahalanobis distance
#across many constrained blocked permutations of an indicator.
#(Here, we're constraining the Mahalanobis distance.)
getConstrainedMDBlockedPerms.md = function(X, indicator, subclass, threshold, perms = 1000){
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  indicator = data$indicator
  subclass = data$subclass

  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #first, collect the indicators that fulfill the balance constraint
  indicator.permutations = replicate(perms,
    permuteData.blocked.constrainedMD(X = X, indicator = indicator, subclass = subclass, threshold = threshold),
    simplify = FALSE)

  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[[i]], covX.inv)
  }
  return(permutations.md)
}
#internal function to associate a color with each assignment mechanism
#(used only for plotting, for consistency)
getAssignmentColor = function(assignment){
  #the different assignment colors are:
  # complete: black
  # blocked: blue
  # constrained diffs: green
  # constrained md: red
  assignment.colors = vector(length = length(assignment))
  for(i in 1:length(assignment.colors)){
    if(assignment[i] == "complete"){
      assignment.colors[i] = "black"
    }
    if(assignment[i] == "blocked"){
      assignment.colors[i] = "blue"
    }
    if(assignment[i] == "constrained diffs"){
      assignment.colors[i] = "green"
    }
    if(assignment[i] == "constrained md"){
      assignment.colors[i] = "red"
    }
    if(assignment[i] == "blocked constrained diffs"){
      assignment.colors[i] = "orange"
    }
    if(assignment[i] == "blocked constrained md"){
      assignment.colors[i] = "purple"
    }
  }
  return(assignment.colors)
}






