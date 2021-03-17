#Returns the Mahalanobis distance (MD)
#according to an indicator (instrument or exposure).
#Specifically, the MD for treatment-versus-control (T or C) is defined as:
#((Nt*(N-Nt))/N)*(\bar{X}_T - \bar{X}_C)%*%[cov(X)]^{-1}%*%(\bar{X}_T - \bar{X}_C)
#Here, Nt is the number of treated units (or indicator = 1 units),
#N is the number of units,
#\bar{X}_T is the vector of covariate means for treated units,
#\bar{X}_C is the vector of covariate means for control units,
#and cov(X) is the covariate-covariance matrix.

#Note that we make cov(X) correspond to
#the covariate-covariance matrix in the full dataset,
#rather than the matched dataset.

#Sometimes it can be computationally efficient to provide
#the inverse of cov(X) to this function; that can be provided
#via the argument covX.inv
getMD = function(X.matched, indicator.matched, covX.inv = NULL,
  X.full = NULL){

  #If the covariate matrix in the full dataset is not specified,
  #just set it equal to X.matched.
  if(is.null(X.full)){
    X.full = X.matched
  }

  data = data.frame(X.matched, indicator.matched = indicator.matched)
  #treatment group
  treatmentData = subset(data, indicator.matched == 1)
  #control group
  controlData = subset(data, indicator.matched == 0)

  #now we can get rid of the indicator variable
  treatmentData = subset(treatmentData, select = -c(indicator.matched))
  controlData = subset(controlData, select = -c(indicator.matched))

  #covariate mean difference
  covMeanDiffs = colMeans(treatmentData) - colMeans(controlData)
  if(is.null(covX.inv)){
    covX.inv = solve(as.matrix(stats::cov(X.full)))
  }

  n = as.numeric(nrow(data))
  n.t = as.numeric(sum(indicator.matched))
  md = ((n.t*(n-n.t))/n)*t(covMeanDiffs)%*%covX.inv%*%covMeanDiffs
  return(md)
}