#returns the standardized covariate mean differences between two groups (1 or 0), i.e.,
#returns (\bar{X}_T - \bar{X}_C)/sd(X), where
#\bar{X}_T is the vector of covariate means for treated units,
#\bar{X}_C is the vector of covariate means for control units,
#and sd(X) is the square root of the pooled variance between groups.
getStandardizedCovMeanDiffs = function (X.matched, indicator.matched,
  X.full = NULL, indicator.full = NULL) {
    #Here we standardize by the pooled variance
    # *in the full, pre-matched dataset*
    #which is the recommendation in the literature.
    #If you specify just the matched data
    #(i.e., don't specify X.full and indicator.full),
    #then we'll use the pooled variance in the matched data:
    if(is.null(X.full) | is.null(indicator.full)){
      X.full = X.matched
      indicator.full = indicator.matched
    }
    
    #put together the matched and full data
    matchedData = data.frame(X.matched, indicator.matched = indicator.matched)
    fullData = data.frame(X.full, indicator.matched = indicator.full)

    #the treatment and control data for the matched data
    treatmentData.matched = subset(matchedData, indicator.matched == 1)
    controlData.matched = subset(matchedData, indicator.matched == 0)
    treatmentData.matched = subset(treatmentData.matched, select = -c(indicator.matched))
    controlData.matched = subset(controlData.matched, select = -c(indicator.matched))

    #the treatment and control data for the full data
    treatmentData.full = subset(fullData, indicator.matched == 1)
    controlData.full = subset(fullData, indicator.matched == 0)
    treatmentData.full = subset(treatmentData.full, select = -c(indicator.matched))
    controlData.full = subset(controlData.full, select = -c(indicator.matched))

    #standardization before matching
    cov.variances.treatment = apply(treatmentData.full,
      MARGIN = 2, FUN = stats::var)
    cov.variances.control = apply(controlData.full,
      MARGIN = 2, FUN = stats::var)
    pooled.cov.variance = (cov.variances.treatment + cov.variances.control)/2
    #mean difference in matched data
    covMeanDiffs = colMeans(treatmentData.matched) - colMeans(controlData.matched)
    #standardized covariate mean differences
    covMeanDiffs.standardized = covMeanDiffs/sqrt(pooled.cov.variance)
    return(covMeanDiffs.standardized)
}