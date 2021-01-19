#create a Love plot of the standardized covariate mean differences
#across an indicator (usually a treatment) for two different datasets.
#X1 and X2 denote the covariate matrices for these two datasets
#indicator1 and indicator2 denote the indicators for these two datasets.
#dataNames denotes the names the user wants to give these datasets.
lovePlotCompare = function(X1, indicator1, X2, indicator2, dataNames = c("Dataset1", "Dataset2")){
  #must provide covariate matrices
  if(is.null(X1) | is.null(X2)){
    stop("Error: Must provide a covariate matrix X1 and a covariate matrix X2.")
  }
  #must provide indicators
  if(is.null(indicator1) | is.null(indicator2)){
    stop("Error: Must provide two indicators of 1s and 0s: indicator1 and indicator2.")
  }
  #the column dimensions of the covariate matrices must be the same
  if(ncol(X1) != ncol(X2)){
    stop("Error: X1 and X2 must have the same number of columns (and correspond to the same covariates).")
  }

  #compute standardized covariate mean differences
  covMeanDiffs1 = getStandardizedCovMeanDiffs(X = X1, indicator = indicator1)
  covMeanDiffs2 = getStandardizedCovMeanDiffs(X = X2, indicator = indicator2)

  #number of covariates
  K = ncol(X1)

  #get range of covMeanDiffs for plot limits
  plot.min = min( c(covMeanDiffs1, covMeanDiffs2) )
  plot.max = max( c(covMeanDiffs1, covMeanDiffs2) )

  #be sure that plot settings are reset after function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  #layout for plot
  m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)

  graphics::layout(mat = m, heights = c(1,0.25))
  graphics::par(mar=c(4.1,4.1,4.1,1.1))

  #make love plot
  #covariate mean differences for X1
  graphics::plot(covMeanDiffs1, 1:K,
    xlim = c(plot.min, plot.max), xlab = "Standardized Covariate Mean Differences", ylab = "", main = "", yaxt = "n",
    pch = 16)
  #covariate mean differences for X2
  graphics::points(covMeanDiffs2, 1:K,
    pch = 17, col = "red")
  graphics::abline(v = 0, col = "gray")
  graphics::axis(side=2, at=1:K, labels = names(covMeanDiffs1), las = 1, cex.axis = 0.5)
  #legend
  graphics::par(mar=c(1.1,1.1,1.1,1.1))
  graphics::plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  graphics::legend("top",inset = 0,
  legend = dataNames,
  col = c("black", "red"),
  pch = c(16, 17),
  horiz = TRUE)
}