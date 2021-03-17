#Display the randomization distribution of the Mahalanobis distance
#across an indicator for different assignment mechanisms.
#Currently, this function supports the following assignment mechanisms:
#-complete randomization ("complete")
#-block randomization ("blocked")
#-constrained MD randomization ("constrained md")
#-constrained standardized diffs randomization ("constrained diffs")
#-blocked constrained MD randomization ("blocked constrained md")
#-blocked constrained standardized diffs randomization ("blocked constrained diffs")
asIfRandPlot = function(X.matched,
  indicator.matched,
  assignment = c("complete"),
  subclass = NULL, threshold = NULL,
  perms = 1000,
  X.full = NULL, indicator.full = NULL){
  #must provide a covariate matrix
  if(is.null(X.matched)){
    stop("Error: Must provide a covariate matrix X.matched.")
  }
  #must provide an indicator
  if(is.null(indicator.matched)){
    stop("Error: Must provide an indicator (indicator.matched) of 1s and 0s.")
  }
  #must provide an assignment mechanism
  if(is.null(assignment)){
    stop("Error: Must provide at least one assignment mechanism. 
      Mechanisms available: complete, blocked, constrained md, constrained diffs, blocked constrained md, blocked constrained diffs")
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

  #vector of MDs across different assignment mechanisms
  compPerms = vector()
  blockPerms = vector()
  constrainedDiffsPerms = vector()
  constrainedMDPerms = vector()
  blockedConstrainedDiffsPerms = vector()
  blockedConstrainedMDPerms = vector()

  #get the xlim and ylim for the plot
  plot.x.vec = vector()
  plot.y.vec = vector()
  if("complete" %in% assignment){
    compPerms = getCompletePerms.md(
      X.matched = X.matched, indicator.matched = indicator.matched,
      perms = perms, X.full = X.full)
    plot.x.vec = append(plot.x.vec, compPerms)
    plot.y.vec = append(plot.y.vec, stats::density(compPerms)$y)
  }
  if("blocked" %in% assignment){
    blockPerms = getBlockPerms.md(
      X.matched = X.matched, indicator.matched = indicator.matched,
      subclass = subclass, perms = perms,
      X.full = X.full)
    plot.x.vec = append(plot.x.vec, blockPerms)
    plot.y.vec = append(plot.y.vec, stats::density(blockPerms)$y)
  }
  if("constrained diffs" %in% assignment){
    constrainedDiffsPerms = getConstrainedDiffsPerms.md(
      X.matched = X.matched, indicator.matched = indicator.matched,
      threshold = threshold, perms = perms,
      X.full = X.full, indicator.full = indicator.full)
    plot.x.vec = append(plot.x.vec, constrainedDiffsPerms)
    plot.y.vec = append(plot.y.vec, stats::density(constrainedDiffsPerms)$y)
  }
  if("constrained md" %in% assignment){
    constrainedMDPerms = getConstrainedMDPerms.md(
      X.matched = X.matched, indicator.matched = indicator.matched,
      threshold = threshold, perms = perms,
      X.full = X.full)
    plot.x.vec = append(plot.x.vec, constrainedMDPerms)
    plot.y.vec = append(plot.y.vec, stats::density(constrainedMDPerms)$y)
  }
  if("blocked constrained diffs" %in% assignment){
    blockedConstrainedDiffsPerms = getConstrainedDiffsBlockedPerms.md(
      X.matched = X.matched, indicator.matched = indicator.matched,
      subclass = subclass, threshold = threshold, perms = perms,
      X.full = X.full, indicator.full = indicator.full)
    plot.x.vec = append(plot.x.vec, blockedConstrainedDiffsPerms)
    plot.y.vec = append(plot.y.vec, stats::density(blockedConstrainedDiffsPerms)$y)
  }
  if("blocked constrained md" %in% assignment){
    blockedConstrainedMDPerms = getConstrainedMDBlockedPerms.md(
      X.matched = X.matched, indicator.matched = indicator.matched,
      subclass = subclass, threshold = threshold, perms = perms,
      X.full)
    plot.x.vec = append(plot.x.vec, blockedConstrainedMDPerms)
    plot.y.vec = append(plot.y.vec, stats::density(blockedConstrainedMDPerms)$y)
  }
  #to get the xlim and ylim for the plot,
  #also consider the observed MD
  md.obs = as.numeric(getMD(X.matched = X.matched, indicator.matched = indicator.matched, X.full = X.full))
  plot.x.max = max(c(plot.x.vec, md.obs))
  plot.y.max = max(plot.y.vec)
  if(plot.y.max == -Inf){plot.y.max = 1}

  #be sure that plot settings are reset after function
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  #layout for plot
  m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)

  graphics::layout(mat = m, heights = c(1,0.25))
  graphics::par(mar=c(4.1,4.1,4.1,1.1))

  #make a blank plot
  graphics::plot(0,0, col = "white", xlim = c(0, plot.x.max), ylim = c(0,plot.y.max),
    xlab = "Mahalanobis Distance", ylab = "Density", main = "")
  #add various densities
  if("complete" %in% assignment){
    graphics::lines(stats::density(compPerms))
  }
  if("blocked" %in% assignment){
    graphics::lines(stats::density(blockPerms), col = "blue")
  }
  if("constrained diffs" %in% assignment){
    graphics::lines(stats::density(constrainedDiffsPerms), col = "green")
  }
  if("constrained md" %in% assignment){
    graphics::lines(stats::density(constrainedMDPerms), col = "red")
  }
  if("blocked constrained diffs" %in% assignment){
    graphics::lines(stats::density(blockedConstrainedDiffsPerms), col = "orange")
  }
  if("blocked constrained md" %in% assignment){
    graphics::lines(stats::density(blockedConstrainedMDPerms), col = "purple")
  }
  #add the observed MD
  graphics::abline(v = md.obs, col = "gray", lty = 2)
  graphics::mtext(expression({MD^{obs}}), side = 1, line = 0.1, at = md.obs)
  #add legend
  #the colors for each assignment are
  assignment.colors = getAssignmentColor(assignment = assignment)

  graphics::par(mar=c(1.1,1.1,1.1,1.1))
  graphics::plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  graphics::legend("top",inset = 0,
  legend = assignment,
  col = assignment.colors,
  lwd = 2, lty = 1, horiz = TRUE)

  #return p-values for different assignment mechanisms
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