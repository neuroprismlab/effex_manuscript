##########################################################################
# Given the mean and standard deviation of a normal distribution, obtain the 
# proportion of that distribution that falls under a certain value in magnitude
#   > Essentially a wrapper for pnorm that is convenient for characterizing the 
#     cross-brain effect size distribution mean=mu and sd=tau)
#
# example:
# tau_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)
# proportions <- sapply(tau_list, function(tau) prop_abs_below(sd=tau))
#
#######################################################################


prop_abs_below <- function(mean=0, sd, threshold=0.2){
  
  if (threshold < 0) {
    stop('Threshold must be positive for two-tailed test - proportion of magnitudes less than this amount (i.e., between -thresh and +thresh).')
  }
  
  # proportion with absolute value below threshold
  return(pnorm(abs(threshold), mean=mean, sd=sd) - pnorm(-abs(threshold), mean=mean, sd=sd))
  
}
