#' SIS algorithm implementation:
#' observations = (y_1, \dots, y_n) = data; nparticles = number of thetaparticles in simulation; rprior = sample from prior function
#' target = (dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, parameters = target_parameters)
#' 
SIS <- function(observations, tuning_parameters, target) {
  datalength <- length(observations)
  nparticles <- tuning_parameters$nparticles
  
  # draw samples from prior
  thetaparticles <- target$rprior(nparticles, target$parameters)
  
  # vector of log-IS weights 
  logweights <- rep(0, nparticles)
  
  # vector of IS weights, normalized to sum to one 
  normalized_weights <- rep(1/nparticles, nparticles)
  
  # list of ESS at each step
  ESS <- rep(0, datalength)
  
  # vector of log evidence, that is, at step i, log p(Y_{1}, Y_{2}, ... , Y_{i})
  IS_log_evidences <- rep(0, datalength)
  
  # go through each observation
  for (idata in 1:datalength){
    # update the (log) importance weight with the new log-likelihood 
    logweights <- logweights + target$loglikelihood(thetaparticles, observations[idata,], target$parameters)
    
    # remove maximum, for numerical stability
    maxlogweights <- max(logweights)
    
    # exponentiate the log-weights
    weights <- exp(logweights - maxlogweights)
    
    # compute the log evidence, using a Monte Carlo approximation
    IS_log_evidences[idata] <- maxlogweights + log(mean(weights))
    normalized_weights <- weights / sum(weights)
    
    # compute the ESS
    ESS[idata] <- 1/sum(normalized_weights^2)
  }
  
  return(list(log_evidences = IS_log_evidences, normalized_weights = normalized_weights, ESS = ESS, theta = thetaparticles))
}
