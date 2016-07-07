#' Computing posterior moments for Normal-Normal model
#' observations = (y_1, \dots, y_n) = data
#' target = (dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, parameters = target_parameters)
#' 
#' 
Stirling <- function(n) ((2*pi)^.5)*(n/exp(1))^n
logStirling <- function(n) n*log(n)-n+0.5*log(2*pi*n)

posterior_moments <- function(observations, target) {
  datalength <- nrow(observations)
  prior_precision <- 1/(target$parameters$prior_sd^2)
  posterior_means <- rep(0, datalength)
  posterior_precisions <- rep(0, datalength)
  log_incremental_evidences <- rep(0, datalength)
  
  # the posterior precisions and means are going to be recursively updated
  posterior_precision <- 1/(target$parameters$prior_sd^2)
  posterior_mean <- target$parameters$prior_mean
  
  # go through each observation 
  for (idata in 1:datalength){
    y <- observations[idata,1]
    
    # for the incremental evidence we use the above formula
    log_incremental_evidences[idata] <- dnorm(posterior_mean - y, 
                                              mean = 0, sd = sqrt(1 / posterior_precision + target$parameters$likelihood_sd^2),
                                              log = TRUE)
    
    # for the posterior moments, we use the standard conjugacy formula
    posterior_mean <- (1/(1/(target$parameters$likelihood_sd^2) + posterior_precision)) *
      (posterior_precision*posterior_mean + y/(target$parameters$likelihood_sd^2))
    posterior_precision <- 1/(target$parameters$likelihood_sd^2) + posterior_precision
    posterior_means[idata] <- posterior_mean
    posterior_precisions[idata] <- posterior_precision
  }
  return(list(mean = posterior_means, precision = posterior_precisions, 
              log_incremental_evidences = log_incremental_evidences))
}


#' Computing posterior moments for Gamma-Poisson model
#' observations = (y_1, \dots, y_n) = data
#' target = (dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, parameters = target_parameters)
#' 
posterior_moments_gam_pois <- function(observations, target) {
  datalength <- nrow(observations)
  #log_incremental_evidences <- rep(0, datalength)
  prior_shape <- target$parameters$prior_shape
  prior_rate <- target$parameters$prior_rate
  
  #for (idata in 1:datalength) {
  #  y <- observations[idata,1]
    
    # independent log increments
    #log_incremental_evidences[idata] <- prior_shape * log(prior_rate) - log(gamma(prior_shape) * factorial(y)) + log(gamma(prior_shape+y)) - (prior_shape+y) * log(prior_rate)
  #}
  
  log_evidence <- prior_shape*log(prior_rate) + log(gamma(prior_shape + sum(observations))) - log(gamma(prior_shape)) - sum(log(factorial(observations))) - (prior_shape+sum(observations))*log(prior_rate+datalength)
  #log_evidence <- prior_shape * log(prior_rate) + logStirling(prior_shape + sum(observations) - 1) - log(gamma(prior_shape)) - sum(logStirling(observations)) - (prior_shape + sum(observations)) * log(prior_rate + datalength)
  return(log_evidence)
  #return(list(log_evidence = sum(log_incremental_evidences), log_incremental_evidences = log_incremental_evidences))
}



#' Computing exact evidence for simple 2-Gaussian mixture model
#' observations = (y_1, \dots, y_n)
#' target = ...
#'
exact_evidence_2g_mix <- function(observations, target) {
  datalength <- nrow(observations)
  prior_prob <- target$parameters$prior_prob
  means <- target$parameters$means
  sds <- target$parameters$sds
  
  log_evidence <- log(prior_prob) + log(exp(sum(dnorm(observations, mean = means[1], sd = sds[1], log = TRUE))) + exp(sum(dnorm(observations, mean = means[2], sd = sds[2], log = TRUE))))
}
