#' Implementation of Iterated Batch Importance Sampling (IBIS)
#' Uses MCMC steps with MVN proposals when ESS < c * N
#' 
#' tuning_parameters = (nparticles, nmoves, ESSthreshold = 0.9)
#' target = (dimension = 1, rprior, dprior, loglikelihood, parameters = (prior_mean, prior_sd, likelihood_sd))
#' observations = (y_1, \dots, y_n)
#'

library(mvtnorm)
library(foreach)
library(doRNG)
library(dplyr)
source("usefulfunctions.R")

IBIS <- function(tuning_parameters, target, observations){
  nparticles <- tuning_parameters$nparticles
  ESSthreshold <- tuning_parameters$ESSthreshold
  nmoves <- tuning_parameters$nmoves
  datalength <- nrow(observations)
  # initialize particle system
  currentparticles <- target$rprior(nparticles, parameters = target$parameters)
  logtargetvalues <- target$dprior(currentparticles, target$parameters)
  logweights <- rep(0, nparticles)
  normweights <- rep(1/nparticles, nparticles)
  ESSarray <- rep(0, datalength)
  resamplingtimes <- c()
  log_evidence_incr <- rep(0, datalength)
  for (idata in 1:datalength){
    ## importance sampling step
    # incremental log-likelihood {log u_t^i} for i = 1, ..., N
    loglikelihood_incr <- 
      target$loglikelihood(currentparticles, observations[idata,], 
                           target$parameters)
    # compute the incremental evidence
    # that is
    # sum_{i=1}^N w_{t-1}^i * u_t^i
    # = sum_{i=1}^N w_{t-1}^i * exp( log u_t^i )
    # = exp(C) sum_{i=1}^N w_{t-1}^i * exp( log u_t^i - C)
    # where C = max(log u_t^i)
    max_log_incr <- max(loglikelihood_incr)
    log_evidence_incr[idata] <- 
      max_log_incr +
      log(sum(normweights * exp(loglikelihood_incr - max_log_incr)))
    # update the weights to w_t^i for i=1,...,N
    logweights <- logweights + loglikelihood_incr
    # update the target density values at this point (useful for the MCMC steps later)
    logtargetvalues <- logtargetvalues + loglikelihood_incr
    # normalize the weights
    max_log_weights <- max(logweights)
    w <- exp(logweights - max_log_weights)
    normweights <- w / sum(w)
    # compute ESS
    currentESS <- ESSfunction(normweights)
    ESSarray[idata] <- currentESS
    # if ESS is smaller than the threshold, perform resample move steps
    if (currentESS < ESSthreshold * nparticles){
      cat("At step", idata, ", ESS = ", currentESS / nparticles * 100, "% \n")
      cat("-> resample-move\n")
      ## compute mean and covariance of current particles
      mean_particles <- wmean(currentparticles, normweights)
      cov_particles <- wcovariance(currentparticles, normweights, mean_particles)
      ## resampling step
      resamplingtimes <- c(resamplingtimes, idata)
      ancestors <- systematic_resampling(normweights, nparticles, runif(1))
      currentparticles <- matrix(currentparticles[ancestors,], 
                                 ncol = target$dimension)
      logtargetvalues <- logtargetvalues[ancestors]
      logweights <- rep(0, nparticles)
      ## move step (possibly multiple ones, according to 'nmoves')
      if (nmoves > 0){
        for (imove in 1:nmoves){
          # propose parameter values according to 
          # an independent Normal distribution fitted on the current particles
          proposals <- rmvnorm(n = nparticles, mean = mean_particles, 
                               sigma = cov_particles)
          proposals_targetd <- target$dprior(proposals, target$parameters)
          for (jdata in 1:idata){
            proposals_targetd <- proposals_targetd + 
              target$loglikelihood(proposals, observations[jdata,], 
                                   target$parameters)
          }
          # compute numerator and denominator of Metropolis-Hastings ratio
          numerator <- proposals_targetd + 
            dmvnorm(x = currentparticles, mean = mean_particles, sigma = cov_particles, 
                    log = TRUE)
          denominator <- logtargetvalues + 
            dmvnorm(x = proposals, mean = mean_particles, sigma = cov_particles, 
                    log = TRUE)
          # acceptance ratio, on the log-scale
          acceptance_ratio <- numerator - denominator
          uniforms <- runif(nparticles)
          accepts <- (log(uniforms) < acceptance_ratio)
          cat("Mean acceptance rate = ", mean(accepts) * 100, "% \n")
          currentparticles[accepts,] <- matrix(proposals[accepts,], 
                                               ncol = target$dimension)
          logtargetvalues[accepts] <- proposals_targetd[accepts]
        }
      }
    }
    # step <- step + 1
  }
  return(list(currentparticles = currentparticles, normweights = normweights, 
              ESSarray = ESSarray,
              resamplingtimes = resamplingtimes, log_evidence_incr = log_evidence_incr))
}