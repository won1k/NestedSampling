#' ---
#' title: "Normal Normal model: SMC sampler"
#' author: "--"
#' date: "February 2016"
#' output:
#'  pdf_document:
#'    fig_width: 15
#'    fig_height: 5
#' ---
#' We consider a Normal prior: $\theta\sim\mathcal{N}(\mu,\tau^2)$, and $n$ observations $Y_1,\ldots,Y_n$
#' independently distributed as $Y_i\sim\mathcal{N}(\theta,\sigma^2)$.
#' 
#+ presets, echo = FALSE, warning = FALSE, message = FALSE
rm(list = ls())
set.seed(17)
library(doMC)
library(foreach)
library(doRNG)
library(dplyr)
library(reshape2)
registerDoMC(cores = 10)

#' We define a "target" list that contains functions relating to the prior and the likelihood.
#+
# Gaussian prior
# generate theta ~ Normal(mu, tau^2), in the form of a matrix
rprior <- function(nparticles, parameters){
  particles <- rnorm(nparticles, mean = parameters$prior_mean, 
                     sd = parameters$prior_sd)
  return(matrix(particles, ncol = 1))
}
# evaluate the log-density of the prior, for each particle
dprior <- function(thetaparticles, parameters){
  logdensities <- dnorm(thetaparticles[,1], mean = parameters$prior_mean, 
                        sd = parameters$prior_sd, log = TRUE)
  return(logdensities)
}
# Gaussian log-likelihood for one observation, and a matrix of particles
# Y_i ~ Normal(theta, sigma^2) where sigma is denoted likelihood_sd
loglikelihood <- function(thetaparticles, observation, parameters){
  logdensities <- dnorm(observation, mean = thetaparticles[,1], 
                        sd = parameters$likelihood_sd, log = TRUE)
  return(logdensities)
}
# target parameters
target_parameters <- list(prior_mean = 0, prior_sd = 1, likelihood_sd = 1)
# group all the objects relating to the problem in one list
target <- list(dimension = 1, rprior = rprior, dprior = dprior, 
               loglikelihood = loglikelihood, parameters = target_parameters)

#+ compute_exact, echo = FALSE
posterior_moments <- function(observations, target){
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
#' Now let us generate some data, and compute the exact posterior moments, and the evidence.
#+ plot_data, warning = FALSE
# generate data
datalength <- 1000
observations <- matrix(rnorm(datalength), ncol = 1)
# compute exact posterior moments
posterior_quantities <- posterior_moments(observations, target)
# we cumulatively sum the incremental log evidences to retrieve the log evidences
exact_log_evidences <- cumsum(posterior_quantities$log_incremental_evidences)
library(ggplot2)
qplot(x = 1:datalength, y = observations, geom = "line") + xlab("# data") + theme_bw()
#' Now let us approximate the posterior moments and the evidence using a Sequential Monte Carlo Sampler (Del Moral, Doucet, Jasra 2006).
#' We assimilate the data one by one, the algorithm is thus also called Iterated Batch Importance Sampling (Chopin 2002).
# import some basic functions for the ESS, 
# for the resampling step
# and for the computation of moments from weighted samples
source("usefulfunctions.R")
smcsampler <- function(tuning_parameters, target, observations){
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

# now run and compare with SIS
nparticles <- 1000
tuning_parameters <- list(nparticles = nparticles, nmoves = 5, ESSthreshold = 0.9)
source("SIS.R")
comparison.df <- data.frame()
nrep <- 100
comparison.df <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
  SISresults <- SIS(observations, nparticles, target)
  SMCresults <- smcsampler(tuning_parameters, target, observations)
  data.frame(irep = irep, idata = 1:datalength, 
             SIS_log_evidences = SISresults$log_evidences, 
             SMC_log_evidences = cumsum(SMCresults$log_evidence_incr))
}


# now compare with exact evidence and compute root mean squared error
exact.df <- data.frame(idata = 1:datalength, exact_log_evidences = exact_log_evidences)
comparison.df <- merge(comparison.df, exact.df, by = "idata")
comparison.df <- comparison.df %>% 
  mutate(SIS_error = (SIS_log_evidences - exact_log_evidences)^2,
         SMC_error = (SMC_log_evidences - exact_log_evidences)^2)
#' We plot the error of each of 100 runs of SIS and SMC.
comparison.melt.df <- melt(comparison.df %>% 
               select(irep, idata, SIS_error, SMC_error), c("irep", "idata"))
gerror <- ggplot(comparison.melt.df, 
                 aes(x = idata, y = value, group = interaction(variable, irep), 
                     colour = variable))
gerror <- gerror + geom_line()  + theme_bw()
gerror <- gerror + xlab("# data") + ylab(expression(paste(L[2], " error"))) 
gerror <- gerror + scale_colour_discrete(name = "method:") + theme(legend.position = "bottom")
print(gerror)
#' We plot the root mean squared error of SIS and SMC, computed from the 100 runs.
rmse.df <- comparison.df %>% group_by(idata) %>% 
  summarize(SIS_rmse = sqrt(mean(SIS_error)), SMC_rmse = sqrt(mean(SMC_error)))

rmse.df <- melt(rmse.df, "idata")
grmse <- ggplot(rmse.df, aes(x = idata, y = value, group = variable, colour = variable)) +
  geom_line() + theme_bw()
grmse <- grmse + xlab("# data") + ylab("RMSE") +
  scale_colour_discrete(name = "method:") + theme(legend.position = "bottom")
print(grmse)