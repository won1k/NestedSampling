#' ---
#' title: "Normal Normal model: exact computation and Importance Sampling"
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
#+ presets, echo = FALSE
rm(list = ls())
set.seed(17)
#' We define a "target" list that contains functions relating to the prior and the likelihood.

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

#' The prior and the likelihood are conjugate so that we can compute exactly the moments of the posterior distributions.
#' Note also that we can compute the normalizing constant (also called the evidence) exactly. Indeed, the normalizing constant
#' $p(Y_{1}, Y_{2}, ... , Y_{n})$ satisfies
#' \[
#' p(Y_{1}, Y_{2}, ... , Y_{n}) = p(Y_1) \prod_{i=2}^n p(Y_i| Y_1, \ldots, Y_{i-1}),
#' \]
#' and each incremental evidence satisfies
#' \[
#' p(Y_i | Y_{1}, Y_{2}, ... , Y_{i-1}) = \int p(Y_i | \theta) p(\theta | Y_{1}, Y_{2}, ... , Y_{i-1}) d\theta.
#' \]
#' This integral has a closed-form solution in the Normal case, given by:
#' \[
#' \int\varphi_{\mathcal{N}}(x;\mu_{1},\lambda_{1})\varphi_{\mathcal{N}}(x;\mu_{2},\lambda_{2})dx=\frac{1}{\sqrt{2\pi\left(\lambda_{1}^{-1}+\lambda_{2}^{-1}\right)}}\exp\left(-\frac{1}{2\left(\lambda_{1}^{-1}+\lambda_{2}^{-1}\right)}\left(\mu_{1}-\mu_{2}\right)^{2}\right).
#' \]
#' where $\varphi_{\mathcal{N}}(x;\mu,\lambda)$ denotes the Normal density evaluated at $x$,
#' with mean $\mu$ and precision $\lambda$ (that is, one over the variance). 
#' Let's test this, to be sure.
mu1 <- 0.234
mu2 <- 0.4
lambda1 <- 0.74
lambda2 <- 2
f <- function(x) dnorm(x, mu1, 1/sqrt(lambda1)) * dnorm(x, mu2, 1/sqrt(lambda2))
integrate(f = f, lower = -10, upper = 10)
dnorm(mu1 - mu2, mean = 0, sd = sqrt(1/lambda1 + 1/lambda2))

#' Now let's compute the posterior moments and the log-incremental evidences.
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
#' Now let us generate some data, and compute the exact posterior moments, and the evidence.
#+ plot_data, warning = FALSE

# generate data
datalength <- 100
observations <- matrix(rnorm(datalength), ncol = 1)

# compute exact posterior moments
posterior_quantities <- posterior_moments(observations, target)

# we cumulatively sum the incremental log evidences to retrieve the log evidences
exact_log_evidences <- cumsum(posterior_quantities$log_incremental_evidences)
library(ggplot2)
qplot(x = 1:datalength, y = observations, geom = "line") + xlab("# data") + theme_bw()


#' Now let us approximate the posterior moments and the evidence using a Sequential Importance Sampling
#' scheme. This is most basic, and it is not expected to work well in general, but in this simple setting it should be
#' fine. Furthermore, it is in some sense the foundation on which sequential Monte Carlo methods are built.
#'

# number of samples
nparticles <- 10^5

# draw samples from prior
thetaparticles <- rprior(nparticles, target$parameters)

# vector of log-IS weights 
logweights <- rep(0, nparticles)

# vector of IS weights, normalized to sum to one 
normalized_weights <- rep(1/nparticles, nparticles)

# vector of log evidence, that is, at step i, log p(Y_{1}, Y_{2}, ... , Y_{i})
IS_log_evidences <- rep(0, datalength)

# go through each observation
for (idata in 1:datalength){
  # update the (log) importance weight with the new log-likelihood 
  logweights <- logweights + 
    target$loglikelihood(thetaparticles, observations[idata,], target$parameters)
  # remove maximum, for numerical stability
  maxlogweights <- max(logweights)
  # exponentiate the log-weights
  weights <- exp(logweights - maxlogweights)
  # compute the log evidence, using a Monte Carlo approximation of the formula:
  # p(Y_{1}, Y_{2}, ... , Y_{i}) 
  # = integral p(Y_{1}, Y_{2}, ... , Y_{i} | theta) p(theta) dtheta
  # (don't forget to add the maximum log-weight back)
  IS_log_evidences[idata] <- maxlogweights + log(mean(weights))
  normalized_weights <- weights / sum(weights)
}


#' SIS algorithm implementation:
#' observations = (y_1, \dots, y_n) = data; nparticles = number of thetaparticles in simulation; rprior = sample from prior function
#' target = (dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, parameters = target_parameters)
#' 
SIS <- function(observations, nparticles, rprior, target) {
  datalength <- length(observations)
  
  # draw samples from prior
  thetaparticles <- rprior(nparticles, target$parameters)
  
  # vector of log-IS weights 
  logweights <- rep(0, nparticles)
  
  # vector of IS weights, normalized to sum to one 
  normalized_weights <- rep(1/nparticles, nparticles)
  
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
  }
  
  return(list(log_evidences = IS_log_evidences, normalized_weights = normalized_weights, theta = thetaparticles))
}

IS_ex <- SIS(observations, nparticles, rprior, target)



#
#' We now plot the evidence calculated exactly, with the evidence approximated by Monte Carlo.
g <- qplot(x = 1:datalength, y = IS_log_evidences, geom = "blank") + geom_line(aes(colour = "importance sampling"))
g <- g + xlab("# data") + ylab("log evidence") + theme_bw()
g <- g + geom_line(aes(y = exact_log_evidences, colour = "exact"), linetype = 2)
g <- g + scale_color_discrete(name = "method:")
print(g)
#' The approximation seems very good (we used a lot of samples!).
#' 




#' As a sanity check, we also plot the sample approximation of the final posterior,
#' in a histogram. We overlay the density of the exact posterior distribution, as a red curve.
#' 
g <- qplot(x = thetaparticles[,1], weight = normalized_weights, geom = "blank") +
  geom_histogram(aes(y = ..density..), bins = 1000) + theme_bw()
g <- g + stat_function(fun = function(x) dnorm(x, mean = posterior_quantities$mean[datalength], 
                                        sd = 1/sqrt(posterior_quantities$precision[datalength])),
                n = 1000, colour = "red") 
g <- g + xlab(expression(theta))
print(g)


