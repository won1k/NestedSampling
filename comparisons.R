#' Comparisons of evidence computation by MC methods on various examples
#' Nested Sampling, SIS, SMC
#' 
#' 

#'
#' Setup files
#' 
setwd("~/R/NS")
source("NS_det.R")
source("NS_rand.R")
source("SIS.R")
source("posterior_moments.R")


#'
#' Setup of variables
#' 
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
target <- list(dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, parameters = target_parameters)



#' 
#' Setup of simulation parameters and tuning parameters
#' 
N <- 100 # number of simulation runs
datalength <- 100
f <- 0.01
tuning_parameters <- list(f = f, nparticles = 10^3, nestimates = 100)


#' Importance sampling vs. Nested sampling comparison:
#' 1. Generate 100 random points (y_1, \dots, y_n) \sim N(0,1) [same data every run]
#' 2. Compute exact log evidence based on conjugacy
#' 3. Compute IS log evidence
#' 4. Compute NS log evidence
#' 

# 1. Generate 100 data points
observations <- matrix(rnorm(datalength), ncol = 1)
tuning_parameters$nparticles <- 10^2

IS.log.evidences <- rep(0,N)
NS.det.log.evidences <- rep(0,N)
NS.det.iterations <- rep(0,N)
NS.rand.log.evidences <- rep(0,N)
NS.rand.var.evidences <- rep(0,N)
NS.rand.iterations <- rep(0,N)

# 2. Compute posterior moments, extract exact log evidence
posterior_quantities <- posterior_moments(observations, target)
exact_log_evidences <- cumsum(posterior_quantities$log_incremental_evidences)[datalength]

for (i in 1:N) {
  print(i)

  # 3. IS log evidence
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.log.evidences[i] <- IS.results$log_evidences[datalength]
  cat("IS step",i, "complete\n")
  
  # 4a. NS log evidence
  NS.results <- NS_det(observations, tuning_parameters, target)
  NS.det.log.evidences[i] <- log(NS.results$evidence)
  NS.det.iterations[i] <- NS.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand(observations, tuning_parameters, target)
  NS.rand.log.evidences[i] <- NS.rand.results$log.evidence
  NS.rand.var.evidences[i] <- NS.rand.results$var.log.evidence
  NS.rand.iterations[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

plot(IS_log_evidences, type = "l", col = 2, lty = 2, main = "IS vs NS Evidence Estimation, Normal-Normal, 10^2 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS_log_evidences, lty = 5, col = 3)
lines(NS.rand.log.evidences, lty = 4, col = 4)
abline(h = exact_log_evidences)
legend(x = 0, y = -130.95, c("IS","NS-det","NS-rand"), lty = c(2,5,4), col = c(2,3,4))
var(IS_log_evidences) # 0.0609
var(NS_log_evidences) # 0.0213
var(NS.rand.log.evidences) # 0.0160
exact_log_evidences # -130.480
mean(IS_log_evidences) # -130.524
mean(NS_log_evidences) # -130.463 
mean(NS.rand.log.evidences) # -130.483

# RMSE
IS.rmse <- sqrt(mean((IS_log_evidences - exact_log_evidences)^2)) # 0.248
NS.det.rmse <- sqrt(mean((NS_log_evidences - exact_log_evidences)^2)) # 0.146
NS.rand.rmse <- sqrt(mean((NS.rand.log.evidences - exact_log_evidences)^2)) # 0.127

# Iterations
mean(NS.det.iterations) # 690.8
mean(NS.rand.iterations) # 692.6


#' 
#' Try for nparticles = 10^3
#' 
tuning_parameters$nparticles <- 10^3
tuning_parameters$f <- 0.01

IS.log.evidences.10e3 <- rep(0,N)
NS.det.log.evidences.10e3 <- rep(0,N)
NS.det.niter.10e3 <- rep(0,N)
NS.rand.log.evidences.10e3 <- rep(0,N)
NS.rand.var.evidences.10e3 <- rep(0,N)
NS.rand.niter.10e3 <- rep(0,N)

# 2. Compute posterior moments, extract exact log evidence
posterior.quantities <- posterior_moments(observations, target)
exact.log.evidence.10e3 <- cumsum(posterior.quantities$log_incremental_evidences)[datalength]

for (i in 1:N) {
  print(i)
  
  # 3. IS log evidence
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.log.evidences.10e3[i] <- IS.results$log_evidences[datalength]
  cat("IS step",i, "complete\n")
  
  # 4. NS-det log evidence
  NS.det.results <- NS_det(observations, tuning_parameters, target)
  NS.det.log.evidences.10e3[i] <- log(NS.det.results$evidence)
  NS.det.niter.10e3[i] <- NS.det.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand(observations, tuning_parameters, target)
  NS.rand.log.evidences.10e3[i] <- NS.rand.results$log.evidence
  print(NS.rand.results$log.evidence)
  NS.rand.var.evidences.10e3[i] <- NS.rand.results$var.log.evidence
  NS.rand.niter.10e3[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

plot(IS.log.evidences.10e3, type = "l", lty = 2, col = 2, main = "IS vs NS Evidence Estimation, Normal-Normal, 10^3 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS.det.log.evidences.10e3, lty = 5, col = 3)
lines(NS.rand.log.evidences.10e3, lty = 4, col = 4)
abline(h = exact.log.evidence.10e3)
legend(x = 78.5, y = -130.23, c("IS","NS-det","NS-rand"), lty = c(3,2,3), col = c(2,3,4))
var(IS.log.evidences.10e3) # 0.00576
var(NS.det.log.evidences.10e3) # 0.00203
var(NS.rand.log.evidences.10e3) # 0.00179
exact.log.evidence.10e3 # -130.480
mean(IS.log.evidences.10e3) # -130.474
mean(NS.det.log.evidences.10e3) # -130.488
mean(NS.rand.log.evidences.10e3) # -130.488

# RMSE
IS.rmse.10e3 <- sqrt(mean((IS.log.evidences.10e3 - exact.log.evidence.10e3)^2)) # 0.0802
NS.det.rmse.10e3 <- sqrt(mean((NS.det.log.evidences.10e3 - exact.log.evidence.10e3)^2)) # 0.0401
NS.rand.rmse.10e3 <- sqrt(mean((NS.rand.log.evidences.10e3 - exact.log.evidence.10e3)^2)) # 0.00937

# Variations on evidences
mean(NS.rand.var.evidences) # 0.0189
mean(NS.rand.var.evidences.10e3) # 0.00195

# Number of iterations
mean(NS.det.niter.10e3) # 6921.68
mean(NS.rand.niter.10e3) # 6921.59

#' 
#' Try for nparticles = 10^3
#' With f = 0.05 instead; does this result in lower variance (i.e. bias-variance tradeoff)?
#' 
tuning_parameters$nparticles <- 10^3
tuning_parameters$f <- 0.05

IS.log.evidences.1000.005 <- rep(0,N)
NS.log.evidences.1000.005 <- rep(0,N)

# 2. Compute posterior moments, extract exact log evidence
posterior.quantities <- posterior.moments(observations, target)
exact.log.evidence.1000.005 <- cumsum(posterior.quantities$log.incremental.evidences)[datalength]

for (i in 1:N) {
  print(i)
  
  # 3. IS log evidence
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.log.evidences.1000.005[i] <- IS.results$log.evidences[datalength]
  
  # 4. NS log evidence
  NS.results <- NS_det(observations, tuning_parameters, target)
  NS.log.evidences.1000.005[i] <- log(NS.results$evidence)
}

plot(IS.log.evidences.1000.005, type = "l", lty = 3, main = "IS vs NS Evidence Estimation, Normal-Normal, 10^3 Particles, f = 0.05", xlab = "", ylab = "Log-Evidence")
lines(NS.log.evidences.1000.005, lty = 2)
abline(h = exact.log.evidences)
legend(x = 80, y = -158.2, c("IS","NS"), lty = c(3,2))
var(IS.log.evidences.1000.005) # 0.006493
var(NS.log.evidences.1000.005) # 0.002205
exact.log.evidence.1000.005 # -157.974
mean(IS.log.evidences.1000.005) # -157.988
mean(NS.log.evidences.1000.005) # -158.029



#' 
#' Try for nparticles = 10^4, f = 0.01
#' 
tuning_parameters$nparticles <- 10^4
tuning_parameters$f <- 0.01

IS.log.evidences.10e4 <- rep(0,N)
NS.det.log.evidences.10e4 <- rep(0,N)
NS.det.niter.10e4 <- rep(0,N)
NS.rand.log.evidences.10e4 <- rep(0,N)
NS.rand.var.evidences.10e4 <- rep(0,N)
NS.rand.niter.10e4 <- rep(0,N)

# 2. Compute posterior moments, extract exact log evidence
posterior.quantities <- posterior_moments(observations, target)
exact.log.evidence.10e4 <- cumsum(posterior.quantities$log_incremental_evidences)[datalength]

for (i in 1:N) {
  print(i)
  
  # 3. IS log evidence
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.log.evidences.10e4[i] <- IS.results$log_evidences[datalength]
  cat("IS step",i, "complete\n")
  
  # 4. NS-det log evidence
  NS.det.results <- NS_det(observations, tuning_parameters, target)
  NS.det.log.evidences.10e4[i] <- log(NS.det.results$evidence)
  NS.det.niter.10e4[i] <- NS.det.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand(observations, tuning_parameters, target)
  NS.rand.log.evidences.10e4[i] <- NS.rand.results$log.evidence
  print(NS.rand.results$log.evidence)
  NS.rand.var.evidences.10e4[i] <- NS.rand.results$var.log.evidence
  NS.rand.niter.10e4[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

plot(IS.log.evidences.10e4, type = "l", col = 2, lty = 2, main = "IS vs NS Evidence Estimation, Normal-Normal, 10^4 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS.det.log.evidences.10e4, col = 3, lty = 5)
lines(NS.rand.log.evidences.10e4, col = 4, lty = 4)
abline(h = exact.log.evidence.10e4)
legend(x = 0, y = -130.512, c("IS","NS-det","NS-rand"), lty = c(2,5,4), col = c(2,3,4))
var(IS.log.evidences.10e4) # 0.0004888
var(NS.det.log.evidences.10e4) # 0.0001399
var(NS.rand.log.evidences.10e4) # 0.0001795
exact.log.evidence.10e4 # -130.480
mean(IS.log.evidences.10e4) # -130.4785
mean(NS.det.log.evidences.10e4) # -130.4916
mean(NS.rand.log.evidences.10e4) # -130.4893

# RMSE
IS.rmse.10e4 <- sqrt(mean((IS.log.evidences.10e4 - exact.log.evidence.10e4)^2)) # 0.02204
NS.det.rmse.10e4 <- sqrt(mean((NS.det.log.evidences.10e4 - exact.log.evidence.10e4)^2)) # 0.01663
NS.rand.rmse.10e4 <- sqrt(mean((NS.rand.log.evidences.10e4 - exact.log.evidence.10e4)^2)) # 0.01633

# Variations on evidences
mean(NS.rand.var.evidences.10e4) # 0.0001883

# Number of iterations
mean(NS.det.niter.10e4) # 69247.87
mean(NS.rand.niter.10e4) # 69224.68


