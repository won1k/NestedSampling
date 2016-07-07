#' Comparisons of evidence computation by MC methods on Gamma-Poisson conjugate model
#' Nested Sampling, SIS, SMC
#' Mismatched prior [i.e. expectation too small]
#' 
#' \lambda \sim Gamma(\alpha, \beta) [\alpha = 3, \beta = 1]
#' y_1, \dots, y_n | \lambda \sim Pois(\lambda) (iid) [true \lambda = 10]
#' 

#'
#' Setup files
#' 
setwd("~/R/NS")
source("NS_det_pois.R")
source("NS_rand_pois.R")
source("SIS.R")
source("posterior_moments.R")


#'
#' Setup of variables
#' 
rm(list = ls())
set.seed(17)
#' We define a "target" list that contains functions relating to the prior and the likelihood.

# Gamma prior
# generate theta ~ Normal(mu, tau^2), in the form of a matrix
rprior <- function(nparticles, parameters){
  particles <- rgamma(nparticles, shape = parameters$prior_shape, rate = parameters$prior_rate)
  return(matrix(particles, ncol = 1))
}


# evaluate the log-density of the prior, for each particle
dprior <- function(thetaparticles, parameters){
  logdensities <- dgamma(thetaparticles[,1], shape = parameters$prior_shape, rate = parameters$prior_rate, log = TRUE)
  return(logdensities)
}


# Poisson log-likelihood for one observation, and a matrix of particles
# Y_i ~ Normal(theta, sigma^2) where sigma is denoted likelihood_sd
loglikelihood <- function(thetaparticles, observation, parameters){
  logdensities <- dpois(observation, lambda = thetaparticles[,1], log = TRUE)
  return(logdensities)
}

# target parameters
target_parameters <- list(prior_shape = 3, prior_rate = 1)
# group all the objects relating to the problem in one list
target <- list(dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, parameters = target_parameters)



#' 
#' Setup of simulation parameters and tuning parameters
#' 
N <- 100 # number of simulation runs
datalength <- 100
tuning_parameters <- list(f = 0.01, nparticles = 10^2, nestimates = 100)


#' Importance sampling vs. Nested sampling comparison:
#' 1. Generate 100 random points (y_1, \dots, y_n) \sim N(0,1) [same data every run]
#' 2. Compute exact log evidence based on conjugacy
#' 3. Compute IS log evidence
#' 4. Compute NS log evidence
#' 

# 1. Generate 100 data points
observations <- matrix(rpois(datalength, lambda = 10), ncol = 1)

IS.pois.log.evidences <- rep(0,N)
NS.det.pois.log.evidences <- rep(0,N)
NS.det.pois.iterations <- rep(0,N)
NS.rand.pois.log.evidences <- rep(0,N)
NS.rand.pois.var.evidences <- rep(0,N)
NS.rand.pois.iterations <- rep(0,N)

# 2. Compute posterior moments, extract exact log evidence
#posterior_moments <- posterior_moments_gam_pois(observations, target)
exact_log_evidence <- posterior_moments_gam_pois(observations, target)

for (i in 1:N) {
  print(i)
  
  # 3. IS log evidence
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.pois.log.evidences[i] <- IS.results$log_evidences[datalength]
  cat("IS step",i, "complete\n")
  
  # 4a. NS log evidence
  NS.results <- NS_det_pois(observations, tuning_parameters, target)
  NS.det.pois.log.evidences[i] <- log(NS.results$evidence)
  NS.det.pois.iterations[i] <- NS.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand_pois(observations, tuning_parameters, target)
  NS.rand.pois.log.evidences[i] <- NS.rand.results$log.evidence
  NS.rand.pois.var.evidences[i] <- NS.rand.results$var.log.evidence
  NS.rand.pois.iterations[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

# Save simulations
write(IS.pois.log.evidences, "IS_pois_log_evidences.txt")
write(NS.det.pois.log.evidences, "NS_det_pois_log_evidences.txt")
write(NS.rand.pois.log.evidences, "NS_rand_pois_log_evidences.txt")

# Plot simulations
plot(IS.pois.log.evidences, type = "l", col = 2, lty = 2, main = "IS vs NS Evidence Estimation, Gamma-Poisson, 10^2 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS.det.pois.log.evidences, lty = 5, col = 3)
lines(NS.rand.pois.log.evidences, lty = 4, col = 4)
abline(h = exact_log_evidence)
legend(x = 0, y = -129.9, c("IS","NS-det","NS-rand"), lty = c(2,5,4), col = c(2,3,4))
var(IS.pois.log.evidences) # 137.466
var(NS.det.pois.log.evidences) # 0.0558
var(NS.rand.pois.log.evidences) # 0.0574
exact_log_evidence # -249.668
mean(IS.pois.log.evidences) # -259.552
mean(NS.det.pois.log.evidences) # -250.643 
mean(NS.rand.pois.log.evidences) # -250.569

# RMSE
IS.rmse <- sqrt(mean((IS.pois.log.evidences - exact_log_evidence)^2)) # 15.290
NS.det.rmse <- sqrt(mean((NS.det.pois.log.evidences - exact_log_evidence)^2)) # 1.003
NS.rand.rmse <- sqrt(mean((NS.rand.pois.log.evidences - exact_log_evidence)^2)) # 0.931

# Iterations
mean(NS.det.pois.iterations) # 700.83
mean(NS.rand.pois.iterations) # 700.48





#' 
#' Mismatch prior (too high; actual lambda is low)
#' \lambda \sim Gamma(\alpha, \beta) [\alpha = 2, \beta = 1]
#' y_1, \dots, y_n | \lambda \sim Pois(\lambda) (iid) [true \lambda = 0.5]
#' 

# target parameters
target_parameters <- list(prior_shape = 2, prior_rate = 1)
# group all the objects relating to the problem in one list
target <- list(dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, parameters = target_parameters)



#' 
#' Setup of simulation parameters and tuning parameters
#' 
N <- 100 # number of simulation runs
datalength <- 100
tuning_parameters <- list(f = 0.01, nparticles = 10^2, nestimates = 100)


#' Importance sampling vs. Nested sampling comparison:
#' 1. Generate 100 random points (y_1, \dots, y_n) \sim N(0,1) [same data every run]
#' 2. Compute exact log evidence based on conjugacy
#' 3. Compute IS log evidence
#' 4. Compute NS log evidence
#' 

# 1. Generate 100 data points
observations.l <- matrix(rpois(datalength, lambda = 0.5), ncol = 1)

IS.pois.log.evidences.l <- rep(0,N)
NS.det.pois.log.evidences.l <- rep(0,N)
NS.det.pois.iterations.l <- rep(0,N)
NS.rand.pois.log.evidences.l <- rep(0,N)
NS.rand.pois.var.evidences.l <- rep(0,N)
NS.rand.pois.iterations.l <- rep(0,N)

# 2. Compute posterior moments, extract exact log evidence
#posterior_moments <- posterior_moments_gam_pois(observations, target)
exact_log_evidence <- posterior_moments_gam_pois(observations.l, target)

for (i in 1:N) {
  print(i)
  
  # 3. IS log evidence
  IS.results <- SIS(observations.l, tuning_parameters, target)
  IS.pois.log.evidences.l[i] <- IS.results$log_evidences[datalength]
  cat("IS step",i, "complete\n")
  
  # 4a. NS log evidence
  NS.results <- NS_det_pois(observations.l, tuning_parameters, target)
  NS.det.pois.log.evidences.l[i] <- log(NS.results$evidence)
  NS.det.pois.iterations.l[i] <- NS.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand_pois(observations.l, tuning_parameters, target)
  NS.rand.pois.log.evidences.l[i] <- NS.rand.results$log.evidence
  NS.rand.pois.var.evidences.l[i] <- NS.rand.results$var.log.evidence
  NS.rand.pois.iterations.l[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

plot(IS.pois.log.evidences.l, type = "l", col = 2, lty = 2, main = "IS vs NS Evidence Estimation, Gamma-Poisson, 10^2 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS.det.pois.log.evidences.l, lty = 5, col = 3)
lines(NS.rand.pois.log.evidences.l, lty = 4, col = 4)
abline(h = exact_log_evidence)
legend(x = 20, y = -100, c("IS","NS-det","NS-rand"), lty = c(2,5,4), col = c(2,3,4))
var(IS.pois.log.evidences.l) # 137.466
var(NS.det.pois.log.evidences.l) # 0.0558
var(NS.rand.pois.log.evidences.l) # 0.0574
exact_log_evidence # -249.668
mean(IS.pois.log.evidences.l) # -259.552
mean(NS.det.pois.log.evidences.l) # -250.643 
mean(NS.rand.pois.log.evidences.l) # -250.569

# RMSE
IS.rmse.l <- sqrt(mean((IS.pois.log.evidences.l - exact_log_evidence)^2)) # 0.414
NS.det.rmse.l <- sqrt(mean((NS.det.pois.log.evidences.l - exact_log_evidence)^2)) # 0.
NS.rand.rmse.l <- sqrt(mean((NS.rand.pois.log.evidences.l - exact_log_evidence)^2)) # 0.931

# Iterations
mean(NS.det.pois.iterations) # 700.83
mean(NS.rand.pois.iterations) # 700.48


