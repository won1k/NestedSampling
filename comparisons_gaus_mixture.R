#' Comparisons of evidence computation by MC methods on Gaussian mixture model
#' Nested Sampling, SIS, SMC
#' 2 Gaussians (multimodal)
#' 
#' \lambda ~ Bern(p)
#' y | \lambda = i ~ N(\mu_i, \sigma^2_i)
#' 
#' Note: Start off with \mu_i, \sigma_i^2 known; then proceed to unknown (general) parameters)
#' p = 1/2 (prior); actual p = 0.3

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

# Bernoulli prior on \lambda
# generate lambda ~ Bern(1/2), in the form of a matrix
rprior <- function(nparticles, parameters){
  particles <- rbinom(nparticles, 1, parameters$prior_prob) + 1
  return(matrix(particles, ncol = 1))
}


# evaluate the log-density for each particle
dprior <- function(thetaparticles, parameters){
  logdensities <- dbinom(thetaparticles[,1], 1, parameters$prior_prob, log = TRUE)
  return(logdensities)
}


# Normal log-likelihood for one observation, and a matrix of particles
# Y_i ~ Normal(theta, sigma^2) where sigma is denoted likelihood_sd
loglikelihood <- function(thetaparticles, observation, parameters){
  nparticles <- length(thetaparticles[,1])
  logdensities <- rep(0, nparticles)
  for (i in 1:nparticles) logdensities[i] <- dnorm(observation, mean = parameters$means[thetaparticles[i,1]], sd = parameters$sds[thetaparticles[i,1]], log = TRUE)
  return(logdensities)
}

# Compute log-likelihood for all observations, and one particle
# Y_i ~ Normal(mu_i, sigma^2) where sigma is denoted obs_sd
loglikelihood_all <- function(thetaparticle, observations, parameters) {
  logdensities <- dnorm(observations, mean = parameters$means[thetaparticle], sd = parameters$sds[thetaparticle], log = TRUE)
  return(logdensities)
}

# target parameters
target_parameters <- list(prior_prob = 1/2, means = c(0, 2), sds = c(1,1), nsamples = 100000)
# group all the objects relating to the problem in one list
target <- list(dimension = 1, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, loglikelihood_all = loglikelihood_all, parameters = target_parameters)



#' 
#' Setup of simulation parameters and tuning parameters
#' 
N <- 100 # number of simulation runs
datalength <- 100
tuning_parameters <- list(f = 0.01, nparticles = 10^2)


#' 
#' Generate observations + compute evidence
#' 
lambda <- rbinom(datalength, 1, 0.3) # actual p = 0.3
observations <- matrix(rnorm(datalength, mean = target$parameters$means[lambda+1], sd = target$parameters$sds[lambda+1]), ncol = 1)

# Exact evidence
exact_MC_evidence_known_means <- MC_evidence_known_means(observations, target$parameters)$log_evidence

# IS, NS-det, NS-rand computation
IS.2g.log.evidences <- rep(0, N)
NS.det.2g.log.evidences <- rep(0, N)
NS.det.2g.iterations <- rep(0, N)
NS.rand.2g.log.evidences <- rep(0, N)
NS.rand.2g.iterations <- rep(0, N)

for (i in 1:N) {
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.2g.log.evidences[i] <- IS.results$log_evidences[datalength]
  cat("IS step",i, "complete\n")
  
  # 4a. NS log evidence
  NS.results <- NS_det(observations, tuning_parameters, target)
  NS.det.2g.log.evidences[i] <- log(NS.results$evidence)
  NS.det.2g.iterations[i] <- NS.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand(observations, tuning_parameters, target)
  NS.rand.2g.log.evidences[i] <- NS.rand.results$log.evidence
  NS.rand.2g.iterations[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

# Plot simulations
plot(IS.2g.log.evidences, type = "l", col = 2, lty = 2, ylim = c(-195.1, -194.1), main = "IS vs NS Evidence Estimation, 2-Gaussian Mixture, 10^2 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS.det.2g.log.evidences, lty = 5, col = 3)
lines(NS.rand.2g.log.evidences, lty = 4, col = 4)
abline(h = exact_MC_evidence_known_means)
legend(x = 0, y = -194.9, c("IS","NS-det","NS-rand"), lty = c(2,5,4), col = c(2,3,4))
var(IS.2g.log.evidences) # 0.0139
var(NS.det.2g.log.evidences) # 0.0254
var(NS.rand.2g.log.evidences) # 0.0189
exact_MC_evidence_known_means # -194.67
mean(IS.2g.log.evidences) # -194.35
mean(NS.det.2g.log.evidences) # -194.66 
mean(NS.rand.2g.log.evidences) # -194.63

# RMSE
sqrt(mean((IS.2g.log.evidences - exact_MC_evidence_known_means)^2)) # 0.344
sqrt(mean((NS.det.2g.log.evidences - exact_MC_evidence_known_means)^2)) # 0.159
sqrt(mean((NS.rand.2g.log.evidences - exact_MC_evidence_known_means)^2)) # 0.142




#'
#' 2-Normal Mixture
#' Unknown means, known variances
#' 
#' \lambda ~ Bern(p)
#' \mu_i ~ N(\theta, \tau^2)
#' y | \lambda = i, \mu_i ~ N(\mu_i, \sigma^2_i)
#' 
#' \mu_i unknown, \sigma_i^2 known
#' p = 1/2 (prior); actual p = 0.3
#' 

# Independent prior on \lambda, \mu_1, \mu_2
# generate in the form of a matrix
rprior <- function(nparticles, parameters) {
  particles <- matrix(nrow = nparticles, ncol = 3)
  particles[,1] <- rbinom(nparticles, 1, parameters$prior_prob)
  particles[,2] <- rnorm(nparticles, mean = parameters$hyper_mean, sd = parameters$hyper_sd)
  particles[,3] <- rnorm(nparticles, mean = parameters$hyper_mean, sd = parameters$hyper_sd)
  return(particles)
}


# evaluate the log-density for each particle
dprior <- function(thetaparticles, parameters) {
  logdensities <- dbinom(thetaparticles[,1], 1, parameters$prior_prob, log = TRUE) +
    dnorm(thetaparticles[,2], mean = parameters$hyper_mean, sd = parameters$hyper_sd, log = TRUE) +
    dnorm(thetaparticles[,3], mean = parameters$hyper_mean, sd = parameters$hyper_sd, log = TRUE) 
  return(logdensities)
}


# Compute log-likelihood for one observation, and a matrix of particles
# Y_i ~ Normal(mu_i, sigma^2) where sigma is denoted obs_sd
loglikelihood <- function(thetaparticles, observation, parameters) {
  nparticles <- length(thetaparticles[,1])
  logdensities <- rep(0, nparticles)
  for (i in 1:nparticles) logdensities[i] <- dnorm(observation, mean = thetaparticles[i, thetaparticles[i,1] + 2], sd = parameters$obs_sd, log = TRUE)
  return(logdensities)
}

# Compute log-likelihood for all observations, and one particle
# Y_i ~ Normal(mu_i, sigma^2) where sigma is denoted obs_sd
loglikelihood_all <- function(thetaparticle, observations, parameters) {
  logdensities <- dnorm(observations, mean = thetaparticle[thetaparticle[1] + 2], sd = parameters$obs_sd, log = TRUE)
  return(logdensities)
}

# target parameters
target_parameters <- list(prior_prob = 1/2, hyper_mean = 0, hyper_sd = 10, obs_sd = 1, nsamples = 100)
# group all the objects relating to the problem in one list
target <- list(dimension = 3, rprior = rprior, dprior = dprior, loglikelihood = loglikelihood, loglikelihood_all = loglikelihood_all, parameters = target_parameters)

#' 
#' Setup of simulation parameters and tuning parameters
#' 
N <- 100 # number of simulation runs
datalength <- 100
tuning_parameters <- list(f = 0.01, nparticles = 10^2, nestimates = 100)

#' 
#' Generate observations + compute evidence
#' 
lambda <- rbinom(datalength, 1, target$parameters$prior_prob) + 1
means <- rnorm(2, mean = 2, sd = 5) # 3.55, -0.4
observations <- matrix(rnorm(datalength, mean = means[lambda], sd = 1), ncol = 1)

# Exact evidence
exact_MC_evidence <- MC_evidence(observations, target$parameters)

# SIS/NS evidence computations
IS.2g.means.log.evidences <- rep(0, N)
NS.det.2g.means.log.evidences <- rep(0, N)
NS.det.2g.means.iterations <- rep(0, N)
NS.rand.2g.means.log.evidences <- rep(0, N)
NS.rand.2g.means.iterations <- rep(0, N)
IS.2g.means.ESS <- matrix(rep(0, N*datalength), ncol = datalength)

for (i in 1:N) {
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.2g.means.log.evidences[i] <- IS.results$log_evidences[datalength]
  IS.2g.means.ESS[i,] <- IS.results$ESS
  cat("IS step",i, "complete\n")
  
  # 4a. NS log evidence
  NS.results <- NS_det(observations, tuning_parameters, target)
  NS.det.2g.means.log.evidences[i] <- log(NS.results$evidence)
  NS.det.2g.means.iterations[i] <- NS.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand(observations, tuning_parameters, target)
  NS.rand.2g.means.log.evidences[i] <- NS.rand.results$log.evidence
  NS.rand.2g.means.iterations[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

# Save simulations
write(IS.2g.means.log.evidences, "IS_2g_means_log_evidences.txt")
write(NS.det.2g.means.log.evidences, "NS_det_2g_means_log_evidences.txt")
write(NS.rand.2g.means.log.evidences, "NS_rand_2g_means_log_evidences.txt")

# Plot log evidences
plot(IS.2g.means.log.evidences, type = "l", lty = 5, col = 2, ylim = c(-365, -335), main = "IS vs NS Evidence Estimation, Two-Normal Mixture, 10^2 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS.det.2g.means.log.evidences, lty = 2, col = 3)
lines(NS.rand.2g.means.log.evidences, lty = 3, col = 4)
abline(h = exact_MC_evidence$log_evidence)
legend(x = 75, y = -358, c("IS","NS-det","NS-rand"), lty = c(2,5,4), col = c(2,3,4))

# Plot ESS
plot(IS.2g.means.ESS[1,], type = "l", ylim = c(0,20), xlab = "", ylab = "Effective Sample Size", main = "ESS for SIS Sampler, Two-Normal, 10^2 Particles")
for (i in 2:N) {
  lines(IS.2g.means.ESS[i,])
}

# Metrics
exact_MC_evidence$log_evidence # -339.06
mean(IS.2g.means.log.evidences) # -340.210
mean(NS.det.2g.means.log.evidences) # -338.708
mean(NS.rand.2g.means.log.evidences) # -338.695
var(IS.2g.means.log.evidences) # 16.754
var(NS.det.2g.means.log.evidences) # 0.0439
var(NS.rand.2g.means.log.evidences) # 0.0394

# RMSE
sqrt(mean((IS.2g.means.log.evidences - exact_MC_evidence$log_evidence)^2)) # 4.476
sqrt(mean((NS.det.2g.means.log.evidences - exact_MC_evidence$log_evidence)^2)) # 0.411
sqrt(mean((NS.rand.2g.means.log.evidences - exact_MC_evidence$log_evidence)^2)) # 0.395


#'
#' 2-Normal Mixture
#' Unknown means, known variances
#' 
#' nparticles = 10^3 Particles
#' 

#' 
#' Setup of simulation parameters and tuning parameters
#' 
N <- 50 # number of simulation runs; 50 for computational costs
datalength <- 100
tuning_parameters <- list(f = 0.01, nparticles = 10^3, nestimates = 100)

# SIS/NS evidence computations
IS.2g.means.log.evidences.10e3 <- rep(0, N)
NS.det.2g.means.log.evidences.10e3 <- rep(0, N)
NS.det.2g.means.iterations.10e3 <- rep(0, N)
NS.rand.2g.means.log.evidences.10e3 <- rep(0, N)
NS.rand.2g.means.iterations.10e3 <- rep(0, N)

for (i in 1:N) {
  IS.results <- SIS(observations, tuning_parameters, target)
  IS.2g.means.log.evidences.10e3[i] <- IS.results$log_evidences[datalength]
  cat("IS step",i, "complete\n")
  
  # 4a. NS log evidence
  NS.results <- NS_det(observations, tuning_parameters, target)
  NS.det.2g.means.log.evidences.10e3[i] <- log(NS.results$evidence)
  NS.det.2g.means.iterations.10e3[i] <- NS.results$iterations
  cat("NS-det step",i, "complete\n")
  
  # 4b. NS-rand log evidence
  NS.rand.results <- NS_rand(observations, tuning_parameters, target)
  NS.rand.2g.means.log.evidences.10e3[i] <- NS.rand.results$log.evidence
  NS.rand.2g.means.iterations.10e3[i] <- NS.rand.results$iterations
  cat("NS-rand step",i, "complete\n")
}

# Save simulations
write(IS.2g.means.log.evidences.10e3, "IS_2g_means_log_evidences_10e3.txt")
write(NS.det.2g.means.log.evidences.10e3, "NS_det_2g_means_log_evidences_10e3.txt")
write(NS.rand.2g.means.log.evidences.10e3, "NS_rand_2g_means_log_evidences_10e3.txt")

# Plot log evidences
plot(IS.2g.means.log.evidences.10e3, type = "l", lty = 5, col = 2, main = "IS vs NS Evidence Estimation, Two-Normal Mixture, 10^3 Particles, f = 0.01", xlab = "", ylab = "Log-Evidence")
lines(NS.det.2g.means.log.evidences.10e3, lty = 2, col = 3)
lines(NS.rand.2g.means.log.evidences.10e3, lty = 3, col = 4)
abline(h = exact_MC_evidence$log_evidence)
legend(x = 0, y = -339.1, c("IS","NS-det","NS-rand"), lty = c(2,5,4), col = c(2,3,4))

# Metrics
exact_MC_evidence$log_evidence # -338.715
mean(NS.det.2g.means.log.evidences.10e3) # -338.711
mean(NS.rand.2g.means.log.evidences.10e3) # -338.703
mean(IS.2g.means.log.evidences.10e3) # -338.725
var(NS.det.2g.means.log.evidences.10e3) # 0.00447
var(NS.rand.2g.means.log.evidences.10e3) # 0.00428
var(IS.2g.means.log.evidences.10e3) # 0.0770

# RMSE
sqrt(mean((IS.2g.means.log.evidences.10e3 - exact_MC_evidence$log_evidence)^2)) # 0.275
sqrt(mean((NS.det.2g.means.log.evidences.10e3 - exact_MC_evidence$log_evidence)^2)) # 0.0663
sqrt(mean((NS.rand.2g.means.log.evidences.10e3 - exact_MC_evidence$log_evidence)^2)) # 0.0659
