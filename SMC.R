#' Generic SMC implementation with 
#' 
#' 
#' 

library(distr)
library(stats)

# Gaussian log-likelihood for one observation, and a matrix of particles
# Y_i ~ Normal(theta, sigma^2) where sigma is denoted likelihood_sd
loglikelihood <- function(thetaparticles, observation, parameters){
  logdensities <- dnorm(observation, mean = thetaparticles[,1], 
                        sd = parameters$likelihood_sd, log = TRUE)
  return(logdensities)
}

#' IBIS implementation; use N(\hat{\mu}_k, \hat{\Sigma}_k)
#' pi0 = initial proposal; y = (y_1, \dots, y_n) = data; posteriors = (\pi_1, \dots, \pi_n) = sequential posteriors;
#'  N = # theta particles; c = fraction of ESS at which to resample

IBIS = function(pi0, y, N, c) {
  s_0 = 10
  s_1 = 1
  n <- length(y)
  thetaparticles <- array(rep(0, (n+1)*N), dim = c(n+1,N))
  prior <- AbscontDistribution(d = pi0)
  rprior <- r(prior)
  thetaparticles[1,] <- rprior(N)
  ess <- N
  W <- array(rep(0, (n+1)*N), dim = c(n+1,N))
  W[1,] <- 1/N
  
  for (k in 1:n) {
    if (ess < c*N) {
      # Note that k is actually the previous step since k+1 is current step
      offspring <- rmultinom(1, N, W[k,])
      current <- 0
      # Resample \tilde(\theta)_k
      theta <- rep(0, N)
      for (i in 1:N) {
        theta[(current+1):(current+offspring[i])] <- thetaparticles[k,i]
      }
      # Move \tilde(\theta)_k to \theta_{k+1} using empirical mean/var Normal proposal
      empmean <- mean(theta); empvar <- var(theta)
      proposal <- rnorm(N, mean = empmean, sd = sqrt(empvar))
      propmean <- mean(proposal); propvar <- var(proposal)
      # Compute pi_k(\theta) (posterior) ~ N(...) (i.e. incremental posterior up to y_k)
      posterior <- function(theta) dnorm(theta, mean = (k*mean(y[1:k])/s_1^2) / (1/s_0^2 + k/s_1^2), sd = 1/(1/s_0^2 + k/s_1^2))
      alpha <- (posteriors[k](proposal) * dnorm(thetaparticles[k,], mean = propmean, sd = sqrt(propvar))) / 
        (posteriors[k](thetaparticles[k,]) * dnorm(proposal, mean = empmean, sd = sqrt(empvar)))
      for (i in 1:N) {
        if (runif(1) < alpha[i]) {
          thetaparticles[k+1,i] <- proposal[i]
        } else {
          thetaparticles[k+1,i] <- thetaparticles[k,i]
        }
      }
    } else {
      thetaparticles[k+1,] <- thetaparticles[k,]
    }
    w <- W[k,] * dnorm(y[k+1], mean = thetaparticles[k,])
    W[k+1,] <- w / sum(w)
  }
  
  return(list(weights = W, theta = thetaparticles))
}

#' IBIS on Normal-Normal example
#' 

s_0 = 10
s_1 = 1
datalength <- 100
observations <- rnorm(datalength)

nll <- function(theta, x = x_ex, sd = s_1) {
  sum(dnorm(x, mean = theta, sd = sd, log = T))
}

nprior <- function(theta, mean = 0, sd = s_0) {
  dnorm(theta, mean = mean, sd = sd)
}

results <- IBIS(nprior, observations, 100, 0.5)


#' SMC implementation
#' pi0 = initial proposal; c = fraction of ESS at which to resample; 

SMC = function(pi0, c, ) {
  
}
