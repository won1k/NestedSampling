library(distr)

#' Deterministic implementation (\Delta x = e^{-i/nparticles}) of nested sampling for nparticlesormal-nparticlesormal example
#' rprior = function to sample from prior; nparticles is number of alive particles; f is precision of estimate
#' target = list of parameters and functions for problem
#' 
NS_det_2g <- function(observations, tuning_parameters, target) {
  # Setup parameters
  datalength <- length(observations)
  n <- 100
  lf <- log(tuning_parameters$f)
  nparticles <- tuning_parameters$nparticles
  # Minimum log-likelihood of alive particles
  lphi <- rep(0, n)
  theta_samp <- rep(0, n)
  
  # Sample alive theta particles
  theta <- target$rprior(nparticles, target$parameters)
  
  t = 1
  # Log-likelihood vector for alive thetas
  ll <- rep(0, nparticles)
  for(j in 1:nparticles) ll[j] <- sum(target$loglikelihood(theta[j], observations, target$parameters))
  # select minimum log-likelihood
  lphi[t] <- min(ll)
  i <- which(ll == lphi[t])[1]
  theta_samp[t] <- theta[i]
  
  # Evaluate evidence
  Z <- exp(lphi[t]) * (1-exp(-t/nparticles))
  
  while( max(ll[-i]) - t/nparticles > lf + log(Z)) {
    t <- t+1
    # Lengthen arrays as necessary
    if(t %% n == 0) {
      m = length(lphi)
      lphi <- lphi[1:(m+n)]
      theta_samp <- theta_samp[1:(m+n)]
    }
    
    # Sample new alive theta; MC until within desired region
    theta[i] <- rprior(1, target$parameters)
    while(sum(target$loglikelihood(theta[i], observations, target$parameters)) < lphi[t-1]) {
      theta[i] <- rprior(1, target$parameters)
    }
    
    # Compute likelihood and select minimum
    ll[i] <- sum(target$loglikelihood(theta[i], observations, target$parameters))
    lphi[t] <- min(ll)
    i <- which(ll == lphi[t])[1]
    theta_samp[t] <- theta[i]
    
    # Incrementally update evidence Z
    Z <- Z + exp(lphi[t]) * (exp(-(t-1)/nparticles) - exp(-t/nparticles))
  }
  
  # Add remainder to Z
  mean_remainder <- mean(exp(ll[-i]))
  Z <- Z + mean_remainder * exp(-(t+1)/N)
  
  # Compile output
  output <- list()
  output$evidence <- Z
  output$iterations <- t
  output$log_incremental_evidence <- lphi
  output$theta <- theta_samp
  
  return(output)
}