library(distr)

#' Stochastic implementation (\Delta x = e^{-i/nparticles}) of nested sampling for Normal-Normal example
#' observations = (y_1, \dots, y_n)
#' tuning_parameters = (nparticles, nestimates, f = 0.01, nmoves, ESSthreshold = 0.9)
#' target = (dimension = 1, rprior, dprior, loglikelihood, parameters = (prior_mean, prior_sd, likelihood_sd))
#' 
#' 
NS_rand_2g <- function(observations, tuning_parameters, target) {
  # Setup parameters
  datalength <- length(observations)
  lf <- log(tuning_parameters$f)
  # nparticles = number of theta particles
  nparticles <- tuning_parameters$nparticles
  # nestimates = number of parallel streams x_t = x_{t-1} * b_t for b_t \sim Beta(nparticles, 1)
  nestimates <- tuning_parameters$nestimates
  
  # Sample alive theta particles
  theta <- target$rprior(nparticles, target$parameters)
  t <- 1
  
  # Log-likelihood vector for alive thetas
  ll <- rep(0, nparticles)
  for(j in 1:nparticles) ll[j] <- sum(target$loglikelihood(theta[j], observations, target$parameters))
  # select minimum log-likelihood
  lphi <- min(ll)
  i <- which(ll == lphi)[1]
  theta_samp <- theta[i]
  
  # Draw quadrature \sim Beta(nparticles, 1) [r_{t,k} are factors x_{t,k} are quadrature points]
  r <- rbeta(nestimates, nparticles, 1)
  x <- matrix(r, ncol = nestimates)
  # Evaluate evidence (Z_k)
  Z <- exp(lphi) * (1 - x)
  
  while( max(ll[-i]) - t/nparticles > lf + mean(log(Z)) ) {
    t <- t + 1
    
    # Sample new alive theta; MC until within desired region
    theta[i] <- rprior(1, target$parameters)
    while( sum(target$loglikelihood(theta[i], observations, target$parameters)) < lphi[t-1] ) {
      theta[i] <- rprior(1, target$parameters)
    }
    
    # Compute likelihood and select minimum
    ll[i] <- sum(target$loglikelihood(theta[i], observations, target$parameters))
    lphi <- c(lphi, min(ll))
    i <- which(ll == lphi[t])[1]
    theta_samp <- c(theta_samp, theta[i])
    
    # Draw quadrature points (Note: for speed, only keep current and immediate past x_t)
    r <- rbeta(nestimates, nparticles, 1)
    x <- matrix(c(x[nrow(x),], x[nrow(x),] * r), ncol = nestimates, byrow = T)
    # x <- matrix(c(as.vector(t(x)), x[t-1,] * r), ncol = nestimates, byrow = T)
    # Incrementally update evidence Z
    Z <- Z + exp(lphi[t]) * (x[1,] - x[2,])
  }
  
  # Add remainder to Z
  mean_remainder <- mean(exp(ll[-i]))
  Z <- Z + mean_remainder * exp(-(t+1)/(N-1))
  
  # Average log(Z_k) values
  avg.log.Z <- mean(log(Z))
  var.log.Z <- var(as.vector(log(Z)))
  
  # Compile output
  output <- list()
  output$log.evidence <- avg.log.Z
  output$all.log.evidence <- log(Z)
  output$log.incremental.evidence <- lphi
  output$var.log.evidence <- var.log.Z
  output$iterations <- t
  #output$quadrature <- x
  output$theta <- theta_samp
  
  return(output)
}