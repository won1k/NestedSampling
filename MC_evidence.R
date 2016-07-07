#'
#' MC computation of exact evidence for mixture Gaussians
#'
#'

#' Known means, known variances
#' \lambda ~ Bern(p)
#' y|\lambda = i ~ N(\mu_i, \sigma_i^2) = N(target$parameters$means[i], target$parameters$sds[i])
#' 
#' p(y) = \sum_{\lambda} p(y|\lambda)p(\lambda)
#' 
#' Sample: \lambda ~ Bern(p) 
#' Compute: Z_j = p(y|\lambda) p(\lambda) = p(y|\lambda) * 2^{-100}
#' Estimate: Z \approx (1/N) * sum(Z_j)

MC_evidence_known_means <- function(observations, parameters) {
  # Set up parameters
  nsamples <- parameters$nsamples
  prior_prob <- parameters$prior_prob
  means <- parameters$means
  sds <- parameters$sds
  n <- length(observations)
  
  # Array of inner integrals
  inner.sums <- rep(0, nsamples)
  
  # Compute inner integrals for MC simulated lambda
  for (i in 1:nsamples) {
    lambda <- rbinom(n, 1, prior_prob) + 1
    n.1 <- length(lambda[lambda == 1])
    n.2 <- length(lambda[lambda == 2])
    rows.1 <- which(lambda == 1)
    rows.2 <- which(lambda == 2)
    mean.obs.1 <- mean(observations[rows.1,1])
    mean.obs.2 <- mean(observations[rows.2,1])
    var.obs.1 <- var(observations[rows.1,1])
    var.obs.2 <- var(observations[rows.2,1])
    
    # Compute MC sum
    inner.sums[i] <- (2*pi*sds[1]^2)^(-n.1/2) * (2*pi*sds[2]^2)^(-n.2/2) * exp(-1/(2*sds[1]^2)*((n.1-1)*var.obs.1 + 
                      n.1*(mean.obs.1-means[1])^2) - 1/(2*sds[2]^2)*((n.2-1)*var.obs.2 + n.2*(mean.obs.2-means[2])^2))
  }
  
  log.inner.sums <- log(inner.sums)
  
  # Remove extreme outliers
  #qnt <- quantile(log.inner.integrals, probs = c(.25, .75), na.rm = TRUE)
  #H <- 1.5 * IQR(x, na.rm = TRUE)
  #log.inner.integrals <- log.inner.integrals[log.inner.integrals > (qnt[1] - H) & log.inner.integrals < (qnt[2] + H)]
  
  return(list(log_inner_sums = log.inner.sums, log_evidence = log(mean(inner.sums))))
}





# Unknown means, known variances
# \lambda ~ Bern(p)
# \mu_i ~ N(\theta, \tau^2) = N(hyper_mean, hyper_sd^2)
# y|\lambda, \mu_1, \mu_2 ~ N(\mu_i, \sigma_i^2) = N(\mu_i, obs_sd^2)
#
# p(y) = \int p(y|\lambda,\mu_1,\mu_2) p(\mu_1,\mu_2) p(\lambda)
#
# Sample \lambda ~ Bern(p) [N times]
# Compute: Z_j = p(y|\lambda,\mu_1,\mu_2) p(\mu_1) p(\mu_2)
# MC estimate: Z \approx (1/N) * sum(Z_j)

MC_evidence <- function(observations, parameters) {
  # Set up parameters
  nsamples <- parameters$nsamples
  prior_prob <- parameters$prior_prob
  hyper_mean <- parameters$hyper_mean
  hyper_sd <- parameters$hyper_sd
  obs_sd <- parameters$obs_sd
  n <- length(observations)
  
  # Array of inner integrals
  inner.integrals <- rep(0, nsamples)
  
  # Compute inner integrals for MC simulated lambda
  for (i in 1:nsamples) {
    lambda <- rbinom(n, 1, prior_prob) + 1
    n.1 <- length(lambda[lambda == 1])
    n.2 <- length(lambda[lambda == 2])
    rows.1 <- which(lambda == 1)
    rows.2 <- which(lambda == 2)
    mean.obs.1 <- mean(observations[rows.1,1])
    mean.obs.2 <- mean(observations[rows.2,1])
    mu.t1 <- ((n.1*mean.obs.1)/(2*obs_sd^2) + hyper_mean/(2*hyper_sd^2))/(n.1/(2*obs_sd^2) + 1/(2*hyper_sd^2))
    mu.t2 <- ((n.2*mean.obs.2)/(2*obs_sd^2) + hyper_mean/(2*hyper_sd^2))/(n.2/(2*obs_sd^2) + 1/(2*hyper_sd^2))
    
    inner.integrals[i] <- (2*pi*obs_sd^2)^(-n/2) * (2*pi*hyper_sd^2)^(-1) * exp(-1/(2*obs_sd^2)*(sum(observations[rows.1]^2)-n.1*mean.obs.1^2 + sum(observations[rows.2]^2) - n.1*mean.obs.2^2) - (n.1/(2*obs_sd^2)*mean.obs.1^2 + hyper_mean^2/(2*hyper_sd^2) + n.2/(2*obs_sd^2)*mean.obs.2^2 + hyper_mean^2/(2*hyper_sd^2)) + (n.1/(2*obs_sd^2) + 1/(2*hyper_sd^2)) * mu.t1^2 + (n.2/(2*obs_sd^2) + 1/(2*hyper_sd^2)) * mu.t2^2) * (2*pi/(n.1/obs_sd^2 + 1/hyper_sd^2))^(-1/2) * (2*pi/(n.2/obs_sd^2 + 1/hyper_sd^2))^(-1/2)
  }
  
  log.inner.integrals <- log(inner.integrals)
  
  # Remove extreme outliers
  #qnt <- quantile(log.inner.integrals, probs = c(.25, .75), na.rm = TRUE)
  #H <- 1.5 * IQR(x, na.rm = TRUE)
  #log.inner.integrals <- log.inner.integrals[log.inner.integrals > (qnt[1] - H) & log.inner.integrals < (qnt[2] + H)]
  
  return(list(log_inner_integrals = log.inner.integrals, log_evidence = mean(log.inner.integrals)))
}