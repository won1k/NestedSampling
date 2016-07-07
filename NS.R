#' Nested Sampling Comparisons
#' Won I. Lee
#' 

library(distr)

#' Deterministic implementation (\Delta x = e^{-i/N}) of nested sampling
#' pi0 is prior; l is log-likelihood function; N is number of alive particles; f is precision of estimate
#' 
NS_det <- function(pi0, l, N, f = 0.001) {
  n <- 100
  lf <- log(f)
  lphi <- rep(NA, n)
  theta_samp <- rep(NA, n)
  
  prior <- AbscontDistribution(d = pi0)
  rprior <- r(prior)
  theta <- rprior(N)
  
  t = 1
  ll <- rep(NA, N)
  for(j in 1:N) ll[j] <- l(theta[j])
  lphi[t] <- min(ll)
  i <- which(ll == lphi[t])[1]
  #print(i)
  theta_samp[t] <- theta[i]
  
  Z <- exp(lphi[t]) * (1-exp(-t/N))
  
  while( max(ll[-i]) - t/N > lf + log(Z)) {
    t <- t+1
    if(t %% n == 0) {
      m = length(phi)
      lphi <- lphi[1:(m+n)]
      theta_samp <- theta_samp[1:(m+n)]
    }
    
    theta[i] <- rprior(1)
    while(l(theta[i]) < lphi[t-1]) {
      theta[i] <- rprior(1)
    }
    
    ll[i] <- l(theta[i])
    lphi[t] <- min(ll)
    i <- which(ll == lphi[t])[1]
    theta_samp[t] <- theta[i]
    
    Z <- Z + exp(lphi[t]) * (exp(-(t-1)/N) - exp(-t/N))
  }
  
  return(list(evidence = Z, log_incremental_evidence = lphi, theta = theta_samp))
}

#NS_rand


#' Examples & Comparisons
#' (1) Normal
#' 

s_0 = 10
s_1 = 1
x_ex <- rnorm(1, mean = 0, sd = sqrt(s_0^2+s_1^2))

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

# Using N = 5
NS_norm_ex5 <- rep(NA, 100)
for (i in 1:100) {
  print(i)
  NS_norm_ex5[i] <- NS_det(nprior, nll, 5)$evidence
}
plot(NS_norm_ex5, xlab = "", ylab = "Z Estimates", main = "Nested Sampling, N = 5", type = "l", lty = 2)
abline(h = dnorm(x_ex, mean = 0, sd = sqrt(s_0^2+s_1^2)))

# Using N = 100
NS_norm_ex100 <- rep(NA, 100)
for (i in 1:100) {
  print(i)
  NS_norm_ex100[i] <- NS_det(nprior, nlike, 100, 0.05)
}
plot(NS_norm_ex100, xlab = "", ylab = "Z Estimates", main = "Nested Sampling, N = 100", ylim = c(0.005, 0.015))
abline(h = dnorm(x_ex, mean = 0, sd = sqrt(s_0^2+s_1^2)))

# Using N = 1000
NS_norm_ex100 <- rep(NA, 100)
for (i in 1:100) {
  print(i)
  NS_norm_ex100[i] <- NS_det(nprior, nlike, 1000, 0.05)
}
plot(NS_norm_ex100, xlab = "", ylab = "Z Estimates", main = "Nested Sampling, N = 1000", ylim = c(0.005, 0.015))
abline(h = dnorm(x_ex, mean = 0, sd = sqrt(s_0^2+s_1^2)))

