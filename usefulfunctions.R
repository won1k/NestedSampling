library(Rcpp)
### 
# useful functions for the SMC sampler
# compute the effective sample size
ESSfunction <- function(normw){
  1. / (sum(normw^2))
}

# systematic resampling
# takes a vector of normalized weights, a number of desired draws, and a uniform random variable
# returns a vector of ancestors a_1, ..., a_N
# such that, marginally, for all k,j in {1,...,M},  Probability(a_k = j) = w_j
systematic_resampling <- function(normalized_weights, ndraws, u){
  N <- length(normalized_weights)
  indices <- rep(0, N)
  normalized_weights <- N * normalized_weights
  j <- 1
  csw <- normalized_weights[1]
  u <- runif(1, min = 0, max = 1)
  for (k in 1:N){
    while (csw < u){
      j <- j + 1
      csw <- csw + normalized_weights[j]
    }
    indices[k] <- j
    u <- u + 1
  }
  return(indices)
}

## functions to compute the mean and the covariance of a weighted sample
## wmean_ takes a matrix X (size N x d) and a vector of not-necessarily normalized weights (size N)
## and returns the weighted mean of the columns of X (of size d)
cppFunction('
NumericVector wmean(const NumericMatrix & x, const NumericVector & unnormalized_w){
  int nrows = x.rows();
  int ncols = x.cols();
  double sumw = sum(unnormalized_w);
  NumericVector result(ncols);
  double cumsumxw;
  for (int icol = 0; icol < ncols; icol++){
    cumsumxw = 0.;
    for (int irow = 0; irow < nrows ; irow++){
      cumsumxw += unnormalized_w(irow) * x(irow, icol);
    }
    result(icol) = cumsumxw / sumw;
  }
  return result;
}
')

## wcovariance_ takes a matrix X (size N x d) and a vector of not-necessarily normalized weights (size N),
## as well as the mean vector (presumably computed with the above function, wmean_)
## and returns the covariance of the columns of X (of size d x d)
cppFunction('
NumericMatrix wcovariance(const NumericMatrix & x, const NumericVector & unnormalized_w, const NumericVector & xbar){
  int nrows = x.rows();
  int ncols = x.cols();
  double sumw = sum(unnormalized_w);
  double sumsqw = sum(unnormalized_w * unnormalized_w);
  NumericMatrix result(ncols, ncols);
  std::fill(result.begin(), result.end(), 0);
  for (int i = 0; i < ncols; i++){
    for (int j = 0; j < ncols; j++){
      for (int irow = 0; irow < nrows ; irow++){
        result(i, j) += unnormalized_w(irow) * (x(irow, i) - xbar(i)) * (x(irow, j) - xbar(j));
      }
      result(i,j) /=  (sumw - sumsqw / sumw);
    }
  }
  return result ;
}
')


## functions to deal with multivariate normal distributions
# we can start by using the package mvtnorm which is slow but does the job
library(mvtnorm)
# we will use the functions dmvnorm and rmvnorm 
# the reason why mvtnorm is slow is that its functions check that your given covariance
# is semi-definite positive

# see faster functions in the comment below.
# NumericMatrix rmvnorm(int nsamples, const NumericVector & mean, const NumericMatrix & covariance){
#   RNGScope scope;
#   int ncols = covariance.cols();
#   const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
#   Eigen::MatrixXd cholesky_covariance(covariance_.llt().matrixU());
#   Eigen::MatrixXd Y(nsamples, ncols);
#   for(int i = 0; i < ncols; i++){
#     Y.col(i) = as<Eigen::ArrayXd>(rnorm(nsamples));
#   }
#   Y = Y * cholesky_covariance;
#   for(int j = 0; j < ncols; j++){
#     for(int i = 0; i < nsamples; i++){
#       Y(i,j) = Y(i,j) + mean(j);
#     }
#   }
#   return wrap(Y);
# }
# 
# 
# NumericMatrix centered_rmvnorm(int nsamples, const Eigen::MatrixXd & cholesky_covariance){
#   // sample centered gaussian variates given the upper triangular factor in the cholesky decomposition of the covariance
#   RNGScope scope;
#   int ncols = cholesky_covariance.cols();
#   Eigen::MatrixXd Y(nsamples, ncols);
#   for(int i = 0; i < ncols; i++){
#     Y.col(i) = as<Eigen::ArrayXd>(rnorm(nsamples));
#   }
#   Y = Y * cholesky_covariance;
#   return wrap(Y);
# }
# 
# // [[Rcpp::export]]
# NumericVector dmvnorm(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance){
#   const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
#   const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
#   Eigen::LLT<Eigen::MatrixXd> lltofcov(covariance_);
#   Eigen::MatrixXd lower = lltofcov.matrixL();
#   Eigen::MatrixXd xcentered(x_);
#   double halflogdeterminant = lower.diagonal().array().log().sum();;
#   double cst = - (halflogdeterminant) - (x.cols() * 0.9189385);
#   for(int j = 0; j < x.cols(); j++){
#     for(int i = 0; i < x.rows(); i++){
#       xcentered(i,j) = xcentered(i,j) - mean(j);
#     }
#   }
#   Eigen::VectorXd results = -0.5 * lower.triangularView<Eigen::Lower>().solve(xcentered.transpose()).colwise().squaredNorm();
#   for (int i = 0; i < results.size(); i++){
#     results(i) = results(i) + cst;
#   }
#   return wrap(results);
# }


