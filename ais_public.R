rm(list = ls())
set.seed(13)
# creating data
A <- matrix(0, nrow = 4, ncol = 36)
feature_count <- dim(A)[1]
dim_count <- dim(A)[2]
#
A[0+1, 1 + c(1 , 6 , 7 , 8 , 13)] <- 1
A[1+1, 1 + c(3 , 4 , 5 , 9 , 11 , 15 , 16 , 17)] <- 1
A[2+1, 1 + c(18 , 24 , 25 , 30 , 31 , 32) ] <- 1
A[3+1, 1 + c(21 , 22 , 23 , 28 , 34) ] <- 1
# setting up some parameters for later
Ndata <- 5
sigma_n <- 0.1
sigma_a <- 5
alpha <- 2
# data-generating Z
Zstar <- matrix(rbinom(n = Ndata * feature_count, size = 1, prob = 0.25), nrow = Ndata, ncol = feature_count)
Zstar
# data itself
X <- Zstar %*% A + sigma_n * matrix(rnorm(Ndata*feature_count, mean = 0, sd = 1), nrow = Ndata, ncol = dim_count)
# store the transpose of X too
tX <- t(X)
# pre-compute some boring matrices
diagNdata <- diag(1, Ndata, Ndata)
diagfeature_count <- diag(1, feature_count, feature_count)
# and some other stuff...
log2pi <- log(2*pi)
# collapsed log-likelihood
cllikelihood <- function(Z, sigma_a, sigma_n){
  tZ <- t(Z)
  inner_term <- tZ %*% Z + (sigma_n^2) / (sigma_a^2) * diagfeature_count
  tr_term <- tX %*% (diagNdata - Z %*% solve(inner_term) %*% tZ) %*% X
  ll <- (- Ndata * dim_count / 2) * log2pi +
    ( -1. * ( Ndata - feature_count ) * dim_count ) * log( sigma_n ) +
    ( -1. * ( feature_count * dim_count ) ) * log( sigma_a ) +
    ( -1. * ( dim_count / 2. ) ) * log(det(inner_term)) +
    ( -1. / 2.0 / ( sigma_n**2 ) ) * sum(diag(tr_term))
  return(ll)
}
# test it
cllikelihood(Zstar, sigma_a, sigma_n)
cllikelihood(Zstar, sigma_a, 0.01)
cllikelihood(matrix(rbinom(Ndata*feature_count, 1, prob = 0.25), ncol = feature_count, nrow = Ndata), sigma_a, sigma_n)

# generate Z from finite feature allocation model (i.e. the prior on Z)
generate_Z <- function(alpha, N, K){
  pi <- rbeta(n = K, shape1 = alpha/K, shape2 = 1)
  Z <- matrix(runif(N*K), ncol = K, nrow = N)
  Z <- t(apply(X = Z, MARGIN = 1, FUN = function(v) as.numeric(v < pi)))
  return(Z)
}

# one sweep of Gibbs updates, where the target is
# (prior(Z) times likelihood(Z)^beta
# for beta in [0,1];
# takes the current log-likelihood as an argument, so that only one new likelihood evaluation
# is required at each step (instead of two)
resample_Z_weighted <- function(Z, current_cll, alpha, sigma_a, sigma_n, beta){
  for (k in 1:feature_count){
    for (n in 1:Ndata){
      cll1 <- 0
      cll0 <- 0
      if (Z[n,k] == 1){
        cll1 <- current_cll
        Z[n,k] <- 0
        cll0 <- cllikelihood(Z, sigma_a, sigma_n)
      } else {
        cll0 <- current_cll
        Z[n,k] <- 1
        cll1 <- cllikelihood(Z, sigma_a, sigma_n)
      }
      mk <- sum(Z[,k]) - Z[n,k]
      p_z1 <- (mk + alpha/feature_count)/(Ndata + alpha/feature_count)
      logp_z0 <- beta * cll0 + log(1-p_z1)
      logp_z1 <- beta * cll1 + log(p_z1)
      m <- max(c(logp_z0, logp_z1))
      logprob1 <- logp_z1 - m - log(exp(logp_z0 - m) + exp(logp_z1 - m))
      z <- (log(runif(1)) < logprob1)
      Z[n,k] <- z
      if (z == 1){
        current_cll <- cll1
      } else {
        current_cll <- cll0
      }
    }
  }
  return(list(Z = Z, current_cll = current_cll))
}

## Annealed importance sampling
AIS <- function(Nparticles, beta_count, beta_step, alpha, sigma_a, sigma_n){
  # store final Z's
  Z_set <- list()
  # store log(R hat)
  logRhat <- rep(0, Nparticles)
  # for each particle, independently...
  for (iparticle in 1:Nparticles){
    # initialize from prior
    Z <- generate_Z(alpha, Ndata, feature_count)
    #
    lweight <- 0
    beta <- 0
    ll_set <- rep(0, beta_count)
    cll <- cllikelihood(Z, sigma_a, sigma_n)
    for (beta_index in 1:beta_count){
      beta <- beta + beta_step
      cll <- cllikelihood(Z, sigma_a, sigma_n)
      lweight <- lweight + beta_step * cll
      resample_res <- resample_Z_weighted(Z, cll, alpha, sigma_a, sigma_n, beta)
      cll <- resample_res$current_cll
      Z <- resample_res$Z
      ll_set[beta_index] <- cll
    }
    logRhat[iparticle] <- lweight
    Z_set[[iparticle]] <- Z
  }
  return(list(logRhat = logRhat, Z_set = Z_set))
}

# number of particles (or chains)
Nparticles <- 50
# number of temperature's beta
beta_count = 100
# space between two beta's
beta_step = 1.0 / beta_count
# NOTE: it's highly sub-optimal to use a linear grid of temperatures
# the temperatures can instead be chosen adaptively based on some criterion
# such as the effective sample size
AIS_results <- AIS(Nparticles, beta_count, beta_step, alpha, sigma_a, sigma_n)
AIS_results$logRhat
# compute log(mean(Rhat))
maxlogevid <- max(AIS_results$logRhat)
maxlogevid + log(mean(exp(AIS_results$logRhat - maxlogevid)))

