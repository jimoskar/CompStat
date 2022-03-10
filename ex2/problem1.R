### Problem 1 ----
library(ggplot2)
load("rain.rda")

## a) ----
head(rain)
ggplot(rain, aes(x = day, y = n.rain)) + 
  geom_jitter() + geom_smooth() + 
  theme_minimal()

## e) ----
T <- 366 # days in a year
n <- rep(39, T) # n in binom. distr.
n[60] = 10 # corrigate for feb 29th
y <- rain$n.rain # response

# Sigmoid function
sigm <- function(tau){
  return(exp(tau)/(1 + exp(tau)))
}
log.sigm <- function(tau){
  return(-log(exp(-tau) + 1))
}
log.1m.sigm <- function(tau){
  return(-log(exp(tau) + 1))
}

# Funciton to sample from inverse gamma
rigamma <- function(n, a, b){
  # uses shape and scale
  return(1/rgamma(n, shape = a, rate = 1/b))
}

# Function to calculate right part of accept. prob.
calc.alpha <- function(t, tau.prop, tau.old){
  pi.prop <- sigm(tau.prop)
  pi.old <- sigm(tau.old)
  log.frac <- y[t]*(log.sigm(tau.prop) - log.sigm(tau.old)) + 
                      (n[t] - y[t])*(log.1m.sigm(tau.prop) - log.1m.sigm(tau.old))
  if(t == 60){
    log.frac = 0
  }
  return(exp(log.frac))
}

# Function for returning qt
calc.q.t <- function(t){
  q <- rep(0, T - 1)
  if(t == 1){
    q[1] = -1
  } else if(t == T){
    q[T - 1] <= -1
  } else{
    q[t - 1] = q[t] = -1
  }
  return(q)
}


# Construct Q without the 1/sigma factor
Q <- diag(rep(2, T))
Q[row(Q) - col(Q) == 1] <-  Q[row(Q) - col(Q) == -1] <- -1
Q[1,1] <- Q[T,T] <- 1

# MCMC with iterative conditioning
mcmc.iterative <- function(iter, sigma0, tau0){
  tau.mat <- matrix(NA, nrow = iter, ncol = T)
  tau.mat[1, ] <- tau0
  sigma.vec <- rep(NA, iter)
  sigma.vec[1] <- sigma0
  for(i in 2:iter){
    # Sample tau
    for(j in 1:T){
      # Generate proposal
      q.tt <- 1/sigma.vec[i-1] * Q[j,j]
      q.t <- calc.q.t(j)
      
      Q.AA <- 1/sigma.vec[i - 1] * q.tt
      Q.AB <- 1/sigma.vec[i - 1] * t(q.t)
      mu.AgB <- 1/Q.AA * Q.AB %*% (tau.mat[i -1, 1:T != j])
      Q.AgB <- q.tt
      tau.prop <- rnorm(1, as.numeric(mu.AgB), solve(Q.AgB))
      
      # Calculate acceptance prob.
      alpha <- min(1, calc.alpha(j, tau.prop, tau.mat[i - 1, j]))
      if(is.na(alpha)){
        browser()
      }
      u <- runif(1)
      if(u < alpha){
        tau.mat[i, j] = tau.prop
      } else{
        tau.mat[i, j] = tau.mat[i-1, j]
      }
    }
    # Generate Gibbs sample for sigma
    shape <- alpha + 0.5
    scale <- 0.5*t(tau.mat[i, ]) %*% Q %*% tau.mat[i, ] + beta
    sigma.vec[i] <- rigamma(1, shape, scale)
  }
  return(list(tau.mat = tau.mat, sigma.vec = sigma.vec))
}

# Parameters for sigma prior
alpha <- 10
beta <- 20

# Run MCMC
iter <- 10000
set.seed(4300)
mcmc <- mcmc.iterative(iter, 10, runif(T))
plot(1:iter, mcmc$sigma.vec, type = "l", ylim = c(0,0.1))
