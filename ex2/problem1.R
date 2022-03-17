### Problem 1 ----
library(ggplot2)
library(MASS)
load("rain.rda")

# sample autocorrelation
source('sacf.R')

## a) ----
ggplot(rain, aes(x = day, y = n.rain)) + 
  geom_point() + geom_smooth() + xlab("Day in year") + ylab("Number of days rain exceeds 1mm") +
  theme_minimal()

## e) ----
T <- 366 # days in a year
n <- rep(39, T) # n in binom. distr.
n[60] = 10 # corrigate for feb 29th
y <- rain$n.rain # response

# logit function
logit <- function(x){
  return(x/(1 - x))
}

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

# Function to sample from inverse gamma
rigamma <- function(n, a, b){
  # Uses shape and scale
  return(1/rgamma(n, shape = a, rate = 1/b))
}

# Function to calculate right part of accept. prob.
calc.alpha <- function(t, tau.prop, tau.old){
  log.frac <- y[t]*(log.sigm(tau.prop) - log.sigm(tau.old)) + 
                      (n[t] - y[t])*(log.1m.sigm(tau.prop) - log.1m.sigm(tau.old))
  return(exp(log.frac))
}

acc.prob <- function(t, tau.prop, tau.old){
  log.frac <- y[t]*(log(exp(-tau.old) + 1) - log(exp(-tau.prop) + 1)) + (n[t] - y[t])*(log(exp(tau.old) + 1) - log(exp(tau.prop) + 1))
  return(min(1, exp(log.frac)))
}

alpha.c <- function(t, tau.prop, tau.old){
  right <- y[t]*(tau.prop - tau.old) + n[t]*log((1 + exp(tau.old))/(1 + exp(tau.prop)))
  return(min(1, exp(right)))
}

# Function for calculating 100(1-alpha) credible intervals
CI <- function(samples, burn = 0, alpha = 0.05, by = 1000){
  samples <- tail(samples, length(samples) - burn) # Remove burn-in samples
  
  iter <- length(samples) %/% by
  ci.mat <- data.frame(lower = rep(NA, iter), upper = rep(NA, iter), idx = rep(NA, iter))
  for(i in 1:iter){
    ci <- quantile(samples[1:(i*by)], probs = c(alpha/2, 1 - alpha/2))
    ci.mat$lower[i] <- ci[1]
    ci.mat$upper[i] <- ci[2]
    ci.mat$idx[i] <- i*by
  }
  return(ci.mat)
}

plot.pi <- function(t, samples, burn = 0, alpha = 0.05, ci.by = 1000){
  pi <- y[t]/n[t]
  pi.samples <- sigm(samples)
  pi.mean <- cumsum(pi.samples)/1:length(samples)
  CI <- CI(pi.samples, burn, alpha, ci.by)
  
  plot(1:length(samples), pi.mean, type = "l", ylim = c(0,1))
  lines(CI$idx, CI$upper, col = "red")
  lines(CI$idx, CI$lower, col = "red")
  abline(a = pi, b = 0, col = "blue")
  
  # par(mfrow = c(2, 1))
  # plot(1:length(samples), tau.mean, type = "l", ylim = c(-2, 2))
  # lines(CI$idx, CI$upper, col = "red")
  # lines(CI$idx, CI$lower, col = "red")
}

# Construct Q without the 1/sigma factor
Q <- diag(rep(2, T))
Q[row(Q) - col(Q) == 1] <-  Q[row(Q) - col(Q) == -1] <- -1
Q[1,1] <- Q[T,T] <- 1

# MCMC with iterative conditioning
mcmc.iterative <- function(num.iter, sigma0, tau0){
  tau.mat <- matrix(NA, nrow = num.iter, ncol = T)
  tau.mat[1, ] <- tau0
  sigma.vec <- rep(NA, num.iter)
  sigma.vec[1] <- sigma0
  count <- 0 # Count of accepted tau-samples
  alpha.vec <- rep(NA, num.iter - 1)
  
  for(i in 2:num.iter){
    # tau.mat[i, ] = tau.mat[i - 1, ]
    # Sample tau
    for(j in 1:T){
      # Generate proposal
      mu.cond <- -1/Q[j,j] * Q[j, 1:T != j] %*% (tau.mat[i - 1, 1:T != j])
      Q.cond <-  1/sigma.vec[i-1] * Q[j,j]
      tau.prop <- rnorm(1, mu.cond, sqrt(1/Q.cond))
      
      # Calculate acceptance prob.
      accept.prob <- alpha.c(j, tau.prop, tau.mat[i - 1, j])
      alpha.vec[i-1] <- accept.prob
      u <- runif(1)
      if(u < accept.prob){
        tau.mat[i, j] = tau.prop
        count <- count + 1
      } else{
        tau.mat[i, j] = tau.mat[i - 1, j]
      }
    }
    # Generate Gibbs sample for sigma
    tQt <- sum((head(tau.mat[i, ], -1) - tail(tau.mat[i, ], -1))^2)
    shape <- alpha + (T-1)/2
    scale <- 0.5*tQt + beta
    #sigma.vec[i] <- rinvgamma(1, shape = shape, scale = scale)
    sigma.vec[i] <- 1/rgamma(1, shape = shape, rate = scale)
  }
  return(list(tau.mat = tau.mat, sigma.vec = sigma.vec, count = count, alpha = alpha.vec))
}

# Parameters for sigma prior
alpha <- 2
beta <- 0.05

# testing
mcmc.1k <- mcmc.iterative(50000, sigma0 =  0.02, tau0 = rnorm(T))
mean(mcmc.1k$alpha)
plot(tail(mcmc.1k$sigma.vec, 50000), type = "l")
plot(tail(mcmc.1k$tau.mat[,201], 50000), type = "l")
hist(tail(1/mcmc.1k$sigma.vec, 9000), breaks = 50)
hist(tail(mcmc.1k$tau.mat[,201],40000), freq = FALSE, breaks = 100, add = TRUE)
mcmc.1k$sigma.vec
plot(sigm(mcmc.1k$tau.mat[,1]), type = "l")
abline( a = y[1]/n[1], b = 0, col = "red")

CI(mcmc.1k$tau.mat[,1], burn = 0)
plot.pi(1, mcmc.1k$tau.mat[,5])


# Run MCMC
set.seed(4300)
num.iter <- 50000
ptm <- proc.time() # For computation time
mcmc <- mcmc.iterative(num.iter, sigma0 =  0.02, tau0 = runif(T))
(proc.time() - ptm)[3] # Computation time of MCMC


# Some plotting
plot(1:num.iter, mcmc$sigma.vec, type = "l")
plot(1:num.iter, mcmc$tau.mat[,201], type = "l")
hist(mcmc$tau.mat[,1], breaks = 100, freq = FALSE)
hist(tail(1/mcmc$sigma.vec, 40000), breaks = 100, freq = FALSE)




mcmc.df <- data.frame(x = 1:num.iter, tau1 = mcmc$tau.mat[, 1],
                      tau201 = mcmc$tau.mat[, 201], tau366 = mcmc$tau.mat[, 366], sigma = mcmc$sigma.vec)
# Trace plots
ggplot(mcmc.df) + 
  geom_line(aes(x = x, y = tau1)) + theme_minimal()
ggplot(mcmc.df) + 
  geom_line(aes(x = x, y = sigma)) + theme_minimal()













# Martin's code
library(ggplot2)
library(MASS)
load("rain.rda")

# sample autocorrelation
source('sacf.R')

# "expit"/"sigmoid" function
expit <- function(x){
  1/(exp(-x)+1)
}

# Negative logarithmic "expit"
neg.log.expit <- function(x){
  log(exp(-x)+1)
}

# Function to sample from inverse gamma
rigamma <- function(n, shape, scale){
  # Uses shape and scale
  return(1/rgamma(n, shape = shape, rate = 1/scale))
}

tau.accept <- function(t, tau.prop, tau.curr, y, n){
  log.acc <- y[t] * ( neg.log.expit(tau.curr) - neg.log.expit(tau.prop) ) +
      (y[t] - n[t]) * ( neg.log.expit(-tau.curr) - neg.log.expit(-tau.prop) )
  
  min(1, exp(log.acc))
}

tau.accept.block <- function(I, tau.prop, tau.curr){
  log.acc <- 0
  for(t in I){
    log.acc <- log.acc + y[t] * ( neg.log.expit(tau.curr) - neg.log.expit(tau.prop) ) +
      (y[t] - n[t]) * ( neg.log.expit(-tau.curr) - neg.log.expit(-tau.prop) )
  }
  
  exp(log.acc)
}


MCMCMC <- function(n.iter, tau, sigma2, Q, y, n){
  # length of tau
  T <- length(tau)
  
  # store tau and sigma2 for all the iterations
  tau.mat       <- matrix(NA, nrow = length(tau), ncol = n.iter+1)
  sigma2.vec    <- rep(NA, n.iter+1)
  tau.mat[,1]   = tau
  sigma2.vec[1] = sigma2
  
  for(i in 1:n.iter){
    # MH steps for tau
    for(t in 1:T){
      mt <- (1:T)[-t]
      Q.AA <- Q[t,t]
      Q.AB <- Q[t,mt]
      
      # tau.A conditioned on tau.B
      mu.cond <- -1/Q.AA * Q.AB %*% tau.mat[mt, i]
      Q.cond <- Q.AA
      #print(paste("mu: ", mu.cond))
      #print(paste("Q: ", Q.cond))
      
      # Generate proposal
      tau.proposal <- rnorm(1, mean = mu.cond, sd = sqrt(sigma2.vec[i]/Q.cond))
      
      #print('Iteration:')
      #print(tau.mat[t,i])
      #print(tau.proposal)
      
      # Calculate acceptance probability
      acc <- tau.accept(t, tau.proposal, tau.mat[t,i], y, n)
      if(is.na(acc)){browser()}
      #print(acc)
      
      # Draw from uniform distribution
      u <- runif(1)
      if(u < acc){
        tau.mat[t, i+1] = tau.proposal
      } else{
        tau.mat[t, i+1] = tau.mat[t, i]
      }
    }
    # Gibbs step for sigma2
    shape <- alpha + (T-1)/2
    scale <- 0.5*t(tau.mat[,i]) %*% Q %*% tau.mat[, i] + beta
    sigma2.vec[i+1] <- rigamma(1, shape = shape, scale = scale)
  }
  
  list("tau" = tau.mat, "sigma2" = sigma2.vec)
}


# Test if the acceptance rate is wrong.
do.iterations <- function(n.iter = 50000){
  T <- 366
  n <- rep(39, T)
  n[60] = 10
  #y = rain$n.rain
  y <- rain$n.rain
  
  alpha <- 2
  beta <- 0.05
  
  # Construct Q (without the 1/sigma2 factor)
  Q <- diag(rep(2, T))
  Q[row(Q) - col(Q) == 1] <-  Q[row(Q) - col(Q) == -1] <- -1
  Q[1,1] <- Q[T,T] <- 1
  
  # initial values
  tau.0 <- runif(366, min=-3, max=0)
  sigma2.0 <- 0.2
  
  # do the MCMC calculations
  time <- proc.time()
  run <- MCMCMC(n.iter, tau.0, sigma2.0, Q, y, n)
  print(proc.time()-time)
  
  run
}

run <- do.iterations(n.iter=2)

hist(tail(run$sigma2, 100), freq = FALSE)
plot(run$tau[100,], type = 'l')







