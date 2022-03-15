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
  if(t == 60){
    log.frac = 0
  }
  return(exp(log.frac))
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
  
  for(i in 2:num.iter){
    # Sample tau
    for(j in 1:T){
      # Generate proposal
      Q.AA <- 1/sigma.vec[i-1] * Q[j,j]
      Q.AB <- 1/sigma.vec[i-1] * Q[j, ]
      Q.AB <- Q.AB[1:T != j]
      
      mu.AgB <- -1/Q.AA * Q.AB %*% (tau.mat[i - 1, 1:T != j])
      Q.AgB <- Q.AA
      tau.prop <- rnorm(1, as.numeric(mu.AgB), solve(Q.AgB))
      
      # Calculate acceptance prob.
      alpha <- min(1, calc.alpha(j, tau.prop, tau.mat[i - 1, j]))
      u <- runif(1)
      if(u < alpha){
        tau.mat[i, j] = tau.prop
        count <- count + 1
      } else{
        tau.mat[i, j] = tau.mat[i-1, j]
      }
    }
    # Generate Gibbs sample for sigma
    shape <- alpha + T/2
    scale <- 0.5*t(tau.mat[i, ]) %*% Q %*% tau.mat[i, ] + beta
    sigma.vec[i] <- rinvgamma(1, shape = shape, scale = scale)
  }
  return(list(tau.mat = tau.mat, sigma.vec = sigma.vec, count = count))
}

# Parameters for sigma prior
alpha <- 2
beta <- 0.05

# Run MCMC
set.seed(4300)
num.iter <- 50000
ptm <- proc.time() # For computation time
mcmc <- mcmc.iterative(num.iter, sigma0 =  0.02, tau0 = runif(T))
(proc.time() - ptm)[3] # Computation time of MCMC

# Some plotting
plot(1:num.iter, mcmc$sigma.vec, type = "l")
plot(1:num.iter, mcmc$tau.mat[,300], type = "l")
hist(mcmc$tau.mat[,366], breaks = 100, freq = FALSE, add = TRUE)
hist(mcmc$sigma.vec, breaks = 100, freq = FALSE)


mcmc.df <- data.frame(x = 1:num.iter, tau1 = mcmc$tau.mat[, 1],
                      tau201 = mcmc$tau.mat[, 201], tau366 = mcmc$tau.mat[, 366], sigma = mcmc$sigma.vec)
# Trace plots
ggplot(mcmc.df) + 
  geom_line(aes(x = x, y = tau1)) + theme_minimal()
ggplot(mcmc.df) + 
  geom_line(aes(x = x, y = sigma)) + theme_minimal()

mcmc.1k <- mcmc.iterative(1000, sigma0 =  0.02, tau0 = runif(T))
plot(501:1000, tail(mcmc.1k$sigma.vec, 500), type = "l")
hist(tail(mcmc.1k$sigma.vec, 500), breaks = 20)
mcmc.1k$sigma.vec











# Martin's code
library(ggplot2)
library(MASS)
load("rain.rda")

# sample autocorrelation
source('sacf.R')

T <- 366
n <- rep(39, T)
n[60] = 10
alpha <- 2
beta <- 0.05

# response
y <- rain$n.rain
# Construct Q (without the 1/sigma2 factor)
Q <- diag(rep(2, T))
Q[row(Q) - col(Q) == 1] <-  Q[row(Q) - col(Q) == -1] <- -1
Q[1,1] <- Q[T,T] <- 1

num.iter <- 500 # 50000 will take about 30 min on my system


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

tau.accept <- function(t, tau.prop, tau.curr){
  log.acc <- y[t] * ( neg.log.expit(tau.curr) - neg.log.expit(tau.prop) ) +
      (y[t] - n[t]) * ( neg.log.expit(-tau.curr) - neg.log.expit(-tau.prop) )
  
  exp(log.acc)
}


MCMCMC <- function(n.iter, tau, sigma2){
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
      
      # Generate proposal
      tau.proposal <- rnorm(1, as.numeric(mu.cond), sigma2.vec[i]/Q.cond)
      
      # Calculate acceptance probability
      alpha <- tau.accept(t, tau.proposal, tau.mat[t,i])
      
      # Draw from uniform distribution
      u <- runif(1)
      if(u < alpha){
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

# set some initial values
tau.0 <- runif(366, min=-3, max=0)
pi.0 <- expit(tau.0)
sigma2.0 <- 0.2

# test the MCMC algorithm
time <- proc.time()
run <- MCMCMC(1000, tau.0, sigma2.0)
print(proc.time()-time)


plot(run$sigma, type = 'l')




