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

# Function to calculate acceptance probability, works for single and multiple indecies I
acceptance.probability <- function(I, tau.proposal, tau.current){
  log.acc <- y[I]*(tau.proposal - tau.current) +
    n[I]*log((1 + exp(tau.current))/(1 + exp(tau.proposal)))
  return(min(1, exp(log.acc)))
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

# maybe useless function
plot.pi <- function(t, samples, burn = 0, alpha = 0.05, ci.by = 1000){
  pi <- y[t]/n[t]
  pi.samples <- sigm(samples)
  pi.mean <- cumsum(pi.samples)/1:length(samples)
  CI <- CI(pi.samples, burn, alpha, ci.by)
  
  plot(1:length(samples), pi.mean, type = "l", ylim = c(0,1))
  lines(CI$idx, CI$upper, col = "red")
  lines(CI$idx, CI$lower, col = "red")
  abline(a = pi, b = 0, col = "blue")
}

# Function for plotting mcmc preds or making dataset. should maybe split up
plot.preds <- function(tau.mat, burn = 0, alpha = 0.05, plot = TRUE){
  tau.mat <- tau.mat[burn:nrow(tau.mat), ]
  tau.mean <- apply(tau.mat, 2, mean)
  lower <- apply(tau.mat, 2, quantile, probs = c(alpha/2))
  upper <- apply(tau.mat, 2, quantile, probs = c(1 - alpha/2))
  pi.mcmc <- sigm(tau.mean)
  gg.df <- data.frame(x = 1:T, p = y/n, pi = pi.mcmc, 
                      lower = sigm(lower), upper = sigm(upper))
  if(plot){
    ggplot(gg.df, aes(x = x)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(y = pi, color = "pi")) +
      geom_line(aes(y = p, color = "yt/nt"), alpha = 0.1) +
      xlab("Day in year, t") + ylab(expression(pi)) +
      scale_color_manual(name = " ", values = c("yt/nt" = "blue", "MCMC-predictions" = "black")) +
      theme_minimal()
  } else{
    return(gg.df)
  }
}

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
      accept.prob <- acceptance.probability(j, tau.prop, tau.mat[i - 1, j])
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
    tQt <- sum(diff(tau.mat[i, ])^2)
    shape <- alpha + (T-1)/2
    scale <- 0.5*tQt + beta
    #sigma.vec[i] <- rinvgamma(1, shape = shape, scale = scale)
    sigma.vec[i] <- 1/rgamma(1, shape = shape, rate = scale)
  }
  return(list(tau.mat = tau.mat, sigma.vec = sigma.vec, count = count, alpha = alpha.vec))
}


problem.e <- function(){
  # Parameters
  T <- 366          # days in a year
  n <- rep(39, T)   # n in binomial distribution
  n[60] = 10        # corrigate for feb 29th
  y <- rain$n.rain  # response
  
  # Construct Q without the 1/sigma factor
  Q <- diag(rep(2, T))
  Q[row(Q) - col(Q) == 1] <-  Q[row(Q) - col(Q) == -1] <- -1
  Q[1,1] <- Q[T,T] <- 1
  
  # Parameters for sigma prior
  alpha <- 2
  beta <- 0.05
  
  # Run MCMC
  num.iter <- 5000
  set.seed(4300)
  ptm <- proc.time()
  mcmc <- mcmc.iterative(num.iter, sigma0 =  0.02, tau0 = rnorm(T))
  elapsed.time <- (proc.time() - ptm)[3]
  print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds." ))
  
  # Calculate acceptance rate
  acceptance.rate <- mcmc$count/(num.iter * T)
  print(paste("The acceptance rate is", round(acceptance.rate,8) ))
  
  # Plot predictions of pi
  plot.preds(mcmc$tau.mat)
  
  # Calculate statistics for 1, 201, 366 and sigma
  tau.idx <- c(1, 201, 366)
  pi.df <- plot.preds(mcmc$tau.mat, plot = FALSE)
  
  tau.table <- data.frame(idx = tau.idx,
                          pi = pi.df$pi[tau.idx], 
                          lower = pi.df$lower[tau.idx],
                          upper = pi.df$upper[tau.idx])
  tau.table
  
  mean(mcmc$sigma.vec)
  quantile(mcmc$sigma.vec, probs = c(0.025, 0.975))
}



## f) ----

mcmc.block <- function(num.iter, sigma0, tau0, M){
  # Blocking is now possible for M>1
  if(M==1){
    return( mcmc.iterative(num.iter, sigma0, tau0) )
  }
  ## Initialising
  tau.mat <- matrix(NA, nrow = num.iter, ncol = T)
  tau.mat[1, ] <- tau0
  sigma.vec <- rep(NA, num.iter)
  sigma.vec[1] <- sigma0
  count <- 0 # Count of accepted tau-samples
  alpha.vec <- rep(NA, num.iter - 1)
  
  
  ## Precomputing, (assuming M<T)
  n.block.total <- ceiling(T/M)
  n.blocks <- list(1, ceiling(T/M)-2, 1)
  M.3 <- T-(ceiling(T/M)-1)*M
  # Q.AA for the three different blocks
  Q.AA <- list( Q[1:M, 1:M], Q[2:(M+1), 2:(M+1)], Q[(T-M.3+1):T, (T-M.3+1):T] )
  # Inverse of Q.AA for all blocks
  Q.AA.inv <- list( solve(Q.AA[[1]]) )
  if(n.blocks[[2]] > 0){
    Q.AA.inv2 <- solve(Q.AA[[2]])
    for(n in 1:n.blocks[[2]]){
      Q.AA.inv <- append( Q.AA.inv, list(Q.AA.inv2) )
    }  
  } 
  Q.AA.inv <- append( Q.AA.inv, list(solve(Q.AA[[3]])) )
  # Q.AB for all blocks
  Q.AB <- list(Q[1:M, (M+1):T])
  if(n.blocks[[2]] > 0){
    for(n in 1:n.blocks[[2]]){
      I <- (M*n + 1):(M*(n+1))
      Q.AB <- append( Q.AB, list(Q[I, (1:T)[-I]]) )
    }   
  }
  Q.AB <- append( Q.AB, list(Q[(T-M.3+1):T, 1:(T-M.3)]) )
  # inv(Q.AA) * Q.AB
  S <- list()
  for(n in 1:n.block.total){
    S <- append( S, list(Q.AA.inv[[n]] %*% Q.AB[[n]]) )
  }
  
  
  for(i in 2:num.iter){
    # MH samples for tau
    idx.block <- 0
    for(j in seq(1,T,M)){
      # Indecies
      idx.block <- idx.block + 1
      a <- j
      b <- min(j+M-1, T)
      I = a:b
      
      # Generate proposal
      mu.cond <- -S[[idx.block]] %*% matrix(tau.mat[i-1, (1:T)[-I]], ncol=1)
      Q.cond.inv <-  sigma.vec[i-1] * Q.AA.inv[[idx.block]]
      tau.prop <- mvrnorm(1, mu = mu.cond, Sigma = Q.cond.inv)
      
      # Calculate acceptance prob.
      accept.prob <- acceptance.probability(I, tau.prop, tau.mat[i-1, I])
      alpha.vec[i-1] <- accept.prob
      
      u <- runif(1)
      if(u < accept.prob){
        tau.mat[i, I] = tau.prop
        count <- count + 1
      } else{
        tau.mat[i, I] = tau.mat[i - 1, I]
      }
    }
    # Generate Gibbs sample for sigma
    tQt <- sum(diff(tau.mat[i, ])^2)
    shape <- alpha + (T-1)/2
    scale <- 0.5*tQt + beta
    sigma.vec[i] <- 1/rgamma(1, shape = shape, rate = scale)
  }
  return(list(tau.mat = tau.mat, sigma.vec = sigma.vec, count = count, alpha = alpha.vec))
}


problem.f <- function(){
  # Parameters
  T <- 366          # days in a year
  n <- rep(39, T)   # n in binomial distribution
  n[60] = 10        # corrigate for feb 29th
  y <- rain$n.rain  # response
  
  # Construct Q without the 1/sigma factor
  Q <- diag(rep(2, T))
  Q[row(Q) - col(Q) == 1] <-  Q[row(Q) - col(Q) == -1] <- -1
  Q[1,1] <- Q[T,T] <- 1
  
  # Parameters for sigma prior
  alpha <- 2
  beta <- 0.05
  
  # Blocking interval
  M <- 10
  
  # Run MCMC
  num.iter <- 5000
  set.seed(4300)
  ptm <- proc.time()
  mcmc <- mcmc.block(num.iter, sigma0 =  0.02, tau0 = rnorm(T), M)
  elapsed.time <- (proc.time() - ptm)[3]
  print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds." ))
  
  # Calculate acceptance rate
  acceptance.rate <- mcmc$count/(num.iter * T)
  print(paste("The acceptance rate is", round(acceptance.rate,8) ))
  
  # Plot predictions of pi
  plot.preds(mcmc$tau.mat)
  
  # Calculate statistics for 1, 201, 366 and sigma
  tau.idx <- c(1, 201, 366)
  pi.df <- plot.preds(mcmc$tau.mat, plot = FALSE)
  
  tau.table <- data.frame(idx = tau.idx,
                          pi = pi.df$pi[tau.idx], 
                          lower = pi.df$lower[tau.idx],
                          upper = pi.df$upper[tau.idx])
  tau.table
  
  mean(mcmc$sigma.vec)
  quantile(mcmc$sigma.vec, probs = c(0.025, 0.975))
}










# Martin's old code
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





