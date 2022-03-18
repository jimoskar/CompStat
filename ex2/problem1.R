### Problem 1 ----
library(ggplot2)
library(MASS)
load("rain.rda")
source('sacf.R')

## a) ----

# Applying a monthly moving average filter
rain$smooth.rain <- raster::movingFun(rain$n.rain, n=31, fun=mean, circular=TRUE)

rainfall <- ggplot(rain, aes(x = day)) + geom_point(aes(y = n.rain)) +
  geom_line(aes(y=smooth.rain), colour="cadetblue", size=1) + xlab("Day in year") +
  ylab("Number of days rain exceeds 1mm") + theme_minimal()
rainfall
ggsave("./figures/rainfall.pdf", plot = rainfall, height = 4.0, width = 8.0)

## e) ----

# logit function
logit <- function(x){
  return(x/(1 - x))
}

# Sigmoid function
sigm <- function(tau){
  return(exp(tau)/(1 + exp(tau)))
}

# Function to calculate acceptance probability, works for single and multiple indecies I
acceptance.probability <- function(I, tau.proposal, tau.current){
  log.acc <- y[I]*(tau.proposal - tau.current) +
    n[I]*log((1 + exp(tau.current))/(1 + exp(tau.proposal)))
  return( min(1, exp(sum(log.acc))) )
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
    plot <- ggplot(gg.df, aes(x = x)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(y = pi, color = "pi")) +
      geom_line(aes(y = p, color = "yt/nt"), alpha = 0.1) +
      xlab("Day in year") + ylab(expression(pi)) +
      scale_color_manual(name = " ", values = c("yt/nt" = "blue", "MCMC-predictions" = "black")) +
      theme_minimal()
    return(plot)
  } else{
    return(gg.df)
  }
}

# Construct Q without the 1/sigma factor
Q <- diag(rep(2, T))
Q[abs(row(Q) - col(Q)) == 1] <- -1
Q[1,1] <- Q[T,T] <- 1

# MCMC with iterative conditioning
mcmc.iterative <- function(num.iter, sigma0, tau0){
  tau.mat      <- matrix(NA, nrow = num.iter, ncol = T)
  tau.mat[1, ] <- tau0
  sigma.vec    <- rep(NA, num.iter)
  sigma.vec[1] <- sigma0
  count        <- 0 # Count of accepted tau-samples
  alpha.vec    <- rep(NA, num.iter - 1)
  
  # Iterations
  for(i in 2:num.iter){
    # Sample tau
    for(t in 1:T){
      # Generate proposal
      mu.cond <- -1/Q[t,t] * Q[t, 1:T != t] %*% (tau.mat[i-1, 1:T != t])
      Q.cond <-  1/sigma.vec[i-1] * Q[t,t]
      tau.prop <- rnorm(1, mu.cond, sqrt(1/Q.cond))
      
      # Calculate acceptance probability
      accept.prob <- acceptance.probability(t, tau.prop, tau.mat[i-1, t])
      alpha.vec[i-1] <- accept.prob
      u <- runif(1)
      if(u < accept.prob){
        tau.mat[i, t] = tau.prop
        count <- count + 1
      } else{
        tau.mat[i, t] = tau.mat[i-1, t]
      }
    }
    # Generate Gibbs sample for sigma
    tQt <- sum(diff(tau.mat[i, ])^2)
    shape <- alpha + (T-1)/2
    scale <- 0.5*tQt + beta
    sigma.vec[i] <- 1/rgamma(1, shape = shape, rate = scale)
  }
  
  return( list(tau.mat = tau.mat, sigma.vec = sigma.vec, count = count, alpha = alpha.vec) )

}

problem.e <- function(){
  y <- rain$n.rain  # response
  n <- rain$n.years # number of years
  T <- length(y)    # days in a year (366)
  
  # Construct the precision matrix Q (without the 1/sigma^2_u factor)
  Q <- diag(c(1, rep(2, T-2), 1))
  Q[abs(row(Q) - col(Q)) == 1] <- -1

  # Parameters for the prior of sigma^2
  alpha <- 2
  beta <- 0.05
  
  # Run the MCMC
  num.iter <- 50000
  set.seed(4300)
  ptm <- proc.time()
  mcmc <- mcmc.iterative(num.iter, sigma0 =  0.02, tau0 = rnorm(T))
  elapsed.time <- (proc.time() - ptm)[3]
  
  # Elapsed time
  print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds." ))
  
  # Acceptance rate
  acceptance.rate <- mcmc$count/(num.iter * T)
  print(paste("The acceptance rate is", round(acceptance.rate,8) ))
  
  ## Traceplots, histograms, and sample autocorrelation functions
  mcmc.data <- data.frame("x"         = 1:num.iter,
                          "sigma.vec" = mcmc$sigma.vec,
                          "tau_1"     = mcmc$tau.mat[,1],
                          "tau_201"   = mcmc$tau.mat[,201],
                          "tau_366"   = mcmc$tau.mat[,366]
  )
  max.lag <- 20
  mcmc.corr <- data.frame("lag"       = 0:max.lag,
                          "sigma.vec" = sacf(mcmc$sigma.vec)$rho.hat,
                          "tau_1"     = sacf(mcmc$tau.mat[,1])$rho.hat,
                          "tau_201"   = sacf(mcmc$tau.mat[,201])$rho.hat,
                          "tau_366"   = sacf(mcmc$tau.mat[,366])$rho.hat
                          )
  
  # Traceplots, histograms, and sample autocorrelation for sigma^2
  traceplot.sigma <- ggplot(mcmc.data, aes(x=x,y=sigma.vec)) + 
    geom_line() + xlab("Iterations") + ylab("Variance") +
    theme_minimal()
  ggsave("./figures/traceplot_sigma2.pdf", plot = traceplot.sigma, height = 4.0, width = 8.0)
  
  histogram.sigma <- ggplot(mcmc.data, aes(x=sigma.vec)) + 
    geom_histogram(binwidth=0.01) + xlab("Variance") + ylab("Count") +
    theme_minimal()
  ggsave("./figures/histogram_sigma2.pdf", plot = histogram.sigma, height = 4.0, width = 8.0)
  
  correlation.sigma <- ggplot(mcmc.corr, aes(x=lag, y=sigma.vec)) + 
    geom_line() + xlab("Iteration lag") + ylab("Autocorrelation") +
    theme_minimal()
  ggsave("./figures/correlation_sigma2.pdf", plot = correlation.sigma, height = 4.0, width = 8.0)
  
  
  # tau_1, tau_201, tau_366
  traceplot.tau <- ggplot(mcmc.data, aes(x=x)) +
    geom_line(aes(y=tau_1, colour="tau_1"),size=0.25, alpha=0.6) +
    geom_line(aes(y=tau_201, colour="tau_201"),size=0.25, alpha=0.6) +
    geom_line(aes(y=tau_366, colour="tau_366"),size=0.25, alpha=0.6) +
    scale_color_manual(name="", values=c("tau_1"="red", "tau_201"="blue", "tau_366"="green")) +
    xlab("Iterations") + ylab("Tau") + xlim(0, 1000) + theme_minimal()
  traceplot.tau
  ggsave("./figures/traceplot_tau.pdf", plot = traceplot.tau, height = 4.0, width = 8.0)
  
  histogram.tau <- ggplot(mcmc.data) + 
    geom_histogram(aes(x=tau_1, fill="tau_1"), binwidth=0.05, alpha=0.6) + 
    geom_histogram(aes(x=tau_201, fill="tau_201"), binwidth=0.05, alpha=0.6) + 
    geom_histogram(aes(x=tau_366, fill="tau_366"), binwidth=0.05, alpha=0.6) + 
    scale_color_manual(name = "", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green")) +
    xlab("Tau") + ylab("Count") + theme_minimal()
  ggsave("./figures/histogram_tau.pdf", plot = histogram.tau, height = 4.0, width = 8.0)
  
  correlation.tau <- ggplot(mcmc.corr, aes(x=lag)) + 
    geom_line(aes(y=tau_1, colour="tau_1"),size=0.25) +
    geom_line(aes(y=tau_201, colour="tau_201"),size=0.25) +
    geom_line(aes(y=tau_366, colour="tau_366"),size=0.25) +
    scale_color_manual(name = "", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green")) +
    xlab("Iteration lag") + ylab("Autcorrelation") + theme_minimal()
  ggsave("./figures/correlation_tau.pdf", plot = correlation.tau, height = 4.0, width = 8.0)
  
  
  # Plot predictions of pi
  pi.preds <- plot.preds(mcmc$tau.mat)
  ggsave("./figures/pi_preds.pdf", plot = pi.preds, height = 4.0, width = 8.0)
  
  # Calculate statistics for tau_1, tau_201, tau_366 and sigma
  tau.idx <- c(1, 201, 366)
  pi.df <- plot.preds(mcmc$tau.mat, plot = FALSE)
  
  tau.table <- data.frame(idx = tau.idx,
                          pi = pi.df$pi[tau.idx], 
                          lower = pi.df$lower[tau.idx],
                          upper = pi.df$upper[tau.idx])
  round(tau.table,4)
  
  mean(mcmc$sigma.vec)
  quantile(mcmc$sigma.vec, probs = c(0.025, 0.975))
}


## f) ----

mcmc.block <- function(num.iter, sigma0, tau0, M){
  if(M==1) return( mcmc.iterative(num.iter, sigma0, tau0) )
  
  tau.mat      <- matrix(NA, nrow = num.iter, ncol = T)
  tau.mat[1, ] <- tau0
  sigma.vec    <- rep(NA, num.iter)
  sigma.vec[1] <- sigma0
  count        <- 0 # Count of accepted tau-samples
  alpha.vec    <- rep(NA, num.iter - 1)
  n.blocks     <- ceiling(T/M) # Total number of blocks in one iteration
  n.blocks.mid <- ceiling(T/M)-2
  M.last       <- T-(ceiling(T/M)-1)*M
  
  
  # Precomputing (assuming M < T)
  # Q.AA for the three different blocks
  Q.AA <- list( Q[1:M, 1:M], Q[2:(M+1), 2:(M+1)], Q[(T-M.last+1):T, (T-M.last+1):T] )
  # Inverse of Q.AA for all blocks
  Q.AA.inv <- list( solve(Q.AA[[1]]) )
  if(n.blocks.mid){
    Q.AA.inv2 <- solve(Q.AA[[2]])
    for(n in 1:n.blocks.mid){
      Q.AA.inv <- append( Q.AA.inv, list(Q.AA.inv2) )
    }  
  } 
  Q.AA.inv <- append( Q.AA.inv, list(solve(Q.AA[[3]])) )
  # Q.AB for all blocks
  Q.AB <- list(Q[1:M, (M+1):T])
  if(n.blocks.mid){
    for(n in 1:n.blocks.mid){
      I <- (M*n + 1):(M*(n+1))
      Q.AB <- append( Q.AB, list(Q[I, (1:T)[-I]]) )
    }   
  }
  Q.AB <- append( Q.AB, list(Q[(T-M.last+1):T, 1:(T-M.last)]) )
  # S = -inv(Q.AA) * Q.AB for all blocks
  S <- list()
  for(n in 1:n.blocks){
    S <- append( S, list(Q.AA.inv[[n]] %*% Q.AB[[n]]) )
  }
  
  # Iterations
  for(i in 2:num.iter){
    
    # MH samples for tau
    idx.block <- 0
    for(j in seq(1,T,M)){
      idx.block <- idx.block + 1
      I = j:min(j+M-1, T)
      
      # Generate proposal
      mu.cond <- -S[[idx.block]] %*% matrix(tau.mat[i-1, (1:T)[-I]], ncol=1)
      Q.cond.inv <-  sigma.vec[i-1] * Q.AA.inv[[idx.block]]
      tau.prop <- mvrnorm(1, mu = mu.cond, Sigma = Q.cond.inv)
      
      # Calculate acceptance prob.
      accept.prob <- acceptance.probability(I, tau.prop, tau.mat[i-1, I])
      
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
  
  return( list(tau.mat = tau.mat, sigma.vec = sigma.vec, count = count, alpha = alpha.vec) )



problem.f <- function(){
  y <- rain$n.rain  # response
  n <- rain$n.years # number of years
  T <- length(y)    # days in a year (366)
  
  # Construct the precision matrix Q (without the 1/sigma^2_u factor)
  Q <- diag(c(1, rep(2, T-2), 1))
  Q[abs(row(Q) - col(Q)) == 1] <- -1
  
  # Parameters for the prior of sigma^2
  alpha <- 2
  beta <- 0.05
  
  # Blocking interval length
  M <- 10
  
  # Run the MCMC
  num.iter <- 50000
  set.seed(4300)
  ptm <- proc.time()
  mcmc <- mcmc.block(num.iter, sigma0 =  0.02, tau0 = rnorm(T), M=10)
  elapsed.time <- (proc.time() - ptm)[3]
  
  # Elapsed time
  print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds." ))
  
  # Acceptance rate
  acceptance.rate <- mcmc$count/(num.iter * T)
  print(paste("The acceptance rate is", round(acceptance.rate,8) ))
  
  ## Traceplots, histograms, and estimated autocorrelation functions
  mcmc.data.all <- data.frame(x         = 1:num.iter,
                          sigma.vec = mcmc$sigma.vec[],
                          tau_1     = mcmc$tau.mat[,1],
                          tau_201   = mcmc$tau.mat[,201],
                          tau_366   = mcmc$tau.mat[,366]
  )
  burn.in = 100
  if(burn.in){
    mcmc.data <- mcmc.data.all[-(1:burn.in),]
  } else{
    mcmc.data <- mcmc.data.all
  }
  
  max.lag <- 20
  mcmc.corr <- data.frame(lag       = 0:max.lag,
                          sigma.vec = sacf(mcmc.data$sigma.vec)$rho.hat,
                          tau_1     = sacf(mcmc.data$tau_1)$rho.hat,
                          tau_201   = sacf(mcmc.data$tau_201)$rho.hat,
                          tau_366   = sacf(mcmc.data$tau_366)$rho.hat
  )
  
  # sigma^2
  traceplot.sigma <- ggplot(mcmc.data, aes(x=x,y=sigma.vec)) + 
    geom_line() + xlab("Iterations") + ylab("Variance") +
    theme_minimal()
  #traceplot.sigma
  ggsave("./figures/traceplot_sigma2_block.pdf", plot = traceplot.sigma, height = 4.0, width = 8.0)
  
  histogram.sigma <- ggplot(mcmc.data, aes(x=sigma.vec)) + 
    geom_histogram(binwidth=0.002) + xlab("Variance") + ylab("Count") +
    theme_minimal()
  #histogram.sigma
  ggsave("./figures/histogram_sigma2_block.pdf", plot = histogram.sigma, height = 4.0, width = 8.0)
  
  correlation.sigma <- ggplot(mcmc.corr, aes(x=lag, y=sigma.vec)) + 
    geom_line() + xlab("Lag") + ylab("Correlation") +
    theme_minimal()
  #correlation.sigma
  ggsave("./figures/correlation_sigma2_block.pdf", plot = correlation.sigma, height = 4.0, width = 8.0)
  
  
  # tau_1, tau_201, tau_366
  traceplot.tau <- ggplot(mcmc.data, aes(x=x)) +
    geom_line(aes(y=tau_1, colour="tau_1"),size=0.25, alpha=0.6) +
    geom_line(aes(y=tau_201, colour="tau_201"),size=0.25, alpha=0.6) +
    geom_line(aes(y=tau_366, colour="tau_366"),size=0.25, alpha=0.6) +
    scale_color_manual(name = "", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green")) +
    xlab("Iterations") + ylab("Tau") + xlim(0, 1000) + theme_minimal()
  #traceplot.tau
  ggsave("./figures/traceplot_tau_block.pdf", plot = traceplot.tau, height = 4.0, width = 8.0)
  
  histogram.tau <- ggplot(mcmc.data) + 
    geom_histogram(aes(x=tau_1, fill="tau_1"), binwidth=0.05, alpha=0.6) + 
    geom_histogram(aes(x=tau_201, fill="tau_201"), binwidth=0.05, alpha=0.6) + 
    geom_histogram(aes(x=tau_366, fill="tau_366"), binwidth=0.05, alpha=0.6) + 
    scale_color_manual(name=" ", values=c("tau_1"="red", "tau_201"="blue", "tau_366"="green")) +
    xlab("Tau") + ylab("Count") + theme_minimal()
  #histogram.tau
  ggsave("./figures/histogram_tau_block2.pdf", plot = histogram.tau, height = 4.0, width = 8.0)
  
  correlation.tau <- ggplot(mcmc.corr, aes(x=lag)) + 
    geom_line(aes(y=tau_1, colour="tau_1"),size=0.25) +
    geom_line(aes(y=tau_201, colour="tau_201"),size=0.25) +
    geom_line(aes(y=tau_366, colour="tau_366"),size=0.25) +
    scale_color_manual(name = "", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green")) +
    xlab("Lag") + ylab("Correlation") + theme_minimal()
  #correlation.tau
  ggsave("./figures/correlation_tau_block.pdf", plot = correlation.tau, height = 4.0, width = 8.0)
  
  
  # Estimate computation time/acceptance rate as function of M
  Ms <- 1:20
  elapsed.times <- numeric(20)
  acceptance.probabilities <- numeric(20)
  
  for(M in Ms){
    # Run MCMC
    num.iter <- 1000
    set.seed(4300)
    ptm <- proc.time()
    mcmc <- mcmc.block(num.iter, sigma0 =  0.02, tau0 = rnorm(T), M = M)
    elapsed.time <- (proc.time() - ptm)[3]
    print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds." ))
    
    # Calculate acceptance rate
    acceptance.rate <- mcmc$count/(num.iter * T)
    print(paste("The acceptance rate is", round(acceptance.rate,8) ))
    
    elapsed.times[M] <- elapsed.time
    acceptance.probabilities[M] <- acceptance.rate
  }
  
  plot(elapsed.times[-1], type='l')
  plot(acceptance.probabilities[-1], type='l')
  
  # plot relative number of accepted values per calculation time
  rel.acceptance.rate.per.time <- acceptance.probabilities / elapsed.times
  rel.acceptance.rate.per.time <- rel.acceptance.rate.per.time / max(rel.acceptance.rate.per.time)
  
  plot(rel.acceptance.rate.per.time, type='l')
  
  # Plot predictions of pi
  plot.preds(mcmc$tau.mat)
  
  # Calculate statistics for 1, 201, 366 and sigma
  tau.idx <- c(1, 201, 366)
  pi.df <- plot.preds(mcmc$tau.mat, plot = FALSE)
  
  tau.table <- data.frame(idx   = tau.idx,
                          pi    = pi.df$pi[tau.idx], 
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





