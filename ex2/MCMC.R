### Problem 1 ----

library(pracma)


# sigmoid / expit function
sigm <- function(tau){
  return(exp(tau)/(1 + exp(tau)))
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
      scale_color_manual(name = " ", values = c("yt/nt" = "blue", "MCMC-predictions" = "black"),
                         labels = c(expression(y[t]/n[t]), "MCMC-predictions")) +
      theme_minimal()
    return(plot)
  } else{
    return(gg.df)
  }
}


# Acceptance probability
acceptance.probability <- function(I, tau.proposal, tau.current){
  log.acc <- y[I]*(tau.proposal - tau.current) +
    n[I]*( log(1 + exp(tau.current)) - log(1 + exp(tau.proposal)) )
  
  return( min(1, exp(sum(log.acc))) )
}


mcmc.single <- function(num.iterations, initial.tau, initial.sigma2, y, n, alpha=2.00, beta=0.05){
  ### Initialization
  T            <- length(initial.tau)
  tau          <- matrix(NA, nrow = (num.iterations+1), ncol = T)
  tau[1, ]     <- initial.tau
  sigma2       <- c(initial.sigma2, rep(NA, num.iterations))
  num.accepted <- 0
  alpha.star   <- alpha + (T-1)/2
  root2        <- sqrt(2)
  
  ### Iterations
  for(i in 2:(num.iterations+1)){
    
    ### Metropolis-Hastings steps for tau
    tau[i, ] <- tau[i-1, ]
    # Precompute probabilities
    U <- runif(T)
    Z.sigma <- sqrt(sigma2[i-1]) * rnorm(T)
    # t = 1
    tau.proposal <- tau[i,2] + Z.sigma[1]
    tau.current  <- tau[i,1]
    log.acceptance <- y[1]*(tau.proposal - tau.current) +
      n[1]*log( (1+exp(tau.current))/(1+exp(tau.proposal)) )
    if( U[1] < exp(log.acceptance) ){
      tau[i,1] = tau.proposal
      num.accepted <- num.accepted+1
    }
    # t = 2:T-1
    for(t in 2:(T-1)){
      tau.proposal <- (tau[i,t-1] + tau[i,t+1])/2 + Z.sigma[t]/root2
      tau.current  <- tau[i,t]
      log.acceptance <- y[t]*(tau.proposal - tau.current) +
        n[t]*log( (1+exp(tau.current))/(1+exp(tau.proposal)) )
      if( U[t] < exp(log.acceptance) ){
        tau[i, t] = tau.proposal
        num.accepted <- num.accepted+1
      }
    }
    # t = T
    tau.proposal <- tau[i,T-1] + Z.sigma[T]
    tau.current  <- tau[i,T]
    log.acceptance <- y[T]*(tau.proposal - tau.current) +
      n[T]*log( (1+exp(tau.current))/(1+exp(tau.proposal)) )
    if(U[T] < exp(log.acceptance) ){
      tau[i, T] = tau.proposal
      num.accepted <- num.accepted+1
    }
    
    ### Gibbs step for sigma^2
    sigma2[i] <- 1/rgamma(1, shape = alpha.star, rate = 0.5*sum(diff(tau[i,])^2) + beta)
  }
  
  ### Return tau, sigma^2 and acceptance rates for tau
  list(tau = tau, sigma2 = sigma2, acceptance.rate = num.accepted/(num.iterations*T))
}


# Returns a 2D array of (start, end) indecies for each block
get.indecies <- function(T, M){
  num.blocks <- ceiling(T/M)
  indecies <- matrix(c(seq(1, T, M), seq(M, T+M-1, M)), nrow=num.blocks, ncol=2 )
  indecies[num.blocks, 2] <- T
  return(indecies)
}

mcmc.block <- function(num.iterations, initial.tau, initial.sigma2, y, n, M, alpha=2.00, beta=0.05){
  if(M==1) return( mcmc.single(num.iterations, initial.tau, initial.sigma2, y, n, alpha=2.00, beta=0.05) )
  
  ### Initialization
  T            <- length(initial.tau)
  tau          <- matrix(NA, nrow = (num.iterations+1), ncol = T)
  tau[1, ]     <- initial.tau
  sigma2       <- c(initial.sigma2, rep(NA, num.iterations))
  num.accepted <- 0
  alpha.star   <- alpha + (T-1)/2
  root2        <- sqrt(2)
  
  
  ### Precomputing (assuming M < T)
  indecies <- get.indecies(T,M)
  num.blocks <- dim(indecies)[1]
  M.last <- T - indecies[num.blocks,1]+1
  
  # Create Q matrix
  Q <- diag(c(1, rep(2,T-2),1))
  Q[abs(row(Q) - col(Q))==1] <- -1
  
  
  # Q.AA for the three different blocks
  Q.AA.inv1 <- solve( Q[1:M, 1:M] )
  Q.AA.inv2 <- solve( Q[2:(M+1), 2:(M+1)] )
  Q.AA.inv3 <- solve( Q[indecies[num.blocks,1]:T, indecies[num.blocks,1]:T] )

  
  ### Iterations
  for(i in 2:(num.iter+1)){
    
    ### Metropolis-Hastings steps for tau
    tau[i, ] <- tau[i-1, ]
    sigma2[i] <- sigma2[i-1]
    Q.AA.inv1.sigma2 <- sigma2[i]*Q.AA.inv1
    Q.AA.inv2.sigma2 <- sigma2[i]*Q.AA.inv2
    Q.AA.inv3.sigma2 <- sigma2[i]*Q.AA.inv3
    
    # precompute probabilities
    U <- runif(num.blocks)
    
    # I = (1:M)
    I = 1:M
    tau.proposal <- mvrnorm(1, mu=rep(tau[i,M+1],M), Sigma=Q.AA.inv1.sigma2)
    tau.current  <- tau[i,I]
    log.acceptance <- y[I]*(tau.proposal - tau.current) +
      n[I]*log( (1+exp(tau.current))/(1+exp(tau.proposal)) )
    if( U[1] < exp(sum(log.acceptance)) ){
      tau[i,I] = tau.proposal
      num.accepted <- num.accepted+M
    }
    # I = (M+1:2M), ...
    for(j in 2:(num.blocks-1)){
      
      a <- indecies[j,1]
      b <- indecies[j,2]
      I = a:b
      tau.proposal <- mvrnorm(1, mu=linspace(tau[i,a-1], tau[i,b+1], M+2)[2:(M+1)], Sigma=Q.AA.inv2.sigma2)
      tau.current  <- tau[i,I]
      log.acceptance <- y[I]*(tau.proposal - tau.current) +
        n[I]*log( (1+exp(tau.current))/(1+exp(tau.proposal)) )
      if( U[j] < exp(sum(log.acceptance)) ){
        tau[i,I] = tau.proposal
        num.accepted <- num.accepted+M
      }
    }
    # I = ((T-M.last+1):T)
    I = indecies[num.blocks,1]:T
    tau.proposal <- mvrnorm(1, mu=rep(tau[i,T-M], M.last), Sigma=Q.AA.inv3.sigma2)
    tau.current  <- tau[i,I]
    log.acceptance <- y[I]*(tau.proposal - tau.current) +
      n[I]*log( (1+exp(tau.current))/(1+exp(tau.proposal)) )
    if( U[num.blocks] < exp(sum(log.acceptance)) ){
      tau[i,I] = tau.proposal
      num.accepted <- num.accepted+M.last
    }
    
    ### Gibbs step for sigma^2
    sigma2[i] <- 1/rgamma(1, shape = alpha.star, rate = 0.5*sum(diff(tau[i,])^2) + beta)
  }
  
  ### Return tau, sigma^2 and acceptance rates for tau
  list(tau = tau, sigma2 = sigma2, acceptance.rate = num.accepted/(num.iterations*T))
}


mcmc.block.old <- function(num.iterations, initial.tau, initial.sigma2, y, n, M){
  if(M==1) return( mcmc.iterative(num.iter, sigma0, tau0) )
  
  tau.mat      <- matrix(NA, nrow = num.iterations, ncol = T)
  tau.mat[1, ] <- initial.tau
  sigma.vec    <- rep(NA, num.iterations)
  sigma.vec[1] <- initial.sigma2
  count        <- 0 # Count of accepted tau-samples
  alpha.vec    <- rep(NA, num.iterations - 1)
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
  for(i in 2:num.iterations){
    
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
}



# Code for checking computation times with profvis.
if(0){
  library(profvis)
  load("rain.rda")
  
  y <- rain$n.rain
  n <- rain$n.years
  T <- length(y)
  num.iter <- 50000
  set.seed(4300)
  init.tau    <- runif(T,-3,0)
  init.sigma2 <- 0.01
  Q <- diag(c(1, rep(2,T-2),1))
  Q[abs(row(Q) - col(Q))==1] <- -1
  M <- 10
  
  profvis( mcmc.single(num.iter, init.tau, init.sigma2, y, n) ) # 3 sec
  #profvis( mcmc.block(num.iter, init.tau, init.sigma2, y, n, M) ) # 36 sec
  #rofvis( mcmc.block.old(num.iter, init.tau, init.sigma2, y, n, M) ) # 41 sec
}



