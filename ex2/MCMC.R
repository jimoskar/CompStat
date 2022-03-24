### Problem 1 ----

# logit function
logit <- function(x){
  return(x/(1 - x))
}

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



if(0){
  load("rain.rda")
  
  y <- rain$n.rain
  n <- rain$n.years
  T <- length(y)
  num.iter <- 50000
  set.seed(4300)
  init.tau    <- runif(T,-3,0)
  init.sigma2 <- 0.01
  
  #library(profvis)
  #profvis( mcmc(num.iter, init.tau, init.sigma2) )
  
  mcmc <- mcmc.single(num.iter, init.tau, init.sigma2)
  
  # Acceptance rate
  print(paste("The acceptance rate is", round(mcmc$acceptance.rate,8) ))
  
  ## Traceplots, histograms, and sample autocorrelation functions
  mcmc.data.all <- data.frame("x"         = 1:(num.iter+1),
                              "sigma.vec" = mcmc$sigma2,
                              "pi_1"     = sigm(mcmc$tau[,1]),
                              "pi_201"   = sigm(mcmc$tau[,201]),
                              "pi_366"   = sigm(mcmc$tau[,366])
  )
  
  burn.in = 50
  if(burn.in){
    mcmc.data <- mcmc.data.all[-(1:burn.in),]
  } else{
    mcmc.data <- mcmc.data.all
  }
  
  max.lag <- 50
  mcmc.corr <- data.frame("lag"       = 0:max.lag,
                          "sigma.vec" = sacf(mcmc$sigma2)$rho.hat,
                          "pi_1"     = sacf(sigm(mcmc$tau[,1]))$rho.hat,
                          "pi_201"   = sacf(sigm(mcmc$tau[,201]))$rho.hat,
                          "pi_366"   = sacf(sigm(mcmc$tau[,366]))$rho.hat
  )
  
  # Traceplots, histograms, and sample autocorrelation for sigma^2
  traceplot.sigma <- ggplot(mcmc.data, aes(x=x,y=sigma.vec)) + 
    geom_line() + xlab("Iterations") + ylab(expression(sigma[u]^2)) +
    theme_minimal()
  traceplot.sigma
  ggsave("./figures/traceplot_sigma2.pdf", plot = traceplot.sigma, height = 4.0, width = 8.0)
  
  histogram.sigma <- ggplot(mcmc.data, aes(x=sigma.vec)) + 
    geom_histogram(binwidth=0.0005) + xlab(expression(sigma[u]^2)) + ylab("Count") +
    theme_minimal()
  histogram.sigma
  ggsave("./figures/histogram_sigma2.pdf", plot = histogram.sigma, height = 4.0, width = 8.0)
  
  correlation.sigma <- ggplot(mcmc.corr, aes(x=lag, y=sigma.vec)) + 
    geom_line() + xlab("Iteration lag") + ylab("Autocorrelation") +
    theme_minimal()
  correlation.sigma
  ggsave("./figures/correlation_sigma2.pdf", plot = correlation.sigma, height = 4.0, width = 8.0)
  
  
  # pi_1, pi_201, pi_366
  traceplot.tau <- ggplot(mcmc.data, aes(x=x)) +
    geom_line(aes(y=pi_1, colour="tau_1"),size=0.25, alpha=0.4) +
    geom_line(aes(y=pi_201, colour="tau_201"),size=0.25, alpha=0.4) +
    geom_line(aes(y=pi_366, colour="tau_366"),size=0.25, alpha=0.4) +
    scale_color_manual(name="", values=c("tau_1"="red", "tau_201"="blue", "tau_366"="green"),
                       labels = expression(pi(tau[1]),pi(tau[201]),pi(tau[366]))) +
    xlab("Iterations") + ylab(" ") + theme_minimal()
  traceplot.tau
  ggsave("./figures/traceplot_tau.pdf", plot = traceplot.tau, height = 4.0, width = 8.0)
  
  histogram.tau <- ggplot(mcmc.data) + 
    geom_histogram(aes(x=pi_1, fill="tau_1"), binwidth=0.005, alpha=0.6) + 
    geom_histogram(aes(x=pi_201, fill="tau_201"), binwidth=0.005, alpha=0.6) + 
    geom_histogram(aes(x=pi_366, fill="tau_366"), binwidth=0.005, alpha=0.6) + 
    scale_fill_manual(name = " ", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green"),
                      labels = expression(pi(tau[1]),pi(tau[201]),pi(tau[366]))) +
    xlab(" ") + ylab("Count") + theme_minimal()
  histogram.tau
  ggsave("./figures/histogram_tau.pdf", plot = histogram.tau, height = 4.0, width = 8.0)
  
  correlation.tau <- ggplot(mcmc.corr, aes(x=lag)) + 
    geom_line(aes(y=pi_1, colour="tau_1"),size=0.25) +
    geom_line(aes(y=pi_201, colour="tau_201"),size=0.25) +
    geom_line(aes(y=pi_366, colour="tau_366"),size=0.25) +
    scale_color_manual(name = "", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green"),
                       labels = expression(pi(tau[1]),pi(tau[201]),pi(tau[366]))) +
    xlab("Iteration lag") + ylab("Autcorrelation") + theme_minimal()
  correlation.tau
  ggsave("./figures/correlation_tau.pdf", plot = correlation.tau, height = 4.0, width = 8.0)
  
  
  # Plot predictions of pi
  pi.preds <- plot.preds(mcmc$tau)
  pi.preds
  ggsave("./figures/pi_preds.pdf", plot = pi.preds, height = 4.0, width = 8.0)
}