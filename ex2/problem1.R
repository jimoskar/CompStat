### Problem 1 ----
library(ggplot2)
library(MASS)
library(coda)
source('MCMC.R')
source('sacf.R')


## a) ----
load("rain.rda")

# Applying a monthly moving average filter
rain$smooth.rain <- raster::movingFun(rain$n.rain, n=31, fun=mean, circular=TRUE)

rainfall <- ggplot(rain, aes(x = day)) + geom_point(aes(y = n.rain)) +
  geom_line(aes(y=smooth.rain), colour="blue", size=1) + xlab("Day in year") +
  ylab("Number of days rain exceeds 1mm") + theme_minimal()
rainfall
ggsave("./figures/rainfall.pdf", plot = rainfall, height = 4.0, width = 8.0)

## e) ----
y <- rain$n.rain  # response
n <- rain$n.years # number of years
T <- length(y)    # days in a year (366)

# Parameters for the prior of sigma^2
alpha <- 2.00
beta <- 0.05


# Run the MCMC for 50000 iterations
set.seed(4300)
num.iter    <- 50000
init.tau    <- runif(T,-3,0)
init.sigma2 <- 0.01
ptm <- proc.time()
mcmc <- mcmc.single(num.iter, init.tau, init.sigma2, y, n, alpha, beta)
elapsed.time <- (proc.time() - ptm)[3]


# Elapsed time
print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds."))
# Acceptance rate
print(paste("The acceptance rate is", round(mcmc$acceptance.rate,8) ))

## Traceplots, histograms, and sample autocorrelation functions
mcmc.data.all <- data.frame("x"         = 1:(num.iter+1),
                            "sigma.vec" = mcmc$sigma2,
                            "pi_1"     = sigm(mcmc$tau[,1]),
                            "pi_201"   = sigm(mcmc$tau[,201]),
                            "pi_366"   = sigm(mcmc$tau[,366])
                            )

burn.in = 500
if(burn.in){
  mcmc.data <- mcmc.data.all[-(1:burn.in),]
} else{
  mcmc.data <- mcmc.data.all
}

max.lag <- 50
mcmc.corr <- data.frame("lag"       = 0:max.lag,
                        "sigma.vec" = sacf(mcmc$sigma2, max.lag)$rho.hat,
                        "pi_1"     = sacf(sigm(mcmc$tau[,1]), max.lag)$rho.hat,
                        "pi_201"   = sacf(sigm(mcmc$tau[,201]), max.lag)$rho.hat,
                        "pi_366"   = sacf(sigm(mcmc$tau[,366]), max.lag)$rho.hat)


# Traceplots, histograms, and sample autocorrelation for sigma^2
traceplot.sigma <- ggplot(mcmc.data.all, aes(x=x,y=sigma.vec)) + 
  geom_line() + xlab("Iterations") + ylab(expression(sigma[u]^2)) +
  ylim(0,0.05) + theme_minimal()
traceplot.sigma
ggsave("./figures/traceplot_sigma2.pdf", plot = traceplot.sigma, height = 4.0, width = 8.0)

histogram.sigma <- ggplot(mcmc.data, aes(x=sigma.vec)) + 
  geom_histogram(binwidth=0.001) + xlab(expression(sigma[u]^2)) + ylab("Count") +
  theme_minimal()
histogram.sigma
ggsave("./figures/histogram_sigma2.pdf", plot = histogram.sigma, height = 4.0, width = 8.0)

correlation.sigma <- ggplot(mcmc.corr, aes(x=lag, y=sigma.vec)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  xlab("Lag") + ylab("Correlation") + theme_minimal()
correlation.sigma
ggsave("./figures/correlation_sigma2.pdf", plot = correlation.sigma, height = 4.0, width = 8.0)


# pi_1, pi_201, pi_366
traceplot.tau <- ggplot(mcmc.data, aes(x=x)) +
  geom_line(aes(y=pi_1, colour="tau_1"),size=0.25, alpha=0.4) +
  geom_line(aes(y=pi_201, colour="tau_201"),size=0.25, alpha=0.4) +
  geom_line(aes(y=pi_366, colour="tau_366"),size=0.25, alpha=0.4) +
  scale_color_manual(name="", values=c("tau_1"="red", "tau_201"="blue", "tau_366"="green"),
                     
                      l= expression(pi(tau[1]),pi(tau[201]),pi(tau[366]))) +
  xlab("Iterations") + ylab(" ") + theme_minimal()
traceplot.tau
ggsave("./figures/traceplot_tau.pdf", plot = traceplot.tau, height = 4.0, width = 8.0)

histogram.tau <- ggplot(mcmc.data) + 
  geom_histogram(aes(x=pi_1, fill="tau_1"), binwidth=0.005, alpha=0.6) + 
  geom_histogram(aes(x=pi_201, fill="tau_201"), binwidth=0.005, alpha=0.6) + 
  geom_histogram(aes(x=pi_366, fill="tau_366"), binwidth=0.005, alpha=0.6) + 
  scale_fill_manual(name = " ", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green"),
                    label = expression(pi(tau[1]),pi(tau[201]),pi(tau[366]))) +
  xlab(" ") + ylab("Count") + theme_minimal()
histogram.tau
ggsave("./figures/histogram_tau.pdf", plot = histogram.tau, height = 4.0, width = 8.0)

correlation.tau <- ggplot(mcmc.corr, aes(x=lag)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(y=pi_1, colour="tau_1"),size=0.25) + 
  geom_point(aes(y=pi_1, colour="tau_1"), size =0.3) + 
  geom_line(aes(y=pi_201, colour="tau_201"),size=0.25) +
  geom_point(aes(y=pi_201, colour="tau_201"), size =0.3) +
  geom_line(aes(y=pi_366, colour="tau_366"),size=0.25) +
  geom_point(aes(y=pi_366, colour="tau_366"), size =0.3) + 
  xlab("Lag") + ylab("Correlation") +
  scale_color_manual(name = "", values = c("tau_1" = "red", "tau_201" = "blue", "tau_366" = "green"),
                     label = expression(pi(tau[1]),pi(tau[201]),pi(tau[366]))) +
  xlab("Iteration lag") + ylab("Autcorrelation") + theme_minimal()
correlation.tau
ggsave("./figures/correlation_tau.pdf", plot = correlation.tau, height = 4.0, width = 8.0)


# Plot predictions of pi
pi.preds <- plot.preds(mcmc$tau)
pi.preds
ggsave("./figures/pi_preds.pdf", plot = pi.preds, height = 4.0, width = 8.0)


### Calculate statistics for tau_1, tau_201, tau_366 and sigma
tau.idx <- c(1, 201, 366)
pi.df <- plot.preds(mcmc$tau, plot = FALSE)


tau.table <- data.frame(idx = tau.idx,
                        pi = pi.df$pi[tau.idx], 
                        lower = pi.df$lower[tau.idx],
                        upper = pi.df$upper[tau.idx])
print(tau.table)


sigma.q <- quantile(mcmc$sigma2, probs = c(0.025, 0.975))
sigma.table = c(pred = mean(mcmc$sigma2), lower = sigma.q[1], upper = sigma.q[2])
print(sigma.table)

# Estimated ESS
effectiveSize(as.mcmc(mcmc.data))


## f) ----
y <- rain$n.rain  # response
n <- rain$n.years # number of years
T <- length(y)    # days in a year (366)

# Parameters for the prior of sigma^2
alpha <- 2.00
beta <- 0.05

# Blocking parameter
M <- 10

# Run the MCMC for 50000 iterations
set.seed(4300)
num.iter    <- 50000
init.tau    <- runif(T,-3,0)
init.sigma2 <- 0.01
ptm <- proc.time()
mcmc <- mcmc.block(num.iter, init.tau, init.sigma2, y, n, M, alpha, beta)
elapsed.time <- (proc.time() - ptm)[3]

# Elapsed time
print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds."))
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
                        "sigma.vec" = sacf(mcmc$sigma2, max.lag)$rho.hat,
                        "pi_1"     = sacf(sigm(mcmc$tau[,1]), max.lag)$rho.hat,
                        "pi_201"   = sacf(sigm(mcmc$tau[,201]), max.lag)$rho.hat,
                        "pi_366"   = sacf(sigm(mcmc$tau[,366]), max.lag)$rho.hat
)

# Traceplots, histograms, and sample autocorrelation for sigma^2
traceplot.sigma <- ggplot(mcmc.data.all, aes(x=x,y=sigma.vec)) + 
  geom_line() + xlab("Iterations") + ylab(expression(sigma[u]^2)) + 
  ylim(0, 0.05) + theme_minimal()
traceplot.sigma
ggsave("./figures/traceplot_sigma2_block.pdf", plot = traceplot.sigma, height = 4.0, width = 8.0)

histogram.sigma <- ggplot(mcmc.data, aes(x=sigma.vec)) + 
  geom_histogram(binwidth=0.001) + xlab(expression(sigma[u]^2)) + ylab("Count") +
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


### Calculate statistics for tau_1, tau_201, tau_366 and sigma
tau.idx <- c(1, 201, 366)
pi.df <- plot.preds(mcmc$tau, plot = FALSE)

tau.table <- data.frame(idx = tau.idx,
                        pi = pi.df$pi[tau.idx], 
                        lower = pi.df$lower[tau.idx],
                        upper = pi.df$upper[tau.idx])
print(tau.table)


sigma.q <- quantile(mcmc$sigma2, probs = c(0.025, 0.975))
sigma.table = c(pred = mean(mcmc$sigma2), lower = sigma.q[1], upper = sigma.q[2])
print(sigma.table)  



problem1f <- function(){
  # Estimate computation time/acceptance rate as function of M
  Ms               <- 1:100
  elapsed.times    <- numeric(100)
  acceptance.rates <- numeric(100)
  
  for(M in Ms){
    # Run MCMC
    num.iter <- 1000
    set.seed(4300)
    ptm <- proc.time()
    mcmc <- mcmc.block(num.iter, sigma0 =  0.02, tau0 = runif(T, min=-3, max=0), M = M)
    elapsed.time <- (proc.time() - ptm)[3]
    print(paste("Time elapsed for", num.iter, "iterations is", round(elapsed.time,8), "seconds." ))
    ?runif
    # Calculate acceptance rate
    acceptance.rate <- mcmc$count/((num.iter-1) * ceiling(T/M))
    print(paste("The acceptance rate for M=", M, "is", round(acceptance.rate,8) ))
    
    elapsed.times[M] <- elapsed.time
    acceptance.rates[M] <- acceptance.rate
  }
  
  data.M <- data.frame("M" = Ms,
                       "elapsed.times" = elapsed.times,
                       "acceptance.rates" = acceptance.rates)
  
  plot(elapsed.times, type='l')
  plot(acceptance.rates, type='l')
  
  elapsed.times.plot <- ggplot(data.M, aes(x=M, y=elapsed.times)) + 
    geom_line() + xlab("Block size M") + ylab("Computation time [seconds]") +
    theme_minimal()
  elapsed.times.plot
  ggsave("./figures/elapsed.times.plot.pdf", plot = elapsed.times.plot, height = 2.5, width = 5.0)
  
  acceptance.rates.plot <- ggplot(data.M, aes(x=M, y=acceptance.rates)) + 
    geom_line() + xlab("Block size M") + ylab("Acceptance rates") +
    theme_minimal()
  acceptance.rates.plot
  ggsave("./figures/acceptance.rates.plot.pdf", plot = acceptance.rates.plot, height = 2.5, width = 5.0)
  
  
  # plot relative number of accepted values per calculation time
  rel.acceptance.rate.per.time <- ceiling(T/Ms)/T * acceptance.rates / elapsed.times
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