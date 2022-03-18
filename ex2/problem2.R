### Problem 2 ----
library(INLA)
library(pals)

load("rain.rda")
alpha <- 2
beta <- 0.05
T <- 366 # days in a year
n <- rep(39, T) # n in binom. distr.
n[60] = 10 # corrigate for feb 29th
y <- rain$n.rain # response
  
## a) ----
problem.2a <- function(){
  control.inla = list(strategy="simplified.laplace", int.strategy="ccd")
  ptm <- proc.time() # To calculate computation time
  mod <- inla(n.rain ~ -1 + f(day, model="rw1", 
                              hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))), 
                              constr=FALSE),
              data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
              family="binomial", verbose=TRUE, control.inla=control.inla)
  (proc.time() - ptm)[3] # Computation time
  
  # Plot y/n and INLA-predictions
  pi.df <- data.frame(pi = sigm(mod$summary.random$day$mean),
                      lower = sigm(mod$summary.random$day$`0.025quant`),
                      upper = sigm(mod$summary.random$day$`0.975quant`),
                      p = y/n, x = 1:T)
                      
  ggplot(pi.df, aes(x = x)) + geom_line(aes(y = pi, color = "pi")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(aes(y = p, color = "yt/nt"), alpha = 0.1) +
    xlab("Day in year, t") + ylab(expression(pi)) +
    scale_color_manual(name = " ", values = c("yt/nt" = "blue", "INLA-predictions" = "black")) +
    theme_minimal()
}

mod <- inla(n.rain ~ -1 + f(day, model="rw1", 
                            hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))), 
                            constr=FALSE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)
sqrt(1/mod$summary.hyperpar$mean)

problem.2a()
# Plot marginals
plot(inla.smarginal(mod$marginals.random$day$index.201), type = "l")
plot(inla.smarginal(mod$marginals.hyperpar$`Precision for day`), type = "l")
  
plot(inla.smarginal(m))

## b) ----
# Plot prediction for 12 different control input combinations
problem.2b <- function(){
  strategy = c('gaussian', 'simplified.laplace', 'laplace', 'adaptive')
  int.strategy = c('ccd', 'grid', 'eb')
  controls <- expand.grid(strategy, int.strategy)
  
  pal <- unname(tol()) # Color palette
  times <- rep(NA, 12) # Computation times
  sc <- list() # For ggplot
  p <- ggplot() # Init. plot
  for(i in 1:12){
    # Set inputs
    control.inla = list(strategy=paste(controls[i,1]), 
                        int.strategy=paste(controls[i,2]))
    # Fit model and time it
    ptm <- proc.time()
    mod <- inla(n.rain ~ -1 + f(day, model="rw1", 
                                hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))), 
                                constr=FALSE),
                data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
                family="binomial", verbose=TRUE, control.inla=control.inla)
    times[i] <- (proc.time() - ptm)[3]
    
    p <- p + geom_line(data = data.frame(x = 1:T, y = sigm(mod$summary.random$day$mean)),
                       aes(x, y, color = paste(i)), size = 0.1)
    sc[[i]] = pal[i]
  }
  
  names(sc) <- 1:12
  p + scale_color_manual(name = " ",  values = unlist(sc))+
    xlab("t")  + ylab(expression(pi)) + theme_minimal()
  
  comp.times <- data.frame(strategy = controls[, 1], 
                     int.strategy = controls[ ,2 ], 
                     time = times)
}
xtable(comp.times)

plot(inla.smarginal(mod$marginals.random$day$index.366), type = "l")
plot(inla.smarginal(mod$marginals.random$day$index.366), type = "l")

## c) ----

# Model from a)
mod.a <- inla(n.rain ~ -1 + f(day, model="rw1", 
                            hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))), 
                            constr=FALSE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)

# Model with intercept and constraint 
mod.intercept <- inla(n.rain ~ f(day, model="rw1", 
                        hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))),
                        constr=TRUE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)

# Plot predictions together
df <- data.frame(x = 1:T, pred.a =  sigm(mod.a$summary.random$day$mean),
                 pred.int =  sigm(mod.intercept$summary.random$day$mean + mod.intercept$summary.fixed$mean))
ggplot(df, aes(x = x)) + geom_line(aes(y = pred.a, color = "Without intercept")) + 
  geom_line(aes(y = pred.int, color = "With intercept")) + 
  scale_color_manual(name = " ", values = c("Without intercept" = "red", "With intercept" = "blue")) + 
  xlab("Day in year, t") + ylab(expression("Predictions of "~pi)) + theme_minimal()

# Plot values of tau
df <- data.frame(x = 1:T, pred.a =  mod.a$summary.random$day$mean,
                 pred.int =  mod.intercept$summary.random$day$mean)
ggplot(df, aes(x = x)) + geom_line(aes(y = pred.a, color = "Without intercept")) + 
  geom_line(aes(y = pred.int, color = "With intercept")) + 
  scale_color_manual(name = " ", values = c("Without intercept" = "red", "With intercept" = "blue")) + 
  xlab("Day in year, t") + ylab(expression("Predictions of "~tau)) + theme_minimal()

