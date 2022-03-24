### Problem 2 ----
library(INLA)
library(pals)
library(xtable)

load("rain.rda")
alpha <- 2
beta <- 0.05
T <- 366 # days in a year
n <- rain$n.years
y <- rain$n.rain # response
  
## a) ----

control.inla = list(strategy="simplified.laplace", int.strategy="ccd")
ptm <- proc.time() # To calculate computation time
mod <- inla(n.rain ~ -1 + f(day, model="rw1", 
                            hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))), 
                            constr=FALSE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)
(proc.time() - ptm)[3] # Computation time

# Create plot of y/n and INLA-predictions and save
pi.df <- data.frame(pi = sigm(mod$summary.random$day$mean),
                    lower = sigm(mod$summary.random$day$`0.025quant`),
                    upper = sigm(mod$summary.random$day$`0.975quant`),
                    p = y/n, x = 1:T)
                    
pi.preds <- ggplot(pi.df, aes(x = x)) + geom_line(aes(y = pi, color = "pi")) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = p, color = "yt/nt"), alpha = 0.1) +
  xlab("Day in year, t") + ylab(expression(pi)) +
  scale_color_manual(name = " ", values = c("yt/nt" = "blue", "INLA-predictions" = "black"),
                     labels = c(expression(y[t]/n[t]),"INLA-predictions")) + theme_minimal()

ggsave("./figures/pi_preds_inla.pdf", plot = pi.preds, height = 4.0, width = 8.0)


plot(inla.smarginal(mod$marginals.hyperpar$`Precision for day`), type = "l")

# Plot marginals
plot(inla.smarginal(mod$marginals.random$day$index.201), type = "l")
plot(inla.smarginal(mod$marginals.hyperpar$`Precision for day`), type = "l")
  
plot(inla.smarginal(m))

## b) ----
# Plot prediction for 12 different control input combinations
strategy = c('gaussian', 'simplified.laplace', 'laplace', 'adaptive')
int.strategy = c('ccd', 'grid', 'eb')
controls <- expand.grid(strategy, int.strategy)

pal <- unname(tol()) # Color palette
times <- rep(NA, 12) # Computation times
sc <- list() # For ggplot
p <- ggplot() # Init. plot
max.vec <- rep(0, 12) # Vector for maximal value in preds
combs <- 1:12 # Combination index
for(i in combs){
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
  
  preds <- sigm(mod$summary.random$day$mean)
  max.vec[i] <- max(preds)
  p <- p + geom_line(data = data.frame(x = 1:T, y = preds),
                     aes(x, y, color = paste(i)), size = 0.3)
  sc[[i]] = pal[i]
}

# Plotting
names(sc) <- combs
plot.control <- p + scale_color_manual(name = " ",  values = unlist(sc))+
  xlab("Day in year, t")  + ylab(expression(pi)) + theme_minimal()
ggsave("./figures/inla_control.pdf", plot = plot.control, height = 4.0, width = 8.0)


# Create dataframe of information
summary.df <- data.frame(strategy = controls[, 1], 
                   int.strategy = controls[ ,2], 
                   time = times, max = max.vec)

# Create table in LaTeX
xtable(comp.times)

## c) ----

# Model with intercept and constraint 
mod.intercept <- inla(n.rain ~ f(day, model="rw1", 
                        hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))),
                        constr=TRUE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)

# Plot predictions together
df <- data.frame(x = 1:T, pred.a =  sigm(mod$summary.random$day$mean),
                 pred.int =  sigm(mod.intercept$summary.random$day$mean + mod.intercept$summary.fixed$mean))
plot.together <- ggplot(df, aes(x = x)) + geom_line(aes(y = pred.a, color = "Without intercept")) + 
  geom_line(aes(y = pred.int, color = "With intercept")) + 
  scale_color_manual(name = " ", values = c("Without intercept" = "red", "With intercept" = "blue")) + 
  xlab("Day in year, t") + ylab(expression("Predictions of "~pi)) + theme_minimal()
ggsave("./figures/inla_together.pdf", plot = plot.together, height = 4.0, width = 8.0)

# Plot values of tau
df <- data.frame(x = 1:T, pred.a =  mod$summary.random$day$mean,
                 pred.int =  mod.intercept$summary.random$day$mean)
plot.tau <- ggplot(df, aes(x = x)) + geom_line(aes(y = pred.a, color = "Without intercept")) + 
  geom_line(aes(y = pred.int, color = "With intercept")) + 
  scale_color_manual(name = " ", values = c("Without intercept" = "red", "With intercept" = "blue")) + 
  xlab("Day in year, t") + ylab(expression("Predictions of "~tau)) + theme_minimal()
ggsave("./figures/inla_tau.pdf", plot = plot.tau, height = 4.0, width = 8.0)

## Testing ----

plot(inla.smarginal(mod$marginals.random$day$index.366), type = "l")
plot(inla.smarginal(mod$marginals.random$day$index.366), type = "l")
