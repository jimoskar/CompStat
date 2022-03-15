### Problem 2 ----

library(INLA)
load("rain.rda")
alpha <- 2
beta <- 0.05
  
  
control.inla = list(strategy="simplified.laplace", int.strategy="ccd")
ptm <- proc.time() # To calculate computation time
mod <- inla(n.rain ~ -1 + f(day, model="rw1", hyper = list(theta = list(prior="loggamma", param=c(alpha, beta))), constr=FALSE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)
(proc.time() - ptm)[3] # Computation time

mod$summary.random

plot(inla.smarginal(mod$marginals.random$day$index.366), type = "l")
plot(inla.smarginal(mod$marginals.hyperpar$`Precision for day`), type = "l")
  
plot(inla.smarginal(m))

#Run inla.doc("rw1") for documentation provided by INLA on its built-in RW(1) model
inla.doc("rw1")


## a) ----


## b) ----


## c) ----

# We consider the following model in INLA:
mod <- inla(n.rain ~ f(day, model="rw1", constr=TRUE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)

summary(mod)

