### Problem 1 ----
library(ggplot2)
load("rain.rda")

## a) ----
head(rain)
ggplot(rain, aes(x = day, y = n.rain)) + 
  geom_jitter() + geom_smooth() + 
  theme_minimal()

## e) ----
T <- 366 # days in a year
n <- rep(39, T) # n in binom. distr.
n[60] = 10 # corrigate for feb 29th

# Construct Q without the 1/sigma factor
Q <- diag(rep(2, T))
Q[row(Q) - col(Q) == 1] <-  Q[row(Q) - col(Q) == -1] <- -1
Q[1,1] <- Q[T,T] <- 1

mcmc.iterative <- function(iter, sigma0, tau0){
  tau.mat <- matrix(NA, nrow = iter, ncol = T)
  tau.mat[0, ] <- tau0
  sigma.vec <- rep(NA, iter)
  sigma.vec[0] <- sigma0
  for(i in 1:iter){
    # sample tau
    for(j in 1:T){
      Q.AA <- 1/sigma * Q[j,j]
      temp <- rep(0, T-1)
      temp[j] <- -1
      Q.AB<- 1/sigma * t(temp) 
      mu.AgB <- 1/Q.AA*Q.AB(tau[1:T != j])
      Q.AgB <- Q.AA
      tau[j] <- rnorm(mu.AgB, solve(Q.AgB))
    }
    sigma <- 
  }
}

