### Problem A ----
library(ggplot2)
library(MASS)
library(rstan)
source('files/probAhelp.R')
source('files/probAdata.R')

## 2. ----

# implementation of the EM algorithm
EM.alg = function(z, u, lambda0, lambda1, tol=1e-8) {
  n <- length(z)
  print(n)
  err <- Inf
  lambda0.iter <- c(lambda0)
  lambda1.iter <- c(lambda1)
  
  t <- 0
  while(err > tol){
    t <- t+1
    
    # add new lambda
    lambda0.iter <- c( lambda0.iter, n/sum(u*z + (1-u)*(1/lambda0.iter[t] - z/(exp(lambda0.iter[t]*z) - 1))) )
    lambda1.iter <- c( lambda1.iter, n/sum((1-u)*z + u*(1/lambda1.iter[t] - z/(exp(lambda1.iter[t]*z) - 1))) )
    
    # calculate maximum error
    err <- max(abs(lambda0.iter[t] - lambda0.iter[t+1]), abs(lambda1.iter[t] - lambda1.iter[t+1]))
  }
  
  return(list(lambda0=lambda0.iter, lambda1=lambda1.iter))
}

# read in the data
z = read.table('files/z.txt')[,1]
u = read.table('files/u.txt')[,1]

# find the MLE of lambda_0 and lambda_1
MLE = EM.alg(z, u, 1, 1)



# TODO: create convergence plots of the MLE estimates

