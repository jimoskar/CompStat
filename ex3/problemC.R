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

MLE$lambda0
MLE$lambda1


library(latex2exp)


# convergence plots of the MLE estimates
convergence_lambda0 <- ggplot(data.frame(x=1:length(MLE$lambda0)-1, lam0=MLE$lambda0, lam1=MLE$lambda1), aes(x=x)) +
  geom_line(aes(y=lam0), lwd=0.5, color='blue') + xlab(expression(Iteration~t)) + ylab(expression(lambda[0]^{(t)})) + theme_bw()
ggsave("./figures/convergence_lambda0.pdf", plot = convergence_lambda0, width = 5, height = 3)

convergence_lambda1 <- ggplot(data.frame(x=1:length(MLE$lambda0)-1, lam0=MLE$lambda0, lam1=MLE$lambda1), aes(x=x)) +
  geom_line(aes(y=lam1), lwd=0.5, color='blue') + xlab(expression(Iteration~t)) + ylab(expression(lambda[1]^{(t)})) + theme_bw()
ggsave("./figures/convergence_lambda1.pdf", plot = convergence_lambda1, width = 5, height = 3)






