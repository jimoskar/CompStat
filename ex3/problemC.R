### Problem A ----
library(ggplot2)
library(MASS)
library(rstan)
library(xtable)
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

print(tail(MLE$lambda0,1)) # 3.465735
print(tail(MLE$lambda1,1)) # 9.353215

df.lambda <- data.frame(x=1:length(MLE$lambda0)-1, lam0=MLE$lambda0, lam1=MLE$lambda1)

# convergence plots of the MLE estimates
convergence_lambda0 <- ggplot(df.lambda, aes(x=x)) +
  geom_line(aes(y=lam0), lwd=0.5, color='blue') + xlab(expression(Iteration~t)) +
  ylab(expression(lambda[0]^{(t)})) + theme_bw()
ggsave("./figures/convergence_lambda0.pdf", plot = convergence_lambda0, width = 5, height = 3)

convergence_lambda1 <- ggplot(df.lambda, aes(x=x)) +
  geom_line(aes(y=lam1), lwd=0.5, color='blue') + xlab(expression(Iteration~t)) +
  ylab(expression(lambda[1]^{(t)})) + theme_bw()
ggsave("./figures/convergence_lambda1.pdf", plot = convergence_lambda1, width = 5, height = 3)

## 3. ----

# Bootstrapping algo. for lambda0 and lambda1
boot.lambda <- function(B, seed = 4300){
  set.seed(seed)
  n <- length(z)
  lambda0.b <- lambda1.b <- rep(NA, B)
  for(i in 1:B){
    sample.idx <- sample(1:n, replace = TRUE)
    zb <- z[sample.idx]
    ub <- u[sample.idx]
    lambdas <- EM.alg(zb, ub, 1, 1)
    lambda0.b[i] <- tail(lambdas$lambda0, 1)
    lambda1.b[i] <- tail(lambdas$lambda1, 1)
  }
  return(list(lambda0 = lambda0.b, lambda1 = lambda1.b))
}

# Run bootstrap
lambdas.b <- boot.lambda(30000)
lambda0.b <- lambdas.b$lambda0
lambda1.b <- lambdas.b$lambda1

# Find standard deviation and bias
sdev <- c(sd(lambda0.b), sd(lambda1.b))
bias <- c(mean(lambda0.b) - tail(MLE$lambda0, 1),
          mean(lambda1.b) - tail(MLE$lambda1, 1))
boot.df <- data.frame(sdev, bias)
xtable(boot.df)

# Find correlation
cor(lambdas.b$lambda0, lambdas.b$lambda1)

## 4. ----

log.lik <- function(lambdas){
  l0 <- lambdas[1]
  l1 <- lambdas[2]
  
  M0 <- which(u == 0)
  M1 <- which(u == 1)
  
  n <- length(z)
  e0 <- exp(-l0*z)
  e1 <- exp(-l1*z)
  return(n*(log(l0*l1) - log(l0 + l1)) + sum(log(e0*(1 - e1) + e1*(1 - e0))))
}

log.lik2 <- function(lambdas){
  idx_u_0 = which(u == 0)
  idx_u_1 = which(u == 1)
  n0 = length(idx_u_0)
  n1 = length(idx_u_1)
  ll = n0 * log(lambdas[2])+n1*log(lambdas[1]) +
    sum(log(1-exp(-lambdas[1]*z[idx_u_0])) - lambdas[2]*z[idx_u_0]) +
    sum(log(1-exp(-lambdas[2]*z[idx_u_1])) - lambdas[1]*z[idx_u_1])
}

mle <- optim(par = c(3,9), fn = log.lik, control = list(fnscale = -1))
mle
