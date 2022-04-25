### Problem A ----
library(ggplot2)
library(MASS)
source('files/probAhelp.R')
source('files/probAdata.R')

## 1. ----

## Multiple residual bootstraping
ARp.beta.bootstrap <- function(x, p, B){
  # Fit the regression model
  beta <- ARp.beta.est(x, p)
  T <- length(x)
  
  # Find the residuals
  resid.LS <- ARp.resid(x, beta$LS)
  resid.LA <- ARp.resid(x, beta$LA)
  
  # Initialize matrices to store the bootstrap estimates of beta
  betas.LS <- matrix(NA, nrow=B, ncol=p)
  betas.LA <- matrix(NA, nrow=B, ncol=p)
  
  for(b in 1:B){
    # Find consequetive starting values for the resampling
    idx <- sample(1:(T-p+1), 1)
    x0 <- x[idx:(idx+p-1)]
    
    # Resample residuals and x
    resample.resid.LS <- sample(resid.LS, replace = TRUE)
    resample.resid.LA <- sample(resid.LA, replace = TRUE)
    bootstrap.x.LS <- ARp.filter(x0, beta$LS, resample.resid.LS)
    bootstrap.x.LA <- ARp.filter(x0, beta$LA, resample.resid.LA)
    
    # Store new estimates for beta
    betas.LS[b,] <- ARp.beta.est(bootstrap.x.LS, p)$LS
    betas.LA[b,] <- ARp.beta.est(bootstrap.x.LA, p)$LA
  }
  return(list(LS=betas.LS, LA=betas.LA))
}

# Time series of interest
x <- data3A$x
p <- 2
beta <- ARp.beta.est(x, p)

# Plot x
ts <- ggplot(data.frame(x=1:100, y=x), aes(x=x, y=y)) +
  geom_line(lwd=0.5) + labs(x="x", y="y") + theme_bw()
ggsave("./figures/x.pdf", plot = ts, width = 4, height = 2)


B <- 10000
set.seed(4200)
betas <- ARp.beta.bootstrap(x, p, B)

## Estimate variance of the two estimators
var.LS <- apply(betas$LS, 2, var)
var.LA <- apply(betas$LA, 2, var)

# Estimated bias of the two samples using the bootstrap expectation
bias.LS <- apply(betas$LS, 2, mean) - beta$LS
bias.LA <- apply(betas$LA, 2, mean) - beta$LA 



## 2. ----

ARp.pred.bootstrap <- function(x, p, B){
  # Fit the regression model
  beta <- ARp.beta.est(x, p)
  T <- length(x)
  
  # Find the residuals
  resid.LS <- ARp.resid(x, beta$LS)
  resid.LA <- ARp.resid(x, beta$LA)
  
  # Initialize matrices to store the bootstrap estimates of x_{T+1}
  pred.LS <- matrix(NA, nrow=B, ncol=1)
  pred.LA <- matrix(NA, nrow=B, ncol=1)
  
  for(b in 1:B){
    # Find consequetive starting points for the resampling
    idx <- sample(1:(T-p+1), 1, replace = TRUE)
    x0 <- x[idx:(idx+p-1)]
  
    # Resample residuals and x
    resample.resid.LS <- sample(resid.LS, replace = TRUE)
    resample.resid.LA <- sample(resid.LA, replace = TRUE)
    bootstrap.x.LS <- ARp.filter(x0, beta$LS, resample.resid.LS)
    bootstrap.x.LA <- ARp.filter(x0, beta$LA, resample.resid.LA)
    
    # Estimate bootstrap betas
    bootstrap.beta.LS <- ARp.beta.est(bootstrap.x.LS, p)$LS
    bootstrap.beta.LA <- ARp.beta.est(bootstrap.x.LA, p)$LA
    
    # Estimate new distribution for residuals
    bootstrap.resis.LS <- ARp.resid(bootstrap.x.LS, bootstrap.beta.LS)
    bootstrap.resis.LA <- ARp.resid(bootstrap.x.LA, bootstrap.beta.LA)
    
    # Store new estimates for x_{T+1}
    pred.LS[b,] <- bootstrap.beta.LS %*% rev(x[(T-p+1):T]) + sample(bootstrap.resis.LS, 1)
    pred.LA[b,] <- bootstrap.beta.LA %*% rev(x[(T-p+1):T]) + sample(bootstrap.resis.LA, 1)
  }
  return(list(LS=pred.LS, LA=pred.LA))
}

B <- 10000
set.seed(4200)
preds <- ARp.pred.bootstrap(x, p, B)

quantile(preds$LS, probs = c(0.025, 0.975))
quantile(preds$LA, probs = c(0.025, 0.975))
