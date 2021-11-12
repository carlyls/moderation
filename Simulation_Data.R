### Simulation of data for moderation project ###

## See Simulation Plan in Moderation Dropbox folder for documentation

# smaller difference
n.trial <- 1000
n.obs <- 2000
reg.coefs <- c(-.5,.7,.4,.2,.4,.4,.2,.2,.2) #check this
mu.trial <- c(1,1)
mu.obs <- c(0,.5)
Sigma.trial <- Sigma.obs <- cbind(c(.5,0),c(0,.5))
p.trial <- 0.5
p.obs <- 0.2

DataGen <- function(n.trial, n.obs, reg.coefs, mu.trial, mu.obs,
                    Sigma.trial, Sigma.obs, p.trial, p.obs) {
  
  ## Generate X1 and X2 for each sample
  tmp.trial <- rmvnorm(n.trial, mean = mu.trial, sigma = Sigma.trial)
  tmp.obs <- rmvnorm(n.obs, mean = mu.obs, sigma = Sigma.obs)

  ## Get X1 and X2 for each sample
  X1.trial <- tmp.trial[,1]
  X1.obs <- tmp.obs[,1]
  
  X2.trial <- tmp.trial[,2]
  X2.obs <- tmp.obs[,2]
  
  ## Generate X3 for each sample
  X3.trial <- sample(0:1, size=n.trial, prob=c(1-p.trial,p.trial), replace=TRUE)
  X3.obs <- sample(0:1, size=n.obs, prob=c(1-p.obs,p.obs), replace=TRUE)
  
  ## Compute treatment assignment probabilities for the observational sample
  lin_pred <- -1/2 + X1.obs/2 + .4*X2.obs + .5X1.obs^2 + .4*X2.obs^2
  pp <- exp(lin_pred)/(1 + exp(lin_pred))
  
  ## Generate treatment assignments for each sample
  A.trial <- sample(0:1, size=n.trial, prob=c(.5,.5), replace=TRUE)
  A.obs <- rep(0, n.obs)
  uu <- runif(n.obs)
  A.obs[uu < pp] <- 1
  
  
  ## Compute outcome probabilities for each person in each trial
  outcome.prob.trial <- rep(0.0, n.trial)
  outcome.prob.obs <- rep(0.0, n.obs)
  for(k in 1:n.trial) {
    dd <- c(1, A.trial[k], X1.trial[k], X2.trial[k], X3.trial[k], (X1.trial[k])^2, (X2.trial[k])^2,
            A.trial[k]*X2.trial[k], A.trial[k]*X3.trial[k])
    lin_pred <- sum(reg.coefs*dd)
    outcome.prob.trial[k] <- exp(lin_pred)/(1 + exp(lin_pred))
  }
  for(k in 1:n.obs) {
    dd <- c(1, A.obs[k], X1.obs[k], X2.obs[k], X3.obs[k], (X1.obs[k])^2, (X2.obs[k])^2, 
            A.obs[k]*X2.obs[k], A.obs[k]*X3.obs[k])
    lin_pred <- sum(reg.coefs*dd)
    outcome.prob.obs[k] <- exp(lin_pred)/(1 + exp(lin_pred))
  }
  
  ## Generate outcomes for each person in each trial
  uu.trial <- runif(n.trial)
  uu.obs <- runif(n.obs)
  y.trial <- rep(0, n.trial)
  y.obs <- rep(0, n.obs)
  
  y.trial[uu.trial < outcome.prob.trial] <- 1
  y.obs[uu.obs < outcome.prob.obs] <- 1
  
  ## Compute true baseline risk for each trial
  # dd <- reg.coefs[1] + reg.coefs[3]*X1.trial + reg.coefs[4]*X2.trial + reg.coefs[7]*U.trial
  # risk.trial <- exp(dd)/(1 + exp(dd))
  # dd <- reg.coefs[1] + reg.coefs[3]*X1.obs + reg.coefs[4]*X2.obs + reg.coefs[7]*U.obs
  # risk.obs <- exp(dd)/(1 + exp(dd))
  # dd <- reg.coefs[1] + reg.coefs[3]*X1.target + reg.coefs[4]*X2.target + reg.coefs[7]*U.target
  # risk.target <- exp(dd)/(1 + exp(dd))
  
  ## Store data for each sample in data frames
  dat.trial <- data.frame(y=y.trial, A=A.trial, X1=X1.trial, X2=X2.trial) #removed risk column
  dat.obs <- data.frame(y=y.obs, A=A.obs, X1=X1.obs, X2=X2.obs)
  dat.target <- data.frame(y=y.target, A=A.target, X1=X1.target, X2=X2.target)
  
  ans <- list()
  ans$dat.trial <- dat.trial
  ans$dat.obs <- dat.obs
  ans$dat.target <- dat.target
  return(ans)
}
