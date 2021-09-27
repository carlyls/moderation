GenerateSimData <- function(n.trial, n.obs, n.target, reg.coefs, mu.trial, mu.obs, mu.target,
                         Sigma.trial, Sigma.obs, Sigma.target, p.trial, p.obs, p.target) {
     ####################################################################################
     ##  Function to generate data according to the simulation setup described in the methodology paper
     ##
     ##    Input:
     ##      n.trial, n.obs, n.target - the sample sizes in each of the samples      
     ##      reg.coefs - the regression coefficients used in the outcome model
     ##      mu.trial, mu.obs, mu.target - mean vectors for the distribution of (X1, U), each one has length 2
     ##      Sigma.trial, Sigma.obs, Sigma.target - covariance matrices for the distribution of (X1, U), each one is 2x2
     ##      p.trial, p.obs, p.target - Pr(X_2 = 1) for each sample
     ##
     ##    Output:
     ##      A list with the following components
     ##         dat.trial - data.frame containing simulated trial data
     ##         dat.obs - data.frame containing simulated observational data
     ##         dat.target - data.frame containing simulated target sample
     #####################################################################################
     
     ## Generate X1 and U for each sample
     tmp.trial <- rmvnorm(n.trial, mean = mu.trial, sigma = Sigma.trial)
     tmp.obs <- rmvnorm(n.obs, mean = mu.obs, sigma = Sigma.obs)
     tmp.target <- rmvnorm(n.target, mean = mu.target, sigma = Sigma.target)

     ## Get X1 and U for each sample
     X1.trial <- tmp.trial[,1]
     X1.obs <- tmp.obs[,1]
     X1.target <- tmp.target[,1]

     U.trial <- tmp.trial[,2]
     U.obs <- tmp.obs[,2]
     U.target <- tmp.target[,2]

     ## Generate X2 for each sample
     X2.trial <- sample(0:1, size=n.trial, prob=c(p.trial,1 - p.trial), replace=TRUE)
     X2.obs <- sample(0:1, size=n.obs, prob=c(p.obs, 1 - p.obs), replace=TRUE)
     X2.target <- sample(0:1, size=n.target, prob=c(p.target, 1 - p.target), replace=TRUE)

     ## Compute treatment assignment probabilities for the observational sample
     lin_pred <- -1/2 + X1.obs/2 + .4*X2.obs + .4*U.obs
     pp <- exp(lin_pred)/(1 + exp(lin_pred))

     ## Generate treatment assignments for each sample
     A.trial <- sample(0:1, size=n.trial, prob=c(.5,.5), replace=TRUE)
     A.obs <- rep(0, n.obs)
     uu <- runif(n.obs)
     A.obs[uu < pp] <- 1
     A.target <- sample(0:1, size=n.target, prob=c(.5,.5), replace=TRUE)


     ## Compute outcome probabilities for each person in each trial
     outcome.prob.trial <- rep(0.0, n.trial)
     outcome.prob.obs <- rep(0.0, n.obs)
     outcome.prob.target <- rep(0.0, n.target)
     for(k in 1:n.trial) {
         dd <- c(1, A.trial[k], X1.trial[k], X2.trial[k], A.trial[k]*X1.trial[k], A.trial[k]*X2.trial[k], U.trial[k], A.trial[k]*U.trial[k])
         lin_pred <- sum(reg.coefs*dd)
         outcome.prob.trial[k] <- exp(lin_pred)/(1 + exp(lin_pred))
     }
     for(k in 1:n.obs) {
         dd <- c(1, A.obs[k], X1.obs[k], X2.obs[k], A.obs[k]*X1.obs[k], A.obs[k]*X2.obs[k], U.obs[k], A.obs[k]*U.obs[k])
         lin_pred <- sum(reg.coefs*dd)
         outcome.prob.obs[k] <- exp(lin_pred)/(1 + exp(lin_pred))
     }
     for(k in 1:n.target) {
         dd <- c(1, A.target[k], X1.target[k], X2.target[k], A.target[k]*X1.target[k], A.target[k]*X2.target[k], U.target[k], A.target[k]*U.target[k])
         lin_pred <- sum(reg.coefs*dd)
         outcome.prob.target[k] <- exp(lin_pred)/(1 + exp(lin_pred))
     }

     ## Generate outcomes for each person in each trial
     uu.trial <- runif(n.trial)
     uu.obs <- runif(n.obs)
     uu.target <- runif(n.target)
     y.trial <- rep(0, n.trial)
     y.obs <- rep(0, n.obs)
     y.target <- rep(0, n.target)

     y.trial[uu.trial < outcome.prob.trial] <- 1
     y.obs[uu.obs < outcome.prob.obs] <- 1
     y.target[uu.target < outcome.prob.target] <- 1

     ## Compute true baseline risk for each trial
     dd <- reg.coefs[1] + reg.coefs[3]*X1.trial + reg.coefs[4]*X2.trial + reg.coefs[7]*U.trial
     risk.trial <- exp(dd)/(1 + exp(dd))
     dd <- reg.coefs[1] + reg.coefs[3]*X1.obs + reg.coefs[4]*X2.obs + reg.coefs[7]*U.obs
     risk.obs <- exp(dd)/(1 + exp(dd))
     dd <- reg.coefs[1] + reg.coefs[3]*X1.target + reg.coefs[4]*X2.target + reg.coefs[7]*U.target
     risk.target <- exp(dd)/(1 + exp(dd))

     ## Store data for each sample in data frames
     dat.trial <- data.frame(y=y.trial, A=A.trial, X1=X1.trial, X2=X2.trial, risk=risk.trial)
     dat.obs <- data.frame(y=y.obs, A=A.obs, X1=X1.obs, X2=X2.obs, risk=risk.obs)
     dat.target <- data.frame(y=y.target, A=A.target, X1=X1.target, X2=X2.target, risk=risk.target)

     ans <- list()
     ans$dat.trial <- dat.trial
     ans$dat.obs <- dat.obs
     ans$dat.target <- dat.target
     return(ans)
}


