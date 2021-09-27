
#################################
## Code to perform simulation study 1 and to produce the results presented in
## the Table describing the results of simulation study 1
###############################

library(crampaper)
library(xtable)

### Define parameters used in Simulation Study 1:
Sigma.trial <- Sigma.obs <- cbind(c(.5,0),c(0,.5))
Sigma.target <- cbind(c(1,0),c(0,1.5))
mu.trial <- c(1,1)
mu.obs <- c(1,.5)
mu.target <- c(0,0.5)
p.trial <- 1/2
p.obs <- .8
p.target <- .2
reg.coefs <- c(-5,7,4,4,2,2,2,1)/10
n.trial <- n.obs <- n.target <- 1000

set.seed(4156)
## Number of simulation replications ... 1,000 used for paper results - 1,000 reps may take a few hours to finish.
nreps <- 1000

## setup storage
naive_est <- true_est <- rs_est <- rs_est_trim <- css_est <- css_est_trim <- cram_est <- cram_est2 <- cram_est_trim <- rep(0.0, nreps)

for(k in 1:nreps) {
  ## Generate simulated data sets (trial, observational, and target)
  sim_dat <- GenerateSimData(n.trial, n.obs, n.target, reg.coefs, mu.trial,
                             mu.obs, mu.target, Sigma.trial, Sigma.obs,
                             Sigma.target, p.trial, p.obs, p.target)

  ## Compute Naive estimate
  naive_model <- glm(y ~ A, data = sim_dat$dat.trial, family=binomial)
  naive_est[k] <- naive_model$coefficients[2]

  ## Compute "true" estimate
  true_model <- glm(y ~ A, data = sim_dat$dat.target, family=binomial)
  true_est[k] <- true_model$coefficients[2]

  ## Compute baseline risk estimates and add them to the simulated data
  br.est <- BaselineRiskEstimate(sim_dat$dat.trial, sim_dat$dat.obs,
                                 sim_dat$dat.target)
  sim_dat$dat.trial$est.base.risk <- br.est$brisk.trial
  sim_dat$dat.obs$est.base.risk <- br.est$brisk.obs
  sim_dat$dat.target$est.base.risk <- br.est$brisk.target

  ## Compute risk-standardized estimate of treatment effect
  rs_est[k] <- RiskStdEstimate(sim_dat$dat.trial, sim_dat$dat.target)
  rs_est_trim[k] <- RiskStdEstimate(sim_dat$dat.trial, sim_dat$dat.target, trim.weights=TRUE)
  css_est[k] <- RiskStdEstimate(sim_dat$dat.trial, sim_dat$dat.target, method="CS")
  css_est_trim[k] <- RiskStdEstimate(sim_dat$dat.trial, sim_dat$dat.target, trim.weights=TRUE, method="CS")

  ## Compute targets from the trial data
  sim_dat$dat.trial$y <- as.integer(sim_dat$dat.trial$y)
  targs <- ComputeTargets(sim_dat$dat.trial, sim_dat$dat.obs)

  ## Calibration: ompute optimal U-parameters from the observational data
  vy <- var(sim_dat$dat.obs$y)
  bdds <- FindInterval(sigma=1/2, vy=vy, cor.target=.7)
  opt_val <- optimize(CRAMU, interval=bdds, UPar2=1/2, data=sim_dat$dat.obs, target=targs)

  ## Finally, compute the cram projection estimates of trt. effect in target sample
  cram_est[k] <- CramProject(b.y=opt_val$minimum, sigma=1/2, dat.obs=sim_dat$dat.obs, dat.target=sim_dat$dat.target)
  cram_est_trim[k] <- CramProject(b.y=opt_val$minimum, sigma=1/2, dat.obs=sim_dat$dat.obs, dat.target=sim_dat$dat.target, trim.weights=TRUE)

  print(k)
}

## remove indices where cram_est returned a missing value
ind <- which(is.na(cram_est))


### Record mean and sd results for the three methods
ResultMat <- matrix(NA, nrow=8, ncol=3)
rownames(ResultMat) <- c("Naive","Risk Standardized", "Trimmed Risk Standardized", "Covariate Standardized","Trimmed Covariate Standardized","CRAM", "Trimmed CRAM", "Truth")
colnames(ResultMat) <- c("Mean", "SD", "MSE")
if(length(ind) > 0) {
   mean_true_est <- mean(true_est[-ind])
   mse_naive <- mean((naive_est[-ind] - mean_true_est)*(naive_est[-ind] - mean_true_est))
   mse_rsest <- mean((rs_est[-ind] - mean_true_est)*(rs_est[-ind] - mean_true_est))
   mse_cssest <- mean((css_est[-ind] - mean_true_est)*(css_est[-ind] - mean_true_est))
   mse_cssest_trim <- mean((css_est_trim[-ind] - mean_true_est)*(css_est_trim[-ind] - mean_true_est))
   mse_rs_trim <- mean((rs_est_trim[-ind] - mean_true_est)*(rs_est_trim[-ind] - mean_true_est))
   mse_cram <- mean((cram_est[-ind] - mean_true_est)*(cram_est[-ind] - mean_true_est))
   mse_cram_trim <- mean((cram_est_trim[-ind] - mean_true_est)*(cram_est_trim[-ind] - mean_true_est))
   mse_true <- mean((true_est[-ind] - mean_true_est)*(true_est[-ind] - mean_true_est))

   ResultMat[1,] <- c(mean(naive_est[-ind]), sd(naive_est[-ind]), mse_naive)
   ResultMat[2,] <- c(mean(rs_est[-ind]), sd(rs_est[-ind]), mse_rsest)
   ResultMat[3,] <- c(mean(rs_est_trim[-ind]), sd(rs_est_trim[-ind]), mse_rs_trim)
   ResultMat[4,] <- c(mean(css_est[-ind]), sd(css_est[-ind]), mse_cssest)
   ResultMat[5,] <- c(mean(css_est_trim[-ind]), sd(css_est_trim[-ind]), mse_cssest_trim)
   ResultMat[6,] <- c(mean(cram_est[-ind]), sd(cram_est[-ind]),mse_cram)
   ResultMat[7,] <- c(mean(cram_est_trim[-ind]), sd(cram_est_trim[-ind]), mse_cram_trim)
   ResultMat[8,] <- c(mean(true_est[-ind]), sd(true_est[-ind]), mse_true)
} else{
  mean_true_est <- mean(true_est)
  mse_naive <- mean((naive_est - mean_true_est)*(naive_est - mean_true_est))
  mse_rsest <- mean((rs_est - mean_true_est)*(rs_est - mean_true_est))
  mse_cssest <- mean((css_est - mean_true_est)*(css_est - mean_true_est))
  mse_cssest_trim <- mean((css_est_trim - mean_true_est)*(css_est_trim - mean_true_est))
  mse_rs_trim <- mean((rs_est_trim - mean_true_est)*(rs_est_trim - mean_true_est))
  mse_cram <- mean((cram_est - mean_true_est)*(cram_est - mean_true_est))
  mse_cram_trim <- mean((cram_est_trim - mean_true_est)*(cram_est_trim - mean_true_est))
  mse_true <- mean((true_est - mean_true_est)*(true_est - mean_true_est))

  ResultMat[1,] <- c(mean(naive_est), sd(naive_est), mse_naive)
  ResultMat[2,] <- c(mean(rs_est), sd(rs_est), mse_rsest)
  ResultMat[3,] <- c(mean(rs_est_trim), sd(rs_est_trim), mse_rs_trim)
  ResultMat[4,] <- c(mean(css_est), sd(css_est), mse_cssest)
  ResultMat[5,] <- c(mean(css_est_trim), sd(css_est_trim), mse_cssest_trim)
  ResultMat[6,] <- c(mean(cram_est), sd(cram_est),mse_cram)
  ResultMat[7,] <- c(mean(cram_est_trim), sd(cram_est_trim), mse_cram_trim)
  ResultMat[8,] <- c(mean(true_est), sd(true_est), mse_true)
}

round(ResultMat, 3)

cap=c("Results from simulation study 1. Treatment effect estimates for the target sample and their standard deviations for
different approaches (1,000 simulations).")
xtable(ResultMat, digits=3, caption=cap)








