############################################################
## Bayesian Longitudinal Mixed Effects Model Simulation
## JAGS Version
############################################################

rm(list=ls())

library(rjags)
library(coda)
library(MASS)
library(dplyr)

set.seed(12345)

############################################################
## True Parameters
############################################################

true.par <- c(
  
  beta0 = 2,
  
  beta1 = 0.5,
  
  sigma_b = 1,
  
  sigma = 1
  
)

############################################################
## Simulation Setting
############################################################

n.subject <- 100

n.time <- 5

times <- 0:(n.time-1)

nsim <- 100

############################################################
## Storage
############################################################

result <- matrix(NA,
                 nsim,
                 8)

colnames(result) <- c(
  
  "beta0",
  
  "beta1",
  
  "sigma_b",
  
  "sigma",
  
  "sd_beta0",
  
  "sd_beta1",
  
  "sd_sigma_b",
  
  "sd_sigma"
  
)

coverage <- matrix(0,
                   nsim,
                   4)

############################################################
## JAGS model
############################################################

modelString <- "

model{

for(i in 1:N){

    b[i] ~ dnorm(0,tau_b)

    for(j in 1:nt){

        mu[i,j] <- beta0 + beta1*time[j] + b[i]

        y[i,j] ~ dnorm(mu[i,j],tau)

    }

}

#############################

beta0 ~ dnorm(0,0.0001)

beta1 ~ dnorm(0,0.0001)

sigma ~ dunif(0,20)

sigma_b ~ dunif(0,20)

tau <- pow(sigma,-2)

tau_b <- pow(sigma_b,-2)

}

"

writeLines(modelString,"lme_jags.txt")

############################################################
## Simulation
############################################################

for(sim in 1:nsim){
  
  cat("Simulation:",sim,"\n")
  
  #############################
  ## Generate data
  #############################
  
  b <- rnorm(n.subject,
             0,
             true.par["sigma_b"])
  
  Y <- matrix(NA,
              n.subject,
              n.time)
  
  for(i in 1:n.subject){
    
    for(j in 1:n.time){
      
      Y[i,j] <-
        
        true.par["beta0"] +
        
        true.par["beta1"]*times[j] +
        
        b[i] +
        
        rnorm(1,0,true.par["sigma"])
      
    }
    
  }
  
  #############################
  ## Data for JAGS
  #############################
  
  data.jags <- list(
    
    y = Y,
    
    N = n.subject,
    
    nt = n.time,
    
    time = times
    
  )
  
  #############################
  ## Initial values
  #############################
  
  init <- function(){
    
    list(
      
      beta0=rnorm(1),
      
      beta1=rnorm(1),
      
      sigma=runif(1,0.5,2),
      
      sigma_b=runif(1,0.5,2)
      
    )
    
  }
  
  #############################
  ## Fit model
  #############################
  
  jm <- jags.model(
    
    file="lme_jags.txt",
    
    data=data.jags,
    
    inits=init,
    
    n.chains=3,
    
    quiet=TRUE
    
  )
  
  update(jm,3000)
  
  post <- coda.samples(
    
    jm,
    
    variable.names=c(
      
      "beta0",
      
      "beta1",
      
      "sigma",
      
      "sigma_b"
      
    ),
    
    n.iter=5000
    
  )
  
  sum.post <- summary(post)
  
  est <- sum.post$statistics[,"Mean"]
  
  sd.est <- sum.post$statistics[,"SD"]
  
  result[sim,] <-
    
    c(
      
      est["beta0"],
      
      est["beta1"],
      
      est["sigma_b"],
      
      est["sigma"],
      
      sd.est["beta0"],
      
      sd.est["beta1"],
      
      sd.est["sigma_b"],
      
      sd.est["sigma"]
      
    )
  
  #############################
  ## Coverage
  #############################
  
  CI <- sum.post$quantiles[,c("2.5%","97.5%")]
  
  truth <- c(
    
    true.par["beta0"],
    
    true.par["beta1"],
    
    true.par["sigma_b"],
    
    true.par["sigma"]
    
  )
  
  for(k in 1:4){
    
    coverage[sim,k] <-
      
      (CI[k,1]<=truth[k]) &
      
      (CI[k,2]>=truth[k])
    
  }
  
}

############################################################
## Performance Evaluation
############################################################

truth <- c(
  
  true.par["beta0"],
  
  true.par["beta1"],
  
  true.par["sigma_b"],
  
  true.par["sigma"]
  
)

Estimate <- colMeans(result[,1:4])

Bias <- Estimate-truth

RB <- Bias/truth*100

RMSE <- sqrt(
  
  colMeans(
    
    (result[,1:4]-matrix(truth,
                         
                         nsim,
                         
                         4,
                         
                         byrow=TRUE))^2
    
  )
  
)

MeanSD <- colMeans(result[,5:8])

Coverage <- colMeans(coverage)

############################################################
## Summary Table
############################################################

summary.table <-
  
  data.frame(
    
    Parameter=
      
      c(
        
        "beta0",
        
        "beta1",
        
        "sigma_b",
        
        "sigma"
        
      ),
    
    True=truth,
    
    Estimate=round(Estimate,3),
    
    Bias=round(Bias,4),
    
    RB=round(RB,2),
    
    RMSE=round(RMSE,4),
    
    MeanSD=round(MeanSD,4),
    
    Coverage=round(Coverage,3)
    
  )

print(summary.table)