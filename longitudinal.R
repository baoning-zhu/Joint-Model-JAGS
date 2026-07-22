install.packages(c("lme4","dplyr"))

########################################################
# Longitudinal Mixed Effects Model Simulation Study
# Author: PhD Simulation Study
########################################################

rm(list = ls())

set.seed(123)

library(lme4)
library(dplyr)

########################################################
# 1. True Parameters
########################################################

true_par <- c(
  
  beta0 = 2.0,
  
  beta1 = 0.5,
  
  sd_b0 = 1.0,
  
  sd_b1 = 0.4,
  
  sigma = 1.0
  
)

########################################################
# 2. Simulate Longitudinal Dataset
########################################################

simulate_lme_data <- function(
    
  N = 100,
  
  n_time = 5,
  
  beta0 = 2,
  
  beta1 = 0.5,
  
  sd_b0 = 1,
  
  sd_b1 = 0.4,
  
  sigma = 1
  
){
  
  id <- rep(
    1:N,
    each = n_time
  )
  
  time <- rep(
    0:(n_time - 1),
    times = N
  )
  
  b0 <- rnorm(
    N,
    mean = 0,
    sd = sd_b0
  )
  
  b1 <- rnorm(
    N,
    mean = 0,
    sd = sd_b1
  )
  
  y <- numeric(
    N * n_time
  )
  
  for(i in 1:N){
    
    idx <- which(id == i)
    
    y[idx] <-
      
      beta0 +
      
      beta1 * time[idx] +
      
      b0[i] +
      
      b1[i] * time[idx] +
      
      rnorm(
        n_time,
        mean = 0,
        sd = sigma
      )
  }
  
  data.frame(
    
    id = factor(id),
    
    time = time,
    
    y = y
    
  )
}

########################################################
# 3. Fit Mixed Effects Model
########################################################

fit_one_lme <- function(dat){
  
  fit <- lmer(
    
    y ~ time +
      
      (time | id),
    
    data = dat,
    
    REML = FALSE
    
  )
  
  fit
}

########################################################
# 4. Extract Estimates and Confidence Intervals
########################################################

extract_results <- function(fit){
  
  est <- fixef(fit)
  
  se <- sqrt(
    diag(vcov(fit))
  )
  
  ci <- confint(
    
    fit,
    
    parm = c(
      "(Intercept)",
      "time"
    ),
    
    method = "Wald"
    
  )
  
  data.frame(
    
    Parameter = c(
      "beta0",
      "beta1"
    ),
    
    Estimate = c(
      est[1],
      est[2]
    ),
    
    SE = c(
      se[1],
      se[2]
    ),
    
    Lower = c(
      ci[1,1],
      ci[2,1]
    ),
    
    Upper = c(
      ci[1,2],
      ci[2,2]
    )
    
  )
}

########################################################
# 5. One Monte Carlo Replicate
########################################################

one_replicate <- function(
    
  rep_id,
  
  N,
  
  true_par
  
){
  
  dat <- simulate_lme_data(
    
    N = N,
    
    beta0 = true_par["beta0"],
    
    beta1 = true_par["beta1"],
    
    sd_b0 = true_par["sd_b0"],
    
    sd_b1 = true_par["sd_b1"],
    
    sigma = true_par["sigma"]
    
  )
  
  fit <- fit_one_lme(dat)
  
  res <- extract_results(fit)
  
  res$True <- c(
    
    true_par["beta0"],
    
    true_par["beta1"]
    
  )
  
  res$Cover <-
    
    as.integer(
      
      res$Lower <= res$True &
        
        res$Upper >= res$True
      
    )
  
  res$Replicate <- rep_id
  
  res
}

########################################################
# 6. Monte Carlo Simulation
########################################################

run_mc <- function(
    
  n_rep = 100,
  
  N = 100,
  
  true_par
  
){
  
  results <- vector(
    
    "list",
    
    n_rep
    
  )
  
  for(r in 1:n_rep){
    
    cat(
      "Replication",
      r,
      "of",
      n_rep,
      "\n"
    )
    
    results[[r]] <-
      
      tryCatch(
        
        one_replicate(
          
          rep_id = r,
          
          N = N,
          
          true_par = true_par
          
        ),
        
        error = function(e){
          
          cat(
            "Failed:",
            r,
            "\n"
          )
          
          NULL
        }
        
      )
  }
  
  results <- results[
    !sapply(results,is.null)
  ]
  
  do.call(
    rbind,
    results
  )
}

########################################################
# 7. Monte Carlo Summary
########################################################

summarize_mc <- function(res){
  
  by_param <- split(
    
    res,
    
    res$Parameter
    
  )
  
  out <- lapply(
    
    by_param,
    
    function(df){
      
      true_val <- unique(df$True)
      
      mean_est <- mean(df$Estimate)
      
      mean_sd <- mean(df$SE)
      
      bias <- mean_est - true_val
      
      rbias <-
        
        100 *
        
        bias /
        
        true_val
      
      mse <- mean(
        
        (df$Estimate -
           true_val)^2
        
      )
      
      rmse <- sqrt(mse)
      
      coverage <- mean(df$Cover)
      
      mcse <-
        
        sd(df$Estimate) /
        
        sqrt(nrow(df))
      
      data.frame(
        
        Parameter =
          unique(df$Parameter),
        
        True_Value =
          round(true_val,4),
        
        Mean_Estimate =
          round(mean_est,4),
        
        Mean_SD =
          round(mean_sd,4),
        
        Bias =
          round(bias,4),
        
        Relative_Bias =
          round(rbias,2),
        
        MSE =
          round(mse,4),
        
        RMSE =
          round(rmse,4),
        
        Coverage =
          round(coverage,4),
        
        MCSE =
          round(mcse,5)
        
      )
    }
    
  )
  
  do.call(
    rbind,
    out
  )
}

########################################################
# 8. Test Run
########################################################

mc_test <- run_mc(
  
  n_rep = 20,
  
  N = 100,
  
  true_par = true_par
  
)

test_summary <- summarize_mc(
  
  mc_test
  
)

print(test_summary)

write.csv(
  
  mc_test,
  
  "LME_Test_Detail.csv",
  
  row.names = FALSE
  
)

write.csv(
  
  test_summary,
  
  "LME_Test_Summary.csv",
  
  row.names = FALSE
  
)

########################################################
# 9. Official Simulation Study
########################################################

sample_sizes <- c(
  
  100,
  
  300,
  
  500
  
)

all_tables <- list()

for(n in sample_sizes){
  
  cat("\n")
  
  cat("========================\n")
  
  cat("Running Sample Size =",n,"\n")
  
  cat("========================\n")
  
  mc <- run_mc(
    
    n_rep = 100,
    
    N = n,
    
    true_par = true_par
    
  )
  
  tab <- summarize_mc(mc)
  
  tab$Sample_Size <- n
  
  all_tables[[as.character(n)]] <- tab
  
  write.csv(
    
    tab,
    
    paste0(
      
      "LME_Summary_N",
      
      n,
      
      ".csv"
      
    ),
    
    row.names = FALSE
    
  )
  
  print(tab)
}

########################################################
# 10. Final Results
########################################################

final_results <- do.call(
  
  rbind,
  
  all_tables
  
)

write.csv(
  
  final_results,
  
  "LME_Simulation_All_Results.csv",
  
  row.names = FALSE
  
)

cat("\n")

cat("Simulation Study Finished\n")

print(final_results)