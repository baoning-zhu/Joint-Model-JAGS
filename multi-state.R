# Clear all objects in R environment to ensure a clean running environment
rm(list = ls())
# Set random seed to 42 for reproducible results
set.seed(42)

# Load JAGS interface package for Bayesian modeling
library(rjags)
# Load package for analyzing MCMC results
library(coda)
# Load parallel computing package to improve computation speed
library(parallel)

#==================================================
# 0. Write JAGS model file
#==================================================
# Define JAGS model code as a character string
jags_code <- "
model{

  # Baseline hazard parameters: log-normal prior (ensures positivity)
  lambda12 ~ dlnorm(0, 1)  # Baseline hazard for state 1→2
  lambda13 ~ dlnorm(0, 1)  # Baseline hazard for state 1→3
  lambda23 ~ dlnorm(0, 1)  # Baseline hazard for state 2→3

  # Covariate coefficients: normal prior (mean 0, variance 4)
  beta12 ~ dnorm(0, 1.0/4)  # Covariate coefficient for 1→2
  beta13 ~ dnorm(0, 1.0/4)  # Covariate coefficient for 1→3
  beta23 ~ dnorm(0, 1.0/4)  # Covariate coefficient for 2→3

  # Precision of random effects = inverse of variance
  tau_b <- pow(sigma_b, -2)
  # Standard deviation of random effects: uniform prior (0 to 5)
  sigma_b ~ dunif(0, 5)

  # Loop over individual random effects: one per subject
  for(i in 1:N){
    b[i] ~ dnorm(0, tau_b)  # Individual random effect, normal mean 0, precision tau_b
  }

  # Loop over transition records: each row is one transition interval
  for(n in 1:n_trans){

    dt[n] <- end_time[n] - start_time[n]  # Time interval = end - start
    zi[n] <- z[id[n]]                     # Covariate for current individual
    bi[n] <- b[id[n]]                     # Random effect for current individual

    # Hazard functions in Cox model form
    h12[n] <- lambda12 * exp(beta12 * zi[n] + bi[n])  # Hazard 1→2
    h13[n] <- lambda13 * exp(beta13 * zi[n] + bi[n])  # Hazard 1→3
    h23[n] <- lambda23 * exp(beta23 * zi[n] + bi[n])  # Hazard 2→3

    # Indicators for starting state
    is_from1[n] <- equals(from[n], 1)  # 1 if from state 1, else 0
    is_from2[n] <- equals(from[n], 2)  # 1 if from state 2, else 0

    # Indicators for ending state
    is_to0[n] <- equals(to[n], 0)    # 1 if censored, else 0
    is_to2[n] <- equals(to[n], 2)    # 1 if to state 2, else 0
    is_to3[n] <- equals(to[n], 3)    # 1 if to state 3, else 0

    # Log-likelihood when starting from state 1
    loglik1[n] <-
      is_from1[n] * (
        is_to2[n] * (log(h12[n]) - (h12[n] + h13[n]) * dt[n]) +  # Likelihood for 1→2
        is_to3[n] * (log(h13[n]) - (h12[n] + h13[n]) * dt[n]) +  # Likelihood for 1→3
        is_to0[n] * (-(h12[n] + h13[n]) * dt[n])                 # Likelihood for censoring from 1
      )

    # Log-likelihood when starting from state 2
    loglik2[n] <-
      is_from2[n] * (
        is_to3[n] * (log(h23[n]) - h23[n] * dt[n]) +  # Likelihood for 2→3
        is_to0[n] * (-h23[n] * dt[n])                 # Likelihood for censoring from 2
      )

    loglik[n] <- loglik1[n] + loglik2[n]  # Total log-likelihood
    phi[n] <- -loglik[n] + C              # Poisson parameter transformation (zero-trick)
    zeros[n] ~ dpois(phi[n])              # Zero count for constructing likelihood in JAGS
  }
}
"
# Write JAGS model code to a local file for later use
writeLines(jags_code, "multi_state_model_cov.jags")

#==================================================
# 1. Define true parameter values (for data simulation)
#==================================================
true_par <- c(
  lambda12 = 1,   # True: baseline hazard 1→2
  lambda13 = 1,   # True: baseline hazard 1→3
  lambda23 = 1,   # True: baseline hazard 2→3
  beta12   = -1,  # True: covariate coefficient 1→2
  beta13   =  2,  # True: covariate coefficient 1→3
  beta23   =  1,  # True: covariate coefficient 2→3
  sigma_b  =  0.5 # True: SD of individual random effects
)

#==================================================
# 2. Simulate data for a single multi-state model
#==================================================
# Define data simulation function
simulate_multistate_data <- function(
    N = 1000,            # Sample size: default 80 individuals
    lambda12 = 1,      # Baseline hazard 1→2
    lambda13 = 1,      # Baseline hazard 1→3
    lambda23 = 1,      # Baseline hazard 2→3
    beta12   = -1,     # Coefficient 1→2
    beta13   = 2,      # Coefficient 1→3
    beta23   = 1,      # Coefficient 2→3
    sigma_b  = 0.5,    # SD of random effects
    censor_rate = 0.08,# Censoring rate: default 0.15
    Tmax = 15          # Maximum observation time: default 10
){
  z <- rnorm(N)                          # Individual covariate (standard normal)
  b <- rnorm(N, 0, sigma_b)              # Individual random effects
  censor_time <- rexp(N, rate = censor_rate)  # Censoring time (exponential)
  
  multi_state_list <- list()             # Initialize list to store transitions
  row_id <- 1                            # Initialize row counter
  
  # Simulate state transitions for each individual
  for (i in 1:N) {
    t_current <- 0          # Current time: start at 0
    state_current <- 1      # Current state: start at 1
    obs_end <- min(censor_time[i], Tmax)  # Actual observation end time
    
    # Continue simulating until observation ends or absorbing state is reached
    while (t_current < obs_end && state_current != 0) {
      
      if (state_current == 1) {  # Current state is 1
        # Hazards for 1→2 and 1→3
        h12 <- lambda12 * exp(beta12 * z[i] + b[i])
        h13 <- lambda13 * exp(beta13 * z[i] + b[i])
        rate <- h12 + h13  # Total hazard
        
        # Simulate next transition time (exponential)
        t_next <- t_current - log(runif(1)) / rate
        
        # If next transition exceeds observation time: censor
        if (t_next > obs_end) {
          multi_state_list[[row_id]] <- data.frame(
            id = i,          # Individual ID
            from = 1,        # Starting state
            to = 0,          # Ending state: censored
            start = t_current,  # Start time
            end = obs_end,      # End time
            status = 0       # Status: censored
          )
          row_id <- row_id + 1
          break
        }
        
        # Probabilities of each transition
        probs <- c(h12, h13) / (h12 + h13)
        next_state <- sample(c(2, 3), size = 1, prob = probs)
        
      } else if (state_current == 2) {  # Current state is 2
        # Hazard for 2→3
        h23 <- lambda23 * exp(beta23 * z[i] + b[i])
        rate <- h23
        
        # Simulate next transition time
        t_next <- t_current - log(runif(1)) / rate
        
        # If next transition exceeds observation time: censor
        if (t_next > obs_end) {
          multi_state_list[[row_id]] <- data.frame(
            id = i,
            from = 2,
            to = 0,
            start = t_current,
            end = obs_end,
            status = 0
          )
          row_id <- row_id + 1
          break
        }
        
        next_state <- 3  # Only possible transition from 2 is to 3
      } else {  # Absorbing state: stop
        break
      }
      
      # Record valid transition
      multi_state_list[[row_id]] <- data.frame(
        id = i,
        from = state_current,
        to = next_state,
        start = t_current,
        end = t_next,
        status = 1  # Status: event occurred
      )
      row_id <- row_id + 1
      
      # Update current time and state
      t_current <- t_next
      state_current <- next_state
    }
  }
  
  # Combine list into data frame
  multi_state_data <- do.call(rbind, multi_state_list)
  
  # Return simulated data
  list(
    N = N,
    z = z,
    multi_state_data = multi_state_data
  )
}

#==================================================
# 3. Fit JAGS model to one simulated dataset
#==================================================
# Define model fitting function
fit_one_jags <- function(simdat,
                         n_adapt = 100,   # Adaptive iterations
                         n_burn  = 300,   # Burn-in iterations
                         n_iter  = 800,   # Posterior sampling iterations
                         C_const = 10000){ # Constant for likelihood transformation
  
  multi_state_data <- simdat$multi_state_data
  N <- simdat$N
  z <- simdat$z
  
  # Prepare data list for JAGS
  jags_data <- list(
    N = N,
    n_trans = nrow(multi_state_data),
    from = multi_state_data$from,
    to = multi_state_data$to,
    id = multi_state_data$id,
    start_time = multi_state_data$start,
    end_time = multi_state_data$end,
    z = z,
    zeros = rep(0, nrow(multi_state_data)),
    C = C_const
  )
  
  # Initial values
  inits1 <- list(
    lambda12 = 0.8,
    lambda13 = 0.8,
    lambda23 = 0.8,
    beta12 = 0,
    beta13 = 0,
    beta23 = 0,
    sigma_b = 1,
    b = rnorm(N, 0, 0.1)
  )
  
  # Parameters to monitor
  params <- c("lambda12","lambda13","lambda23",
              "beta12","beta13","beta23","sigma_b")
  
  # Compile and initialize JAGS model
  jm <- jags.model(
    file = "multi_state_model_cov.jags",
    data = jags_data,
    inits = inits1,
    n.chains = 3,
    n.adapt = n_adapt
  )
  
  # Burn-in phase
  update(jm, n_burn)
  
  # Posterior sampling
  samp <- coda.samples(
    jm,
    variable.names = params,
    n.iter = n_iter
  )
  
  return(samp)
}

#==================================================
# 4. Extract posterior summary statistics (credible intervals)
#==================================================
# Function to extract posterior intervals
extract_post_ci <- function(samp, probs = c(0.025, 0.5, 0.975)){
  mat <- as.matrix(samp)
  qs <- apply(mat, 2, quantile, probs = probs)
  data.frame(
    param = colnames(mat),
    lwr = qs[1, ],
    med = qs[2, ],
    upr = qs[3, ],
    row.names = NULL
  )
}

#==================================================
# 5. Single Monte Carlo replicate (simulate + fit + extract)
#==================================================
# Function for one replication
one_replicate <- function(rep_id,
                          N = 80,
                          true_par,
                          n_adapt = 100,
                          n_burn  = 300,
                          n_iter  = 800){
  
  
  # Simulate data
  simdat <- simulate_multistate_data(
    N = N,
    lambda12 = true_par["lambda12"],
    lambda13 = true_par["lambda13"],
    lambda23 = true_par["lambda23"],
    beta12   = true_par["beta12"],
    beta13   = true_par["beta13"],
    beta23   = true_par["beta23"],
    sigma_b  = true_par["sigma_b"]
  )
  
  # Fit JAGS model
  samp <- fit_one_jags(
    simdat,
    n_adapt = n_adapt,
    n_burn  = n_burn,
    n_iter  = n_iter
  )
  
  # Extract posterior intervals
  ci_tab <- extract_post_ci(samp)
  
  # Posterior samples matrix
  mat <- as.matrix(samp)
  
  # Posterior means
  post_mean <- colMeans(mat)
  
  # Posterior SDs
  post_sd <- apply(mat, 2, sd)
  
  # Add information
  ci_tab$true <- true_par[ci_tab$param]
  
  ci_tab$cover <- as.integer(
    ci_tab$lwr <= ci_tab$true &
      ci_tab$true <= ci_tab$upr
  )
  
  ci_tab$estimate <- post_mean[ci_tab$param]
  
  ci_tab$post_sd <- post_sd[ci_tab$param]
  
  ci_tab$replicate <- rep_id
  
  return(ci_tab)
}

#==================================================
# 6. Monte Carlo simulation (Serial Version)
#==================================================

run_mc_serial <- function(
    n_rep = 20,
    N = 100,
    true_par,
    n_adapt = 100,
    n_burn = 300,
    n_iter = 800
){
  
  res_list <- vector("list", n_rep)
  
  for(r in 1:n_rep){
    
    cat("Replication", r, "of", n_rep, "\n")
    
    res_list[[r]] <- tryCatch(
      
      one_replicate(
        rep_id = r,
        N = N,
        true_par = true_par,
        n_adapt = n_adapt,
        n_burn = n_burn,
        n_iter = n_iter
      ),
      
      error = function(e){
        
        cat(
          "Replication",
          r,
          "failed:",
          e$message,
          "\n"
        )
        
        NULL
      }
    )
  }
  
  ok <- !sapply(res_list,is.null)
  
  res <- do.call(
    rbind,
    res_list[ok]
  )
  
  list(
    detail = res,
    ok_replicates = sum(ok),
    failed_replicates = sum(!ok)
  )
}

summarize_bias_rmse <- function(mc_obj){
  
  res <- mc_obj$detail
  
  if(is.null(res)){
    stop("No successful replicates.")
  }
  
  by_param <- split(res,res$param)
  
  out <- lapply(by_param,function(df){
    
    mean_est <- mean(df$estimate)
    
    true_val <- unique(df$true)
    
    bias <- mean_est - true_val
    
    rbias <- (mean_est - true_val)/true_val
    
    mse <- mean(
      (df$estimate - true_val)^2
    )
    
    rmse <- sqrt(mse)
    
    coverage <- mean(df$cover)
    
    mean_sd <- mean(df$post_sd)
    
    mcse <- sd(df$estimate)/sqrt(nrow(df))
    
    data.frame(
      
      Parameter = unique(df$param),
      
      True_Value = round(true_val,4),
      
      Mean_Estimate = round(mean_est,4),
      
      Posterior_SD = round(mean_sd,4),
      
      Bias = round(bias,4),
      
      Relative_Bias = round(rbias,4),
      
      MSE = round(mse,4),
      
      RMSE = round(rmse,4),
      
      Coverage = round(coverage,4),
      
      MCSE = round(mcse,5)
      
    )
  })
  
  do.call(rbind,out)
}

mc_test <- run_mc_serial(
  
  n_rep = 10,
  
  N = 100,
  
  true_par = true_par,
  
  n_adapt = 100,
  
  n_burn = 300,
  
  n_iter = 800
)

cat(
  "Successful replicates:",
  mc_test$ok_replicates,
  "\n"
)

cat(
  "Failed replicates:",
  mc_test$failed_replicates,
  "\n"
)

coverage_table_test <- summarize_bias_rmse(mc_test)

print(coverage_table_test)

#==================================================
# 9. Official Simulation Study
#==================================================

sample_sizes <- c(100,300,500)

all_tables <- list()

for(n in sample_sizes){
  
  cat("\n")
  cat("=====================\n")
  cat("Running N =",n,"\n")
  cat("=====================\n")
  
  mc <- run_mc_serial(
    
    n_rep = 100,
    
    N = n,
    
    true_par = true_par,
    
    n_adapt = 100,
    
    n_burn = 300,
    
    n_iter = 800
  )
  
  cat(
    "Successful replicates:",
    mc$ok_replicates,
    "\n"
  )
  
  cat(
    "Failed replicates:",
    mc$failed_replicates,
    "\n"
  )
  
  tab <- summarize_bias_rmse(mc)
  
  tab$Sample_Size <- n
  
  all_tables[[as.character(n)]] <- tab
  
  write.csv(
    tab,
    paste0(
      "Simulation_Summary_N",
      n,
      ".csv"
    ),
    row.names = FALSE
  )
  
  print(tab)
}

final_results <- do.call(
  rbind,
  all_tables
)

write.csv(
  final_results,
  "Simulation_All_Results.csv",
  row.names = FALSE
)

cat("\n")
cat("Simulation Finished\n")

print(final_results)