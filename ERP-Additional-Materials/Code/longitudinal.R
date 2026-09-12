
############################################################
## Bayesian longitudinal simulation study
## Random-intercept and random-slope linear mixed model
############################################################
## 0. Clear workspace
############################################################

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

set.seed(2026)


############################################################
## 1. Required packages
############################################################

required_packages <- c(
  "MASS",
  "rjags",
  "coda",
  "dplyr"
)

not_installed <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(not_installed) > 0) {
  install.packages(
    not_installed,
    dependencies = TRUE
  )
}

library(MASS)
library(rjags)
library(coda)
library(dplyr)


############################################################
## 2. Output directory
############################################################

output_directory <- "longitudinal_simulation_results"

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 3. True parameter values
############################################################

true_values <- c(
  beta0_y = 2,
  beta1_y = 0.5,
  beta2_y = 1,
  sd_b0   = 1,
  sd_b1   = 0.6,
  rho     = 0.3,
  sigma_y = 1
)

sample_sizes <- c(
  500
)

n_time <- 8

R_simulation <- 100


############################################################
## 4. MCMC settings
##
## These settings improve convergence without changing
## the statistical model or prior distributions.
############################################################

n_chains <- 3
n_adapt <- 2000
n_burnin <- 5000
n_iter <- 10000
thin <- 1


############################################################
## 5. Data-generating function
############################################################

simulate_longitudinal_data <- function(
    N,
    n_time = 8,
    beta0_y = 2,
    beta1_y = 0.5,
    beta2_y = 1,
    sd_b0 = 1,
    sd_b1 = 0.6,
    rho = 0.3,
    sigma_y = 1
) {
  
  ##########################################################
  ## Input checks
  ##########################################################
  
  if (N <= 0) {
    stop("N must be a positive integer.")
  }
  
  if (n_time <= 1) {
    stop("n_time must be greater than 1.")
  }
  
  if (sd_b0 <= 0 || sd_b1 <= 0 || sigma_y <= 0) {
    stop("All standard deviations must be positive.")
  }
  
  if (abs(rho) >= 1) {
    stop("rho must lie strictly between -1 and 1.")
  }
  
  
  ##########################################################
  ## Baseline covariate
  ##
  ## x_i ~ N(0, 1)
  ##########################################################
  
  x <- rnorm(
    n = N,
    mean = 0,
    sd = 1
  )
  
  
  ##########################################################
  ## True covariance matrix of random effects
  ##########################################################
  
  D_true <- matrix(
    c(
      sd_b0^2,
      rho * sd_b0 * sd_b1,
      rho * sd_b0 * sd_b1,
      sd_b1^2
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
  
  
  ##########################################################
  ## Subject-specific random effects
  ##
  ## (b0_i, b1_i)' ~ N_2(0, D)
  ##########################################################
  
  random_effects <- MASS::mvrnorm(
    n = N,
    mu = c(0, 0),
    Sigma = D_true
  )
  
  b0 <- random_effects[, 1]
  b1 <- random_effects[, 2]
  
  
  ##########################################################
  ## Irregular observation times
  ##
  ## Eight times independently generated from Uniform(0, 8)
  ## and sorted within subject
  ##########################################################
  
  time_matrix <- matrix(
    NA_real_,
    nrow = N,
    ncol = n_time
  )
  
  for (i in seq_len(N)) {
    
    time_matrix[i, ] <- sort(
      runif(
        n = n_time,
        min = 0,
        max = 8
      )
    )
  }
  
  
  ##########################################################
  ## Longitudinal outcome generation
  ##
  ## m_i(s_ij)
  ## = beta0_y
  ## + beta1_y s_ij
  ## + beta2_y x_i
  ## + b0_i
  ## + b1_i s_ij
  ##
  ## Y_ij ~ N(m_i(s_ij), sigma_y^2)
  ##########################################################
  
  mu_matrix <- matrix(
    NA_real_,
    nrow = N,
    ncol = n_time
  )
  
  y_matrix <- matrix(
    NA_real_,
    nrow = N,
    ncol = n_time
  )
  
  for (i in seq_len(N)) {
    
    for (j in seq_len(n_time)) {
      
      mu_matrix[i, j] <-
        beta0_y +
        beta1_y * time_matrix[i, j] +
        beta2_y * x[i] +
        b0[i] +
        b1[i] * time_matrix[i, j]
      
      y_matrix[i, j] <- rnorm(
        n = 1,
        mean = mu_matrix[i, j],
        sd = sigma_y
      )
    }
  }
  
  
  ##########################################################
  ## Long-format data for checking
  ##########################################################
  
  long_data <- data.frame(
    id = rep(
      seq_len(N),
      each = n_time
    ),
    
    occasion = rep(
      seq_len(n_time),
      times = N
    ),
    
    time = as.vector(
      t(time_matrix)
    ),
    
    x = rep(
      x,
      each = n_time
    ),
    
    y = as.vector(
      t(y_matrix)
    ),
    
    true_mu = as.vector(
      t(mu_matrix)
    ),
    
    true_b0 = rep(
      b0,
      each = n_time
    ),
    
    true_b1 = rep(
      b1,
      each = n_time
    )
  )
  
  
  ##########################################################
  ## Return generated objects
  ##########################################################
  
  return(
    list(
      N = N,
      n_time = n_time,
      x = x,
      time = time_matrix,
      y = y_matrix,
      mu = mu_matrix,
      b0 = b0,
      b1 = b1,
      D_true = D_true,
      long_data = long_data
    )
  )
}


############################################################
## 6. JAGS model
############################################################

jags_model_string <- "
model {

  ##########################################################
  ## Longitudinal likelihood
  ##########################################################

  for (i in 1:N) {

    for (j in 1:n_time) {

      y[i, j] ~ dnorm(mu[i, j], tau_y)

      mu[i, j] <-
        beta0_y +
        beta1_y * time[i, j] +
        beta2_y * x[i] +
        b[i, 1] +
        b[i, 2] * time[i, j]
    }


    ########################################################
    ## Random intercept and random slope
    ########################################################

    b[i, 1:2] ~ dmnorm(
      zero[1:2],
      Tau_b[1:2, 1:2]
    )
  }


  ##########################################################
  ## Fixed-effect priors
  ##
  ## beta_k ~ N(0, 10^4)
  ## JAGS uses precision = 1 / variance = 0.0001
  ##########################################################

  beta0_y ~ dnorm(0, 0.0001)
  beta1_y ~ dnorm(0, 0.0001)
  beta2_y ~ dnorm(0, 0.0001)


  ##########################################################
  ## Random-effects covariance prior
  ##
  ## D ~ IW(nu0, S0)
  ## implemented through:
  ## Tau_b = D^{-1} ~ Wishart(S0^{-1}, nu0)
  ##########################################################

  Tau_b[1:2, 1:2] ~ dwish(
    R_b[1:2, 1:2],
    nu0
  )

  D[1:2, 1:2] <- inverse(
    Tau_b[1:2, 1:2]
  )


  ##########################################################
  ## Derived random-effects parameters
  ##########################################################

  var_b0 <- D[1, 1]
  var_b1 <- D[2, 2]
  cov_b01 <- D[1, 2]

  sd_b0 <- sqrt(var_b0)
  sd_b1 <- sqrt(var_b1)

  rho <- cov_b01 /
    sqrt(var_b0 * var_b1)


  ##########################################################
  ## Residual variance prior
  ##
  ## sigma_y^2 ~ IG(0.001, 0.001)
  ## equivalent to:
  ## tau_y ~ Gamma(0.001, 0.001)
  ##########################################################

  tau_y ~ dgamma(
    0.001,
    0.001
  )

  sigma2_y <- 1 / tau_y
  sigma_y <- sqrt(sigma2_y)
}
"


############################################################
## 7. Fit one simulated dataset
############################################################

fit_one_dataset <- function(
    simulated_data,
    model_string = jags_model_string,
    n_chains = 3,
    n_adapt = 2000,
    n_burnin = 5000,
    n_iter = 10000,
    thin = 1,
    seed = 2026
) {
  
  N <- simulated_data$N
  n_time <- simulated_data$n_time
  
  
  ##########################################################
  ## Prior settings
  ##########################################################
  
  S0 <- matrix(
    c(
      1, 0,
      0, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  R_b <- solve(S0)
  nu0 <- 4
  
  
  ##########################################################
  ## JAGS data
  ##########################################################
  
  jags_data <- list(
    N = N,
    n_time = n_time,
    y = simulated_data$y,
    time = simulated_data$time,
    x = simulated_data$x,
    zero = c(0, 0),
    R_b = R_b,
    nu0 = nu0
  )
  
  
  ##########################################################
  ## Stable but dispersed initial values
  ##
  ## These do not change the model or true parameter values.
  ##########################################################
  
  initial_values <- list(
    
    list(
      beta0_y = 1.0,
      beta1_y = 0.2,
      beta2_y = 0.5,
      tau_y = 0.8,
      Tau_b = matrix(
        c(
          0.8, 0,
          0, 2.0
        ),
        nrow = 2,
        byrow = TRUE
      ),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = seed + 1000
    ),
    
    list(
      beta0_y = 2.0,
      beta1_y = 0.5,
      beta2_y = 1.0,
      tau_y = 1.0,
      Tau_b = matrix(
        c(
          1.0, 0,
          0, 2.5
        ),
        nrow = 2,
        byrow = TRUE
      ),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = seed + 2000
    ),
    
    list(
      beta0_y = 3.0,
      beta1_y = 0.8,
      beta2_y = 1.5,
      tau_y = 1.2,
      Tau_b = matrix(
        c(
          1.2, 0,
          0, 3.0
        ),
        nrow = 2,
        byrow = TRUE
      ),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = seed + 3000
    )
  )
  
  if (n_chains != 3) {
    stop(
      "This code currently defines initial values for exactly 3 chains."
    )
  }
  
  
  ##########################################################
  ## Compile model
  ##########################################################
  
  model_connection <- textConnection(
    model_string
  )
  
  on.exit(
    close(model_connection),
    add = TRUE
  )
  
  jags_fit <- jags.model(
    file = model_connection,
    data = jags_data,
    inits = initial_values,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  
  
  ##########################################################
  ## Burn-in
  ##########################################################
  
  update(
    object = jags_fit,
    n.iter = n_burnin,
    progress.bar = "none"
  )
  
  
  ##########################################################
  ## Monitored parameters
  ##########################################################
  
  monitored_parameters <- c(
    "beta0_y",
    "beta1_y",
    "beta2_y",
    "sd_b0",
    "sd_b1",
    "rho",
    "sigma_y"
  )
  
  
  ##########################################################
  ## Posterior sampling
  ##########################################################
  
  posterior_samples <- coda.samples(
    model = jags_fit,
    variable.names = monitored_parameters,
    n.iter = n_iter,
    thin = thin,
    progress.bar = "none"
  )
  
  return(posterior_samples)
}


############################################################
## 8. Extract posterior estimates and diagnostics
############################################################

extract_posterior_results <- function(
    posterior_samples,
    replication,
    sample_size
) {
  
  combined_samples <- as.matrix(
    posterior_samples
  )
  
  target_parameters <- c(
    "beta0_y",
    "beta1_y",
    "beta2_y",
    "sd_b0",
    "sd_b1",
    "rho",
    "sigma_y"
  )
  
  
  ##########################################################
  ## Check parameter availability
  ##########################################################
  
  missing_parameters <- setdiff(
    target_parameters,
    colnames(combined_samples)
  )
  
  if (length(missing_parameters) > 0) {
    
    stop(
      paste(
        "Missing parameters:",
        paste(
          missing_parameters,
          collapse = ", "
        )
      )
    )
  }
  
  
  ##########################################################
  ## Gelman-Rubin diagnostics
  ##########################################################
  
  gelman_output <- tryCatch(
    
    gelman.diag(
      posterior_samples[, target_parameters],
      multivariate = FALSE,
      autoburnin = FALSE
    )$psrf,
    
    error = function(e) {
      
      matrix(
        NA_real_,
        nrow = length(target_parameters),
        ncol = 2,
        dimnames = list(
          target_parameters,
          c(
            "Point est.",
            "Upper C.I."
          )
        )
      )
    }
  )
  
  
  ##########################################################
  ## Effective sample size
  ##########################################################
  
  ess_output <- tryCatch(
    
    effectiveSize(
      posterior_samples[, target_parameters]
    ),
    
    error = function(e) {
      
      setNames(
        rep(
          NA_real_,
          length(target_parameters)
        ),
        target_parameters
      )
    }
  )
  
  
  ##########################################################
  ## Extract posterior summaries
  ##########################################################
  
  result_list <- lapply(
    target_parameters,
    function(parameter_name) {
      
      parameter_samples <-
        combined_samples[, parameter_name]
      
      data.frame(
        Replication = replication,
        Sample_Size = sample_size,
        Parameter = parameter_name,
        
        Estimate = mean(
          parameter_samples
        ),
        
        Posterior_SD = sd(
          parameter_samples
        ),
        
        Lower_95 = unname(
          quantile(
            parameter_samples,
            probs = 0.025
          )
        ),
        
        Upper_95 = unname(
          quantile(
            parameter_samples,
            probs = 0.975
          )
        ),
        
        Rhat = unname(
          gelman_output[
            parameter_name,
            "Point est."
          ]
        ),
        
        Rhat_Upper = unname(
          gelman_output[
            parameter_name,
            "Upper C.I."
          ]
        ),
        
        ESS = unname(
          ess_output[
            parameter_name
          ]
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  bind_rows(
    result_list
  )
}


############################################################
## 9. Run repeated simulation for one sample size
############################################################

run_simulation_setting <- function(
    N,
    R = 100,
    n_time = 8,
    n_chains = 3,
    n_adapt = 2000,
    n_burnin = 5000,
    n_iter = 10000,
    thin = 1,
    base_seed = 2026,
    save_directory = NULL
) {
  
  all_results <- vector(
    mode = "list",
    length = R
  )
  
  if (!is.null(save_directory)) {
    
    dir.create(
      save_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  for (r in seq_len(R)) {
    
    cat(
      "\nSample size:",
      N,
      "| Replication:",
      r,
      "of",
      R,
      "\n"
    )
    
    replication_seed <-
      base_seed +
      N * 10000 +
      r * 100
    
    set.seed(
      replication_seed
    )
    
    
    ########################################################
    ## Generate one dataset
    ########################################################
    
    simulated_data <- simulate_longitudinal_data(
      N = N,
      n_time = n_time,
      beta0_y = true_values["beta0_y"],
      beta1_y = true_values["beta1_y"],
      beta2_y = true_values["beta2_y"],
      sd_b0 = true_values["sd_b0"],
      sd_b1 = true_values["sd_b1"],
      rho = true_values["rho"],
      sigma_y = true_values["sigma_y"]
    )
    
    
    ########################################################
    ## Fit model
    ########################################################
    
    fit_result <- tryCatch(
      
      {
        
        posterior_samples <- fit_one_dataset(
          simulated_data = simulated_data,
          n_chains = n_chains,
          n_adapt = n_adapt,
          n_burnin = n_burnin,
          n_iter = n_iter,
          thin = thin,
          seed = replication_seed
        )
        
        extract_posterior_results(
          posterior_samples = posterior_samples,
          replication = r,
          sample_size = N
        )
      },
      
      error = function(e) {
        
        warning(
          paste(
            "Sample size",
            N,
            "replication",
            r,
            "failed:",
            conditionMessage(e)
          )
        )
        
        data.frame(
          Replication = r,
          Sample_Size = N,
          Parameter = names(true_values),
          Estimate = NA_real_,
          Posterior_SD = NA_real_,
          Lower_95 = NA_real_,
          Upper_95 = NA_real_,
          Rhat = NA_real_,
          Rhat_Upper = NA_real_,
          ESS = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    )
    
    all_results[[r]] <- fit_result
    
    
    ########################################################
    ## Save intermediate results after each replication
    ########################################################
    
    if (!is.null(save_directory)) {
      
      current_results <- bind_rows(
        all_results[
          seq_len(r)
        ]
      )
      
      output_file <- file.path(
        save_directory,
        paste0(
          "longitudinal_raw_results_N",
          N,
          ".csv"
        )
      )
      
      write.csv(
        current_results,
        file = output_file,
        row.names = FALSE
      )
    }
    
    
    ########################################################
    ## Clean memory
    ########################################################
    
    rm(
      simulated_data,
      fit_result
    )
    
    gc()
  }
  
  bind_rows(
    all_results
  )
}


############################################################
## 10. Summarise simulation performance
############################################################

summarise_simulation_results <- function(
    raw_results,
    true_parameter_values = true_values
) {
  
  true_value_table <- data.frame(
    Parameter = names(
      true_parameter_values
    ),
    
    True_Value = as.numeric(
      true_parameter_values
    ),
    
    stringsAsFactors = FALSE
  )
  
  results_with_truth <- raw_results %>%
    
    left_join(
      true_value_table,
      by = "Parameter"
    ) %>%
    
    mutate(
      
      Covered = case_when(
        
        is.na(Lower_95) |
          is.na(Upper_95) ~ NA_real_,
        
        Lower_95 <= True_Value &
          Upper_95 >= True_Value ~ 1,
        
        TRUE ~ 0
      ),
      
      Converged = case_when(
        
        is.na(Rhat) ~ NA_real_,
        
        Rhat < 1.05 ~ 1,
        
        TRUE ~ 0
      )
    )
  
  summary_table <- results_with_truth %>%
    
    group_by(
      Sample_Size,
      Parameter,
      True_Value
    ) %>%
    
    summarise(
      
      Number_Successful = sum(
        !is.na(Estimate)
      ),
      
      Mean_Estimate = mean(
        Estimate,
        na.rm = TRUE
      ),
      
      Empirical_SD = sd(
        Estimate,
        na.rm = TRUE
      ),
      
      Mean_Posterior_SD = mean(
        Posterior_SD,
        na.rm = TRUE
      ),
      
      Bias = mean(
        Estimate,
        na.rm = TRUE
      ) - first(True_Value),
      
      MSE = mean(
        (Estimate - True_Value)^2,
        na.rm = TRUE
      ),
      
      RMSE = sqrt(
        mean(
          (Estimate - True_Value)^2,
          na.rm = TRUE
        )
      ),
      
      Coverage_95 = mean(
        Covered,
        na.rm = TRUE
      ),
      
      MCSE_Mean_Estimate = if (
        sum(!is.na(Estimate)) > 1
      ) {
        
        sd(
          Estimate,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              !is.na(Estimate)
            )
          )
        
      } else {
        
        NA_real_
      },
      
      Mean_Rhat = mean(
        Rhat,
        na.rm = TRUE
      ),
      
      Maximum_Rhat = if (
        all(
          is.na(Rhat)
        )
      ) {
        
        NA_real_
        
      } else {
        
        max(
          Rhat,
          na.rm = TRUE
        )
      },
      
      Mean_ESS = mean(
        ESS,
        na.rm = TRUE
      ),
      
      Convergence_Rate = mean(
        Converged,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      Relative_Bias = case_when(
        
        is.na(Bias) ~ NA_real_,
        
        abs(True_Value) <= 1e-12 ~ NA_real_,
        
        TRUE ~ Bias / True_Value
      )
    ) %>%
    
    select(
      Sample_Size,
      Parameter,
      True_Value,
      Number_Successful,
      Mean_Estimate,
      Empirical_SD,
      Mean_Posterior_SD,
      Bias,
      Relative_Bias,
      MSE,
      RMSE,
      Coverage_95,
      MCSE_Mean_Estimate,
      Mean_Rhat,
      Maximum_Rhat,
      Mean_ESS,
      Convergence_Rate
    )
  
  return(
    list(
      individual_results =
        results_with_truth,
      
      summary_table =
        summary_table
    )
  )
}


############################################################
## 11. Optional small test
##
## Run this before the full simulation.
## This section does not alter the formal simulation settings.
############################################################

run_small_test <- FALSE

if (run_small_test) {
  
  test_results <- run_simulation_setting(
    N = 100,
    R = 3,
    n_time = n_time,
    n_chains = n_chains,
    n_adapt = 1000,
    n_burnin = 2000,
    n_iter = 3000,
    thin = thin,
    base_seed = 2026,
    save_directory = file.path(
      output_directory,
      "small_test"
    )
  )
  
  test_summary <- summarise_simulation_results(
    raw_results = test_results,
    true_parameter_values = true_values
  )
  
  print(
    test_summary$summary_table
  )
  
  write.csv(
    test_summary$summary_table,
    file = file.path(
      output_directory,
      "small_test_summary.csv"
    ),
    row.names = FALSE
  )
}


############################################################
## 12. Full simulation study
############################################################

full_results_list <- lapply(
  sample_sizes,
  function(current_N) {
    
    run_simulation_setting(
      N = current_N,
      R = R_simulation,
      n_time = n_time,
      n_chains = n_chains,
      n_adapt = n_adapt,
      n_burnin = n_burnin,
      n_iter = n_iter,
      thin = thin,
      base_seed = 2026,
      save_directory = file.path(
        output_directory,
        "full_results"
      )
    )
  }
)


############################################################
## 13. Combine all sample-size settings
############################################################

full_raw_results <- bind_rows(
  full_results_list
)


############################################################
## 14. Calculate simulation summaries
############################################################

full_output <- summarise_simulation_results(
  raw_results = full_raw_results,
  true_parameter_values = true_values
)

individual_results <-
  full_output$individual_results

final_summary_table <-
  full_output$summary_table


############################################################
## 15. Display final summary
############################################################

print(
  final_summary_table,
  n = Inf
)


############################################################
## 16. Save complete results
############################################################

write.csv(
  individual_results,
  file = file.path(
    output_directory,
    "longitudinal_individual_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_summary_table,
  file = file.path(
    output_directory,
    "longitudinal_complete_summary.csv"
  ),
  row.names = FALSE
)


############################################################
## 17. Save five core performance measures
############################################################

five_core_metrics <- final_summary_table %>%
  
  select(
    Sample_Size,
    Parameter,
    True_Value,
    Mean_Estimate,
    Empirical_SD,
    Relative_Bias,
    Coverage_95,
    MSE
  )

write.csv(
  five_core_metrics,
  file = file.path(
    output_directory,
    "longitudinal_five_core_metrics.csv"
  ),
  row.names = FALSE
)

print(
  five_core_metrics,
  n = Inf
)


############################################################
## 18. Rounded tables for thesis presentation
############################################################

final_summary_rounded <- final_summary_table %>%
  
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 4)
    )
  )

five_core_metrics_rounded <- five_core_metrics %>%
  
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 4)
    )
  )

write.csv(
  final_summary_rounded,
  file = file.path(
    output_directory,
    "longitudinal_complete_summary_rounded.csv"
  ),
  row.names = FALSE
)

write.csv(
  five_core_metrics_rounded,
  file = file.path(
    output_directory,
    "longitudinal_five_core_metrics_rounded.csv"
  ),
  row.names = FALSE
)


############################################################
## 19. Final message
############################################################

cat(
  "\nLongitudinal simulation study completed.\n"
)

cat(
  "Output directory:",
  normalizePath(
    output_directory,
    mustWork = FALSE
  ),
  "\n"
)
