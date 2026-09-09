rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

set.seed(20260810)

############################################################
## 0. PACKAGES
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
## 1. BASIC SETTINGS
############################################################

N <- 500
R <- 100

Tmax <- 15

## Planned number of longitudinal observations
n_planned <- 8

## Random censoring rate
lambda_c <- 0.02

############################################################
## 2. TRUE LONGITUDINAL PARAMETERS
############################################################

beta0_y_true <- 2
beta1_y_true <- 0.5
beta2_y_true <- 1

sigma_b0_true <- 1
sigma_b1_true <- 0.6
rho_true <- 0.3

sigma_y_true <- 1

############################################################
## Random-effects covariance matrix D
############################################################

cov_b01_true <-
  rho_true *
  sigma_b0_true *
  sigma_b1_true

D_true <- matrix(
  c(
    sigma_b0_true^2,
    cov_b01_true,
    cov_b01_true,
    sigma_b1_true^2
  ),
  nrow = 2,
  byrow = TRUE
)

print(D_true)

############################################################
## 3. TRUE MULTI-STATE PARAMETERS
############################################################

## Baseline transition intensities

lambda012_true <- 0.08
lambda013_true <- 0.04
lambda023_true <- 0.10

## Baseline covariate effects

gamma12_true <- -0.5
gamma13_true <- 0.8
gamma23_true <- 0.5

############################################################
## 4. TRUE ASSOCIATION PARAMETERS
############################################################

## Random intercept associations

alpha0_12_true <- 0.4
alpha0_13_true <- 0.8
alpha0_23_true <- 0.6

## Random slope associations

alpha1_12_true <- 0.15
alpha1_13_true <- 0.30
alpha1_23_true <- 0.25


############################################################
## 5. VECTOR OF TRUE PARAMETERS
############################################################

true_values <- c(
  
  beta0_y = beta0_y_true,
  beta1_y = beta1_y_true,
  beta2_y = beta2_y_true,
  
  sigma_b0 = sigma_b0_true,
  sigma_b1 = sigma_b1_true,
  rho = rho_true,
  sigma_y = sigma_y_true,
  
  lambda012 = lambda012_true,
  lambda013 = lambda013_true,
  lambda023 = lambda023_true,
  
  gamma12 = gamma12_true,
  gamma13 = gamma13_true,
  gamma23 = gamma23_true,
  
  alpha0_12 = alpha0_12_true,
  alpha0_13 = alpha0_13_true,
  alpha0_23 = alpha0_23_true,
  
  alpha1_12 = alpha1_12_true,
  alpha1_13 = alpha1_13_true,
  alpha1_23 = alpha1_23_true
)


############################################################
## 6. DATA-GENERATING FUNCTION
############################################################

simulate_joint_data <- function(
    N = 500,
    Tmax = 15,
    n_planned = 8,
    lambda_c = 0.02
) {
  
  ##########################################################
  ## 6.1 Baseline covariate
  ##########################################################
  
  x <- rnorm(
    N,
    mean = 0,
    sd = 1
  )
  
  ##########################################################
  ## 6.2 Generate random intercept and random slope
  ##########################################################
  
  b <- MASS::mvrnorm(
    n = N,
    mu = c(0, 0),
    Sigma = D_true
  )
  
  b0 <- b[, 1]
  b1 <- b[, 2]
  
  ##########################################################
  ## 6.3 Subject-specific transition intensities
  ##########################################################
  
  lambda12 <- lambda012_true *
    exp(
      gamma12_true * x +
        alpha0_12_true * b0 +
        alpha1_12_true * b1
    )
  
  lambda13 <- lambda013_true *
    exp(
      gamma13_true * x +
        alpha0_13_true * b0 +
        alpha1_13_true * b1
    )
  
  lambda23 <- lambda023_true *
    exp(
      gamma23_true * x +
        alpha0_23_true * b0 +
        alpha1_23_true * b1
    )
  
  ##########################################################
  ## 6.4 Random censoring
  ##########################################################
  
  C_star <- rexp(
    N,
    rate = lambda_c
  )
  
  C <- pmin(
    C_star,
    Tmax
  )
  
  ##########################################################
  ## 6.5 Objects for multi-state history
  ##########################################################
  
  ## Time spent in state 1
  time1 <- numeric(N)
  
  ## Time spent in state 2
  time2 <- numeric(N)
  
  ## Transition indicators
  d12 <- integer(N)
  d13 <- integer(N)
  d23 <- integer(N)
  
  ## Absorbing-state time
  T_absorb <- rep(Inf, N)
  
  ## Observed follow-up time
  follow_end <- numeric(N)
  
  ## Final observed state
  final_state <- integer(N)
  
  ## Transition times for descriptive purposes
  transition12_time <- rep(NA_real_, N)
  transition13_time <- rep(NA_real_, N)
  transition23_time <- rep(NA_real_, N)
  
  ##########################################################
  ## 6.6 Generate multi-state event histories
  ##########################################################
  
  for (i in seq_len(N)) {
    
    ########################################################
    ## Potential competing transitions from state 1
    ########################################################
    
    W12 <- rexp(
      1,
      rate = lambda12[i]
    )
    
    W13 <- rexp(
      1,
      rate = lambda13[i]
    )
    
    first_event_time <- min(
      W12,
      W13
    )
    
    ########################################################
    ## Case A:
    ## No transition before censoring
    ########################################################
    
    if (first_event_time >= C[i]) {
      
      time1[i] <- C[i]
      time2[i] <- 0
      
      d12[i] <- 0
      d13[i] <- 0
      d23[i] <- 0
      
      follow_end[i] <- C[i]
      final_state[i] <- 1
      
    } else {
      
      ######################################################
      ## Case B:
      ## Direct transition 1 -> 3
      ######################################################
      
      if (W13 < W12) {
        
        d13[i] <- 1
        
        transition13_time[i] <- W13
        
        time1[i] <- W13
        time2[i] <- 0
        
        T_absorb[i] <- W13
        follow_end[i] <- W13
        
        final_state[i] <- 3
        
      } else {
        
        ####################################################
        ## Case C:
        ## Transition 1 -> 2
        ####################################################
        
        d12[i] <- 1
        
        transition12_time[i] <- W12
        
        time1[i] <- W12
        
        ####################################################
        ## Generate waiting time for 2 -> 3
        ####################################################
        
        W23 <- rexp(
          1,
          rate = lambda23[i]
        )
        
        T23_absolute <- W12 + W23
        
        ####################################################
        ## 2 -> 3 before censoring
        ####################################################
        
        if (T23_absolute < C[i]) {
          
          d23[i] <- 1
          
          transition23_time[i] <- T23_absolute
          
          time2[i] <- W23
          
          T_absorb[i] <- T23_absolute
          follow_end[i] <- T23_absolute
          
          final_state[i] <- 3
          
        } else {
          
          ##################################################
          ## Censored while in state 2
          ##################################################
          
          d23[i] <- 0
          
          time2[i] <- C[i] - W12
          
          follow_end[i] <- C[i]
          
          final_state[i] <- 2
        }
      }
    }
  }
  
  ##########################################################
  ## 6.7 Generate planned longitudinal observation times
  ##########################################################
  
  longitudinal_list <- vector(
    "list",
    N
  )
  
  for (i in seq_len(N)) {
    
    ########################################################
    ## 8 planned observation times
    ########################################################
    
    planned_times <- sort(
      runif(
        n_planned,
        min = 0,
        max = Tmax
      )
    )
    
    ########################################################
    ## Measurements after follow-up end are unavailable
    ########################################################
    
    observed_times <- planned_times[
      planned_times <= follow_end[i]
    ]
    
    n_obs_i <- length(observed_times)
    
    if (n_obs_i > 0) {
      
      ######################################################
      ## True latent longitudinal trajectory
      ######################################################
      
      mu_y <-
        beta0_y_true +
        beta1_y_true * observed_times +
        beta2_y_true * x[i] +
        b0[i] +
        b1[i] * observed_times
      
      ######################################################
      ## Observed longitudinal outcome
      ######################################################
      
      y_i <- rnorm(
        n_obs_i,
        mean = mu_y,
        sd = sigma_y_true
      )
      
      longitudinal_list[[i]] <- data.frame(
        id = i,
        time = observed_times,
        y = y_i
      )
    }
  }
  
  ##########################################################
  ## Remove subjects with no observed longitudinal measure
  ## from longitudinal dataframe ONLY.
  ##
  ## They remain in the multi-state model.
  ##########################################################
  
  longitudinal_data <- do.call(
    rbind,
    longitudinal_list
  )
  
  if (is.null(longitudinal_data)) {
    stop("No longitudinal observations were generated.")
  }
  
  rownames(longitudinal_data) <- NULL
  
  ##########################################################
  ## Multi-state dataset
  ##########################################################
  
  multistate_data <- data.frame(
    
    id = seq_len(N),
    
    x = x,
    
    b0_true = b0,
    b1_true = b1,
    
    lambda12_true = lambda12,
    lambda13_true = lambda13,
    lambda23_true = lambda23,
    
    C_star = C_star,
    C = C,
    
    time1 = time1,
    time2 = time2,
    
    d12 = d12,
    d13 = d13,
    d23 = d23,
    
    transition12_time = transition12_time,
    transition13_time = transition13_time,
    transition23_time = transition23_time,
    
    absorbing_time = T_absorb,
    follow_end = follow_end,
    final_state = final_state
  )
  
  ##########################################################
  ## Return
  ##########################################################
  
  list(
    longitudinal = longitudinal_data,
    multistate = multistate_data
  )
}


############################################################
## 7. TEST ONE SIMULATED DATASET
############################################################

test_data <- simulate_joint_data(
  N = N,
  Tmax = Tmax,
  n_planned = n_planned,
  lambda_c = lambda_c
)

cat("\n============================================\n")
cat("SIMULATED DATA CHECK\n")
cat("============================================\n")

cat(
  "\nTotal longitudinal observations:",
  nrow(test_data$longitudinal),
  "\n"
)

cat(
  "Mean longitudinal observations per subject:",
  nrow(test_data$longitudinal) / N,
  "\n"
)

cat(
  "\nTransitions 1 -> 2:",
  sum(test_data$multistate$d12),
  "\n"
)

cat(
  "Transitions 1 -> 3:",
  sum(test_data$multistate$d13),
  "\n"
)

cat(
  "Transitions 2 -> 3:",
  sum(test_data$multistate$d23),
  "\n"
)

cat("\nFinal states:\n")

print(
  table(test_data$multistate$final_state)
)


############################################################
## 8. JAGS MODEL
############################################################

jags_model_string <- "
model {

############################################################
## A. LONGITUDINAL SUBMODEL
############################################################

for (j in 1:Nobs) {

    mu_y[j] <-
        beta0_y +
        beta1_y * obs_time[j] +
        beta2_y * x[obs_id[j]] +
        b[obs_id[j],1] +
        b[obs_id[j],2] * obs_time[j]

    y[j] ~ dnorm(
        mu_y[j],
        tau_y
    )
}

############################################################
## B. RANDOM EFFECTS
############################################################

for (i in 1:N) {

    b[i,1:2] ~ dmnorm(
        zero_b[],
        Tau_b[,]
    )
}

############################################################
## C. MULTI-STATE SUBMODEL
############################################################

for (i in 1:N) {

    ########################################################
    ## log transition intensities
    ########################################################

    log_lambda12[i] <-
        log_lambda012 +
        gamma12 * x[i] +
        alpha0_12 * b[i,1] +
        alpha1_12 * b[i,2]

    log_lambda13[i] <-
        log_lambda013 +
        gamma13 * x[i] +
        alpha0_13 * b[i,1] +
        alpha1_13 * b[i,2]

    log_lambda23[i] <-
        log_lambda023 +
        gamma23 * x[i] +
        alpha0_23 * b[i,1] +
        alpha1_23 * b[i,2]

    ########################################################
    ## Transition intensities
    ########################################################

    lambda12[i] <- exp(log_lambda12[i])
    lambda13[i] <- exp(log_lambda13[i])
    lambda23[i] <- exp(log_lambda23[i])

    ########################################################
    ## Individual multi-state log likelihood
    ##
    ## State 1:
    ## d12 log(lambda12)
    ## + d13 log(lambda13)
    ## - (lambda12 + lambda13) time1
    ##
    ## State 2:
    ## d23 log(lambda23)
    ## - lambda23 time2
    ########################################################

    loglik_ms[i] <-

        d12[i] * log_lambda12[i] +

        d13[i] * log_lambda13[i] -

        (lambda12[i] + lambda13[i]) *
        time1[i] +

        d23[i] * log_lambda23[i] -

        lambda23[i] *
        time2[i]

    ########################################################
    ## Zero trick
    ########################################################

    phi[i] <- -loglik_ms[i] + Czero

    zeros[i] ~ dpois(phi[i])
}

############################################################
## D. LONGITUDINAL FIXED-EFFECT PRIORS
############################################################

beta0_y ~ dnorm(0, 0.0001)
beta1_y ~ dnorm(0, 0.0001)
beta2_y ~ dnorm(0, 0.0001)

############################################################
## E. RESIDUAL VARIANCE
##
## sigma_y^2 ~ IG(0.001, 0.001)
## equivalent:
## tau_y ~ Gamma(0.001, 0.001)
############################################################

tau_y ~ dgamma(0.001, 0.001)

sigma2_y <- 1 / tau_y
sigma_y <- sqrt(sigma2_y)

############################################################
## F. RANDOM-EFFECT COVARIANCE PRIOR
##
## D ~ IW(S0, nu0)
##
## Tau_b = inverse(D)
############################################################

Tau_b[1:2,1:2] ~ dwish(
    S0[,],
    nu0
)

D[1:2,1:2] <- inverse(
    Tau_b[,]
)

sigma_b0 <- sqrt(D[1,1])
sigma_b1 <- sqrt(D[2,2])

rho <-
    D[1,2] /
    sqrt(
        D[1,1] *
        D[2,2]
    )

############################################################
## G. BASELINE TRANSITION INTENSITIES
############################################################

lambda012 ~ dlnorm(0, 1)
lambda013 ~ dlnorm(0, 1)
lambda023 ~ dlnorm(0, 1)

log_lambda012 <- log(lambda012)
log_lambda013 <- log(lambda013)
log_lambda023 <- log(lambda023)

############################################################
## H. TRANSITION-SPECIFIC REGRESSION COEFFICIENTS
##
## gamma ~ N(0,4)
##
## variance = 4
## precision = 1/4 = 0.25
############################################################

gamma12 ~ dnorm(0, 0.25)
gamma13 ~ dnorm(0, 0.25)
gamma23 ~ dnorm(0, 0.25)

############################################################
## I. ASSOCIATION PARAMETERS
##
## alpha ~ N(0,100)
##
## variance = 100
## precision = 0.01
############################################################

alpha0_12 ~ dnorm(0, 0.01)
alpha0_13 ~ dnorm(0, 0.01)
alpha0_23 ~ dnorm(0, 0.01)

alpha1_12 ~ dnorm(0, 0.01)
alpha1_13 ~ dnorm(0, 0.01)
alpha1_23 ~ dnorm(0, 0.01)

}
"


############################################################
## 9. WRITE JAGS MODEL TO FILE
############################################################

model_file <- "joint_shared_RE_model.jags"

writeLines(
  jags_model_string,
  con = model_file
)


############################################################
## 10. PREPARE JAGS DATA
############################################################

prepare_jags_data <- function(sim_data) {
  
  long <- sim_data$longitudinal
  ms <- sim_data$multistate
  
  ##########################################################
  ## Inverse-Wishart hyperparameters
  ##########################################################
  
  S0 <- matrix(
    c(
      1, 0,
      0, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  list(
    
    N = nrow(ms),
    
    ########################################################
    ## Longitudinal data
    ########################################################
    
    Nobs = nrow(long),
    
    y = long$y,
    
    obs_time = long$time,
    
    obs_id = as.integer(long$id),
    
    ########################################################
    ## Baseline covariate
    ########################################################
    
    x = ms$x,
    
    ########################################################
    ## Multi-state data
    ########################################################
    
    time1 = ms$time1,
    
    time2 = ms$time2,
    
    d12 = ms$d12,
    
    d13 = ms$d13,
    
    d23 = ms$d23,
    
    ########################################################
    ## Zero trick
    ########################################################
    
    zeros = rep(
      0,
      nrow(ms)
    ),
    
    Czero = 1000,
    
    ########################################################
    ## Random-effects prior
    ########################################################
    
    zero_b = c(0, 0),
    
    S0 = S0,
    
    nu0 = 4
  )
}


############################################################
## 11. PARAMETERS TO MONITOR
############################################################

parameters_to_monitor <- c(
  
  ##########################################################
  ## Longitudinal
  ##########################################################
  
  "beta0_y",
  "beta1_y",
  "beta2_y",
  
  "sigma_b0",
  "sigma_b1",
  "rho",
  
  "sigma_y",
  
  ##########################################################
  ## Multi-state baseline intensities
  ##########################################################
  
  "lambda012",
  "lambda013",
  "lambda023",
  
  ##########################################################
  ## Covariate effects
  ##########################################################
  
  "gamma12",
  "gamma13",
  "gamma23",
  
  ##########################################################
  ## Association parameters
  ##########################################################
  
  "alpha0_12",
  "alpha0_13",
  "alpha0_23",
  
  "alpha1_12",
  "alpha1_13",
  "alpha1_23"
)


############################################################
## 12. INITIAL VALUES
############################################################

make_inits <- function(chain_id) {
  
  list(
    
    beta0_y = rnorm(1, 2, 0.1),
    beta1_y = rnorm(1, 0.5, 0.05),
    beta2_y = rnorm(1, 1, 0.1),
    
    tau_y = rgamma(
      1,
      shape = 2,
      rate = 2
    ),
    
    lambda012 = exp(
      rnorm(
        1,
        log(0.08),
        0.1
      )
    ),
    
    lambda013 = exp(
      rnorm(
        1,
        log(0.04),
        0.1
      )
    ),
    
    lambda023 = exp(
      rnorm(
        1,
        log(0.10),
        0.1
      )
    ),
    
    gamma12 = rnorm(
      1,
      -0.5,
      0.1
    ),
    
    gamma13 = rnorm(
      1,
      0.8,
      0.1
    ),
    
    gamma23 = rnorm(
      1,
      0.5,
      0.1
    ),
    
    alpha0_12 = rnorm(
      1,
      0.4,
      0.05
    ),
    
    alpha0_13 = rnorm(
      1,
      0.8,
      0.05
    ),
    
    alpha0_23 = rnorm(
      1,
      0.6,
      0.05
    ),
    
    alpha1_12 = rnorm(
      1,
      0.15,
      0.05
    ),
    
    alpha1_13 = rnorm(
      1,
      0.30,
      0.05
    ),
    
    alpha1_23 = rnorm(
      1,
      0.25,
      0.05
    ),
    
    .RNG.name = "base::Wichmann-Hill",
    
    .RNG.seed =
      10000 +
      chain_id * 100 +
      sample(
        1:99,
        1
      )
  )
}


############################################################
## 13. FIT ONE JAGS MODEL
############################################################

fit_joint_model <- function(
    sim_data,
    n_adapt = 2000,
    n_burn = 5000,
    n_iter = 10000,
    n_chains = 2
) {
  
  jags_data <- prepare_jags_data(
    sim_data
  )
  
  inits_list <- lapply(
    seq_len(n_chains),
    make_inits
  )
  
  ##########################################################
  ## Create model
  ##########################################################
  
  model <- rjags::jags.model(
    
    file = model_file,
    
    data = jags_data,
    
    inits = inits_list,
    
    n.chains = n_chains,
    
    n.adapt = n_adapt,
    
    quiet = TRUE
  )
  
  ##########################################################
  ## Burn-in
  ##########################################################
  
  update(
    model,
    n.iter = n_burn,
    progress.bar = "none"
  )
  
  ##########################################################
  ## Posterior sampling
  ##########################################################
  
  samples <- rjags::coda.samples(
    
    model = model,
    
    variable.names =
      parameters_to_monitor,
    
    n.iter = n_iter,
    
    thin = 1,
    
    progress.bar = "none"
  )
  
  samples
}


############################################################
## 14. EXTRACT POSTERIOR ESTIMATES
############################################################

extract_estimates <- function(samples) {
  
  sample_matrix <- as.matrix(
    samples
  )
  
  estimates <- data.frame(
    
    Parameter =
      parameters_to_monitor,
    
    Estimate = NA_real_,
    
    Posterior_SD = NA_real_,
    
    Lower95 = NA_real_,
    
    Upper95 = NA_real_,
    
    Rhat = NA_real_,
    
    ESS = NA_real_,
    
    stringsAsFactors = FALSE
  )
  
  ##########################################################
  ## Posterior summaries
  ##########################################################
  
  for (p in parameters_to_monitor) {
    
    values <- sample_matrix[, p]
    
    estimates[
      estimates$Parameter == p,
      "Estimate"
    ] <- mean(values)
    
    estimates[
      estimates$Parameter == p,
      "Posterior_SD"
    ] <- sd(values)
    
    estimates[
      estimates$Parameter == p,
      "Lower95"
    ] <- quantile(
      values,
      0.025
    )
    
    estimates[
      estimates$Parameter == p,
      "Upper95"
    ] <- quantile(
      values,
      0.975
    )
  }
  
  ##########################################################
  ## R-hat
  ##########################################################
  
  rhat_result <- tryCatch(
    
    coda::gelman.diag(
      samples,
      multivariate = FALSE,
      autoburnin = FALSE
    )$psrf[, 1],
    
    error = function(e) {
      rep(
        NA_real_,
        length(parameters_to_monitor)
      )
    }
  )
  
  if (
    length(rhat_result) ==
    length(parameters_to_monitor)
  ) {
    
    estimates$Rhat <-
      rhat_result[
        estimates$Parameter
      ]
  }
  
  ##########################################################
  ## Effective sample size
  ##########################################################
  
  ess_result <- tryCatch(
    
    coda::effectiveSize(
      samples
    ),
    
    error = function(e) {
      rep(
        NA_real_,
        length(parameters_to_monitor)
      )
    }
  )
  
  if (
    length(ess_result) > 0
  ) {
    
    estimates$ESS <-
      ess_result[
        estimates$Parameter
      ]
  }
  
  estimates
}


############################################################
## 15. RUN ONE REPLICATION
############################################################

run_one_replication <- function(
    replicate_id,
    N = 500
) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "N =", N,
    "| replicate =", replicate_id,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  ##########################################################
  ## Reproducibility
  ##########################################################
  
  set.seed(
    100000 +
      N * 1000 +
      replicate_id
  )
  
  ##########################################################
  ## Generate data
  ##########################################################
  
  sim_data <- simulate_joint_data(
    
    N = N,
    
    Tmax = Tmax,
    
    n_planned = n_planned,
    
    lambda_c = lambda_c
  )
  
  ##########################################################
  ## Descriptive information
  ##########################################################
  
  n_long <- nrow(
    sim_data$longitudinal
  )
  
  n12 <- sum(
    sim_data$multistate$d12
  )
  
  n13 <- sum(
    sim_data$multistate$d13
  )
  
  n23 <- sum(
    sim_data$multistate$d23
  )
  
  cat(
    "Longitudinal observations:",
    n_long,
    "\n"
  )
  
  cat(
    "Transitions:",
    "12 =", n12,
    "| 13 =", n13,
    "| 23 =", n23,
    "\n"
  )
  
  ##########################################################
  ## Fit model
  ##########################################################
  
  fit_result <- tryCatch(
    
    {
      
      samples <- fit_joint_model(
        
        sim_data = sim_data,
        
        n_adapt = 30000,
        
        n_burn = 20000,
        
        n_iter = 50000,
        
        n_chains = 2
      )
      
      estimates <- extract_estimates(
        samples
      )
      
      ######################################################
      ## Check estimates
      ######################################################
      
      valid <-
        all(
          is.finite(
            estimates$Estimate
          )
        ) &&
        all(
          is.finite(
            estimates$Lower95
          )
        ) &&
        all(
          is.finite(
            estimates$Upper95
          )
        )
      
      if (!valid) {
        stop(
          "Non-finite posterior estimates detected."
        )
      }
      
      ######################################################
      ## Additional diagnostics
      ######################################################
      
      max_rhat <- suppressWarnings(
        max(
          estimates$Rhat,
          na.rm = TRUE
        )
      )
      
      min_ess <- suppressWarnings(
        min(
          estimates$ESS,
          na.rm = TRUE
        )
      )
      
      estimates$Replicate <-
        replicate_id
      
      estimates$N <-
        N
      
      estimates$N_long <-
        n_long
      
      estimates$N12 <-
        n12
      
      estimates$N13 <-
        n13
      
      estimates$N23 <-
        n23
      
      estimates$Max_Rhat <-
        max_rhat
      
      estimates$Min_ESS <-
        min_ess
      
      estimates$Fit_Status <-
        "success"
      
      estimates
    },
    
    error = function(e) {
      
      cat(
        "\nFit failed:",
        "N =", N,
        "replicate =", replicate_id,
        "error =", conditionMessage(e),
        "\n"
      )
      
      data.frame(
        
        Parameter =
          parameters_to_monitor,
        
        Estimate =
          NA_real_,
        
        Posterior_SD =
          NA_real_,
        
        Lower95 =
          NA_real_,
        
        Upper95 =
          NA_real_,
        
        Rhat =
          NA_real_,
        
        ESS =
          NA_real_,
        
        Replicate =
          replicate_id,
        
        N =
          N,
        
        N_long =
          n_long,
        
        N12 =
          n12,
        
        N13 =
          n13,
        
        N23 =
          n23,
        
        Max_Rhat =
          NA_real_,
        
        Min_ESS =
          NA_real_,
        
        Fit_Status =
          paste0(
            "failed: ",
            conditionMessage(e)
          )
      )
    }
  )
  
  ##########################################################
  ## Save replication immediately
  ##
  ## Important for CSF:
  ## results survive if job stops halfway.
  ##########################################################
  
  dir.create(
    "replicate_results",
    showWarnings = FALSE
  )
  
  write.csv(
    
    fit_result,
    
    file = file.path(
      "replicate_results",
      paste0(
        "N",
        N,
        "_rep",
        sprintf(
          "%03d",
          replicate_id
        ),
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  gc()
  
  fit_result
}


############################################################
## 16. RUN ALL MONTE CARLO REPLICATIONS
############################################################

all_results_list <- vector(
  "list",
  R
)

for (r in seq_len(R)) {
  
  all_results_list[[r]] <-
    run_one_replication(
      
      replicate_id = r,
      
      N = N
    )
  
  ##########################################################
  ## Save cumulative results after every replication
  ##########################################################
  
  current_results <- do.call(
    rbind,
    all_results_list[
      seq_len(r)
    ]
  )
  
  write.csv(
    
    current_results,
    
    file =
      paste0(
        "joint_raw_results_N",
        N,
        ".csv"
      ),
    
    row.names = FALSE
  )
}


############################################################
## 17. COMBINE RESULTS
############################################################

all_results <- do.call(
  rbind,
  all_results_list
)

write.csv(
  
  all_results,
  
  file =
    paste0(
      "joint_raw_results_N",
      N,
      ".csv"
    ),
  
  row.names = FALSE
)


############################################################
## 18. FAILED-FIT SUMMARY
############################################################

fit_status <- all_results %>%
  
  select(
    Replicate,
    Fit_Status
  ) %>%
  
  distinct()


cat("\n============================================\n")
cat("FIT STATUS\n")
cat("============================================\n")

print(
  table(
    fit_status$Fit_Status
  )
)

n_success <- sum(
  fit_status$Fit_Status ==
    "success"
)

n_failed <-
  R -
  n_success

cat(
  "\nSuccessful fits:",
  n_success,
  "\n"
)

cat(
  "Failed fits:",
  n_failed,
  "\n"
)


############################################################
## 19. PARAMETER NAME MATCHING
############################################################

parameter_name_map <- c(
  
  beta0_y = "beta0_y",
  beta1_y = "beta1_y",
  beta2_y = "beta2_y",
  
  sigma_b0 = "sigma_b0",
  sigma_b1 = "sigma_b1",
  rho = "rho",
  sigma_y = "sigma_y",
  
  lambda012 = "lambda012",
  lambda013 = "lambda013",
  lambda023 = "lambda023",
  
  gamma12 = "gamma12",
  gamma13 = "gamma13",
  gamma23 = "gamma23",
  
  alpha0_12 = "alpha0_12",
  alpha0_13 = "alpha0_13",
  alpha0_23 = "alpha0_23",
  
  alpha1_12 = "alpha1_12",
  alpha1_13 = "alpha1_13",
  alpha1_23 = "alpha1_23"
)


############################################################
## 20. PERFORMANCE MEASURES
############################################################

calculate_performance <- function(
    all_results,
    true_values
) {
  
  ##########################################################
  ## Only successful fits
  ##########################################################
  
  successful <- all_results %>%
    filter(
      Fit_Status == "success"
    )
  
  performance_list <- list()
  
  counter <- 1
  
  for (
    parameter_name
    in names(true_values)
  ) {
    
    true_value <-
      true_values[
        parameter_name
      ]
    
    parameter_data <-
      successful %>%
      filter(
        Parameter ==
          parameter_name
      )
    
    M <- nrow(
      parameter_data
    )
    
    ########################################################
    ## No valid fits
    ########################################################
    
    if (M == 0) {
      
      performance_list[[counter]] <-
        data.frame(
          
          Parameter =
            parameter_name,
          
          True_Value =
            true_value,
          
          M =
            0,
          
          Mean_Estimate =
            NA_real_,
          
          Empirical_SD =
            NA_real_,
          
          Mean_Posterior_SD =
            NA_real_,
          
          Bias =
            NA_real_,
          
          Relative_Bias =
            NA_real_,
          
          MSE =
            NA_real_,
          
          RMSE =
            NA_real_,
          
          Coverage =
            NA_real_,
          
          MCSE =
            NA_real_
        )
      
      counter <- counter + 1
      
      next
    }
    
    ########################################################
    ## Estimates
    ########################################################
    
    estimates <-
      parameter_data$Estimate
    
    ########################################################
    ## Mean estimate
    ########################################################
    
    mean_estimate <-
      mean(
        estimates
      )
    
    ########################################################
    ## Empirical standard deviation
    ########################################################
    
    empirical_sd <-
      sd(
        estimates
      )
    
    ########################################################
    ## Mean posterior standard deviation
    ########################################################
    
    mean_posterior_sd <-
      mean(
        parameter_data$Posterior_SD
      )
    
    ########################################################
    ## Bias
    ########################################################
    
    bias <-
      mean_estimate -
      true_value
    
    ########################################################
    ## Relative bias
    ########################################################
    
    if (
      abs(true_value) <
      .Machine$double.eps
    ) {
      
      relative_bias <-
        NA_real_
      
    } else {
      
      relative_bias <-
        100 *
        bias /
        true_value
    }
    
    ########################################################
    ## MSE
    ########################################################
    
    mse <-
      mean(
        (
          estimates -
            true_value
        )^2
      )
    
    ########################################################
    ## RMSE
    ########################################################
    
    rmse <-
      sqrt(
        mse
      )
    
    ########################################################
    ## Coverage
    ########################################################
    
    coverage <-
      mean(
        
        parameter_data$Lower95 <=
          true_value &
          
          parameter_data$Upper95 >=
          true_value
      )
    
    ########################################################
    ## Monte Carlo standard error of mean estimate
    ########################################################
    
    mcse <-
      empirical_sd /
      sqrt(M)
    
    ########################################################
    ## Store
    ########################################################
    
    performance_list[[counter]] <-
      data.frame(
        
        Parameter =
          parameter_name,
        
        True_Value =
          true_value,
        
        M =
          M,
        
        Mean_Estimate =
          mean_estimate,
        
        Empirical_SD =
          empirical_sd,
        
        Mean_Posterior_SD =
          mean_posterior_sd,
        
        Bias =
          bias,
        
        Relative_Bias =
          relative_bias,
        
        MSE =
          mse,
        
        RMSE =
          rmse,
        
        Coverage =
          coverage,
        
        MCSE =
          mcse
      )
    
    counter <- counter + 1
  }
  
  do.call(
    rbind,
    performance_list
  )
}


############################################################
## 21. FINAL PERFORMANCE TABLE
############################################################

performance_results <-
  calculate_performance(
    
    all_results =
      all_results,
    
    true_values =
      true_values
  )


############################################################
## 22. ADD SAMPLE SIZE AND FIT INFORMATION
############################################################

performance_results$Sample_Size <-
  N

performance_results$Successful_Fits <-
  n_success

performance_results$Failed_Fits <-
  n_failed

performance_results <-
  performance_results %>%
  
  select(
    
    Sample_Size,
    
    Parameter,
    
    True_Value,
    
    Mean_Estimate,
    
    Mean_Posterior_SD,
    
    Empirical_SD,
    
    Bias,
    
    Relative_Bias,
    
    MSE,
    
    RMSE,
    
    Coverage,
    
    MCSE,
    
    M,
    
    Successful_Fits,
    
    Failed_Fits
  )


############################################################
## 23. SAVE FINAL RESULTS
############################################################

write.csv(
  
  performance_results,
  
  file =
    paste0(
      "joint_performance_N",
      N,
      ".csv"
    ),
  
  row.names = FALSE
)


############################################################
## 24. CONVERGENCE SUMMARY
############################################################

convergence_summary <-
  
  all_results %>%
  
  filter(
    Fit_Status ==
      "success"
  ) %>%
  
  group_by(
    Parameter
  ) %>%
  
  summarise(
    
    Mean_Rhat =
      mean(
        Rhat,
        na.rm = TRUE
      ),
    
    Max_Rhat =
      max(
        Rhat,
        na.rm = TRUE
      ),
    
    Mean_ESS =
      mean(
        ESS,
        na.rm = TRUE
      ),
    
    Min_ESS =
      min(
        ESS,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


write.csv(
  
  convergence_summary,
  
  file =
    paste0(
      "joint_convergence_N",
      N,
      ".csv"
    ),
  
  row.names = FALSE
)


############################################################
## 25. FINAL PRINT
############################################################

cat("\n\n============================================\n")
cat("FINAL PERFORMANCE RESULTS\n")
cat("============================================\n\n")

print(
  performance_results
)

cat("\n\n============================================\n")
cat("CONVERGENCE SUMMARY\n")
cat("============================================\n\n")

print(
  convergence_summary
)

cat("\n\nSimulation completed.\n")