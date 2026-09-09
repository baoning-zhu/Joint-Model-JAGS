############################################################
## Conventional multi-state simulation under MODERATE shared RE DGM
## States: 1 -> 2, 1 -> 3, 2 -> 3
##
## Data-generating model contains shared random intercept and
## random slope. The fitted conventional multi-state model
## deliberately omits them and uses only x.
##
## Moderate-association scenario:
## alpha0_rs = 0.4, alpha1_rs = 0.2,
## sigma_b0 = 0.7, sigma_b1 = 0.3, rho = 0.3.
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260808)

############################################################
## 0. Packages
############################################################
required_packages <- c("MASS", "rjags", "coda")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(MASS)
library(rjags)
library(coda)

############################################################
## 1. Simulation settings
############################################################
true_par <- c(
  lambda012 = 0.08,
  lambda013 = 0.04,
  lambda023 = 0.10,
  gamma12   = -0.50,
  gamma13   =  0.80,
  gamma23   =  0.50,
  alpha012  =  0.40,
  alpha013  =  0.80,
  alpha023  =  0.60,
  alpha112  =  0.15,
  alpha113  =  0.30,
  alpha123  =  0.25,
  sigma_b0  =  0.70,
  sigma_b1  =  0.30,
  rho       =  0.30
)

sample_sizes <- c(500L)
R_rep        <- 100L
Tmax         <- 15
censor_rate  <- 0.08

## MCMC settings in the written plan
n_chains <- 3L
n_adapt  <- 1000L
n_burn   <- 2000L
n_iter   <- 5000L
n_thin   <- 1L

## For a quick code test, set TRUE. Full study uses FALSE.
quick_test <- FALSE
if (quick_test) {
  sample_sizes <- 100L
  R_rep  <- 2L
  n_adapt <- 200L
  n_burn  <- 300L
  n_iter  <- 500L
}

output_dir <- "multistate_simulation_results_moderate"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. Utility functions
############################################################
make_D <- function(sigma_b0, sigma_b1, rho) {
  matrix(
    c(
      sigma_b0^2, rho * sigma_b0 * sigma_b1,
      rho * sigma_b0 * sigma_b1, sigma_b1^2
    ),
    nrow = 2L,
    byrow = TRUE
  )
}

## Generate one common dataset from the true shared-RE mechanism.
simulate_multistate_data <- function(
    N,
    pars = true_par,
    Tmax = 15,
    censor_rate = 0.08
) {
  stopifnot(N >= 1L, Tmax > 0, censor_rate > 0)

  x <- rnorm(N, mean = 0, sd = 1)

  D <- make_D(
    sigma_b0 = pars["sigma_b0"],
    sigma_b1 = pars["sigma_b1"],
    rho      = pars["rho"]
  )

  b <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = D)
  b0 <- b[, 1]
  b1 <- b[, 2]

  lambda12 <- unname(pars["lambda012"]) * exp(
    unname(pars["gamma12"])  * x +
      unname(pars["alpha012"]) * b0 +
      unname(pars["alpha112"]) * b1
  )
  lambda13 <- unname(pars["lambda013"]) * exp(
    unname(pars["gamma13"])  * x +
      unname(pars["alpha013"]) * b0 +
      unname(pars["alpha113"]) * b1
  )
  lambda23 <- unname(pars["lambda023"]) * exp(
    unname(pars["gamma23"])  * x +
      unname(pars["alpha023"]) * b0 +
      unname(pars["alpha123"]) * b1
  )

  W1 <- rexp(N, rate = lambda12 + lambda13)
  p12 <- lambda12 / (lambda12 + lambda13)
  first_destination <- ifelse(runif(N) < p12, 2L, 3L)
  W2 <- rexp(N, rate = lambda23)

  C_random <- rexp(N, rate = censor_rate)
  follow_end <- pmin(C_random, Tmax)

  rows <- vector("list", 2L * N)
  row_id <- 0L

  subject <- data.frame(
    id = seq_len(N), x = x, b0 = b0, b1 = b1,
    lambda12 = lambda12, lambda13 = lambda13, lambda23 = lambda23,
    W1 = W1, first_destination = first_destination, W2 = W2,
    censor_time = C_random, follow_end = follow_end
  )

  for (i in seq_len(N)) {
    ## No first transition observed: censored in state 1.
    if (W1[i] >= follow_end[i]) {
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        id = i,
        start = 0,
        stop = follow_end[i],
        duration = follow_end[i],
        from = 1L,
        to = NA_integer_,
        event12 = 0L,
        event13 = 0L,
        event23 = 0L,
        x = x[i]
      )
      next
    }

    ## First transition is observed.
    if (first_destination[i] == 3L) {
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        id = i,
        start = 0,
        stop = W1[i],
        duration = W1[i],
        from = 1L,
        to = 3L,
        event12 = 0L,
        event13 = 1L,
        event23 = 0L,
        x = x[i]
      )
      next
    }

    ## Observed 1 -> 2 interval.
    row_id <- row_id + 1L
    rows[[row_id]] <- data.frame(
      id = i,
      start = 0,
      stop = W1[i],
      duration = W1[i],
      from = 1L,
      to = 2L,
      event12 = 1L,
      event13 = 0L,
      event23 = 0L,
      x = x[i]
    )

    ## Potential 2 -> 3 transition at absolute time W1 + W2.
    T3 <- W1[i] + W2[i]
    row_id <- row_id + 1L

    if (T3 < follow_end[i]) {
      rows[[row_id]] <- data.frame(
        id = i,
        start = W1[i],
        stop = T3,
        duration = W2[i],
        from = 2L,
        to = 3L,
        event12 = 0L,
        event13 = 0L,
        event23 = 1L,
        x = x[i]
      )
    } else {
      rows[[row_id]] <- data.frame(
        id = i,
        start = W1[i],
        stop = follow_end[i],
        duration = follow_end[i] - W1[i],
        from = 2L,
        to = NA_integer_,
        event12 = 0L,
        event13 = 0L,
        event23 = 0L,
        x = x[i]
      )
    }
  }

  start_stop <- do.call(rbind, rows[seq_len(row_id)])
  rownames(start_stop) <- NULL

  if (any(start_stop$duration <= 0)) {
    stop("Non-positive risk interval generated.")
  }

  list(
    subject = subject,
    start_stop = start_stop,
    D = D,
    diagnostics = c(
      proportion_randomly_or_admin_censored_before_first_transition =
        mean(W1 >= follow_end),
      observed_12 = sum(start_stop$event12),
      observed_13 = sum(start_stop$event13),
      observed_23 = sum(start_stop$event23)
    )
  )
}

############################################################
## 3. JAGS model: conventional multi-state model
############################################################
## The model omits b0 and b1 by design.
##
## For an interval in state 1:
## log L = d12*log(lambda12) + d13*log(lambda13)
##         - (lambda12 + lambda13)*duration
## For an interval in state 2:
## log L = d23*log(lambda23) - lambda23*duration
##
## C_zero is an additive constant and does not affect the posterior.
## A moderate value is safer numerically than 10,000.
jags_model_string <- "
model {
  for (k in 1:K) {
    lambda12[k] <- lambda012 * exp(gamma12 * x[k])
    lambda13[k] <- lambda013 * exp(gamma13 * x[k])
    lambda23[k] <- lambda023 * exp(gamma23 * x[k])

    loglik_state1[k] <-
      event12[k] * log(lambda12[k]) +
      event13[k] * log(lambda13[k]) -
      (lambda12[k] + lambda13[k]) * duration[k]

    loglik_state2[k] <-
      event23[k] * log(lambda23[k]) -
      lambda23[k] * duration[k]

    loglik[k] <- is_state1[k] * loglik_state1[k] +
                 is_state2[k] * loglik_state2[k]

    phi[k] <- -loglik[k] + C_zero
    zeros[k] ~ dpois(phi[k])
  }

  log_lambda012 ~ dnorm(0, 1)
  log_lambda013 ~ dnorm(0, 1)
  log_lambda023 ~ dnorm(0, 1)

  lambda012 <- exp(log_lambda012)
  lambda013 <- exp(log_lambda013)
  lambda023 <- exp(log_lambda023)

  gamma12 ~ dnorm(0, 0.25)
  gamma13 ~ dnorm(0, 0.25)
  gamma23 ~ dnorm(0, 0.25)
}
"

model_file <- file.path(output_dir, "conventional_multistate_model.jags")
writeLines(jags_model_string, con = model_file)

prepare_jags_data <- function(start_stop, C_zero = 10000) {
  K <- nrow(start_stop)

  ## Check that phi = -loglik + C_zero should remain positive near
  ## plausible values. The model itself still requires C_zero > max(loglik).
  list(
    K = K,
    x = start_stop$x,
    duration = start_stop$duration,
    event12 = start_stop$event12,
    event13 = start_stop$event13,
    event23 = start_stop$event23,
    is_state1 = as.integer(start_stop$from == 1L),
    is_state2 = as.integer(start_stop$from == 2L),
    zeros = rep(0L, K),
    C_zero = C_zero
  )
}

make_inits <- function(chain_id) {
  list(
    log_lambda012 = rnorm(1, log(0.10), 0.15),
    log_lambda013 = rnorm(1, log(0.05), 0.15),
    log_lambda023 = rnorm(1, log(0.15), 0.15),
    gamma12 = rnorm(1, -0.3, 0.15),
    gamma13 = rnorm(1,  0.3, 0.15),
    gamma23 = rnorm(1,  0.3, 0.15),
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = 10000L + chain_id * 1000L + sample.int(999L, 1L)
  )
}

fit_conventional_model <- function(start_stop) {
  data_jags <- prepare_jags_data(start_stop, C_zero = 100)
  inits <- lapply(seq_len(n_chains), make_inits)

  jm <- rjags::jags.model(
    file = model_file,
    data = data_jags,
    inits = inits,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )

  update(jm, n.iter = n_burn, progress.bar = "none")

  monitored <- c(
    "lambda012", "lambda013", "lambda023",
    "gamma12", "gamma13", "gamma23"
  )

  samples <- rjags::coda.samples(
    model = jm,
    variable.names = monitored,
    n.iter = n_iter,
    thin = n_thin,
    progress.bar = "none"
  )

  sm <- summary(samples)
  stat <- sm$statistics
  quant <- sm$quantiles

  ## Ensure a fixed parameter order.
  param_order <- monitored
  stat <- stat[param_order, , drop = FALSE]
  quant <- quant[param_order, , drop = FALSE]

  gelman <- tryCatch(
    coda::gelman.diag(samples, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1],
    error = function(e) setNames(rep(NA_real_, length(param_order)), param_order)
  )
  gelman <- gelman[param_order]

  ess <- coda::effectiveSize(samples)
  ess <- ess[param_order]

  data.frame(
    parameter = param_order,
    estimate = stat[, "Mean"],
    posterior_sd = stat[, "SD"],
    lower = quant[, "2.5%"],
    median = quant[, "50%"],
    upper = quant[, "97.5%"],
    rhat = as.numeric(gelman),
    ess = as.numeric(ess),
    row.names = NULL
  )
}

############################################################
## 4. One-replicate wrapper with error handling
############################################################
run_one_replicate <- function(N, replicate_id) {
  dat <- simulate_multistate_data(
    N = N,
    pars = true_par,
    Tmax = Tmax,
    censor_rate = censor_rate
  )

  fit <- tryCatch(
    fit_conventional_model(dat$start_stop),
    error = function(e) {
      message(
        "Fit failed: N=", N,
        ", replicate=", replicate_id,
        ", error=", conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(data.frame(
      N = N,
      replicate = replicate_id,
      parameter = NA_character_,
      estimate = NA_real_, posterior_sd = NA_real_,
      lower = NA_real_, median = NA_real_, upper = NA_real_,
      rhat = NA_real_, ess = NA_real_,
      censor_before_first = unname(dat$diagnostics[1]),
      n12 = unname(dat$diagnostics[2]),
      n13 = unname(dat$diagnostics[3]),
      n23 = unname(dat$diagnostics[4]),
      converged = FALSE
    ))
  }

  fit$N <- N
  fit$replicate <- replicate_id
  fit$censor_before_first <- unname(dat$diagnostics[1])
  fit$n12 <- unname(dat$diagnostics[2])
  fit$n13 <- unname(dat$diagnostics[3])
  fit$n23 <- unname(dat$diagnostics[4])
  fit$converged <- is.finite(fit$rhat) & fit$rhat < 1.05 & fit$ess >= 200

  fit[, c(
    "N", "replicate", "parameter", "estimate", "posterior_sd",
    "lower", "median", "upper", "rhat", "ess",
    "censor_before_first", "n12", "n13", "n23", "converged"
  )]
}

############################################################
## 5. Run simulation
############################################################
all_results <- list()
result_counter <- 0L

for (N in sample_sizes) {
  message("Starting N = ", N)

  for (r in seq_len(R_rep)) {
    message("  replicate ", r, " / ", R_rep)

    result_counter <- result_counter + 1L
    all_results[[result_counter]] <- run_one_replicate(N, r)

    ## Save progress after every replicate.
    current_results <- do.call(rbind, all_results)
    saveRDS(
      current_results,
      file = file.path(output_dir, "raw_results_progress.rds")
    )
  }
}

raw_results <- do.call(rbind, all_results)
write.csv(
  raw_results,
  file = file.path(output_dir, "raw_posterior_results.csv"),
  row.names = FALSE
)
saveRDS(raw_results, file.path(output_dir, "raw_posterior_results.rds"))

############################################################
## 6. Performance summaries
############################################################
## Important interpretation:
## The fitted conventional model is misspecified relative to the conditional
## data-generating model because it omits b0 and b1. Therefore, comparisons
## with the conditional generating values quantify discrepancy from those
## generating values; they are not a pure measure of finite-sample estimator
## bias under a correctly specified model.
comparison_true <- c(
  lambda012 = true_par["lambda012"],
  lambda013 = true_par["lambda013"],
  lambda023 = true_par["lambda023"],
  gamma12   = true_par["gamma12"],
  gamma13   = true_par["gamma13"],
  gamma23   = true_par["gamma23"]
)
comparison_true <- as.numeric(comparison_true)
names(comparison_true) <- c(
  "lambda012", "lambda013", "lambda023",
  "gamma12", "gamma13", "gamma23"
)

summarise_one_parameter <- function(df) {
  p <- unique(df$parameter)
  tv <- comparison_true[p]

  valid <- is.finite(df$estimate) & is.finite(df$lower) & is.finite(df$upper)
  d <- df[valid, , drop = FALSE]

  if (nrow(d) == 0L) {
    return(data.frame(
      N = unique(df$N), parameter = p, true_value = tv,
      successful_replicates = 0L, mean_estimate = NA_real_,
      empirical_sd = NA_real_, mean_posterior_sd = NA_real_,
      bias = NA_real_, relative_bias = NA_real_, mse = NA_real_,
      rmse = NA_real_, coverage_95 = NA_real_, mcse_mean = NA_real_,
      mean_rhat = NA_real_, max_rhat = NA_real_, mean_ess = NA_real_,
      convergence_rate = NA_real_
    ))
  }

  bias_values <- d$estimate - tv

  data.frame(
    N = unique(df$N),
    parameter = p,
    true_value = tv,
    successful_replicates = nrow(d),
    mean_estimate = mean(d$estimate),
    empirical_sd = sd(d$estimate),
    mean_posterior_sd = mean(d$posterior_sd),
    bias = mean(bias_values),
    relative_bias = if (abs(tv) > .Machine$double.eps) {
      mean(bias_values) / tv
    } else {
      NA_real_
    },
    mse = mean(bias_values^2),
    rmse = sqrt(mean(bias_values^2)),
    coverage_95 = mean(d$lower <= tv & d$upper >= tv),
    mcse_mean = sd(d$estimate) / sqrt(nrow(d)),
    mean_rhat = mean(d$rhat, na.rm = TRUE),
    max_rhat = max(d$rhat, na.rm = TRUE),
    mean_ess = mean(d$ess, na.rm = TRUE),
    convergence_rate = mean(d$converged, na.rm = TRUE)
  )
}

valid_rows <- raw_results[!is.na(raw_results$parameter), , drop = FALSE]
split_results <- split(
  valid_rows,
  interaction(valid_rows$N, valid_rows$parameter, drop = TRUE)
)
performance_summary <- do.call(
  rbind,
  lapply(split_results, summarise_one_parameter)
)
rownames(performance_summary) <- NULL
performance_summary <- performance_summary[
  order(performance_summary$N, match(
    performance_summary$parameter,
    c("lambda012", "lambda013", "lambda023", "gamma12", "gamma13", "gamma23")
  )),
]

write.csv(
  performance_summary,
  file = file.path(output_dir, "performance_summary.csv"),
  row.names = FALSE
)

############################################################
## 7. Dataset-level transition and censoring summaries
############################################################
dataset_summary <- unique(raw_results[, c(
  "N", "replicate", "censor_before_first", "n12", "n13", "n23"
)])

event_summary <- do.call(
  rbind,
  lapply(split(dataset_summary, dataset_summary$N), function(d) {
    data.frame(
      N = unique(d$N),
      mean_censor_before_first = mean(d$censor_before_first),
      mean_n12 = mean(d$n12),
      mean_n13 = mean(d$n13),
      mean_n23 = mean(d$n23)
    )
  })
)
rownames(event_summary) <- NULL

write.csv(
  event_summary,
  file = file.path(output_dir, "event_censoring_summary.csv"),
  row.names = FALSE
)

############################################################
## 8. Session information
############################################################
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

message("Simulation completed. Results saved in: ", normalizePath(output_dir))
print(performance_summary)
print(event_summary)
