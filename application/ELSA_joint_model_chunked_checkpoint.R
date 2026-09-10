############################################################
# ELSA FOUR-STATE BAYESIAN JOINT MODEL
#
# Longitudinal outcome:
#   cflisd = delayed recall
#
# LONGITUDINAL FIXED EFFECTS:
#   time
#   age
#   cardiovascular condition
#   high blood pressure
#   diabetes
#   stroke
#   cancer
#   arthritis
#   visual acuity / eyesight
#
# Longitudinal random effects:
#   random intercept + random slope for time
#
# MULTISTATE PROCESS:
#   0 = Non-frail
#   1 = Pre-frail
#   2 = Frail
#   3 = Death (absorbing)
#
# All 9 existing transitions are retained:
#   0->1, 0->2, 0->3,
#   1->0, 1->2, 1->3,
#   2->0, 2->1, 2->3
#
# MULTISTATE BASELINE COVARIATES (8):
#   age
#   cardiovascular condition
#   high blood pressure
#   diabetes
#   stroke
#   cancer
#   arthritis
#   visual acuity / eyesight
#
# IMPORTANT:
#   baseline delayed recall is REMOVED from the multistate X.
#   Repeated cflisd is modelled longitudinally and linked to the
#   9 transition intensities through shared random intercept/slope:
#
#   log lambda_ik =
#       log lambda0_k
#       + gamma_k' X_i
#       + alpha0_k*b0_i
#       + alpha1_k*b1_i
#
# The existing multistate exposure/count likelihood, 4 states,
# 9 transitions, zero trick, lambda0 priors and gamma priors
# are otherwise unchanged.
############################################################


############################################################
# 0. SETTINGS
############################################################

LONG_FILE <- "17_ELSA_JOINT_LONGITUDINAL.csv"
MS_FILE   <- "19_ELSA_JOINT_MULTISTATE_TRANSITIONS.csv"

# ==========================================================
# RUN MODE
# ==========================================================
# "TEST" = local smoke test on a small participant subset.
# "FULL" = full analysis.
RUN_MODE <- "FULL"

OUT_DIR_BASE <- "ELSA_JOINT_9TRANSITION_8COV_RESULTS"
OUT_DIR <- if (RUN_MODE == "TEST") {
  paste0(OUT_DIR_BASE, "_TEST")
} else {
  OUT_DIR_BASE
}

MODEL_FILE <- file.path(
  OUT_DIR,
  "elsa_joint_9transition_8cov_shared_RE.jags"
)

CHECKPOINT_DIR <- file.path(
  OUT_DIR,
  "MCMC_CHECKPOINTS"
)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  CHECKPOINT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

SEED <- 20260903
set.seed(SEED)

# ==========================================================
# TEST SETTINGS
# ==========================================================
# The test mode is ONLY for checking that the whole workflow:
# compile -> burn-in -> sample -> save -> summarize
# works correctly. Do not use TEST estimates in the thesis.
TEST_N_SUBJECTS <- 200L
TEST_N_CHAINS   <- 2L
TEST_N_ADAPT    <- 100L
TEST_N_BURN     <- 200L
TEST_N_ITER     <- 400L
TEST_N_THIN     <- 1L
TEST_CHUNK_ITER <- 100L

# ==========================================================
# FULL SETTINGS
# ==========================================================
FULL_N_CHAINS   <- 2L
FULL_N_ADAPT    <- 1000L
FULL_N_BURN     <- 2000L
FULL_N_ITER     <- 20000L
FULL_N_THIN     <- 5L

# Posterior sampling is saved every FULL_CHUNK_ITER raw
# iterations. If the run is interrupted, at most the current
# unfinished chunk is lost.
FULL_CHUNK_ITER <- 250L

# Optional one-time extension after the initial target.
AUTO_EXTEND <- FALSE
EXTEND_BURN <- 0L
EXTEND_ITER <- 20000L
EXTEND_THIN <- 10L
EXTEND_CHUNK_ITER <- 1000L

# DIC can be very expensive. Keep FALSE while testing/running.
# Change to TRUE only after you have a satisfactory final MCMC.
RUN_DIC <- FALSE

if (RUN_MODE == "TEST") {
  N_CHAINS <- TEST_N_CHAINS
  N_ADAPT  <- TEST_N_ADAPT
  N_BURN   <- TEST_N_BURN
  N_ITER   <- TEST_N_ITER
  N_THIN   <- TEST_N_THIN
  CHUNK_ITER <- TEST_CHUNK_ITER
  AUTO_EXTEND <- FALSE
  RUN_DIC <- FALSE
} else {
  N_CHAINS <- FULL_N_CHAINS
  N_ADAPT  <- FULL_N_ADAPT
  N_BURN   <- FULL_N_BURN
  N_ITER   <- FULL_N_ITER
  N_THIN   <- FULL_N_THIN
  CHUNK_ITER <- FULL_CHUNK_ITER
}

if (N_ITER %% CHUNK_ITER != 0) {
  stop("N_ITER must be an exact multiple of CHUNK_ITER.")
}

if (AUTO_EXTEND && EXTEND_ITER %% EXTEND_CHUNK_ITER != 0) {
  stop("EXTEND_ITER must be an exact multiple of EXTEND_CHUNK_ITER.")
}

cat(
  "\n================ RUN SETTINGS ================\n",
  "RUN_MODE: ", RUN_MODE, "\n",
  "OUT_DIR: ", OUT_DIR, "\n",
  "N_CHAINS: ", N_CHAINS, "\n",
  "N_ADAPT: ", N_ADAPT, "\n",
  "N_BURN: ", N_BURN, "\n",
  "N_ITER: ", N_ITER, "\n",
  "N_THIN: ", N_THIN, "\n",
  "CHUNK_ITER: ", CHUNK_ITER, "\n",
  "AUTO_EXTEND: ", AUTO_EXTEND, "\n",
  "RUN_DIC: ", RUN_DIC, "\n",
  "================================================\n"
)


############################################################
# 1. PACKAGES
############################################################

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "stringr",
  "rjags",
  "coda",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(rjags)
library(coda)
library(ggplot2)


############################################################
# 2. READ DATA
############################################################

long_raw <- read_csv(
  LONG_FILE,
  show_col_types = FALSE
)

ms_raw <- read_csv(
  MS_FILE,
  show_col_types = FALSE
)

cat(
  "\nLongitudinal data dimensions:",
  nrow(long_raw),
  "x",
  ncol(long_raw),
  "\n"
)

cat(
  "Multistate data dimensions:",
  nrow(ms_raw),
  "x",
  ncol(ms_raw),
  "\n"
)


############################################################
# 2A. OPTIONAL SMALL-SAMPLE TEST MODE
############################################################

if (RUN_MODE == "TEST") {

  long_ids <- unique(
    as.character(long_raw$idauniq)
  )

  ms_ids <- unique(
    as.character(ms_raw$idauniq)
  )

  eligible_test_ids <- intersect(
    long_ids,
    ms_ids
  )

  if (length(eligible_test_ids) == 0) {
    stop("No common participant IDs were found for TEST mode.")
  }

  set.seed(SEED)

  n_test <- min(
    TEST_N_SUBJECTS,
    length(eligible_test_ids)
  )

  test_ids <- sample(
    eligible_test_ids,
    size = n_test,
    replace = FALSE
  )

  long_raw <- long_raw %>%
    dplyr::filter(
      as.character(idauniq) %in% test_ids
    )

  ms_raw <- ms_raw %>%
    dplyr::filter(
      as.character(idauniq) %in% test_ids
    )

  cat(
    "\n================ TEST SUBSET ================\n",
    "Participants requested: ", TEST_N_SUBJECTS, "\n",
    "Participants sampled: ", n_test, "\n",
    "Longitudinal rows: ", nrow(long_raw), "\n",
    "Multistate rows: ", nrow(ms_raw), "\n",
    "=============================================\n"
  )
}


############################################################
# 3. REQUIRED VARIABLE CHECKS
############################################################

required_long <- c(
  "idauniq",
  "year",
  "age_final",
  "hediman",
  "hedimmi",
  "hedimhf",
  "hedimar",
  "hedimbp",
  "hedimdi",
  "hedimst",
  "hedibca",
  "hedibar",
  "heeye",
  "cflisd"
)

required_ms <- c(
  "idauniq",
  "from",
  "to",
  "start",
  "stop",
  "calendar_start",
  "calendar_stop",
  "interval"
)

missing_long <- setdiff(
  required_long,
  names(long_raw)
)

missing_ms <- setdiff(
  required_ms,
  names(ms_raw)
)

if (length(missing_long) > 0) {
  stop(
    "Longitudinal file is missing required variable(s): ",
    paste(missing_long, collapse = ", ")
  )
}

if (length(missing_ms) > 0) {
  stop(
    "Multistate file is missing required variable(s): ",
    paste(missing_ms, collapse = ", ")
  )
}


############################################################
# 4. CLEAN AND VALIDATE MULTISTATE DATA
#    EXISTING SETTINGS RETAINED
############################################################

ms <- ms_raw %>%
  dplyr::transmute(
    idauniq = as.character(idauniq),
    from = as.integer(from),
    to = as.integer(to),
    start = as.numeric(start),
    stop = as.numeric(stop),
    calendar_start = as.numeric(calendar_start),
    calendar_stop = as.numeric(calendar_stop),
    interval = as.numeric(interval),
    transition_source =
      if ("transition_source" %in% names(ms_raw)) {
        as.character(transition_source)
      } else {
        NA_character_
      },
    transition_type =
      if ("transition_type" %in% names(ms_raw)) {
        as.character(transition_type)
      } else {
        NA_character_
      }
  )

if (
  anyNA(
    ms[
      ,
      c(
        "idauniq",
        "from",
        "to",
        "start",
        "stop",
        "interval"
      )
    ]
  )
) {
  stop(
    "Missing ID/state/time values were found in the multistate file."
  )
}

if (
  any(
    ms$from < 0 |
      ms$from > 3 |
      ms$to < 0 |
      ms$to > 3
  )
) {
  stop(
    "State values outside {0,1,2,3} were found."
  )
}

if (any(ms$from == 3)) {
  stop(
    "Outgoing transitions from State 3 were found, but State 3 must be absorbing."
  )
}

if (any(ms$interval <= 0)) {
  stop(
    "All multistate intervals must be > 0."
  )
}

interval_difference <- abs(
  ms$interval -
    (ms$stop - ms$start)
)

if (any(interval_difference > 1e-8)) {
  warning(
    sum(interval_difference > 1e-8),
    " row(s) have interval != stop-start. ",
    "The supplied 'interval' variable is used."
  )
}

allowed_pairs <- tribble(
  ~from, ~to,
  0L, 0L,
  0L, 1L,
  0L, 2L,
  0L, 3L,
  1L, 0L,
  1L, 1L,
  1L, 2L,
  1L, 3L,
  2L, 0L,
  2L, 1L,
  2L, 2L,
  2L, 3L
)

invalid_pairs <- ms %>%
  dplyr::distinct(
    from,
    to
  ) %>%
  dplyr::anti_join(
    allowed_pairs,
    by = c(
      "from",
      "to"
    )
  )

if (nrow(invalid_pairs) > 0) {
  print(invalid_pairs)
  stop(
    "Observed transition(s) fall outside the four-state model."
  )
}


############################################################
# 5. OBSERVED TRANSITION AUDIT
############################################################

transition_audit <- ms %>%
  dplyr::count(
    from,
    to,
    name = "n"
  ) %>%
  dplyr::arrange(
    from,
    to
  )

cat(
  "\n================ OBSERVED TRANSITIONS ================\n"
)

print(
  transition_audit,
  n = Inf
)

write_csv(
  transition_audit,
  file.path(
    OUT_DIR,
    "00_observed_transition_counts.csv"
  )
)

n_rows_before <- nrow(ms)


############################################################
# 6. DEFINE 9 TRANSITIONS
#    UNCHANGED
############################################################

transition_map <- tribble(
  ~k, ~from, ~to, ~transition,
  1L, 0L, 1L, "0->1",
  2L, 0L, 2L, "0->2",
  3L, 0L, 3L, "0->3",
  4L, 1L, 0L, "1->0",
  5L, 1L, 2L, "1->2",
  6L, 1L, 3L, "1->3",
  7L, 2L, 0L, "2->0",
  8L, 2L, 1L, "2->1",
  9L, 2L, 3L, "2->3"
)


############################################################
# 7. SUBJECT-SPECIFIC MULTISTATE BASELINE YEAR
############################################################

ms_baseline <- ms %>%
  dplyr::group_by(
    idauniq
  ) %>%
  dplyr::summarise(
    ms_baseline_year =
      min(
        calendar_start,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


############################################################
# 8. PREPARE SOURCE DATA FOR THE 8 BASELINE COVARIATES
############################################################

long_cov <- long_raw %>%
  dplyr::transmute(
    idauniq = as.character(idauniq),
    year = as.numeric(year),
    age_final = as.numeric(age_final),

    # cardiovascular components
    hediman = as.numeric(hediman),
    hedimmi = as.numeric(hedimmi),
    hedimhf = as.numeric(hedimhf),
    hedimar = as.numeric(hedimar),

    # remaining health covariates
    hedimbp = as.numeric(hedimbp),
    hedimdi = as.numeric(hedimdi),
    hedimst = as.numeric(hedimst),
    hedibca = as.numeric(hedibca),
    hedibar = as.numeric(hedibar),
    heeye = as.numeric(heeye)
  ) %>%
  dplyr::mutate(
    cardio_any = pmax(
      hediman,
      hedimmi,
      hedimhf,
      hedimar,
      na.rm = TRUE
    ),
    cardio_any = ifelse(
      is.infinite(cardio_any),
      NA_real_,
      cardio_any
    )
  ) %>%
  dplyr::filter(
    idauniq %in%
      ms_baseline$idauniq
  )


############################################################
# 9. EXTRACT SUBJECT-LEVEL BASELINE COVARIATES
#
# Same rule as existing multistate code:
# 1. most recent nonmissing value at/before MS baseline
# 2. if unavailable, nearest subsequent nonmissing value
############################################################

extract_baseline_var <- function(
  long_df,
  baseline_df,
  varname
) {

  tmp <- data.frame(
    idauniq = long_df[["idauniq"]],
    year = long_df[["year"]],
    value = long_df[[varname]],
    stringsAsFactors = FALSE
  )

  tmp <- tmp[
    !is.na(tmp$value) &
      !is.na(tmp$year),
    ,
    drop = FALSE
  ]

  base_tmp <- data.frame(
    idauniq = baseline_df[["idauniq"]],
    ms_baseline_year = baseline_df[["ms_baseline_year"]],
    stringsAsFactors = FALSE
  )

  tmp <- merge(
    tmp,
    base_tmp,
    by = "idauniq",
    all = FALSE,
    sort = FALSE
  )

  tmp$side <- ifelse(
    tmp$year <=
      tmp$ms_baseline_year,
    0L,
    1L
  )

  tmp$distance <- abs(
    tmp$year -
      tmp$ms_baseline_year
  )

  tmp <- tmp[
    order(
      tmp$idauniq,
      tmp$side,
      tmp$distance,
      -tmp$year
    ),
    ,
    drop = FALSE
  ]

  tmp <- tmp[
    !duplicated(tmp$idauniq),
    ,
    drop = FALSE
  ]

  tmp$source_relation <- ifelse(
    tmp$year ==
      tmp$ms_baseline_year,
    "exact_baseline_year",
    ifelse(
      tmp$year <
        tmp$ms_baseline_year,
      "previous_available_year",
      "subsequent_available_year"
    )
  )

  result <- data.frame(
    idauniq = tmp$idauniq,
    value = tmp$value,
    source_year = tmp$year,
    source_relation = tmp$source_relation,
    stringsAsFactors = FALSE
  )

  names(result)[
    names(result) ==
      "value"
  ] <- paste0(
    "base_",
    varname
  )

  names(result)[
    names(result) ==
      "source_year"
  ] <- paste0(
    "source_year_",
    varname
  )

  names(result)[
    names(result) ==
      "source_relation"
  ] <- paste0(
    "source_relation_",
    varname
  )

  result
}


base_age <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "age_final"
)

base_cardio <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "cardio_any"
)

base_bp <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "hedimbp"
)

base_dm <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "hedimdi"
)

base_stroke <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "hedimst"
)

base_cancer <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "hedibca"
)

base_arth <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "hedibar"
)

base_eye <- extract_baseline_var(
  long_cov,
  ms_baseline,
  "heeye"
)


baseline_cov <- ms_baseline %>%
  dplyr::left_join(
    base_age,
    by = "idauniq"
  ) %>%
  dplyr::left_join(
    base_cardio,
    by = "idauniq"
  ) %>%
  dplyr::left_join(
    base_bp,
    by = "idauniq"
  ) %>%
  dplyr::left_join(
    base_dm,
    by = "idauniq"
  ) %>%
  dplyr::left_join(
    base_stroke,
    by = "idauniq"
  ) %>%
  dplyr::left_join(
    base_cancer,
    by = "idauniq"
  ) %>%
  dplyr::left_join(
    base_arth,
    by = "idauniq"
  ) %>%
  dplyr::left_join(
    base_eye,
    by = "idauniq"
  )


mode_numeric <- function(x) {

  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0) {
    return(NA_real_)
  }

  ux <- unique(x)

  ux[
    which.max(
      tabulate(
        match(
          x,
          ux
        )
      )
    )
  ]
}


fallback_values <- list(

  age_final =
    median(
      baseline_cov$base_age_final,
      na.rm = TRUE
    ),

  cardio_any =
    mode_numeric(
      baseline_cov$base_cardio_any
    ),

  hedimbp =
    mode_numeric(
      baseline_cov$base_hedimbp
    ),

  hedimdi =
    mode_numeric(
      baseline_cov$base_hedimdi
    ),

  hedimst =
    mode_numeric(
      baseline_cov$base_hedimst
    ),

  hedibca =
    mode_numeric(
      baseline_cov$base_hedibca
    ),

  hedibar =
    mode_numeric(
      baseline_cov$base_hedibar
    ),

  heeye =
    median(
      baseline_cov$base_heeye,
      na.rm = TRUE
    )
)


baseline_cov <- baseline_cov %>%
  dplyr::mutate(

    imputed_age_final =
      is.na(base_age_final),

    imputed_cardio_any =
      is.na(base_cardio_any),

    imputed_hedimbp =
      is.na(base_hedimbp),

    imputed_hedimdi =
      is.na(base_hedimdi),

    imputed_hedimst =
      is.na(base_hedimst),

    imputed_hedibca =
      is.na(base_hedibca),

    imputed_hedibar =
      is.na(base_hedibar),

    imputed_heeye =
      is.na(base_heeye),

    base_age_final =
      replace_na(
        base_age_final,
        fallback_values$age_final
      ),

    base_cardio_any =
      replace_na(
        base_cardio_any,
        fallback_values$cardio_any
      ),

    base_hedimbp =
      replace_na(
        base_hedimbp,
        fallback_values$hedimbp
      ),

    base_hedimdi =
      replace_na(
        base_hedimdi,
        fallback_values$hedimdi
      ),

    base_hedimst =
      replace_na(
        base_hedimst,
        fallback_values$hedimst
      ),

    base_hedibca =
      replace_na(
        base_hedibca,
        fallback_values$hedibca
      ),

    base_hedibar =
      replace_na(
        base_hedibar,
        fallback_values$hedibar
      ),

    base_heeye =
      replace_na(
        base_heeye,
        fallback_values$heeye
      )
  )


binary_baseline_vars <- c(
  "base_cardio_any",
  "base_hedimbp",
  "base_hedimdi",
  "base_hedimst",
  "base_hedibca",
  "base_hedibar"
)

for (v in binary_baseline_vars) {

  if (
    any(
      !baseline_cov[[v]] %in%
        c(
          0,
          1
        )
    )
  ) {
    stop(
      "Baseline ",
      v,
      " contains values other than 0/1."
    )
  }
}


cat(
  "\n================ BASELINE COVARIATE EXTRACTION ================\n"
)

print(
  baseline_cov %>%
    dplyr::summarise(
      age_final =
        sum(imputed_age_final),

      cardio_any =
        sum(imputed_cardio_any),

      hedimbp =
        sum(imputed_hedimbp),

      hedimdi =
        sum(imputed_hedimdi),

      hedimst =
        sum(imputed_hedimst),

      hedibca =
        sum(imputed_hedibca),

      hedibar =
        sum(imputed_hedibar),

      heeye =
        sum(imputed_heeye)
    )
)


write_csv(
  baseline_cov,
  file.path(
    OUT_DIR,
    "01_subject_baseline_covariates_8cov.csv"
  )
)


############################################################
# 10. AGGREGATE MULTISTATE EXPOSURE AND REPEATED COUNTS
#     UNCHANGED
############################################################

ms_subject <- ms %>%
  dplyr::group_by(
    idauniq
  ) %>%
  dplyr::summarise(

    time0 =
      sum(
        interval[from == 0],
        na.rm = TRUE
      ),

    time1 =
      sum(
        interval[from == 1],
        na.rm = TRUE
      ),

    time2 =
      sum(
        interval[from == 2],
        na.rm = TRUE
      ),

    d01 =
      sum(from == 0 & to == 1),

    d02 =
      sum(from == 0 & to == 2),

    d03 =
      sum(from == 0 & to == 3),

    d10 =
      sum(from == 1 & to == 0),

    d12 =
      sum(from == 1 & to == 2),

    d13 =
      sum(from == 1 & to == 3),

    d20 =
      sum(from == 2 & to == 0),

    d21 =
      sum(from == 2 & to == 1),

    d23 =
      sum(from == 2 & to == 3),

    n_intervals = n(),

    .groups = "drop"
  ) %>%
  dplyr::left_join(
    baseline_cov,
    by = "idauniq"
  )


if (
  nrow(ms_subject) !=
    n_distinct(ms$idauniq)
) {
  stop(
    "Subject aggregation failed: participant count changed."
  )
}


required_baseline_model_vars <- c(
  "base_age_final",
  "base_cardio_any",
  "base_hedimbp",
  "base_hedimdi",
  "base_hedimst",
  "base_hedibca",
  "base_hedibar",
  "base_heeye"
)


if (
  anyNA(
    ms_subject[
      ,
      required_baseline_model_vars
    ]
  )
) {
  stop(
    "Missing baseline covariates remain after extraction/fallback."
  )
}


############################################################
# 11. VERIFY NO MULTISTATE TRANSITIONS WERE LOST
############################################################

aggregate_check <- tibble(

  transition =
    c(
      "0->1",
      "0->2",
      "0->3",
      "1->0",
      "1->2",
      "1->3",
      "2->0",
      "2->1",
      "2->3"
    ),

  original_n =
    c(
      sum(ms$from == 0 & ms$to == 1),
      sum(ms$from == 0 & ms$to == 2),
      sum(ms$from == 0 & ms$to == 3),
      sum(ms$from == 1 & ms$to == 0),
      sum(ms$from == 1 & ms$to == 2),
      sum(ms$from == 1 & ms$to == 3),
      sum(ms$from == 2 & ms$to == 0),
      sum(ms$from == 2 & ms$to == 1),
      sum(ms$from == 2 & ms$to == 3)
    ),

  model_n =
    c(
      sum(ms_subject$d01),
      sum(ms_subject$d02),
      sum(ms_subject$d03),
      sum(ms_subject$d10),
      sum(ms_subject$d12),
      sum(ms_subject$d13),
      sum(ms_subject$d20),
      sum(ms_subject$d21),
      sum(ms_subject$d23)
    )
) %>%
  dplyr::mutate(
    retained =
      original_n ==
        model_n
  )


cat(
  "\n================ TRANSITION RETENTION CHECK ================\n"
)

print(
  aggregate_check,
  n = Inf
)


if (!all(aggregate_check$retained)) {
  stop(
    "At least one observed transition was lost during aggregation."
  )
}


if (nrow(ms) != n_rows_before) {
  stop(
    "Multistate interval rows were unexpectedly dropped."
  )
}


write_csv(
  aggregate_check,
  file.path(
    OUT_DIR,
    "02_transition_retention_check.csv"
  )
)


############################################################
# 12. SCALE/CODE THE 8 COVARIATES
#
# Same coding for longitudinal and multistate submodels.
#
# age:    per 10-year increase
# heeye:  per 0.2-unit increase
# binary health variables: 0/1
############################################################

AGE_CENTER <- mean(
  ms_subject$base_age_final
)

EYE_CENTER <- mean(
  ms_subject$base_heeye
)


ms_subject <- ms_subject %>%
  dplyr::mutate(

    x_age10 =
      (
        base_age_final -
          AGE_CENTER
      ) / 10,

    x_cardio =
      base_cardio_any,

    x_hedimbp =
      base_hedimbp,

    x_hedimdi =
      base_hedimdi,

    x_hedimst =
      base_hedimst,

    x_hedibca =
      base_hedibca,

    x_hedibar =
      base_hedibar,

    x_heeye02 =
      (
        base_heeye -
          EYE_CENTER
      ) / 0.2
  )


covariate_columns <- c(
  "x_age10",
  "x_cardio",
  "x_hedimbp",
  "x_hedimdi",
  "x_hedimst",
  "x_hedibca",
  "x_hedibar",
  "x_heeye02"
)


X_ms <- as.matrix(
  ms_subject[
    ,
    covariate_columns
  ]
)

storage.mode(X_ms) <- "double"


# The longitudinal model uses the SAME 8 subject-level
# baseline covariates.
X_long <- X_ms


covariate_map <- tibble(

  p = 1:8,

  covariate =
    c(
      "Age (per 10-year increase)",
      "Cardiovascular condition (any; 1 vs 0)",
      "High blood pressure (1 vs 0)",
      "Diabetes (1 vs 0)",
      "Stroke (1 vs 0)",
      "Cancer (1 vs 0)",
      "Arthritis (1 vs 0)",
      "Eyesight score (per 0.2-unit increase)"
    )
)


write_csv(
  tibble(
    age_center_years = AGE_CENTER,
    heeye_center = EYE_CENTER
  ),
  file.path(
    OUT_DIR,
    "03_covariate_centering_constants.csv"
  )
)


############################################################
# 13. SUBJECT INDEX
#     SAME INDEX USED BY BOTH SUBMODELS
############################################################

subject_index <- ms_subject %>%
  dplyr::transmute(
    idauniq,
    subject_id =
      dplyr::row_number()
  )


############################################################
# 14. PREPARE LONGITUDINAL DATA
#
# Outcome:
#   repeated cflisd
#
# Time origin:
#   subject-specific multistate baseline year
#
# Covariates:
#   supplied through subject-level X_long
#
# No multistate rows/transitions are removed here.
############################################################

longitudinal_data <- long_raw %>%
  dplyr::transmute(
    idauniq =
      as.character(idauniq),

    year =
      as.numeric(year),

    cflisd =
      as.numeric(cflisd)
  ) %>%
  dplyr::inner_join(
    ms_baseline,
    by = "idauniq"
  ) %>%
  dplyr::inner_join(
    subject_index,
    by = "idauniq"
  ) %>%
  dplyr::filter(
    !is.na(year),
    !is.na(cflisd),
    year >=
      ms_baseline_year
  ) %>%
  dplyr::mutate(
    obs_time =
      year -
        ms_baseline_year
  ) %>%
  dplyr::arrange(
    subject_id,
    obs_time,
    year
  )


if (nrow(longitudinal_data) == 0) {
  stop(
    "No usable longitudinal cflisd observations remain after alignment."
  )
}


if (any(longitudinal_data$obs_time < 0)) {
  stop(
    "Negative longitudinal follow-up times were found."
  )
}


longitudinal_subject_audit <- subject_index %>%
  dplyr::left_join(
    longitudinal_data %>%
      dplyr::count(
        idauniq,
        name = "n_long_obs"
      ),
    by = "idauniq"
  ) %>%
  dplyr::mutate(
    n_long_obs =
      replace_na(
        n_long_obs,
        0L
      )
  )


cat(
  "\n================ LONGITUDINAL DATA ================\n"
)

cat(
  "Longitudinal observations:",
  nrow(longitudinal_data),
  "\n"
)

cat(
  "Subjects with >=1 longitudinal observation:",
  sum(
    longitudinal_subject_audit$n_long_obs >= 1
  ),
  "of",
  nrow(subject_index),
  "\n"
)

cat(
  "Subjects with >=2 longitudinal observations:",
  sum(
    longitudinal_subject_audit$n_long_obs >= 2
  ),
  "\n"
)


write_csv(
  longitudinal_data,
  file.path(
    OUT_DIR,
    "04_longitudinal_model_data.csv"
  )
)


write_csv(
  longitudinal_subject_audit,
  file.path(
    OUT_DIR,
    "05_longitudinal_subject_audit.csv"
  )
)


############################################################
# 15. CRUDE MULTISTATE RATES
#     EXISTING MULTISTATE SETTING RETAINED
############################################################

exposure_by_source <- c(

  `0` =
    sum(ms_subject$time0),

  `1` =
    sum(ms_subject$time1),

  `2` =
    sum(ms_subject$time2)
)


observed_events <- aggregate_check$original_n
source_for_k <- transition_map$from


crude_rates <- vapply(

  seq_len(
    nrow(transition_map)
  ),

  function(k) {

    observed_events[k] /
      exposure_by_source[
        as.character(
          source_for_k[k]
        )
      ]
  },

  numeric(1)
)


crude_rate_table <- transition_map %>%
  dplyr::mutate(

    events =
      observed_events,

    source_exposure_years =
      exposure_by_source[
        as.character(from)
      ],

    crude_rate_per_person_year =
      crude_rates,

    sparse_transition =
      events < 100
  )


write_csv(
  crude_rate_table,
  file.path(
    OUT_DIR,
    "06_crude_transition_rates.csv"
  )
)


############################################################
# 16. JOINT JAGS MODEL
############################################################

jags_model_string <- "
model {

  ############################################################
  # A. LONGITUDINAL SUBMODEL
  ############################################################

  for (j in 1:Nobs) {

    mu_y[j] <-
      beta0_y +
      beta_time_y * obs_time[j] +
      inprod(beta_long[1:P], X_long[obs_id[j],1:P]) +
      b[obs_id[j],1] +
      b[obs_id[j],2] * obs_time[j]

    y[j] ~ dnorm(
      mu_y[j],
      tau_y
    )
  }

  ############################################################
  # B. SHARED RANDOM INTERCEPT + SLOPE
  ############################################################

  for (i in 1:N) {

    b[i,1:2] ~ dmnorm(
      zero_b[],
      Tau_b[,]
    )
  }

  ############################################################
  # C. MULTISTATE SUBMODEL
  #    4 states / all 9 transitions retained
  ############################################################

  for (i in 1:N) {

    for (k in 1:K) {

      log_lambda[i,k] <-
        log_lambda0[k] +
        inprod(
          gamma[k,1:P],
          X_ms[i,1:P]
        ) +
        alpha0[k] * b[i,1] +
        alpha1[k] * b[i,2]

      lambda[i,k] <-
        exp(
          log_lambda[i,k]
        )
    }

    loglik_ms[i] <-

      d01[i] * log_lambda[i,1] +
      d02[i] * log_lambda[i,2] +
      d03[i] * log_lambda[i,3] -
      time0[i] * (
        lambda[i,1] +
        lambda[i,2] +
        lambda[i,3]
      ) +

      d10[i] * log_lambda[i,4] +
      d12[i] * log_lambda[i,5] +
      d13[i] * log_lambda[i,6] -
      time1[i] * (
        lambda[i,4] +
        lambda[i,5] +
        lambda[i,6]
      ) +

      d20[i] * log_lambda[i,7] +
      d21[i] * log_lambda[i,8] +
      d23[i] * log_lambda[i,9] -
      time2[i] * (
        lambda[i,7] +
        lambda[i,8] +
        lambda[i,9]
      )

    phi[i] <-
      -loglik_ms[i] +
      Czero

    zeros[i] ~ dpois(
      phi[i]
    )
  }

  ############################################################
  # D. LONGITUDINAL PRIORS
  ############################################################

  beta0_y ~ dnorm(
    0,
    0.0001
  )

  beta_time_y ~ dnorm(
    0,
    0.25
  )

  for (p in 1:P) {

    beta_long[p] ~ dnorm(
      0,
      0.0001
    )
  }

  tau_y ~ dgamma(
    0.001,
    0.001
  )

  sigma2_y <-
    1 /
    tau_y

  sigma_y <-
    sqrt(
      sigma2_y
    )

  ############################################################
  # E. RANDOM-EFFECT COVARIANCE
  ############################################################

  Tau_b[1:2,1:2] ~ dwish(
    S0[,],
    nu0
  )

  D[1:2,1:2] <-
    inverse(
      Tau_b[,]
    )

  sigma_b0 <-
    sqrt(
      D[1,1]
    )

  sigma_b1 <-
    sqrt(
      D[2,2]
    )

  rho <-
    D[1,2] /
    sqrt(
      D[1,1] *
      D[2,2]
    )

  ############################################################
  # F. BASELINE TRANSITION INTENSITIES
  #    SAME AS CURRENT MULTISTATE MODEL
  ############################################################

  for (k in 1:K) {

    lambda0[k] ~ dlnorm(
      0,
      1
    )

    log_lambda0[k] <-
      log(
        lambda0[k]
      )
  }

  ############################################################
  # G. MULTISTATE COVARIATE EFFECTS
  #    SAME PRIOR AS CURRENT MULTISTATE MODEL
  ############################################################

  for (k in 1:K) {

    for (p in 1:P) {

      gamma[k,p] ~ dnorm(
        0,
        0.25
      )
    }
  }

  ############################################################
  # H. JOINT ASSOCIATION PARAMETERS
  #    SAME PRIOR FAMILY AS SIMULATION MODEL
  ############################################################

  for (k in 1:K) {

    alpha0[k] ~ dnorm(
      0,
      0.01
    )

    alpha1[k] ~ dnorm(
      0,
      0.01
    )
  }
}
"


writeLines(
  jags_model_string,
  MODEL_FILE
)

cat(
  "\nJAGS model written to:\n",
  MODEL_FILE,
  "\n"
)


############################################################
# 17. JAGS DATA
############################################################

S0 <- diag(2)
nu0 <- 2


jags_data <- list(

  N =
    nrow(ms_subject),

  Nobs =
    nrow(longitudinal_data),

  K = 9L,

  # IMPORTANT:
  # delayed recall was removed from multistate covariates,
  # therefore P is now 8, not 9.
  P = 8L,

  ##########################################################
  # longitudinal
  ##########################################################

  y =
    as.numeric(
      longitudinal_data$cflisd
    ),

  obs_time =
    as.numeric(
      longitudinal_data$obs_time
    ),

  obs_id =
    as.integer(
      longitudinal_data$subject_id
    ),

  X_long =
    X_long,

  ##########################################################
  # multistate
  ##########################################################

  X_ms =
    X_ms,

  time0 =
    as.numeric(
      ms_subject$time0
    ),

  time1 =
    as.numeric(
      ms_subject$time1
    ),

  time2 =
    as.numeric(
      ms_subject$time2
    ),

  d01 =
    as.integer(
      ms_subject$d01
    ),

  d02 =
    as.integer(
      ms_subject$d02
    ),

  d03 =
    as.integer(
      ms_subject$d03
    ),

  d10 =
    as.integer(
      ms_subject$d10
    ),

  d12 =
    as.integer(
      ms_subject$d12
    ),

  d13 =
    as.integer(
      ms_subject$d13
    ),

  d20 =
    as.integer(
      ms_subject$d20
    ),

  d21 =
    as.integer(
      ms_subject$d21
    ),

  d23 =
    as.integer(
      ms_subject$d23
    ),

  zeros =
    rep(
      0L,
      nrow(ms_subject)
    ),

  Czero = 10000,

  zero_b =
    c(
      0,
      0
    ),

  S0 =
    S0,

  nu0 =
    nu0
)


############################################################
# 18. INITIAL VALUES
############################################################

init_rate <- pmax(
  crude_rates,
  1e-5
)


lm_start <- try(

  lm(
    cflisd ~ obs_time,
    data =
      longitudinal_data
  ),

  silent = TRUE
)


if (
  inherits(
    lm_start,
    "try-error"
  )
) {

  beta0_start <- mean(
    longitudinal_data$cflisd,
    na.rm = TRUE
  )

  beta_time_start <- 0

} else {

  beta0_start <- unname(
    coef(lm_start)[1]
  )

  beta_time_start <- unname(
    coef(lm_start)[2]
  )

  if (
    !is.finite(
      beta_time_start
    )
  ) {
    beta_time_start <- 0
  }
}


make_inits <- function(
  chain_id
) {

  rng_names <- c(
    "base::Wichmann-Hill",
    "base::Marsaglia-Multicarry",
    "base::Super-Duper",
    "base::Mersenne-Twister"
  )

  list(

    beta0_y =
      beta0_start +
        rnorm(
          1,
          0,
          0.05
        ),

    beta_time_y =
      beta_time_start +
        rnorm(
          1,
          0,
          0.02
        ),

    beta_long =
      rnorm(
        8,
        0,
        0.05
      ),

    tau_y = 1,

    Tau_b =
      diag(2),

    lambda0 =
      init_rate *
        exp(
          rnorm(
            9,
            0,
            0.10
          )
        ),

    gamma =
      matrix(
        rnorm(
          9 * 8,
          0,
          0.05
        ),
        nrow = 9,
        ncol = 8
      ),

    alpha0 =
      rnorm(
        9,
        0,
        0.05
      ),

    alpha1 =
      rnorm(
        9,
        0,
        0.05
      ),

    .RNG.name =
      rng_names[
        chain_id
      ],

    .RNG.seed =
      SEED +
        chain_id *
        100
  )
}


inits <- lapply(
  seq_len(N_CHAINS),
  make_inits
)


############################################################
# 19. COMPILE AND RUN JAGS
#     CHUNKED SAMPLING + CHECKPOINTS + RESTART SUPPORT
############################################################

parameters_to_monitor <- c(

  # longitudinal
  "beta0_y",
  "beta_time_y",
  "beta_long",
  "sigma_y",

  # random effects distribution
  "sigma_b0",
  "sigma_b1",
  "rho",

  # multistate
  "lambda0",
  "gamma",

  # joint association
  "alpha0",
  "alpha1"
)


# ----------------------------------------------------------
# Helper: safely combine saved MCMC chunks chain-by-chain.
# Each chunk is an mcmc.list with the same number of chains.
# ----------------------------------------------------------
combine_mcmc_chunks <- function(chunk_list) {

  if (length(chunk_list) == 0) {
    stop("No completed MCMC chunks are available.")
  }

  n_chain_list <- vapply(
    chunk_list,
    length,
    integer(1)
  )

  if (any(n_chain_list != n_chain_list[1])) {
    stop("Saved MCMC chunks have inconsistent chain counts.")
  }

  nch <- n_chain_list[1]

  combined_chains <- lapply(
    seq_len(nch),
    function(ch) {

      mats <- lapply(
        chunk_list,
        function(x) {
          as.matrix(x[[ch]])
        }
      )

      coda::mcmc(
        do.call(
          rbind,
          mats
        ),
        start = 1,
        thin = 1
      )
    }
  )

  coda::mcmc.list(
    combined_chains
  )
}


# ----------------------------------------------------------
# Helper: list and read completed chunks for a phase.
# ----------------------------------------------------------
get_chunk_files <- function(prefix) {

  files <- list.files(
    CHECKPOINT_DIR,
    pattern = paste0(
      "^",
      prefix,
      "_chunk_[0-9]{3}\\.rds$"
    ),
    full.names = TRUE
  )

  sort(files)
}


load_all_completed_chunks <- function() {

  files <- c(
    get_chunk_files("initial"),
    get_chunk_files("extend")
  )

  if (length(files) == 0) {
    return(list())
  }

  lapply(
    files,
    readRDS
  )
}


# ----------------------------------------------------------
# Helper: write progress immediately after every chunk.
# ----------------------------------------------------------
write_checkpoint_status <- function(
  phase,
  chunk_id,
  raw_iter_completed_this_run,
  total_completed_chunk_files,
  combined_samples
) {

  retained_per_chain <- nrow(
    as.matrix(
      combined_samples[[1]]
    )
  )

  status <- tibble(
    timestamp = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    run_mode = RUN_MODE,
    phase = phase,
    last_completed_chunk = chunk_id,
    raw_iter_completed_this_run =
      raw_iter_completed_this_run,
    total_completed_chunk_files =
      total_completed_chunk_files,
    chains = length(combined_samples),
    retained_samples_per_chain =
      retained_per_chain,
    total_retained_draws =
      retained_per_chain *
        length(combined_samples)
  )

  write_csv(
    status,
    file.path(
      OUT_DIR,
      "MCMC_PROGRESS.csv"
    )
  )

  saveRDS(
    combined_samples,
    file.path(
      OUT_DIR,
      "18_joint_mcmc_samples_PARTIAL.rds"
    )
  )

  invisible(status)
}


# ----------------------------------------------------------
# Helper: convergence diagnostics.
# ----------------------------------------------------------
get_convergence <- function(
  samples_obj
) {

  if (length(samples_obj) < 2) {
    stop("Rhat requires at least two chains.")
  }

  rhat <- gelman.diag(
    samples_obj,
    multivariate = FALSE,
    autoburnin = FALSE
  )$psrf[
    ,
    "Point est."
  ]

  ess <- effectiveSize(
    samples_obj
  )

  tibble(
    parameter =
      names(rhat),

    Rhat =
      as.numeric(rhat),

    ESS =
      as.numeric(
        ess[
          names(rhat)
        ]
      ),

    Rhat_ok =
      Rhat < 1.05,

    ESS_ok =
      ESS >= 400
  )
}


# ----------------------------------------------------------
# Determine what has already been saved.
# ----------------------------------------------------------
initial_chunk_files <- get_chunk_files(
  "initial"
)

extend_chunk_files <- get_chunk_files(
  "extend"
)

N_INITIAL_CHUNKS <- N_ITER %/% CHUNK_ITER
N_EXTEND_CHUNKS <- if (AUTO_EXTEND) {
  EXTEND_ITER %/% EXTEND_CHUNK_ITER
} else {
  0L
}

cat(
  "\n================ CHECKPOINT STATUS ================\n",
  "Completed initial chunks found: ",
  length(initial_chunk_files),
  " / ",
  N_INITIAL_CHUNKS,
  "\n",
  "Completed extension chunks found: ",
  length(extend_chunk_files),
  " / ",
  N_EXTEND_CHUNKS,
  "\n",
  "===================================================\n"
)


# ----------------------------------------------------------
# Compile a fresh JAGS model only if more sampling is needed.
#
# Important restart behavior:
# - Completed chunks are never overwritten.
# - On a new R session, JAGS itself cannot resume an external
#   pointer exactly where the process was killed.
# - Therefore the new session recompiles, burns in again, and
#   adds NEW posterior chunks to the already saved chunks.
# - This preserves all completed posterior draws.
# ----------------------------------------------------------
need_initial_sampling <-
  length(initial_chunk_files) <
    N_INITIAL_CHUNKS

need_possible_extension <-
  AUTO_EXTEND &&
  length(extend_chunk_files) <
    N_EXTEND_CHUNKS

if (
  need_initial_sampling ||
    need_possible_extension
) {

  cat(
    "\n================ COMPILING JOINT MODEL ================\n"
  )

  jm <- jags.model(
    file =
      MODEL_FILE,

    data =
      jags_data,

    inits =
      inits,

    n.chains =
      N_CHAINS,

    n.adapt =
      N_ADAPT
  )

  cat(
    "\n================ BURN-IN ================\n"
  )

  update(
    jm,
    n.iter =
      N_BURN
  )
}


# ----------------------------------------------------------
# INITIAL POSTERIOR: run only missing chunks.
# ----------------------------------------------------------
if (need_initial_sampling) {

  start_chunk <-
    length(initial_chunk_files) + 1L

  cat(
    "\n================ CHUNKED POSTERIOR SAMPLE ================\n",
    "Starting at initial chunk ",
    start_chunk,
    " of ",
    N_INITIAL_CHUNKS,
    "\n"
  )

  raw_iter_completed_this_run <- 0L

  for (
    chunk_id in seq.int(
      start_chunk,
      N_INITIAL_CHUNKS
    )
  ) {

    cat(
      "\n--- Initial chunk ",
      chunk_id,
      "/",
      N_INITIAL_CHUNKS,
      " : ",
      CHUNK_ITER,
      " raw iterations, thin=",
      N_THIN,
      " ---\n"
    )

    chunk_samples <- coda.samples(
      model =
        jm,

      variable.names =
        parameters_to_monitor,

      n.iter =
        CHUNK_ITER,

      thin =
        N_THIN
    )

    chunk_file <- file.path(
      CHECKPOINT_DIR,
      sprintf(
        "initial_chunk_%03d.rds",
        chunk_id
      )
    )

    # Save the completed chunk FIRST.
    saveRDS(
      chunk_samples,
      chunk_file
    )

    raw_iter_completed_this_run <-
      raw_iter_completed_this_run +
      CHUNK_ITER

    # Re-load all completed chunks so the combined partial file
    # is always reconstructable even after an interruption.
    all_chunks_now <-
      load_all_completed_chunks()

    samples_partial <-
      combine_mcmc_chunks(
        all_chunks_now
      )

    write_checkpoint_status(
      phase = "initial",
      chunk_id = chunk_id,
      raw_iter_completed_this_run =
        raw_iter_completed_this_run,
      total_completed_chunk_files =
        length(all_chunks_now),
      combined_samples =
        samples_partial
    )

    cat(
      "Saved: ",
      chunk_file,
      "\n",
      "Current retained samples per chain: ",
      nrow(
        as.matrix(
          samples_partial[[1]]
        )
      ),
      "\n"
    )
  }
}


# ----------------------------------------------------------
# Load all INITIAL chunks and check convergence.
# ----------------------------------------------------------
initial_chunk_files <- get_chunk_files(
  "initial"
)

if (
  length(initial_chunk_files) <
    N_INITIAL_CHUNKS
) {
  stop(
    "Initial sampling is incomplete. ",
    "Completed chunks remain safely saved in ",
    CHECKPOINT_DIR,
    ". Re-run the same script to continue."
  )
}

initial_samples <- combine_mcmc_chunks(
  lapply(
    initial_chunk_files,
    readRDS
  )
)

convergence <- get_convergence(
  initial_samples
)

cat(
  "\n================ CONVERGENCE AFTER INITIAL TARGET ================\n"
)

cat(
  "Maximum Rhat:",
  max(
    convergence$Rhat,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum ESS:",
  min(
    convergence$ESS,
    na.rm = TRUE
  ),
  "\n"
)


# ----------------------------------------------------------
# OPTIONAL EXTENSION.
# It is also chunked and every completed extension chunk is
# immediately saved.
# ----------------------------------------------------------
needs_extension_now <-
  AUTO_EXTEND &&
  any(
    !convergence$Rhat_ok |
      !convergence$ESS_ok
  )

if (needs_extension_now) {

  cat(
    "\nConvergence target not reached after initial sampling.\n"
  )

  extend_chunk_files <- get_chunk_files(
    "extend"
  )

  if (
    length(extend_chunk_files) <
      N_EXTEND_CHUNKS
  ) {

    # If this script was started with all initial chunks already
    # present, jm may not exist yet. Compile + burn before adding
    # fresh extension chunks.
    if (!exists("jm")) {

      cat(
        "\n================ RECOMPILING FOR EXTENSION ================\n"
      )

      jm <- jags.model(
        file =
          MODEL_FILE,

        data =
          jags_data,

        inits =
          inits,

        n.chains =
          N_CHAINS,

        n.adapt =
          N_ADAPT
      )

      update(
        jm,
        n.iter =
          N_BURN
      )
    }

    if (EXTEND_BURN > 0) {
      update(
        jm,
        n.iter =
          EXTEND_BURN
      )
    }

    start_extend <-
      length(extend_chunk_files) + 1L

    raw_iter_completed_this_run <- 0L

    for (
      chunk_id in seq.int(
        start_extend,
        N_EXTEND_CHUNKS
      )
    ) {

      cat(
        "\n--- Extension chunk ",
        chunk_id,
        "/",
        N_EXTEND_CHUNKS,
        " : ",
        EXTEND_CHUNK_ITER,
        " raw iterations, thin=",
        EXTEND_THIN,
        " ---\n"
      )

      chunk_samples <- coda.samples(
        model =
          jm,

        variable.names =
          parameters_to_monitor,

        n.iter =
          EXTEND_CHUNK_ITER,

        thin =
          EXTEND_THIN
      )

      chunk_file <- file.path(
        CHECKPOINT_DIR,
        sprintf(
          "extend_chunk_%03d.rds",
          chunk_id
        )
      )

      saveRDS(
        chunk_samples,
        chunk_file
      )

      raw_iter_completed_this_run <-
        raw_iter_completed_this_run +
        EXTEND_CHUNK_ITER

      all_chunks_now <-
        load_all_completed_chunks()

      samples_partial <-
        combine_mcmc_chunks(
          all_chunks_now
        )

      write_checkpoint_status(
        phase = "extend",
        chunk_id = chunk_id,
        raw_iter_completed_this_run =
          raw_iter_completed_this_run,
        total_completed_chunk_files =
          length(all_chunks_now),
        combined_samples =
          samples_partial
      )

      cat(
        "Saved: ",
        chunk_file,
        "\n",
        "Current retained samples per chain: ",
        nrow(
          as.matrix(
            samples_partial[[1]]
          )
        ),
        "\n"
      )
    }
  }
}


# ----------------------------------------------------------
# FINAL AVAILABLE SAMPLE = every completed saved chunk.
# Even if no extension was needed, this gives the initial set.
# ----------------------------------------------------------
all_completed_chunks <-
  load_all_completed_chunks()

samples <- combine_mcmc_chunks(
  all_completed_chunks
)

saveRDS(
  samples,
  file.path(
    OUT_DIR,
    "18_joint_mcmc_samples.rds"
  )
)

convergence <- get_convergence(
  samples
)

cat(
  "\n================ FINAL AVAILABLE CONVERGENCE ================\n"
)

cat(
  "Completed checkpoint chunks:",
  length(all_completed_chunks),
  "\n"
)

cat(
  "Retained samples per chain:",
  nrow(
    as.matrix(
      samples[[1]]
    )
  ),
  "\n"
)

cat(
  "Maximum Rhat:",
  max(
    convergence$Rhat,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum ESS:",
  min(
    convergence$ESS,
    na.rm = TRUE
  ),
  "\n"
)


write_csv(
  convergence,
  file.path(
    OUT_DIR,
    "07_convergence_Rhat_ESS.csv"
  )
)


############################################################
# 21. TRACE / DENSITY / GELMAN PLOTS
############################################################

pdf(
  file.path(
    OUT_DIR,
    "08_trace_density.pdf"
  ),
  width = 12,
  height = 8
)

plot(samples)

dev.off()


pdf(
  file.path(
    OUT_DIR,
    "09_gelman_plots.pdf"
  ),
  width = 12,
  height = 8
)

gelman.plot(
  samples,
  autoburnin = FALSE
)

dev.off()


############################################################
# 22. POSTERIOR MATRIX
############################################################

post <- as.matrix(
  samples
)


posterior_summary_one <- function(x) {

  qs <- quantile(
    x,
    probs =
      c(
        0.025,
        0.5,
        0.975
      ),
    names = FALSE
  )

  p_gt0 <- mean(
    x > 0
  )

  c(
    Mean = mean(x),
    SD = sd(x),
    Median = qs[2],
    CrI_2.5 = qs[1],
    CrI_97.5 = qs[3],
    Pr_gt_0 = p_gt0,
    Direction_probability =
      max(
        p_gt0,
        1 - p_gt0
      )
  )
}


############################################################
# 23. LONGITUDINAL FIXED-EFFECT RESULTS
############################################################

longitudinal_covariate_results <- map_dfr(

  1:8,

  function(p) {

    parname <- sprintf(
      "beta_long[%d]",
      p
    )

    x <- post[
      ,
      parname
    ]

    s <- posterior_summary_one(
      x
    )

    tibble(
      p = p,
      parameter = parname,
      beta_mean =
        unname(s["Mean"]),
      beta_sd =
        unname(s["SD"]),
      beta_median =
        unname(s["Median"]),
      beta_CrI_2.5 =
        unname(s["CrI_2.5"]),
      beta_CrI_97.5 =
        unname(s["CrI_97.5"]),
      Significant_95CrI =
        (
          unname(s["CrI_2.5"]) > 0 |
            unname(s["CrI_97.5"]) < 0
        )
    )
  }
) %>%
  dplyr::left_join(
    covariate_map,
    by = "p"
  ) %>%
  dplyr::select(
    covariate,
    beta_mean,
    beta_sd,
    beta_median,
    beta_CrI_2.5,
    beta_CrI_97.5,
    Significant_95CrI,
    parameter
  )


longitudinal_core_results <- map_dfr(

  c(
    "beta0_y",
    "beta_time_y",
    "sigma_y",
    "sigma_b0",
    "sigma_b1",
    "rho"
  ),

  function(parname) {

    x <- post[
      ,
      parname
    ]

    qs <- quantile(
      x,
      c(
        0.025,
        0.5,
        0.975
      ),
      names = FALSE
    )

    tibble(
      parameter =
        parname,

      mean =
        mean(x),

      sd =
        sd(x),

      median =
        qs[2],

      CrI_2.5 =
        qs[1],

      CrI_97.5 =
        qs[3]
    )
  }
)


cat(
  "\n================ LONGITUDINAL CORE RESULTS ================\n"
)

print(
  longitudinal_core_results,
  n = Inf
)


cat(
  "\n================ LONGITUDINAL COVARIATE RESULTS ================\n"
)

print(
  longitudinal_covariate_results,
  n = Inf,
  width = Inf
)


write_csv(
  longitudinal_core_results,
  file.path(
    OUT_DIR,
    "10_longitudinal_core_results.csv"
  )
)


write_csv(
  longitudinal_covariate_results,
  file.path(
    OUT_DIR,
    "11_longitudinal_8cov_results.csv"
  )
)


############################################################
# 24. MULTISTATE GAMMA / HR RESULTS
############################################################

regression_results <- map_dfr(

  1:9,

  function(k) {

    map_dfr(

      1:8,

      function(p) {

        parname <- sprintf(
          "gamma[%d,%d]",
          k,
          p
        )

        x <- post[
          ,
          parname
        ]

        s <- posterior_summary_one(
          x
        )

        low <- unname(
          s["CrI_2.5"]
        )

        high <- unname(
          s["CrI_97.5"]
        )

        tibble(
          k = k,
          p = p,
          parameter = parname,

          beta_mean =
            unname(
              s["Mean"]
            ),

          beta_sd =
            unname(
              s["SD"]
            ),

          beta_median =
            unname(
              s["Median"]
            ),

          beta_CrI_2.5 =
            low,

          beta_CrI_97.5 =
            high,

          HR =
            exp(
              unname(
                s["Mean"]
              )
            ),

          HR_CrI_2.5 =
            exp(low),

          HR_CrI_97.5 =
            exp(high),

          Significant_95CrI =
            (
              low > 0 |
                high < 0
            )
        )
      }
    )
  }
) %>%
  dplyr::left_join(
    transition_map,
    by = "k"
  ) %>%
  dplyr::left_join(
    covariate_map,
    by = "p"
  ) %>%
  dplyr::select(
    transition,
    from,
    to,
    covariate,
    beta_mean,
    beta_sd,
    beta_median,
    beta_CrI_2.5,
    beta_CrI_97.5,
    HR,
    HR_CrI_2.5,
    HR_CrI_97.5,
    Significant_95CrI,
    parameter
  )


cat(
  "\n================ MULTISTATE 8-COVARIATE RESULTS ================\n"
)

print(
  regression_results,
  n = Inf,
  width = Inf
)


write_csv(
  regression_results,
  file.path(
    OUT_DIR,
    "12_multistate_8cov_HR_results.csv"
  )
)


############################################################
# 25. JOINT ASSOCIATION RESULTS
############################################################

association_results <- map_dfr(

  1:9,

  function(k) {

    map_dfr(

      c(
        "alpha0",
        "alpha1"
      ),

      function(type) {

        parname <- sprintf(
          "%s[%d]",
          type,
          k
        )

        x <- post[
          ,
          parname
        ]

        s <- posterior_summary_one(
          x
        )

        low <- unname(
          s["CrI_2.5"]
        )

        high <- unname(
          s["CrI_97.5"]
        )

        tibble(
          k = k,

          association =
            ifelse(
              type == "alpha0",
              "Shared random intercept",
              "Shared random slope"
            ),

          parameter =
            parname,

          mean =
            unname(
              s["Mean"]
            ),

          sd =
            unname(
              s["SD"]
            ),

          median =
            unname(
              s["Median"]
            ),

          CrI_2.5 =
            low,

          CrI_97.5 =
            high,

          HR_per_1unit_RE =
            exp(
              unname(
                s["Mean"]
              )
            ),

          HR_CrI_2.5 =
            exp(low),

          HR_CrI_97.5 =
            exp(high),

          Significant_95CrI =
            (
              low > 0 |
                high < 0
            )
        )
      }
    )
  }
) %>%
  dplyr::left_join(
    transition_map,
    by = "k"
  ) %>%
  dplyr::select(
    transition,
    from,
    to,
    association,
    mean,
    sd,
    median,
    CrI_2.5,
    CrI_97.5,
    HR_per_1unit_RE,
    HR_CrI_2.5,
    HR_CrI_97.5,
    Significant_95CrI,
    parameter
  )


cat(
  "\n================ JOINT ASSOCIATION RESULTS ================\n"
)

print(
  association_results,
  n = Inf,
  width = Inf
)


write_csv(
  association_results,
  file.path(
    OUT_DIR,
    "13_joint_association_alpha_results.csv"
  )
)


############################################################
# 26. BASELINE TRANSITION INTENSITIES
############################################################

baseline_results <- map_dfr(

  1:9,

  function(k) {

    parname <- sprintf(
      "lambda0[%d]",
      k
    )

    x <- post[
      ,
      parname
    ]

    qs <- quantile(
      x,
      c(
        0.025,
        0.5,
        0.975
      ),
      names = FALSE
    )

    tibble(
      k = k,

      lambda0_mean =
        mean(x),

      lambda0_sd =
        sd(x),

      lambda0_median =
        qs[2],

      lambda0_CrI_2.5 =
        qs[1],

      lambda0_CrI_97.5 =
        qs[3]
    )
  }
) %>%
  dplyr::left_join(
    transition_map,
    by = "k"
  ) %>%
  dplyr::select(
    transition,
    from,
    to,
    lambda0_mean,
    lambda0_sd,
    lambda0_median,
    lambda0_CrI_2.5,
    lambda0_CrI_97.5
  )


write_csv(
  baseline_results,
  file.path(
    OUT_DIR,
    "14_baseline_transition_intensities.csv"
  )
)


############################################################
# 27. DIC
############################################################

if (RUN_DIC) {

  cat(
    "\n================ DIC ================\n"
  )

  if (!exists("jm")) {

    jm <- jags.model(
      file = MODEL_FILE,
      data = jags_data,
      inits = inits,
      n.chains = N_CHAINS,
      n.adapt = N_ADAPT
    )

    update(
      jm,
      n.iter = N_BURN
    )
  }

  dic_obj <- dic.samples(
    jm,
    n.iter = 5000,
    thin = 5,
    type = "pD"
  )

  print(dic_obj)

  capture.output(
    print(dic_obj),
    file =
      file.path(
        OUT_DIR,
        "15_DIC.txt"
      )
  )
}


############################################################
# 28. FOREST PLOT FOR MULTISTATE COVARIATE EFFECTS
############################################################

plot_data <- regression_results %>%
  dplyr::mutate(

    transition =
      factor(
        transition,
        levels =
          transition_map$transition
      ),

    covariate =
      factor(
        covariate,
        levels =
          rev(
            covariate_map$covariate
          )
      )
  )


p_forest <- ggplot(
  plot_data,
  aes(
    x = HR,
    y = covariate
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2
  ) +
  geom_errorbarh(
    aes(
      xmin =
        HR_CrI_2.5,
      xmax =
        HR_CrI_97.5
    ),
    height = 0.15
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    ~ transition,
    scales = "free_x",
    ncol = 3
  ) +
  scale_x_log10() +
  labs(
    x =
      "Hazard ratio (posterior mean; 95% credible interval)",
    y = NULL,
    title =
      "Joint model: transition-specific effects of 8 baseline covariates"
  ) +
  theme_bw(
    base_size = 11
  )


ggsave(
  filename =
    file.path(
      OUT_DIR,
      "16_multistate_8cov_HR_forest_plot.pdf"
    ),
  plot = p_forest,
  width = 13,
  height = 10
)


ggsave(
  filename =
    file.path(
      OUT_DIR,
      "16_multistate_8cov_HR_forest_plot.png"
    ),
  plot = p_forest,
  width = 13,
  height = 10,
  dpi = 300
)


############################################################
# 29. SAVE DATA / MCMC / COMPLETE RESULTS
############################################################

write_csv(
  ms_subject,
  file.path(
    OUT_DIR,
    "17_subject_level_joint_model_data.csv"
  )
)


saveRDS(
  samples,
  file.path(
    OUT_DIR,
    "18_joint_mcmc_samples.rds"
  )
)


saveRDS(
  list(

    transition_map =
      transition_map,

    covariate_map =
      covariate_map,

    age_center =
      AGE_CENTER,

    heeye_center =
      EYE_CENTER,

    baseline_covariates =
      baseline_cov,

    multistate_subject_data =
      ms_subject,

    longitudinal_data =
      longitudinal_data,

    convergence =
      convergence,

    longitudinal_core_results =
      longitudinal_core_results,

    longitudinal_covariate_results =
      longitudinal_covariate_results,

    multistate_results =
      regression_results,

    association_results =
      association_results,

    baseline_intensity_results =
      baseline_results
  ),

  file.path(
    OUT_DIR,
    "19_complete_joint_analysis_results.rds"
  )
)


############################################################
# 30. FINAL MODEL CHECK
############################################################

cat(
  "\n============================================================\n"
)

cat(
  "FINAL JOINT MODEL CHECK SUMMARY\n"
)

cat(
  "============================================================\n"
)

cat(
  "Multistate intervals retained:",
  nrow(ms),
  "of",
  n_rows_before,
  "\n"
)

cat(
  "Participants:",
  nrow(ms_subject),
  "\n"
)

cat(
  "All 9 transition totals retained:",
  all(
    aggregate_check$retained
  ),
  "\n"
)

cat(
  "Number of multistate baseline covariates:",
  ncol(X_ms),
  "\n"
)

cat(
  "Number of longitudinal baseline covariates:",
  ncol(X_long),
  "\n"
)

cat(
  "Baseline delayed recall included in multistate X:",
  FALSE,
  "\n"
)

cat(
  "Longitudinal observations:",
  nrow(longitudinal_data),
  "\n"
)

cat(
  "Subjects with longitudinal observations:",
  sum(
    longitudinal_subject_audit$n_long_obs > 0
  ),
  "\n"
)

cat(
  "Maximum Rhat:",
  round(
    max(
      convergence$Rhat,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Minimum ESS:",
  round(
    min(
      convergence$ESS,
      na.rm = TRUE
    ),
    1
  ),
  "\n"
)

cat(
  "\nFINAL MODEL STRUCTURE:\n"
)

cat(
  "Longitudinal:\n"
)

cat(
  "  cflisd ~ time + age + cardiovascular + BP + diabetes + ",
  "stroke + cancer + arthritis + eyesight + random intercept + random slope\n"
)

cat(
  "Multistate:\n"
)

cat(
  "  9 transition intensities ~ age + cardiovascular + BP + diabetes + ",
  "stroke + cancer + arthritis + eyesight\n"
)

cat(
  "Joint association:\n"
)

cat(
  "  each transition has alpha0[k]*random_intercept + ",
  "alpha1[k]*random_slope\n"
)

cat(
  "\nAll outputs saved in:",
  OUT_DIR,
  "\n"
)

cat(
  "============================================================\n"
)


cat(
  "Checkpoint directory:",
  CHECKPOINT_DIR,
  "\n"
)

cat(
  "To run the full dataset, change only: RUN_MODE <- \"FULL\"\n"
)

cat(
  "If interrupted, re-run the same script. Completed chunk files are reused.\n"
)
