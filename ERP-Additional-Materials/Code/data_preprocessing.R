############################################################
# 0. PACKAGES
############################################################
required_packages <- c(
  "haven", "dplyr", "tidyr", "purrr", "stringr", "readr",
  "survival", "ggplot2", "survminer"
)

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages) > 0) install.packages(new_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

############################################################
# 1. SETTINGS
############################################################
# Use the actual input file name. Change this only if your raw
# Stata file has a different name.
DATA_FILE <- "new_dataset_noL 1.dta"
OUT_DIR <- "ELSA_JOINT_MODEL_DATA_FINAL"

ID_VAR <- "idauniq"
WAVE_VAR <- "wave"
YEAR_VAR <- "year"
AGE_VAR <- "age"
LONGITUDINAL_OUTCOME <- "cflisd"
MULTISTATE_STATE <- "frailgr"
DEATH_YEAR_VAR <- "raxyear"
LAST_ALIVE_YEAR_VAR <- "ralstcorey"
EOL_INTERVIEW_YEAR_VAR <- "raxtiwy"

# Primary covariates for the joint analysis.
# These variables are retained throughout the processing pipeline.
PRIMARY_COVARIATES <- c("age", "hedibar", "heeye", "hedimbp")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

############################################################
# 2. READ RAW DATA
############################################################
cat("\n====================================================\n")
cat("READING RAW ELSA DATA\n")
cat("====================================================\n")
raw <- read_dta(DATA_FILE)
cat("Rows:", nrow(raw), "\n")
cat("Columns:", ncol(raw), "\n")

raw_names <- names(raw)
empty_name_index <- which(is.na(raw_names) | raw_names == "")
if (length(empty_name_index) > 0) {
  names(raw)[empty_name_index] <- paste0("unnamed_variable_", empty_name_index)
}

required_variables <- c(
  ID_VAR, WAVE_VAR, YEAR_VAR, AGE_VAR,
  LONGITUDINAL_OUTCOME, MULTISTATE_STATE
)
missing_required <- setdiff(required_variables, names(raw))
if (length(missing_required) > 0) {
  stop(paste("Missing required variables:", paste(missing_required, collapse = ", ")))
}

death_variable_available <- DEATH_YEAR_VAR %in% names(raw)
last_alive_available <- LAST_ALIVE_YEAR_VAR %in% names(raw)
eol_year_available <- EOL_INTERVIEW_YEAR_VAR %in% names(raw)

cat("Death variable:", DEATH_YEAR_VAR, "available =", death_variable_available, "\n")
cat("Last alive core year:", LAST_ALIVE_YEAR_VAR, "available =", last_alive_available, "\n")
cat("EOL interview year:", EOL_INTERVIEW_YEAR_VAR, "available =", eol_year_available, "\n")

############################################################
# 3. CORE VARIABLE CONVERSION AND CLEANING
############################################################
elsa <- raw %>%
  mutate(
    idauniq = as.numeric(.data[[ID_VAR]]),
    wave = as.numeric(.data[[WAVE_VAR]]),
    year = as.numeric(.data[[YEAR_VAR]]),
    age_original = as.numeric(.data[[AGE_VAR]]),
    cflisd = as.numeric(.data[[LONGITUDINAL_OUTCOME]]),
    frailgr = as.numeric(.data[[MULTISTATE_STATE]])
  ) %>%
  mutate(
    year_valid = !is.na(year) & is.finite(year) & year >= 1900 & year <= 2100,
    cflisd = if_else(!is.na(cflisd) & cflisd >= 0, cflisd, NA_real_),
    frailgr = if_else(frailgr %in% c(0, 1, 2), frailgr, NA_real_),
    age_original = if_else(!is.na(age_original) & age_original >= 50 & age_original <= 110,
                           age_original, NA_real_)
  ) %>%
  arrange(idauniq, year, wave)

############################################################
# 4. FRAILTY DISTRIBUTION
############################################################
frailty_distribution <- elsa %>%
  count(frailgr, name = "n") %>%
  mutate(percent = n / sum(n) * 100)
print(frailty_distribution)
write.csv(frailty_distribution,
          file.path(OUT_DIR, "00_frailty_distribution.csv"), row.names = FALSE)

############################################################
# 5. AGE RECONSTRUCTION
############################################################
age_information <- elsa %>%
  filter(!is.na(idauniq), year_valid, !is.na(age_original)) %>%
  mutate(
    birth_year_estimate = year - age_original,
    # 新增：过滤出生年份不合理的记录（早于1900年或晚于访谈年份）
    birth_year_valid = birth_year_estimate >= 1900 & birth_year_estimate <= year
  ) %>%
  filter(birth_year_valid) # 仅保留出生年份合理的记录

birth_year_data <- age_information %>%
  group_by(idauniq) %>%
  summarise(
    n_age_observed = n(),
    birth_year_median = median(birth_year_estimate, na.rm = TRUE),
    birth_year_min = min(birth_year_estimate, na.rm = TRUE),
    birth_year_max = max(birth_year_estimate, na.rm = TRUE),
    birth_year_range = birth_year_max - birth_year_min,
    .groups = "drop"
  )

elsa <- elsa %>% left_join(birth_year_data, by = "idauniq") %>%
  mutate(
    age_reconstructed = case_when(
      is.na(age_original) & !is.na(birth_year_median) & year_valid &
        birth_year_range <= 2 ~ year - birth_year_median,
      TRUE ~ NA_real_
    ),
    age_final = coalesce(age_original, age_reconstructed),
    age_source = case_when(
      !is.na(age_original) ~ "observed",
      !is.na(age_reconstructed) ~ "reconstructed",
      TRUE ~ "missing"
    )
  )

age_outlier_check <- elsa %>%
  filter(!is.na(age_final), age_final < 40 | age_final > 110) %>%
  select(idauniq, wave, year, age_original, age_reconstructed, age_final, age_source)
write.csv(age_outlier_check,
          file.path(OUT_DIR, "01_age_outlier_check.csv"), row.names = FALSE)

age_audit <- elsa %>% count(age_source, name = "n") %>%
  mutate(percent = n / sum(n) * 100)
write.csv(age_audit,
          file.path(OUT_DIR, "02_age_reconstruction_audit.csv"), row.names = FALSE)

############################################################
# 6. DUPLICATE ID‑WAVE AUDIT
############################################################
duplicate_id_wave <- elsa %>%
  filter(!is.na(idauniq), !is.na(wave)) %>%
  count(idauniq, wave, name = "n") %>%
  filter(n > 1)
write.csv(duplicate_id_wave,
          file.path(OUT_DIR, "03_duplicate_id_wave.csv"), row.names = FALSE)

############################################################
# 7. PERSON‑LEVEL DEATH DATA
############################################################
if (death_variable_available) {
  death_data <- elsa %>%
    transmute(
      idauniq,
      death_year_raw = as.numeric(.data[[DEATH_YEAR_VAR]]),
      last_alive_year_raw = if (last_alive_available)
        as.numeric(.data[[LAST_ALIVE_YEAR_VAR]]) else NA_real_,
      eol_interview_year_raw = if (eol_year_available)
        as.numeric(.data[[EOL_INTERVIEW_YEAR_VAR]]) else NA_real_
    ) %>%
    filter(!is.na(idauniq)) %>%
    mutate(
      death_year = if_else(death_year_raw >= 1900 & death_year_raw <= 2100,
                           death_year_raw, NA_real_),
      last_alive_year = if_else(last_alive_year_raw >= 1900 &
                                  last_alive_year_raw <= 2100,
                                last_alive_year_raw, NA_real_),
      eol_interview_year = if_else(eol_interview_year_raw >= 1900 &
                                     eol_interview_year_raw <= 2100,
                                   eol_interview_year_raw, NA_real_)
    ) %>%
    group_by(idauniq) %>%
    summarise(
      death_year = if (all(is.na(death_year))) NA_real_ else min(death_year, na.rm = TRUE),
      last_alive_year = if (all(is.na(last_alive_year))) NA_real_ else max(last_alive_year, na.rm = TRUE),
      eol_interview_year = if (all(is.na(eol_interview_year))) NA_real_ else max(eol_interview_year, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      death_observed = !is.na(death_year),
      death_time_precision = if_else(death_observed, "year‑level", "none")
    )
} else {
  death_data <- tibble(
    idauniq = numeric(0), death_year = numeric(0), last_alive_year = numeric(0),
    eol_interview_year = numeric(0), death_observed = logical(0),
    death_time_precision = character(0)
  )
}

write.csv(death_data,
          file.path(OUT_DIR, "04_death_data_person_level.csv"), row.names = FALSE)

death_audit <- death_data %>% summarise(
  participants_with_death = sum(death_observed),
  participants_without_death = sum(!death_observed),
  total_participants = n()
)
write.csv(death_audit,
          file.path(OUT_DIR, "05_death_audit.csv"), row.names = FALSE)

############################################################
# 8. CLEAN PERSON‑WAVE MASTER DATA
############################################################
# IMPORTANT:
# Do NOT select only the core analysis variables here.
# The raw ELSA data contain additional covariates such as age,
# hedibar, heeye and hedimbp.  If select(...) is restricted to
# idauniq/wave/year/age/cflisd/frailgr, all other variables are
# silently dropped and cannot be recovered in downstream files.
#
# Therefore the master person-wave dataset retains ALL variables
# from elsa, together with the person-level death information.
person_wave <- elsa %>%
  left_join(death_data, by = "idauniq") %>%
  arrange(idauniq, year, wave)

# Check that the intended primary covariates are still present.
missing_primary_covariates <- setdiff(PRIMARY_COVARIATES, names(person_wave))
if (length(missing_primary_covariates) > 0) {
  stop(
    paste(
      "Primary covariates missing from person_wave:",
      paste(missing_primary_covariates, collapse = ", ")
    )
  )
}

write.csv(
  person_wave,
  file.path(OUT_DIR, "06_ELSA_CLEAN_PERSON_WAVE.csv"),
  row.names = FALSE
)

# A compact analysis-ready person-wave file containing the variables
# needed directly for the joint analysis.  The full master file above
# remains available for additional analyses and auditing.
joint_core_vars <- c(
  "idauniq", "wave", "year", "time_from_baseline",
  "age", "age_original", "age_reconstructed", "age_final", "age_source",
  "hedibar", "heeye", "hedimbp",
  "cflisd", "frailgr",
  "death_year", "last_alive_year", "eol_interview_year"
)

# time_from_baseline is created after the baseline step, so this
# compact file is generated later in Section 19 together with the
# joint datasets.

############################################################
# 9. GLOBAL BASELINE
# PATCH‑1: add filter: is.na(death_year)|year < death_year
# Baseline = earliest valid observation with valid age and frailty.
# cflisd is not required for baseline.
############################################################
global_baseline <- person_wave %>%
  filter(
    !is.na(idauniq), year_valid, !is.na(age_final),
    frailgr %in% c(0, 1, 2),
    age_final >= 50,
    ## PATCH‑1 修复：禁止死亡当年观测选为基线
    is.na(death_year) | year < death_year
  ) %>%
  group_by(idauniq) %>%
  arrange(year, wave, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    idauniq,
    baseline_year = year,
    baseline_wave = wave,
    baseline_age = age_final,
    baseline_frailty = frailgr,
    baseline_death_year = death_year,
    baseline_death_same_year = !is.na(death_year) & death_year == year,
    baseline_death_before = !is.na(death_year) & death_year < year
  )


person_wave <- person_wave %>%
  left_join(global_baseline, by = "idauniq") %>%
  mutate(time_from_baseline = year - baseline_year)

############################################################
# 10. BASELINE DEATH AUDIT
############################################################
baseline_death_audit <- global_baseline %>%
  mutate(
    status = case_when(
      baseline_death_before ~ "death_before_baseline",
      baseline_death_same_year ~ "death_same_year_as_baseline",
      !is.na(baseline_death_year) & baseline_death_year > baseline_year ~ "death_after_baseline",
      TRUE ~ "no_known_death"
    )
  ) %>%
  count(status, name = "n") %>%
  mutate(percent = n / sum(n) * 100)
write.csv(baseline_death_audit,
          file.path(OUT_DIR, "06A_baseline_death_audit.csv"), row.names = FALSE)

############################################################
# 11. LONGITUDINAL DATA
# PATCH‑2: add filter: is.na(death_year)|year < death_year
############################################################
longitudinal_data <- person_wave %>%
  filter(
    !is.na(idauniq), year_valid, !is.na(baseline_year), !is.na(cflisd),
    !is.na(time_from_baseline), time_from_baseline >= 0,
    ## PATCH‑2 修复：纵向观测与多状态统一，死亡当年不使用
    is.na(death_year) | year < death_year
  ) %>%
  arrange(idauniq, time_from_baseline, wave)

longitudinal_summary <- longitudinal_data %>%
  group_by(idauniq) %>%
  summarise(
    n_longitudinal_observations = n(),
    first_longitudinal_year = min(year),
    last_longitudinal_year = max(year),
    longitudinal_followup_years = max(year) - min(year),
    .groups = "drop"
  )
longitudinal_ids <- longitudinal_data %>% distinct(idauniq)

write.csv(longitudinal_data,
          file.path(OUT_DIR, "07_ELSA_LONGITUDINAL.csv"), row.names = FALSE)
write.csv(longitudinal_summary,
          file.path(OUT_DIR, "08_longitudinal_participant_summary.csv"), row.names = FALSE)

############################################################
# 12. MULTISTATE PERSON‑WAVE DATA
#
# Only observations strictly before a known death year are treated
# as alive frailty states. An observation in the death year is not
# automatically treated as pre‑death because death timing within
# the year is not released as an exact date.
############################################################
multistate_person_wave <- person_wave %>%
  filter(
    !is.na(idauniq), year_valid, !is.na(baseline_year),
    frailgr %in% c(0, 1, 2), time_from_baseline >= 0
  ) %>%
  filter(is.na(death_year) | year < death_year) %>%
  arrange(idauniq, time_from_baseline, wave)

multistate_summary <- multistate_person_wave %>%
  group_by(idauniq) %>%
  summarise(
    n_multistate_observations = n(),
    first_multistate_year = min(year),
    last_multistate_year = max(year),
    .groups = "drop"
  )
multistate_ids <- multistate_person_wave %>% distinct(idauniq)

write.csv(multistate_person_wave,
          file.path(OUT_DIR, "09_ELSA_MULTISTATE_PERSON_WAVE.csv"), row.names = FALSE)
write.csv(multistate_summary,
          file.path(OUT_DIR, "10_multistate_participant_summary.csv"), row.names = FALSE)

############################################################
# 13. NEXT OBSERVATION
############################################################
transition_candidates <- multistate_person_wave %>%
  group_by(idauniq) %>%
  arrange(year, wave, .by_group = TRUE) %>%
  mutate(
    next_year = lead(year),
    next_wave = lead(wave),
    next_state = lead(frailgr),
    next_time = lead(time_from_baseline)
  ) %>%
  ungroup()

############################################################
# 14. OBSERVED FRAILTY TRANSITIONS
#
# A transition is retained only if the next frailty observation
# occurs before any known death year.
############################################################
observed_transitions <- transition_candidates %>%
  filter(
    !is.na(next_year), !is.na(next_state),
    next_year > year,
    frailgr %in% c(0, 1, 2), next_state %in% c(0, 1, 2),
    is.na(death_year) | next_year < death_year
  ) %>%
  transmute(
    idauniq,
    from = frailgr,
    to = next_state,
    start = time_from_baseline,
    stop = next_time,
    calendar_start = year,
    calendar_stop = next_year,
    transition_source = "observed_frailty",
    transition_time_precision = "interval_between_waves"
  )

############################################################
# 15. SAME‑YEAR DEATH / NEXT‑WAVE AMBIGUITY AUDIT
############################################################
# These cases are not assigned a death transition because the order
# of the frailty observation and death is unknown within the year.
same_year_death_audit <- person_wave %>%
  filter(!is.na(death_year)) %>%
  filter(year == death_year) %>%
  select(idauniq, wave, year, frailgr, death_year) %>%
  distinct()

write.csv(same_year_death_audit,
          file.path(OUT_DIR, "11_same_year_death_audit.csv"), row.names = FALSE)

ambiguous_death_between_waves <- transition_candidates %>%
  filter(!is.na(death_year), !is.na(next_year), death_year == next_year) %>%
  select(idauniq, wave, year, frailgr, next_wave, next_year, next_state, death_year) %>%
  distinct()
write.csv(ambiguous_death_between_waves,
          file.path(OUT_DIR, "11A_death_same_year_as_next_wave_audit.csv"), row.names = FALSE)

############################################################
# 16. DEATH TRANSITIONS
#
# Death is assigned from the last observed alive frailty state only
# when death occurs strictly after that observation and before the
# next observed frailty wave. If death occurs in the same calendar
# year as the next wave, no exact ordering is invented.
############################################################
death_transition_candidates <- transition_candidates %>%
  filter(!is.na(death_year), death_year > year) %>%
  mutate(
    death_before_next = is.na(next_year) | death_year < next_year,
    death_same_year_as_next = !is.na(next_year) & death_year == next_year
  )

death_transitions <- death_transition_candidates %>%
  filter(death_before_next) %>%
  group_by(idauniq) %>%
  arrange(desc(year), desc(wave), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    idauniq,
    from = frailgr,
    to = 3,
    start = time_from_baseline,
    stop = death_year - baseline_year,
    calendar_start = year,
    calendar_stop = death_year,
    transition_source = "death",
    transition_time_precision = "year_level"
  ) %>%
  filter(!is.na(start), !is.na(stop), stop > start)

## PATCH‑4：死亡转移ID唯一性校验
dup_death_id <- death_transitions$idauniq[duplicated(death_transitions$idauniq)]
if(length(dup_death_id) > 0){
  warning("WARNING: Duplicated id found in death_transitions!")
  write.csv(death_transitions[death_transitions$idauniq %in% dup_death_id, ],
            file.path(OUT_DIR,"WARNING_dup_death_trans.csv"),row.names = FALSE)
}

############################################################
# 17. COMBINE TRANSITIONS
############################################################
multistate_transitions <- bind_rows(observed_transitions, death_transitions) %>%
  mutate(
    transition = paste0(from, "_to_", to),
    transition_type = case_when(
      to == 3 ~ "death",
      from == to ~ "stay",
      TRUE ~ "frailty_change"
    ),
    from_label = factor(from, levels = c(0, 1, 2),
                        labels = c("Non‑frail", "Pre‑frail", "Frail")),
    to_label = factor(to, levels = c(0, 1, 2, 3),
                      labels = c("Non‑frail", "Pre‑frail", "Frail", "Death")),
    interval = stop - start
  ) %>%
  filter(!is.na(start), !is.na(stop), stop > start) %>%
  arrange(idauniq, start, stop)

write.csv(multistate_transitions,
          file.path(OUT_DIR, "12_ELSA_MULTISTATE_TRANSITIONS.csv"), row.names = FALSE)

transition_summary <- multistate_transitions %>%
  count(from, to, transition, transition_type, name = "n") %>%
  arrange(from, to)
write.csv(transition_summary,
          file.path(OUT_DIR, "13_transition_summary.csv"), row.names = FALSE)

frailty_transition_matrix <- observed_transitions %>%
  count(from, to, name = "n") %>%
  complete(from = c(0, 1, 2), to = c(0, 1, 2), fill = list(n = 0)) %>%
  arrange(from, to)
write.csv(frailty_transition_matrix,
          file.path(OUT_DIR, "14_frailty_transition_matrix.csv"), row.names = FALSE)

death_transition_summary <- death_transitions %>%
  count(from, to, name = "n") %>% arrange(from)
write.csv(death_transition_summary,
          file.path(OUT_DIR, "15_death_transition_summary.csv"), row.names = FALSE)

############################################################
# 18. JOINT OVERLAP
############################################################
joint_ids <- inner_join(longitudinal_ids, multistate_ids, by = "idauniq")

joint_participant_audit <- tibble(
  analysis_group = c("Raw ELSA", "Longitudinal", "Multistate", "Joint overlap"),
  participants = c(
    n_distinct(person_wave$idauniq),
    n_distinct(longitudinal_ids$idauniq),
    n_distinct(multistate_ids$idauniq),
    n_distinct(joint_ids$idauniq)
  )
)
write.csv(joint_participant_audit,
          file.path(OUT_DIR, "16_joint_participant_audit.csv"), row.names = FALSE)

joint_longitudinal_data <- longitudinal_data %>% semi_join(joint_ids, by = "idauniq")
joint_multistate_person_wave <- multistate_person_wave %>% semi_join(joint_ids, by = "idauniq")
joint_multistate_transitions <- multistate_transitions %>% semi_join(joint_ids, by = "idauniq")

write.csv(joint_longitudinal_data,
          file.path(OUT_DIR, "17_ELSA_JOINT_LONGITUDINAL.csv"), row.names = FALSE)
write.csv(joint_multistate_person_wave,
          file.path(OUT_DIR, "18_ELSA_JOINT_MULTISTATE_PERSON_WAVE.csv"), row.names = FALSE)
write.csv(joint_multistate_transitions,
          file.path(OUT_DIR, "19_ELSA_JOINT_MULTISTATE_TRANSITIONS.csv"), row.names = FALSE)

############################################################
# 19. JOINT BASELINE DATASET
############################################################
joint_baseline <- global_baseline %>%
  semi_join(joint_ids, by = "idauniq") %>%
  left_join(
    longitudinal_data %>%
      group_by(idauniq) %>%
      arrange(time_from_baseline, wave, .by_group = TRUE) %>%
      filter(!is.na(cflisd)) %>%
      slice(1) %>%
      ungroup() %>%
      select(idauniq,
             baseline_cflisd = cflisd,
             cflisd_baseline_year = year,
             cflisd_baseline_wave = wave),
    by = "idauniq"
  )
write.csv(joint_baseline,
          file.path(OUT_DIR, "20_ELSA_JOINT_BASELINE.csv"), row.names = FALSE)

# Baseline covariates for reporting/model specification.
# These are taken from the same baseline record used to define
# baseline_frailty, so covariates and baseline state are aligned.
joint_baseline_covariates <- person_wave %>%
  semi_join(joint_ids, by = "idauniq") %>%
  filter(year == baseline_year, wave == baseline_wave) %>%
  select(
    idauniq, baseline_year, baseline_wave,
    age, age_final, hedibar, heeye, hedimbp,
    baseline_age, baseline_frailty
  ) %>%
  distinct(idauniq, .keep_all = TRUE)

write.csv(
  joint_baseline_covariates,
  file.path(OUT_DIR, "20A_ELSA_JOINT_BASELINE_COVARIATES.csv"),
  row.names = FALSE
)

############################################################
# 20. KAPLAN‑MEIER DATA
#
# Baseline exposure: baseline_frailty
# Event: death_year strictly after baseline_year
# Censoring: ralstcorey when available; otherwise last observed
# alive frailty year.
############################################################
last_observed_alive_year <- multistate_person_wave %>%
  group_by(idauniq) %>%
  summarise(last_observed_alive_year = max(year), .groups = "drop")

km_data_raw <- global_baseline %>%
  left_join(death_data, by = "idauniq", suffix = c("", "_deathdata")) %>%
  left_join(last_observed_alive_year, by = "idauniq") %>%
  mutate(
    last_alive_year_final = coalesce(last_alive_year, last_observed_alive_year),
    baseline_status = case_when(
      !is.na(death_year) & death_year < baseline_year ~ "death_before_baseline",
      !is.na(death_year) & death_year == baseline_year ~ "death_same_year_as_baseline",
      !is.na(death_year) & death_year > baseline_year ~ "death_after_baseline",
      TRUE ~ "no_known_death"
    ),
    death_event = if_else(baseline_status == "death_after_baseline", 1L, 0L),
    followup_end = case_when(
      death_event == 1L ~ death_year,
      !is.na(last_alive_year_final) ~ last_alive_year_final,
      TRUE ~ NA_real_
    ),
    survival_time = followup_end - baseline_year,
    ## PATCH‑3：标记生存时间等于0样本
    surv_time_zero = (survival_time == 0)
  )

# Save all baseline survival cases, including excluded ambiguity cases.
write.csv(km_data_raw,
          file.path(OUT_DIR, "21A_ELSA_KM_SURVIVAL_RAW_WITH_AUDIT.csv"), row.names = FALSE)

km_exclusion_audit <- km_data_raw %>%
  count(baseline_status, surv_time_zero, name = "n") %>%
  mutate(percent = n / sum(n) * 100)
write.csv(km_exclusion_audit,
          file.path(OUT_DIR, "21B_KM_baseline_exclusion_audit.csv"), row.names = FALSE)

############################################################
# 21. FINAL KM ANALYSIS DATA
############################################################
km_data <- km_data_raw %>%
  filter(
    baseline_status %in% c("death_after_baseline", "no_known_death"),
    baseline_frailty %in% c(0, 1, 2),
    !is.na(survival_time),
    survival_time >= 0
  ) %>%
  select(
    idauniq, baseline_year, baseline_wave, baseline_age, baseline_frailty,
    death_year, death_event, last_alive_year_final, followup_end,
    survival_time, baseline_status, surv_time_zero
  )

write.csv(km_data,
          file.path(OUT_DIR, "21_ELSA_KM_SURVIVAL.csv"), row.names = FALSE)

############################################################
# 22. KM ANALYSIS
############################################################
km_analysis_data <- km_data %>%
  mutate(
    baseline_frailty = factor(
      baseline_frailty,
      levels = c(0, 1, 2),
      labels = c("Non‑frail", "Pre‑frail", "Frail")
    )
  ) %>%
  filter(!is.na(idauniq), !is.na(survival_time), !is.na(death_event))

km_event_summary <- km_analysis_data %>%
  count(baseline_frailty, death_event, name = "n") %>%
  group_by(baseline_frailty) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()
write.csv(km_event_summary,
          file.path(OUT_DIR, "26_KM_event_summary.csv"), row.names = FALSE)

km_surv_object <- Surv(km_analysis_data$survival_time,
                       km_analysis_data$death_event)
km_fit <- survfit(km_surv_object ~ baseline_frailty,
                  data = km_analysis_data)
print(km_fit)

km_logrank <- survdiff(km_surv_object ~ baseline_frailty,
                       data = km_analysis_data)
logrank_p <- 1 - pchisq(km_logrank$chisq,
                        df = length(km_logrank$n) - 1)

logrank_result <- tibble(
  chisq = km_logrank$chisq,
  df = length(km_logrank$n) - 1,
  p_value = logrank_p
)
write.csv(logrank_result,
          file.path(OUT_DIR, "27_KM_logrank_test.csv"), row.names = FALSE)

km_median <- survminer::surv_median(km_fit)
write.csv(km_median,
          file.path(OUT_DIR, "28_KM_median_survival.csv"), row.names = FALSE)

############################################################
# 23. KM CURVE WITH RISK TABLE
############################################################
km_plot <- ggsurvplot(
  km_fit,
  data = km_analysis_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  censor = TRUE,
  surv.median.line = "none",
  xlab = "Years since baseline",
  ylab = "Overall survival probability",
  title = "Kaplan‑Meier Survival Curves by Baseline Frailty",
  legend.title = "Baseline frailty",
  legend.labs = c("Non‑frail", "Pre‑frail", "Frail"),
  risk.table.title = "Number at risk",
  risk.table.height = 0.25,
  ggtheme = theme_classic(base_size = 12)
)

print(km_plot)

ggsave(
  file.path(OUT_DIR, "Figure_KM_baseline_frailty.png"),
  km_plot$plot,
  width = 8, height = 6, dpi = 600
)

png(
  file.path(OUT_DIR, "Figure_KM_baseline_frailty_with_risktable.png"),
  width = 2400, height = 2100, res = 300
)
print(km_plot)
dev.off()

############################################################
# 24. FULL MISSINGNESS AUDIT
############################################################
missing_audit <- tibble(
  variable = names(person_wave),
  n_missing = sapply(person_wave, function(x) sum(is.na(x))),
  percent_missing = sapply(person_wave, function(x) mean(is.na(x)) * 100)
) %>%
  arrange(desc(percent_missing))
write.csv(missing_audit,
          file.path(OUT_DIR, "22_full_missingness_audit.csv"), row.names = FALSE)

############################################################
# 25. PARTICIPANT RETENTION AUDIT
############################################################
raw_n <- n_distinct(person_wave$idauniq)
longitudinal_n <- n_distinct(longitudinal_data$idauniq)
multistate_n <- n_distinct(multistate_person_wave$idauniq)
joint_n <- n_distinct(joint_ids$idauniq)
transition_n <- n_distinct(multistate_transitions$idauniq)
km_n <- n_distinct(km_data$idauniq)

participant_retention <- tibble(
  stage = c("Raw ELSA", "Longitudinal", "Multistate", "Joint overlap",
            "Multistate transition", "KM"),
  participants = c(raw_n, longitudinal_n, multistate_n, joint_n,
                   transition_n, km_n)
) %>%
  mutate(percent_of_raw = participants / raw_n * 100)
write.csv(participant_retention,
          file.path(OUT_DIR, "23_participant_retention.csv"), row.names = FALSE)

############################################################
# 26. FINAL SAMPLE AUDIT
############################################################
final_sample_audit <- tibble(
  dataset = c(
    "Raw", "Clean person‑wave", "Longitudinal", "Multistate person‑wave",
    "Observed transitions", "Death transitions", "All multistate transitions",
    "Joint longitudinal", "Joint multistate", "Joint transitions", "KM"
  ),
  rows = c(
    nrow(raw), nrow(person_wave), nrow(longitudinal_data),
    nrow(multistate_person_wave), nrow(observed_transitions),
    nrow(death_transitions), nrow(multistate_transitions),
    nrow(joint_longitudinal_data), nrow(joint_multistate_person_wave),
    nrow(joint_multistate_transitions), nrow(km_data)
  ),
  participants = c(
    n_distinct(raw$idauniq), n_distinct(person_wave$idauniq),
    n_distinct(longitudinal_data$idauniq), n_distinct(multistate_person_wave$idauniq),
    n_distinct(observed_transitions$idauniq), n_distinct(death_transitions$idauniq),
    n_distinct(multistate_transitions$idauniq), n_distinct(joint_longitudinal_data$idauniq),
    n_distinct(joint_multistate_person_wave$idauniq),
    n_distinct(joint_multistate_transitions$idauniq), n_distinct(km_data$idauniq)
  )
)
write.csv(final_sample_audit,
          file.path(OUT_DIR, "24_final_sample_audit.csv"), row.names = FALSE)

############################################################
# 27. TRANSITION SANITY CHECK
# PATCH‑5：控制台打印无效转移行数
############################################################
invalid_transitions <- multistate_transitions %>%
  filter(
    !(from %in% c(0, 1, 2)) |
      !(to %in% c(0, 1, 2, 3)) |
      is.na(start) | is.na(stop) | stop <= start
  )

write.csv(invalid_transitions,
          file.path(OUT_DIR, "25_INVALID_TRANSITIONS.csv"), row.names = FALSE)

n_invalid_trans <- nrow(invalid_transitions)
if(n_invalid_trans > 0){
  warning(sprintf("WARNING: %d invalid transitions found, check 25_INVALID_TRANSITIONS.csv", n_invalid_trans))
}

############################################################
# 28. ADDITIONAL DEATH CONSISTENCY AUDIT
############################################################
death_consistency_audit <- death_data %>%
  mutate(
    death_before_last_alive = !is.na(death_year) &
      !is.na(last_alive_year) & death_year < last_alive_year,
    death_same_year_last_alive = !is.na(death_year) &
      !is.na(last_alive_year) & death_year == last_alive_year,
    death_after_last_alive = !is.na(death_year) &
      !is.na(last_alive_year) & death_year > last_alive_year
  ) %>%
  summarise(
    n_with_death = sum(death_observed),
    n_death_before_last_alive = sum(death_before_last_alive, na.rm = TRUE),
    n_death_same_year_last_alive = sum(death_same_year_last_alive, na.rm = TRUE),
    n_death_after_last_alive = sum(death_after_last_alive, na.rm = TRUE)
  )
write.csv(death_consistency_audit,
          file.path(OUT_DIR, "29_death_consistency_audit.csv"), row.names = FALSE)

############################################################
# 29. FINAL CONSOLE SUMMARY
# PATCH‑6：打印survival_time=0样本数
############################################################
n_surv_zero <- sum(km_data$surv_time_zero, na.rm = TRUE)

cat("\n====================================================\n")
cat("ELSA JOINT MODEL DATA CONSTRUCTION COMPLETE\n")
cat("====================================================\n")
cat("Output directory:", OUT_DIR, "\n\n")
cat("Raw participants:", raw_n, "\n")
cat("Longitudinal participants:", longitudinal_n, "\n")
cat("Multistate participants:", multistate_n, "\n")
cat("Joint participants:", joint_n, "\n")
cat("Observed frailty transitions:", nrow(observed_transitions), "\n")
cat("Death transitions:", nrow(death_transitions), "\n")
cat("Total multistate transitions:", nrow(multistate_transitions), "\n")
cat("KM participants:", km_n, "\n")
cat("KM survival_time ==0 count:", n_surv_zero,"\n")
cat("KM log‑rank P:", format.pval(logrank_p, digits = 4), "\n")
cat("\nKey files to inspect:\n")
cat("1. 06A_baseline_death_audit.csv\n")
cat("2. 11_same_year_death_audit.csv\n")
cat("3. 11A_death_same_year_as_next_wave_audit.csv\n")
cat("4. 14_frailty_transition_matrix.csv\n")
cat("5. 15_death_transition_summary.csv\n")
cat("6. 21B_KM_baseline_exclusion_audit.csv\n")
cat("7. 23_participant_retention.csv\n")
cat("8. 25_INVALID_TRANSITIONS.csv\n")
cat("9. Figure_KM_baseline_frailty_with_risktable.png\n")
cat("====================================================\n")

############################################################
# END
############################################################
