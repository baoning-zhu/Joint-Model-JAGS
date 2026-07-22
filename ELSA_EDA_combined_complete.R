############################################################
## ELSA DATASET
## Complete exploratory data analysis reproduction
##
## Reproduce figures and results in elsa(3).pdf, pages 15–23
############################################################

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

set.seed(42)

############################################################
## 0. Package installation and loading
############################################################

required_packages <- c(
  "haven",
  "dplyr",
  "tidyr",
  "ggplot2",
  "survival",
  "survminer",
  "scales",
  "VIM",
  "naniar"
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

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(scales)
library(VIM)
library(naniar)

############################################################
## 1. File path and output directory
############################################################

## 修改成你电脑中数据文件的实际位置
data_file <- "new_dataset_noL 1.dta"

## 所有结果保存到该文件夹
output_dir <- "ELSA_EDA_results"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

############################################################
## 2. Read data
############################################################

if (!file.exists(data_file)) {
  stop(
    paste0(
      "Cannot find the data file: ",
      data_file,
      "\nPlease modify data_file to the correct path."
    )
  )
}

elsa_raw <- haven::read_dta(data_file)

cat("\n========================================\n")
cat("Original dataset dimensions\n")
cat("========================================\n")

cat("Rows:", nrow(elsa_raw), "\n")
cat("Columns:", ncol(elsa_raw), "\n")

############################################################
## 3. Keep only variables used in the EDA
############################################################

core_vars <- c(
  "idauniq",
  "wave",
  "year",
  "age",
  "count"
)

health_vars <- c(
  "hehelf",
  "hehelfi",
  "heeye",
  "heeyei",
  "hehear",
  "heheari"
)

disease_vars <- c(
  "hedimbp",
  "hediman",
  "hedimmi",
  "hedimhf",
  "hedimar",
  "hedimdi",
  "hedimst",
  "hediblu",
  "hedibas",
  "hedibar",
  "hedibos",
  "hedibca",
  "hedibpd",
  "hedibps",
  "hedibad",
  "hedibde"
)

mobility_vars <- c(
  "hemobwa",
  "hemobsi",
  "hemobch",
  "hemobcs",
  "hemobcl",
  "hemobst",
  "hemobre",
  "hemobpu",
  "hemobli",
  "hemobpi"
)

adl_vars <- c(
  "headldr",
  "headlwa",
  "headlba",
  "headlea",
  "headlbe",
  "headlwc"
)

iadl_vars <- c(
  "headlma",
  "headlpr",
  "headlsh",
  "headlph",
  "headlme",
  "headlho",
  "headlmo"
)

cognition_vars <- c(
  "cfdatd",
  "cfdatm",
  "cfdaty",
  "cfday",
  "cflisen",
  "cflisd",
  "cfliseni",
  "cflisdi"
)

cesd_vars <- c(
  "psceda",
  "pscedb",
  "pscedc",
  "pscedd",
  "pscede",
  "pscedf",
  "pscedg",
  "pscedh"
)

frailty_vars <- c(
  "denominator",
  "fraill",
  "frailgr",
  "nonfrail",
  "prefrail",
  "frail"
)

death_vars <- c(
  "raxseason",
  "raxyear",
  "radage",
  "radagef",
  "radtoivwm",
  "radtoivwy",
  "radtoivwf",
  "radloc_e",
  "radloc",
  "racod_e",
  "ragcod",
  "radexpec",
  "raddur",
  "raxtiwy",
  "ralstcorey",
  "radmarrp",
  "radlivnh"
)

all_used_vars <- unique(
  c(
    core_vars,
    health_vars,
    disease_vars,
    mobility_vars,
    adl_vars,
    iadl_vars,
    cognition_vars,
    cesd_vars,
    frailty_vars,
    death_vars
  )
)

existing_vars <- intersect(
  all_used_vars,
  names(elsa_raw)
)

missing_vars <- setdiff(
  all_used_vars,
  names(elsa_raw)
)

if (length(missing_vars) > 0) {
  cat(
    "\nVariables not found and therefore skipped:\n"
  )
  print(missing_vars)
}

elsa <- elsa_raw %>%
  select(all_of(existing_vars))

rm(elsa_raw)
gc()

############################################################
## 4. Convert labelled variables to numeric
############################################################

elsa <- elsa %>%
  mutate(
    across(
      where(haven::is.labelled),
      as.numeric
    )
  )

elsa <- elsa %>%
  arrange(
    idauniq,
    wave
  )

############################################################
## 5. Basic checks
############################################################

cat("\n========================================\n")
cat("Basic dataset information\n")
cat("========================================\n")

cat("Number of rows:", nrow(elsa), "\n")
cat("Number of columns:", ncol(elsa), "\n")

n_unique <- dplyr::n_distinct(
  elsa$idauniq,
  na.rm = TRUE
)

cat(
  "Number of unique individuals:",
  n_unique,
  "\n"
)

cat(
  "Available waves:",
  paste(
    sort(unique(elsa$wave)),
    collapse = ", "
  ),
  "\n"
)

duplicate_check <- elsa %>%
  count(
    idauniq,
    wave,
    name = "n"
  ) %>%
  filter(n > 1)

cat(
  "Duplicated ID-wave combinations:",
  nrow(duplicate_check),
  "\n"
)

############################################################
## 6. Wave summary
############################################################

wave_summary <- elsa %>%
  filter(!is.na(wave)) %>%
  count(
    wave,
    name = "Number_of_observations"
  ) %>%
  arrange(wave) %>%
  mutate(
    Percentage = 100 *
      Number_of_observations /
      sum(Number_of_observations)
  )

print(wave_summary)

write.csv(
  wave_summary,
  file.path(
    output_dir,
    "Table_01_wave_summary.csv"
  ),
  row.names = FALSE
)

person_wave_summary <- elsa %>%
  group_by(idauniq) %>%
  summarise(
    Number_of_waves = n_distinct(
      wave,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  count(
    Number_of_waves,
    name = "Number_of_individuals"
  ) %>%
  arrange(Number_of_waves)

write.csv(
  person_wave_summary,
  file.path(
    output_dir,
    "Table_02_participation_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## 7. Age distribution
############################################################

age_summary <- elsa %>%
  summarise(
    N = sum(!is.na(age)),
    Missing = sum(is.na(age)),
    Mean = mean(age, na.rm = TRUE),
    SD = sd(age, na.rm = TRUE),
    Minimum = min(age, na.rm = TRUE),
    Q1 = as.numeric(
      quantile(
        age,
        0.25,
        na.rm = TRUE
      )
    ),
    Median = median(
      age,
      na.rm = TRUE
    ),
    Q3 = as.numeric(
      quantile(
        age,
        0.75,
        na.rm = TRUE
      )
    ),
    Maximum = max(
      age,
      na.rm = TRUE
    )
  )

print(age_summary)

write.csv(
  age_summary,
  file.path(
    output_dir,
    "Table_03_age_summary.csv"
  ),
  row.names = FALSE
)

age_plot <- ggplot(
  elsa %>%
    filter(!is.na(age)),
  aes(x = age)
) +
  geom_histogram(
    binwidth = 5,
    boundary = 0,
    colour = "black",
    fill = "steelblue",
    alpha = 0.85
  ) +
  scale_x_continuous(
    breaks = seq(
      20,
      100,
      by = 20
    )
  ) +
  labs(
    title = "Age distribution",
    x = "Age",
    y = NULL
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(age_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_01_age_distribution.png"
  ),
  plot = age_plot,
  width = 6,
  height = 6,
  dpi = 300
)

############################################################
## 8. Create CES-D score
############################################################

available_cesd <- intersect(
  cesd_vars,
  names(elsa)
)

if (length(available_cesd) == 8) {
  
  elsa <- elsa %>%
    mutate(
      cesd_valid_items = rowSums(
        !is.na(
          across(
            all_of(available_cesd)
          )
        )
      ),
      cesd = if_else(
        cesd_valid_items == 8,
        rowSums(
          across(
            all_of(available_cesd)
          ),
          na.rm = FALSE
        ),
        NA_real_
      )
    )
  
} else {
  
  elsa$cesd <- NA_real_
  
  warning(
    "Some CES-D variables are absent. CES-D score was not calculated."
  )
}

############################################################
## 9. Create mobility, ADL and IADL scores
############################################################

available_mobility <- intersect(
  mobility_vars,
  names(elsa)
)

available_adl <- intersect(
  adl_vars,
  names(elsa)
)

available_iadl <- intersect(
  iadl_vars,
  names(elsa)
)

if (length(available_mobility) > 0) {
  
  elsa <- elsa %>%
    mutate(
      mobility_valid = rowSums(
        !is.na(
          across(
            all_of(available_mobility)
          )
        )
      ),
      mobility_score = if_else(
        mobility_valid ==
          length(available_mobility),
        rowSums(
          across(
            all_of(available_mobility)
          ),
          na.rm = FALSE
        ),
        NA_real_
      )
    )
  
} else {
  elsa$mobility_score <- NA_real_
}

if (length(available_adl) > 0) {
  
  elsa <- elsa %>%
    mutate(
      adl_valid = rowSums(
        !is.na(
          across(
            all_of(available_adl)
          )
        )
      ),
      adl_score = if_else(
        adl_valid ==
          length(available_adl),
        rowSums(
          across(
            all_of(available_adl)
          ),
          na.rm = FALSE
        ),
        NA_real_
      )
    )
  
} else {
  elsa$adl_score <- NA_real_
}

if (length(available_iadl) > 0) {
  
  elsa <- elsa %>%
    mutate(
      iadl_valid = rowSums(
        !is.na(
          across(
            all_of(available_iadl)
          )
        )
      ),
      iadl_score = if_else(
        iadl_valid ==
          length(available_iadl),
        rowSums(
          across(
            all_of(available_iadl)
          ),
          na.rm = FALSE
        ),
        NA_real_
      )
    )
  
} else {
  elsa$iadl_score <- NA_real_
}

############################################################
## 10. Frailty group coding
############################################################

elsa <- elsa %>%
  mutate(
    frailgr_factor = factor(
      frailgr,
      levels = c(0, 1, 2),
      labels = c(
        "Nonfrail",
        "Prefrail",
        "Frail"
      )
    )
  )

frailty_summary <- elsa %>%
  filter(!is.na(frailgr_factor)) %>%
  count(
    frailgr_factor,
    name = "Count"
  ) %>%
  mutate(
    Percent = 100 *
      Count /
      sum(Count)
  )

print(frailty_summary)

write.csv(
  frailty_summary,
  file.path(
    output_dir,
    "Table_04_frailty_group_summary.csv"
  ),
  row.names = FALSE
)


############################################################
## 11. Correct wave-based survival dataset
############################################################

first_valid <- function(
    x,
    lower = -Inf,
    upper = Inf
) {
  
  x <- as.numeric(x)
  
  valid <- x[
    !is.na(x) &
      is.finite(x) &
      x >= lower &
      x <= upper
  ]
  
  if (length(valid) == 0) {
    return(NA_real_)
  }
  
  valid[1]
}

last_valid <- function(
    x,
    lower = -Inf,
    upper = Inf
) {
  
  x <- as.numeric(x)
  
  valid <- x[
    !is.na(x) &
      is.finite(x) &
      x >= lower &
      x <= upper
  ]
  
  if (length(valid) == 0) {
    return(NA_real_)
  }
  
  valid[length(valid)]
}

############################################################
## 11.1 Identify actual interview rows
############################################################

elsa <- elsa %>%
  mutate(
    real_interview =
      !is.na(year) &
      year >= 1900 &
      year <= 2100
  )

cat("\nActual interview observations by wave:\n")

print(
  elsa %>%
    filter(real_interview) %>%
    count(wave)
)

############################################################
## 11.2 Build a wave-year lookup table
############################################################

## Derive the typical calendar year corresponding to each wave.
## For this dataset it should be close to:
## wave 1=2002, 2=2004, ..., 8=2016.

wave_year_lookup <- elsa %>%
  filter(
    real_interview,
    !is.na(wave)
  ) %>%
  group_by(wave) %>%
  summarise(
    survey_year = round(
      median(year, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  arrange(wave)

print(wave_year_lookup)

write.csv(
  wave_year_lookup,
  file.path(
    output_dir,
    "Table_wave_year_lookup.csv"
  ),
  row.names = FALSE
)

############################################################
## 11.3 Convert death year to study wave
############################################################

## Map a calendar death year to a wave.
##
## A death occurring between two survey years is assigned to
## the next scheduled wave. For example:
##
## death in 2003 -> wave 2
## death in 2005 -> wave 3

year_to_wave <- function(death_year) {
  
  if (
    length(death_year) == 0 ||
    is.na(death_year) ||
    !is.finite(death_year)
  ) {
    return(NA_real_)
  }
  
  eligible <- wave_year_lookup$wave[
    wave_year_lookup$survey_year >= death_year
  ]
  
  if (length(eligible) > 0) {
    return(min(eligible))
  }
  
  ## Death after the last scheduled interview
  max(wave_year_lookup$wave)
}

############################################################
## 11.4 One row per participant
############################################################

surv_person <- elsa %>%
  arrange(
    idauniq,
    wave
  ) %>%
  group_by(idauniq) %>%
  summarise(
    
    first_wave = if (
      any(real_interview)
    ) {
      min(
        wave[real_interview],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    
    last_wave = if (
      any(real_interview)
    ) {
      max(
        wave[real_interview],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    
    entry_year = first_valid(
      year[real_interview],
      lower = 1900,
      upper = 2100
    ),
    
    last_year = last_valid(
      year[real_interview],
      lower = 1900,
      upper = 2100
    ),
    
    death_year = first_valid(
      raxyear,
      lower = 1900,
      upper = 2100
    ),
    
    baseline_age = first_valid(
      age[real_interview],
      lower = 20,
      upper = 120
    ),
    
    baseline_fraill = first_valid(
      fraill[real_interview],
      lower = 0,
      upper = 1
    ),
    
    baseline_frailgr = first_valid(
      frailgr[real_interview],
      lower = 0,
      upper = 2
    ),
    
    baseline_cognition = first_valid(
      cflisen[real_interview],
      lower = 0,
      upper = 10
    ),
    
    .groups = "drop"
  )

############################################################
## 11.5 Create death wave and follow-up time
############################################################

surv_person$death_wave <- vapply(
  surv_person$death_year,
  year_to_wave,
  numeric(1)
)

surv_data_wave <- surv_person %>%
  mutate(
    
    died = as.integer(
      !is.na(death_year)
    ),
    
    ## Death event must occur no earlier than entry
    death_wave = case_when(
      is.na(death_wave) ~ NA_real_,
      death_wave < first_wave ~ first_wave,
      TRUE ~ death_wave
    ),
    
    ## Absolute end wave
    end_wave = case_when(
      died == 1 ~ death_wave,
      died == 0 ~ last_wave,
      TRUE ~ NA_real_
    ),
    
    ## Follow-up starts at zero for every participant
    survtime = end_wave - first_wave,
    
    ## Death occurring during the entry interval
    survtime = case_when(
      is.na(survtime) ~ NA_real_,
      survtime < 0 ~ NA_real_,
      died == 1 & survtime == 0 ~ 0.5,
      TRUE ~ survtime
    ),
    
    frailgr = factor(
      baseline_frailgr,
      levels = c(0, 1, 2),
      labels = c(
        "Nonfrail",
        "Prefrail",
        "Frail"
      )
    )
  ) %>%
  filter(
    !is.na(survtime),
    survtime >= 0,
    survtime <= 7
  )


############################################################
## 11.6 Create year-based survival data for formal analysis
############################################################

surv_data_year <- surv_person %>%
  mutate(

    died = as.integer(
      !is.na(death_year)
    ),

    ## End of follow-up:
    ## deceased participants use death year;
    ## censored participants use their last observed interview year.
    end_year = case_when(
      died == 1 ~ death_year,
      died == 0 ~ last_year,
      TRUE ~ NA_real_
    ),

    survtime = end_year - entry_year,

    ## Remove impossible follow-up values.
    ## A death in the entry year is assigned 0.5 years so that
    ## the event remains in the survival analysis.
    survtime = case_when(
      is.na(survtime) ~ NA_real_,
      survtime < 0 ~ NA_real_,
      died == 1 & survtime == 0 ~ 0.5,
      TRUE ~ survtime
    ),

    frailgr = factor(
      baseline_frailgr,
      levels = c(0, 1, 2),
      labels = c(
        "Nonfrail",
        "Prefrail",
        "Frail"
      )
    )
  ) %>%
  filter(
    !is.na(survtime),
    is.finite(survtime),
    survtime >= 0
  )


############################################################
## 12. Essential checks
############################################################

cat("\nNumber of participants:\n")
print(nrow(surv_data_wave))

cat("\nDeath status:\n")
print(
  table(
    surv_data_wave$died,
    useNA = "ifany"
  )
)

cat("\nDeath percentage:\n")
print(
  round(
    prop.table(
      table(surv_data_wave$died)
    ) * 100,
    3
  )
)

cat("\nFollow-up time distribution:\n")
print(
  table(
    surv_data_wave$survtime,
    useNA = "ifany"
  )
)

cat("\nDeath events by follow-up time:\n")
print(
  with(
    surv_data_wave,
    table(
      survtime,
      died
    )
  )
)

cat("\nFrailty group and death status:\n")
print(
  with(
    surv_data_wave,
    table(
      frailgr,
      died,
      useNA = "ifany"
    )
  )
)

############################################################
## 13. Overall Kaplan-Meier curve
############################################################

km_all_data <- surv_data_wave %>%
  filter(
    !is.na(survtime),
    !is.na(died),
    is.finite(survtime)
  )

km_fit <- survival::survfit(
  survival::Surv(
    time = survtime,
    event = died
  ) ~ 1,
  data = km_all_data
)

print(km_fit)
print(summary(km_fit))

max_followup <- ceiling(
  max(
    km_all_data$survtime,
    na.rm = TRUE
  )
)

km_all_plot <- survminer::ggsurvplot(
  fit = km_fit,
  data = km_all_data,
  
  conf.int = FALSE,
  risk.table = FALSE,
  
  censor = TRUE,
  censor.shape = "+",
  censor.size = 3,
  
  xlab = "Wave since entry",
  ylab = "Survival probability",
  
  xlim = c(0, max_followup),
  ylim = c(0, 1.05),
  
  break.time.by = 1,
  
  legend = "top",
  legend.title = "Strata",
  legend.labs = "All",
  
  palette = "#4C78A8",
  
  ggtheme = ggplot2::theme_classic(
    base_size = 12
  )
)

print(km_all_plot)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Figure_02_overall_KM_corrected.png"
  ),
  plot = km_all_plot$plot,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
## 14. KM curves by frailty group
############################################################

km_frail_data <- surv_data_wave %>%
  filter(
    !is.na(survtime),
    !is.na(died),
    !is.na(frailgr),
    is.finite(survtime)
  )

cat("\nParticipants by frailty group:\n")

print(
  table(
    km_frail_data$frailgr
  )
)

cat("\nDeaths by frailty group:\n")

print(
  with(
    km_frail_data,
    table(
      frailgr,
      died
    )
  )
)

km_frail <- survival::survfit(
  survival::Surv(
    time = survtime,
    event = died
  ) ~ frailgr,
  data = km_frail_data
)

print(km_frail)

############################################################
## Log-rank test
############################################################

logrank_frail <- survival::survdiff(
  survival::Surv(
    survtime,
    died
  ) ~ frailgr,
  data = km_frail_data
)

logrank_p <- 1 -
  pchisq(
    logrank_frail$chisq,
    df = length(logrank_frail$n) - 1
  )

cat(
  "\nLog-rank p-value:",
  format.pval(
    logrank_p,
    digits = 4,
    eps = 0.0001
  ),
  "\n"
)

############################################################
## Plot without theme_cleantable()
############################################################

km_frail_plot <- survminer::ggsurvplot(
  fit = km_frail,
  data = km_frail_data,
  
  conf.int = FALSE,
  
  risk.table = TRUE,
  risk.table.height = 0.28,
  risk.table.title = "Number at risk",
  
  risk.table.y.text = TRUE,
  risk.table.y.text.col = TRUE,
  
  censor = TRUE,
  censor.shape = "+",
  censor.size = 3,
  
  pval = TRUE,
  pval.method = FALSE,
  
  xlab = "Wave since entry",
  ylab = "Survival probability",
  
  xlim = c(
    0,
    max_followup
  ),
  
  ylim = c(
    0,
    1.05
  ),
  
  break.time.by = 1,
  
  legend = "top",
  legend.title = "Strata",
  legend.labs = c(
    "Nonfrail",
    "Prefrail",
    "Frail"
  ),
  
  palette = c(
    "#35B779",
    "#F5A623",
    "#F04E3E"
  ),
  
  ggtheme = ggplot2::theme_classic(
    base_size = 12
  ),
  
  tables.theme = ggplot2::theme_classic(
    base_size = 10
  )
)

print(km_frail_plot)

png(
  filename = file.path(
    output_dir,
    "Figure_03_KM_by_frailty_corrected.png"
  ),
  width = 2200,
  height = 1700,
  res = 250
)

print(km_frail_plot)

dev.off()


############################################################
## 16. Disease prevalence over waves
############################################################

available_disease_vars <- intersect(
  disease_vars,
  names(elsa)
)

disease_labels <- c(
  hedimbp = "Hypertension",
  hediman = "Angina",
  hedimmi = "Heart attack",
  hedimhf = "Heart failure",
  hedimar = "Abnormal rhythm",
  hedimdi = "Diabetes",
  hedimst = "Stroke",
  hediblu = "Lung disease",
  hedibas = "Asthma",
  hedibar = "Arthritis",
  hedibos = "Osteoporosis",
  hedibca = "Cancer",
  hedibpd = "Parkinson's",
  hedibps = "Psychiatric condition",
  hedibad = "Alzheimer's",
  hedibde = "Dementia"
)

disease_long <- elsa %>%
  select(
    wave,
    all_of(available_disease_vars)
  ) %>%
  pivot_longer(
    cols = all_of(
      available_disease_vars
    ),
    names_to = "Variable",
    values_to = "Disease"
  ) %>%
  filter(
    !is.na(wave),
    !is.na(Disease),
    Disease %in% c(0, 1)
  )

disease_prevalence <- disease_long %>%
  group_by(
    wave,
    Variable
  ) %>%
  summarise(
    Number_nonmissing = n(),
    Number_cases = sum(
      Disease == 1
    ),
    Prevalence = 100 *
      mean(
        Disease == 1
      ),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = unname(
      disease_labels[Variable]
    )
  )

print(disease_prevalence)

write.csv(
  disease_prevalence,
  file.path(
    output_dir,
    "Table_07_disease_prevalence_by_wave.csv"
  ),
  row.names = FALSE
)

############################################################
## 17. Disease prevalence plot matching PDF
############################################################

selected_conditions <- c(
  "Arthritis",
  "Asthma",
  "Diabetes",
  "Heart attack",
  "Lung disease"
)

disease_plot_data <- disease_prevalence %>%
  filter(
    Condition %in%
      selected_conditions
  )

disease_plot <- ggplot(
  disease_plot_data,
  aes(
    x = wave,
    y = Prevalence,
    colour = Condition,
    group = Condition
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2
  ) +
  scale_x_continuous(
    breaks = sort(
      unique(
        disease_plot_data$wave
      )
    )
  ) +
  labs(
    title = "Disease prevalence over waves",
    x = "Wave",
    y = "Prevalence (%)",
    colour = "Condition"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(disease_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_04_disease_prevalence_over_waves.png"
  ),
  plot = disease_plot,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
## 18. Death season table
############################################################

season_labels <- c(
  `1` = "Winter",
  `2` = "Spring",
  `3` = "Summer",
  `4` = "Autumn"
)

season_df <- elsa %>%
  filter(
    !is.na(raxseason),
    raxseason %in% 1:4
  ) %>%
  distinct(
    idauniq,
    .keep_all = TRUE
  ) %>%
  mutate(
    Season = factor(
      season_labels[
        as.character(raxseason)
      ],
      levels = c(
        "Winter",
        "Spring",
        "Summer",
        "Autumn"
      )
    )
  ) %>%
  count(
    Season,
    name = "Count"
  ) %>%
  tidyr::complete(
    Season,
    fill = list(
      Count = 0
    )
  ) %>%
  mutate(
    Percent = round(
      100 *
        Count /
        sum(Count),
      1
    )
  )

print(season_df)

write.csv(
  season_df,
  file.path(
    output_dir,
    "Table_08_death_season.csv"
  ),
  row.names = FALSE
)

season_plot <- ggplot(
  season_df,
  aes(
    x = Season,
    y = Count,
    fill = Season
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        Count,
        "\n",
        Percent,
        "%"
      )
    ),
    vjust = -0.2,
    size = 4
  ) +
  labs(
    title = "Season of death",
    x = NULL,
    y = "Count"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5
    )
  )

print(season_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_05_death_season.png"
  ),
  plot = season_plot,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
## 19. Person-level baseline missing-data dataset
############################################################

missing_variable_names <- c(
  "fraill",
  "age",
  "cflisen",
  "cflisd",
  "hedimdi",
  "hedibas",
  "hedibar",
  "hediblu",
  "hedimbp",
  "hedimmi",
  "hedimst",
  "hedibca"
)

available_missing_vars <- intersect(
  missing_variable_names,
  names(elsa)
)

## Return the first observed value for each variable
first_observed_or_na <- function(x) {

  x <- as.numeric(x)

  valid <- x[
    !is.na(x) &
      is.finite(x)
  ]

  if (length(valid) == 0) {
    return(NA_real_)
  }

  valid[1]
}

missing_person_with_id <- elsa %>%
  filter(real_interview) %>%
  arrange(
    idauniq,
    wave
  ) %>%
  group_by(idauniq) %>%
  summarise(
    across(
      all_of(available_missing_vars),
      first_observed_or_na
    ),
    .groups = "drop"
  )

## VIM::aggr must receive analysis variables only.
## Keeping idauniq in the data would create an unnamed column after renaming,
## which causes:
## "row names contain missing values".
missing_person <- missing_person_with_id %>%
  select(
    -idauniq
  )

## Stop early with a clear message if no variables are available.
if (ncol(missing_person) == 0) {
  stop(
    "None of the requested missing-data variables were found in the dataset."
  )
}

############################################################
## 20. Rename variables and summarize missingness
############################################################

display_names <- c(
  fraill = "Frailty Index",
  age = "Age",
  cflisen = "Immediate recall",
  cflisd = "Delayed recall",
  hedimdi = "Diabetes",
  hedibas = "Asthma",
  hedibar = "Arthritis",
  hediblu = "Lung disease",
  hedimbp = "Hypertension",
  hedimmi = "Heart attack",
  hedimst = "Stroke",
  hedibca = "Cancer"
)

new_missing_names <- unname(
  display_names[
    names(missing_person)
  ]
)

if (any(is.na(new_missing_names))) {
  stop(
    paste0(
      "The following variables do not have display names: ",
      paste(
        names(missing_person)[is.na(new_missing_names)],
        collapse = ", "
      )
    )
  )
}

names(missing_person) <- new_missing_names

missing_summary <- data.frame(
  Variable = names(missing_person),

  Missing_count = sapply(
    missing_person,
    function(x) {
      sum(is.na(x))
    }
  ),

  Missing_proportion = sapply(
    missing_person,
    function(x) {
      mean(is.na(x))
    }
  ),

  row.names = NULL
) %>%
  mutate(
    Missing_percent =
      100 * Missing_proportion
  )

print(missing_summary)

write.csv(
  missing_summary,
  file.path(
    output_dir,
    "Table_09_person_level_missingness.csv"
  ),
  row.names = FALSE
)

############################################################
## 21. Missingness bar plot
############################################################

missing_bar_plot <- ggplot(
  missing_summary,
  aes(
    x = reorder(
      Variable,
      Missing_proportion
    ),
    y = Missing_proportion
  )
) +
  geom_col(
    fill = "tomato"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    title = "Person-level baseline missing-data proportions",
    x = NULL,
    y = "Proportion missing"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    )
  )

print(missing_bar_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_06_person_level_missingness_barplot.png"
  ),
  plot = missing_bar_plot,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
## 22. Missing-data aggregation plot
############################################################

## This object is person-level and contains only about one row
## per participant, so sampling is not required.


## Validate the object before plotting.
stopifnot(
  is.data.frame(missing_person),
  nrow(missing_person) > 0,
  ncol(missing_person) > 0,
  !anyNA(names(missing_person)),
  all(nzchar(names(missing_person)))
)

cat(
  "\nMissing-pattern data dimensions:",
  nrow(missing_person),
  "rows x",
  ncol(missing_person),
  "variables\n"
)


png(
  filename = file.path(
    output_dir,
    "Figure_07_missing_pattern_person_level.png"
  ),
  width = 2400,
  height = 1600,
  res = 220
)

VIM::aggr(
  missing_person,

  col = c(
    "lightblue",
    "red"
  ),

  numbers = FALSE,

  sortVars = FALSE,

  sortCombs = TRUE,

  prop = TRUE,

  combined = TRUE,

  gap = 3,

  cex.axis = 0.68,

  ylab = c(
    "Proportion of missings",
    "Combinations"
  )
)

dev.off()

## Display in the R graphics device.
## try() prevents a graphics-device problem from stopping all later exports.
try(
  VIM::aggr(
    missing_person,

    col = c(
      "lightblue",
      "red"
    ),

    numbers = FALSE,

    sortVars = FALSE,

    sortCombs = TRUE,

    prop = TRUE,

    combined = TRUE,

    gap = 3,

    cex.axis = 0.68,

    ylab = c(
      "Proportion of missings",
      "Combinations"
    )
  ),
  silent = FALSE
)

############################################################
## Optional: longitudinal-row missing-data pattern
############################################################

missing_long <- elsa %>%
  filter(real_interview) %>%
  select(
    all_of(
      available_missing_vars
    )
  )

new_missing_long_names <- unname(
  display_names[
    names(missing_long)
  ]
)

if (any(is.na(new_missing_long_names))) {
  stop(
    paste0(
      "The following longitudinal missing-data variables do not have display names: ",
      paste(
        names(missing_long)[is.na(new_missing_long_names)],
        collapse = ", "
      )
    )
  )
}

names(missing_long) <- new_missing_long_names

png(
  filename = file.path(
    output_dir,
    "Figure_07B_missing_pattern_longitudinal.png"
  ),
  width = 2400,
  height = 1600,
  res = 220
)

VIM::aggr(
  missing_long,
  col = c(
    "lightblue",
    "red"
  ),
  numbers = FALSE,
  sortVars = FALSE,
  sortCombs = TRUE,
  prop = TRUE,
  combined = TRUE,
  gap = 3,
  cex.axis = 0.68,
  ylab = c(
    "Proportion of missings",
    "Combinations"
  )
)

dev.off()

gc()

############################################################
## 23. Candidate longitudinal outcomes summary
############################################################

candidate_vars <- c(
  "fraill",
  "cflisen",
  "cflisd",
  "hehelf",
  "hedimdi",
  "hedibas",
  "hedibar",
  "hediblu",
  "cesd"
)

candidate_vars <- intersect(
  candidate_vars,
  names(elsa)
)

candidate_summary <- lapply(
  candidate_vars,
  function(variable_name) {
    
    x <- elsa[[variable_name]]
    
    data.frame(
      Variable = variable_name,
      N = sum(!is.na(x)),
      Missing = sum(is.na(x)),
      Missing_percent =
        100 * mean(is.na(x)),
      Mean = mean(
        x,
        na.rm = TRUE
      ),
      SD = sd(
        x,
        na.rm = TRUE
      ),
      Median = median(
        x,
        na.rm = TRUE
      ),
      Minimum = min(
        x,
        na.rm = TRUE
      ),
      Maximum = max(
        x,
        na.rm = TRUE
      )
    )
  }
) %>%
  bind_rows()

print(candidate_summary)

write.csv(
  candidate_summary,
  file.path(
    output_dir,
    "Table_10_candidate_longitudinal_variables.csv"
  ),
  row.names = FALSE
)


############################################################
## 23A. Additional comprehensive EDA modules
############################################################

## These analyses supplement the original reproduction code.
## They distinguish person-wave summaries from participant-level summaries.

############################################################
## 23A.1 Complete longitudinal continuous-variable summary
############################################################

safe_summary <- function(data, variables) {

  variables <- intersect(
    variables,
    names(data)
  )

  if (length(variables) == 0) {
    return(data.frame())
  }

  lapply(
    variables,
    function(v) {

      x <- suppressWarnings(
        as.numeric(data[[v]])
      )

      valid <- x[
        !is.na(x) &
          is.finite(x)
      ]

      if (length(valid) == 0) {
        return(
          data.frame(
            Variable = v,
            N = 0,
            Missing = length(x),
            Missing_percent = 100,
            Mean = NA_real_,
            SD = NA_real_,
            Median = NA_real_,
            Q1 = NA_real_,
            Q3 = NA_real_,
            Minimum = NA_real_,
            Maximum = NA_real_
          )
        )
      }

      data.frame(
        Variable = v,
        N = length(valid),
        Missing = sum(
          is.na(x) |
            !is.finite(x)
        ),
        Missing_percent =
          100 *
            mean(
              is.na(x) |
                !is.finite(x)
            ),
        Mean = mean(valid),
        SD = if (
          length(valid) > 1
        ) {
          sd(valid)
        } else {
          NA_real_
        },
        Median = median(valid),
        Q1 = as.numeric(
          quantile(
            valid,
            0.25
          )
        ),
        Q3 = as.numeric(
          quantile(
            valid,
            0.75
          )
        ),
        Minimum = min(valid),
        Maximum = max(valid)
      )
    }
  ) %>%
    bind_rows()
}

complete_continuous_vars <- c(
  "age",
  "hehelf",
  "heeye",
  "hehear",
  "cflisen",
  "cflisd",
  "fraill",
  "cesd",
  "mobility_score",
  "adl_score",
  "iadl_score"
)

complete_longitudinal_summary <- safe_summary(
  elsa %>%
    filter(real_interview),
  complete_continuous_vars
)

print(complete_longitudinal_summary)

write.csv(
  complete_longitudinal_summary,
  file.path(
    output_dir,
    "Table_11_complete_longitudinal_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## 23A.2 One-row-per-participant baseline dataset
############################################################

baseline_data <- elsa %>%
  filter(real_interview) %>%
  arrange(
    idauniq,
    wave
  ) %>%
  group_by(idauniq) %>%
  summarise(
    baseline_wave = first_valid(
      wave,
      lower = 1
    ),
    baseline_year = first_valid(
      year,
      lower = 1900,
      upper = 2100
    ),
    age = first_valid(
      age,
      lower = 20,
      upper = 120
    ),
    hehelf = if (
      "hehelf" %in% names(cur_data())
    ) {
      first_valid(hehelf)
    } else {
      NA_real_
    },
    heeye = if (
      "heeye" %in% names(cur_data())
    ) {
      first_valid(heeye)
    } else {
      NA_real_
    },
    hehear = if (
      "hehear" %in% names(cur_data())
    ) {
      first_valid(hehear)
    } else {
      NA_real_
    },
    cflisen = if (
      "cflisen" %in% names(cur_data())
    ) {
      first_valid(
        cflisen,
        0,
        10
      )
    } else {
      NA_real_
    },
    cflisd = if (
      "cflisd" %in% names(cur_data())
    ) {
      first_valid(
        cflisd,
        0,
        10
      )
    } else {
      NA_real_
    },
    fraill = if (
      "fraill" %in% names(cur_data())
    ) {
      first_valid(
        fraill,
        0,
        1
      )
    } else {
      NA_real_
    },
    frailgr = if (
      "frailgr" %in% names(cur_data())
    ) {
      first_valid(
        frailgr,
        0,
        2
      )
    } else {
      NA_real_
    },
    cesd = if (
      "cesd" %in% names(cur_data())
    ) {
      first_valid(
        cesd,
        0,
        8
      )
    } else {
      NA_real_
    },
    mobility_score = if (
      "mobility_score" %in% names(cur_data())
    ) {
      first_valid(
        mobility_score
      )
    } else {
      NA_real_
    },
    adl_score = if (
      "adl_score" %in% names(cur_data())
    ) {
      first_valid(
        adl_score
      )
    } else {
      NA_real_
    },
    iadl_score = if (
      "iadl_score" %in% names(cur_data())
    ) {
      first_valid(
        iadl_score
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    frailgr_factor = factor(
      frailgr,
      levels = c(
        0,
        1,
        2
      ),
      labels = c(
        "Nonfrail",
        "Prefrail",
        "Frail"
      )
    )
  )

baseline_summary <- safe_summary(
  baseline_data,
  complete_continuous_vars
)

print(baseline_summary)

write.csv(
  baseline_summary,
  file.path(
    output_dir,
    "Table_12_baseline_summary.csv"
  ),
  row.names = FALSE
)

baseline_age_plot <- ggplot(
  baseline_data %>%
    filter(
      !is.na(age),
      is.finite(age)
    ),
  aes(x = age)
) +
  geom_histogram(
    binwidth = 5,
    boundary = 0,
    colour = "black",
    fill = "steelblue",
    alpha = 0.85
  ) +
  labs(
    title = "Baseline age distribution",
    subtitle = "Analysis unit: participant",
    x = "Baseline age",
    y = "Number of participants"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

print(baseline_age_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_08_baseline_age_distribution.png"
  ),
  plot = baseline_age_plot,
  width = 7,
  height = 5,
  dpi = 300
)

baseline_frailty_summary <- baseline_data %>%
  filter(
    !is.na(frailgr_factor)
  ) %>%
  count(
    frailgr_factor,
    name = "Count"
  ) %>%
  mutate(
    Percent =
      100 *
        Count /
        sum(Count)
  )

write.csv(
  baseline_frailty_summary,
  file.path(
    output_dir,
    "Table_13_baseline_frailty_groups.csv"
  ),
  row.names = FALSE
)

baseline_frailty_plot <- ggplot(
  baseline_frailty_summary,
  aes(
    x = frailgr_factor,
    y = Count,
    fill = frailgr_factor
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        Count,
        "\n",
        round(
          Percent,
          1
        ),
        "%"
      )
    ),
    vjust = -0.2
  ) +
  labs(
    title = "Baseline frailty-group distribution",
    subtitle = "Analysis unit: participant",
    x = "Frailty group",
    y = "Number of participants"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

print(baseline_frailty_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_09_baseline_frailty_groups.png"
  ),
  plot = baseline_frailty_plot,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
## 23A.3 Binary-variable frequency summary
############################################################

binary_vars <- intersect(
  c(
    "hehelfi",
    "heeyei",
    "heheari",
    disease_vars
  ),
  names(elsa)
)

binary_summary <- lapply(
  binary_vars,
  function(v) {

    x <- as.numeric(
      elsa[[v]]
    )

    valid_binary <-
      !is.na(x) &
        x %in% c(0, 1)

    data.frame(
      Variable = v,
      Valid_binary_N =
        sum(valid_binary),
      No_0 =
        sum(
          x == 0,
          na.rm = TRUE
        ),
      Yes_1 =
        sum(
          x == 1,
          na.rm = TRUE
        ),
      Prevalence_percent =
        if (
          sum(valid_binary) > 0
        ) {
          100 *
            mean(
              x[valid_binary] == 1
            )
        } else {
          NA_real_
        },
      Missing =
        sum(is.na(x)),
      Missing_percent =
        100 *
          mean(is.na(x)),
      Invalid_nonmissing =
        sum(
          !is.na(x) &
            !x %in% c(0, 1)
        )
    )
  }
) %>%
  bind_rows()

print(binary_summary)

write.csv(
  binary_summary,
  file.path(
    output_dir,
    "Table_14_binary_variable_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## 23A.4 Frailty trajectories over waves
############################################################

if ("fraill" %in% names(elsa)) {

  frailty_wave_summary <- elsa %>%
    filter(
      real_interview,
      !is.na(fraill),
      is.finite(fraill)
    ) %>%
    group_by(wave) %>%
    summarise(
      N = n(),
      Mean_frailty =
        mean(fraill),
      SD_frailty =
        sd(fraill),
      SE_frailty =
        SD_frailty /
          sqrt(N),
      Lower = pmax(
        Mean_frailty -
          1.96 *
            SE_frailty,
        0
      ),
      Upper = pmin(
        Mean_frailty +
          1.96 *
            SE_frailty,
        1
      ),
      .groups = "drop"
    )

  write.csv(
    frailty_wave_summary,
    file.path(
      output_dir,
      "Table_15_frailty_index_by_wave.csv"
    ),
    row.names = FALSE
  )

  frailty_trajectory_plot <- ggplot(
    frailty_wave_summary,
    aes(
      x = wave,
      y = Mean_frailty
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Lower,
        ymax = Upper
      ),
      alpha = 0.2
    ) +
    geom_line(
      linewidth = 0.9
    ) +
    geom_point(
      size = 2
    ) +
    scale_x_continuous(
      breaks = sort(
        unique(
          frailty_wave_summary$wave
        )
      )
    ) +
    labs(
      title = "Mean frailty index over waves",
      subtitle = "Ribbon shows approximate 95% confidence intervals",
      x = "Wave",
      y = "Mean frailty index"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )

  print(frailty_trajectory_plot)

  ggsave(
    filename = file.path(
      output_dir,
      "Figure_10_frailty_index_trajectory.png"
    ),
    plot = frailty_trajectory_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
}

frailty_group_wave <- elsa %>%
  filter(
    real_interview,
    !is.na(frailgr_factor)
  ) %>%
  count(
    wave,
    frailgr_factor,
    name = "Count"
  ) %>%
  group_by(wave) %>%
  mutate(
    Percent =
      100 *
        Count /
        sum(Count)
  ) %>%
  ungroup()

write.csv(
  frailty_group_wave,
  file.path(
    output_dir,
    "Table_16_frailty_group_percent_by_wave.csv"
  ),
  row.names = FALSE
)

frailty_group_wave_plot <- ggplot(
  frailty_group_wave,
  aes(
    x = wave,
    y = Percent,
    colour = frailgr_factor,
    group = frailgr_factor
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2
  ) +
  scale_x_continuous(
    breaks = sort(
      unique(
        frailty_group_wave$wave
      )
    )
  ) +
  labs(
    title = "Frailty-group proportions over waves",
    x = "Wave",
    y = "Percentage",
    colour = "Frailty group"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

print(frailty_group_wave_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_11_frailty_group_percent_by_wave.png"
  ),
  plot = frailty_group_wave_plot,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
## 23A.5 Cognitive-score distributions
############################################################

available_memory_vars <- intersect(
  c(
    "cflisen",
    "cflisd"
  ),
  names(elsa)
)

if (length(available_memory_vars) > 0) {

  cognition_long <- elsa %>%
    filter(real_interview) %>%
    select(
      all_of(
        available_memory_vars
      )
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Measure",
      values_to = "Score"
    ) %>%
    mutate(
      Measure = recode(
        Measure,
        cflisen =
          "Immediate recall",
        cflisd =
          "Delayed recall"
      )
    )

  cognition_plot <- ggplot(
    cognition_long %>%
      filter(
        !is.na(Score),
        is.finite(Score)
      ),
    aes(
      x = Score,
      fill = Measure
    )
  ) +
    geom_histogram(
      binwidth = 1,
      boundary = -0.5,
      colour = "black",
      alpha = 0.75
    ) +
    facet_wrap(
      ~ Measure,
      ncol = 1,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = 0:10
    ) +
    labs(
      title = "Distribution of memory-recall scores",
      subtitle = "Analysis unit: person-wave observation",
      x = "Number of words recalled",
      y = "Frequency"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )

  print(cognition_plot)

  ggsave(
    filename = file.path(
      output_dir,
      "Figure_12_cognitive_distributions.png"
    ),
    plot = cognition_plot,
    width = 7,
    height = 7,
    dpi = 300
  )
}

############################################################
## 23A.6 Plot prevalence of all chronic diseases
############################################################

if (
  exists("disease_prevalence") &&
    nrow(disease_prevalence) > 0
) {

  all_disease_plot <- ggplot(
    disease_prevalence,
    aes(
      x = wave,
      y = Prevalence,
      colour = Condition,
      group = Condition
    )
  ) +
    geom_line(
      linewidth = 0.65
    ) +
    geom_point(
      size = 1.3
    ) +
    scale_x_continuous(
      breaks = sort(
        unique(
          disease_prevalence$wave
        )
      )
    ) +
    labs(
      title = "Prevalence of all chronic conditions over waves",
      x = "Wave",
      y = "Prevalence (%)",
      colour = "Condition"
    ) +
    theme_classic(
      base_size = 10
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "right"
    )

  print(all_disease_plot)

  ggsave(
    filename = file.path(
      output_dir,
      "Figure_13_all_disease_prevalence.png"
    ),
    plot = all_disease_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
}

############################################################
## 23A.7 Year-based Kaplan-Meier curves
############################################################

km_year_data <- surv_data_year %>%
  filter(
    !is.na(survtime),
    !is.na(died),
    is.finite(survtime)
  )

if (nrow(km_year_data) > 0) {

  km_year_fit <- survival::survfit(
    survival::Surv(
      survtime,
      died
    ) ~ 1,
    data = km_year_data
  )

  max_year_followup <- ceiling(
    max(
      km_year_data$survtime,
      na.rm = TRUE
    )
  )

  km_year_plot <- survminer::ggsurvplot(
    fit = km_year_fit,
    data = km_year_data,
    conf.int = TRUE,
    risk.table = TRUE,
    censor = TRUE,
    xlab = "Years since entry",
    ylab = "Survival probability",
    title = "Overall Kaplan-Meier curve: year time scale",
    xlim = c(
      0,
      max_year_followup
    ),
    ylim = c(
      0,
      1.05
    ),
    break.time.by = 2,
    risk.table.height = 0.25,
    ggtheme = theme_classic(
      base_size = 12
    ),
    tables.theme = theme_classic(
      base_size = 9
    )
  )

  print(km_year_plot)

  png(
    filename = file.path(
      output_dir,
      "Figure_14_overall_KM_year.png"
    ),
    width = 2200,
    height = 1700,
    res = 250
  )

  print(km_year_plot)

  dev.off()
}

km_year_frail_data <- surv_data_year %>%
  filter(
    !is.na(survtime),
    !is.na(died),
    !is.na(frailgr),
    is.finite(survtime)
  )

if (
  nrow(km_year_frail_data) > 0 &&
    n_distinct(
      km_year_frail_data$frailgr
    ) >= 2
) {

  km_year_frail <- survival::survfit(
    survival::Surv(
      survtime,
      died
    ) ~ frailgr,
    data = km_year_frail_data
  )

  logrank_year <- survival::survdiff(
    survival::Surv(
      survtime,
      died
    ) ~ frailgr,
    data = km_year_frail_data
  )

  logrank_year_p <- 1 -
    pchisq(
      logrank_year$chisq,
      df =
        length(logrank_year$n) - 1
    )

  logrank_year_table <- data.frame(
    Time_scale = "Year",
    Chi_square =
      logrank_year$chisq,
    Degrees_of_freedom =
      length(logrank_year$n) - 1,
    P_value =
      logrank_year_p
  )

  write.csv(
    logrank_year_table,
    file.path(
      output_dir,
      "Table_17_logrank_year.csv"
    ),
    row.names = FALSE
  )

  km_year_frail_plot <- survminer::ggsurvplot(
    fit = km_year_frail,
    data = km_year_frail_data,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.28,
    censor = TRUE,
    pval = TRUE,
    xlab = "Years since entry",
    ylab = "Survival probability",
    title = "Kaplan-Meier curves by baseline frailty group",
    xlim = c(
      0,
      ceiling(
        max(
          km_year_frail_data$survtime,
          na.rm = TRUE
        )
      )
    ),
    ylim = c(
      0,
      1.05
    ),
    break.time.by = 2,
    legend = "top",
    legend.title =
      "Baseline frailty group",
    legend.labs = c(
      "Nonfrail",
      "Prefrail",
      "Frail"
    ),
    ggtheme = theme_classic(
      base_size = 12
    ),
    tables.theme = theme_classic(
      base_size = 9
    )
  )

  print(km_year_frail_plot)

  png(
    filename = file.path(
      output_dir,
      "Figure_15_KM_by_frailty_year.png"
    ),
    width = 2200,
    height = 1750,
    res = 250
  )

  print(km_year_frail_plot)

  dev.off()
}

############################################################
## 23A.8 Missingness by wave
############################################################

missing_wave_vars <- intersect(
  c(
    "fraill",
    "age",
    "hehelf",
    "cflisen",
    "cflisd",
    "hedimdi",
    "hedibas",
    "hedibar",
    "hediblu",
    "hedimbp",
    "hedimmi",
    "hedimst",
    "hedibca",
    "cesd"
  ),
  names(elsa)
)

missing_by_wave <- elsa %>%
  filter(real_interview) %>%
  select(
    wave,
    all_of(
      missing_wave_vars
    )
  ) %>%
  pivot_longer(
    cols = -wave,
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  group_by(
    wave,
    Variable
  ) %>%
  summarise(
    N = n(),
    Missing_count =
      sum(is.na(Value)),
    Missing_percent =
      100 *
        mean(is.na(Value)),
    .groups = "drop"
  )

print(missing_by_wave)

write.csv(
  missing_by_wave,
  file.path(
    output_dir,
    "Table_18_missingness_by_wave.csv"
  ),
  row.names = FALSE
)

missing_by_wave_plot <- ggplot(
  missing_by_wave,
  aes(
    x = wave,
    y = Missing_percent,
    colour = Variable,
    group = Variable
  )
) +
  geom_line(
    linewidth = 0.75
  ) +
  geom_point(
    size = 1.5
  ) +
  scale_x_continuous(
    breaks = sort(
      unique(
        missing_by_wave$wave
      )
    )
  ) +
  labs(
    title = "Missingness by wave",
    x = "Wave",
    y = "Missing values (%)",
    colour = "Variable"
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid.minor =
      element_blank(),
    legend.position = "right"
  )

print(missing_by_wave_plot)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_16_missingness_by_wave.png"
  ),
  plot = missing_by_wave_plot,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
## 23A.9 Longitudinal missingness bar plot and heatmap
############################################################

long_missing_summary <- data.frame(
  Variable = names(missing_long),
  Missing_count = sapply(
    missing_long,
    function(x) {
      sum(is.na(x))
    }
  ),
  Missing_proportion = sapply(
    missing_long,
    function(x) {
      mean(is.na(x))
    }
  ),
  row.names = NULL
) %>%
  mutate(
    Missing_percent =
      100 *
        Missing_proportion
  ) %>%
  arrange(
    desc(
      Missing_proportion
    )
  )

write.csv(
  long_missing_summary,
  file.path(
    output_dir,
    "Table_19_longitudinal_missingness.csv"
  ),
  row.names = FALSE
)

long_missing_bar <- ggplot(
  long_missing_summary,
  aes(
    x = reorder(
      Variable,
      Missing_proportion
    ),
    y = Missing_proportion
  )
) +
  geom_col(
    fill = "steelblue"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    expand =
      expansion(
        mult = c(
          0,
          0.08
        )
      )
  ) +
  labs(
    title = "Longitudinal missing-data proportions",
    subtitle = "Analysis unit: person-wave observation",
    x = NULL,
    y = "Proportion missing"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

print(long_missing_bar)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_17_longitudinal_missingness.png"
  ),
  plot = long_missing_bar,
  width = 8,
  height = 6,
  dpi = 300
)

set.seed(42)

missing_heatmap_sample <- missing_long %>%
  slice_sample(
    n = min(
      3000,
      nrow(missing_long)
    )
  )

missing_heatmap <- naniar::vis_miss(
  missing_heatmap_sample,
  cluster = TRUE,
  sort_miss = TRUE,
  show_perc = TRUE,
  warn_large_data = FALSE
) +
  labs(
    title = "Longitudinal missing-data heatmap",
    subtitle = "Reproducible sample of up to 3,000 records"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(missing_heatmap)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_18_longitudinal_missing_heatmap.png"
  ),
  plot = missing_heatmap,
  width = 10,
  height = 6,
  dpi = 300
)

############################################################
## 23A.10 Participant-flow summary
############################################################

participant_sets <- elsa %>%
  filter(real_interview) %>%
  distinct(
    idauniq,
    wave
  )

available_waves <- sort(
  unique(
    participant_sets$wave
  )
)

participant_flow <- lapply(
  available_waves,
  function(w) {

    current_ids <- participant_sets %>%
      filter(wave == w) %>%
      pull(idauniq) %>%
      unique()

    previous_ids <- participant_sets %>%
      filter(wave < w) %>%
      pull(idauniq) %>%
      unique()

    future_ids <- participant_sets %>%
      filter(wave > w) %>%
      pull(idauniq) %>%
      unique()

    data.frame(
      Wave = w,
      Participants =
        length(current_ids),
      New_participants =
        sum(
          !current_ids %in%
            previous_ids
        ),
      Previously_seen =
        sum(
          current_ids %in%
            previous_ids
        ),
      Observed_again_later =
        sum(
          current_ids %in%
            future_ids
        ),
      Not_observed_later =
        sum(
          !current_ids %in%
            future_ids
        )
    )
  }
) %>%
  bind_rows()

print(participant_flow)

write.csv(
  participant_flow,
  file.path(
    output_dir,
    "Table_20_participant_flow_by_wave.csv"
  ),
  row.names = FALSE
)

############################################################
## 23A.11 Manual survival-construction check
############################################################

set.seed(42)

check_ids <- sample(
  unique(
    elsa$idauniq[
      !is.na(
        elsa$idauniq
      )
    ]
  ),
  size = min(
    20,
    n_distinct(
      elsa$idauniq,
      na.rm = TRUE
    )
  )
)

manual_check_long <- elsa %>%
  filter(
    idauniq %in%
      check_ids
  ) %>%
  select(
    any_of(
      c(
        "idauniq",
        "wave",
        "year",
        "age",
        "fraill",
        "frailgr",
        "raxyear",
        "raxseason"
      )
    )
  ) %>%
  arrange(
    idauniq,
    wave
  )

manual_check_survival <- surv_person %>%
  filter(
    idauniq %in%
      check_ids
  ) %>%
  arrange(idauniq)

write.csv(
  manual_check_long,
  file.path(
    output_dir,
    "Table_21_manual_check_long.csv"
  ),
  row.names = FALSE
)

write.csv(
  manual_check_survival,
  file.path(
    output_dir,
    "Table_22_manual_check_survival.csv"
  ),
  row.names = FALSE
)


############################################################
## 24. Save cleaned data
############################################################

saveRDS(
  elsa,
  file.path(
    output_dir,
    "ELSA_EDA_cleaned_long_data.rds"
  )
)

saveRDS(
  baseline_data,
  file.path(
    output_dir,
    "ELSA_EDA_baseline_data.rds"
  )
)

saveRDS(
  surv_data_wave,
  file.path(
    output_dir,
    "ELSA_EDA_wave_survival_data.rds"
  )
)

saveRDS(
  surv_data_year,
  file.path(
    output_dir,
    "ELSA_EDA_year_survival_data.rds"
  )
)

saveRDS(
  missing_person,
  file.path(
    output_dir,
    "ELSA_EDA_person_missing_data.rds"
  )
)


saveRDS(
  missing_person_with_id,
  file.path(
    output_dir,
    "ELSA_EDA_person_missing_data_with_id.rds"
  )
)

############################################################
## 25. Save session information
############################################################

sink(
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

print(sessionInfo())

sink()

############################################################
## 26. Final report to R console
############################################################

cat("\n")
cat("============================================================\n")
cat("ELSA EDA COMPLETED\n")
cat("============================================================\n")

cat(
  "Rows in longitudinal dataset:",
  nrow(elsa),
  "\n"
)

cat(
  "Unique individuals:",
  n_distinct(elsa$idauniq),
  "\n"
)

cat(
  "Available waves:",
  paste(
    sort(
      unique(
        elsa$wave
      )
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Number of deaths:",
  sum(
    surv_data_year$died == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Death percentage:",
  round(
    100 *
      mean(
        surv_data_year$died == 1,
        na.rm = TRUE
      ),
    2
  ),
  "%\n"
)

cat(
  "Mean age:",
  round(
    mean(
      elsa$age,
      na.rm = TRUE
    ),
    2
  ),
  "\n"
)

cat(
  "Mean frailty index:",
  round(
    mean(
      elsa$fraill,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Results saved in:",
  normalizePath(output_dir),
  "\n"
)

cat("============================================================\n")