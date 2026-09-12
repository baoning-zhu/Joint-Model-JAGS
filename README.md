# Additional Materials for the Extended Research Project

## Bayesian Joint Modelling of Longitudinal Outcomes and Multi-State Survival Processes

This directory contains the computational materials accompanying the Extended Research Project on Bayesian joint modelling of longitudinal outcomes and multi-state survival processes.

The materials include the R and JAGS code used for data preprocessing, simulation studies, empirical model fitting, posterior inference, convergence assessment, and generation of the main numerical and graphical results.

The empirical application uses data from the **English Longitudinal Study of Ageing (ELSA)**. The original ELSA data are not included in this repository because access and redistribution are subject to the relevant ELSA data-use conditions.

---

## 1. Repository structure

```text
ERP-Additional-Materials/
│
├── README.md
│
├── Code/
│   ├── data_preprocessing.R
│   ├── ELSA_joint_model application.R
│   ├── joint_shared_RE_model.jags
│   ├── longitudinal.R
│   ├── multistate_simulation_moderate_complete.R
│   └── simulation_study.R
│
├── Results/
│   ├── 00_observed_transition_counts.csv
│   ├── 01_subject_baseline_covariates_8cov.csv
│   ├── 02_transition_retention_check.csv
│   ├── 03_covariate_centering_constants.csv
│   ├── 04_longitudinal_model_data.csv
│   ├── 05_longitudinal_subject_audit.csv
│   ├── 06_crude_transition_rates.csv
│   ├── 07_convergence_Rhat_ESS.csv
│   ├── 10_longitudinal_core_results.csv
│   ├── 11_longitudinal_8cov_results.csv
│   ├── 12_multistate_8cov_HR_results.csv
│   ├── 13_joint_association_alpha_results.csv
│   └── 14_baseline_transition_intensities.csv
│
├── Figures/
│   ├── trace_density.pdf
│   ├── gelman_plots.pdf
│   └── multistate_HR_forest_plot.png
│
└── Simulation/
    ├── simulation_results.csv
    ├── simulation_summary.csv
    ├── joint_convergence_N500.csv
    ├── longitudinal_simulation_results/
    └── multistate_simulation_results_moderate/
```

The numbering of files in the `Results/` directory follows the numbering used in the original analysis pipeline. Intermediate files that are not required for interpretation or reproduction of the reported results are not included.

---

## 2. Software requirements

The analyses were performed in **R** and **JAGS**.

The main R packages used across the scripts include:

* `MASS`
* `rjags`
* `coda`
* `haven`
* `readr`
* `dplyr`
* `tidyr`
* `purrr`
* `tibble`
* `stringr`
* `ggplot2`
* `survival`
* `survminer`

A separate installation of **JAGS** is required for the Bayesian analyses because `rjags` provides the interface between R and JAGS.

Some scripts check for missing R packages and install them automatically. JAGS itself must be installed separately.

---

## 3. ELSA empirical application

The empirical analysis jointly investigates longitudinal cognitive performance and transitions between frailty states.

### Longitudinal outcome

The longitudinal outcome is delayed word recall:

```text
cflisd
```

The longitudinal model contains fixed effects for time and baseline covariates, together with a subject-specific random intercept and random slope.

### Multi-state process

The four states are:

```text
0 = Non-frail
1 = Pre-frail
2 = Frail
3 = Death
```

Death is treated as an absorbing state.

The final empirical model retains nine possible transitions:

```text
0 -> 1
0 -> 2
0 -> 3

1 -> 0
1 -> 2
1 -> 3

2 -> 0
2 -> 1
2 -> 3
```

The final joint model includes the following eight baseline covariates:

* age
* cardiovascular condition
* high blood pressure
* diabetes
* stroke
* cancer
* arthritis
* eyesight

Repeated delayed recall measurements are modelled longitudinally rather than being included as a baseline covariate in the multi-state component.

The two processes are linked through transition-specific shared random effects representing the subject-specific cognitive level and rate of change.

---

## 4. Data availability and preprocessing

The original ELSA dataset is **not redistributed in this repository**.

Researchers wishing to reproduce the empirical analysis should first obtain authorised access to the required ELSA data.

The data preprocessing workflow is implemented in:

```text
Code/data_preprocessing.R
```

The script currently expects a local Stata dataset named:

```text
new_dataset_noL 1(4).dta
```

If the locally available dataset has a different name or location, the `DATA_FILE` setting at the beginning of the script should be changed accordingly.

The preprocessing script performs the main data preparation steps required for the longitudinal and multi-state analyses, including:

1. importing and cleaning the ELSA data;
2. checking and cleaning age, year, delayed recall, and frailty variables;
3. constructing longitudinal delayed-recall observations;
4. constructing individual frailty histories;
5. incorporating mortality information;
6. removing inappropriate post-death observations;
7. identifying and auditing ambiguous same-year death and frailty observations;
8. constructing multi-state transition intervals;
9. preparing baseline covariates; and
10. generating the datasets used for joint-model estimation.

### Treatment of mortality information

The mortality variable used in preprocessing is:

```text
raxyear
```

This records the **year of death rather than an exact date of death**.

Therefore, the preprocessing procedure does not treat the recorded year as an exact event time. In particular, frailty and mortality observations occurring within the same year are treated as potentially ambiguous rather than imposing an artificial ordering.

---

## 5. Processed data required by the final joint model

The final empirical fitting script expects the following two processed datasets:

```text
17_ELSA_JOINT_LONGITUDINAL.csv
19_ELSA_JOINT_MULTISTATE_TRANSITIONS.csv
```

These files are produced as part of the ELSA data-processing workflow and contain participant-level information. They are not included in this public repository.

After obtaining authorised access to the underlying ELSA data, the supplied preprocessing code can be used to reconstruct the analysis datasets.

---

## 6. Final ELSA joint model

The main empirical analysis is implemented in:

```text
Code/ELSA_joint_model application.R
```

This script fits the four-state Bayesian joint longitudinal and multi-state model.

The longitudinal component is a linear mixed-effects model with:

* fixed effects for time and baseline covariates;
* a subject-specific random intercept; and
* a subject-specific random slope for time.

The multi-state component contains transition-specific baseline intensities and transition-specific covariate effects.

The longitudinal and multi-state components are linked through shared random effects. For each transition, separate association parameters link the subject-specific longitudinal random intercept and random slope to the corresponding transition intensity.

### Final MCMC settings

The full empirical analysis uses:

```text
Number of chains:       2
Adaptation iterations:  1,000
Burn-in iterations:     2,000
Posterior iterations:   20,000 per chain
Thinning interval:      5
Chunk size:              250 iterations
```

The random seed is:

```text
20260903
```

The full model is computationally intensive. The script also contains a smaller `TEST` mode for checking whether the model compilation, burn-in, sampling, saving, and summarisation workflow operates correctly. Results from the `TEST` mode should not be used for substantive inference.

---

## 7. JAGS model specification

The repository also contains:

```text
Code/joint_shared_RE_model.jags
```

This provides the JAGS specification of the shared-random-effects joint model.

The model contains:

* longitudinal fixed effects;
* subject-specific random intercepts and slopes;
* the random-effects covariance structure;
* transition-specific baseline intensities;
* transition-specific covariate effects;
* shared random-effect association parameters; and
* Bayesian prior distributions.

---

## 8. Simulation studies

The repository contains three main simulation components.

### 8.1 Joint-model simulation

The principal simulation study is implemented in:

```text
Code/simulation_study.R
```

The study uses:

```text
Sample size:                 N = 500
Number of replications:      100
Maximum follow-up:           15
Planned longitudinal visits: 8
Random censoring rate:       0.02
```

The longitudinal data-generating model uses:

```text
beta0 = 2
beta1 = 0.5
beta2 = 1

sigma_b0 = 1
sigma_b1 = 0.6
rho = 0.3
sigma_y = 1
```

The multi-state component uses the baseline transition intensities:

```text
lambda12 = 0.08
lambda13 = 0.04
lambda23 = 0.10
```

with covariate effects:

```text
gamma12 = -0.5
gamma13 =  0.8
gamma23 =  0.5
```

and shared random-effect association parameters:

```text
alpha0_12 = 0.40
alpha0_13 = 0.80
alpha0_23 = 0.60

alpha1_12 = 0.15
alpha1_13 = 0.30
alpha1_23 = 0.25
```

The simulation evaluates parameter recovery using quantities including:

* empirical bias;
* root mean squared error (RMSE); and
* empirical coverage probability of 95% credible intervals.

The replication-level and summary outputs are provided in:

```text
Simulation/simulation_results.csv
Simulation/simulation_summary.csv
Simulation/joint_convergence_N500(1).csv
```

---

### 8.2 Longitudinal simulation

The standalone longitudinal simulation is implemented in:

```text
Code/longitudinal.R
```

This examines the Bayesian longitudinal random-intercept and random-slope model separately.

The simulation uses:

```text
Sample size:            500
Number of replications: 100
Repeated measurements:  8
```

The corresponding outputs are stored in:

```text
Simulation/longitudinal_simulation_results/
```

---

### 8.3 Conventional multi-state simulation

The standalone multi-state simulation is implemented in:

```text
Code/multistate_simulation_moderate_complete.R
```

This simulation considers a three-state transition structure:

```text
1 -> 2
1 -> 3
2 -> 3
```

The data-generating mechanism includes shared random-intercept and random-slope effects, while the fitted conventional multi-state model deliberately omits these shared random effects.

The simulation therefore provides a comparison with the joint modelling framework under a moderate shared-random-effects data-generating mechanism.

The corresponding outputs are stored in:

```text
Simulation/multistate_simulation_results_moderate/
```

---

## 9. Empirical results

The `Results/` directory contains the main numerical outputs used to verify the ELSA analysis.

### Data and transition checks

```text
00_observed_transition_counts.csv
```

contains observed counts for the retained multi-state transitions.

```text
01_subject_baseline_covariates_8cov.csv
```

contains the baseline covariate information used for the eight-covariate joint analysis.

```text
02_transition_retention_check.csv
```

records checks on the transition structure retained for modelling.

```text
03_covariate_centering_constants.csv
```

contains the centring constants applied to continuous covariates.

```text
04_longitudinal_model_data.csv
```

contains the prepared longitudinal model data used in the empirical analysis.

```text
05_longitudinal_subject_audit.csv
```

contains participant-level audit information for the longitudinal component.

```text
06_crude_transition_rates.csv
```

contains descriptive crude transition rates.

---

## 10. Posterior results

### MCMC diagnostics

```text
07_convergence_Rhat_ESS.csv
```

contains numerical convergence diagnostics, including Gelman-Rubin statistics and effective sample size measures.

### Longitudinal model

```text
10_longitudinal_core_results.csv
```

contains posterior summaries for the main longitudinal parameters, including the time effect and random-effects parameters.

```text
11_longitudinal_8cov_results.csv
```

contains posterior summaries for the eight baseline covariates in the longitudinal sub-model.

### Multi-state model

```text
12_multistate_8cov_HR_results.csv
```

contains transition-specific hazard ratios and posterior credible intervals for the baseline covariates.

### Joint association parameters

```text
13_joint_association_alpha_results.csv
```

contains the transition-specific association parameters linking the longitudinal random intercept and random slope to the multi-state transition intensities.

### Baseline transition intensities

```text
14_baseline_transition_intensities.csv
```

contains posterior summaries for the transition-specific baseline intensities.

---

## 11. Graphical outputs

The `Figures/` directory contains graphical outputs used to assess model convergence and present the multi-state results.

```text
Figures/trace_density.pdf
```

contains MCMC trace plots and posterior density plots for monitored parameters.

```text
Figures/gelman_plots.pdf
```

contains graphical Gelman-Rubin convergence diagnostics.

```text
Figures/multistate_HR_forest_plot.png
```

presents transition-specific hazard-ratio estimates and uncertainty intervals for the baseline covariates in the multi-state model.

---

## 12. Suggested reproduction workflow

A researcher wishing to reproduce the analyses should use the following general workflow.

### Step 1 — Obtain the ELSA data

Obtain authorised access to the ELSA data required for the empirical application.

### Step 2 — Prepare the ELSA analysis datasets

Place the source dataset in the working directory and run:

```r
source("Code/data_preprocessing.R")
```

If necessary, first modify the `DATA_FILE` path at the beginning of the script.

The preprocessing workflow constructs the datasets required for the empirical joint model.

### Step 3 — Fit the final ELSA joint model

Ensure that:

```text
17_ELSA_JOINT_LONGITUDINAL.csv
19_ELSA_JOINT_MULTISTATE_TRANSITIONS.csv
```

are available in the working directory expected by the model script.

Then run:

```r
source("Code/ELSA_joint_model application.R")
```

The full MCMC analysis can require substantial computation time.

### Step 4 — Check empirical results

Compare the resulting posterior summaries and diagnostics with the supplied files in:

```text
Results/
```

and:

```text
Figures/
```

### Step 5 — Run the joint simulation study

Run:

```r
source("Code/simulation_study.R")
```

Compare the resulting estimates with:

```text
Simulation/simulation_results.csv
Simulation/simulation_summary.csv
```

### Step 6 — Run the component simulations

For the standalone longitudinal simulation, run:

```r
source("Code/longitudinal.R")
```

For the conventional multi-state simulation, run:

```r
source("Code/multistate_simulation_moderate_complete.R")
```

Their outputs can be compared with the corresponding subdirectories within `Simulation/`.

---

## 13. Interpretation of the supplied materials

The supplied outputs are intended to make it possible to verify whether the analyses have been reproduced successfully.

For the Bayesian analyses, exact equality of individual MCMC draws is not necessary for substantive replication. Posterior summaries, credible intervals, convergence measures, and substantive conclusions should nevertheless be comparable when the same data, model specification, priors, random seeds, and MCMC settings are used.

Some parameters, particularly those corresponding to less frequent transitions, may have lower effective sample sizes or greater posterior uncertainty. The convergence diagnostics supplied in the repository should therefore be considered alongside the posterior estimates.

---

## 14. Reproducibility

Random seeds are specified within the simulation and empirical analysis scripts.

Small numerical differences may arise across computing environments because of differences in:

* R versions;
* JAGS versions;
* R package versions;
* operating systems;
* random-number generation; and
* numerical libraries.

The code, simulation outputs, posterior summaries, convergence diagnostics, and figures included here are intended to provide sufficient information for an authorised researcher to reconstruct and verify the analyses reported in the project.

---

## 15. Summary of materials

The repository provides:

* ELSA data preprocessing code;
* the final four-state Bayesian joint-model application;
* the JAGS shared-random-effects model specification;
* longitudinal, multi-state, and joint-model simulation code;
* simulation replication and summary results;
* processed-data audit outputs;
* posterior summaries from the empirical analysis;
* transition-specific hazard ratios;
* shared random-effect association estimates;
* MCMC convergence diagnostics; and
* major graphical outputs.

The original ELSA data and participant-level model input datasets are not publicly redistributed because they are subject to the relevant data-access conditions.
