# ERP Additional Materials

Supporting code, results, figures, and simulation outputs for Joint-Model-JAGS.

## Contents

- [Code/](Code/): data preprocessing, ELSA joint model, joint simulation study, and supplementary model scripts.
- [Results/](Results/): existing ELSA analysis CSV outputs.
- [Figures/](Figures/): trace/density and Gelman diagnostic PDFs, and the existing forest plot PNG.
- [Simulation/](Simulation/): joint simulation results and summary, convergence diagnostics, and supplementary longitudinal/multistate simulation outputs.

## Missing requested files

The following files were not present in the repository at reorganization time and have not been fabricated:

- `Technical_Appendix.pdf`
- `Code/generate_results.R` (result generation is included in the existing ELSA analysis script; no separate script was found)
- `Results/01_subject_baseline_covariates_8cov.csv`
- `Figures/multistate_HR_forest_plot.pdf` (only the PNG version is available)

This README was newly created because the repository did not contain a README.

## Execution notes

Existing script contents and all scientific outputs are preserved byte-for-byte. This is a directory/file-name reorganization, not a change to the analysis or a rerun.

Scripts retain their original input/output settings and generated output filenames. Review these settings and choose a suitable R working directory before execution; rerunning scripts does not automatically reproduce this publication directory layout.

The preprocessing script expects `new_dataset_noL 1(4).dta`. The ELSA analysis script expects `17_ELSA_JOINT_LONGITUDINAL.csv` and `19_ELSA_JOINT_MULTISTATE_TRANSITIONS.csv`. These input files are not included in the repository.

`Simulation/simulation_results.csv` and `Simulation/simulation_summary.csv` are the existing joint-model N=500 raw results and performance summary, respectively. Supplementary longitudinal and multistate simulation outputs are retained in their own subdirectories.

## Original-to-current file mapping

| Original path | Current path (relative to this folder) |
| --- | --- |
| `application/00_observed_transition_counts.csv` | `Results/00_observed_transition_counts.csv` |
| `application/02_transition_retention_check.csv` | `Results/02_transition_retention_check.csv` |
| `application/03_covariate_centering_constants.csv` | `Results/03_covariate_centering_constants.csv` |
| `application/06_crude_transition_rates.csv` | `Results/06_crude_transition_rates.csv` |
| `application/07_convergence_Rhat_ESS.csv` | `Results/07_convergence_Rhat_ESS.csv` |
| `application/08_trace_density.pdf` | `Figures/trace_density.pdf` |
| `application/09_gelman_plots.pdf` | `Figures/gelman_plots.pdf` |
| `application/10_longitudinal_core_results.csv` | `Results/10_longitudinal_core_results.csv` |
| `application/11_longitudinal_8cov_results.csv` | `Results/11_longitudinal_8cov_results.csv` |
| `application/12_multistate_8cov_HR_results.csv` | `Results/12_multistate_8cov_HR_results.csv` |
| `application/13_joint_association_alpha_results.csv` | `Results/13_joint_association_alpha_results.csv` |
| `application/14_baseline_transition_intensities.csv` | `Results/14_baseline_transition_intensities.csv` |
| `application/16_multistate_8cov_HR_forest_plot.png` | `Figures/multistate_HR_forest_plot.png` |
| `application/ELSA_joint_model_chunked_checkpoint.R` | `Code/ELSA_joint_model.R` |
| `data_generate_corrected.R` | `Code/data_preprocessing.R` |
| `final simulation/joint(4).R` | `Code/simulation_study.R` |
| `final simulation/joint_convergence_N500(1).csv` | `Simulation/joint_convergence_N500(1).csv` |
| `final simulation/joint_performance_N500(1).csv` | `Simulation/simulation_summary.csv` |
| `final simulation/joint_raw_results_N500(1).csv` | `Simulation/simulation_results.csv` |
| `final simulation/joint_shared_RE_model.jags` | `Code/joint_shared_RE_model.jags` |
| `longitudinal.R` | `Code/longitudinal.R` |
| `longitudinal_simulation_results/full_results/longitudinal_raw_results_N500.csv` | `Simulation/longitudinal_simulation_results/full_results/longitudinal_raw_results_N500.csv` |
| `longitudinal_simulation_results/longitudinal_complete_summary.csv` | `Simulation/longitudinal_simulation_results/longitudinal_complete_summary.csv` |
| `longitudinal_simulation_results/longitudinal_complete_summary_rounded.csv` | `Simulation/longitudinal_simulation_results/longitudinal_complete_summary_rounded.csv` |
| `longitudinal_simulation_results/longitudinal_five_core_metrics.csv` | `Simulation/longitudinal_simulation_results/longitudinal_five_core_metrics.csv` |
| `longitudinal_simulation_results/longitudinal_five_core_metrics_rounded.csv` | `Simulation/longitudinal_simulation_results/longitudinal_five_core_metrics_rounded.csv` |
| `longitudinal_simulation_results/longitudinal_individual_results.csv` | `Simulation/longitudinal_simulation_results/longitudinal_individual_results.csv` |
| `multistate_simulation_moderate_complete.R` | `Code/multistate_simulation_moderate_complete.R` |
| `multistate_simulation_results_moderate/conventional_multistate_model.jags` | `Simulation/multistate_simulation_results_moderate/conventional_multistate_model.jags` |
| `multistate_simulation_results_moderate/event_censoring_summary.csv` | `Simulation/multistate_simulation_results_moderate/event_censoring_summary.csv` |
| `multistate_simulation_results_moderate/performance_summary.csv` | `Simulation/multistate_simulation_results_moderate/performance_summary.csv` |
| `multistate_simulation_results_moderate/raw_posterior_results.csv` | `Simulation/multistate_simulation_results_moderate/raw_posterior_results.csv` |
| `multistate_simulation_results_moderate/raw_posterior_results.rds` | `Simulation/multistate_simulation_results_moderate/raw_posterior_results.rds` |
| `multistate_simulation_results_moderate/raw_results_progress.rds` | `Simulation/multistate_simulation_results_moderate/raw_results_progress.rds` |
| `multistate_simulation_results_moderate/sessionInfo.txt` | `Simulation/multistate_simulation_results_moderate/sessionInfo.txt` |
