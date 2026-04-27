###############################################################################
# Do economic and study planning calculations from cross-brain effect distributions
# Example: get_cost(0.16, 0.1, 2) 
# Requirements: run crossbrain_effect_estimator.R, pwr
###############################################################################

get_cost <- function(tau_sample, tau_pop, n_groups) {
  
  library(pwr)
  
  # params
  sig_level <- 0.05
  power_target <- 0.8
  cost_per_scan <- 500
  n_brain_areas <- 30000
  
  # expected n
  
  # tau_sample <- 0.1
  # tau_pop <- 0.04
  
  n_expected <- pwr.t.test(d=tau_sample,sig.level=sig_level,power=power_target,type="two.sample")$n
  n_actual <- pwr.t.test(d=tau_pop,sig.level=sig_level,power=power_target,type="two.sample")$n
  
  # this is N PER GROUP
  
  n_expected <- n_expected * n_groups
  n_actual <- n_actual * n_groups
  factor_increase = n_actual / n_expected
  
  # actual power
  
  pwr__from_sample <- pwr.t.test(d=tau_pop,sig.level=sig_level,n=n_expected,type="two.sample")
  
  # proportion detectable
  
  prop_detect__from_sample <- proportion_detectable(sig_level, 1-power_target , 0, tau_pop * sqrt(n_expected/2^2),2,2)
  prop_detect__from_pop <- proportion_detectable(sig_level, 1-power_target , 0, tau_pop * sqrt(n_actual/2^2),2,2)
  
  prop_detect_fdr__from_sample <- BH_proportion_detectable(0, sig_level, 1-power_target , 0, tau_pop * sqrt(n_expected/2^2),2,2)
  prop_detect_fdr__from_pop <- BH_proportion_detectable(0, sig_level, 1-power_target , 0, tau_pop * sqrt(n_actual/2^2),2,2)
  
  n_brain_areas_detectable__from_sample <- prop_detect__from_sample * n_brain_areas
  n_brain_areas_detectable__from_pop <- prop_detect__from_pop * n_brain_areas
  
  n_brain_areas_detectable_fdr__from_sample <- prop_detect_fdr__from_sample * n_brain_areas
  n_brain_areas_detectable_fdr__from_pop <- prop_detect_fdr__from_pop * n_brain_areas
  
  # cost
  
  cost__from_sample <- n_expected * cost_per_scan
  cost__from_pop <- n_actual * cost_per_scan
  cost_per_brain_area__from_sample <- cost__from_sample / n_brain_areas_detectable_fdr__from_sample
  cost_per_brain_area__from_pop <- cost__from_pop / n_brain_areas_detectable_fdr__from_pop
  factor_increase_cost_per_brain_area__for_sample <- cost_per_brain_area__from_sample / cost_per_brain_area__from_pop
  
  # save all to csv
  df <- data.frame(
    tau_sample = tau_sample,
    tau_pop = tau_pop,
    n_expected = n_expected,
    n_actual = n_actual,
    factor_increase = factor_increase,
    prop_detect_fdr__from_sample = prop_detect_fdr__from_sample,
    prop_detect_fdr__from_pop = prop_detect_fdr__from_pop,
    prop_detect__from_sample = prop_detect__from_sample,
    prop_detect__from_pop = prop_detect__from_pop,
    n_brain_areas_detectable__from_sample = n_brain_areas_detectable__from_sample,
    n_brain_areas_detectable__from_pop = n_brain_areas_detectable__from_pop,
    n_brain_areas_detectable_fdr__from_sample = n_brain_areas_detectable_fdr__from_sample,
    n_brain_areas_detectable_fdr__from_pop = n_brain_areas_detectable_fdr__from_pop,
    cost__from_sample = cost__from_sample,
    cost__from_pop = cost__from_pop,
    cost_per_brain_area__from_sample = cost_per_brain_area__from_sample,
    cost_per_brain_area__from_pop = cost_per_brain_area__from_pop,
    factor_increase_cost_per_brain_area__for_sample
  )
  
  # print the following summary:
  # Uncorrected: expect n=*n* subjects to detect d=*d* (power = *pwr*)
  # Corrected: n=*n* subjects to detect d=*d* (factor increase = *factor_increase*)
  # Uncorrected: Detect *prop_detect_fdr__from_sample* of the brain (*n_brain_areas_detectable_fdr__from_sample*; FDR)
  # Corrected: Detect *prop_detect_fdr__from_pop* of the brain (*n_brain_areas_detectable_fdr__from_pop*; FDR)
  # Uncorrected: Costs $*cost__from_sample* overall -> $*cost_per_brain_area__from_sample*/brain area
  # Corrected: Costs $*cost__from_pop* overall -> $*cost_per_brain_area__from_pop*/brain area
  
  write.csv(df, "~/Desktop/economics.csv", row.names = FALSE)
  summary_df <- data.frame(
    summary = c(
      sprintf("Uncorrected: expect n=%.0f subjects to detect d=%.2f (power = %.2f)", n_expected, tau_sample, pwr__from_sample$power),
      sprintf("Corrected: n=%.0f subjects to detect d=%.2f (factor increase = %.2f)", n_actual, tau_pop, factor_increase),
      sprintf("Uncorrected: Detect %.2f%% of the brain (%.0f; FDR)", prop_detect_fdr__from_sample * 100, n_brain_areas_detectable_fdr__from_sample),
      sprintf("Corrected: Detect %.2f%% of the brain (%.0f; FDR)", prop_detect_fdr__from_pop * 100, n_brain_areas_detectable_fdr__from_pop),
      sprintf("Uncorrected: Costs $%.2f overall -> $%.2f/brain area (factor increase = %.2f)", cost__from_sample, cost_per_brain_area__from_sample, factor_increase_cost_per_brain_area__for_sample),
      sprintf("Corrected: Costs $%.2f overall -> $%.2f/brain area", cost__from_pop, cost_per_brain_area__from_pop)
    )
  )
  
  write.csv(summary_df, "~/Desktop/economics_summary.csv", row.names = FALSE)

  
 }