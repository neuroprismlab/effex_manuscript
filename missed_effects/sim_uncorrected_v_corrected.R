####################################################
# 
# Simulations for spatial extent of detections
#
####################################################


######### SET PATHS & sim_params #########

# paths
project_dir <- "/Users/stephanienoble/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/"
out_master_dir <- paste0(project_dir, "manuscript/figures/plots/sim/")
crossbrain_effect_estimator_script <- paste0(project_dir, "scripts/crossbrain_effects/crossbrain_effect_estimator.R")
effect_size_estimate_data <- paste0(project_dir, "manuscript/figures/plots/crossbrain_effects/pooling.none.motion.regression.mv.none/point_res.Rdata")

# libraries
library(MASS)
library(ggplot2)
library(reshape2)
library(pwr)
source(crossbrain_effect_estimator_script)

# parameters
#   -for sim
sim_params <- list (
  outcome_category = "psychological", # c("psychological","physical", "task activation", "task connectivity")
  n_subjects__gt = 10000, # ground truth full sample size
  n_regions = 500,
  sample_sizes = c(25, 50, 100, 500, 1000, 5000),
  n_reps = 100,
  alpha_fdr = 0.05,
  targeted_power = 0.8,
  target_effect_type__from_basis = "mean_without_selection" # c("min", "median", "max", "mean_without_selection")
)

#   -cross-brain area and cross-subject distribution parameters (Cohen's)
mu_crossbrain <- 0
# var_crossbrain <- 0.04^2
# var_crosssubject <- 1.29
#   -covariance between brain areas
# rho_values <- seq(0, 0.9, by = 0.1)
# rho <- c(0.5, 0.3)

#   -for parametric methods
n_groups <- 1 # note: this is hard-coded in the simulations
n_sides <- 1 # note: this is hard-coded in the simulations
proportion_tp_target <- 0.05
targeted_d <- 0.2



######### SETUP #########

# load sim_params for this outcome
load(effect_size_estimate_data) # load "res" from cross-brain effect estimation
var_crossbrain <- res[sim_params$outcome_category,"est"]
var_crosssubject <- res[sim_params$outcome_category,"phi2_est"]
sd_crossbrain <- sqrt(var_crossbrain)
sd_crosssubject <- sqrt(var_crosssubject)

# set files for this outcome
out_dir <- paste0(out_master_dir, sim_params$outcome_category, "/",sim_params$target_effect_type__from_basis, "_effect/")
sim_results_file <- paste0(out_dir, "sim_results.Rdata")
sim_params_file <- paste0(out_dir, "sim_params.rds")
summary_file <- paste0(out_dir, "summary_results.rds")

# make out dir for this outcome
if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }

# calculate target based on number of regions
num_tp_target <- 0.05 * sim_params$n_regions



######### FUNCTIONS #########

# generate per-region and per-subject effect sizes
simulate_data <- function(n_subjects, n_regions, mu_crossbrain, sd_crossbrain, sd_crosssubject) {
  ground_truth_true_effects <- rnorm(n_regions, mean = mu_crossbrain, sd = sd_crossbrain)
  sim_data <- matrix(0, nrow = n_subjects, ncol = n_regions)
  for (i in 1:n_subjects) {
    subject_var <- rnorm(n_regions, mean = 0, sd = sd_crosssubject)
    sim_data[i, ] <- as.vector(ground_truth_true_effects * sd_crosssubject + subject_var)
  }
  return(list(
    ground_truth_true_effects = ground_truth_true_effects,
    sim_data = sim_data
  ))
}

# return detections and effect sizes
get_pos_and_effects <- function(X, alpha = 0.05, use_correction = TRUE) {
  
  results <- apply(X, 2, function(region) t.test(region, alternative = "greater"))
  
  pvals <- sapply(results, function(region) region$p.value)
  if (use_correction) {
    pvals_final <- p.adjust(pvals, method = "fdr")
  } else {
    pvals_final <- pvals
  }
  positives <- pvals_final < alpha
  
  est <- sapply(results, function(region) region$estimate)
  sd <- apply(X, 2, sd) # TODO: uncorrect for n-1? depends on whether people usually apply the correction (not a big deal for large N)
  effects <- est / sd
  
  return(list(positives = positives, effects = effects))
}
    


########################### MAIN ###########################



######### SIMULATE DATA #########

# if sim_data exists, see whether current params are the same as previous
sim_data_exists <- FALSE
if (file.exists(sim_params_file)) {
# if (exists("sim_data") && nrow(sim_data) == sim_params$n_subjects__gt && ncol(sim_data) == sim_params$n_regions && all(names(positives__replication) == sim_params$sample_sizes) && sim_params$outcome_category == outcome_category__previous ) {
  old_sim_params <- readRDS(sim_params_file)
  if (identical(names(old_sim_params), names(sim_params)) &&
      all(vapply(names(sim_params), function(param) {
        identical(old_sim_params[[param]], sim_params[[param]])
      }, logical(1)))) {
    sim_data_exists <- TRUE
  }
}

# if it exists, ask whether to replace
run_replication_only <- FALSE
if (sim_data_exists) {
  run_sim <- readline(prompt = "Simulations using same params already exist. Replace? (yes/no; select yes if sim code has changed): ")
  run_sim <- tolower(run_sim) == "yes"
  if (run_sim) {
    run_replication_only <- readline(prompt = "Replace all? Otherwise will just re-run replication (yes/no): ")
    run_replication_only <- tolower(run_replication_only) == "no"
  }
} else {
  run_sim <- TRUE
}


# Start Simulation

if (run_sim) {
if (!run_replication_only) {

# simulate master dataset
sim_info <- simulate_data(sim_params$n_subjects__gt, sim_params$n_regions, mu_crossbrain, sd_crossbrain, sd_crosssubject)
ground_truth_true_effects <- sim_info$ground_truth_true_effects  #sim_info[1:sim_params$n_regions]
sim_data <- sim_info$sim_data #matrix(sim_info[-(1:sim_params$n_regions)], nrow = sim_params$n_subjects__gt, ncol = sim_params$n_regions)

if (sim_params$n_regions %% 2 == 0) { # drop the min effect before mean to make things easier
  target_effect__actual_mean <- mean(ground_truth_true_effects[ground_truth_true_effects > min(ground_truth_true_effects > 0)], na.rm = TRUE)
} else {
  target_effect__actual_mean <- mean(ground_truth_true_effects[ground_truth_true_effects > 0])
}
# get number that is closest to but bigger than mean so it's a liberal comparison (easier for inflation-based planning to hit)
target_effect_idx__actual_mean <- which(ground_truth_true_effects==sort(ground_truth_true_effects[which(ground_truth_true_effects >= target_effect__actual_mean)], decreasing=FALSE)[1])


######### 1A. BASIS STUDY #########

# initialize
sig_mask__basis <- list()
num_tp__basis <- list()
n_pos_above_target_effect_type__basis <- c()
num_tp_expect__from_basis <- list()
num_fp_expect__from_basis <- list()
num_tp_expect__in_basis <- list()

target_effect_idx__basis <- list()
expected_n_to_replicate_basis_effect <- c()

for (sample_size in sim_params$sample_sizes) {
  
  # preallocate for this sample size
  sig_mask__basis[[as.character(sample_size)]] <- matrix(FALSE, nrow = sim_params$n_reps, ncol = sim_params$n_regions)
  num_tp__basis[[as.character(sample_size)]] <- matrix(0, nrow = sim_params$n_reps, ncol = 1)
  n_pos_above_target_effect_type__basis[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
  num_tp_expect__from_basis[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
  num_fp_expect__from_basis[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
  num_tp_expect__in_basis[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
  
  target_effect_idx__basis[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
  expected_n_to_replicate_basis_effect[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
  
  n_reps_successful <- 0
  for (rep in 1:sim_params$n_reps) {
    
    # Basis Study: estimate the mean effect size (Cohen's d) for each region and use the significant results for power analysis
    
    subsample_indices__study1 <- sample(1:sim_params$n_subjects__gt, sample_size)
    data_subsample1 <- sim_data[subsample_indices__study1, ]
    
    this_res <- get_pos_and_effects(data_subsample1)
    this_sig_mask <- this_res$positives
    sig_mask__basis[[as.character(sample_size)]][rep, ] <- this_sig_mask
    
    
    # Store parameters for replication study: general or significant effect sizes and associated n
    # Goal: detect the same effects
    
    # setup
    if (sum(this_sig_mask) == 0) { # skip if no sig results from study 1
      next
    }
    n_reps_successful <- n_reps_successful + 1
    num_tp__basis[[as.character(sample_size)]][rep] <- sum(this_sig_mask[ground_truth_true_effects > 0])
    
    sig_effects_basis <- as.numeric(this_res$effects[this_sig_mask])
    
    ### "Same Sample" Case: let's say we just want to be powered to detect the same effects we got in the basis study, without any adjustment
    if (sim_params$target_effect_type__from_basis == "same_sample") {
      target_effect_size__from_basis <- min(sig_effects_basis, na.rm = TRUE) # there's an expectation that we'll be able to detect the same minimum effect size previously detected
      n_pos_above_target_effect_type__basis[[as.character(sample_size)]][rep] <- sum(this_sig_mask)
      
      num_tp_expect__from_basis[[as.character(sample_size)]][rep] <- sum(this_sig_mask) * sim_params$targeted_power # we're planning to detect all the original effects, even if in our heart of hearts we know some are false positives
      num_fp_expect__from_basis[[as.character(sample_size)]][rep] <- (sim_params$alpha_fdr/(1-sim_params$alpha_fdr)) * num_tp_expect__from_basis[[as.character(sample_size)]][rep] # 5 fp for every 95 tp
      
      num_tp_expect__in_basis[[as.character(sample_size)]][rep] <- sum(this_sig_mask) * (1-sim_params$alpha_fdr) # expect 5% of sig are FP
    
    ### "Min" Case: alright, let's say we're doing our best and planning for the lowest effect size reported
    } else if (sim_params$target_effect_type__from_basis == "min") {
      target_effect_size__from_basis <- min(sig_effects_basis, na.rm = TRUE)
      n_pos_above_target_effect_type__basis[[as.character(sample_size)]][rep] <- sum(this_sig_mask)
     
      num_tp_expect__from_basis[[as.character(sample_size)]][rep] <- sum(this_sig_mask) * sim_params$targeted_power # conservative estimate when planning for target power to detect this minimum effect size
      num_fp_expect__from_basis[[as.character(sample_size)]][rep] <- (sim_params$alpha_fdr/(1-sim_params$alpha_fdr)) * num_tp_expect__from_basis[[as.character(sample_size)]][rep] # 5 fp for every 95 tp
      
      num_tp_expect__in_basis[[as.character(sample_size)]][rep] <- sum(this_sig_mask) * (1-sim_params$alpha_fdr) # expect 5% of sig are FP
    
    ### "Median" Case: we're more optimistic and planning for the average effect size reported
    } else if (sim_params$target_effect_type__from_basis == "median") {
      if (sum(this_sig_mask) == 1) { # if only 1, choose that one
        target_effect_size__from_basis <- sig_effects_basis
      } else if (sum(this_sig_mask) %% 2 == 0) {
        # if even, remove the highest and take median (conservative approach)
        sig_effects_basis_tmp <- sig_effects_basis[sig_effects_basis < max(sig_effects_basis, na.rm = TRUE)]
        target_effect_size__from_basis <- median(sig_effects_basis_tmp, na.rm = TRUE)
      } else {
        target_effect_size__from_basis <- median(sig_effects_basis, na.rm = TRUE)
      }
      n_pos_above_target_effect_type__basis[[as.character(sample_size)]][rep] <- sum(sig_effects_basis >= target_effect_size__from_basis)
      
      num_tp_expect__from_basis[[as.character(sample_size)]][rep] <- sum(sig_effects_basis >= target_effect_size__from_basis) * sim_params$targeted_power # conservative estimate when planning for target power to detect this median effect size
      num_fp_expect__from_basis[[as.character(sample_size)]][rep] <- (sim_params$alpha_fdr/(1-sim_params$alpha_fdr)) * num_tp_expect__from_basis[[as.character(sample_size)]][rep] # 5 fp for every 95 tp
      
      num_tp_expect__in_basis[[as.character(sample_size)]][rep] <- sum(sig_effects_basis >= target_effect_size__from_basis) * (1-sim_params$alpha_fdr) # expect 5% of sig are FP (probably conservative with these larger-than-median effects)
        
    ### "Max" Case: we're really optimistic and planning for the highest effect size reported (e.g., peak activation)
    } else if (sim_params$target_effect_type__from_basis == "max") {
      target_effect_size__from_basis <- max(sig_effects_basis, na.rm = TRUE)
      n_pos_above_target_effect_type__basis[[as.character(sample_size)]][rep] <- 1
      
      num_tp_expect__from_basis[[as.character(sample_size)]][rep] <- sim_params$targeted_power # conservative estimate when planning to for target power to detect this maximum effect size
      num_fp_expect__from_basis[[as.character(sample_size)]][rep] <- (sim_params$alpha_fdr/(1-sim_params$alpha_fdr)) * sim_params$targeted_power # 5 fp for every 95 tp
      
      num_tp_expect__in_basis[[as.character(sample_size)]][rep] <- 1-sim_params$alpha_fdr # expect 5% of sig are FP
    
    ### "Mean without selection" Case: we're doing our best to plan for the average effect size in the brain, but we don't want to be biased by selection of significant effects
    } else if (sim_params$target_effect_type__from_basis == "mean_without_selection") {
      this_res <- get_pos_and_effects(data_subsample1, alpha=1, use_correction = FALSE) # only for getting effect sizes for all regions, not sig mask
      target_effect_size__from_basis <- mean(as.numeric(this_res$effects[this_res$effects > 0]), na.rm = TRUE)
      n_pos_above_target_effect_type__basis[[as.character(sample_size)]][rep] <- sum(this_res$effects >= target_effect_size__from_basis)
      
      num_tp_expect__from_basis[[as.character(sample_size)]][rep] <- sum(this_res$effects >= target_effect_size__from_basis) * sim_params$targeted_power # conservative estimate when planning for target power to detect this mean effect size
      num_fp_expect__from_basis[[as.character(sample_size)]][rep] <- (sim_params$alpha_fdr/(1-sim_params$alpha_fdr)) * num_tp_expect__from_basis[[as.character(sample_size)]][rep] # 5 fp for every 95 tp
      
      num_tp_expect__in_basis[[as.character(sample_size)]][rep] <- sum(this_res$effects >= target_effect_size__from_basis) * (1-sim_params$alpha_fdr) # expect 5% of sig are FP
    
    } else {
      stop("Invalid sim_params$target_effect_type__from_basis value.")
    }
    
    # store exactly where the closest univariate effect occurred
    if (sim_params$target_effect_type__from_basis == "mean_without_selection") {
      # using mean means there's no true effect for comparison, so we'll pick the closest mean idx of ground truth>0 since to stand in for the fact that there is always a true effect for ground truth>0
      target_effect_idx__basis[[as.character(sample_size)]][rep] <- target_effect_idx__actual_mean
    } else {
      target_effect_idx__basis[[as.character(sample_size)]][rep] <- which(as.numeric(this_res$effects) == target_effect_size__from_basis)
    }
    
    # store planned sample size for replication
    expected_n_to_replicate_basis_effect[[as.character(sample_size)]][rep] <- ceiling(pwr.t.test(power=sim_params$targeted_power, d = target_effect_size__from_basis, sig.level = sim_params$alpha_fdr, type = "one.sample", alternative = "greater")$n)
    
  }
}



######### 1B. REPLICATION STUDY #########

} else {
  message("Using existing sim_data and basis study results for replication study.")
  load(sim_results_file) # note this is only to get pre-calculated basis study and we will replace the replication
}

# initialize
positives__replication <- list()
effect_sizes__replication <- list()
positives_uncorr__replication <- list()
effect_sizes_uncorr__replication <- list()
proportion_overlap_with_basis__replication <- list()

for (sample_size in sim_params$sample_sizes) {
  
  # preallocate for this sample size
  positives__replication[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = sim_params$n_regions)
  effect_sizes__replication[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = sim_params$n_regions)
  positives_uncorr__replication[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = sim_params$n_regions)
  effect_sizes_uncorr__replication[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = sim_params$n_regions)
  proportion_overlap_with_basis__replication[[as.character(sample_size)]] <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
  
  for (rep in 1:sim_params$n_reps) {
    
    if (is.na(expected_n_to_replicate_basis_effect[[as.character(sample_size)]][rep])) {
      next
    }
    if (expected_n_to_replicate_basis_effect[[as.character(sample_size)]][rep] > sim_params$n_subjects__gt) {
      warning(paste0("Planned sample size for replication (", expected_n_to_replicate_basis_effect[[as.character(sample_size)]][rep], ") exceeds available subjects (", sim_params$n_subjects__gt, "). Skipping this replication."))
      next
    }
    
    # sample sim_data
    subsample_indices <- sample(1:sim_params$n_subjects__gt, expected_n_to_replicate_basis_effect[[as.character(sample_size)]][rep]) # TODO: some overlap
    data_subsample <- sim_data[subsample_indices,]
    
    # results
    this_res <- get_pos_and_effects(data_subsample)
    positives__replication[[as.character(sample_size)]][rep, ] <- this_res$positives
    effect_sizes__replication[[as.character(sample_size)]][rep, ] <- this_res$effects
    
    this_sig_mask__basis <- sig_mask__basis[[as.character(sample_size)]][rep, ]
    proportion_overlap_with_basis__replication[[as.character(sample_size)]][rep] <- sum(which(this_res$positives) %in% which(this_sig_mask__basis)) / sum(this_sig_mask__basis)
    
    # repeat for uncorrected - used for single-region analysis
    this_res_uncorr <- get_pos_and_effects(data_subsample, use_correction = FALSE)
    positives_uncorr__replication[[as.character(sample_size)]][rep, ] <- this_res_uncorr$positives
    effect_sizes_uncorr__replication[[as.character(sample_size)]][rep, ] <- this_res_uncorr$effects
    
  }
}


  
######### 2. BENCHMARK-BASED STUDY #########

# Study 2B: Corrected approach: use ground truth effect sizes for study planning

# Goal: be adequately powered to detect a target proportion of effects 
# Goal: be adequately powered to detect effects above a target size
# Subgoal: understand how many effects were missed in study 1

expected_tpr__corrected <- c()
expected_num_tp__corrected <- c()

# reflecting on study 1, estimate number of detections & avg power using our knowledge of the ground truth
for (sample_size in sim_params$sample_sizes) {
  # expected_prop_detectable__corrected[[as.character(sample_size)]] <- BH_proportion_detectable(0, sim_params$alpha_fdr, 1-sim_params$targeted_power, 0, sd_crossbrain*sqrt(sample_size/n_groups^2), n_groups, n_sides)
  expected_tpr__corrected[[as.character(sample_size)]] <- BHpower(pi0=0, alphaFDR=sim_params$alpha_fdr, 0, sd_crossbrain*sqrt(sample_size/n_groups^2), n_groups, n_sides)
  expected_num_tp__corrected[[as.character(sample_size)]] <- sim_params$n_regions * expected_tpr__corrected[[as.character(sample_size)]]
}

# Let's instead plan for a desired number of detections
for (n in seq(10, 10000, by = 10)) {
  this_num_tp <- sim_params$n_regions * BHpower(0, sim_params$alpha_fdr, 0, sd_crossbrain*sqrt(n/n_groups^2), n_groups, n_sides)
  if (this_num_tp >= num_tp_target) {
    expected_n_for_target_num_tp__corrected <- n
    break
  }
}

# Save results

save(
  list = c("num_tp__basis", "n_pos_above_target_effect_type__basis", "num_tp_expect__from_basis", "num_fp_expect__from_basis", "num_tp_expect__in_basis",
           "target_effect_idx__basis", "expected_n_to_replicate_basis_effect",
           "positives__replication", "effect_sizes__replication", "positives_uncorr__replication", "effect_sizes_uncorr__replication", "proportion_overlap_with_basis__replication",
           "expected_tpr__corrected", "expected_num_tp__corrected", "expected_n_for_target_num_tp__corrected"),
  file = sim_results_file
)
saveRDS(sim_params, file = sim_params_file)
  
} else { # If skip sim, load previously saved results
  load(sim_results_file)
}


# test with resampling 
# num_tp__actual_from_corrected <- matrix(NA, nrow = sim_params$n_reps, ncol = 1)
# for (rep in 1:sim_params$n_reps) {
#   subsample_indices <- sample(1:sim_params$n_subjects__gt, expected_n_for_target_num_tp__corrected)
#   data_subsample <- sim_data[subsample_indices,]
#   this_res <- get_pos_and_effects(data_subsample)
#   num_tp__actual_from_corrected[rep] <- sum(this_res$positives & ground_truth_true_effects > 0) / sim_params$n_regions
# }
# mean_num_tp__actual_from_corrected <- mean(num_tp__actual_from_corrected, na.rm = TRUE)

# We can also see what happens if we plan to be adequately powered to detect everything above a certain effect size
# how much of ground truth distribution > 0.2

# prop_gt_dist_above_targeted_d <- 1 - pnorm(targeted_d, mean = mu_crossbrain, sd = sd_crossbrain) # TODO: will have to update
# if (n_sides == 2) { # take advantage of symmetry
#   prop_gt_dist_above_targeted_d <- prop_gt_dist_above_targeted_d * 2
# }
# # and now we can again plan for a targeted number of detections
# num_tp_target <- prop_gt_dist_above_targeted_d
# for (n in seq(10, 10000, by = 10)) {
#   expected_num_tp__corrected <- sim_params$n_regions * BH_proportion_detectable(targeted_d, sim_params$alpha_fdr, 1-sim_params$targeted_power, 0, sd_crossbrain*sqrt(n/n_groups^2), n_groups, n_sides)
#   if (expected_num_tp__corrected >= num_tp_target) {
#     planned_n_for_d__corrected <- n
#     break
#   }
# }
# 
# # test with resampling 
# true_effect_magnitude <- matrix(NA, nrow = sim_params$n_reps, ncol = sim_params$n_regions)
# for (rep in 1:sim_params$n_reps) {
#   subsample_indices <- sample(1:sim_params$n_subjects__gt, planned_n_for_d__corrected)
#   data_subsample <- sim_data[subsample_indices,]
#   this_res <- get_pos_and_effects(data_subsample)
#   true_effect_magnitude[rep,this_res$positives] <- ground_truth_true_effects[this_res$positives]
# }
# mean_true_effect_magnitude <- mean(true_effect_magnitude, na.rm = TRUE)



### ESTIMATE ERROR RATES ###

# Estimate error rates based on ground truth crossbrain distribution.
#
# specs drives three things:
#   1. which columns appear in results_by_region / results_by_rep / results_for_basis_study_region
#   2. how the summary sim_data frame is initialised and filled (region/rep get _mean+_sd; target mean only)
#   3. which plots are generated and with what labels / y-limits
#
# Fields per entry:
#   source  "region"  → per-brain-region long table; summary gets <col>_mean + <col>_sd
#           "rep"     → per-rep long table;           summary gets <col>_mean + <col>_sd
#           "target"  → target-effect table;          summary gets <col> (mean only, no sd)
#           "overall"    → calculated over all regions and reps; summary gets "fwer"
#   col     column name written into the results table
#   label   plot axis / title label
#   ylim    y-axis limits for the plot
#
# To add a new metric: add one list() entry here. No other changes needed.
#
# Misc notes:
# Type M Error = E( |est| / actual )
# Type S Error = E( sign(est) != sign(actual) )
#

specs <- list(
  list(source = "overall",col = "num_tp__expected_corr",label = "Expected Number of True Positives (Corrected)",    ylim = c(0, sim_params$n_regions/2), plot = TRUE),
  list(source = "overall",col = "fwer",                 label = "Familywise Error Rate",                            ylim = c(0, 1), plot = FALSE),
  list(source = "rep",    col = "n_pos_above_target_effect_type__basis", label = "Number of Positives in Basis Study",  ylim = c(0, sim_params$n_regions/2), plot = FALSE),
  list(source = "rep",    col = "num_tp__replication",  label = "Number of True Positives in Replication Study",    ylim = c(0, sim_params$n_regions/2), plot = TRUE),
  list(source = "rep",    col = "num_tp__basis",        label = "Number of True Positives in Basis Study",          ylim = c(0, sim_params$n_regions/2), plot = TRUE),
  list(source = "rep",    col = "num_tp_expect__from_basis",  label = "Expected Number of TP From Basis Study (Adj if Min Sig Basis)",   ylim = c(0, sim_params$n_regions/2), plot = TRUE),
  list(source = "rep",    col = "expect_v_actual_n_tp__based_on_basis",  label = "Proportion of Expected TPs Detected (Uncorr)",          ylim = c(0, 1.5), plot = TRUE),
  list(source = "rep",    col = "expect_v_actual_n_tp__based_on_corr",   label = "Proportion of Expected TPs Detected (Corr)",            ylim = c(0, 1.5), plot = TRUE),
  list(source = "rep",    col = "prop_tp__replication", label = "Replication Proportion of True Positives",         ylim = c(0, 0.5), plot = TRUE),
  list(source = "rep",    col = "prop_above_target_effect_type__basis",  label = "Expected Proportion of True Positives (FDR Adj)",       ylim = c(0, 0.5), plot = FALSE),
  list(source = "rep",    col = "overlap",              label = "Proportion of Original Findings Replicated",  ylim = c(0, 1), plot = TRUE),
  list(source = "rep",    col = "expected_n__replication", label = "Planned Sample Size for Target Effect (Basis)",      ylim = c(0, max(sim_params$sample_sizes)), plot = TRUE),
  list(source = "rep",    col = "fdr",                  label = "False Discovery Rate",                             ylim = c(0, 0.15), plot = TRUE),
  list(source = "region", col = "tpr",                  label = "Power (Mean Over Regions)",                        ylim = c(0, 0.5), plot = TRUE),
  list(source = "region", col = "type_m",               label = "Type M Error (Mean Over Regions)",                 ylim = c(0, 10), plot = TRUE),
  list(source = "region", col = "type_s",               label = "Type S Error (Mean Over Regions)",                 ylim = c(0, 0.5), plot = TRUE),
  list(source = "target", col = "tpr_target_univ",      label = "Power for Target (Univariate)",                    ylim = c(0, 1), plot = TRUE),
  list(source = "target", col = "type_m_target_univ",   label = "Type M Error for Target (Univariate)",             ylim = c(0, 10), plot = TRUE),
  list(source = "target", col = "type_s_target_univ",   label = "Type S Error for Target (Univariate)",             ylim = c(0, 0.5), plot = TRUE),
  list(source = "target", col = "target_true",     label = "Rate Target Effect is a True Effect",              ylim = c(0, 1), plot = TRUE),
  list(source = "target", col = "target_detected", label = "Rate Target Effect is Detected",                   ylim = c(0, 1), plot = FALSE),
  list(source = "target", col = "target_tp",       label = "True Positive Rate for Target Effect",             ylim = c(0, 1), plot = TRUE),
  list(source = "target", col = "target_fn",       label = "False Negative Rate for Target Effect",            ylim = c(0, 1), plot = TRUE),
  list(source = "target", col = "target_fp",       label = "False Positive Rate for Target Effect",            ylim = c(0, 0.1), plot = TRUE),
  list(source = "target", col = "target_tn",       label = "True Negative Rate for Target Effect",             ylim = c(0, 1), plot = FALSE)
)

# helper: extract col names for a given source from specs
spec_cols <- function(src) sapply(Filter(function(s) s$source == src, specs), function(s) s$col)

true_effect_regions  <- which(ground_truth_true_effects > 0)
false_effect_regions <- which(ground_truth_true_effects <= 0)

results_by_rep            <- list()
results_by_region         <- list()
results_for_basis_study_region <- list()

for (sample_size in sim_params$sample_sizes) {
  ss      <- as.character(sample_size)
  pos_mat <- positives__replication[[ss]]
  eff_mat <- effect_sizes__replication[[ss]]
  pos_mat_univ <- positives_uncorr__replication[[ss]]
  eff_mat_univ <- effect_sizes_uncorr__replication[[ss]]

  # --- per-rep metrics: columns initialised from spec_cols("rep") ---
  rep_df <- data.frame(
    sample_size = sample_size,
    rep         = seq_len(sim_params$n_reps),
    setNames(replicate(length(spec_cols("rep")), rep(NA_real_, sim_params$n_reps), simplify = FALSE),
             spec_cols("rep"))
  )
  
  rep_df$num_tp__replication    <- rowSums(pos_mat[, true_effect_regions, drop = FALSE], na.rm = TRUE)
  rep_df$prop_tp__replication   <- rep_df$num_tp__replication / sim_params$n_regions
  rep_df$num_tp__basis     <- num_tp__basis[[ss]]
  rep_df$n_pos_above_target_effect_type__basis <- as.numeric(n_pos_above_target_effect_type__basis[[ss]])
  # if (sim_params$target_effect_type__from_basis == "min") {
  #   rep_df$num_tp_expect__from_basis <- rep_df$n_pos_above_target_effect_type__basis * (1-sim_params$alpha_fdr) # if significance thresholded, expect 5% of pos are FP
  # } else {
  #   rep_df$num_tp_expect__from_basis <- rep_df$n_pos_above_target_effect_type__basis * sim_params$targeted_power # no theory exists for adjustment, so researcher likely expects all to be TP  # TODO: would the average researcher adjust this expectation further?
  # }
  rep_df$num_tp_expect__from_basis <- num_tp_expect__from_basis[[ss]]
  rep_df$prop_above_target_effect_type__basis <- rep_df$n_pos_above_target_effect_type__basis / sim_params$n_regions
  rep_df$expect_v_actual_n_tp__based_on_basis <- rep_df$num_tp__replication / rep_df$num_tp_expect__from_basis
  this_idx__zero_not_na <- rep_df$n_pos_above_target_effect_type__basis == 0 & !is.na(rep_df$n_pos_above_target_effect_type__basis)
  rep_df$expect_v_actual_n_tp__based_on_basis[this_idx__zero_not_na] <- ifelse(rep_df$num_tp__replication[this_idx__zero_not_na]==0,1,0) # TODO: see below comment
  rep_df$expect_v_actual_n_tp__based_on_corr <- rep_df$num_tp__basis /expected_num_tp__corrected[[ss]]
  
  if (expected_num_tp__corrected[[ss]] == 0) {
    rep_df$expect_v_actual_n_tp__based_on_corr <- ifelse(num_tp__basis[[ss]]==0,1,0) # if match, set to 1 - TODO: hack to deal with studies having no actual effects; let's just take a simple difference between extents (divide by sim_params$n_regions)
  }
  rep_df$overlap         <- as.numeric(proportion_overlap_with_basis__replication[[ss]])
  rep_df$expected_n__replication <- expected_n_to_replicate_basis_effect[[ss]]
  total_pos_by_rep <- rowSums(pos_mat, na.rm = TRUE)
  false_pos_by_rep <- rowSums(pos_mat[, false_effect_regions, drop = FALSE], na.rm = TRUE)
  rep_df$fdr              <- ifelse(total_pos_by_rep > 0, false_pos_by_rep / total_pos_by_rep, NA_real_)
  
  results_by_rep[[ss]]   <- rep_df
  
  
  # --- per-region metrics: columns initialised from spec_cols("region") ---
  region_df <- data.frame(
    sample_size = sample_size,
    region      = seq_len(sim_params$n_regions),
    setNames(replicate(length(spec_cols("region")), rep(NA_real_, sim_params$n_regions), simplify = FALSE),
             spec_cols("region"))
  )
  positive_rate <- colMeans(pos_mat, na.rm = TRUE)
  region_df$tpr[true_effect_regions]    <- positive_rate[true_effect_regions] # TODO: adjust so it's not counting tpr for false_effect_regions
  region_df$type_m[true_effect_regions] <- sapply(seq_along(true_effect_regions), function(j) {
    mean(abs(eff_mat[, true_effect_regions[j]]) / ground_truth_true_effects[true_effect_regions[j]], na.rm = TRUE)
  })
  region_df$type_s[true_effect_regions] <- sapply(seq_along(true_effect_regions), function(j) {
    mean(sign(eff_mat[, true_effect_regions[j]]) != sign(ground_truth_true_effects[true_effect_regions[j]]), na.rm = TRUE)
  })
  results_by_region[[ss]] <- region_df

  # --- target-region metrics: columns initialized from spec_cols("target") ---
  basis_study_df <- data.frame(
    sample_size = sample_size,
    rep         = seq_len(sim_params$n_reps),
    setNames(replicate(length(spec_cols("target")), rep(NA, sim_params$n_reps), simplify = FALSE),
             spec_cols("target"))
  )
  target_idx_vec <- as.integer(target_effect_idx__basis[[ss]])
  for (this_rep in seq_len(sim_params$n_reps)) {
    this_idx <- target_idx_vec[this_rep]
    if (is.na(this_idx) || this_idx < 1 || this_idx > sim_params$n_regions) next

    is_true     <- ground_truth_true_effects[this_idx] > 0
    is_detected <- pos_mat[this_rep, this_idx]

    basis_study_df$target_true[this_rep]     <- is_true
    basis_study_df$target_detected[this_rep] <- is_detected
    basis_study_df$target_tp[this_rep]       <-  is_true &&  is_detected
    basis_study_df$target_fp[this_rep]       <- !is_true &&  is_detected
    basis_study_df$target_fn[this_rep]       <-  is_true && !is_detected
    basis_study_df$target_tn[this_rep]       <- !is_true && !is_detected
    
    positive_rate_univ <- colMeans(pos_mat_univ, na.rm = TRUE)
    basis_study_df$tpr_target_univ[this_rep]    <- positive_rate_univ[target_effect_idx__actual_mean]
    basis_study_df$type_m_target_univ[this_rep] <- abs(eff_mat_univ[this_rep, target_effect_idx__actual_mean]) / ground_truth_true_effects[target_effect_idx__actual_mean]
    basis_study_df$type_s_target_univ[this_rep] <-
      (eff_mat_univ[this_rep, target_effect_idx__actual_mean] < 0 & ground_truth_true_effects[target_effect_idx__actual_mean] > 0) |
        (eff_mat_univ[this_rep, target_effect_idx__actual_mean] > 0 & ground_truth_true_effects[target_effect_idx__actual_mean] < 0)
    
  }
  results_for_basis_study_region[[ss]] <- basis_study_df
}

region_metrics <- do.call(rbind, results_by_region)
rep_metrics <- do.call(rbind, results_by_rep)
basis_study_metrics <- do.call(rbind, results_for_basis_study_region)

# derive the summary column name from a spec entry:
#   region/rep → "<col>_mean"  (paired with "<col>_sd")
#   target     → "<col>"       (mean only)
#   fwer       → "fwer"        (special case)
summary_col_name <- function(spec) {
  # if (spec$source == "fwer")   return("fwer")
  if (spec$source == "target" || spec$source == "overall") return(spec$col)
  paste0(spec$col, "_mean")
}

# initialise summary sim_data frame — columns are derived from specs so they stay in sync automatically
summary_col_names <- c("sample_size", unlist(lapply(specs, function(spec) {
  mean_col <- summary_col_name(spec)
  # region and rep sources also get a paired _sd column
  if (spec$source %in% c("region", "rep")) c(mean_col, sub("_mean$", "_sd", mean_col)) else mean_col
})))
summary <- setNames(
  data.frame(matrix(NA_real_, nrow = length(sim_params$sample_sizes), ncol = length(summary_col_names))),
  summary_col_names
)
summary$sample_size <- sim_params$sample_sizes

# fill summary: loop over sample sizes, then over specs
for (i in seq_along(sim_params$sample_sizes)) {
  sample_size <- sim_params$sample_sizes[i]

  # subset each long table to this sample size
  rm_ss   <- region_metrics[region_metrics$sample_size == sample_size, ]
  rep_ss  <- rep_metrics[rep_metrics$sample_size == sample_size, ]
  targ_ss <- basis_study_metrics[basis_study_metrics$sample_size == sample_size, ]
  src     <- list(region = rm_ss, rep = rep_ss, target = targ_ss)

  for (spec in specs) {
    mean_col <- summary_col_name(spec)

    if (spec$source == "overall") {
      # special setup
      # fwer = fraction of reps where at least one false positive occurred
      pos_mat <- positives__replication[[as.character(sample_size)]]
      if (spec$col == "fwer") {
        summary$fwer[i] <- mean(rowSums(pos_mat[, false_effect_regions, drop = FALSE], na.rm = TRUE) > 0)
      } else if (spec$col == "num_tp__expected_corr") {
        summary$num_tp__expected_corr[i] <- expected_num_tp__corrected[[as.character(sample_size)]]
      }
        
    } else {
      vals <- as.numeric(src[[spec$source]][[spec$col]])
      summary[[mean_col]][i] <- mean(vals, na.rm = TRUE)
      # region and rep metrics also store sd across their respective units
      if (spec$source %in% c("region", "rep")) {
        summary[[sub("_mean$", "_sd", mean_col)]][i] <- sd(vals, na.rm = TRUE)
      }
    }
  }
}

plot_summary <- function(summary, metric, metric_label, y_limits, out_dir = "~/Desktop/sim/") {
  
  str_mean <- metric
  if (grepl("_mean$", metric)) {
    str_sd <- sub("_mean", "_sd", metric)
  } else {
    str_sd <- NULL
  }

  # preallocate df_plot
  df_plot <- data.frame(
    sample_size = summary$sample_size, ub = NA_real_, lb = NA_real_, point = NA_real_
  )
  
  # preprocess to avoid ghosting artifacts from too big data
  padding <- 2
  df_plot$point <- summary[, str_mean]
  df_plot$point[df_plot$point > (y_limits[2]+padding)] <- y_limits[2]+padding
  if (!is.null(str_sd)) {
    df_plot$ub <- summary[, str_mean] + summary[, str_sd]
    df_plot$lb <- summary[, str_mean] - summary[, str_sd]
    df_plot$ub[df_plot$ub > (y_limits[2]+padding)] <- y_limits[2]+padding
    df_plot$lb[df_plot$lb > (y_limits[2]+padding)] <- y_limits[2]+padding
    df_plot$lb[df_plot$lb < (y_limits[1]-padding)] <- y_limits[1]-padding
  }
  
  df_plot$x_num <- as.integer(factor(df_plot$sample_size))
  if (grepl("Corr", metric_label)) {
    x_label <- "Sample Size"
  } else {
    x_label <- "Sample Size of Basis Study"
  }

  p <- ggplot(df_plot, aes(x = x_num, y = point)) + # nolint
    geom_line(size = 0.5) +
    scale_x_continuous(breaks = df_plot$x_num, labels = df_plot$sample_size)
  
  # add ribbon if sd exists
  if (!is.null(str_sd)) {
    p <- p + geom_ribbon(aes(x = x_num, ymin = lb, ymax = ub), alpha = 0.2, inherit.aes = FALSE) # nolint
  }
    
  p <- p +
    coord_cartesian(ylim = y_limits) +
    labs(title = metric_label, x = x_label, y = metric_label) +
    theme_bw()

  ggsave(filename = paste0(out_dir, metric, ".png"), plot = p, width = 6, height = 4)
}

# plots — driven directly from specs, so adding a spec entry automatically adds a plot
for (spec in specs) {
  if (isTRUE(spec$plot)) {
    plot_summary(summary, summary_col_name(spec), spec$label, spec$ylim, out_dir = out_dir)
  }
}
# save summary variable if doesn't exist
# if (!file.exists(summary_file)) {
  saveRDS(summary, file = summary_file)
# }


