do_simulation <- function(
  sim_params,
  effect_size_estimate_data,
  out_master_dir,
  ask_to_overwrite_if_exists = FALSE,
  default_overwrite = FALSE
) {

######### PER-ITERATION SETUP #########

# cross-brain area and cross-subject distribution parameters (Cohen's)
mu_crossbrain <- 0

# for parametric methods
n_groups <- 1 # note: this is hard-coded in the simulations
n_sides <- 1 # note: this is hard-coded in the simulations

# load sim_params for this outcome
load(effect_size_estimate_data) # load "res" from cross-brain effect estimation
var_crossbrain <- res[sim_params$outcome_category, "est"]
var_crosssubject <- res[sim_params$outcome_category, "phi2_est"]
sd_crossbrain <- sqrt(var_crossbrain)
sd_crosssubject <- sqrt(var_crosssubject)

# set files for this outcome
out_dir <- paste0(out_master_dir, sim_params$outcome_category, "/", sim_params$target_effect_type__from_basis, "_effect/")
sim_results_file <- paste0(out_dir, "sim_results.Rdata")
sim_params_file <- paste0(out_dir, "sim_params.rds")
summary_file <- paste0(out_dir, "summary_results.rds")

# make out dir for this outcome
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# calculate target based on number of regions
num_tp_target <- 0.05 * sim_params$n_regions


######### HELPER FUNCTIONS #########

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

simulate_data_ctrl <- function(n_subjects, n_regions, mu_crossbrain, sd_crossbrain, sd_crosssubject) {
  ground_truth_true_effects <- replicate(0, n_regions) # all effects truly 0 so difference with treatment gives specified effect size
  sim_data <- matrix(0, nrow = n_subjects, ncol = n_regions)
  for (i in 1:n_subjects) {
    subject_var <- rnorm(n_regions, mean = 0, sd = sd_crosssubject)
    sim_data[i, ] <- as.vector(ground_truth_true_effects + subject_var)
  }
  return(list(
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

print(paste0("Running: ", sim_params$outcome_category, " with target effect type from basis: ", sim_params$target_effect_type__from_basis))

######### SIMULATE DATA #########

# if sim_data exists, see whether current params are the same as previous
sim_data_exists <- FALSE
if (file.exists(sim_params_file)) {
  old_sim_params <- readRDS(sim_params_file)
  if (identical(names(old_sim_params), names(sim_params)) &&
      all(vapply(names(sim_params), function(param) {
        identical(old_sim_params[[param]], sim_params[[param]])
      }, logical(1)))) {
    sim_data_exists <- TRUE
  }
}

# if it exists, ask whether to replace
if (sim_data_exists) {
  if (ask_to_overwrite_if_exists) {
    run_sim <- readline(prompt = "Simulations using same params already exist. Replace? (yes/no; select yes if sim code has changed): ")
    run_sim <- tolower(run_sim) == "yes"
  } else {
    run_sim <- default_overwrite
  }
} else {
  run_sim <- TRUE
}


# Start Simulation

if (run_sim) {
  
print("Simulating data...")

# simulate master dataset
sim_info <- simulate_data(sim_params$n_subjects__gt, sim_params$n_regions, mu_crossbrain, sd_crossbrain, sd_crosssubject)
ground_truth_true_effects <- sim_info$ground_truth_true_effects  #sim_info[1:sim_params$n_regions]
sim_data <- sim_info$sim_data #matrix(sim_info[-(1:sim_params$n_regions)], nrow = sim_params$n_subjects__gt, ncol = sim_params$n_regions)

# if 2-sample test: simulate control dataset  - TODO: in progress

if (sim_params$n_regions %% 2 == 0) { # drop the min effect before mean to make things easier
  target_effect__actual_mean <- mean(ground_truth_true_effects[ground_truth_true_effects > min(ground_truth_true_effects > 0)], na.rm = TRUE)
} else {
  target_effect__actual_mean <- mean(ground_truth_true_effects[ground_truth_true_effects > 0])
}
# get number that is closest to but bigger than mean so it's a conservative comparison (easier for inflation-based planning to hit)
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
    
    # if 2-sample, subsample ctrl group also - TODO
    
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
        
    ### "Mean" Case: we're more optimistic and planning for the average effect size reported, even if it's biased by selection of significant effects
    } else if (sim_params$target_effect_type__from_basis == "mean") {
      
      if (sum(this_sig_mask) == 1) { # if only 1, choose that one
        target_effect_size__from_basis <- sig_effects_basis
      } else {
        target_effect_size__from_basis <- mean(sig_effects_basis, na.rm = TRUE)
      }

      n_pos_above_target_effect_type__basis[[as.character(sample_size)]][rep] <- sum(sig_effects_basis >= target_effect_size__from_basis)
      
      num_tp_expect__from_basis[[as.character(sample_size)]][rep] <- sum(sig_effects_basis >= target_effect_size__from_basis) * sim_params$targeted_power # conservative estimate when planning for target power to detect this mean effect size
      num_fp_expect__from_basis[[as.character(sample_size)]][rep] <- (sim_params$alpha_fdr/(1-sim_params$alpha_fdr)) * num_tp_expect__from_basis[[as.character(sample_size)]][rep] # 5 fp for every 95 tp
      
      num_tp_expect__in_basis[[as.character(sample_size)]][rep] <- sum(sig_effects_basis >= target_effect_size__from_basis) * (1-sim_params$alpha_fdr) # expect 5% of sig are FP (probably conservative with these larger-than-mean effects)
    
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
      effects_num <- as.numeric(this_res$effects)
      if (!is.finite(target_effect_size__from_basis) || all(!is.finite(effects_num))) {
        target_effect_idx__basis[[as.character(sample_size)]][rep] <- NA_integer_
      } else {
        target_effect_idx__basis[[as.character(sample_size)]][rep] <- which.min(abs(effects_num - target_effect_size__from_basis))
      }
    }
    
    # store planned sample size for replication
    if (sim_params$target_effect_type__from_basis == "same_sample") {
      expected_n_to_replicate_basis_effect[[as.character(sample_size)]][rep] <- sample_size
    } else {
      expected_n_to_replicate_basis_effect[[as.character(sample_size)]][rep] <- ceiling(pwr.t.test(power=sim_params$targeted_power, d = target_effect_size__from_basis, sig.level = sim_params$alpha_fdr, type = "one.sample", alternative = "greater")$n)
    }
    
  }
}



######### 1B. REPLICATION STUDY #########

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
  list = c("ground_truth_true_effects",  # must be saved so true_effect_regions / false_effect_regions are always consistent with the stored positives/effects
           "num_tp__basis", "n_pos_above_target_effect_type__basis", "num_tp_expect__from_basis", "num_fp_expect__from_basis", "num_tp_expect__in_basis",
           "target_effect_idx__basis", "expected_n_to_replicate_basis_effect",
           "positives__replication", "effect_sizes__replication", "positives_uncorr__replication", "effect_sizes_uncorr__replication", "proportion_overlap_with_basis__replication",
           "expected_tpr__corrected", "expected_num_tp__corrected", "expected_n_for_target_num_tp__corrected"),
  file = sim_results_file
)
saveRDS(sim_params, file = sim_params_file)
  
} else { # If skip sim, load previously saved results
  print("Simulations already exist - loading saved results...")
  load(sim_results_file)
  readRDS(sim_params_file)
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

  return(list(
    out_dir = out_dir,
    summary_file = summary_file,
    sim_results_file = sim_results_file,
    sim_params_file = sim_params_file,
    ground_truth_true_effects = ground_truth_true_effects,
    num_tp__basis = num_tp__basis,
    n_pos_above_target_effect_type__basis = n_pos_above_target_effect_type__basis,
    num_tp_expect__from_basis = num_tp_expect__from_basis,
    num_fp_expect__from_basis = num_fp_expect__from_basis,
    num_tp_expect__in_basis = num_tp_expect__in_basis,
    target_effect_idx__basis = target_effect_idx__basis,
    expected_n_to_replicate_basis_effect = expected_n_to_replicate_basis_effect,
    positives__replication = positives__replication,
    effect_sizes__replication = effect_sizes__replication,
    positives_uncorr__replication = positives_uncorr__replication,
    effect_sizes_uncorr__replication = effect_sizes_uncorr__replication,
    proportion_overlap_with_basis__replication = proportion_overlap_with_basis__replication,
    expected_tpr__corrected = expected_tpr__corrected,
    expected_num_tp__corrected = expected_num_tp__corrected,
    expected_n_for_target_num_tp__corrected = expected_n_for_target_num_tp__corrected
  ))
}



