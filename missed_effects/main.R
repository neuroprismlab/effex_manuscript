####################################################
#
# Simulations for spatial extent of detections
#
####################################################


######### SET PATHS & GLOBAL CONFIG #########

# paths
# project_dir <- "/Users/steph/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/"
project_dir <- "/Users/stephanienoble/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/"
scripts_dir <- paste0(project_dir, "scripts/crossbrain_effects/missed_effects/")
out_master_dir <- paste0(project_dir, "manuscript/figures/plots/missed_effects/")
crossbrain_effect_estimator_script <- paste0(project_dir, "scripts/crossbrain_effects/crossbrain_effect_estimator.R")
effect_size_estimate_data <- paste0(project_dir, "manuscript/figures/plots/crossbrain_effects/pooling.none.motion.regression.mv.none/point_res.Rdata")

# libraries
library(MASS)
library(ggplot2)
library(reshape2)
library(pwr)
source(crossbrain_effect_estimator_script)
source(paste0(scripts_dir, "do_simulation.R"))
source(paste0(scripts_dir, "plot_results_single.R"))
source(paste0(scripts_dir, "plot_results.R"))

all_outcome_categories <- c("psychological", "physical", "task activation", "task connectivity")
all_target_effect_type__from_basis <- c("same_sample", "max", "mean_without_selection")

run_outcome_categories <- all_outcome_categories
run_target_effect_type__from_basis <- all_target_effect_type__from_basis
do_test <- FALSE
if (do_test) {
  run_outcome_categories <- run_outcome_categories[1]
  run_target_effect_type__from_basis <- run_target_effect_type__from_basis[1]
}

# overwrite params
ask_to_overwrite_if_exists <- FALSE
default_overwrite <- FALSE

# simulation params (global defaults; outcome/effect type are injected inside the loop)
run_simulation <- TRUE
sim_params <- list(
  n_subjects__gt = 10000,
  n_regions = 500,
  sample_sizes = c(25, 50, 100, 500, 1000, 5000),
  n_reps = 200,
  alpha_fdr = 0.05,
  targeted_power = 0.8
)

# plotting params
cats <- c("psychological", "physical", "task activation", "task connectivity")
cat_colors <- RColorBrewer::brewer.pal(length(cats), "Set1")
cat_colors[c(1, 2)] <- cat_colors[c(2, 1)]
names(cat_colors) <- cats
transparency_main <- 0.6
transparency_overlay <- 0.5
text_size <- 20
ticks_size <- 24


######### SIMULATION + PLOTS PER OUTCOME CATEGORY / BASIS TYPE #########

if (run_simulation) {
  for (outcome_category in run_outcome_categories) {
    for (target_effect_type__from_basis in run_target_effect_type__from_basis) {
    
      # Per-combination params
      sim_params$outcome_category <- outcome_category
      sim_params$target_effect_type__from_basis <- target_effect_type__from_basis
    
      sim_output <- do_simulation(sim_params = sim_params, effect_size_estimate_data = effect_size_estimate_data,
        out_master_dir = out_master_dir, ask_to_overwrite_if_exists = ask_to_overwrite_if_exists,
        default_overwrite = default_overwrite
      )
    
      plot_results_single(
        sim_params = sim_params, sim_output = sim_output, cat_colors = cat_colors,
        out_master_dir = out_master_dir, text_size = text_size, ticks_size = ticks_size,
        transparency_main = transparency_main, transparency_overlay = transparency_overlay
      )
    
    }
  }
}


######### MULTI-CATEGORY PLOTS #########

plot_results(
  all_target_effect_type__from_basis = all_target_effect_type__from_basis,
  all_outcome_categories = all_outcome_categories,
  cat_colors = cat_colors, out_master_dir = out_master_dir,
  text_size = text_size, ticks_size = ticks_size,
  transparency_main = transparency_main, transparency_overlay = transparency_overlay
)
