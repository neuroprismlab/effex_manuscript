# Demo for results manuscript figures

# source plot_results.R library
scripts_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(scripts_dir,"/plot_results.R"))

# set params
estimate = 'd' # only works for d currently, need to fix calculate_effex/effect_size/scripts/checker.R to work for r_sq
pooling_methods = c('none', 'net')
motion_method_t = "threshold" # motion method for one-sample t-tests # NOTE: "threshold" refers to 0.1mm threshold, use "threshold2" for 0.2mm threshold
motion_method_not_t = "regression" # motion method for other tests
get_data_from_OSF <- FALSE
save_plots = TRUE # whether to save the plots

# setup dirs - can also define both directly
base_dir <- list.files('/Users', pattern = '^ste.*', full.names = TRUE, all.files = TRUE)
base_dir <- paste0(base_dir,'/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/')
if (!get_data_from_OSF) {
  data_dir <- paste0(base_dir, 'scripts/BrainEffeX_utils/inst/meta/v.RData') # user-specified - only needed if not 
}

# load data
if (get_data_from_OSF) {
  library(osfr)
  analysis_plan <- osf_retrieve_file("g84tk") %>% osf_download()
  data_path <- "braineffex_data_09-24-25.RData"
  load(data_path)
} else {
  load(data_dir)
}

# generate plots
for (pooling_method in pooling_methods) {
  print(paste0('Doing pooling method: ', pooling_method))
  combo_name <- paste0('pooling.', pooling_method, '.motion.',motion_method_not_t)
  output_basedir <- paste0(base_dir, 'manuscript/figures/plots/results/curve_fit/',combo_name,'/')  # user-defined path to save the plots
  plot_results(estimate, output_basedir, v, motion_method_t, motion_method_not_t, save_plots)
}

