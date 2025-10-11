# Demo for generating Effect Size Project results manuscript figures

# set paths
source('/Users/stephanienoble/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/scripts/smn_pathsetter.R')
if (!exists('project_dir')) { # temp trigger to use pre-set paths
  print('Assuming script is in working directory and creating results in working directory.')
  script_path <- 'plot_results.R'
  data_path <- 'braineffex_data_09-24-25.RData' # name of file we will d/l from OSF
  output_dir <- 'results/'
}

# source library
source(script_path)

# set params
estimate = 'd' # only works for d currently, need to fix calculate_effex/effect_size/scripts/checker.R to work for r_sq
pooling_methods = c('none', 'net')
motion_method_t = "threshold" # motion method for one-sample t-tests # NOTE: "threshold" refers to 0.1mm threshold, use "threshold2" for 0.2mm threshold
# TODO: ^critical - this should not be used except for multivariate one-sample t-test, and n should be updated accordingly
motion_method_not_t = "regression" # motion method for other tests
save_plots = TRUE # whether to save the plots
get_data_from_OSF = FALSE
osf_file_id <- 'g84tk'

# load data
if (get_data_from_OSF) {
  library(osfr)
  analysis_plan <- osf_retrieve_file("g84tk") %>% osf_download()
}
load(data_path)

# generate plots
for (pooling_method in pooling_methods) {
  print(paste0('Doing pooling method: ', pooling_method))
  combo_name <- paste0('pooling.', pooling_method, '.motion.',motion_method_not_t,'.mv.none')
  output_basedir <- paste0(output_dir, combo_name,'/')  # user-defined path to save the plots
  plot_results(estimate, output_basedir, v, motion_method_t, motion_method_not_t, combo_name, save_plots)
}

