# Demo for generating Effect Size Project results manuscript figures

# set paths
source('/Users/steph/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/scripts/smn_pathsetter.R')
if (!exists('project_dir')) { # temp trigger to use pre-set paths
  print('Assuming script is in working directory and creating results in working directory.')
  script_path <- 'plot_results.R'
  data_path <- 'braineffex_data_09-24-25.RData' # name of file we will d/l from OSF
  output_dir <- 'crossbrain_effects/'
}

# source library
source(script_path)

# set params
estimate <- 'd' # only works for d currently, need to fix calculate_effex/effect_size/scripts/checker.R to work for r_sq
pooling_methods <- c('none', 'net')
motion_methods <- c('none','regression','threshold') # motion method for other tests
save_plots <- TRUE # whether to save the plots
get_data_from_OSF <- FALSE
osf_file_id <- 'g84tk'

# load data
if (get_data_from_OSF) {
  library(osfr)
  analysis_plan <- osf_retrieve_file("g84tk") %>% osf_download()
}
load(data_path)

# generate plots
for (pooling_method in pooling_methods) {
  for (motion_method in motion_methods) {
    print(paste0('Doing pooling method: ', pooling_method, ', motion: ', motion_method))
    combo_name <- paste0('pooling.', pooling_method, '.motion.',motion_method,'.mv.none') # note: automatically does mv.none->mv.multi in the background
    output_basedir <- paste0(output_dir, combo_name,'/')  # user-defined path to save the plots
    estimate_xb_effects(estimate, output_basedir, v, combo_name, save_plots)
  }
}

