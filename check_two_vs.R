# check a new v file against an old v file
# I loaded my (Hallee's) v file, and named it v_h
# loaded Steph's v file, and named it v_s

# first make sure all study names and brain mask names are lower case for both vs
names(v_h$data) = tolower(names(v_h$data))
names(v_s$data) = tolower(names(v_s$data))
names(v_h$brain_masks) = tolower(names(v_h$brain_masks))
names(v_s$brain_masks) = tolower(names(v_s$brain_masks))
v_h$study$name = tolower(v_h$study$name)
v_s$study$name = tolower(v_s$study$name)

motion_methods = c("none", "regression", "threshold")
pooling_methods = c("none", "net")

studies_in_both = intersect(names(v_h$data), names(v_s$data))

estimates = c("d", "r_sq", "sim_ci_lb", "sim_ci_ub", "r_sq_sim_ci_lb", "r_sq_sim_ci_ub")

list_of_bad_studies = c()
list_of_bad_combos = c()

# check the first few values of a given study
for (study in studies_in_both) {
  for (motion in motion_methods) {
    for (pooling in pooling_methods) {
      this_motion_idx = grepl(paste0("motion\\.",motion, "\\."), names(v_h$data[[study]]))
      this_pooling_idx = grepl(paste0("pooling\\.", pooling, "\\."), names(v_h$data[[study]]))
      these_combos_idx = this_motion_idx & this_pooling_idx
      these_combo_names = names(v_h$data[[study]])[these_combos_idx]
      
      for (combo in these_combo_names) {
        for (estimate in estimates) {
          
          if(!isTRUE(all.equal(as.vector(v_h$data[[study]][[combo]][[estimate]]), as.vector(v_s$data[[study]][[combo]][[estimate]]), tolerance = 0.001))) {
            if (!any(grepl(study, list_of_bad_studies))) {
              list_of_bad_studies = c(list_of_bad_studies, study)
            }
            if (!any(grepl(combo, list_of_bad_combos))) {
              list_of_bad_combos = c(list_of_bad_combos, combo)
            }
            print(paste0("study: ", study, "combo: ", combo))
            print(paste0(estimate, " values do not match for this study!"))
            print(v_h$data[[study]][[combo]][[estimate]][v_h$data[[study]][[combo]][[estimate]] == v_s$data[[study]][[combo]][[estimate]]])
            print(v_h$data[[study]][[combo]][[estimate]][v_h$data[[study]][[combo]][[estimate]] == v_s$data[[study]][[combo]][[estimate]]])
          } 
        }
      }
    }
  }
}


# repeat for meta-analysis results

motion_methods = c("none", "regression", "threshold")
pooling_methods = c("none", "net")

studies_in_both = intersect(names(v_h$meta_category$data), names(v_s$meta_category$data))

estimates = c("d", "r_sq", "sim_ci_lb", "sim_ci_ub", "r_sq_sim_ci_lb", "r_sq_sim_ci_ub")

list_of_bad_studies = c()
list_of_bad_combos = c()

# check the first few values of a given study
for (study in studies_in_both) {
  for (motion in motion_methods) {
    for (pooling in pooling_methods) {
      this_motion_idx = grepl(paste0("motion\\.",motion, "\\."), names(v_h$data[[study]]))
      this_pooling_idx = grepl(paste0("pooling\\.", pooling, "\\."), names(v_h$data[[study]]))
      these_combos_idx = this_motion_idx & this_pooling_idx
      these_combo_names = names(v_h$meta_category$data[[study]])[these_combos_idx]
      print(study)
      for (combo in these_combo_names) {
        for (estimate in estimates) {
          
          if(!isTRUE(all.equal(as.vector(v_h$meta_category$data[[study]][[combo]][[estimate]]), as.vector(v_s$meta_category$data[[study]][[combo]][[estimate]]), tolerance = 0.001))) {
            if (!any(grepl(study, list_of_bad_studies))) {
              list_of_bad_studies = c(list_of_bad_studies, study)
            }
            if (!any(grepl(combo, list_of_bad_combos))) {
              list_of_bad_combos = c(list_of_bad_combos, combo)
            }
            print(paste0("study: ", study, "combo: ", combo))
            print(paste0(estimate, " values do not match for this study!"))
            print(v_h$meta_category$data[[study]][[combo]][[estimate]][v_h$data[[study]][[combo]][[estimate]] == v_s$meta_category$data[[study]][[combo]][[estimate]]])
            print(v_h$meta_category$data[[study]][[combo]][[estimate]][v_h$data[[study]][[combo]][[estimate]] == v_s$meta_category$data[[study]][[combo]][[estimate]]])
          } 
        }
      }
    }
  }
}

