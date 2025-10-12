#' Create results manuscript plots
#' Authors: Stephanie Noble & Hallee Shearer
#'
#' This function is the master plotter for making effect size plots for the results manuscript.
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @import lme4
#' @import pwr
#' @param estimate
#' @ param fn_basedir Directory to save results
#' @param v_data list of studies in the BrainEffeX format (see OSF link: https://doi.org/10.17605/OSF.IO/CWNJD)
#' @param motion_method_t motion correction method to use for studies with one sample t-tests (default: "threshold")
#' @param motion_method_not_t motion correction method to use for studies without one sample t-tests (default: "regression")
#' @param save_plots Whether to save plots to files (default: TRUE)
#'
#' @return Saves plots to specified output directory
#'
#' @examples
#' # Example usage
#' \dontrun{
#' plot_results(estimate = 'd', fn_basedir = 'results/', v_data = v)
#' }

# ------------- MAIN -------------------------

plot_results <- function(estimate = 'd', fn_basedir, v_data, motion_method_t = "threshold", motion_method_not_t = "regression", combo_name, save_plots = TRUE) {

## Setup
  
# libraries
library(dplyr) # TODO: don't think we need to re-import here etc
library(ggplot2)
library(ggrepel)
library(lme4)
library(pwr) # for sample size calculations
library(tidyr)
library(RColorBrewer)
library(metafor)
  
# set add'l params
n_large_threshold <- 900
cats <- c("psychological", "physical", "task activation", "task connectivity")
cat_colors <- RColorBrewer::brewer.pal(length(cats), "Set1")
n_pts <- 1000

# make output directory if it doesn't exist
if (!dir.exists(fn_basedir)) {
  print(paste0('Creating output directory: ', fn_basedir))
  dir.create(fn_basedir, recursive = TRUE)
}

## Get summaries
# - also making separate data frame with conservative estimates (facilitates reuse of later functions on cons est)

# 1. using individual studies
summary_data <- get_study_summaries(v_data$data, v_data$study, estimate, combo_name)
summary_data_cons <- summary_data
summary_data_cons$mean <- summary_data_cons$mean_cons
summary_data_cons$sd <- summary_data_cons$sd_cons
summary_data_cons$char_mag <- summary_data_cons$char_mag_cons

# 2. using meta-analysis results

## Extract additional info from meta-analysis before can summarize

# first, for meta, category is inexplicably group_level - rename 
v_data$meta_category$study$category <- v_data$meta_category$study$group_level

# Assign all relevant datasets, n's, and study names to meta_category$study
v_data$meta_category$study$dataset <- vector("list", length(v_data$meta_category$study$name))
v_data$meta_category$study$each_n <- vector("list", length(v_data$meta_category$study$name))
v_data$meta_category$study$n <- vector("list", length(v_data$meta_category$study$name)) # number of unique subjects
v_data$meta_category$study$n1 <- vector("list", length(v_data$meta_category$study$name)) # number of unique subjects
v_data$meta_category$study$n2 <- vector("list", length(v_data$meta_category$study$name)) # number of unique subjects
v_data$meta_category$study$included_study_names <- vector("list", length(v_data$meta_category$study$name))
v_data$meta_category$study$n_studies <- integer(length(v_data$meta_category$study$name))

# For each meta-analysis, get all studies included, their associated datasets, and sample sizes
for (i in seq_along(v_data$meta_category$study$name)) {
  
  # get category and ref from meta_name
  meta_name <- v_data$meta_category$study$name[i]
  parts <- strsplit(meta_name, "_reference_")[[1]]
  meta_cat <- parts[1]
  meta_ref <- parts[2]
  
  matches <- which(v_data$study$category == meta_cat & v_data$study$ref == meta_ref)
  
  included_study_names <- v_data$study$name[matches]
  v_data$meta_category$study$included_study_names[[i]] <- included_study_names
  v_data$meta_category$study$dataset[[i]] <- v_data$study$dataset[matches]
  v_data$meta_category$study$overarching_category[[i]] <- unique(summary_data$overarching_category[matches])
  v_data$meta_category$study$each_n[[i]] <- summary_data$n[match(included_study_names, summary_data$name)]
  v_data$meta_category$study$n_studies[[i]] <- length(included_study_names)
  
  v_data$meta_category$study$n <- vector("list", length(v_data$meta_category$study$name))
  # for (i in seq_along(v_data$meta_category$study$name)) {
    meta_datasets <- v_data$meta_category$study$dataset[[i]]
    unique_datasets <- unique(meta_datasets)
    unique_n <- 0
    for (ds in unique_datasets) {
      n_vals <- unlist(v_data$meta_category$study$each_n[[i]][v_data$meta_category$study$dataset[[i]] == ds])
      if (length(n_vals) > 0) {
        unique_n <- unique_n + max(n_vals, na.rm = TRUE)
      }
      # else do nothing (skip if no n)
    }
    
    # go through each field in v_data$meta_category$data
    for (field in names(v_data$meta_category$data[[i]])) {
      if (is.list(v_data$meta_category$data[[i]][[field]])) {
        # if list, subset to only those included in this meta-analysis
        v_data$meta_category$data[[i]][[field]]$n <- unique_n
      }
    }
    
    # add orig stat type
    v_data$meta_category$study$orig_stat_type[[i]] <- unique(v_data$study$orig_stat_type[matches])
    
    # collect common stat_type for this meta-analysis (same as datasets collation)
    stat_types_in_meta <- v_data$study$orig_stat_type[matches]
    orig_stat_type <- unique(stat_types_in_meta)
    if (length(orig_stat_type) == 1) {
      v_data$meta_category$study$orig_stat_type[[i]] <- orig_stat_type
    } else {
      # If multiple stat types, use the most common one or a combination label
      stat_type_counts <- table(stat_types_in_meta)
      v_data$meta_category$study$orig_stat_type[[i]] <- names(which.max(stat_type_counts))
    }
    
  # }
}
names(v_data$meta_category$study$dataset) <- v_data$meta_category$study$name
names(v_data$meta_category$study$each_n) <- v_data$meta_category$study$name
names(v_data$meta_category$study$included_study_names) <- v_data$meta_category$study$name
names(v_data$meta_category$study$n_studies) <- v_data$meta_category$study$name

# summary - meta
summary_data__meta <- get_study_summaries(v_data$meta_category$data, v_data$meta_category$study,estimate, combo_name)
summary_data_cons__meta <- summary_data__meta
summary_data_cons__meta$mean <- summary_data_cons__meta$mean_cons
summary_data_cons__meta$sd <- summary_data_cons__meta$sd_cons
summary_data_cons__meta$char_mag <- summary_data_cons__meta$char_mag_cons

## Extra info

# get total unique subjects across all studies
all_datasets <- unique(summary_data$dataset)
all_datasets <- all_datasets[all_datasets != "hcp_voxel"] # remove hcp_voxel, which is a subset of hcp_shen_268
all_ns <- numeric(length(all_datasets))
for (j in seq_along(all_datasets)) {
  ds <- all_datasets[j]
  ns <- summary_data$n[summary_data$dataset == ds]
  if (length(ns) > 0) {
    all_ns[j] <- max(ns, na.rm = TRUE)
  } else {
    all_ns[j] <- 0
  }
}
total_n <- sum(all_ns, na.rm = TRUE)


## Plot
  
plot_extra <- FALSE # TODO - tmp

# get point est and conservative est
res <- estimate_params(summary_data, summary_data__meta,  n_pts, "Point Estimates", paste0(fn_basedir,'point'))
# get multivariate effect sizes
res_mv <- estimate_params(summary_data, summary_data__meta,  n_pts, "Multivariate Effect Sizes", paste0(fn_basedir,'mv_est'), plot_type = "mv")

# Make density plots
sigmas_master <- plot_densities(res, res_mv,  n_pts, fn_basedir, cats, cat_colors, save_plots)

# Req'd n plots:
# - original
required_n_df <- make_required_n_df(n_pts, sigmas_master, do_mv = FALSE)
plot_required_n_panel(required_n_df, do_mv = FALSE, cat_colors,fn_basedir)

# - multivariate
required_n_df_mv <- make_required_n_df(n_pts, sigmas_master, res_mv = res_mv, do_mv = TRUE)
plot_required_n_panel(required_n_df_mv, do_mv = TRUE, cat_colors,fn_basedir)



if (plot_extra) {
  # conservative and large n
  res_cons <- estimate_params(summary_data_cons, summary_data_cons__meta,  n_pts, "Conservative Estimates", paste0(fn_basedir,'extra/cons'))
  res_large <- estimate_params(summary_data[summary_data$n > n_large_threshold,], summary_data__meta[summary_data__meta$n > n_large_threshold,], n_pts, "Point Estimates (n > 900)", paste0(fn_basedir,'extra/point_n900'))
}


} # function






# ------------- HELPER FUNCTIONS -------------------------


get_characteristic_magnitude <- function(mean, sd) {
  # search over values to find characteristic magnitude where 68% of values are < |mag|
  mag <- 0
  # first pass
  while ((pnorm(mag, mean = mean, sd = sd) - pnorm(-mag, mean = mean, sd = sd)) < 0.68) {
    mag <- mag + 0.01
  }
  # refine estimate
  mag <- mag - 0.01
  while ((pnorm(mag, mean = mean, sd = sd) - pnorm(-mag, mean = mean, sd = sd)) < 0.68) {
    mag <- mag + 0.001
  }
  return(mag)
}


# helper
d_se <- function(d, n1, n2 = NULL) {
  if (is.null(n2)) { # one-sample
    se <- sqrt(1 / n1 + (d^2 / (2 * n1)))
  } else { # two-sample
    se <- sqrt((n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2))))
  }
  return(se)
}

# helper
r_sq_se <- function(r_sq, n) {
  r <- sqrt(r_sq)
  se_r <- sqrt((1 - r^2) / (n - 2));
  se <- se_r^2
  return(se)
}



# ------------- ORGANIZING, FITTING, AND PLOTTING -------------------------


########### Summarize studies ########### 

get_study_summaries <- function(data, study, estimate, combo_name) {
  
  # preallocate
  d_name <- numeric(length(data))
  d_mean <- numeric(length(data))
  d_sd <- numeric(length(data))
  d_characteristic_mag <- numeric(length(data))
  d_cons_mean <- numeric(length(data))
  d_cons_sd <- numeric(length(data))
  d_characteristic_mag_cons <- numeric(length(data))
  d_n <- numeric(length(data))
  d_mv <- numeric(length(data))
  d_mv_lb <- numeric(length(data))
  d_mv_ub <- numeric(length(data))
  vi <- numeric(length(data))
  vi_char_mag <- numeric(length(data))
  vi_mv <- numeric(length(data))
  
  if (estimate == "d") {
    ci_lb <- "sim_ci_lb"
    ci_ub <- "sim_ci_ub"
  } else if (estimate == "r_sq") {
    ci_lb <- "r_sq_sim_ci_lb"
    ci_ub <- "r_sq_sim_ci_ub"
  }
  
  for (i in seq_along(data)) {
    
    # print(paste0('Processing study: ', study$name[i]))
    
    # 1. Get point estimate mean and sd across vars
    
    # to get d, check whether dim is null (nested list)
    if (is.null(dim(data[[i]][[combo_name]][[estimate]]))) {
      d <- data[[i]][[combo_name]][[estimate]]
    } else {
      d <- data[[i]][[combo_name]][[estimate]][1,]
    }
    
    d_mean[i] <- mean(d)
    d_sd[i] <- sd(d)
    d_characteristic_mag[i] <- get_characteristic_magnitude(d_mean[i], d_sd[i])
    
    # 2. Get conservative estimate mean and sd across vars
    # from prep_data_for_plot.R
    # omitting na check in mean and sd and downsampling
    
    # unlist sim CIs if list
    if (is.list(data[[i]][[combo_name]][[ci_lb]])) {
      data[[i]][[combo_name]][[ci_lb]] <- unlist(data[[i]][[combo_name]][[ci_lb]])
      data[[i]][[combo_name]][[ci_ub]] <- unlist(data[[i]][[combo_name]][[ci_ub]])
    }
    
    na_idx <- is.na(data[[i]][[combo_name]][[estimate]]) | is.na(data[[i]][[combo_name]][[ci_lb]]) | is.na(data[[i]][[combo_name]][[ci_ub]])
    data[[i]][[combo_name]][[estimate]] <- data[[i]][[combo_name]][[estimate]][!na_idx]
    data[[i]][[combo_name]][[ci_lb]] <- data[[i]][[combo_name]][[ci_lb]][!na_idx]
    data[[i]][[combo_name]][[ci_ub]] <- data[[i]][[combo_name]][[ci_ub]][!na_idx]
    
    # sort data from smallest to largest effect size
    sorted_indices <- order(data[[i]][[combo_name]][[estimate]])
    sorted_estimate <- data[[i]][[combo_name]][[estimate]][sorted_indices]
    sorted_upper_bounds <- data[[i]][[combo_name]][[ci_ub]][sorted_indices]
    sorted_lower_bounds <- data[[i]][[combo_name]][[ci_lb]][sorted_indices]
    
    sorted_cons_estimate <- ifelse((abs(sorted_lower_bounds) > abs(sorted_upper_bounds)),
                                   ifelse((sorted_upper_bounds < 0),
                                          round(sorted_upper_bounds, 2), 0),
                                   ifelse((sorted_lower_bounds > 0),
                                          round(sorted_lower_bounds, 2), 0))
    
    
    d_cons_mean[i] <- mean(sorted_cons_estimate)
    d_cons_sd[i] <- sd(sorted_cons_estimate)
    d_characteristic_mag_cons[i] <- get_characteristic_magnitude(d_cons_mean[i], d_cons_sd[i])
    
    # 3. Get multivariate effect size and bounds if available
    
    mv_combo_name <- gsub("mv.none", "mv.multi", combo_name) # get multivariate combo name
    all_combos <- names(data[[i]])
    mv_combo_idx <- grep(mv_combo_name, all_combos)
    # if doesn't exist, replace motion.xxx. with motion.threshold.
    if (length(mv_combo_idx) == 0) {
      mv_combo_name_alt <- gsub("motion\\.[^\\.]+\\.", "motion.threshold.", mv_combo_name)
      mv_combo_idx <- grep(mv_combo_name_alt, all_combos)
    }
    mv_combo_name <- all_combos[mv_combo_idx]
    d_mv[i] <- data[[i]][[mv_combo_name]]$d
    d_mv_lb[i] <- data[[i]][[mv_combo_name]]$sim_ci_lb
    d_mv_ub[i] <- data[[i]][[mv_combo_name]]$sim_ci_ub
    
    # 4. Get sample size and variances for meta-regression
    
    if (study$orig_stat_type[[i]] == "t2") {
      if (!is.null(data[[i]][[combo_name]]$n1)) {
        d_n[i] <- data[[i]][[combo_name]]$n1 + data[[i]][[combo_name]]$n2
        # get var for meta
        if (estimate == "d") {
          vi[i] <- d_se(d_sd, data[[i]][[combo_name]]$n1, data[[i]][[combo_name]]$n2)^2
          vi_char_mag[i] <- d_se(d_characteristic_mag[i], data[[i]][[combo_name]]$n1, data[[i]][[combo_name]]$n2)^2
          vi_mv[i] <- d_se(d_mv[i], data[[i]][[combo_name]]$n1, data[[i]][[combo_name]]$n2)^2
        } else if (estimate == "r_sq") {
          vi[i] <- r_sq_se(d_sd, data[[i]][[combo_name]]$n1 + data[[i]][[combo_name]]$n2)^2
          vi_char_mag[i] <- r_sq_se(d_characteristic_mag[i], data[[i]][[combo_name]]$n1 + data[[i]][[combo_name]]$n2)^2
          vi_mv[i] <- r_sq_se(d_mv[i], data[[i]][[combo_name]]$n1 + data[[i]][[combo_name]]$n2)^2
        }
      } else {
        # TODO: for the meta-analytic results, we should really just get the results from the confidence intervals,
        # since these come directly from the meta-analysis. For now, we don't use the variances previously calculated
        # from the meta-analysis for the meta-regression below, so we will ignore
        # (add helper alongside other se estimators below)
        d_n[i] <- data[[i]][[combo_name]]$n
        vi[i] <- NA
        vi_char_mag[i] <- NA
        vi_mv[i] <- NA
      }
    } else if (study$orig_stat_type[[i]] == "t" || study$orig_stat_type[[i]] == "r") {
    # if (!is.null(data[[i]][[combo_name]]$n)) {
      d_n[i] <- data[[i]][[combo_name]]$n
      
      if (estimate == "d") {
        if (study$orig_stat_type[[i]] == "r") { # treat as 2-sample t-test
          vi[i] <- d_se(d_sd, d_n[i]/2, d_n[i]/2)^2
          vi_char_mag[i] <- d_se(d_characteristic_mag[i], d_n[i]/2, d_n[i]/2)^2
          vi_mv[i] <- d_se(d_mv[i], d_n[i]/2, d_n[i]/2)^2
        } else { # normal 1-sample t-test
          vi[i] <- d_se(d_sd, d_n[i])^2
          vi_char_mag[i] <- d_se(d_characteristic_mag[i], d_n[i])^2
          vi_mv[i] <- d_se(d_mv[i], d_n[i])^2
        }
      } else if (estimate == "r_sq") {
        vi[i] <- r_sq_se(d_sd, d_n[i])^2
        vi_char_mag[i] <- r_sq_se(d_characteristic_mag[i], d_n[i])^2
        vi_mv[i] <- r_sq_se(d_mv[i], d_n[i])^2
      }
    } else {
      d_n[i] <- NA
    }
    
    # d_biased_sd[i] <- sd(d) * sqrt((length(d) - 1) / length(d)) # biased sd - doesn'd change results
    # TODO: catch normality test results
    # test for normality (usually light- or heavy-tailed, some approximately normal)
    # d_subset <- d[seq(1, length(d), length.out = 5000)] # max 5000 variables for shapiro test
    # d_normal_test[i] <- shapiro.test(d_subset)$p.value # p-value < 2e-16 -> not normal
    # qqnorm(d, pch=20, cex=0.5); qqline(d) # visualize
    
    
  }
  
  # set up final data frame
  
  df <- data.frame(name = names(data), mean = d_mean, sd = d_sd, char_mag = d_characteristic_mag, mean_cons = d_cons_mean, sd_cons = d_cons_sd, char_mag_cons = d_characteristic_mag_cons, n = d_n, category = as.factor(study$category), dataset = I(study$dataset), ref = study$ref, orig_stat_type = study$orig_stat_type, mv = d_mv, mv_lb = d_mv_lb, mv_ub = d_mv_ub, vi_sd = vi, vi_char_mag = vi_char_mag, vi_mv = vi_mv, stringsAsFactors = FALSE)
  
  # for meta: if exists, add n_studies
  if ("n_studies" %in% colnames(study)) {
    df$n_studies <- study$n_studies
  }
  
  # add reference type to dataset to distinguish data used for act from data used for FC (encompasses substantially different processing that should be nested within ref category)
  df <- df %>%
    mutate(dataset = paste(dataset, ref, sep = "_"))
  
  df <- df %>%
    mutate(overarching_category = case_when(
      category %in% c("biometric", "sex (demographic)", "age (demographic)") ~ "physical",
      category %in% c("cognitive", "psychiatric") ~ "psychological",
      # category == "cognitive (task)" ~ "task (within-sub)",
      category == "cognitive (task)" & !grepl("voxel", ref) ~ "task connectivity",
      category == "cognitive (task)" & grepl("voxel", ref) ~ "task activation",
      TRUE ~ "other"
    ))
  df$overarching_category <- as.factor(df$overarching_category)
  
  return(df)
  
}


########### Estimate Parameters & Plot Fits ########### 

# Function for plotting effect sizes (mean, sd, n) for each study - point est and conservative  
# Add a flag to toggle between sqrt(n) and 1/sqrt(n) for plotting

estimate_params <- function(df, df_meta, n_pts, main_title, fn, plot_type = "sd") {
  
  print(paste0('Fitting lines for ', main_title))
  
  # params
  nmax_plt <- 10000000
  nmax <- 1000000000
  points_only_dir <- 'points_only/'
  plot_fitted_lines <- TRUE
  n_pts <- 1000
  use_char_mag <- FALSE
  plot_with_inv_sqrt_n <- FALSE # FALSE for sqrt(n), TRUE for 1/sqrt(n)
  
  # Determine y variable and settings based on plot_type
  if (plot_type == "sd") {
    y_var <- "sd"
    y_label <- "SD(d) across brain areas"
    fit_lines <- TRUE
    y_limits <- c(-0.05, 0.5)
  } else if (plot_type == "mv") {
    y_var <- "mv"
    y_label <- "Multivariate effect size"
    y_limits <- NULL  # Let ggplot auto-scale
  } else {
    stop("plot_type must be either 'sd' or 'mv'")
  }
  
  # sort by sample size
  df <- df[order(df$n, decreasing = FALSE), ]
  
  for (add_meta in c(TRUE, FALSE)) {
    
    if (add_meta) {
      meta_str <- '__meta'
      nmax_plt <- max(df_meta$n)
      alpha <- 0.15 # make non-meta points more transparent
    } else {
      meta_str <- ''
      nmax_plt <- max(df$n)
      alpha <- 0.7
    }
    
    # Fit params & plot fits
    
    if (plot_fitted_lines) {
      fitlines_str <- '__fits'
      
      # interpolate predictions for all n
      n_seq <- data.frame(n=seq(0, nmax_plt, length.out=n_pts))
      n_seq[length(n_seq$n)+1,] <- nmax # ensure max n is included
      
      # set up for each overarching category
      unique_cats <- unique(df$overarching_category)
      # if 'task activation' exists in unique_cats, make it come right after 'task connectivity'
      if ("task activation" %in% unique_cats) {
        unique_cats <- c(setdiff(unique_cats, "task connectivity"), "task connectivity")
        unique_cats <- c(setdiff(unique_cats, "task activation"), "task activation")
      }
      predicted_y <- vector("list", length(unique_cats))
      names(predicted_y) <- unique_cats
      
      do_by_cat <- FALSE
      
      if (do_by_cat) {
        
        # preallocate results
        res <- vector("list", length(unique_cats))
        
        for (cat in unique_cats) {
          
          # setup df
          df_cat <- df[df$overarching_category == cat, ]
          # df_cat <- df[df$overarching_category == cat & !is.na(df$mv), ]
          n_obs <- nrow(df_cat)
          n_levels <- length(unique(df_cat$dataset))
          
          if (cat != "task activation") { # special procedure for activation
            
            # # estimate fits with nesting by dataset, unless only one unique dataset
            # if (length(unique(df_cat$dataset)) > 2) {
            #   
            #   # fit
            #   if (plot_type == "mv") {
            #     fit_cat <- lmer(mv ~ 1 + (1|dataset), data = df_cat)
            #     d <- as.matrix(model.matrix(~ 1, data = n_seq))
            #   } else {
            #     if (use_char_mag) {
            #       fit_cat <- lmer(char_mag ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            #     } else {
            #       fit_cat <- lmer(sd ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            #     }
            #     # fit_cat <- lmer(char_mag ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            #     d <- model.matrix(~ I(1/sqrt(n)), data = n_seq)
            #   }
            #   
            #   # Use fixed effects for prediction
            #   preds <- as.vector(d %*% fixef(fit_cat) )
            #   preds_se <- sqrt(diag(d %*% vcov(fit_cat) %*% t(d)))
            #   lwr <- preds - 1.96 * preds_se
            #   upr <- preds + 1.96 * preds_se
            #   predicted_y[[cat]] <- cbind(fit = preds, lwr = lwr, upr = upr)
            #   
            # } else {
            
            if (plot_type == "mv") {
              fit_cat <- lm(mv ~ 1, data=df_cat)
            } else {
              if (use_char_mag) {
                fit_cat <- lm(char_mag ~ I(1/sqrt(n)), data=df_cat)
              } else {
                fit_cat <- lm(sd ~ I(1/sqrt(n)), data=df_cat)
              }
            }
            predicted_y[[cat]] <- predict(fit_cat, n_seq, interval="confidence")
            # }
            
          } else {
            
            # use slope from task connectivity to estimate intercept
            
            if (plot_type == "mv") { # same as above
              # fit_cat2 <- lm(mv ~ 1, data=df_cat)
              # # For mv, just use fit_cat2 estimates directly
              # d <- model.matrix(~ 1, data = n_seq)
              # preds <- as.vector(d %*% coef(fit_cat2))
              # preds_se <- sqrt(diag(d %*% vcov(fit_cat2) %*% t(d)))
              
              fit_cat <- lm(mv ~ 1, data=df_cat)
              predicted_y[[cat]] <- predict(fit_cat, n_seq, interval="confidence")
              
            } else {
              
              # Get slope from previous fit
              # slope <- fixef(fit_cat)["I(1/sqrt(n))"]
              slope <- fit_cat$coefficients["I(1/sqrt(n))"]
              
              if (use_char_mag) { # different - fit intercept with slope fixed from connectivity
                fit_cat2 <- lm(I(char_mag - I(slope * 1/sqrt(n))) ~ 1, data=df_cat)
              } else {
                fit_cat2 <- lm(I(sd - I(slope * 1/sqrt(n))) ~ 1, data=df_cat)
              }
              d <- model.matrix(~ I(1/sqrt(n)), data = n_seq)
              
              # Create combined coefficients: intercept from fit_cat2, slope from fit_cat
              # also, corresponding variances for standard errors
              combined_coef <- c(coef(fit_cat2)["(Intercept)"], slope)
              names(combined_coef) <- c("(Intercept)", "I(1/sqrt(n))")
              vcov_combined <- vcov(fit_cat)
              vcov_combined["(Intercept)", "(Intercept)"] <- vcov(fit_cat2)["(Intercept)", "(Intercept)"]
              
              # Use combined coefficients for prediction
              preds <- as.vector(d %*% combined_coef)
              preds_se <- sqrt(diag(d %*% vcov_combined %*% t(d)))
              lwr <- preds - 1.96 * preds_se
              upr <- preds + 1.96 * preds_se
              predicted_y[[cat]] <- cbind(fit = preds, lwr = lwr, upr = upr)
              
              fit_cat <- fit_cat2 # for getting Intercept est, lwr, upr below
            }
          }
          
          est <- summary(fit_cat)$coefficients[,"Estimate"]
          lwr <- summary(fit_cat)$coefficients[,"Estimate"] - 1.96 * summary(fit_cat)$coefficients[,"Std. Error"]
          upr <- summary(fit_cat)$coefficients[,"Estimate"] + 1.96 * summary(fit_cat)$coefficients[,"Std. Error"]
          
          # max_val[[cat]] <- predicted_y[[cat]][length(n_seq$n),]
          # res[[cat]] <- c(beta[[cat]], max_val[[cat]])
          
          if (plot_type == "mv" | cat == "task activation") {
            res[[cat]] <- c(est = est, lwr = lwr, upr = upr)
          } else {
            res[[cat]] <- data.frame(est = est["(Intercept)"], lwr = lwr["(Intercept)"], upr = upr["(Intercept)"], row.names = NULL)
          }
          
        }
        
        # Combine all categories into a tidy data frame
        res <- do.call(rbind, res)
        
        
      } else { # nesting within category, dataset x category
        
        # fit meta-analysis nesting studies by dataset and overarching category
        
        # setup grouping variables
        df$dataset_nested <- interaction(df$overarching_category, df$dataset, drop = TRUE)
        # df$study_nested <- interaction(df$overarching_category, df$dataset, df$name, drop = TRUE)
        
        if (plot_type == "mv") {
          fit_all <- rma.mv(yi = mv, 
                            V = vi_mv,  # approximate variance from CI
                            random = ~ 1 | overarching_category/dataset_nested,
                            # random = ~ 1 | overarching_category/dataset_nested/study_nested,
                            data = df,
                            method = "REML")
        } else {
          if (use_char_mag) {
            fit_all <- rma.mv(yi = char_mag, 
                              V = vi_char_mag,  # approximate variance
                              mods = ~ I(1/sqrt(n)),
                              random = ~ 1 | overarching_category/dataset_nested,
                              data = df,
                              method = "REML")
          } else {
            fit_all <- rma.mv(yi = sd, 
                              V = vi_sd,  # approximate variance
                              mods = ~ I(1/sqrt(n)),
                              random = ~ 1 | overarching_category/dataset_nested,
                              data = df,
                              method = "REML")
          }
        }
        
        # Extract estimates and confidence intervals from meta-analysis
        # Extract category-specific random effects (BLUPs)
        random_effects <- ranef(fit_all)
        category_effects <- random_effects$overarching_category
        
        # Create results with category-specific intercepts
        res <- vector("list", length(unique_cats))
        names(res) <- unique_cats
        predicted_y <- vector("list", length(unique_cats))
        names(predicted_y) <- unique_cats
        
        for (cat in unique_cats) {
          # Get category-specific intercept (overall + random effect)
          if (cat %in% rownames(category_effects)) {
            cat_intercept <- fit_all$beta[1] + category_effects[cat, "intrcpt"]
            cat_intercept_se <- sqrt(fit_all$vb[1,1] + category_effects[cat, "se"]^2)
          } else {
            # If category not found, use overall intercept
            cat_intercept <- fit_all$beta[1]
            cat_intercept_se <- sqrt(fit_all$vb[1,1])
          }
          
          if (plot_type == "mv") {
            # For intercept-only model
            res[[cat]] <- data.frame(
              est = cat_intercept,
              lwr = cat_intercept - 1.96 * cat_intercept_se,
              upr = cat_intercept + 1.96 * cat_intercept_se,
              row.names = paste0(cat, "_intercept")
            )
            
            # Create predicted values (constant line for mv)
            preds <- rep(cat_intercept, length(n_seq$n))
            preds_se <- rep(cat_intercept_se, length(n_seq$n))
            lwr <- preds - 1.96 * preds_se
            upr <- preds + 1.96 * preds_se
            predicted_y[[cat]] <- cbind(fit = preds, lwr = lwr, upr = upr)
            
          } else {
            # For models with moderators (intercept and slope)
            # Note: slope is the same for all categories (fixed effect)
            slope <- fit_all$beta[2]
            slope_se <- sqrt(fit_all$vb[2,2])
            
            # Only include intercept in results (not slope)
            res[[cat]] <- data.frame(
              est = cat_intercept,
              lwr = cat_intercept - 1.96 * cat_intercept_se,
              upr = cat_intercept + 1.96 * cat_intercept_se,
              row.names = paste0(cat, "_intercept")
            )
            
            # Create predicted values with category-specific intercept
            preds <- cat_intercept + slope * (1/sqrt(n_seq$n))
            # Standard errors for predictions (more complex with random effects)
            X_pred <- cbind(1, 1/sqrt(n_seq$n))
            preds_se_fixed <- sqrt(diag(X_pred %*% fit_all$vb %*% t(X_pred)))
            preds_se <- sqrt(preds_se_fixed^2 + cat_intercept_se^2)
            lwr <- preds - 1.96 * preds_se
            upr <- preds + 1.96 * preds_se
            predicted_y[[cat]] <- cbind(fit = preds, lwr = lwr, upr = upr)
          }
        }
        
        # Combine results into single dataframe
        res <- do.call(rbind, res)
        
      }
      
    } else {
      # If not plotting fitlines, change output folder to 'points_only' dir
      fn_plot <- file.path(dirname(fn), points_only_dir, basename(fn))
      if (!dir.exists(dirname(fn_plot))) {
        dir.create(dirname(fn_plot), recursive = TRUE)
      }
      res <- NULL
    }
    
    # plot
    if (plot_with_inv_sqrt_n) {
      df$x_plot <- 1/sqrt(df$n)
      df_meta$x_plot <- 1/sqrt(df_meta$n)
      x_label <- "1 / sqrt(n)"
      x_limits <- c(0, max(df$x_plot, df_meta$x_plot, na.rm=TRUE))
    } else {
      df$x_plot <- sqrt(df$n)
      df_meta$x_plot <- sqrt(df_meta$n)
      x_label <- "sqrt(n)"
      x_limits <- c(0, max(df$x_plot, df_meta$x_plot, na.rm=TRUE))
    }
    cats <- levels(df$overarching_category)
    color_map <- setNames(RColorBrewer::brewer.pal(length(cats), "Set1"), cats)
    
    p <- ggplot(df, aes_string(x = "x_plot", y = y_var, color = "overarching_category")) +
      geom_point(size = 1.5, alpha = alpha, stroke = 0) +
      scale_color_manual(values = color_map) +
      labs(title = main_title,
           x = x_label,
           y = y_label,
           color = "Category") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.key.size = unit(0.7, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3)
    
    if (!is.null(y_limits)) {
      p <- p + scale_y_continuous(limits = y_limits, expand = c(0, 0))
    }
    
    if (plot_fitted_lines) {
      for (cat in unique_cats) {
        pred_mat <- predicted_y[[cat]]
        if (!is.null(pred_mat) && is.matrix(pred_mat) && all(c("fit","lwr","upr") %in% colnames(pred_mat)) && nrow(pred_mat) == length(n_seq$n)) {
          pred_df <- data.frame(
            x_plot = if (plot_with_inv_sqrt_n) 1/sqrt(n_seq$n) else sqrt(n_seq$n),
            fit = pred_mat[,"fit"],
            lwr = pred_mat[,"lwr"],
            upr = pred_mat[,"upr"],
            overarching_category = cat
          )
          p <- p +
            geom_line(data = pred_df, aes(x = x_plot, y = fit, color = overarching_category), linewidth = 1, alpha = alpha) +
            geom_line(data = pred_df, aes(x = x_plot, y = lwr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha) +
            geom_line(data = pred_df, aes(x = x_plot, y = upr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha)
        }
      }
      
      if (add_meta && nrow(df_meta) > 0) {
        df_meta$label <- paste0(gsub("_reference_", " (", df_meta$name), ")\n", df_meta$n_studies, " ", ifelse(df_meta$n_studies == 1, "study", "studies"))
        p <- p +
          geom_point(data = df_meta, aes_string(x = "x_plot", y = y_var, color = "overarching_category"), shape = 8, size = 2) +
          geom_text_repel(data = df_meta, aes_string(x = "x_plot", y = y_var, label = "label"), size = 2.2, segment.size = 0.3, force = 10, max.overlaps = Inf)
      }
      
      fitlines_suffix <- if (plot_fitted_lines) "__fits" else ""
      fn_plot <- fn
      this_fn <- paste0(fn_plot, fitlines_suffix, meta_str, '.pdf')
      if (save_plots) {
        ggsave(this_fn, plot = p, width = 5, height = 4)
      } else {
        show(p)
      }
    }
  }
  
  return(res)
}


###########  Make Estimated Density Plots ########### 

plot_densities <- function(res, res_mv,  n_pts, fn_basedir, cats, cat_colors, save_plots = TRUE) {

  print('Making density plots')
  # cats <- unique(rownames(res))
  # cat_colors <- RColorBrewer::brewer.pal(length(cats), "Set1")
  
  # params
  do_horizontal_panels <- TRUE
  
  for (do_mv in c(FALSE, TRUE)) {
    if (do_mv) {
      mv_suffix <- '_mv'
      xlim <- c(0, 3.5)
    } else {
      mv_suffix <- ''
      xlim <- c(-0.8, 0.8)
    }
    d <- seq(xlim[1], xlim[2], length.out = n_pts)
    
  
  # preallocate
  density_list <- list()
  if (!do_mv) {
    sigmas_master <- NULL
  }
  
  for (i in seq_along(cats)) {
    cat <- cats[i]
    
    if (do_mv) {
      
      for (name in names(res_mv[cat, ])) {
        mu <- res_mv[cat, name]
        
        y <- rep(0, length(d))
        closest_idx <- which.min(abs(d - mu))
        y[closest_idx] <- 10
        
        density_list[[length(density_list)+1]] <- data.frame(
          d = d,
          density = y,
          category = cat,
          sigma_type = name,
          overarching_category = cat
        )
      }
      
    } else {
      sigmas <- c(res[cat,"lwr"], res[cat,"est"], res[cat,"upr"])
      # print(sigmas)
      names(sigmas) <- c('lb','est','ub')
      sigmas_master <- rbind(sigmas_master, sigmas)
      
      # If not mv, make y as dnorm for each sigma (lb, est, ub)
      for (name in names(sigmas)) {
        
        s <- sigmas[name]
        
        if (s >= 0) {
          y <- dnorm(d, mean = 0, sd = s)
          density_list[[length(density_list)+1]] <- data.frame(
            d = d,
            density = y,
            category = cat,
            sigma_type = name,
            overarching_category = cat
          )
        } else {
          # point mass at 0
          y <- rep(0, length(d))
          # closest_idx <- which.min(abs(d))
          y[ceiling(length(d)/2)] <- 1
          density_list[[length(density_list)+1]] <- data.frame(
            d = d,
            density = y,
            category = cat,
            sigma_type = name,
            overarching_category = cat
          )
        }
      }
    }
  }
  
  if (!do_mv) {
    sigmas_master <- as.data.frame(sigmas_master)
    rownames(sigmas_master) <- cats
    # TODO: can probably just use res, instead of recreating and renaming sigmas_master
  }
  
  
  density_df <- do.call(rbind, density_list)

  y_max <- max(density_df$density)
  
  # Overlapping densities colored by overarching category, legend inset
  p_density <- ggplot(density_df, aes(x = d, y = density, color = overarching_category, linetype = sigma_type)) +
    geom_line(size = 1) +
    scale_color_manual(values = cat_colors) +
    labs(title = "Density Plot by Overarching Category", x = "Cohen's d", y = "Density", color = "Category", linetype = "Sigma Type") +
    theme_bw() +
    coord_cartesian(xlim = xlim, ylim = c(0, y_max)) +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key.size = unit(0.7, "lines")
    )
  if (save_plots) {
    ggsave(paste0(fn_basedir, 'density',mv_suffix,'.pdf'), p_density, width = 5, height = 4)
  } else {
    print(p_density)
  }
  
  # Panel densities
  
  if (do_horizontal_panels) {
    nrow <- 1
    width = 4 * length(cats)
    height = 4
  } else {
    nrow <- length(cats)
    width = 5
    height = 4 * length(cats)
  }
  density_df$fill <- density_df$sigma_type=="est"
  
  # Separate panels per category, colored by overarching category, legend inset
  p_density_panel <- ggplot() +
    # ribbon only for "est" rows
    geom_ribbon(
      data = subset(density_df, sigma_type == "est"),
      aes(x = d, ymin = 0, ymax = density, fill = overarching_category),
      alpha = 0.35,
      colour = NA,
      inherit.aes = FALSE,
      show.legend = FALSE   # hide separate fill legend if you prefer
    ) +
    # lines for all rows (color + linetype)
    geom_line(
      data = density_df,
      aes(x = d, y = density, color = overarching_category, linetype = sigma_type),
      size = 0.9
    ) +
    facet_wrap(~category, nrow = nrow, scales = "free_y") +
    scale_color_manual(values = cat_colors, name = "Category") +
    scale_fill_manual(values = cat_colors, guide = "none") +  # keep colors consistent, hide fill legend
    labs(title = "Density Curves by Category", x = "Cohen's d", y = "Density", color = "Category", linetype = "Sigma Type") +
    theme_bw() +
    coord_cartesian(xlim = xlim) +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key.size = unit(0.7, "lines")
    )
  
  # p_density_panel <- ggplot(density_df, aes(x = d, y = density, color = overarching_category, fill = overarching_category, linetype = sigma_type)) +
  #   geom_ribbon(aes(ymin = 0, ymax = density,
  #                   alpha = ifelse(sigma_type == "est", 0.35, 0)),
  #               colour = NA, inherit.aes = TRUE, show.legend = FALSE) +
  #   geom_line(size = 1) +
  #   facet_wrap(~category, nrow = nrow, scales = "free_y") +
  #   scale_color_manual(values = cat_colors) + # remap line colors
  #   scale_fill_manual(values = cat_colors) + # remap fill colors
  #   labs(title = "Density Curves by Category", x = "Cohen's d", y = "Density", color = "Category", linetype = "Sigma Type") +
  #   theme_bw() +
  #   coord_cartesian(xlim = xlim) +
  #   theme(
  #     legend.position = c(0.02, 0.98),
  #     legend.justification = c("left", "top"),
  #     legend.background = element_rect(fill = "white", color = "grey80"),
  #     legend.key.size = unit(0.7, "lines")
  #   )
  if (save_plots) {
    ggsave(paste0(fn_basedir, 'density__panels',mv_suffix,'.pdf'), p_density_panel, width = width, height = height)
  } else {
    print(p_density_panel)
  }
  
  # Save sigmas_master to file
  if (save_plots) {
    if (do_mv) {
      write.csv(res_mv, file=paste0(fn_basedir, 'param_mus',mv_suffix,'.csv'), row.names=TRUE)
    } else {
      write.csv(sigmas_master, file=paste0(fn_basedir, 'param_sigmas',mv_suffix,'.csv'), row.names=TRUE)
    }
  }
  
  }
  
  return(sigmas_master)
}

##### REQUIRED N PLOTS #####

# get req'd n and bin
make_required_n_df <- function(n_pts, sigmas_master, res_mv = NULL, do_mv = FALSE) {
  
  if (do_mv) {
    mv_suffix <- '_mv'
    xlim <- c(0, 6)
  } else {
    mv_suffix <- ''
    xlim <- c(-0.8, 0.8)
  }
  d <- seq(xlim[1], xlim[2], length.out = n_pts)

  
  # Safe wrapper for pwr.t.test to handle edge cases
  safe_pwr <- function(dd, test_type) {
    if (is.na(dd) || abs(dd) < 1e-6) {
      return(Inf)  # Very small effect sizes require infinite sample size
    }
    tryCatch({
      pwr.t.test(power = 0.8, d = abs(dd), sig.level = 0.05/2, type = test_type)$n
    }, error = function(e) {
      return(Inf)  # Return Inf if calculation fails
    })
  }
  
  n_bins <- c(0, 25, 50, 100, 500, 1000, 5000, 50000, 300000, Inf)
  bin_labels <- paste(head(n_bins, -1), n_bins[-1]-1, sep = "–")
  required_n_df <- NULL
  
  # get cats from res
  cats <- rownames(sigmas_master)
  for (cat in cats) {
    
    # run power
    if (grepl("task", cat)) {
      test_type <- "one.sample"
    } else {
      test_type <- "two.sample"
    }
    
    if (do_mv) {
      y <- rep(0, length(d))
      if (!is.null(res_mv) && cat %in% rownames(res_mv)) {
        closest_idx <- which.min(abs(d - res_mv[cat, "est"])) # TODO: this should use sigmas_master for mv instead; but then again, we should be able to just use res instead of sigmas throughout
        y[closest_idx] <- 1
      }
      
      n_detect <- sapply(d, function(dd) safe_pwr(dd, test_type))
      # n_detect_two <- sapply(d, function(dd) safe_pwr(dd, "two.sample"))
      bin_indices <- cut(n_detect, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
      
    } else {
      sigmas <- sigmas_master[cat, ]
      y <- dnorm(d, mean = 0, sd = as.numeric(sigmas["ub"]))
      
      
      n_detect <- sapply(d, function(dd) safe_pwr(dd, test_type))
      bin_indices <- cut(n_detect, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
    }
    
    binned_sums <- tapply(y, bin_indices, sum, na.rm = TRUE)
    binned_sums <- binned_sums / sum(binned_sums, na.rm = TRUE)
    binned_sums[!is.na(binned_sums)] <- cumsum(binned_sums[!is.na(binned_sums)])
    tmp_df <- data.frame(
      bin = bin_labels,
      cumulative_proportion = as.numeric(binned_sums),
      overarching_category = cat
    )
    required_n_df <- rbind(required_n_df, tmp_df)
  }
  
  required_n_df$bin <- factor(required_n_df$bin, levels = bin_labels, ordered = TRUE)
  return(required_n_df)
}


# plot req'd n

plot_required_n_panel <- function(df, do_mv = FALSE, cat_colors, fn_basedir) {
  
  do_horizontal_panels <- TRUE
  
  title <- if (do_mv) "Required n for 80% Power by Category (Multivariate)" else "Required n for 80% Power by Category"
  filename <- if (do_mv) 'sample_size__panels_mv.pdf' else 'sample_size__panels.pdf'
  # Add line segment from y=0 to y=y[1] at x=1 if y[1] != 0 for each category
  segment_df <- df %>% group_by(overarching_category) %>% filter(row_number() == 1 & cumulative_proportion[1] != 0)
  
  if (do_horizontal_panels) {
    nrow <- 1
    width = 4 * length(cat_colors)
    height = 4
  } else {
    nrow <- length(cat_colors)
    width = 5
    height = 4 * length(cat_colors)
  }
  
  p <- ggplot(df, aes(x = bin, y = cumulative_proportion, group = overarching_category, color = overarching_category)) +
    geom_ribbon(aes(ymin = 0, ymax = cumulative_proportion, fill = overarching_category),
                alpha = 0.35, inherit.aes = TRUE, show.legend = FALSE) +
    geom_line(size = 1) +
    scale_color_manual(values = cat_colors) +
    scale_fill_manual(values = cat_colors) + # remap fill colors
    facet_wrap(~overarching_category, nrow = nrow, scales = "free_y") +
    labs(title = title, x = "Minimum sample size bin", y = "Cumulative Proportion of Effects", color = "Category") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key.size = unit(0.7, "lines")
    )
  if (nrow(segment_df) > 0) {
    p <- p + geom_segment(data = segment_df, aes(x = 1, xend = 1, y = 0, yend = cumulative_proportion, color = overarching_category), inherit.aes = FALSE, size = 1)
  }
  if (save_plots) {
    ggsave(paste0(fn_basedir, filename), p, width = width, height = height)
  } else {
    print(p)
  }
}