#' Estimate cross-brain effects
#' Authors: Stephanie Noble & Hallee Shearer
#'
#' Estimate and plot cross-brain effect size estimates from BrainEffeX data
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @import lme4
#' @import pwr
#' @param estimate
#' @ param fn_basedir Directory to save results
#' @param v_data list of studies in the BrainEffeX format (see OSF link: https://doi.org/10.17605/OSF.IO/CWNJD)
#' @param save_plots Whether to save plots to files (default: TRUE)
#'
#' @return Saves plots to specified output directory
#'
#' @examples
#' # Example usage
#' \dontrun{
#' estimate_xb_effects(estimate = 'd', fn_basedir = 'results/', v_data = v)
#' }

# ------------- MAIN -------------------------

estimate_xb_effects <- function(estimate = 'd', fn_basedir, v_data, combo_name, use_bayesian_fit = FALSE, save_plots = TRUE) {
  
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
  library(patchwork)
  library(bayesmeta) # Bayes regression one
  
  # set add'l params
  n_large_threshold <- 900
  cats <- c("psychological", "physical", "task activation", "task connectivity")
  cat_colors <- RColorBrewer::brewer.pal(length(cats), "Set1")
  cat_colors[c(1,2)] <- cat_colors[c(2,1)]
  names(cat_colors) <- cats
  n_pts <- 10000
  
  # Centralize plotting and export sizes/appearance in one place.
  plot_params <- list(
    estimate_fits = list(
      base_width = 5,
      base_height = 4
    ),
    density = list(
      do_horizontal_panels = TRUE,
      xlim_annotate = c(-0.2, 0.2),
      xbreaks = c(-2, 1, -0.5, -0.2, 0, 0.2, 0.5, 1, 2),
      axis_text_size = 20,
      axis_title_size = 20,
      base_width = 4.2,
      base_height = 3.5
    ),
    power = list(
      do_horizontal_panels = TRUE,
      base_width = 4.2,
      single_col_width = 4.5,
      base_height = 3.5,
      axis_text_size = 20,
      axis_title_size = 20
    )
  )
  
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
  summary_data_cons$var_xv <- summary_data_cons$var_xv_cons
  summary_data_cons$var_xv__emp <- summary_data_cons$var_xv__emp_cons
  
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
  summary_data_cons__meta$var_xv <- summary_data_cons__meta$var_xv_cons
  summary_data_cons__meta$var_xv__emp <- summary_data_cons__meta$var_xv__emp_cons
  
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
  
  
  ## Estimate & Plot
  
  plot_extra <- FALSE # TODO - tmp
  
  # Estimate effect sizes & plot param estimation plot
  # - mass univariate (corrected cross-brain distribution)
  res_fn_basename <- paste0(fn_basedir,'point')
  # res_fn <- paste0(res_fn_basename,'_res.Rdata')
  # if (file.exists(res_fn)) { # first try to load file if exists
  #   load(res_fn)
  # } else { ## df, df_meta, n_pts, main_title, fn, plot_type = "crossvariable", plot_params = NULL
  all_res <- estimate_params(summary_data, summary_data__meta,  n_pts, "Parameter Estimation Plot: Cross-Brain Effects", res_fn_basename, use_bayesian_fit = use_bayesian_fit, plot_params = plot_params)
  #res is only est lwr upr
  res <- all_res[c("est", "lwr", "upr")]
  phi2 <- all_res[c("phi2_est", "phi2_lwr", "phi2_upr")]
  # }
  # - multivariate
  res_fn_mv_basename <- paste0(fn_basedir,'mv_est')
  # res_fn_mv <- paste0(res_fn_mv_basename,'_res.Rdata')
  # if (file.exists(res_fn_mv)) {
  #   res_mv <- get(load(res_fn_mv))
  # } else {
  res_mv <- estimate_params(summary_data, summary_data__meta,  n_pts, "Parameter Estimation Plot: Multivariate Effects", res_fn_mv_basename, plot_type = "mv", use_bayesian_fit = use_bayesian_fit, plot_params = plot_params)
  # }
  
  # Make density plots
  sigmas_master <- plot_densities(res, res_mv,  n_pts, fn_basedir, cat_colors, save_plots, plot_params = plot_params)
  
  # Power plots
  do_other_power_plots <- TRUE # TODO: temporary
  if (do_other_power_plots) {
    # - mass univariate
    results_uv <- get_average_power(sigmas_master, do_mv = FALSE)
    avg_power <- results_uv$avg_power
    proportion_detectable <- results_uv$proportion_detectable
    
    plot_average_power(avg_power, do_mv = FALSE, cat_colors,fn_basedir, save_plots = save_plots, plot_params = plot_params)
    plot_proportion_detectable(proportion_detectable, do_mv = FALSE, cat_colors,fn_basedir, save_plots = save_plots, plot_params = plot_params)
    
    # - multivariate
    results_uv <- get_average_power(sigmas_master, res_mv = res_mv, do_mv = TRUE)
    avg_power_mv <- results_uv$avg_power
    proportion_detectable_mv <- results_uv$proportion_detectable
    
    plot_average_power(avg_power_mv, do_mv = TRUE, cat_colors,fn_basedir, save_plots = save_plots, plot_params = plot_params)
    plot_proportion_detectable(proportion_detectable_mv, do_mv = TRUE, cat_colors,fn_basedir, save_plots = save_plots, plot_params = plot_params)
  }
  
  # Power mismatch plots
  plot_proportion_difference(sigmas_master, phi2, cat_colors, fn_basedir)
  
  # # Req'd n plots:
  # # - mass univariate
  # required_n_df <- make_required_n_df(n_pts, sigmas_master, do_mv = FALSE)
  # plot_required_n_panel(required_n_df, do_mv = FALSE, cat_colors,fn_basedir)
  # 
  # # - multivariate
  # required_n_df_mv <- make_required_n_df(n_pts, sigmas_master, res_mv = res_mv, do_mv = TRUE)
  # plot_required_n_panel(required_n_df_mv, do_mv = TRUE, cat_colors,fn_basedir)
  
  
  
  if (plot_extra) {
    # conservative and large n
    res_cons <- estimate_params(summary_data_cons, summary_data_cons__meta,  n_pts, "Conservative Estimates", paste0(fn_basedir,'extra/cons'), use_bayesian_fit = use_bayesian_fit, plot_params = plot_params)
    res_large <- estimate_params(summary_data[summary_data$n > n_large_threshold,], summary_data__meta[summary_data__meta$n > n_large_threshold,], n_pts, "Point Estimates (n > 900)", paste0(fn_basedir,'extra/point_n900'), use_bayesian_fit = use_bayesian_fit, plot_params = plot_params)
  }
  
  
} # function






# ------------- HELPER FUNCTIONS -------------------------

# empirical sd counterpart

get_characteristic_magnitude_for_sd <- function(mean, sd) {
  # find empirical counterpart to standard deviation in describing distribution spread
  # by searching over values to find characteristic magnitude where 68% of values are < |mag|
  sd_emp <- 0
  # first pass
  while ((pnorm(sd_emp, mean = mean, sd = sd) - pnorm(-sd_emp, mean = mean, sd = sd)) < 0.68) {
    sd_emp <- sd_emp + 0.01
  }
  # refine estimate
  sd_emp <- sd_emp - 0.01
  while ((pnorm(sd_emp, mean = mean, sd = sd) - pnorm(-sd_emp, mean = mean, sd = sd)) < 0.68) {
    sd_emp <- sd_emp + 0.001
  }
  return(sd_emp)
}

# effect size standard error functions

d_se <- function(d, n1, n2 = NULL) {
  if (is.null(n2)) { # one-sample
    se <- sqrt(1 / n1 + (d^2 / (2 * n1)))
  } else { # two-sample
    se <- sqrt((n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2))))
  }
  return(se)
}

r_sq_se <- function(r_sq, n) {
  r <- sqrt(r_sq)
  se_r <- sqrt((1 - r^2) / (n - 2));
  se <- se_r^2
  return(se)
}

# power functions - provided by Thomas Nichols
#                   2-sided extensions added by Steph

# pi0  : proportion of true nulls
# delta: noncentrality parameter = d * sqrt(n)
# sigma_delta: standard deviation of noncentrality parameter across tests - sigma_d * sqrt(n)
# alphaFDR : target FDR level (e.g. 0.05)
# pp is p-value threshold for rejection of test (test-specific Type I error)
#   uncorrected: pp=alpha
#   Bonferroni: pp=alpha/k
#   FDR: pp=BHthresh

# CDF of non-null p-values under H1
# Usage:
#   For FDR, BHpower(pi0, alphaFDR, d*sqrt(n))
#   For uncorrected, power=F1(alphaFDR,d*sqrt(n))
#   For Bonferroni, power=F1(alphaFDR/k,d*sqrt(n))
F1 <- function(pp, delta, sigma_delta=0, n_groups, n_sides)  {
  n_sides * (1 - pnorm((qnorm(1 - pp/2) - delta/n_groups) / sqrt(1 + (sigma_delta/n_groups)^2)))
}

# BH-FDR threshold anticipated by the process approach to FDR
BHthresh <- function(pi0, alphaFDR, delta, sigma_delta=0, n_groups, n_sides) {
  pi1 <- 1 - pi0
  
  # FDR Power Implicit Equation, solved for pp (p')
  #    pi0 * pp + pi1 * F1(t) = pp / alphaFDR
  # pp is the long-run p-value threshold for this setting, and power
  # is then usual power at this threshold.
  f <- function(t) {
    pi0 * t + pi1 * F1(t,delta,sigma_delta,n_groups,n_sides) - t / (alphaFDR)
  }
  
  ## Trivial solution t = 0 always exists.
  ## Look for a nontrivial root in (0, alphaFDR), if it exists.
  lower <- 1e-10
  upper <- min(alphaFDR - 1e-10, 1 - 1e-10)
  
  if (f(upper) * f(lower) > 0) {
    ## No sign change -> only solution is t* = 0
    pp <- 0
  } else {
    pp <- uniroot(f, c(lower, upper))$root
  }
  
  pp
}

# Average BH power (per non-null) for given pi0, delta, alphaFDR
BHpower <- function(pi0, alphaFDR, delta, sigma_delta=0, n_groups, n_sides) {
  if (!is.finite(delta) || !is.finite(sigma_delta)) return(1) # catch case where effects infinitely big so guaranteed to be detected
  pp <- BHthresh(pi0, alphaFDR, delta, sigma_delta, n_groups, n_sides)
  if (pp == 0) return(0)
  F1(pp, delta, sigma_delta, n_groups, n_sides)
}

# additional power functions - Steph

# additional definitions
# beta: Type II error (1-power)
# n_groups (1 or 2): 1- or 2-sample test
# n_sides (1 or 2): 1- or 2-sided test (note: 1-sided assumes positive tail; multiply data by -1 manually for negative tail)

# Proportion of tests with above "adequate" (1-beta) power at significance level alpha
proportion_detectable <- function(pp, beta, delta, sigma_delta=0, n_groups=1, n_sides=2)  {
  if (n_sides == 2) { pp <- pp / 2 } # if 2-sided, significance threshold/2
  d_star <- qnorm(1 - pp) - qnorm(beta)
  n_sides*(1 - pnorm((d_star - delta/n_groups) / sigma_delta/n_groups))
}

BH_proportion_detectable <- function(pi0, alphaFDR, beta, delta, sigma_delta=0, n_groups=1, n_sides=2) {
  if (!is.finite(delta) || !is.finite(sigma_delta)) return(1) # catch case where effects infinitely big so guaranteed to be detected
  pp <- BHthresh(pi0, alphaFDR, delta, sigma_delta, n_groups, n_sides)
  if (n_sides == 2) { pp <- pp * 2 } # workaround for proportion_detectable function: for 2-sample, need to check both sides for beta, but not correct the pp (since already adjusted)
  if (pp == 0) return(0)
  proportion_detectable(pp, beta, delta, sigma_delta, n_groups, n_sides)
}


# ------------- ORGANIZING, FITTING, AND PLOTTING -------------------------


########### Summarize studies ########### 

get_study_summaries <- function(data, study, estimate, combo_name) {
  
  # preallocate
  d_name <- numeric(length(data))
  d_mean <- numeric(length(data))
  d_var_xv <- numeric(length(data))
  d_var_xv__emp <- numeric(length(data))
  d_cons_mean <- numeric(length(data))
  d_cons_var_xv <- numeric(length(data))
  d_var_xv__emp_cons <- numeric(length(data))
  d_n <- numeric(length(data))
  d_k <- numeric(length(data))
  d_mv <- numeric(length(data))
  d_mv_lb <- numeric(length(data))
  d_mv_ub <- numeric(length(data))
  vi <- numeric(length(data))
  vi_var_xv__emp <- numeric(length(data))
  vi_mv <- numeric(length(data))
  shapiro <- numeric(length(data))
  
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
    d_var_xv[i] <- var(d)
    d_var_xv__emp[i] <- (get_characteristic_magnitude_for_sd(d_mean[i], sqrt(d_var_xv[i])))^2
    
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
    
    # get mask - already pre-masked
    # if (combo_name %in% names(masks[[i]])) { # meta-analysis
    #   mask <- masks[[i]][[combo_name]]$mask
    # } else { # individual study
    #   mask <- masks[[i]]$mask
    # }
    # if (sum(mask) < length(data[[i]][[combo_name]][[estimate]])) {
    #   sorted_indices <- which(mask == 1)[order(data[[i]][[combo_name]][[estimate]][mask])]
    # } else {
    sorted_indices <- order(data[[i]][[combo_name]][[estimate]])
    # }
    
    # sort data from smallest to largest effect size
    sorted_estimate <- data[[i]][[combo_name]][[estimate]][sorted_indices]
    sorted_upper_bounds <- data[[i]][[combo_name]][[ci_ub]][sorted_indices]
    sorted_lower_bounds <- data[[i]][[combo_name]][[ci_lb]][sorted_indices]
    
    sorted_cons_estimate <- ifelse((abs(sorted_lower_bounds) > abs(sorted_upper_bounds)),
                                   ifelse((sorted_upper_bounds < 0),
                                          round(sorted_upper_bounds, 2), 0),
                                   ifelse((sorted_lower_bounds > 0),
                                          round(sorted_lower_bounds, 2), 0))
    
    
    d_cons_mean[i] <- mean(sorted_cons_estimate)
    d_cons_var_xv[i] <- var(sorted_cons_estimate)
    d_var_xv__emp_cons[i] <- (get_characteristic_magnitude_for_sd(d_cons_mean[i], sqrt(d_cons_var_xv[i])))^2
    
    # 3. Get multivariate effect size and bounds if available
    
    mv_combo_name <- gsub("mv.none", "mv.multi", combo_name) # get multivariate combo name
    all_combos <- names(data[[i]])
    mv_combo_idx <- grep(mv_combo_name, all_combos)
    # if doesn't exist, replace motion.xxx. with motion.threshold.
    # if (length(mv_combo_idx) == 0) {
    #   mv_combo_name_alt <- gsub("motion\\.[^\\.]+\\.", "motion.threshold.", mv_combo_name)
    #   mv_combo_idx <- grep(mv_combo_name_alt, all_combos)
    # }
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
          vi[i] <- d_se(d_var_xv[i], data[[i]][[combo_name]]$n1, data[[i]][[combo_name]]$n2)^2
          vi_var_xv__emp[i] <- d_se(d_var_xv__emp[i], data[[i]][[combo_name]]$n1, data[[i]][[combo_name]]$n2)^2
          vi_mv[i] <- d_se(d_mv[i], data[[i]][[combo_name]]$n1, data[[i]][[combo_name]]$n2)^2
        } else if (estimate == "r_sq") {
          vi[i] <- r_sq_se(d_var_xv[i], data[[i]][[combo_name]]$n1 + data[[i]][[combo_name]]$n2)^2
          vi_var_xv__emp[i] <- r_sq_se(d_var_xv__emp[i], data[[i]][[combo_name]]$n1 + data[[i]][[combo_name]]$n2)^2
          vi_mv[i] <- r_sq_se(d_mv[i], data[[i]][[combo_name]]$n1 + data[[i]][[combo_name]]$n2)^2
        }
      } else {
        # TODO: for the meta-analytic results, we should really just get the results from the confidence intervals,
        # since these come directly from the meta-analysis. For now, we don't use the variances previously calculated
        # from the meta-analysis for the meta-regression below, so we will ignore
        # (add helper alongside other se estimators below)
        d_n[i] <- data[[i]][[combo_name]]$n
        vi[i] <- NA
        vi_var_xv__emp[i] <- NA
        vi_mv[i] <- NA
      }
      d_k[i] <- 4
    } else if (study$orig_stat_type[[i]] == "t" || study$orig_stat_type[[i]] == "r") {
      # if (!is.null(data[[i]][[combo_name]]$n)) {
      d_n[i] <- data[[i]][[combo_name]]$n
      if (estimate == "d") {
        if (study$orig_stat_type[[i]] == "r") { # treat as 2-sample t-test
          vi[i] <- d_se(d_var_xv[i], d_n[i]/2, d_n[i]/2)^2
          vi_var_xv__emp[i] <- d_se(d_var_xv__emp[i], d_n[i]/2, d_n[i]/2)^2
          vi_mv[i] <- d_se(d_mv[i], d_n[i]/2, d_n[i]/2)^2
          d_k[i] <- 4
        } else { # normal 1-sample t-test
          vi[i] <- d_se(d_var_xv[i], d_n[i])^2
          vi_var_xv__emp[i] <- d_se(d_var_xv__emp[i], d_n[i])^2
          vi_mv[i] <- d_se(d_mv[i], d_n[i])^2
          d_k[i] <- 1
        }
      } else if (estimate == "r_sq") {
        vi[i] <- r_sq_se(d_var_xv[i], d_n[i])^2
        vi_var_xv__emp[i] <- r_sq_se(d_var_xv__emp[i], d_n[i])^2
        vi_mv[i] <- r_sq_se(d_mv[i], d_n[i])^2
        d_k[i] <- 4
      }
    } else {
      d_n[i] <- NA
    }
    
    # d_biased_sd[i] <- sd(d) * sqrt((length(d) - 1) / length(d)) # biased sd - doesn'd change results
    
    # test for normality (usually light- or heavy-tailed, some approximately normal)
    if (length(sorted_estimate) > 1) {
      k <- min(5000, length(sorted_estimate))  # max 5000 variables for shapiro test
      d_subset <- sorted_estimate[seq(1, length(sorted_estimate), length.out = k)]
      shapiro[i] <- shapiro.test(d_subset)$p.value # p-value < 2e-16 -> not normal
      # shapiro[i] <- ks.test(sorted_estimate,'pnorm',mean=mean(sorted_estimate),sd=sd(sorted_estimate))$p.value # p-value < 2e-16 -> not normal
      # qqnorm(d, pch=20, cex=0.5); qqline(d) # visualize
    } else {
      shapiro[i] <- NA # not enough data to test
    }
  }
  
  # set up final data frame
  
  df <- data.frame(name = names(data), mean = d_mean, var_xv = d_var_xv, var_xv__emp = d_var_xv__emp, mean_cons = d_cons_mean, var_xv_cons = d_cons_var_xv, var_xv__emp_cons = d_var_xv__emp_cons, n = d_n, category = as.factor(study$category), dataset = I(study$dataset), ref = study$ref, orig_stat_type = unlist(study$orig_stat_type), mv = d_mv, mv_lb = d_mv_lb, mv_ub = d_mv_ub, vi_var_xv = vi, vi_var_xv__emp = vi_var_xv__emp, vi_mv = vi_mv, k = d_k, shapiro = shapiro, stringsAsFactors = FALSE)
  
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

estimate_params <- function(df, df_meta, n_pts, main_title, fn, plot_type = "crossvariable", use_bayesian_fit = FALSE, plot_params = NULL) {
  
  print(paste0('Fitting lines for ', main_title))
  
  # params
  ndivk_max_extra_padding <- 1000 # so plot xlim extend a bit beyond max n
  use_var_xv__emp <- FALSE
  fit_base_width <- plot_params$estimate_fits$base_width
  fit_base_height <- plot_params$estimate_fits$base_height
  
  # Determine y variable and settings based on plot_type
  if (plot_type == "crossvariable") {
    y_var <- "var_xv"
    y_label <- expression("Observed cross-brain effect variance (var("*hat(theta)*"))")
    y_limits <- c(-0.021, 0.155)
  } else if (plot_type == "mv") { # note: all the same procedures here can also be used for the univariate case, just need to change the "mv" variable
    y_var <- "mv"
    y_label <- expression("Observed multivariate effect ("*hat(theta)*")")
    y_limits <- NULL  # Let ggplot auto-scale
  } else {
    stop("plot_type must be either 'crossvariable' or 'mv'")
  }
  
  # sort by sample size
  df <- df[order(df$n/df$k, decreasing = FALSE), ]
  
  # set ndivk_max_plt to the maximum observed sample size in df (fallback to large default if unavailable)
  ndivk_min_plt <- 10
  ndivk_max_plt <- max(df$n, na.rm = TRUE) + ndivk_max_extra_padding
  
  # interpolate predictions for all n (use plotted max n only)
  # NOTE: this used to be recomputed inside the add_meta loop below (identically,
  # every time) -- moved up here since it doesn't depend on add_meta.
  ndivk__seq <- seq(ndivk_min_plt, ndivk_max_plt, length.out = n_pts)
  
  # set up for each overarching category
  unique_cats <- unique(df$overarching_category)
  # # if 'task activation' exists in unique_cats, make it come right after 'task connectivity'
  # if ("task activation" %in% unique_cats) {
  #   unique_cats <- c(setdiff(unique_cats, "task connectivity"), "task connectivity")
  #   unique_cats <- c(setdiff(unique_cats, "task activation"), "task activation")
  # }
  
  # fit meta-analysis nesting studies by dataset and overarching category
  
  # setup grouping variables (used by the frequentist rma.mv nesting)
  df$dataset_nested <- interaction(df$overarching_category, df$dataset, drop = TRUE)
  
  if (plot_type == "mv") {
    if (use_bayesian_fit) {
      # NOTE: bmr() has no equivalent of the category/dataset_nested nesting (pools heterogeneity into a single tau).
      # Category is entered as a fixed effect (one-hot columns in X), and each coefficients is obtained from posterior summary.
      
      keep <- !is.na(df$mv) & !is.na(df$vi_mv)
      df_keep <- droplevels(df[keep, ])
      X_bayes <- model.matrix(~ 0 + overarching_category, data = df_keep)
      colnames(X_bayes) <- levels(df_keep$overarching_category)
      print("  Fitting bmr() [mv, intercept-only]...")
      t0 <- Sys.time()
      fit_all <- bmr(y = df$mv[keep],
                     sigma = sqrt(df$vi_mv[keep]),
                     X = X_bayes,
                     labels = df$name[keep],
                     tau.prior = "uniform")
      print(paste0("    ...done in ", round(difftime(Sys.time(), t0, units = "secs"), 1), " sec"))
    } else {
      fit_all <- rma.mv(yi = mv, 
                        V = vi_mv,  # approximate variance from CI
                        random = ~ 1 | overarching_category/dataset_nested,
                        data = df,
                        method = "REML")
    }
    
  } else {
    # NOTE: y_var/v_var below used to be overwritten with the actual numeric
    # data vectors here, clobbering the column-name STRING that y_var holds
    # above (used in aes_string() for plotting). That "worked" by accident
    # for the main df layer (same row count as the data), but broke the
    # df_meta overlay layer (different row count) with a check_aesthetics
    # length-mismatch error. Using separate y_vals/v_vals for the fit data
    # keeps y_var as the plotting column name throughout.
    if (use_var_xv__emp) {
      y_vals <- df$var_xv__emp
      v_vals <- df$vi_var_xv__emp
    } else {
      y_vals <- df$var_xv
      v_vals <- df$vi_var_xv
    }
    
    if (use_bayesian_fit) {
      keep <- !is.na(y_vals) & !is.na(v_vals) & !is.na(df$k) & !is.na(df$n)
      df_keep <- droplevels(df[keep, ])
      X_cat <- model.matrix(~ 0 + overarching_category, data = df_keep)
      colnames(X_cat) <- levels(df_keep$overarching_category)
      X_bayes <- cbind(X_cat, "invn" = (df$k / df$n)[keep])  # slope shared across categories
      print("  Fitting bmr() [crossvariable]...")
      t0 <- Sys.time()
      fit_all <- bmr(y = y_vals[keep], sigma = sqrt(v_vals[keep]),
                     X = X_bayes,
                     labels = df$name[keep],
                     tau.prior = "uniform")
      print(paste0("    ...done in ", round(difftime(Sys.time(), t0, units = "secs"), 1), " sec"))
    } else {
      if (use_var_xv__emp) {
        fit_all <- rma.mv(yi = var_xv__emp, 
                          V = vi_var_xv__emp,
                          mods = ~ I(k/n),
                          random = ~ 1 | overarching_category/dataset_nested,
                          data = df,
                          method = "REML")
      } else {
        fit_all <- rma.mv(yi = var_xv, 
                          V = vi_var_xv,
                          mods = ~ I(k/n),
                          random = ~ 1 | overarching_category/dataset_nested,
                          data = df,
                          method = "REML")
      }
    }
  }
  
  # ------------------------------------------------------------------
  # Extract estimates and credible/confidence intervals (also only once).
  # ------------------------------------------------------------------
  res <- vector("list", length(unique_cats))
  names(res) <- unique_cats
  predicted_y <- vector("list", length(unique_cats))
  names(predicted_y) <- unique_cats
  
  if (use_bayesian_fit) {
    
    # bmr(): category is a fixed effect in X (one column per category), so
    # each category's estimate is read straight off the posterior summary.
    # NOTE: confirm these row labels match your installed bayesmeta/bmr
    # version once via: rownames(fit_all$summary); colnames(fit_all$summary)
    post_summary <- fit_all$summary
    est_row <- "mean"
    lwr_row <- "95% lower"
    upr_row <- "95% upper"
    
    # bmr() renames names with spaces internally (somewhere between X_bayes and fit_all$summary)
    # so below we matches exactly first, then fall back to a make.names()-normalized match to be 
    # robust to convention is in play (can confirm with colnames(fit_all$summary); colnames(X_bayes), or X_cat)
    resolve_col <- function(name, choices) {
      if (name %in% choices) return(name)
      hit <- choices[make.names(choices) == make.names(name)]
      if (length(hit) == 1) return(hit)
      return(NA_character_)
    }
    
    if (plot_type == "mv") {
      
      for (cat in unique_cats) {
        cat_col <- resolve_col(as.character(cat), colnames(post_summary))
        
        # Leave res[[cat]]/predicted_y[[cat]] as NULL rather than NA to silently drop those nulls during do.call(rbind, res) rather than tripping over NAs later
        # more detail: unique_cats is derived from full df, but X_bayes (and thus post_summary's columns) only includes categories that had at least one non-missing row for this estimate/plot_type.
        if (is.na(cat_col)) {
          warning(paste0("Category '", as.character(cat), "' has no usable data for this fit (plot_type = '",
                         plot_type, "') -- excluding from results."))
          next
        }
        
        cat_est <- post_summary[est_row, cat_col]
        cat_lwr <- post_summary[lwr_row, cat_col]
        cat_upr <- post_summary[upr_row, cat_col]
        
        res[[cat]] <- data.frame(
          est = cat_est,
          lwr = cat_lwr,
          upr = cat_upr,
          row.names = paste0(cat, "_intercept")
        )
        
        # Constant line for mv (intercept-only model)
        predicted_y[[cat]] <- cbind(
          fit = rep(cat_est, length(ndivk__seq)),
          lwr = rep(cat_lwr, length(ndivk__seq)),
          upr = rep(cat_upr, length(ndivk__seq))
        )
      }
      
    } else {
      
      slope_est <- post_summary[est_row, "invn"]
      slope_lwr <- post_summary[lwr_row, "invn"]
      slope_upr <- post_summary[upr_row, "invn"]
      
      # BEWARE: rposterior() is slow -- samples via numerical inversion, and by default (tau.sample=TRUE) also draws tau,
      # requiring root-finding per draw. TODO: consider setting tau.sample=FALSE
      n_draws <- 500
      print(paste0("  Drawing ", n_draws, " posterior samples via rposterior()..."))
      t0 <- Sys.time()
      draws <- fit_all$rposterior(n_draws)
      print(paste0("    ...done in ", round(difftime(Sys.time(), t0, units = "secs"), 1), " sec"))
      # NOTE: inspect once with str(draws) / colnames(draws) to confirm
      # it returns a draw per X column (+ tau) under these same names.
      
      for (cat in unique_cats) {
        cat_col <- resolve_col(as.character(cat), colnames(post_summary))
        
        # as above, set to NULL rather than NA to catch early
        if (is.na(cat_col)) {
          warning(paste0("Category '", as.character(cat), "' has no usable data for this fit (plot_type = '",
                         plot_type, "') -- excluding from results."))
          next
        }
        
        # draws may use the same (possibly renamed) col names as post_summary, so
        # resolve against draws' own colnames rather than assuming match
        draws_col <- resolve_col(as.character(cat), colnames(draws))
        if (is.na(draws_col)) {
          warning(paste0("Category '", as.character(cat), "' found in post_summary but not in ",
                         "rposterior() draws -- excluding from results. Check colnames(draws) ",
                         "against colnames(fit_all$summary)."))
          next
        }
        
        cat_est <- post_summary[est_row, cat_col]
        cat_lwr <- post_summary[lwr_row, cat_col]
        cat_upr <- post_summary[upr_row, cat_col]
        
        res[[cat]] <- data.frame(
          est = cat_est,
          lwr = cat_lwr,
          upr = cat_upr,
          phi2_est = slope_est,
          phi2_lwr = slope_lwr,
          phi2_upr = slope_upr,
          row.names = paste0(cat, "_intercept")
        )
        
        pred_draws <- outer(draws[, draws_col], rep(1, length(ndivk__seq))) +
          outer(draws[, "invn"], 1 / ndivk__seq)
        predicted_y[[cat]] <- cbind(
          fit = colMeans(pred_draws),
          lwr = apply(pred_draws, 2, quantile, probs = 0.025),
          upr = apply(pred_draws, 2, quantile, probs = 0.975)
        )
      }
    }
    
  } else {
    
    # Frequentist (rma.mv) path -- category-specific intercepts via BLUPs,
    # kept as an alternative/sensitivity estimation procedure.
    random_effects <- ranef(fit_all)
    category_effects <- random_effects$overarching_category
    
    for (cat in unique_cats) {
      # Get category-specific intercept (overall + random effect)
      if (cat %in% rownames(category_effects)) {
        cat_intercept <- fit_all$beta[1] + category_effects[cat, "intrcpt"]
        cat_intercept_se <- sqrt(fit_all$vb[1,1] + category_effects[cat, "se"]^2)
      } else {
        # If category not found, use overall intercept
        cat_intercept <- NA
        cat_intercept_se <- NA
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
        preds <- rep(cat_intercept, length(ndivk__seq))
        preds_se <- rep(cat_intercept_se, length(ndivk__seq))
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
          phi2_est = slope,
          phi2_lwr = slope - 1.96 * slope_se,
          phi2_upr = slope + 1.96 * slope_se,
          row.names = paste0(cat, "_intercept")
        )
        
        # Create predicted values with category-specific intercept
        preds <- cat_intercept + slope * 1/ndivk__seq
        # Standard errors for predictions (more complex with random effects)
        X_pred <- cbind(1, 1/ndivk__seq)
        preds_se_fixed <- sqrt(diag(X_pred %*% fit_all$vb %*% t(X_pred)))
        preds_se <- sqrt(preds_se_fixed^2 + cat_intercept_se^2)
        lwr <- preds - 1.96 * preds_se
        upr <- preds + 1.96 * preds_se
        predicted_y[[cat]] <- cbind(fit = preds, lwr = lwr, upr = upr)
      }
    }
  }
  
  # Combine results into single dataframe
  res <- do.call(rbind, res)
  
  # ------------------------------------------------------------------
  # Plot, optionally adding meta points/labels
  # ------------------------------------------------------------------
  
  df$x_plot <- df$n/df$k
  x_label <- "n/k (log scale)"
  # use ndivk_max_plt (prediction max) to set the upper x limit
  x_limits <- c(ndivk_min_plt, ndivk_max_plt)
  
  cats <- levels(df$overarching_category)
  color_map <- setNames(RColorBrewer::brewer.pal(length(cats), "Set1"), cats)
  
  for (add_meta in c(TRUE, FALSE)) {
    
    if (add_meta) {
      meta_str <- '__meta'
      alpha <- 0.15 # make non-meta points more transparent
    } else {
      meta_str <- ''
      alpha <- 0.7
    }
    
    fitlines_str <- '__fits'
    
    p <- ggplot(df, aes_string(x = "x_plot", y = y_var, color = "overarching_category")) +
      geom_point(size = 1.5, alpha = alpha, stroke = 0) +
      scale_color_manual(values = color_map) +
      labs(title = main_title,
           x = x_label,
           y = y_label,
           color = "Category") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5,face = "bold"),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.key.size = unit(0.7, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(expand = c(0, 0), trans = "log") +
      coord_cartesian(xlim = x_limits) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3)
    
    if (!is.null(y_limits)) {
      # p <- p + scale_y_continuous(limits = y_limits, expand = c(0, 0))
      p <- p + coord_cartesian(ylim = y_limits)
    }
    
    for (cat in unique_cats) {
      pred_mat <- predicted_y[[cat]]
      if (!is.null(pred_mat) && is.matrix(pred_mat) && all(c("fit","lwr","upr") %in% colnames(pred_mat)) && nrow(pred_mat) == length(ndivk__seq)) {
        pred_df <- data.frame(
          x_plot = ndivk__seq,
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
      df_meta$x_plot <- df_meta$n/df_meta$k
      p <- p +
        geom_point(data = df_meta, aes_string(x = "x_plot", y = y_var, color = "overarching_category"), shape = 8, size = 2) +
        geom_text_repel(data = df_meta, aes_string(x = "x_plot", y = y_var, label = "label"), size = 2.2, segment.size = 0.3, force = 10, max.overlaps = Inf)
    }
    
    fn_plot <- fn
    this_fn <- paste0(fn_plot, "__fits", meta_str, '.pdf')
    if (save_plots) {
      ggsave(this_fn, plot = p, width = fit_base_width, height = fit_base_height)
    } else {
      show(p)
    }
    
    # Save normality test results
    if (save_plots  && plot_type != "mv") {
      write.csv(paste0('Proportion significantly non-normal: ', sum(df$shapiro < 0.05)/length(df$shapiro),' (',sum(df$shapiro < 0.05),' studies)'),file=paste0(fn, '_shapiro_proportion_sig.csv'))
      write.csv(df$shapiro, file=paste0(fn, '_shapiro.csv'), row.names=TRUE)
    }
  }
  
  # save res Rdata
  this_fn <- paste0(fn, '_res.Rdata')
  save(res, file = this_fn)
  
  return(res)
}


###########  Make Estimated Density Plots ########### 

plot_densities <- function(res, res_mv,  n_pts, fn_basedir, cat_colors, save_plots = TRUE, plot_params = NULL) {
  
  print('Making density plots')
  cats <- unique(rownames(res))
  # cat_colors <- RColorBrewer::brewer.pal(length(cats), "Set1")
  
  # params
  density_cfg <- plot_params$density
  do_horizontal_panels <- density_cfg$do_horizontal_panels
  xlim_annotate <- density_cfg$xlim_annotate
  xbreaks_base <- density_cfg$xbreaks
  axis_text_size <- density_cfg$axis_text_size
  axis_title_size <- density_cfg$axis_title_size
  base_width <- density_cfg$base_width
  base_height <- density_cfg$base_height
  
  for (do_mv in c(FALSE, TRUE)) {
    if (do_mv) {
      mv_suffix <- '_mv'
      xlim <- c(0, 5)
      xlim_annotate[1] <- 0
      xtick_angle <- 45 # fit stuff if mv
      xbreaks <- sort(unique(xbreaks_base[xbreaks_base >= xlim[1]]))
      xbreaks <- xbreaks[abs(xbreaks - 0.2) > 1e-9]
      vjust <- 1
      hjust <- 1
    } else {
      mv_suffix <- ''
      xlim <- c(-0.8, 0.8)
      xtick_angle <- 45
      xbreaks <- xbreaks_base
      vjust <- 1
      hjust <- 1
    }
    d <- seq(xlim[1], xlim[2], length.out = n_pts)
    width <- base_width
    height <- base_height
    
    
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
          if (!is.na(mu) && mu < 0) {
            mu <- 0
          }
          
          y <- rep(0, length(d))
          closest_idx <- which.min(abs(d - mu))
          y[closest_idx] <- 1
          
          density_list[[length(density_list)+1]] <- data.frame(
            d = d,
            density = y,
            category = cat,
            sigma_type = name,
            overarching_category = cat
          )
        }
        
      } else {
        # vars <- c(res[cat,"lwr"], res[cat,"est"], res[cat,"upr"])
        # vars <- pmax(vars, 0) # no negative variances
        # sigmas <- sqrt(vars)
        # # print(sigmas)
        # names(sigmas) <- c('lwr','est','upr')
        # sigmas_master <- rbind(sigmas_master, sigmas)
        
        res[res < 0] <- 0 # no negative variances
        sigmas_master <- sqrt(res)
        sigmas <- sigmas_master[cat, ]
        
        
        # If not mv, make y as dnorm for each sigma (lwr, est, upr)
        for (name in names(sigmas)) {
          
          s <- sigmas[[name]]
          
          if (s > 0) {
            y <- dnorm(d, mean = 0, sd = s)
          } else {
            # point mass at 0 with exponential taper toward the midpoint so it plots nicely
            # note that these are arbitrary values set so point mass shows up for plots with ymax=~3 - ymax=~8
            y <- rep(0, length(d))
            midpt <- ceiling(length(d)/2)
            max_density <- 30
            n_taper <- 70
            weights <- exp(seq(log(0.1), log(1), length.out = n_taper))
            y[(midpt - n_taper/2 + 1):(midpt + n_taper/2)] <- max_density * weights
            # weights <- exp(seq(log(0.1), log(1), length.out = n_taper))
            # y[(midpt - n_taper + 1):midpt] <- max_density * weights
            # y[(midpt + 1):(midpt + n_taper)] <- max_density * rev(weights)
          }
          
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
    
    # if (!do_mv) {
    # sigmas_master <- as.data.frame(sigmas_master)
    # rownames(sigmas_master) <- cats
    # TODO: can probably just use res, instead of recreating and renaming sigmas_master
    # }
    
    
    density_df <- do.call(rbind, density_list)
    
    y_max <- max(density_df$density)
    
    # Overlapping densities colored by overarching category, legend inset
    p_density <- ggplot(density_df, aes(x = d, y = density, color = overarching_category, linetype = sigma_type)) +
      geom_line(size = 1) +
      scale_color_manual(values = cat_colors) +
      labs(title = "Density Plot by Overarching Category", x = "Cohen's d", y = "Density", color = "Category", linetype = "Sigma Type") +
      theme_classic() +
      coord_cartesian(xlim = xlim, ylim = c(0, y_max)) +
      theme(
        axis.text.x = element_text(angle = xtick_angle, vjust = vjust, hjust = hjust, size = axis_text_size),
        legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.key.size = unit(0.7, "lines")
      )
    if (save_plots) {
      ggsave(paste0(fn_basedir, 'density',mv_suffix,'.pdf'), p_density, width = width, height = height)
    } else {
      print(p_density)
    }
    
    # Panel densities
    
    if (do_horizontal_panels) {
      nrow <- 1
      panel_width <- width * length(cats) #4 * length(cats)
      panel_height <- height #3.5 #3.8
    } else {
      nrow <- length(cats)
      panel_width <- width
      panel_height <- height * length(cats) #3.8 * length(cats)
    }
    density_df$fill <- density_df$sigma_type=="est"
    
    # Get unique categories and create individual plots
    unique_cats <- unique(density_df$category)
    plot_list <- list()
    
    for (i in seq_along(unique_cats)) {
      cat <- unique_cats[i]
      p <- ggplot(density_df %>% filter(category == cat), aes(x = d, y = density, color = overarching_category, linetype = sigma_type)) +
        geom_ribbon(data = subset(density_df, category == cat & sigma_type == "est"),
                    aes(ymin = 0, ymax = density, fill = overarching_category), 
                    alpha = 0.8, colour = NA, show.legend = FALSE) +
        annotate("rect", xmin = xlim_annotate[1], xmax = xlim_annotate[2], ymin = 0, ymax = Inf, 
                 fill = "gold", colour=NA, alpha = 0.6) + # using annotate to avoid drawing multiple->too high opacity
        geom_line(linewidth = 0.8, lineend = "butt") +
        scale_color_manual(values = cat_colors) +
        # use named linetypes so dash patterns remain visible at thicker linewidth
        scale_linetype_manual(values = c(est = "solid", lwr = "dashed", upr = "dotdash")) +
        scale_fill_manual(values = cat_colors, guide = "none") +
        scale_y_continuous(limits = c(0, max(density_df %>% filter(category == cat, sigma_type == "est") %>% pull(density), na.rm = TRUE))) +
        scale_x_continuous(breaks = xbreaks) +
        # labs(title = cat) +
        labs(title = cat, x = "Cohen's d", y = "Density") +
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = xtick_angle, vjust = vjust, hjust = hjust, size = axis_text_size),
              axis.text.y = element_text(size = axis_text_size),
              # axis.title.x = element_text(size = axis_title_size),
              # axis.title.y = element_text(size = axis_title_size)
              axis.title.x = element_blank(),
              axis.title.y = element_blank()
        )
      
      plot_list[[cat]] <- p
    }
    
    # Combine plots with patchwork
    p_density_panel <- Reduce(`+`, plot_list) +
      plot_layout(nrow = nrow, guides = 'collect') &
      theme(legend.position = c(0.02, 0.98), 
            legend.justification = c("left", "top"),
            legend.background = element_rect(fill = "white", color = "grey80"))
    
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
      ggsave(paste0(fn_basedir, 'density__panels',mv_suffix,'.pdf'), p_density_panel, width = panel_width, height = panel_height)
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

##### POWER PLOTS #####

# make average power vs. sample size plots

get_average_power <- function(sigmas_master, res_mv = NULL, do_mv = FALSE) {
  
  if (do_mv) {
    mv_suffix <- '_mv'
    xlim <- c(0, 6)
  } else {
    mv_suffix <- ''
    xlim <- c(-0.8, 0.8)
  }
  
  alpha <- 0.05
  n_sides <- 2 # two-tailed test
  target_power <- 0.8
  
  n_vector <- c(0, 25, 50, 100, 500, 1000, 5000, 50000, 300000, Inf)
  avg_power <- data.frame()
  proportion_detectable <- data.frame()
  
  # get cats from res
  cats <- rownames(sigmas_master)
  for (cat in cats) {
    
    # run power
    
    # set up temporary data frame for this cat
    avg_power_tmp <- data.frame(n = n_vector, overarching_category = cat)
    proportion_detectable_tmp <- data.frame(n = n_vector, overarching_category = cat)
    
    # set test type
    test_type <- if (grepl("task", cat)) "one.sample" else "two.sample"
    ifelse(test_type == "one.sample", n_groups <- 1, n_groups <- 2)
    
    if (do_mv) {
      avg_power_tmp$uncorrected <- sapply(n_vector, function(n) pwr.t.test(n = n, d = res_mv[cat, "est"], sig.level = alpha, type = test_type, alternative = "greater")$power)
      avg_power_tmp$uncorrected[n_vector==0] <- 0 # power=0 at n=0
      avg_power_tmp$bonferroni <- NA
      avg_power_tmp$fdr <- NA
      
      proportion_detectable_tmp$uncorrected <- as.numeric(avg_power_tmp$uncorrected > target_power)
      proportion_detectable_tmp$bonferroni <- NA
      proportion_detectable_tmp$fdr <- NA
      
    } else {
      
      sigmas <- sigmas_master[cat, ]
      this_sigma <- as.numeric(sigmas["est"])
      
      # get average power at each n
      
      # comparisons based on studies available in meta-analysis 02102026
      if (cat=="task activation") {
        k <- 204899 # based on number of common voxels in meta-analysis
      } else {
        k <- 35778 # Shen atlas dimensionality
      }
      warning(paste0("Using hard-coded number of tests (k=", k, ") for power calculations."))
      
      
      avg_power_tmp$uncorrected <- sapply(n_vector, function(n) F1(alpha, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
      avg_power_tmp$bonferroni <- sapply(n_vector, function(n) F1(alpha/k, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
      avg_power_tmp$fdr <- sapply(n_vector, function(n) BHpower(0, alpha, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
      
      proportion_detectable_tmp$uncorrected <- sapply(n_vector, function(n) proportion_detectable(alpha, 1-target_power, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
      proportion_detectable_tmp$bonferroni <- sapply(n_vector, function(n) proportion_detectable(alpha/k, 1-target_power, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
      proportion_detectable_tmp$fdr <- sapply(n_vector, function(n) BH_proportion_detectable(0, alpha, 1-target_power, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
    }
    
    # make long, moving correction type to a new column
    avg_power_tmp <- avg_power_tmp %>%
      pivot_longer(
        cols = c(uncorrected, bonferroni, fdr),
        names_to = "correction_type",
        values_to = "avg_power"
      )
    
    proportion_detectable_tmp <- proportion_detectable_tmp %>%
      pivot_longer(
        cols = c(uncorrected, bonferroni, fdr),
        names_to = "correction_type",
        values_to = "proportion_detectable"
      )
    
    avg_power <- rbind(avg_power, avg_power_tmp)
    proportion_detectable <- rbind(proportion_detectable, proportion_detectable_tmp)
    
  }
  
  # Preserve facet order to match cats vector
  avg_power$overarching_category <- factor(avg_power$overarching_category, levels = cats)
  proportion_detectable$overarching_category <- factor(proportion_detectable$overarching_category, levels = cats)
  return(list(avg_power=avg_power, proportion_detectable=proportion_detectable))
  
}


# plot avg power

plot_average_power <- function(df, do_mv = FALSE, cat_colors, fn_basedir, save_plots = TRUE, plot_params = NULL) {
  
  print("Making power plots")
  power_cfg <- plot_params$power
  do_horizontal_panels <- power_cfg$do_horizontal_panels
  
  title <- if (do_mv) "Average Power by Category (Multivariate)" else "Average Power by Category"
  filename <- if (do_mv) 'power_panels_mv.pdf' else 'power_panels.pdf'
  
  if (do_horizontal_panels) {
    nrow <- 1
    width <- power_cfg$base_width * length(cat_colors)
    height <- power_cfg$base_height
  } else {
    nrow <- length(cat_colors)
    width <- power_cfg$single_col_width
    height <- power_cfg$base_height * length(cat_colors)
  }
  
  n_str <- as.character(sort(unique(df$n)))
  n_str[length(n_str)] <- "∞"
  x100_idx <- match(100, sort(unique(df$n)))
  
  p <- ggplot(df, aes(x = factor(n), y = avg_power, group = interaction(overarching_category,correction_type), color = overarching_category, linetype = correction_type)) +
    geom_line(size = 1.5) +
    scale_x_discrete(labels = n_str) +
    scale_color_manual(values = cat_colors) +
    scale_linetype_manual(values = c("uncorrected" = "solid", "fdr"  = "dashed", "bonferroni" = "dotted")) +
    facet_wrap(~overarching_category, nrow = nrow, scales = "free_y") +
    labs(title = title, x = "Sqrt Sample Size", y = "Average Power", color = "Category", linetype = "Correction Type") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none"
    )
  
  if (!is.na(x100_idx)) {
    p <- p + geom_vline(xintercept = x100_idx, linetype = "dotted", color = "grey60", linewidth = 0.8)
  }
  
  if (save_plots) {
    ggsave(paste0(fn_basedir, filename), p, width = width, height = height)
  } else {
    print(p)
  }
}

# plot proportion detectable
plot_proportion_detectable <- function(df, do_mv = FALSE, cat_colors, fn_basedir, save_plots = TRUE, plot_params = NULL) {
  
  power_cfg <- plot_params$power
  do_horizontal_panels <- power_cfg$do_horizontal_panels
  
  title <- if (do_mv) "Proportion Detectable by Category (Multivariate)" else "Proportion Detectable by Category"
  filename <- if (do_mv) 'proportion_detect_mv.pdf' else 'proportion_detect.pdf'
  
  if (do_horizontal_panels) {
    nrow <- 1
    width <- power_cfg$base_width * length(cat_colors)
    height <- power_cfg$base_height
  } else {
    nrow <- length(cat_colors)
    width <- power_cfg$single_col_width
    height <- power_cfg$base_height * length(cat_colors)
  }
  axis_text_size <- power_cfg$axis_text_size
  axis_title_size <- power_cfg$axis_title_size
  
  n_str <- as.character(sort(unique(df$n)))
  n_str[length(n_str)] <- "∞"
  x100_idx <- match(100, sort(unique(df$n)))
  
  p <- ggplot(df, aes(x = factor(n), y = proportion_detectable, group = interaction(overarching_category,correction_type), color = overarching_category, linetype = correction_type)) +
    geom_line(size = 1.5) +
    scale_x_discrete(labels = n_str) +
    scale_color_manual(values = cat_colors) +
    scale_linetype_manual(values = c("uncorrected" = "solid", "fdr"  = "dashed", "bonferroni" = "dotted")) +
    facet_wrap(~overarching_category, nrow = nrow, scales = "free_y") +
    labs(title = title, x = "Sample Size", y = "Proportion Detectable", color = "Category", linetype = "Correction Type") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      # axis.title.x = element_text(size = axis_title_size),
      # axis.title.y = element_text(size = axis_title_size),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none"
    )
  
  if (!is.na(x100_idx)) {
    p <- p + geom_vline(xintercept = x100_idx, linetype = "dotted", color = "grey60", linewidth = 0.8)
  }
  
  if (save_plots) {
    ggsave(paste0(fn_basedir, filename), p, width = width, height = height)
  } else {
    print(p)
  }
}


# proportion detected under different estimation conditions
plot_proportion_difference <- function(sigmas_master, phi2, cat_colors, fn_basedir) {
  
  alpha <- 0.05
  n_sides <- 2 # two-tailed test
  target_power <- 0.8
  
  n_vector <- c(0, 25, 50, 100, 500, 1000, 5000)
  
  # get categories
  cats <- rownames(sigmas_master)
  
  # preallocate
  proportion_detectable <- data.frame()
  diff_detections <- data.frame()
  
  for (cat in cats) {
    
    # # set number of groups
    if (grepl("task", cat)) {
      n_groups <- 1
    } else {
      n_groups <- 2
    }
    
    # this_sigma_uncorrected <- sqrt(sigmas_master[cat, "est"]^2 + n_groups^2 * phi2 / n)
    # this_sigma <- sigmas_master[cat, "est"]
    
    #preallocate
    proportion_detectable_tmp <- data.frame(n = n_vector, overarching_category = cat)
    proportion_detectable_tmp__uncorrected <- data.frame(n = n_vector, overarching_category = cat)
    diff_detections_tmp <- data.frame(n = n_vector, overarching_category = cat)
    
    # Scenario 1: if you plan for an average effect
    proportion_detectable_tmp$uncorrected <- sapply(n_vector, function(n) proportion_detectable(alpha, 1-target_power, 0, sigmas_master[cat, "est"]*sqrt(n/n_groups^2),n_groups,n_sides))
    proportion_detectable_tmp__uncorrected$uncorrected <- sapply(n_vector, function(n) proportion_detectable(alpha, 1-target_power, 0, sqrt(sigmas_master[cat, "est"]^2 + n_groups^2 * phi2[cat, "phi2_est"] / n)*sqrt(n/n_groups^2),n_groups,n_sides))
    diff_detections_tmp$uncorrected <- proportion_detectable_tmp$uncorrected - proportion_detectable_tmp__uncorrected$uncorrected
    
    # proportion_detectable_tmp$bonferroni <- sapply(n_vector, function(n) proportion_detectable(alpha/k, 1-target_power, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
    # proportion_detectable_tmp$fdr <- sapply(n_vector, function(n) BH_proportion_detectable(0, alpha, 1-target_power, 0, this_sigma*sqrt(n/n_groups^2),n_groups,n_sides))
    
    
    
    # Scenario 2: if you plan for the strongest effect size (take top 10%)
    # prop_detect_est_strong <- proportion_detectable(alpha = 0.05, power = 0.8, mu = 0, sigma = sigma_uncorrected*sqrt(n/n_groups), n_groups = n_groups, n_sides = 2)
    # prop_detect_actual_strong <- proportion_detectable(alpha = 0.05, power = 0.8, mu = 0, sigma = sigma_actual*sqrt(n/n_groups), n_groups = n_groups, n_sides = 2)
    
    diff_detections_tmp <- diff_detections_tmp %>%
      pivot_longer(
        cols = c(uncorrected),
        names_to = "correction_type",
        values_to = "diff_proportion_detectable"
      )
    
    diff_detections <- rbind(diff_detections, diff_detections_tmp)
    
  }
  
  # TODO: plot both
  
  return(diff_detections)
}