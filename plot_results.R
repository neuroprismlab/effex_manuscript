# v_data is the braineffex data downloaded from OSF: a list called v
# OSF link: https://doi.org/10.17605/OSF.IO/CWNJD
# output_path is the directory to save results (should end with /)

#' Create results manuscript plots
#'
#' This function is the master plotter for making effect size plots for the results manuscript.
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @import lme4
#' @import pwr
#' @param estimate
#' @ param output_path Directory to save results
#' @param v_data v data loaded from braineffex_data file that can be downloaded from OSF
#' @param pooling_methods list of pooling methods to consider: 'none', 'net', or both
#' @param motion_method_t motion correction method to use for studies with one sample t-tests (default: "threshold")
#' @param motion_method_not_t motion correction method to use for studies without one sample t-tests (default: "regression")
#' @param save_plots Whether to save plots to files (default: TRUE)
#'
#' @return Saves plots to specified output directory
#'
#' @examples
#' # Example usage
#' \dontrun{
#' plot_results(estimate = 'd', output_path = 'results/', v_data = v, pooling_methods = c('none', 'net'))
#' }

plot_results <- function(estimate = 'd', output_path, v_data, pooling_methods = c('none', 'net'), motion_method_t = "threshold", motion_method_not_t = "regression", save_plots = TRUE) {

# libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
# library(ggridges)
library(lme4)
library(pwr) # for sample size calculations

  
  ###### HELPER FUNCTIONS ############

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

get_study_summaries <- function(data, study, motion_method_t, motion_method_not_t, meta, estimate) {
  
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
  
  for (i in seq_along(data)) {
    # find test type for this study
    this_test_type = study[i, "orig_stat_type"]
    this_motion = ifelse(this_test_type == "t", motion_method_t, motion_method_not_t) #TODO make dynamic with fn inputs
    
    if (meta) {
      # try regression motion first, if null then switch to threshold
      this_motion = "regression"
    }
    this_test_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.none")
    this_test_mv_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.multi")

    # check if estimate for the given study and combo is NULL
    if (anyNA(data[[i]][[this_test_combo]][[estimate]]) | anyNA(data[[i]][[this_test_mv_combo]][[estimate]])) {
      this_motion = "threshold"
      this_test_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.none")
      this_test_mv_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.multi")
    }

    # if estimate is still NULL, throw an error
    if (anyNA(data[[i]][[this_test_combo]][[estimate]])) {
      stop(paste0("Estimate '", estimate, "' is NULL for study ", i, ": ", study$name[i], " with test type ", this_test_type, " even after switching to threshold motion. Please check the data."))
    }

    print(paste0("using motion method ", this_motion, " for study ", i, ": ", study$name[i]))

    # 1. Get point estimate mean and sd across vars
    
    # to get d, check whether dim is null (nested list)
    if (is.null(dim(data[[i]][[this_test_combo]][[estimate]]))) {
      d <- data[[i]][[this_test_combo]][[estimate]]
    } else {
      d <- data[[i]][[this_test_combo]][[estimate]][1,]
    }
    
    d_mean[i] <- mean(d)
    d_sd[i] <- sd(d)
    d_characteristic_mag[i] <- get_characteristic_magnitude(d_mean[i], d_sd[i])
    
    # 2. Get conservative estimate mean and sd across vars
    # from prep_data_for_plot.R
    # omitting na check in mean and sd and downsampling
    
    # unlist sim CIs if list
    if (is.list(data[[i]][[this_test_combo]][[ci_lb]])) {
      data[[i]][[this_test_combo]][[ci_lb]] <- unlist(data[[i]][[this_test_combo]][[ci_lb]])
      data[[i]][[this_test_combo]][[ci_ub]] <- unlist(data[[i]][[this_test_combo]][[ci_ub]])
    }
    
    na_idx <- is.na(data[[i]][[this_test_combo]][[estimate]]) | is.na(data[[i]][[this_test_combo]][[ci_lb]]) | is.na(data[[i]][[this_test_combo]][[ci_ub]])
    data[[i]][[this_test_combo]][[estimate]] <- data[[i]][[this_test_combo]][[estimate]][!na_idx]
    data[[i]][[this_test_combo]][[ci_lb]] <- data[[i]][[this_test_combo]][[ci_lb]][!na_idx]
    data[[i]][[this_test_combo]][[ci_ub]] <- data[[i]][[this_test_combo]][[ci_ub]][!na_idx]
    
    # sort data from smallest to largest effect size
    sorted_indices <- order(data[[i]][[this_test_combo]][[estimate]])
    sorted_estimate <- data[[i]][[this_test_combo]][[estimate]][sorted_indices]
    sorted_upper_bounds <- data[[i]][[this_test_combo]][[ci_ub]][sorted_indices]
    sorted_lower_bounds <- data[[i]][[this_test_combo]][[ci_lb]][sorted_indices]
    
    sorted_cons_estimate <- ifelse((abs(sorted_lower_bounds) > abs(sorted_upper_bounds)),
                                   ifelse((sorted_upper_bounds < 0),
                                          round(sorted_upper_bounds, 2), 0),
                                   ifelse((sorted_lower_bounds > 0),
                                          round(sorted_lower_bounds, 2), 0))
    
    
    d_cons_mean[i] <- mean(sorted_cons_estimate)
    d_cons_sd[i] <- sd(sorted_cons_estimate)
    d_characteristic_mag_cons[i] <- get_characteristic_magnitude(d_cons_mean[i], d_cons_sd[i])
    
    # 3. Get sample size
    
    if (!is.null(data[[i]][[this_test_combo]]$n)) {
      d_n[i] <- data[[i]][[this_test_combo]]$n
    } else if (!is.null(data[[i]][[this_test_combo]]$n1) && !is.null(data[[i]][[this_test_combo]]$n2)) {
      d_n[i] <- data[[i]][[this_test_combo]]$n1 + data[[i]][[this_test_combo]]$n2
    } else {
      d_n[i] <- NA
    }
    
    # d_biased_sd[i] <- sd(d) * sqrt((length(d) - 1) / length(d)) # biased sd - doesn'd change results
    # TODO: catch normality test results
    # test for normality (usually light- or heavy-tailed, some approximately normal)
    # d_subset <- d[seq(1, length(d), length.out = 5000)] # max 5000 variables for shapiro test
    # d_normal_test[i] <- shapiro.test(d_subset)$p.value # p-value < 2e-16 -> not normal
    # qqnorm(d, pch=20, cex=0.5); qqline(d) # visualize
    
    # 4. Get multivariate effect size and bounds if available
    all_combos <- names(data[[i]])
    mv_combo_idx <- grep(this_test_mv_combo, all_combos)
    this_mv_combo <- all_combos[mv_combo_idx]
    d_mv[i] <- data[[i]][[this_mv_combo]]$d
    d_mv_lb[i] <- data[[i]][[this_mv_combo]]$sim_ci_lb
    d_mv_ub[i] <- data[[i]][[this_mv_combo]]$sim_ci_ub
  }
  
  # set up final data frame
  
  df <- data.frame(name = names(data), mean = d_mean, sd = d_sd, char_mag = d_characteristic_mag, mean_cons = d_cons_mean, sd_cons = d_cons_sd, char_mag_cons = d_characteristic_mag_cons, n = d_n, category = as.factor(study$category), dataset = I(study$dataset), ref = study$ref, mv = d_mv, mv_lb = d_mv_lb, mv_ub = d_mv_ub)
  
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
      category == "cognitive (task)" ~ "task (within-sub)",
      # category == "cognitive (task)" & grepl("voxel", ref) ~ "task activation",
      # category == "cognitive (task)" & !grepl("voxel", ref) ~ "task connectivity",
      TRUE ~ "other"
    ))
  df$overarching_category <- as.factor(df$overarching_category)
  
  return(df)
  
}



# Function for plotting effect sizes (mean, sd, n) for each study - point est and conservative  
plot_study <- function(df, df_meta, main_title, fn, plot_type = "sd") {

print(paste0('Fitting lines for ', main_title))

# params
nmax_plt <- 10000000
nmax <- 1000000000
points_only_dir <- 'points_only/'
plot_fitted_lines <- TRUE
n_pts <- 1000
use_char_mag <- FALSE

# Determine y variable and settings based on plot_type
if (plot_type == "sd") {
  y_var <- "sd"
  y_label <- "SD(d) across brain areas"
  fit_lines <- TRUE
  y_limits <- c(-0.05, 0.5)
} else if (plot_type == "mv") {
  y_var <- "mv"
  y_label <- "Multivariate effect size"
  # fit_lines <- FALSE
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
  
  # Only fit lines for sd plots
  # fit_options <- if (fit_lines) c(FALSE, TRUE) else c(FALSE)
  
  # for (plot_fitted_lines in fit_options) {
    
    # if (plot_fitted_lines && fit_lines) {
    if (plot_fitted_lines) {
      fitlines_str <- '__fitlines'
      
      # # interpolate predictions for all n
      n_seq <- data.frame(n=seq(0, nmax_plt, length.out=n_pts))
      n_seq[length(n_seq$n)+1,] <- nmax # ensure max n is included
      # predicted_sd_seq <- predict(fit, n_seq, interval="confidence")
      
      # set up for each overarching category
      unique_cats <- unique(df$overarching_category)
      predicted_y <- vector("list", length(unique_cats))
      names(predicted_y) <- unique_cats
      
      # preallocate beta
      beta <- vector("list", length(unique_cats))
      max_mv <- vector("list", length(unique_cats))
      res <- vector("list", length(unique_cats))
      names(beta) <- unique_cats
      
      for (cat in unique_cats) {
        
        # setup df
        df_cat <- df[df$overarching_category == cat, ]
        # df_cat <- df[df$overarching_category == cat & !is.na(df$mv), ]
        n_obs <- nrow(df_cat)
        n_levels <- length(unique(df_cat$dataset))
        
        # if only one unique dataset, do without random effect
        if (length(unique(df_cat$dataset)) > 1) {
          
          # fit
          if (plot_type == "mv") {
            fit_cat <- lmer(mv ~ 1 + (1|dataset), data = df_cat)
            d <- as.matrix(model.matrix(~ 1, data = n_seq))
          } else {
            if (use_char_mag) {
              fit_cat <- lmer(char_mag ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            } else {
              fit_cat <- lmer(sd ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            }
            # fit_cat <- lmer(char_mag ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            d <- model.matrix(~ I(1/sqrt(n)), data = n_seq)
          }
          
          # Use fixed effects for prediction
          # d <- as.matrix(model.matrix(~ 1, data = n_seq))
          beta[[cat]] <- fixef(fit_cat) #TODO: remove
          preds <- as.vector(d %*% beta[[cat]])
          preds_se <- sqrt(diag(d %*% vcov(fit_cat) %*% t(d)))
          lwr <- preds - 1.96 * preds_se
          upr <- preds + 1.96 * preds_se
          predicted_y[[cat]] <- cbind(fit = preds, lwr = lwr, upr = upr)
          
        } else {
          
          if (plot_type == "mv") {
            fit_cat <- lm(mv ~ 1, data=df_cat)
          } else {
            fit_cat <- lm(char_mag ~ I(1/sqrt(n)), data=df_cat)
          }
          predicted_y[[cat]] <- predict(fit_cat, n_seq, interval="confidence")
          beta[[cat]] <- coef(fit_cat) #TODO: remove
          
        }
        
        est <- summary(fit_cat)$coefficients[,"Estimate"]
        lwr <- summary(fit_cat)$coefficients[,"Estimate"] - 1.96 * summary(fit_cat)$coefficients[,"Std. Error"]
        upr <- summary(fit_cat)$coefficients[,"Estimate"] + 1.96 * summary(fit_cat)$coefficients[,"Std. Error"]
        
        
        # max_val[[cat]] <- predicted_y[[cat]][length(n_seq$n),]
        # res[[cat]] <- c(beta[[cat]], max_val[[cat]])
        
        if (plot_type == "mv") {
          res[[cat]] <- c(est = est, lwr = lwr, upr = upr)
        } else {
          res[[cat]] <- data.frame(est = est["(Intercept)"], lwr = lwr["(Intercept)"], upr = upr["(Intercept)"], row.names = NULL)
        }
        
      }
        
      
  # Combine all categories into a tidy data frame
  res <- do.call(rbind, res)

      fn_plot <- fn
      
    } else {
      # If not plotting fitlines, change output folder to 'points_only'
      fn_plot <- file.path(dirname(fn), points_only_dir, basename(fn))
      if (!dir.exists(dirname(fn_plot))) {
        dir.create(dirname(fn_plot), recursive = TRUE)
      }
      res <- NULL
    }
    
    # plot
    # ggplot2 implementation
    df$sqrt_n <- sqrt(df$n)
    df_meta$sqrt_n <- sqrt(df_meta$n)
    cats <- levels(df$overarching_category)
    color_map <- setNames(RColorBrewer::brewer.pal(length(cats), "Set1"), cats)

    # Create base plot with dynamic y variable
    p <- ggplot(df, aes_string(x = "sqrt_n", y = y_var, color = "overarching_category")) +
      geom_point(size = 1.5, alpha = alpha, stroke = 0) +
      scale_color_manual(values = color_map) +
      labs(title = main_title,
            x = "Sqrt(n)",
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
      scale_x_continuous(limits = c(0, max(sqrt(nmax_plt), max(sqrt(df_meta$n), na.rm=TRUE), max(sqrt(df$n), na.rm=TRUE))), expand = c(0, 0)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3)
    
    # Add y-axis limits only for sd plots
    if (!is.null(y_limits)) {
      p <- p + scale_y_continuous(limits = y_limits, expand = c(0, 0))
    }
    
    # Add mean points only for sd plots
    # if (plot_type == "sd") {
    #   p <- p + geom_point(aes(x = sqrt_n, y = mean), size = 1.5, alpha = alpha, stroke = 0, color = "black")
    # }

    # Add fitted lines only for sd plots when requested
    # if (plot_fitted_lines && fit_lines) {
    if (plot_fitted_lines) {
      
      # if (plot_type == "mv") {
        for (cat in unique_cats) {
          pred_mat <- predicted_y[[cat]]
          if (!is.null(pred_mat) && is.matrix(pred_mat) && all(c("fit","lwr","upr") %in% colnames(pred_mat)) && nrow(pred_mat) == length(n_seq$n)) {
            pred_df <- data.frame(sqrt_n = sqrt(n_seq$n),
                                  fit = pred_mat[,"fit"],
                                  lwr = pred_mat[,"lwr"],
                                  upr = pred_mat[,"upr"],
                                  overarching_category = cat)
            p <- p +
              geom_line(data = pred_df, aes(x = sqrt_n, y = fit, color = overarching_category), linewidth = 1, alpha = alpha) +
              geom_line(data = pred_df, aes(x = sqrt_n, y = lwr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha) +
              geom_line(data = pred_df, aes(x = sqrt_n, y = upr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha)
          }
        }
        
      # } else {
      #   for (cat in unique_cats) {
      #     pred_mat <- predicted_y[[cat]]
      #     if (!is.null(pred_mat) && is.matrix(pred_mat) && all(c("fit","lwr","upr") %in% colnames(pred_mat)) && nrow(pred_mat) == length(n_seq$n)) {
      #       pred_df <- data.frame(sqrt_n = sqrt(n_seq$n),
      #                             fit = pred_mat[,"fit"],
      #                             lwr = pred_mat[,"lwr"],
      #                             upr = pred_mat[,"upr"],
      #                             overarching_category = cat)
      #       p <- p +
      #         geom_line(data = pred_df, aes(x = sqrt_n, y = fit, color = overarching_category), linewidth = 1, alpha = alpha) +
      #         geom_line(data = pred_df, aes(x = sqrt_n, y = lwr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha) +
      #         geom_line(data = pred_df, aes(x = sqrt_n, y = upr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha)
      #     }
      #   }
      # }
    }

    if (add_meta && nrow(df_meta) > 0) {
      # Use gsub to format label with line break and study count
      df_meta$label <- paste0(gsub("_reference_", " (", df_meta$name), ")\n", df_meta$n_studies, " ", ifelse(df_meta$n_studies == 1, "study", "studies"))
      p <- p +
        geom_point(data = df_meta, aes_string(x = "sqrt_n", y = y_var, color = "overarching_category"), shape = 8, size = 2) +
        geom_text_repel(data = df_meta, aes_string(x = "sqrt_n", y = y_var, label = "label"), size = 2.2, segment.size = 0.3, force = 10, max.overlaps = Inf)
    }

    # Create filename suffix based on plot type and fitted lines
    type_suffix <- if (plot_type == "mv") "__mv" else ""
    # fitlines_suffix <- if (plot_fitted_lines && fit_lines) "__fitlinesY" else ""
    fitlines_suffix <- if (plot_fitted_lines) "__fitlines" else ""
    this_fn <- paste0(fn_plot, type_suffix, fitlines_suffix, meta_str, '.pdf')
    if (save_plots) {
      ggsave(this_fn, plot = p, width = 5, height = 4)
    } else {
      show(p)
    }
    
  # }
}

return(res)
}



# get req'd n and bin
make_required_n_df <- function(cats, d, sigmas_master, res_mv = NULL, mv = FALSE) {
  n_bins <- c(0, 25, 50, 100, 500, 1000, 5000, 50000, 300000, Inf)
  bin_labels <- paste(head(n_bins, -1), n_bins[-1]-1, sep = "–")
  required_n_df <- NULL
  for (cat in cats) {
    if (mv) {
      y <- rep(0, length(d))
      if (!is.null(res_mv) && cat %in% rownames(res_mv)) {
        closest_idx <- which.min(abs(d - res_mv[cat, "est"]))
        y[closest_idx] <- 1
      }
      # For mv, use two-sample logic for binning
      n_detect_two <- sapply(d, function(dd) pwr.t.test(power = 0.8, d = dd, sig.level = 0.05/2, type = "two.sample")$n)
      bin_indices <- cut(n_detect_two, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
    } else {
      sigmas <- sigmas_master[cat, ]
      y <- dnorm(d, mean = 0, sd = as.numeric(sigmas["ub"]))
      
      # Calculate both one-sample and two-sample n
      n_detect_one <- sapply(d, function(dd) pwr.t.test(power = 0.8, d = dd, sig.level = 0.05/2, type = "one.sample")$n)
      n_detect_two <- sapply(d, function(dd) pwr.t.test(power = 0.8, d = dd, sig.level = 0.05/2, type = "two.sample")$n)
      
      # Use n_detect_one for 'task (within-sub)', else n_detect_two
      if (cat == "task (within-sub)") {
        bin_indices <- cut(n_detect_one, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
      } else {
        bin_indices <- cut(n_detect_two, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
      }
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

plot_required_n_panel <- function(df, if_mv = FALSE) {
  title <- if (if_mv) "Required n for 80% Power by Category (Multivariate)" else "Required n for 80% Power by Category"
  filename <- if (if_mv) 'required_n_distributions__panels_mv.pdf' else 'required_n_distributions__panels.pdf'
  # Add line segment from y=0 to y=y[1] at x=1 if y[1] != 0 for each category
  segment_df <- df %>% group_by(overarching_category) %>% filter(row_number() == 1 & cumulative_proportion[1] != 0)
  p <- ggplot(df, aes(x = bin, y = cumulative_proportion, group = overarching_category, color = overarching_category)) +
    geom_line(size = 1) +
    scale_color_manual(values = cat_colors) +
    facet_wrap(~overarching_category, ncol = 1, scales = "free_y") +
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
    ggsave(paste0(fn_basedir, filename), p, width = 5, height = 15)
  } else {
    print(p)
  }
}
  
  ####### END OF HELPER FUNCTIONS ##########
  
for (pooling_method in pooling_methods) {
  print(paste0('plotting for pooling method: ', pooling_method))
  plots_name <- paste0('pooling.', pooling_method)

if (estimate == "d") {
  ci_lb <- "sim_ci_lb"
  ci_ub <- "sim_ci_ub"
} else if (estimate == "r_sq") {
  ci_lb <- "r_sq_sim_ci_lb"
  ci_ub <- "r_sq_sim_ci_ub"
}

# make output directory if it doesn't exist
fn_basedir <- paste0(output_path, plots_name,'/')
if (!dir.exists(fn_basedir)) {
  print(paste0('Creating output directory: ', fn_basedir))
  dir.create(fn_basedir, recursive = TRUE)
}




# get summaries for orig
summary_data <- get_study_summaries(v_data$data, v_data$study, motion_method_t, motion_method_not_t, meta = FALSE, estimate)

# make separate data frame with conservative estimates (facilitates reuse of later functions on cons est)
summary_data_cons <- summary_data
summary_data_cons$mean <- summary_data_cons$mean_cons
summary_data_cons$sd <- summary_data_cons$sd_cons
summary_data_cons$char_mag <- summary_data_cons$char_mag_cons

# get summaries for meta

# first, for meta, category is inexplicably group_level - rename 
v_data$meta_category$study$category <- v_data$meta_category$study$group_level

# Assign all relevant datasets, n's, and study names to meta_category$study
v_data$meta_category$study$dataset <- vector("list", length(v_data$meta_category$study$name))
v_data$meta_category$study$each_n <- vector("list", length(v_data$meta_category$study$name))
v_data$meta_category$study$n <- vector("list", length(v_data$meta_category$study$name)) # number of unique subjects
v_data$meta_category$study$included_study_names <- vector("list", length(v_data$meta_category$study$name))
v_data$meta_category$study$n_studies <- integer(length(v_data$meta_category$study$name))

  
# For each meta-analysis, get all studies included, their associated datasets, and sample sizes
for (i in seq_along(v_data$meta_category$study$name)) {
  print(paste0('Processing meta-analysis ', i, ' of ', length(v_data$meta_category$study$name), ': ', v_data$meta_category$study$name[i]))
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
  for (i in seq_along(v_data$meta_category$study$name)) {
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
    
  }
}
  
names(v_data$meta_category$study$dataset) <- v_data$meta_category$study$name
names(v_data$meta_category$study$each_n) <- v_data$meta_category$study$name
names(v_data$meta_category$study$included_study_names) <- v_data$meta_category$study$name
names(v_data$meta_category$study$n_studies) <- v_data$meta_category$study$name

# get summaries
summary_data__meta <- get_study_summaries(v_data$meta_category$data, v_data$meta_category$study, motion_method_t, motion_method_not_t, meta = TRUE, estimate)
summary_data_cons__meta <- summary_data__meta
summary_data_cons__meta$mean <- summary_data_cons__meta$mean_cons
summary_data_cons__meta$sd <- summary_data_cons__meta$sd_cons

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
  
  
plot_non_mv <- TRUE # TODO - tmp

if (plot_non_mv) {
  print('Fitting lines for point and conservative estimates')
# get point est and conservative est
res <- plot_study(summary_data, summary_data__meta, "Point Estimates", paste0(fn_basedir,'point_est'))
res_cons <- plot_study(summary_data_cons, summary_data_cons__meta, "Conservative Estimates", paste0(fn_basedir,'cons_est'))
  
# repeat, but for large N studies
  print('Fitting lines for large N studies only')
n_large <- 900
res_large <- plot_study(summary_data[summary_data$n > n_large,], summary_data__meta[summary_data__meta$n > n_large,],"Point Estimates (n > 900)", paste0(fn_basedir,'point_estimates_n900'))
}

# get multivariate effect sizes
  print('Fitting lines for multivariate effect sizes')
res_mv <- plot_study(summary_data, summary_data__meta, "Multivariate Effect Sizes", paste0(fn_basedir,'mv_est'), plot_type = "mv")


###########  Projected Density Plots ########### 

library(tidyr)
library(RColorBrewer)

n_pts <- 1000
cats <- c("physical", "psychological", "task (within-sub)")
cat_colors <- RColorBrewer::brewer.pal(length(cats), "Set1")
d <- seq(-2, 2, length.out = n_pts)
  
  for (do_mv in c(FALSE, TRUE)) {
  
  if (do_mv) {
    mv_suffix <- '_mv'
  } else {
    mv_suffix <- ''
  }

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
  
  # Set x limits and y max for each category
  if (do_mv) {
    xlim <- c(-2.5, 2.5)
  } else {
    xlim <- c(-0.8, 0.8)
  }
  # density_df <- density_df[density_df$d >= xlim_min & density_df$d <= xlim_max, ]
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
    ggsave(paste0(fn_basedir, 'projected_distributions',mv_suffix,'.pdf'), p_density, width = 5, height = 4)
  } else {
    print(p_density)
  }
  
  # Separate panels per category, colored by overarching category, legend inset
  p_density_panel <- ggplot(density_df, aes(x = d, y = density, color = overarching_category, linetype = sigma_type)) +
    geom_line(size = 1) +
    facet_wrap(~category, nrow = 1, scales = "free_y") +
    scale_color_manual(values = cat_colors) +
    labs(title = "Density Curves by Category", x = "Cohen's d", y = "Density", color = "Category", linetype = "Sigma Type") +
    theme_bw() +
    coord_cartesian(xlim = xlim) +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key.size = unit(0.7, "lines")
    )
  if (save_plots) {
    ggsave(paste0(fn_basedir, 'projected_distributions__panels',mv_suffix,'.pdf'), p_density_panel, width = 12, height = 4)
  } else {
    print(p_density_panel)
  }
  
  # Save sigmas_master to file
  if (save_plots) {
    if (do_mv) {
        write.csv(res_mv, file=paste0(fn_basedir, 'projected_distributions__mus',mv_suffix,'.csv'), row.names=TRUE)
    } else {
      write.csv(sigmas_master, file=paste0(fn_basedir, 'projected_distributions__sigmas',mv_suffix,'.csv'), row.names=TRUE)
    }
  }

}
  
  # Normal setup
required_n_df <- make_required_n_df(cats, d, sigmas_master, mv = FALSE)
plot_required_n_panel(required_n_df, if_mv = FALSE)

# Multivariate setup
required_n_df_mv <- make_required_n_df(cats, d, sigmas_master, res_mv = res_mv, mv = TRUE)
plot_required_n_panel(required_n_df_mv, if_mv = TRUE)
  
  
  
  
} # pooling method loop
} # function
# ------------------------------------------------------------------



# helper functions to source first:

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

get_study_summaries <- function(data, study, motion_method_t, motion_method_not_t, meta, estimate) {
  
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
  
  for (i in seq_along(data)) {
    # find test type for this study
    this_test_type = study[i, "orig_stat_type"]
    this_motion = ifelse(this_test_type == "t", motion_method_t, motion_method_not_t) #TODO make dynamic with fn inputs
    
    if (meta) {
      # try regression motion first, if null then switch to threshold
      this_motion = "regression"
    }
    this_test_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.none")
    this_test_mv_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.multi")

    # check if estimate for the given study and combo is NULL
    if (anyNA(data[[i]][[this_test_combo]][[estimate]]) | anyNA(data[[i]][[this_test_mv_combo]][[estimate]])) {
      this_motion = "threshold"
      this_test_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.none")
      this_test_mv_combo = paste0("pooling.", pooling_method, ".motion.", this_motion, ".mv.multi")
    }

    # if estimate is still NULL, throw an error
    if (anyNA(data[[i]][[this_test_combo]][[estimate]])) {
      stop(paste0("Estimate '", estimate, "' is NULL for study ", i, ": ", study$name[i], " with test type ", this_test_type, " even after switching to threshold motion. Please check the data."))
    }

    print(paste0("using motion method ", this_motion, " for study ", i, ": ", study$name[i]))

    # 1. Get point estimate mean and sd across vars
    
    # to get d, check whether dim is null (nested list)
    if (is.null(dim(data[[i]][[this_test_combo]][[estimate]]))) {
      d <- data[[i]][[this_test_combo]][[estimate]]
    } else {
      d <- data[[i]][[this_test_combo]][[estimate]][1,]
    }
    
    d_mean[i] <- mean(d)
    d_sd[i] <- sd(d)
    d_characteristic_mag[i] <- get_characteristic_magnitude(d_mean[i], d_sd[i])
    
    # 2. Get conservative estimate mean and sd across vars
    # from prep_data_for_plot.R
    # omitting na check in mean and sd and downsampling
    
    # unlist sim CIs if list
    if (is.list(data[[i]][[this_test_combo]][[ci_lb]])) {
      data[[i]][[this_test_combo]][[ci_lb]] <- unlist(data[[i]][[this_test_combo]][[ci_lb]])
      data[[i]][[this_test_combo]][[ci_ub]] <- unlist(data[[i]][[this_test_combo]][[ci_ub]])
    }
    
    na_idx <- is.na(data[[i]][[this_test_combo]][[estimate]]) | is.na(data[[i]][[this_test_combo]][[ci_lb]]) | is.na(data[[i]][[this_test_combo]][[ci_ub]])
    data[[i]][[this_test_combo]][[estimate]] <- data[[i]][[this_test_combo]][[estimate]][!na_idx]
    data[[i]][[this_test_combo]][[ci_lb]] <- data[[i]][[this_test_combo]][[ci_lb]][!na_idx]
    data[[i]][[this_test_combo]][[ci_ub]] <- data[[i]][[this_test_combo]][[ci_ub]][!na_idx]
    
    # sort data from smallest to largest effect size
    sorted_indices <- order(data[[i]][[this_test_combo]][[estimate]])
    sorted_estimate <- data[[i]][[this_test_combo]][[estimate]][sorted_indices]
    sorted_upper_bounds <- data[[i]][[this_test_combo]][[ci_ub]][sorted_indices]
    sorted_lower_bounds <- data[[i]][[this_test_combo]][[ci_lb]][sorted_indices]
    
    sorted_cons_estimate <- ifelse((abs(sorted_lower_bounds) > abs(sorted_upper_bounds)),
                                   ifelse((sorted_upper_bounds < 0),
                                          round(sorted_upper_bounds, 2), 0),
                                   ifelse((sorted_lower_bounds > 0),
                                          round(sorted_lower_bounds, 2), 0))
    
    
    d_cons_mean[i] <- mean(sorted_cons_estimate)
    d_cons_sd[i] <- sd(sorted_cons_estimate)
    d_characteristic_mag_cons[i] <- get_characteristic_magnitude(d_cons_mean[i], d_cons_sd[i])
    
    # 3. Get sample size
    
    if (!is.null(data[[i]][[this_test_combo]]$n)) {
      d_n[i] <- data[[i]][[this_test_combo]]$n
    } else if (!is.null(data[[i]][[this_test_combo]]$n1) && !is.null(data[[i]][[this_test_combo]]$n2)) {
      d_n[i] <- data[[i]][[this_test_combo]]$n1 + data[[i]][[this_test_combo]]$n2
    } else {
      d_n[i] <- NA
    }
    
    # d_biased_sd[i] <- sd(d) * sqrt((length(d) - 1) / length(d)) # biased sd - doesn'd change results
    # TODO: catch normality test results
    # test for normality (usually light- or heavy-tailed, some approximately normal)
    # d_subset <- d[seq(1, length(d), length.out = 5000)] # max 5000 variables for shapiro test
    # d_normal_test[i] <- shapiro.test(d_subset)$p.value # p-value < 2e-16 -> not normal
    # qqnorm(d, pch=20, cex=0.5); qqline(d) # visualize
    
    # 4. Get multivariate effect size and bounds if available
    all_combos <- names(data[[i]])
    mv_combo_idx <- grep(this_test_mv_combo, all_combos)
    this_mv_combo <- all_combos[mv_combo_idx]
    d_mv[i] <- data[[i]][[this_mv_combo]]$d
    d_mv_lb[i] <- data[[i]][[this_mv_combo]]$sim_ci_lb
    d_mv_ub[i] <- data[[i]][[this_mv_combo]]$sim_ci_ub
  }
  
  # set up final data frame
  
  df <- data.frame(name = names(data), mean = d_mean, sd = d_sd, char_mag = d_characteristic_mag, mean_cons = d_cons_mean, sd_cons = d_cons_sd, char_mag_cons = d_characteristic_mag_cons, n = d_n, category = as.factor(study$category), dataset = I(study$dataset), ref = study$ref, mv = d_mv, mv_lb = d_mv_lb, mv_ub = d_mv_ub)
  
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
      category == "cognitive (task)" ~ "task (within-sub)",
      # category == "cognitive (task)" & grepl("voxel", ref) ~ "task activation",
      # category == "cognitive (task)" & !grepl("voxel", ref) ~ "task connectivity",
      TRUE ~ "other"
    ))
  df$overarching_category <- as.factor(df$overarching_category)
  
  return(df)
  
}



# Function for plotting effect sizes (mean, sd, n) for each study - point est and conservative  
plot_study <- function(df, df_meta, main_title, fn, plot_type = "sd") {

print(paste0('Fitting lines for ', main_title))

# params
nmax_plt <- 10000000
nmax <- 1000000000
points_only_dir <- 'points_only/'
plot_fitted_lines <- TRUE
n_pts <- 1000
use_char_mag <- FALSE

# Determine y variable and settings based on plot_type
if (plot_type == "sd") {
  y_var <- "sd"
  y_label <- "SD(d) across brain areas"
  fit_lines <- TRUE
  y_limits <- c(-0.05, 0.5)
} else if (plot_type == "mv") {
  y_var <- "mv"
  y_label <- "Multivariate effect size"
  # fit_lines <- FALSE
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
  
  # Only fit lines for sd plots
  # fit_options <- if (fit_lines) c(FALSE, TRUE) else c(FALSE)
  
  # for (plot_fitted_lines in fit_options) {
    
    # if (plot_fitted_lines && fit_lines) {
    if (plot_fitted_lines) {
      fitlines_str <- '__fitlines'
      
      # # interpolate predictions for all n
      n_seq <- data.frame(n=seq(0, nmax_plt, length.out=n_pts))
      n_seq[length(n_seq$n)+1,] <- nmax # ensure max n is included
      # predicted_sd_seq <- predict(fit, n_seq, interval="confidence")
      
      # set up for each overarching category
      unique_cats <- unique(df$overarching_category)
      predicted_y <- vector("list", length(unique_cats))
      names(predicted_y) <- unique_cats
      
      # preallocate beta
      beta <- vector("list", length(unique_cats))
      max_mv <- vector("list", length(unique_cats))
      res <- vector("list", length(unique_cats))
      names(beta) <- unique_cats
      
      for (cat in unique_cats) {
        
        # setup df
        df_cat <- df[df$overarching_category == cat, ]
        # df_cat <- df[df$overarching_category == cat & !is.na(df$mv), ]
        n_obs <- nrow(df_cat)
        n_levels <- length(unique(df_cat$dataset))
        
        # if only one unique dataset, do without random effect
        if (length(unique(df_cat$dataset)) > 1) {
          
          # fit
          if (plot_type == "mv") {
            fit_cat <- lmer(mv ~ 1 + (1|dataset), data = df_cat)
            d <- as.matrix(model.matrix(~ 1, data = n_seq))
          } else {
            if (use_char_mag) {
              fit_cat <- lmer(char_mag ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            } else {
              fit_cat <- lmer(sd ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            }
            # fit_cat <- lmer(char_mag ~ I(1/sqrt(n)) + (1|dataset), data = df_cat)
            d <- model.matrix(~ I(1/sqrt(n)), data = n_seq)
          }
          
          # Use fixed effects for prediction
          # d <- as.matrix(model.matrix(~ 1, data = n_seq))
          beta[[cat]] <- fixef(fit_cat) #TODO: remove
          preds <- as.vector(d %*% beta[[cat]])
          preds_se <- sqrt(diag(d %*% vcov(fit_cat) %*% t(d)))
          lwr <- preds - 1.96 * preds_se
          upr <- preds + 1.96 * preds_se
          predicted_y[[cat]] <- cbind(fit = preds, lwr = lwr, upr = upr)
          
        } else {
          
          if (plot_type == "mv") {
            fit_cat <- lm(mv ~ 1, data=df_cat)
          } else {
            fit_cat <- lm(char_mag ~ I(1/sqrt(n)), data=df_cat)
          }
          predicted_y[[cat]] <- predict(fit_cat, n_seq, interval="confidence")
          beta[[cat]] <- coef(fit_cat) #TODO: remove
          
        }
        
        est <- summary(fit_cat)$coefficients[,"Estimate"]
        lwr <- summary(fit_cat)$coefficients[,"Estimate"] - 1.96 * summary(fit_cat)$coefficients[,"Std. Error"]
        upr <- summary(fit_cat)$coefficients[,"Estimate"] + 1.96 * summary(fit_cat)$coefficients[,"Std. Error"]
        
        
        # max_val[[cat]] <- predicted_y[[cat]][length(n_seq$n),]
        # res[[cat]] <- c(beta[[cat]], max_val[[cat]])
        
        if (plot_type == "mv") {
          res[[cat]] <- c(est = est, lwr = lwr, upr = upr)
        } else {
          res[[cat]] <- data.frame(est = est["(Intercept)"], lwr = lwr["(Intercept)"], upr = upr["(Intercept)"], row.names = NULL)
        }
        
      }
        
      
  # Combine all categories into a tidy data frame
  res <- do.call(rbind, res)

      fn_plot <- fn
      
    } else {
      # If not plotting fitlines, change output folder to 'points_only'
      fn_plot <- file.path(dirname(fn), points_only_dir, basename(fn))
      if (!dir.exists(dirname(fn_plot))) {
        dir.create(dirname(fn_plot), recursive = TRUE)
      }
      res <- NULL
    }
    
    # plot
    # ggplot2 implementation
    df$sqrt_n <- sqrt(df$n)
    df_meta$sqrt_n <- sqrt(df_meta$n)
    cats <- levels(df$overarching_category)
    color_map <- setNames(RColorBrewer::brewer.pal(length(cats), "Set1"), cats)

    # Create base plot with dynamic y variable
    p <- ggplot(df, aes_string(x = "sqrt_n", y = y_var, color = "overarching_category")) +
      geom_point(size = 1.5, alpha = alpha, stroke = 0) +
      scale_color_manual(values = color_map) +
      labs(title = main_title,
            x = "Sqrt(n)",
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
      scale_x_continuous(limits = c(0, max(sqrt(nmax_plt), max(sqrt(df_meta$n), na.rm=TRUE), max(sqrt(df$n), na.rm=TRUE))), expand = c(0, 0)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3)
    
    # Add y-axis limits only for sd plots
    if (!is.null(y_limits)) {
      p <- p + scale_y_continuous(limits = y_limits, expand = c(0, 0))
    }
    
    # Add mean points only for sd plots
    # if (plot_type == "sd") {
    #   p <- p + geom_point(aes(x = sqrt_n, y = mean), size = 1.5, alpha = alpha, stroke = 0, color = "black")
    # }

    # Add fitted lines only for sd plots when requested
    # if (plot_fitted_lines && fit_lines) {
    if (plot_fitted_lines) {
      
      # if (plot_type == "mv") {
        for (cat in unique_cats) {
          pred_mat <- predicted_y[[cat]]
          if (!is.null(pred_mat) && is.matrix(pred_mat) && all(c("fit","lwr","upr") %in% colnames(pred_mat)) && nrow(pred_mat) == length(n_seq$n)) {
            pred_df <- data.frame(sqrt_n = sqrt(n_seq$n),
                                  fit = pred_mat[,"fit"],
                                  lwr = pred_mat[,"lwr"],
                                  upr = pred_mat[,"upr"],
                                  overarching_category = cat)
            p <- p +
              geom_line(data = pred_df, aes(x = sqrt_n, y = fit, color = overarching_category), linewidth = 1, alpha = alpha) +
              geom_line(data = pred_df, aes(x = sqrt_n, y = lwr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha) +
              geom_line(data = pred_df, aes(x = sqrt_n, y = upr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha)
          }
        }
        
      # } else {
      #   for (cat in unique_cats) {
      #     pred_mat <- predicted_y[[cat]]
      #     if (!is.null(pred_mat) && is.matrix(pred_mat) && all(c("fit","lwr","upr") %in% colnames(pred_mat)) && nrow(pred_mat) == length(n_seq$n)) {
      #       pred_df <- data.frame(sqrt_n = sqrt(n_seq$n),
      #                             fit = pred_mat[,"fit"],
      #                             lwr = pred_mat[,"lwr"],
      #                             upr = pred_mat[,"upr"],
      #                             overarching_category = cat)
      #       p <- p +
      #         geom_line(data = pred_df, aes(x = sqrt_n, y = fit, color = overarching_category), linewidth = 1, alpha = alpha) +
      #         geom_line(data = pred_df, aes(x = sqrt_n, y = lwr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha) +
      #         geom_line(data = pred_df, aes(x = sqrt_n, y = upr, color = overarching_category), linetype = "dotted", linewidth = 0.8, alpha = alpha)
      #     }
      #   }
      # }
    }

    if (add_meta && nrow(df_meta) > 0) {
      # Use gsub to format label with line break and study count
      df_meta$label <- paste0(gsub("_reference_", " (", df_meta$name), ")\n", df_meta$n_studies, " ", ifelse(df_meta$n_studies == 1, "study", "studies"))
      p <- p +
        geom_point(data = df_meta, aes_string(x = "sqrt_n", y = y_var, color = "overarching_category"), shape = 8, size = 2) +
        geom_text_repel(data = df_meta, aes_string(x = "sqrt_n", y = y_var, label = "label"), size = 2.2, segment.size = 0.3, force = 10, max.overlaps = Inf)
    }

    # Create filename suffix based on plot type and fitted lines
    type_suffix <- if (plot_type == "mv") "__mv" else ""
    # fitlines_suffix <- if (plot_fitted_lines && fit_lines) "__fitlinesY" else ""
    fitlines_suffix <- if (plot_fitted_lines) "__fitlines" else ""
    this_fn <- paste0(fn_plot, type_suffix, fitlines_suffix, meta_str, '.pdf')
    if (save_plots) {
      ggsave(this_fn, plot = p, width = 5, height = 4)
    } else {
      show(p)
    }
    
  # }
}

return(res)
}



# get req'd n and bin
make_required_n_df <- function(cats, d, sigmas_master, res_mv = NULL, mv = FALSE) {
  n_bins <- c(0, 25, 50, 100, 500, 1000, 5000, 50000, 300000, Inf)
  bin_labels <- paste(head(n_bins, -1), n_bins[-1]-1, sep = "–")
  required_n_df <- NULL
  for (cat in cats) {
    if (mv) {
      y <- rep(0, length(d))
      if (!is.null(res_mv) && cat %in% rownames(res_mv)) {
        closest_idx <- which.min(abs(d - res_mv[cat, "est"]))
        y[closest_idx] <- 1
      }
      # For mv, use two-sample logic for binning
      n_detect_two <- sapply(d, function(dd) pwr.t.test(power = 0.8, d = dd, sig.level = 0.05/2, type = "two.sample")$n)
      bin_indices <- cut(n_detect_two, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
    } else {
      sigmas <- sigmas_master[cat, ]
      y <- dnorm(d, mean = 0, sd = as.numeric(sigmas["ub"]))
      
      # Calculate both one-sample and two-sample n
      n_detect_one <- sapply(d, function(dd) pwr.t.test(power = 0.8, d = dd, sig.level = 0.05/2, type = "one.sample")$n)
      n_detect_two <- sapply(d, function(dd) pwr.t.test(power = 0.8, d = dd, sig.level = 0.05/2, type = "two.sample")$n)
      
      # Use n_detect_one for 'task (within-sub)', else n_detect_two
      if (cat == "task (within-sub)") {
        bin_indices <- cut(n_detect_one, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
      } else {
        bin_indices <- cut(n_detect_two, breaks = n_bins, include.lowest = TRUE, labels = bin_labels)
      }
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

plot_required_n_panel <- function(df, if_mv = FALSE) {
  title <- if (if_mv) "Required n for 80% Power by Category (Multivariate)" else "Required n for 80% Power by Category"
  filename <- if (if_mv) 'required_n_distributions__panels_mv.pdf' else 'required_n_distributions__panels.pdf'
  # Add line segment from y=0 to y=y[1] at x=1 if y[1] != 0 for each category
  segment_df <- df %>% group_by(overarching_category) %>% filter(row_number() == 1 & cumulative_proportion[1] != 0)
  p <- ggplot(df, aes(x = bin, y = cumulative_proportion, group = overarching_category, color = overarching_category)) +
    geom_line(size = 1) +
    scale_color_manual(values = cat_colors) +
    facet_wrap(~overarching_category, ncol = 1, scales = "free_y") +
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
    ggsave(paste0(fn_basedir, filename), p, width = 5, height = 15)
  } else {
    print(p)
  }
}