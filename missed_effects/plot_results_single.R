plot_results_single <- function(
  sim_params,
  sim_output,
  cat_colors,
  out_master_dir,
  text_size = 20,
  ticks_size = 24,
  transparency_main = 0.6,
  transparency_overlay = 0.5
) {

list2env(sim_output, envir = environment())
cat_color <- cat_colors[sim_params$outcome_category]

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
  # list(source = "target", col = "tpr_target_univ",      label = "Power for Target (Univariate)",                    ylim = c(0, 1), plot = TRUE),
  # list(source = "target", col = "type_m_target_univ",   label = "Type M Error for Target (Univariate)",             ylim = c(0, 10), plot = TRUE),
  # list(source = "target", col = "type_s_target_univ",   label = "Type S Error for Target (Univariate)",             ylim = c(0, 0.5), plot = TRUE),
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
    # basis_study_df$tpr_target_univ[this_rep]    <- positive_rate_univ[target_effect_idx__actual_mean] # this is only for the mean, not max or same sample
    # basis_study_df$type_m_target_univ[this_rep] <- abs(eff_mat_univ[this_rep, target_effect_idx__actual_mean]) / ground_truth_true_effects[target_effect_idx__actual_mean]
    # basis_study_df$type_s_target_univ[this_rep] <-
    #   (eff_mat_univ[this_rep, target_effect_idx__actual_mean] < 0 & ground_truth_true_effects[target_effect_idx__actual_mean] > 0) |
    #     (eff_mat_univ[this_rep, target_effect_idx__actual_mean] > 0 & ground_truth_true_effects[target_effect_idx__actual_mean] < 0)
    
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

plot_summary <- function(summary, metric, metric_label, y_limits, cat_color, text_size, ticks_size, out_dir = "~/Desktop/sim/") {
  
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
    geom_line(size = 0.5, color = cat_color) +
    scale_x_continuous(breaks = df_plot$x_num, labels = df_plot$sample_size)
  
  # add ribbon if sd exists
  if (!is.null(str_sd)) {
    p <- p + geom_ribbon(aes(x = x_num, ymin = lb, ymax = ub), alpha = transparency_main, fill = cat_color, inherit.aes = FALSE) # nolint
  }
    
  p <- p +
    coord_cartesian(ylim = y_limits) +
    labs(title = metric_label, x = x_label, y = metric_label) +
    theme_bw(base_size = ticks_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = text_size),
      axis.title = element_text(size = text_size)
    )

  ggsave(filename = paste0(out_dir, metric, ".png"), plot = p, width = 6, height = 4)
}

# plots — driven directly from specs, so adding a spec entry automatically adds a plot
for (spec in specs) {
  if (isTRUE(spec$plot)) {
    plot_summary(summary, summary_col_name(spec), spec$label, spec$ylim, cat_color = cat_colors[sim_params$outcome_category], text_size, ticks_size, out_dir = out_dir)
  }
}
# save summary variable if doesn't exist
# if (!file.exists(summary_file)) {
  saveRDS(summary, file = summary_file)
# }


######### OVERLAY PLOT: Expected TPs (Uncorr) + Overlap by Category #########

summary__this_cat <- summary
summary__this_cat$category <- sim_params$outcome_category

if (!is.null(summary__this_cat) && nrow(summary__this_cat) > 0) {

  # Map sample_size to integer x positions
  all_sample_sizes <- sort(unique(summary__this_cat$sample_size))
  summary__this_cat$x_num <- match(summary__this_cat$sample_size, all_sample_sizes)

  # Secondary x-axis labels: expected_n__replication for each sample_size
  sec_x_breaks <- seq_along(all_sample_sizes)
  sec_x_labels <- round(summary__this_cat$expected_n__replication_mean[match(all_sample_sizes, summary__this_cat$sample_size)])

  left_max  <- 1.5

  # Styling: right-axis series uses a dark-grey tint of left-axis color (low saturation, reduced brightness)
  left_color <- cat_color
  this_hsv <- rgb2hsv(col2rgb(left_color))
  right_color <- hsv(h = this_hsv[1, 1], s = this_hsv[2, 1] * 0.3, v = this_hsv[3, 1] * 0.6)

  p_overlay <- ggplot(summary__this_cat, aes(x = x_num)) +
    # Left y-axis: proportion of expected TPs detected (uncorr) — solid line with ribbon
    geom_ribbon(
      aes(x = x_num,
          ymin = expect_v_actual_n_tp__based_on_basis_mean - expect_v_actual_n_tp__based_on_basis_sd,
          ymax = expect_v_actual_n_tp__based_on_basis_mean + expect_v_actual_n_tp__based_on_basis_sd),
      fill = left_color, alpha = transparency_main, colour = NA, inherit.aes = FALSE
    ) +
    geom_ribbon(
      aes(x = x_num,
          ymin = overlap_mean - overlap_sd,
          ymax = overlap_mean + overlap_sd),
      fill = right_color, alpha = transparency_overlay, colour = NA, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 1, colour = "grey50", linetype = "dotted", linewidth = 0.5) +

    geom_line(aes(y = expect_v_actual_n_tp__based_on_basis_mean), linewidth = 0.7, colour = left_color) +
    geom_line(aes(y = overlap_mean), linewidth = 0.7, colour = right_color) +
    # Y-axis
    scale_y_continuous(
      name   = "Proportion Detected"
    ) +
    coord_cartesian(ylim = c(0, left_max)) +
    # X-axes: primary (bottom) = basis study sample size; secondary (top) = planned main N
    scale_x_continuous(
      name   = "Planned Sample Size (Main Study)",
      breaks = sec_x_breaks,
      labels = sec_x_labels,
      sec.axis = sec_axis(
        transform = ~ .,
        breaks = sec_x_breaks,
        labels = all_sample_sizes,
        name   = "Sample Size of Basis Study"
      )
    ) +
    theme_bw(base_size = ticks_size) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = text_size),
      axis.title.y.left  = element_text(size = text_size, color = "black"),
      axis.text.y.left   = element_text(color = "black")
    )

  ggsave(
    filename = paste0(out_master_dir, sim_params$outcome_category, "/",sim_params$target_effect_type__from_basis, "_effect/overlay_tp_and_overlap.png"),
    plot = p_overlay, width = 8, height = 5.5
  )

} else {
  message("No category summary files found; skipping overlay plot.")
}

  return(invisible(list(
    summary = summary,
    out_dir = out_dir,
    summary_file = summary_file
  )))
}


