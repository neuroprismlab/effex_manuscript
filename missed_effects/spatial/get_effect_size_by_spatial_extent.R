
### Libraries

library(dplyr)
library(tidyr)
library(ggplot2)
library(metafor)
library(tibble) # for rownames_to_column

### Set params & filenames
pooling_type <- "none"
motion_type <- "threshold"
multi_string <- "multi" # "multi" or "mv\\.none"
use_high_sample_size_only <- FALSE

# /Users/stephanienoble/
data_dir <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My\ Drive/Lab/xMore/Software/scripts/R/myscripts/effect_size/BrainEffeX_utils/inst/meta/'
results_dir_master <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My\ Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/manuscript/figures/plots/crossbrain_effects__spatial_extent_d2/'

# name results_dir based on motion and pooling
results_dir <- paste0(results_dir_master, "pooling.", pooling_type, ".motion.", motion_type, "/")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# for combining models
combine_models <- TRUE
do_r2 <- FALSE
spatial_extents <- c("01","05", "10", "25", "50", "75", "100")
categories <- c("psychological", "physical", "task activation", "task connectivity")

# plotting params
cat_colors <- setNames(RColorBrewer::brewer.pal(length(categories), "Set1"), categories)
cat_colors[c("physical", "psychological")] <- cat_colors[c("psychological", "physical")] # switch to match usual colors
axis_text_size <- 16
axis_title_size <- 16
line_width_main <- 1.4
line_width_ci <- 1.1
line_width_ref <- 0.6


##### LOAD & ORGANIZE ####

repeat_overwrite_for_all <- FALSE
overwrite_all <- NA
load_data <- TRUE

for (this_extent in spatial_extents) {
  
  # TODO: also need to check whether previous extracted data (e.g., pooling=net) matches the present one, or it will look like it's already loaded
  
  load_data <- TRUE
  # if (exists(paste0("study_level_data_", this_extent))) {
  #   if (exists("existing_pooling_param") && (existing_pooling_param != pooling_type || existing_motion_param != motion_type)) {
  #     load_data <- TRUE
  #   } else {
  #     if (isTRUE(repeat_overwrite_for_all) && !is.na(overwrite_all)) {
  #       response <- if (overwrite_all) "y" else "n"
  #       # cat(paste0("Using saved response for ", this_extent, "%: ", response, "\n"))
  #     } else {
  #       response <- readline(prompt = paste0("Data for ", this_extent, "% was already extracted. Do you want to re-load & extract this data? (y/n): "))
  #       response_repeat <- readline(prompt = "Repeat your choice for all remaining extents? (y/n): ")
  #       if (tolower(response_repeat) == "y") {
  #         repeat_overwrite_for_all <- TRUE
  #         overwrite_all <- tolower(response) == "y"
  #       }
  #     }
  # 
  #     if (tolower(response) != "y") {
  #       cat(paste0("Skipping loading for ", this_extent, "%\n"))
  #       load_data <- FALSE
  #     } else {
  #       cat(paste0("Overwriting data for ", this_extent, "%\n"))
  #       load_data <- TRUE
  #     }
  #   }
  # 
  # } else {
  #   load_data <- TRUE
  #   existing_pooling_param <- pooling_type
  #   existing_motion_param <- motion_type
  # }
  
  if (load_data) {
    v_filename <- list.files(data_dir, pattern = paste0("braineffex_data_", this_extent, "_percent_.*\\.RData$"), full.names = TRUE)
    
    # load data
    # v_filename <- paste0(data_dir, "v.RData")
    load(v_filename)

    extract_metric <- function(dat, metric, pooling_type, motion_type) {
      unlist(lapply(dat, function(x) {
        x1 <- x[grepl("multi", names(x))]
        x1 <- x1[grepl(paste0("pooling\\.", pooling_type), names(x1))]
        lapply(x1[grepl(paste0("motion\\.", motion_type), names(x1))], `[[`, metric)
      }), recursive = FALSE)
    }
    
    # extract data to new data frame 
    d_list <- extract_metric(v$data, "d", pooling_type, motion_type)
    n_list <- extract_metric(v$data, "n", pooling_type, motion_type)
    n1_list <- extract_metric(v$data, "n1", pooling_type, motion_type)
    n2_list <- extract_metric(v$data, "n2", pooling_type, motion_type)
    
    study_name <- sub("\\..*$", "", names(d_list))
    df <- v$study
    
    # remove rows of df that don't match study_name df[!match(tolower(df$name), study_name), ]
    df <- df[match(study_name,tolower(df$name)), ]
    
    df$d <- unlist(d_list[match(tolower(df$name), study_name)])
    df$n <- unlist(n_list[match(tolower(df$name), study_name)])
    df$n1 <- unlist(n1_list[match(tolower(df$name), study_name)])
    df$n2 <- unlist(n2_list[match(tolower(df$name), study_name)])
    df$k2 <- ifelse(df$orig_stat_type == "t", 1, 4)
    
    # remove tests and duplicates
    df <- df[!grepl("test",df$name),]
    df <- df[!duplicated(df$name), ]
    
    # fix categories and map to overarching
    df$category[df$category == "clinical"] <- "psychiatric"
    df$category[df$test_component_2 == "bmi"] <- "biometric"
    df$category[grepl("sex", tolower(df$test_component_2))] <- "sex (demographic)"
    df$category[grepl("gender", tolower(df$test_component_2))] <- "sex (demographic)"
    df$category[df$test_component_2 == "cbcl_scr_syn_internal_t_FU1"] <- "psychiatric" # weird random bug
    df$category[grep("age", df$name)] <- "age (demographic)"
    df$category[df$category == "cognitive" & df$orig_stat_type == "t"] <- "cognitive (task)"

    # silly fix when using the old (not "updated" abcd data)
    need_to_flip <- which(grepl("abcd", tolower(df$dataset)) & df$map_type == "fc" & df$orig_stat_type == "t" & tolower(df$test_component_2) != "rest")
    # since multivariate don't need to change sign, but note for univariate will need to do: df$d[need_to_flip] <- -df$d[need_to_flip]
    if (length(need_to_flip) > 0) {
      tmp <- df$test_component_1[need_to_flip]
      df$test_component_1[need_to_flip] <- df$test_component_2[need_to_flip]
      df$test_component_2[need_to_flip] <- tmp
      df$name[need_to_flip] <- sapply(need_to_flip, function(i) {
        new_name <- gsub(df$test_component_1[i], "tmp", df$name[i])
        new_name <- gsub(df$test_component_2[i], df$test_component_1[i], new_name)
        new_name <- gsub("tmp", df$test_component_2[i], new_name)
        new_name
      })
      df$basefile[need_to_flip] <- sapply(need_to_flip, function(i) {
        new_basefile <- gsub(df$test_component_1[i], "tmp", df$basefile[i])
        new_basefile <- gsub(df$test_component_2[i], df$test_component_1[i], new_basefile)
        new_basefile <- gsub("tmp", df$test_component_2[i], new_basefile)
        new_basefile
      })
    }
    
    
    df <- df %>%
      mutate(overarching_category = case_when(
        category %in% c("biometric", "sex (demographic)", "age (demographic)") ~ "physical",
        category %in% c("cognitive", "psychiatric") ~ "psychological",
        # category == "cognitive (task)" ~ "task (within-sub)",
        category == "cognitive (task)" & !grepl("act", map_type) ~ "task connectivity",
        category == "cognitive (task)" & grepl("act", map_type) ~ "task activation",
        TRUE ~ "other"
      ))
    
    df$overarching_category <- as.factor(df$overarching_category)
    
    
    # average within overarching_category as a quick check
    # df_summary <- df %>%
    #   group_by(overarching_category) %>%
    #   summarise(mean_d = mean(d, na.rm = TRUE), .groups = "drop")
    
    assign(paste0("study_level_data_", this_extent), df)
    
  }
}

# get studies that are shared across all study_level_data, and also get missing studies
common_studies <- Reduce(intersect, lapply(spatial_extents, function(this_extent) tolower(get(paste0("study_level_data_", this_extent))$name)))

missing_studies <- lapply(spatial_extents, function(this_extent) {
  df <- get(paste0("study_level_data_", this_extent))
  setdiff(tolower(df$name), common_studies)
})

# summarize for each missing study which extents have or don't have that study
missing_studies_df <- data.frame(
  extent = rep(spatial_extents, sapply(missing_studies, length)),
  missing_study = unlist(missing_studies)
)
missing_studies_df <- missing_studies_df %>%
  group_by(missing_study) %>%
  summarise(missing_extents = paste(extent, collapse = ", "), .groups = "drop") %>%
  rowwise() %>%
  mutate(present_extents = paste(setdiff(spatial_extents, strsplit(missing_extents, ", ")[[1]]), collapse = ", "))

# only keep common studies
for (this_extent in spatial_extents) {
  df <- get(paste0("study_level_data_", this_extent))
  df <- df[tolower(df$name) %in% common_studies, ]
  assign(paste0("study_level_data_", this_extent), df)
}


##### FIT PARAMS ####

for (this_extent in spatial_extents) {

  df <- get(paste0("study_level_data_", this_extent))
  
  # getting variances
  
  d_se <- function(d, n1, n2 = NULL) {
    if (is.null(n2)) { # one-sample
      se <- sqrt(1 / n1 + (d^2 / (2 * n1)))
    } else { # two-sample
      se <- sqrt((n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2))))
    }
    return(se)
  }
  
  # preallocate
  n <- numeric(nrow(df))
  d_var <- numeric(nrow(df))
  
  for (i in 1:nrow(df)) {
  if (df$orig_stat_type[[i]] == "t2") {
    if (!is.null(df$n1[i])) {
      n[i] <- df$n1[i] + df$n2[i]
      d_var[i] <- d_se(df$d[i], df$n1[i], df$n2[i])^2
    } else {
      n[i] <- df$n[i]
      d_var[i] <- NA
    }
  } else if (df$orig_stat_type[i] == "t" || df$orig_stat_type[i] == "r") {
    n[i] <- df$n[i]
    if (df$orig_stat_type[i] == "r") { # treat as 2-sample t-test
      d_var[i] <- d_se(df$d[i], df$n[i]/2, df$n[i]/2)^2
    } else { # normal 1-sample t-test
      d_var[i] <- d_se(df$d[i], df$n[i])^2
    }
  } else {
    n[i] <- NA
  }
  }
  
  # assign higher variance for small sample sizes to effectively exclude them from the shared slope fit
  # note: n>2000 also effectively removes task, which also has insufficient variance in n to contribute usefully to the fit
  # TODO: as supp
  if (use_high_sample_size_only) {
    d_var[is.na(df$n) | df$n <=2000 ] <- 15*d_var[is.na(df$n) | df$n <=2000]
  }
  

  fit_all <- rma.mv(yi = d, 
                    V = d_var,  # approximate variance
                    mods = ~ I(k2/n),
                    random = ~ 1 | overarching_category,
                    data = df,
                    method = "REML")
  
  random_effects <- ranef(fit_all)
  category_effects <- random_effects$overarching_category
  
  # Create results with category-specific intercepts
  unique_cats <- unique(df$overarching_category)
  model_params <- vector("list", length(unique_cats))
  fit_all_vb <- vector("list", length(unique_cats))
  names(model_params) <- unique_cats
  
  slope <- fit_all$beta[2]
  slope_se <- sqrt(fit_all$vb[2,2])
  
  # Only include intercept in results (not slope)
  for (cat in unique_cats) {
    cat_intercept <- fit_all$beta[1] + category_effects[cat, "intrcpt"]
    cat_intercept_se <- sqrt(fit_all$vb[1,1] + category_effects[cat, "se"]^2)
    model_params[[cat]] <- data.frame(
      est = cat_intercept,
      lwr = cat_intercept - 1.96 * cat_intercept_se,
      upr = cat_intercept + 1.96 * cat_intercept_se,
      est_se = cat_intercept_se,
      phi2_est = slope,
      phi2_lwr = slope - 1.96 * slope_se,
      phi2_upr = slope + 1.96 * slope_se,
      phi2_se = slope_se,
      row.names = paste0(cat, "_intercept")
    )
    fit_all_vb[[cat]] <- list(
      phi2_vb = fit_all$vb
    )
  }
  
  ### rename model params to model_params_<extent>
  assign(paste0("model_params_", this_extent), model_params)

}

### PLOT PARAM FITS FOR EACH SPATIAL EXTENT ###

build_predictions_mat <- function(model_params, extent = NULL) {
  ndivk2_seq <- 10^seq(log10(10), log10(20000), length.out = 100)
  X_pred <- cbind(1, 1/ndivk2_seq)
  prediction_list <- list()

  for (cat in unique(model_params$overarching_category)) {
    idx <- model_params$overarching_category == cat
    est_cat <- model_params$est[idx][1]
    phi2_cat <- model_params$phi2_est[idx][1]
    est_se_cat <- model_params$est_se[idx][1]

    preds <- est_cat + phi2_cat * (1 / ndivk2_seq)
    preds_se_fixed <- sqrt(diag(X_pred %*% fit_all_vb[[cat]]$phi2_vb %*% t(X_pred)))
    preds_se <- sqrt(preds_se_fixed^2 + est_se_cat^2)

    pred_df <- data.frame(
      X = ndivk2_seq,
      preds = preds,
      lwr = preds - 1.96 * preds_se,
      upr = preds + 1.96 * preds_se,
      overarching_category = cat
    )

    if (!is.null(extent)) {
      pred_df$extent <- extent
    }

    prediction_list[[cat]] <- pred_df
  }

  dplyr::bind_rows(prediction_list)
}

plot_param_fits <- function(this_study_level_data, predictions_mat, do_overlapping = FALSE, results_dir, plot_string = "") {
  
  invert_x <- FALSE
  if (invert_x) {
    x_pwr_factor <- -1
  } else {
    x_pwr_factor <- 1
  }
  title_txt <- "Parameter fit plot"
  fit_axis_text_size <- axis_text_size + 2
  fit_axis_title_size <- axis_title_size + 2

  if (do_overlapping) {
    predictions_mat <- dplyr::bind_rows(lapply(names(predictions_mat), function(this_extent) {
      build_predictions_mat(predictions_mat[[this_extent]], extent = this_extent)
    }))

    study_level_data_df <- dplyr::bind_rows(lapply(names(this_study_level_data), function(this_extent) {
      dat <- this_study_level_data[[this_extent]]
      dat$extent <- this_extent
      dat
    }))
    
    study_level_data_df$extent_num <- suppressWarnings(as.numeric(as.character(study_level_data_df$extent)))
    study_level_data_df$point_alpha <- scales::rescale(study_level_data_df$extent_num, to = c(0.25, 1), na.rm = TRUE)
    study_level_data_df$extent_legend <- factor(study_level_data_df$extent, levels = spatial_extents)
    predictions_mat$extent_num <- suppressWarnings(as.numeric(as.character(predictions_mat$extent)))
    predictions_mat$line_alpha <- scales::rescale(predictions_mat$extent_num, to = c(0.25, 1), na.rm = TRUE)
    predictions_mat$extent_legend <- factor(predictions_mat$extent, levels = spatial_extents)
    predictions_mat$line_group <- interaction(predictions_mat$overarching_category, predictions_mat$extent)
    
  } else {
    study_level_data_df <- this_study_level_data
    predictions_mat <- build_predictions_mat(predictions_mat)
    study_level_data_df$point_alpha <- 1
    study_level_data_df$extent_legend <- factor("all")
    predictions_mat$line_alpha <- 1
    predictions_mat$extent_legend <- factor("all")
    predictions_mat$line_group <- predictions_mat$overarching_category
  }

  p <- ggplot(data = study_level_data_df, aes(x = (n/k2)^x_pwr_factor, y = d, color = overarching_category)) +
  # p <- ggplot(data = study_level_data_df, aes(x = k²/n, y = d, color = overarching_category)) +
    geom_point(aes(alpha = extent_legend), size=0.1) +
    geom_line(
      data = predictions_mat,
      mapping = aes(x = X^x_pwr_factor, y = preds, color = overarching_category, group = line_group, alpha = extent_legend),
      linewidth = line_width_main,
      inherit.aes = FALSE
    )
  
    if(!do_overlapping) {
    p <- p +
    geom_line(
      data = predictions_mat,
      mapping = aes(x = X^x_pwr_factor, y = lwr, color = overarching_category, group = line_group, alpha = extent_legend),
      linewidth = line_width_ci,
      linetype = "dotted",
      inherit.aes = FALSE
    ) +
    geom_line(
      data = predictions_mat,
      mapping = aes(x = X^x_pwr_factor, y = upr, color = overarching_category, group = line_group, alpha = extent_legend),
      linewidth = line_width_ci,
      linetype = "dotted",
      inherit.aes = FALSE
    )
    }
    
    p <- p +
    labs(title = title_txt,
         x = "n/k²",
         y = "Effect Size (d)") +
        scale_color_manual(values = cat_colors) +
    scale_x_continuous(
      expand = c(0, 0),
      trans = "log",
      labels = function(x) format(round(x), trim = TRUE, scientific = FALSE)
    ) +
    theme_minimal() +
    # theme(legend.title = element_blank()) +
    theme(
      legend.position = "right",
      legend.justification = c("left", "center"),
      axis.text.x = element_text(size = fit_axis_text_size, angle = 0, hjust = 0.5, vjust = 1, margin = margin(t = 0)),
      axis.text.y = element_text(size = fit_axis_text_size),
      axis.title.x = element_text(size = fit_axis_title_size),
      axis.title.y = element_text(size = fit_axis_title_size),
      legend.text = element_text(size = fit_axis_text_size),
      legend.title = element_text(size = fit_axis_title_size)
    ) +
    coord_cartesian(ylim = c(-2.5, 7.5))
    # coord_cartesian(ylim = c(0.3, 1), xlim = c(0,0.05))

  if (do_overlapping) {
    alpha_values <- scales::rescale(as.numeric(spatial_extents), to = c(0.25, 1), na.rm = TRUE)
    names(alpha_values) <- spatial_extents
    p <- p + scale_alpha_manual(values = alpha_values, name = "Extent (%)")
  } else {
    p <- p + scale_alpha_manual(values = c(all = 1), guide = "none")
  }

  # show(p)
  if (!do_overlapping) {
    # make dir
    results_dir <- paste0(results_dir, "by_extent/")
    if (!dir.exists(results_dir)) {
      dir.create(results_dir)
    }
  }
  
  fname <- paste0(results_dir, "param_fit_plot_", ifelse(do_overlapping, "overlapping", "individual"), ifelse(invert_x, "_inverted", ""), plot_string, ".png")
  ggsave(p, filename = fname, width = 8, height = 5)
}

# # make individual plots
for (this_extent in spatial_extents) {
  this_study_level_data <- get(paste0("study_level_data_", this_extent))
  model_params <- get(paste0("model_params_", this_extent))
  # convert to df
  model_params <- do.call(rbind, model_params) %>%
    rownames_to_column("category") %>%
    mutate(overarching_category = sub("_intercept", "", category))
  plot_param_fits(this_study_level_data, model_params, do_overlapping = FALSE, results_dir, paste0("_", this_extent))
}


# make overlapping plots

model_params_list <- list()
study_level_data_list <- list()

for (this_extent in spatial_extents) {
  this_study_level_data <- get(paste0("study_level_data_", this_extent))
  model_params <- get(paste0("model_params_", this_extent))
  # convert to df
  model_params <- do.call(rbind, model_params) %>%
    rownames_to_column("category") %>%
    mutate(overarching_category = sub("_intercept", "", category))
  # add to list
  model_params_list[[this_extent]] <- model_params
  study_level_data_list[[this_extent]] <- this_study_level_data
}

plot_param_fits(study_level_data_list, model_params_list, do_overlapping = TRUE, results_dir)




### COMBINE RESULTS ACROSS SPATIAL EXTENTS ###

to_named_vec <- function(x) {
  if (is.data.frame(x)) {
    if (nrow(x) != 1) stop("Expected 1-row data.frame")
    v <- unlist(x[1, ], use.names = TRUE)
  } else if (is.list(x)) {
    v <- unlist(x, use.names = TRUE)
  } else {
    v <- x
  }
  setNames(as.numeric(v), names(v))
}
  
combine_category_df <- function(spatial_extents, category, prefix = "model_params_", env = .GlobalEnv) { # TODO: gotta get rid of the global inheritance
  models <- mget(paste0(prefix, spatial_extents), envir = env, inherits = TRUE)
  vecs <- lapply(models, function(m) to_named_vec(m[[category]]))
  
  rn <- Reduce(union, lapply(vecs, names))
  mat <- sapply(vecs, function(v) {
    out <- setNames(rep(NA_real_, length(rn)), rn)
    out[names(v)] <- v
    out
  }, simplify = "matrix")
  
  colnames(mat) <- spatial_extents
  as.data.frame(mat, check.names = FALSE)  # keeps "05", "10", etc.
}

model_params_master <- setNames(vector("list", length(categories)), categories)
for (cat in categories) {
  model_params_master[[cat]] <- combine_category_df(spatial_extents, cat)
}
  
  
### PLOT EXPECTED R2 ACROSS SPATIAL EXTENTS ###
  
# Individual plots
  
plot_model_params__individual <- function(df, category, do_r2) {
  
  df <- as.data.frame(t(df[[category]]))
  extent_labels <- setNames(as.character(as.numeric(rownames(df))), rownames(df))
  
  
  # TODO: figure out how want to deal with one-sample r2 - if (grep(category, "task")) {
  if (do_r2) { # convert d to r2
    d2r2 <- function(d) { ifelse(d < 0, 0, d^2 / (d^2 + 4)) } # standard conversion a la hauselin
    df <- df %>%
      mutate(est = d2r2(est), lwr = d2r2(lwr), upr = d2r2(upr))
    esz_str <- "R²"
  }
  else {esz_str <- "d"}
  
  p <- ggplot(df, aes(x=rownames(df), y = est)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, linewidth = line_width_ci) +
    xlim(rownames(df)) +
    scale_x_discrete(labels = extent_labels) +
    ylim(0, ifelse(do_r2, 1, 4)) +
    labs(title = paste0("Model Parameters for ", category), x = "extent of brain included", y = paste0("Effect size (",esz_str,")")) +
    # add a line that's the max across est
    geom_hline(yintercept = max(df$est, na.rm = TRUE), linetype = "dashed", color = "red", linewidth = line_width_ref) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      axis.title.x = element_text(size = axis_title_size),
      axis.title.y = element_text(size = axis_title_size)
    )

  vline_positions <- which(rownames(df) %in% c("05", "50"))
  if (length(vline_positions) > 0) {
    p <- p + geom_vline(xintercept = vline_positions, linetype = "dashed", color = "grey60", linewidth = 2*line_width_ref)
  }
  
  #make dir
  if (!dir.exists(paste0(results_dir, "/cat"))) {
    dir.create(paste0(results_dir, "/cat"))
  }
  ggsave(filename = paste0(results_dir, "/cat/model_params_individual_", category, "_", ifelse(do_r2, "r2", "d"), ".png"), plot = p, width = 5, height = 4)
}

# make individual plots
# for (cat in categories) {
#   plot_model_params__individual(model_params_master, cat, do_r2)
# }


## Overlapping plots

plot_model_params__overlapping <- function(df, do_r2, results_dir) {
  
  # params
  alpha <- 0.7
  main_title <- paste0("Model Parameters")
  x_label <- "extent of brain included"
  x_limits <- c(0, 100)
  categories <- c("psychological", "physical", "task activation", "task connectivity")
  categories <- factor(categories, levels = categories) # fix order
  this_fn <- "~/Desktop/test_plot.png"

  
  # TODO: figure out how want to deal with one-sample r2 - if (grep(category, "task")) {
  if (do_r2) { # convert d to r2
    d2r2 <- function(d) { ifelse(d < 0, 0, d^2 / (d^2 + 4)) } # standard conversion a la hauselin
    esz_str <- "R²"
  } else {esz_str <- "d"}
  y_label <- paste0("Effect size (",esz_str,")")
  
  all_df <- dplyr::bind_rows(lapply(as.character(categories), function(cat) {
    df_cat <- as.data.frame(t(df[[cat]]))
    df_cat$category <- cat
    df_cat$extent <- suppressWarnings(as.numeric(rownames(df_cat)))
    df_cat
  }))
  
  all_df <- all_df[!is.na(all_df$extent), , drop = FALSE]
  
  if (do_r2) {
    all_df$est <- d2r2(all_df$est)
    all_df$lwr <- d2r2(all_df$lwr)
    all_df$upr <- d2r2(all_df$upr)
  }
  
  
  
  p <- ggplot(all_df, aes(x = extent, y = est, color = category)) +
    geom_line(aes(x = extent, y = est, color = category), linewidth = line_width_main, alpha = alpha) +
    geom_line(aes(x = extent, y = lwr, color = category), linetype = "dotted", linewidth = line_width_ci, alpha = alpha) +
    geom_line(aes(x = extent, y = upr, color = category), linetype = "dotted", linewidth = line_width_ci, alpha = alpha) +
    geom_vline(xintercept = c(5, 50), linetype = "dashed", color = "grey60", linewidth = 2 * line_width_ref, alpha = 0.5) +
    scale_color_manual(values = cat_colors) +
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
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      axis.title.x = element_text(size = axis_title_size),
      axis.title.y = element_text(size = axis_title_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # scale_x_continuous(expand = c(0, 0), trans = "log") +
    scale_x_continuous(
      breaks = as.numeric(spatial_extents),
      labels = as.character(as.numeric(spatial_extents))
    ) +
    coord_cartesian(xlim = x_limits) +
    coord_cartesian(ylim = c(0, ifelse(do_r2, 1, 4))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = line_width_ref) +
    # make legend small
    guides(color = guide_legend(override.aes = list(size = 0.5))) +
    # theme(legend.title = element_blank())
    theme(legend.position = "none")
  
  # show(p)
  # ggsave(this_fn, plot = p, width = 5, height = 4)
  ggsave(p, filename = paste0(results_dir, "model_params_overlapping_", ifelse(do_r2, "r2", "d"), ".png"), width = 6, height = 4)
}

# also save difference between 5% and 50% as table
log_differences <- function(df, do_r2, results_dir) {
  all_df <- dplyr::bind_rows(lapply(as.character(categories), function(cat) {
    df_cat <- as.data.frame(t(df[[cat]]))
    df_cat$category <- cat
    df_cat$extent <- suppressWarnings(as.numeric(rownames(df_cat)))
    df_cat
  }))
  
  diff_df <- all_df %>%
    filter(extent %in% c(5, 50)) %>%
    select(category, extent, est) %>%
    pivot_wider(names_from = extent, values_from = est) %>%
    mutate(diff = `50` - `5`) %>%
    mutate(percent_increase = (`50` - `5`) / abs(`5`) * 100)
  
  if (do_r2) {
    d2r2 <- function(d) { ifelse(d < 0, 0, d^2 / (d^2 + 4)) } # standard conversion a la hauselin
    diff_df <- diff_df %>%
      mutate(`5` = d2r2(`5`), `50` = d2r2(`50`), diff = d2r2(diff))
  }
  
  write.csv(diff_df, file = paste0(results_dir, "model_params_diff_50_minus_05_", ifelse(do_r2, "r2", "d"), ".csv"), row.names = FALSE)
}
    
# make overlapping plot
plot_model_params__overlapping(model_params_master, do_r2, results_dir)
log_differences(model_params_master, do_r2, results_dir)


## Automated merge workflow: run both motion types and combine

compute_motion_results <- function(motion_type_local,
                                   pooling_type,
                                   data_dir,
                                   spatial_extents,
                                   categories,
                                   use_high_sample_size_only = FALSE) {

  extract_metric <- function(dat, metric, pooling_type, motion_type) {
    unlist(lapply(dat, function(x) {
      x1 <- x[grepl("multi", names(x))]
      x1 <- x1[grepl(paste0("pooling\\.", pooling_type), names(x1))]
      lapply(x1[grepl(paste0("motion\\.", motion_type), names(x1))], `[[`, metric)
    }), recursive = FALSE)
  }

  d_se <- function(d, n1, n2 = NULL) {
    if (is.null(n2)) {
      sqrt(1 / n1 + (d^2 / (2 * n1)))
    } else {
      sqrt((n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2))))
    }
  }

  study_level_data_list_local <- list()

  for (this_extent in spatial_extents) {
    v_filename <- list.files(
      data_dir,
      pattern = paste0("braineffex_data_", this_extent, "_percent_.*\\.RData$"),
      full.names = TRUE
    )
    load(v_filename)

    d_list <- extract_metric(v$data, "d", pooling_type, motion_type_local)
    n_list <- extract_metric(v$data, "n", pooling_type, motion_type_local)
    n1_list <- extract_metric(v$data, "n1", pooling_type, motion_type_local)
    n2_list <- extract_metric(v$data, "n2", pooling_type, motion_type_local)

    study_name <- sub("\\..*$", "", names(d_list))
    df <- v$study
    df <- df[match(study_name, tolower(df$name)), ]

    df$d <- unlist(d_list[match(tolower(df$name), study_name)])
    df$n <- unlist(n_list[match(tolower(df$name), study_name)])
    df$n1 <- unlist(n1_list[match(tolower(df$name), study_name)])
    df$n2 <- unlist(n2_list[match(tolower(df$name), study_name)])
    df$k2 <- ifelse(df$orig_stat_type == "t", 1, 4)

    df <- df[!grepl("test", df$name), ]
    df <- df[!duplicated(df$name), ]

    df$category[df$category == "clinical"] <- "psychiatric"
    df$category[df$test_component_2 == "bmi"] <- "biometric"
    df$category[grepl("sex", tolower(df$test_component_2))] <- "sex (demographic)"
    df$category[grepl("gender", tolower(df$test_component_2))] <- "sex (demographic)"
    df$category[df$test_component_2 == "cbcl_scr_syn_internal_t_FU1"] <- "psychiatric"
    df$category[grep("age", df$name)] <- "age (demographic)"
    df$category[df$category == "cognitive" & df$orig_stat_type == "t"] <- "cognitive (task)"

    need_to_flip <- which(
      grepl("abcd", tolower(df$dataset)) &
      df$map_type == "fc" &
      df$orig_stat_type == "t" &
      tolower(df$test_component_2) != "rest"
    )
    if (length(need_to_flip) > 0) {
      tmp <- df$test_component_1[need_to_flip]
      df$test_component_1[need_to_flip] <- df$test_component_2[need_to_flip]
      df$test_component_2[need_to_flip] <- tmp
      df$name[need_to_flip] <- sapply(need_to_flip, function(i) {
        new_name <- gsub(df$test_component_1[i], "tmp", df$name[i])
        new_name <- gsub(df$test_component_2[i], df$test_component_1[i], new_name)
        gsub("tmp", df$test_component_2[i], new_name)
      })
      df$basefile[need_to_flip] <- sapply(need_to_flip, function(i) {
        new_basefile <- gsub(df$test_component_1[i], "tmp", df$basefile[i])
        new_basefile <- gsub(df$test_component_2[i], df$test_component_1[i], new_basefile)
        gsub("tmp", df$test_component_2[i], new_basefile)
      })
    }

    df <- df %>%
      mutate(overarching_category = case_when(
        category %in% c("biometric", "sex (demographic)", "age (demographic)") ~ "physical",
        category %in% c("cognitive", "psychiatric") ~ "psychological",
        category == "cognitive (task)" & !grepl("act", map_type) ~ "task connectivity",
        category == "cognitive (task)" & grepl("act", map_type) ~ "task activation",
        TRUE ~ "other"
      ))

    df$overarching_category <- as.factor(df$overarching_category)
    study_level_data_list_local[[this_extent]] <- df
  }

  common_studies_local <- Reduce(
    intersect,
    lapply(spatial_extents, function(this_extent) tolower(study_level_data_list_local[[this_extent]]$name))
  )

  for (this_extent in spatial_extents) {
    df <- study_level_data_list_local[[this_extent]]
    study_level_data_list_local[[this_extent]] <- df[tolower(df$name) %in% common_studies_local, ]
  }

  model_params_by_extent <- list()
  model_params_list_local <- list()

  for (this_extent in spatial_extents) {
    df <- study_level_data_list_local[[this_extent]]

    n <- numeric(nrow(df))
    d_var <- numeric(nrow(df))

    for (i in seq_len(nrow(df))) {
      if (df$orig_stat_type[[i]] == "t2") {
        if (!is.null(df$n1[i])) {
          n[i] <- df$n1[i] + df$n2[i]
          d_var[i] <- d_se(df$d[i], df$n1[i], df$n2[i])^2
        } else {
          n[i] <- df$n[i]
          d_var[i] <- NA
        }
      } else if (df$orig_stat_type[i] == "t" || df$orig_stat_type[i] == "r") {
        n[i] <- df$n[i]
        if (df$orig_stat_type[i] == "r") {
          d_var[i] <- d_se(df$d[i], df$n[i] / 2, df$n[i] / 2)^2
        } else {
          d_var[i] <- d_se(df$d[i], df$n[i])^2
        }
      } else {
        n[i] <- NA
      }
    }

    if (use_high_sample_size_only) {
      d_var[is.na(df$n) | df$n <= 2000] <- 15 * d_var[is.na(df$n) | df$n <= 2000]
    }

    fit_all <- rma.mv(
      yi = d,
      V = d_var,
      mods = ~ I(k2/n),
      random = ~ 1 | overarching_category,
      data = df,
      method = "REML"
    )

    random_effects <- ranef(fit_all)
    category_effects <- random_effects$overarching_category
    unique_cats <- unique(df$overarching_category)

    model_params <- vector("list", length(unique_cats))
    names(model_params) <- unique_cats

    slope <- fit_all$beta[2]
    slope_se <- sqrt(fit_all$vb[2, 2])

    for (cat in unique_cats) {
      cat_intercept <- fit_all$beta[1] + category_effects[cat, "intrcpt"]
      cat_intercept_se <- sqrt(fit_all$vb[1, 1] + category_effects[cat, "se"]^2)
      model_params[[cat]] <- data.frame(
        est = cat_intercept,
        lwr = cat_intercept - 1.96 * cat_intercept_se,
        upr = cat_intercept + 1.96 * cat_intercept_se,
        est_se = cat_intercept_se,
        phi2_est = slope,
        phi2_lwr = slope - 1.96 * slope_se,
        phi2_upr = slope + 1.96 * slope_se,
        phi2_se = slope_se,
        row.names = paste0(cat, "_intercept")
      )
    }

    model_params_by_extent[[this_extent]] <- model_params
    model_params_list_local[[this_extent]] <- do.call(rbind, model_params) %>%
      rownames_to_column("category") %>%
      mutate(overarching_category = sub("_intercept", "", category))
  }

  model_params_master_local <- setNames(vector("list", length(categories)), categories)
  names(model_params_master_local) <- categories

  for (cat in categories) {
    vecs <- lapply(model_params_by_extent, function(m) {
      x <- m[[cat]]
      if (is.null(x)) return(setNames(numeric(0), character(0)))
      setNames(as.numeric(unlist(x[1, ], use.names = TRUE)), names(unlist(x[1, ], use.names = TRUE)))
    })

    rn <- Reduce(union, lapply(vecs, names))
    mat <- sapply(vecs, function(v) {
      out <- setNames(rep(NA_real_, length(rn)), rn)
      if (length(v) > 0) out[names(v)] <- v
      out
    }, simplify = "matrix")

    colnames(mat) <- spatial_extents
    model_params_master_local[[cat]] <- as.data.frame(mat, check.names = FALSE)
  }

  list(
    study_level_data_list = study_level_data_list_local,
    model_params_list = model_params_list_local,
    model_params_master = model_params_master_local
  )
}

motion_results_reg <- compute_motion_results(
  motion_type_local = "regression",
  pooling_type = pooling_type,
  data_dir = data_dir,
  spatial_extents = spatial_extents,
  categories = categories,
  use_high_sample_size_only = use_high_sample_size_only
)

motion_results_threshold <- compute_motion_results(
  motion_type_local = "threshold",
  pooling_type = pooling_type,
  data_dir = data_dir,
  spatial_extents = spatial_extents,
  categories = categories,
  use_high_sample_size_only = use_high_sample_size_only
)

study_level_data_list_reg <- motion_results_reg$study_level_data_list
model_params_list_reg <- motion_results_reg$model_params_list
model_params_master_reg <- motion_results_reg$model_params_master

study_level_data_list_threshold <- motion_results_threshold$study_level_data_list
model_params_list_threshold <- motion_results_threshold$model_params_list
model_params_master_threshold <- motion_results_threshold$model_params_master

study_level_data_list_merge <- study_level_data_list_reg
for (this_extent in spatial_extents) {
  study_level_data_list_merge[[this_extent]] <- study_level_data_list_merge[[this_extent]][
    !study_level_data_list_merge[[this_extent]]$overarching_category %in% c("task activation", "task connectivity"),
  ]

  thresh_data <- study_level_data_list_threshold[[this_extent]]
  task_data <- thresh_data[thresh_data$overarching_category %in% c("task activation", "task connectivity"), ]
  study_level_data_list_merge[[this_extent]] <- rbind(study_level_data_list_merge[[this_extent]], task_data)
}

model_params_list_merge <- model_params_list_reg
for (this_extent in spatial_extents) {
  model_params_list_merge[[this_extent]] <- model_params_list_merge[[this_extent]][
    !model_params_list_merge[[this_extent]]$overarching_category %in% c("task activation", "task connectivity"),
  ]

  thresh_params <- model_params_list_threshold[[this_extent]]
  task_params <- thresh_params[thresh_params$overarching_category %in% c("task activation", "task connectivity"), ]
  model_params_list_merge[[this_extent]] <- rbind(model_params_list_merge[[this_extent]], task_params)
}

model_params_master_merge <- model_params_master_reg
model_params_master_merge$`task activation` <- model_params_master_threshold$`task activation`
model_params_master_merge$`task connectivity` <- model_params_master_threshold$`task connectivity`

results_dir__merge <- paste0(results_dir_master, "pooling.", pooling_type, ".motion.regression_threshold_merge/")
if (!dir.exists(results_dir__merge)) {
  dir.create(results_dir__merge, recursive = TRUE)
}

plot_param_fits(study_level_data_list_merge, model_params_list_merge, do_overlapping = TRUE, results_dir__merge)
plot_model_params__overlapping(model_params_master_merge, do_r2, results_dir__merge)
log_differences(model_params_master_merge, do_r2, results_dir__merge)
save(
  study_level_data_list_reg,
  model_params_list_reg,
  model_params_master_reg,
  study_level_data_list_threshold,
  model_params_list_threshold,
  model_params_master_threshold,
  study_level_data_list_merge,
  model_params_list_merge,
  model_params_master_merge,
  file = paste0(results_dir__merge, "merged_results.RData")
)
