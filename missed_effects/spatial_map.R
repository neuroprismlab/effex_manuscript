### Set params & filenames
motion_type <- "regression"
multi_string <- "mv\\.none" # "multi" or "mv\\.none"
pooling_type <- "none"
spatial_extents <- c("05", "10", "25", "50", "75", "100")

data_dir <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My\ Drive/Lab/xMore/Software/scripts/R/myscripts/effect_size/BrainEffeX_utils/inst/meta/'
results_dir <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-s.noble@northeastern.edu/My\ Drive/Lab/Tasks-Ongoing/-K99/Effect_Size/manuscript/figures/plots/crossbrain_effects__spatial_extent/'

# rename results_dir based on motion and pooling
results_dir <- paste0(results_dir, "pooling.", pooling_type, ".motion.", motion_type, "/spatial_map/")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# categories + plotting params
categories <- c("physical", "psychological", "task activation", "task connectivity")
cat_colors <- RColorBrewer::brewer.pal(length(categories), "Set1")
color_map <- setNames(RColorBrewer::brewer.pal(length(categories), "Set1"), categories)

if (load_data) {
  v_filename <- list.files(data_dir, pattern = paste0("braineffex_data_100_percent_.*\\.RData$"), full.names = TRUE)
  
  # load data
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
  
}

# plot different top percents for each study
for (this_extent in spatial_extents) {

  p <- ggplot(df, aes(x = overarching_category, y = d, color = overarching_category)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    scale_color_manual(values = color_map) +
    labs(title = paste0("Effect Sizes for ", this_extent, "% Spatial Extent"),
         x = "Overarching Category",
         y = "Effect Size (d)") +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(filename = paste0(results_dir, "spatial_map_", this_extent, ".png"), plot = p, width = 8, height = 6)
}

