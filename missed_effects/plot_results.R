plot_results <- function(
  all_target_effect_type__from_basis,
  all_outcome_categories,
  cat_colors,
  out_master_dir,
  text_size = 26,
  ticks_size = 24,
  transparency_main = 0.6,
  transparency_overlay = 0.5
) {

  axis_title_size <- text_size
  axis_tick_size <- max(ticks_size, text_size * 0.9)

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("plot_results_multicat requires the gridExtra package.")
  }

  build_summary_multicat_plot <- function(target_effect_type, metric, metric_label, y_limits, hline_yintercept = NULL, vline_xintercept = 3) {
    plot_list <- lapply(seq_along(all_outcome_categories), function(category_idx) {
      this_category <- all_outcome_categories[category_idx]
      summary_path <- paste0(out_master_dir, this_category, "/", target_effect_type, "_effect/summary_results.rds")
      if (!file.exists(summary_path)) {
        return(NULL)
      }

      this_summary <- readRDS(summary_path)
      cat_color <- unname(cat_colors[this_category])

      str_mean <- metric
      if (grepl("_mean$", metric)) {
        str_sd <- sub("_mean", "_sd", metric)
      } else {
        str_sd <- NULL
      }

      df_plot <- data.frame(
        sample_size = this_summary$sample_size,
        ub = NA_real_,
        lb = NA_real_,
        point = NA_real_
      )

      padding <- 2
      df_plot$point <- this_summary[, str_mean]
      df_plot$point[df_plot$point > (y_limits[2] + padding)] <- y_limits[2] + padding

      if (!is.null(str_sd) && str_sd %in% names(this_summary)) {
        df_plot$ub <- this_summary[, str_mean] + this_summary[, str_sd]
        df_plot$lb <- this_summary[, str_mean] - this_summary[, str_sd]
        df_plot$ub[df_plot$ub > (y_limits[2] + padding)] <- y_limits[2] + padding
        df_plot$lb[df_plot$lb > (y_limits[2] + padding)] <- y_limits[2] + padding
        df_plot$lb[df_plot$lb < (y_limits[1] - padding)] <- y_limits[1] - padding
      }

      df_plot$x_num <- as.integer(factor(df_plot$sample_size))
      if (grepl("Corr", metric_label)) {
        x_label <- "Sample Size"
      } else {
        x_label <- "Sample Size of Basis Study"
      }

      y_axis_title <- if (category_idx == 1) metric_label else NULL
      y_axis_title_element <- if (category_idx == 1) {
        element_text(size = axis_title_size, color = "black")
      } else {
        element_blank()
      }

      p <- ggplot(df_plot, aes(x = x_num, y = point)) +
        geom_line(linewidth = 0.5, color = cat_color) +
        scale_x_continuous(breaks = df_plot$x_num, labels = df_plot$sample_size)

      if (!is.null(hline_yintercept)) {
        p <- p + geom_hline(
          yintercept = hline_yintercept,
          colour = "grey50",
          linetype = "dashed",
          linewidth = 2,
          alpha = 0.75
        )
      }

      if (!is.null(vline_xintercept)) {
        p <- p + geom_vline(
          xintercept = vline_xintercept,
          colour = "grey50",
          linetype = "dotted",
          linewidth = 3,
          alpha = 0.5
        )
      }

      if (!is.null(str_sd) && str_sd %in% names(this_summary)) {
        p <- p + geom_ribbon(
          aes(x = x_num, ymin = lb, ymax = ub),
          alpha = transparency_main,
          fill = cat_color,
          inherit.aes = FALSE
        )
      }

      p +
        coord_cartesian(ylim = y_limits) +
        labs(title = NULL, x = x_label, y = y_axis_title) +
        theme_classic(base_size = axis_tick_size) +
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = axis_title_size),
          axis.title.y.left = y_axis_title_element,
          axis.text.x = element_text(size = axis_tick_size, color = "black"),
          axis.text.y.left = element_text(size = axis_tick_size, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          plot.margin = margin(t = 5, r = 2, b = 5, l = 2),
          plot.title = element_blank()
        )
    })

    plot_list <- Filter(Negate(is.null), plot_list)

    if (length(plot_list) == 0) {
      message(paste0("No saved summary files found for target effect type: ", target_effect_type, ". Skipping multicategory summary plot for ", metric, "."))
      return(invisible(NULL))
    }

    ncol_used <- length(plot_list)
    combined_plot <- gridExtra::arrangeGrob(grobs = plot_list, ncol = ncol_used)
    out_file <- paste0(out_master_dir, metric, "_multicat_", target_effect_type, ".png")
    ggsave(filename = out_file, plot = combined_plot, width = 5.5 * ncol_used, height = 4)
  }

  build_overlay_multicat_plot <- function(target_effect_type) {
    plot_list <- lapply(seq_along(all_outcome_categories), function(category_idx) {
      this_category <- all_outcome_categories[category_idx]
      summary_path <- paste0(out_master_dir, this_category, "/", target_effect_type, "_effect/summary_results.rds")
      if (!file.exists(summary_path)) {
        return(NULL)
      }

      this_summary <- readRDS(summary_path)
      all_sample_sizes <- sort(unique(this_summary$sample_size))
      this_summary$x_num <- match(this_summary$sample_size, all_sample_sizes)

      sec_x_breaks <- seq_along(all_sample_sizes)
      sec_x_labels <- round(this_summary$expected_n__replication_mean[match(all_sample_sizes, this_summary$sample_size)])

      left_max <- 1.5

      left_color <- unname(cat_colors[this_category])
      this_hsv <- rgb2hsv(col2rgb(left_color))
      right_color <- hsv(h = this_hsv[1, 1], s = this_hsv[2, 1] * 0.3, v = this_hsv[3, 1] * 0.6)
      y_axis_title <- if (category_idx == 1) "Proportion Detected" else NULL
      y_axis_title_element <- if (category_idx == 1) {
        element_text(size = axis_title_size, color = "black")
      } else {
        element_blank()
      }

      ggplot(this_summary, aes(x = x_num)) +
        geom_ribbon(
          aes(
            x = x_num,
            ymin = expect_v_actual_n_tp__based_on_basis_mean - expect_v_actual_n_tp__based_on_basis_sd,
            ymax = expect_v_actual_n_tp__based_on_basis_mean + expect_v_actual_n_tp__based_on_basis_sd
          ),
          fill = left_color,
          alpha = transparency_main,
          colour = NA,
          inherit.aes = FALSE
        ) +
        geom_ribbon(
          aes(
            x = x_num,
            ymin = overlap_mean - overlap_sd,
            ymax = overlap_mean + overlap_sd
          ),
          fill = right_color,
          alpha = transparency_overlay,
          colour = NA,
          inherit.aes = FALSE
        ) +
        geom_hline(yintercept = 1, colour = "grey50", linetype = "dashed", linewidth = 2, alpha = 0.75) +
        geom_vline(xintercept = 3, colour = "grey50", linetype = "dotted", linewidth = 3, alpha = 0.5) +
        geom_line(aes(y = expect_v_actual_n_tp__based_on_basis_mean), linewidth = 0.7, colour = left_color) +
        geom_line(aes(y = overlap_mean), linewidth = 0.7, colour = right_color) +
        scale_y_continuous(name = y_axis_title) +
        coord_cartesian(ylim = c(0, left_max)) +
        scale_x_continuous(
          name = "Planned Sample Size (Main Study)",
          breaks = sec_x_breaks,
          labels = sec_x_labels,
          sec.axis = sec_axis(
            transform = ~ .,
            breaks = sec_x_breaks,
            labels = all_sample_sizes,
            name = "Sample Size of Basis Study"
          )
        ) +
        theme_classic(base_size = axis_tick_size) +
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = axis_title_size),
          axis.title.y.left = y_axis_title_element,
          axis.title.x.top = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_tick_size, color = "black"),
          axis.text.y.left = element_text(size = axis_tick_size, color = "black"),
          axis.text.x.top = element_text(size = axis_tick_size, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_blank()
        )
    })

    plot_list <- Filter(Negate(is.null), plot_list)

    if (length(plot_list) == 0) {
      message(paste0("No saved summary files found for target effect type: ", target_effect_type, ". Skipping multi-category overlay plot."))
      return(invisible(NULL))
    }

    ncol_used <- length(plot_list)
    combined_plot <- gridExtra::arrangeGrob(grobs = plot_list, ncol = ncol_used)
    out_file <- paste0(out_master_dir, "overlay_tp_and_overlap_multicat_", target_effect_type, ".png")
    ggsave(filename = out_file, plot = combined_plot, width = 5.5 * ncol_used, height = 4.75)
  }

  for (target_effect_type__from_basis in all_target_effect_type__from_basis) {
    build_overlay_multicat_plot(target_effect_type__from_basis)

    build_summary_multicat_plot(
      target_effect_type = target_effect_type__from_basis,
      metric = "type_m_mean",
      metric_label = "Mean Type M Error",
      y_limits = c(0, 10),
      hline_yintercept = 1
    )
    build_summary_multicat_plot(
      target_effect_type = target_effect_type__from_basis,
      metric = "type_s_mean",
      metric_label = "Mean Type S Error",
      y_limits = c(0, 0.5),
      hline_yintercept = 0
    )
    build_summary_multicat_plot(
      target_effect_type = target_effect_type__from_basis,
      metric = "expect_v_actual_n_tp__based_on_corr_mean",
      metric_label = "Proportion Expected TPs Detected",
      y_limits = c(0, 1.5),
      hline_yintercept = 1
    )
  }

  invisible(NULL)
}
