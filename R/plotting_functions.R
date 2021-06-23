# Plotting functions for HGIS data


#' Plot of HGIS data vs time
#' 
#' Produces a grid of per-sample plots of parameters vs time. An 
#' error term for error bars may be included if available.
#' 
#' @param data A tibble in format of output from get_hgis_data().
#' @param y_var A column to plot against time
#' @param errors A column containing errors for y_var
#'
#' @return A ggplot object with norm ratio and current vs time,
#' faceted for each pos_name.
#' 
#' @export
#'
plot_hgis_time <- function(data, y_var = normFm, errors = NULL, outlier) {
  y_var <- enquo(y_var)
  outlier <- enquo(outlier)
  
  if (!("pos_name" %in% names(data))) {
    data <- data %>% 
      group_by(sample_name, pos) %>% 
      mutate(pos_name = paste(pos, sample_name, sep = " - "))
  }
  
  if (!("cum_acqtime" %in% names(data))) {
    data <- data %>% 
      group_by(sample_name, pos) %>% 
      mutate(cum_acqtime = cumsum(cycles) / 10)
  }
  
  if (missing(errors)) {
    ggplot(data, aes(cum_acqtime, !!y_var)) +
      geom_point(aes(color = !!outlier)) + 
      facet_wrap(vars(pos_name), scales = "free")
  } else {
    errors <- enquo(errors)
    
    ggplot(data, aes(cum_acqtime, !!y_var)) +
      geom_pointrange(aes(ymin = !!y_var - !!errors, ymax = !!y_var + !!errors,
                         color = !!outlier), size = .2) + 
      facet_wrap(vars(pos_name), scales = "free") +
      labs(title = paste(as_label(y_var), "by target"),
           x = "Time (s)",
           y = as_label(y_var))
  }
}


#' Boxplot of samples in HGIS test, colored by dilution
#' 
#' Normalized ratios are per-run, not per-sample. Use `plot_hgis_summary()`
#' for plotting normalized and blank corrected data
#' 
#' hline for modern is mean value of tank standard rec 101730
#'
#' @param data A tibble in format of output from get_hgis_data().
#'
#' @return A ggplot object with Boxplots of data.
#' @export
#'
plot_hgis <- function(data) {
  
  data %>% 
    mutate(Fm = ifelse(normFm > .3, "modern", "dead"),
           Name = factor(pos_name, levels = unique(pos_name[order(Fm, dil_factor)])),
           Fm = ordered(Fm, levels = c("modern", "dead")),
           Dillution_Ratio = ordered(dil_factor, levels = c(0, 1, 3))) %>%  
    ggplot(aes(Name, normFm, fill = Dillution_Ratio)) +
      geom_boxplot() +
      geom_hline(data = data.frame(yint=1.0377,Fm=ordered("modern", levels = c("modern", "dead"))),
                 aes(yintercept = yint)) + 
      geom_hline(data = data.frame(yint=0,Fm=ordered("dead", levels = c("modern", "dead"))), 
                 aes(yintercept = yint)) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      scale_fill_manual(values = c("0" = "slategray1", "1" = "slategray2", "3" = "slategray3")) +
      labs(x = NULL,
           y = "Fraction Modern") +
      facet_grid(Fm~., scales = "free_y")
}


#' Plot normalized and corrected HGIS summary data
#'
#' @param data 
#'
#' @return A ggplot2 object.
#' @export
#'
plot_hgis_summary <- function(data) {
  ggplot(data, aes(sample_name, fm_corr)) +
    geom_pointrange(aes(ymin = fm_corr - sig_fm_corr, ymax = fm_corr + sig_fm_corr),
                    size = 0.1) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      labs(x = NULL,
           y = "Fraction Modern")
    
}

#' Plot difference from consensus
#'
#' @param data A data frame as output from `compare_consensus()`
#'
#' @return A ggplot2 object
#' @export
#'
plot_hgis_consensus <- function(data) {
  ggplot(data, aes(fm_consensus, fm_diff, color = he12C)) +
    geom_hline(yintercept = 0) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_pointrange(aes(ymin = fm_diff - sig_fm_corr, ymax = fm_diff + sig_fm_corr), 
                    position = position_dodge2(width = 0.1),
                    size = .2) +
    labs(subtitle = "Blank corrected",
         x = "Fm expected",
         y = "Fm difference") +
    theme_classic()
}