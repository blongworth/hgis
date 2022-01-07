# Functions for HGIS data reduction

#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter mutate select summarize group_by ungroup n case_when across left_join
#' @importFrom stringr str_starts
NULL

#' Calculate d13C
#'
#' @param he1312 Normalized 13/12 ratio
#'
#' @return Delta 13C of the sample.
#' @export
#'
calc_d13c <- function(he1312) {
	1000 * (he1312 / 1.12372 -1)
}

#' Flag outliers using a list of measurements
#'
#' @param data A data frame of run data, processed by `process_hgis_results`.
#' @param outliers A data frame of pos and meas of outliers 
#'
#' @return The data frame with outliers flagged as FALSE in ok_calc.
#' @export
#'
flag_outliers <- function(data, outliers) {
  outliers <- mutate(outliers, ok_calc = FALSE)
  data %>% 
    select(-ok_calc) %>% 
    left_join(outliers, by = c("pos", "meas")) %>% 
    mutate(ok_calc = is.na(ok_calc))
}

#' Summarize HGIS data per target
#' 
#' Combine runs of targets and compute per-target statistics. 
#'
#' @param data A data frame of run data, processed by `process_hgis_results`.
#' @param remove_outliers If `TRUE`, do not include samples where `ok_calc == FALSE`.
#' @param get_consensus If `TRUE`, get and add consensus values from NOSAMS DB.
#'
#' @return A data frame of per-sample data.
#' @export
#'
sum_hgis_targets <- function(data, remove_outliers = TRUE, get_consensus = TRUE) {
  if (remove_outliers == TRUE) {
    data <- filter(data, ok_calc)
  }
  data_sum <- data %>% 
    group_by(pos, sample_name, sample_type, wheel) %>% 
    summarize(runtime = min(runtime),
              counts = sum(CntTotGT),
              n_runs = n(),
              ext_err = amstools::se(corr_14_12),
              #int_err = 1/sqrt(counts) * mean(corr_14_12),
              # now using error propagated mean
              int_err = sqrt(sum((corr_14_12/sqrt(CntTotGT))^2))/n(),
              max_err = pmax(ext_err, int_err),
              across(c(le12C, le13C, he12C, he13C,
                       he13_12, he14_12, corr_14_12, norm_del13c), 
                     list(mean = mean, sd = sd),
                     .names = "{ifelse(.fn == 'mean', '', 'sig_')}{.col}"),
              ) %>% 
    mutate(rec_num = case_when(str_starts(sample_name, "LiveGas") ~ 101730,
                               str_starts(sample_name, "C-?1") ~ 83028,
                               str_starts(sample_name, "TIRI-F") ~ 2138,
                               str_starts(sample_name, "TIRI-I") ~ 17185,
                               str_starts(sample_name, "C-?2") ~ 1082,
                               str_starts(sample_name, "NOSAMS-?2") ~ 38809,
                               str_starts(sample_name, "DeadGas") ~ 72446)) %>% 
    ungroup()
  if (get_consensus == TRUE) {
    # TODO: fail gracefully if no DB available or use local file
    std <- amstools::getStdTable()
    data_sum %>%  
      left_join(select(std, rec_num, fm_consensus), by = "rec_num") %>% 
      mutate(fm_consensus = case_when(rec_num == 101730 ~ 1.0398,
                                      rec_num == 72446 ~ 0.0013,
                                      TRUE ~ fm_consensus))
  } else {
    data_sum
  }
}


## Normalize

#' Normalize a gas target run using the measured and consensus value of a standard
#'
#' @param sample Ratio to be normalized.
#' @param standard Measured ratio of the standard.
#' @param stdrat Consensus value of the standard.
#'
#' @return The normalized ratio of the sample.
#' @export
#'
norm_gas <- function(sample, standard, stdrat = 1.0398) {
  sample/standard * stdrat
}

#' Calculate normalization error
#'
#' @param sample The normalized Fm of the sample
#' @param standard The normalized Fm of the standard
#' @param sample_err The error of the sample
#' @param standard_err The error of the standard
#' @param stdrat The known ratio of the standard
#' @param stdrat_err The error of the known ratio of standard
#'
#' @return The propagated error of the normalized sample
#' @export
#'
norm_err <- function(sample, standard, sample_err, standard_err, 
                     stdrat = 1.0398, stdrat_err = 0.0006) {
  sqrt(sample_err^2/sample^2 + 
       standard_err^2/standard^2 + 
       stdrat_err^2/stdrat^2)
}

#' Propagate errors
#' 
#' Add errors in quadrature for a vector of errors
#'
#' @param err A vector of errors
#'
#' @return A propagated error for the vector.
#'
prop_err <- function(err) {
  sqrt(sum(err^2))/length(err)
}
  
#' Summarize HGIS standards
#'
#' @param data An HGIS results dataframe, either raw or summarized by target
#'
#' @return A summary of standard means and errors
#'
summarize_standards <- function(data) {
  data %>% 
    filter(sample_type == "S") %>% 
    mutate(max_err = ifelse("max_err" %in% names(.), max_err, sig_14_12)) %>% 
    summarize(across(corr_14_12, list(mean = mean, se = amstools::se)),
              propagated_err = prop_err(max_err),
              norm_std_err = max(corr_14_12_se, propagated_err))
  # Using greater of se of standards or propagated measurement error of stds as error in norm stds.
}

#' Normalize HGIS data
#'
#' @param data A data frame of per-sample data.
#' @param standards A vector of positions to use for standards. Will use target types in `data` if not provided.
#'
#' @return A data frame with normalized per-sample data.
#' @export
#'
norm_hgis <- function(data, standards = NULL) {
  if (!is.null(standards)) {
    data <- data %>% 
      mutate(sample_type = case_when(pos %in% standards ~ "S",
                             sample_type == "S" ~ "U",
                             TRUE ~ sample_type))
  }
  
  stds <- summarize_standards(data)
  
  data %>% 
    mutate(norm_ratio = norm_gas(corr_14_12, stds$corr_14_12_mean),
           sig_norm_ratio = norm_err(corr_14_12, stds$corr_14_12_mean, max_err, stds$norm_std_err) * norm_ratio, 
           d13C = mean(norm_del13c)
          )
}


## Blank correction

#' Blank correct HGIS data
#'
#' @param data A data frame of normalized per-sample data.
#' @param blanks A vector of positions to use for blanks. Will use target types in `data` if not provided.
#' @param fmstd The expected fraction modern of the normalizing standard.
#'
#' @return A data frame with blank corrected fraction modern and error.
#' @export
#'
blank_cor_hgis <- function(data, blanks = NULL, fmstd = 1.0398) {
  if (!is.null(blanks)) {
    data <- data %>% 
      mutate(sample_type = case_when(pos %in% blanks ~ "B",
                             sample_type == "B" ~ "U",
                             TRUE ~ sample_type ))
  }
  
  blanks <- data %>% 
    filter(sample_type == "B")
  
  meanblank <- mean(blanks$norm_ratio)
  
  # get blank error using SNICSer error floor
  blankerr <- amstools::blankErr(blanks$norm_ratio, blanks$sig_norm_ratio) # uses SNICSer error floor method
  # apply blank correction and propagate error
  data %>% 
    mutate(fm_corr = amstools::doLBC(norm_ratio, meanblank, fmstd),
           sig_fm_corr = amstools::doLBCerr(norm_ratio, meanblank, fmstd, sig_norm_ratio, blankerr)
          )
}

# Master data reduction function

#' Reduce HGIS data from an AMS results file
#' 
#'  Steps:
#'  1. Load data
#'  2. Process raw data.
#'  3. Reduce data
#'  4. Normalize
#'  5. Blank Correct
#'
#' @param file A SNICS results file
#' @param date Date sample run if analyzing a specific day
#' @param standards A vector of standard positions.
#' @param blanks A vector of blank positions.
#' @param outliers A dataframe of the position and measurement number of outlier runs.
#' @param remove_outliers Remove outlier runs from analysis if TRUE
#' @param get_consensus Get consensus values for samples if TRUE
#'
#' @return A dataframe of results.
#' @export
#'
reduce_hgis <- function(file, date = NULL, standards = NULL, 
                        blanks = NULL, outliers = NULL,
                        remove_outliers = TRUE, get_consensus = TRUE) {
  df <- process_hgis_results(file, date, standards)
  if (!is.null(outliers)) {
    df <- flag_outliers(df, outliers)
  }
  df_sum <- df %>% 
    sum_hgis_targets(remove_outliers, get_consensus) %>% 
    norm_hgis(standards) %>% 
    blank_cor_hgis(blanks)
  list(raw = df, results = df_sum)
}