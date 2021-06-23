# Functions for HGIS data reduction


#' Calculate d13C
#'
#' @param he1312 
#'
#' @return
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
              int_err = 1/sqrt(counts) * mean(corr_14_12),
              max_err = pmax(ext_err, int_err),
              across(c(le12C, le13C, he12C, he13C,
                       he13_12, he14_12, corr_14_12), 
                     list(mean = mean, sd = sd),
                     .names = "{ifelse(.fn == 'mean', '', 'sig_')}{.col}"),
              ) %>% 
    mutate(rec_num = case_when(str_starts(sample_name, "LiveGas") ~ 101730,
                               str_starts(sample_name, "C-?1") ~ 83028,
                               str_starts(sample_name, "TIRI-F") ~ 2138,
                               str_starts(sample_name, "TIRI-I") ~ 17185,
                               str_starts(sample_name, "C-?2") ~ 1082,
                               str_starts(sample_name, "NOSAMS-?2") ~ 38809,
                               str_starts(sample_name, "DeadGas") ~ 72446))
  if (get_consensus == TRUE) {
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
  meanstd <- data %>% 
    filter(sample_type == "S") %>% 
    pull(corr_14_12) %>% 
    mean()
  # Using se of standards as standard error. Should probably compare to per-standard error
  stderr <- data %>% 
    filter(sample_type == "S") %>% 
    pull(corr_14_12) %>% 
    se()
  data %>% 
    mutate(norm_ratio = norm_gas(corr_14_12, meanstd),
           sig_norm_ratio = norm_err(corr_14_12, meanstd, max_err, stderr) * norm_ratio 
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
#'  3. Insert raw data into DB
#'  4. Reduce data
#'  5. Normalize
#'  6. Blank Correct
#'  7. Insert results into DB
#'
#' @param file A SNICS results file
#' @param date Date sample run if analysing a specific day
#' @param standards A vector of standard positions.
#' @param blanks A vector of blank positions.
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