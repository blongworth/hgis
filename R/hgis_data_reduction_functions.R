# Functions for HGIS data reduction

# Add 13C and deadtime correction
doCor1412 <- function(he1412, ltcorr, he1312) {
	he1412 / ltcorr / he1312 ^ 2
}

# Calc internal error for a measurement
calcSig1412 <- function(CntTotH, CntTotS, CntTotGT, cor1412) {
RelErrSq <- (CntTotH - CntTotS) * CntTotH ^ 2 / CntTotS ^ 4 + 
             CntTotH ^ 2 / CntTotGT / CntTotS ^ 2
cor1412 * RelErrSq ^ 0.5
}

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

## Normalize
# Find standards


# Find mean of stds
normStds <- function(cor1412std, defstd) {
  # both vectors, same length or 1 element vector if all standards are same
  mean(cor1412std/defstd) 
}

# Normalize to mean of standards
norm1412 <- function(cor1412, meanstd) {
  cor1412/meanstd
}

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




## Blank correction

# Apply large blank
doLBC <- function(fmmeas, fmblank, fmstd) {
	fmmeas - fmblank * (fmstd - fmmeas) / fmstd
}

# Propagate large blank error
doLBCerror <- function(fmmeas, fmblank, fmstd, fmmeaserr, fmblankerr) {
	sqrt(fmmeaserr ^ 2 * (1 + fmblank / fmstd) ^ 2 +
	fmblankerr ^ 2 * ((fmmeas - fmstd) / fmstd) ^ 2)
}

# 
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
    group_by(pos, sample_name, sample_type) %>% 
    summarize(ext_err = amstools::se(corr_14_12),
              int_err = 1/sqrt(sum(CntTotGT)),
              max_err = pmax(ext_err, int_err),
              across(c(le12C, le13C, he12C, he13C,
                       he13_12, he14_12, corr_14_12, 
                       norm_ratio, norm_del13c), 
                     list(mean = mean, sd = sd),
                     .names = "{ifelse(.fn == 'mean', '', 'sig_')}{.col}"),
              runtime = min(runtime),
              counts = sum(CntTotGT),
              n_runs = n()
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
norm_hgis <- function(data, standards) {
  
  if (!missing(standards)) {
    data <- data %>% 
      mutate(Num = case_when(pos %in% standards ~ "S",
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


#' Blank correct HGIS data
#'
#' @param data A data frame of normalized per-sample data.
#' @param blanks A vector of positions to use for blanks. Will use target types in `data` if not provided.
#' @param fmstd The expected fraction modern of the normalizing standard.
#'
#' @return A data frame with blank corrected fraction modern and error.
#' @export
#'
blank_cor_hgis <- function(data, blanks, fmstd = 1.0398) {
  
  if (!missing(blanks)) {
    data <- data %>% 
      mutate(sample_type = case_when(pos %in% blanks ~ "B",
                             sample_type == "B" ~ "U",
                             TRUE ~ sample_type ))
  }
  
  meanblank <- data %>% 
    filter(sample_type == "B") %>% 
    pull(norm_ratio) %>% 
    mean()
  
  blankerr <- data %>% 
    filter(sample_type == "B") %>% 
    pull(norm_ratio) %>% 
    amstools::blankErr() # uses SNICSer error floor method
  
  data %>% 
    mutate(fm_corr = amstools::doLBC(norm_ratio, meanblank, fmstd),
           sig_fm_corr = amstools::doLBCerr(norm_ratio, meanblank, fmstd, sig_norm_ratio, blankerr)
          )
}

# Produce normalized data for wheel per target


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
reduce_hgis <- function(file, date, standards, blanks) {
  #Load data
  
  
  
}
