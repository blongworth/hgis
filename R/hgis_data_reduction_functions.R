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

# Calculate d13C
calcd13c <- function(he1312) {
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


# Propagate normalization error

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
norm_err <- function(sample, standard, sample_err, standard_err) {
  sqrt(sample_err^2/sample^2 + standard_err^2/standard^2)
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
#' @param data A data frame of run data, processed by `get_hgis_data`.
#' @param remove_outliers If `TRUE`, do not include samples where `outlier == TRUE`.
#' @param get_consensus If `TRUE`, get and add consensus values from NOSAMS DB.
#'
#' @return A data frame of per-sample data.
#' @export
#'
sum_hgis_targets <- function(data, remove_outliers = TRUE, get_consensus = TRUE) {
  
  if (remove_outliers == TRUE) {
    data <- filter(data, !outlier)
  }
  
  

data_sum <- data %>% 
    group_by(Pos, Sample.Name, Num) %>% 
    summarize(ext_err = amstools::se(cor1412he) / mean(cor1412he),
              int_err = 1/sqrt(sum(CntTotGT)),
              max_err = pmax(ext_err, int_err),
              across(c(le12C, le13C, he12C, he13C, X13.12he, X14.12he, cor1412he, normFm), 
                     list(mean = mean, sd = sd),
                     .names = "{ifelse(.fn == 'mean', '', 'sig_')}{.col}"),
              ts = min(ts),
              counts = sum(CntTotGT),
              n_runs = n()
              ) %>% 
    mutate(rec_num = case_when(str_starts(Sample.Name, "LiveGas") ~ 101730,
                               str_starts(Sample.Name, "C-?1") ~ 83028,
                               str_starts(Sample.Name, "TIRI-F") ~ 2138,
                               str_starts(Sample.Name, "TIRI-I") ~ 17185,
                               str_starts(Sample.Name, "C-?2") ~ 1082,
                               str_starts(Sample.Name, "NOSAMS-?2") ~ 38809,
                               str_starts(Sample.Name, "DeadGas") ~ 72446))
    
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
      mutate(Num = case_when(Pos %in% standards ~ "S",
                             Num == "S" ~ "U",
                             TRUE ~ Num))
  }
  
  meanstd <- data %>% 
    filter(Num == "S") %>% 
    pull(cor1412he) %>% 
    mean()
  
  # Using sd of standards as standard error. Should probably compare to per-standard error
  stderr <- data %>% 
    filter(Num == "S") %>% 
    pull(cor1412he) %>% 
    sd()
  
  data %>% 
    mutate(norm_ratio = norm_gas(cor1412he, meanstd),
           sig_norm_ratio = norm_err(cor1412he, meanstd, max_err, stderr) * norm_ratio # Replace with proper error propagation
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
      mutate(Num = case_when(Pos %in% blanks ~ "B",
                             Num == "B" ~ "U",
                             TRUE ~ Num))
  }
  
  meanblank <- data %>% 
    filter(Num == "B") %>% 
    pull(norm_ratio) %>% 
    mean()
  
  blankerr <- data %>% 
    filter(Num == "B") %>% 
    pull(norm_ratio) %>% 
    amstools::blankErr() # uses SNICSer error floor method
  
  data %>% 
    mutate(fm_corr = amstools::doLBC(norm_ratio, meanblank, fmstd),
           sig_fm_corr = amstools::doLBCerr(norm_ratio, meanblank, fmstd, sig_norm_ratio, blankerr)
          )
}

# Produce normalized data for wheel per target


#' 
#'
#' @param data 
#' @param standards 
#' @param blanks 
#'
#' @return
#' @export
#'
reduce_hgis <- function(data, standards, blanks) {
  #Load data
  
  
}
