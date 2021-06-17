# Functions for Hybrid GIS calculations


## Efficiency

# convert L CO2 to atoms C
lCO2toatC <- function (vol) {
  molC <- vol * 1/22.4 # 22.4 L of gas per mole at stp
  atC <- molC * 6.022E23
  atC
}

# convert current in A to atoms C/sec
CcurtoatC <- function(cur, cs){
  at12C <- cur * 6.25E18 # electrons per coulomb
  at12C <- at12C / cs # divide by charge state
  at12C
}


#' Calculate efficiency from ml/m CO2 and uA C12
#'
#' @param ulm Flow in ulm-1.
#' @param uA Post-accelerator 12C ion current in uA.
#' @param cs Charge state to convert particle-uA to 12C current.
#'
#' @return Efficiency as carbon atoms in vs carbon atoms measured.
#' @export
#'
hgis_eff <- function(ulm, uA, cs = 3) {
  
  # if flow is zero or less, efficiency doesn't make sense
  ulm <- ifelse(ulm <= 0, NA, ulm)
  
  vol <- ulm * 1E-6 # conv to L
  cur <- uA * 1E-6 # conv to A
  atCin <- lCO2toatC(vol) / 60 # gives atoms/s
  atC12out <- CcurtoatC(cur, cs)
  atCout <- atC12out / .99 # 99% of C is 12C, we need total C
  atC12out / atCin
}


#' Read SNICS results file
#'
#' @param file A path to a SNICS results file
#'
#' @return A tibble of results for a wheel
#' @export
#'
read_results_file <- function(file) {
  res_cols <- cols(
    runtime = col_datetime(format = "%a %b %d %H:%M:%S %Y"),
    pos = col_integer(),
    meas = col_integer(),
    sample_name = col_character(),
    sample_type = col_character(),
    cycles = col_double(),
    le12C = col_double(),
    le13C = col_double(),
    he12C = col_double(),
    he13C = col_double(),
    CntTotH = col_integer(),
    CntTotS = col_integer(),
    CntTotGT = col_integer(),
    he13_12 = col_double(),
    he14_12 = col_double()
  )
  
  read_tsv(file, col_names = names(res_cols$cols), col_types = res_cols, skip = 5, comment = "=") %>%
    mutate(corr_14_12 = he14_12/he13_12^2,
	   sig_14_12 = 1/sqrt(CntTotGT),
	   corr_lt = cycles/10)
}


#' Get and process hgis data from results file.
#'
#' Read data file and perform standard munging. 
#' Normalize using mean of samples marked as "S", or a list of standards.
#' 
#' @param file Path to a NOSAMS format AMS results file.
#' @param date A vector of dates to subset from file. All data used if missing.
#' @param standards A vector of the wheel positions containing standards
#'
#' @return
#' @export
#'
process_hgis_results <- function(file, date = NULL, standards = NULL) {
  data <- read_results_file(file)

  if (!is.null(date)) {
    data <- filter(data, as.Date(runtime) %in% date)
  }
  
  if (!is.null(standards)) {
    data <- data %>% 
      mutate(sample_type = case_when(pos %in% standards ~ "S",
                             sample_type == "S" ~ "U",
                             TRUE ~ sample_type))
  }
  
  data <- data %>%
    mutate(pos_name = reorder(paste(pos, sample_name, sep = " - "), pos),
           wheel = str_extract(file, "USAMS\\d{6}"),
           ok_calc = TRUE)
    
  meanstd <- data %>% 
    filter(sample_type == "S") %>% 
    pull(corr_14_12) %>% 
    mean()
  
  mean13cstd <- data %>%
    filter(sample_type == "S") %>% 
    pull(he13_12) %>% 
    mean()

  data %>%
    mutate(norm_ratio = norm_gas(corr_14_12, meanstd),
           sig_norm_ratio = sig_14_12 * norm_ratio,
           norm_del13c = calc_d13c(he13_12))
}
  
	   
#' Get and process hgis data from results file.
#'
#' Read data file and perform standard munging with amstools::mungeResfile(). 
#' Add dilution factor. Normalize using mean of samples marked as "S", or a list of standards.
#' 
#' @param file Path to a NOSAMS format AMS results file.
#' @param date A vector of dates to subset from file. All data used if missing.
#' @param standards A vector of the wheel postions containing standards
#'
#' @return
#' @export
#'
get_hgis_data <- function(file, date = NULL, standards = NULL) {
  data <- readResfile(file) %>% 
    mungeResfile() 
  
  if (!is.null(date)) {
    data <- filter(data, as.Date(ts) %in% date)
  }
  
  if (!is.null(standards)) {
    data <- data %>% 
      mutate(Num = case_when(Pos %in% standards ~ "S",
                             Num == "S" ~ "U",
                             TRUE ~ Num))
  }
  
  data <- data %>%
    mutate(dil_factor = case_when(str_ends(Sample.Name, "1") ~ 1,
                                 str_starts(Sample.Name, "Dil") ~ 3,
                                 TRUE ~ 0),
           pos_name = reorder(paste(Pos, Sample.Name, sep = " - "), as.numeric(Pos)),
           wheel = str_extract(file, "USAMS\\d{6}")) %>% 
    group_by(pos_name) %>% 
    mutate(cum_acqtime = cumsum(Cycles) / 10,
           outlier = ifelse(is.na(removeOutliers(cor1412he)), TRUE, FALSE)) %>% 
    ungroup()
    
  meanstd <- data %>% 
    filter(Num == "S" & !outlier) %>% 
    pull(cor1412he) %>% 
    mean()
  
  data <- data %>%
    mutate(normFm = norm_gas(cor1412he, meanstd))
}


#' Summary table of hgis data
#'
#' @param data A tibble in format of output from get_hgis_data().
#'
#' @return A tibble with summary data.
#' @export
#'
sum_hgis <- function(data) {
  data %>% 
    group_by(wheel, Pos, Sample.Name, dil_factor) %>%
    summarise(Cur = mean(he12C),
              Cur.sd = sd(he12C),
              mean = mean(normFm),
              sd = sd(normFm),
              exterr = se(normFm),
              interr = mean * (1/sqrt(sum(CntTotGT))), # multiplied internal error by fm, not sure if correct.
              merr = pmax(exterr, interr),
              acqtime = sum(Cycles)/10,
              N_acq = n()) 
}


#' Compare agreement of replicate samples
#' 
#' Will use normalized data and error if blank corrected data unavailable.
#'
#' @param data 
#'
#' @return
#' @export
#'
compare_replicates <- function(data) {
  
  # Use normalized data if BC data is not available
  if ("fm_corr" %in% names(data)) {
    cols <- c("he12C", "fm_corr", "sig_fm_corr")
  } else {
    cols <- c("he12C", "norm_ratio", "sig_norm_ratio")
    warning("fm_corr unavailable, using norm_ratio for comparison")
  }
  
  data %>% 
    mutate(he12C = he12C * 1E6) %>% 
    group_by(rec_num) %>% 
    filter(n() > 1) %>% 
    summarize(Name = str_remove(sample_name[1], "_.*$"),
              across(cols,
                     list(mean = mean, sd = sd)),
              N = n()) 
}

#' Compare data to consensus Fm
#'
#' @param data A results dataframe with fm_consensus, fm_corr, and sig_fm_corr
#'
#' @return A data frame showing agreement with consensus for all samples.
#' @export
#'
compare_consensus <- function(data) {
  data %>% 
    filter(!is.na(fm_consensus)) %>% 
    mutate(Name = case_when(rec_num == 101730 ~ "LiveGas",
                                rec_num == 83028 ~ "C-1",
                                rec_num == 2138 ~ "TIRI-F",
                                rec_num == 17185 ~ "TIRI-I",
                                rec_num == 1082 ~ "C-2",
                                rec_num == 38809 ~ "NOSAMS2",
                                rec_num == 72446 ~ "DeadGas"),
           Fm_diff = fm_corr - fm_consensus,
           sigma = amstools::sigma(fm_corr, fm_consensus, sig_fm_corr))
}