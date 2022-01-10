# Functions for HGIS file IO

#' @importFrom magrittr `%>%`
NULL


#' Read SNICS results file
#'
#' @param file A path to a SNICS results file
#'
#' @return A tibble of results for a wheel
#' @export
#'
read_results_file <- function(file) {
  res_cols <- readr::cols(
    runtime = readr::col_datetime(format = "%a %b %d %H:%M:%S %Y"),
    pos = "i",
    meas = "i",
    sample_name = "c",
    sample_type = "c",
    cycles = "d",
    le12C = "d",
    le13C = "d",
    he12C = "d",
    he13C = "d",
    CntTotH = "i",
    CntTotS = "i",
    CntTotGT = "i",
    he13_12 = "d",
    he14_12 = "d"
  )
  
  readr::read_tsv(file, col_names = names(res_cols$cols),
                  col_types = res_cols, skip = 5, comment = "=") %>% 
    dplyr::mutate(corr_14_12 = he14_12/he13_12^2,
	   sig_14_12 = 1/sqrt(CntTotGT),
	   corr_lt = cycles/10)
}


#' Get and process hgis data from results file.
#'
#' Read data file and perform standard munging. 
#' Normalize using mean of samples marked as "S", or a list of standards.
#' This is a rough per-run normalization, not the final Fm of a sample.
#' Use `sum_hgis_targets()` and `norm_hgis()` to produce a per-sample result.
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
    data <- dplyr::filter(data, as.Date(runtime) %in% date)
  }
  
  if (!is.null(standards)) {
    data <- data %>% 
      dplyr::mutate(sample_type = dplyr::case_when(pos %in% standards ~ "S",
                             sample_type == "S" ~ "U",
                             TRUE ~ sample_type))
  }
  
  data <- data %>%
    dplyr::mutate(wheel = stringr::str_extract(file, "USAMS\\d{6}"),
           ok_calc = TRUE)
    
  stds <- summarize_standards(data)
  
  mean13cstd <- data %>%
    dplyr::filter(sample_type == "S") %>% 
    dplyr::pull(he13_12) %>% 
    mean()

  data %>%
    dplyr::mutate(norm_ratio = amsdata::norm_run(corr_14_12, stds$corr_14_12_mean),
           sig_norm_ratio = sig_14_12 * norm_ratio,
           norm_he13_12 = amsdata::norm_run(he13_12, mean13cstd, 1.111618),
           norm_del13c = calc_d13c(norm_he13_12))
}
  

###########
# Older functions below here
###########


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
    mutate(normFm = amsdata::norm_run(cor1412he, meanstd))
}
