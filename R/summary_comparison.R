
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
#' @param data A reduced hgis data object
#'
#' @return A summary of replicate sample agreement
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
           fm_diff = fm_corr - fm_consensus,
           sigma = amsdata::sigma(fm_corr, fm_consensus, sig_fm_corr))
}
