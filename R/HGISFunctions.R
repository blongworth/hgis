# Functions for Hybrid GIS calculations


# Poiseuille's; convert delta pressure in Pa to flow in m^3/s
#' Poiseuille's equation for flow through small tubing.
#'
#' @param dp Pressure differential in Pa.
#' @param r Inner radius of capillary in m.
#' @param u Viscosity in kg/m.
#' @param l Length of capillary in m.
#'
#' @return Predicted capillary flow in m^3s-1.
#' @export
#'
prestoflow <- function(dp, r, u, l) {
  # r = radius of capillary
  # l = length of capillary
  # u = viscosity kg/m
  flow <- dp * pi * r^4 / 8 / u / l
  flow
}

#' Calculate flow through a capillary or tube.
#'
#' @param pres Differential pressure in kPa.
#' @param ... 
#'
#' @return Flow in uL/min.
#' @export
#'
flowcalc <- function(pres, ...) {
  dppa <- pres * 1E3 #convert kPa to Pa
  flowcms <- prestoflow(dppa, ...)
  flowuls <- flowcms * 1E6 * 1E3
  flowulm <- flowuls * 60
  flowulm
}

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

#' Get and process hgis data from results file.
#'
#' Read data file and perform standard munging with amstools::mungeResfile(). 
#' Add dilution factor. Normalize using mean of samples marked as "S".
#' 
#' @param file Path to a NOSAMS format AMS results file.
#' @param date A vector of dates to subset from file. All data used if missing.
#'
#' @return
#' @export
#'
get_hgis_data <- function(file, date) {
  data <- readResfile(file) %>% 
    mungeResfile()
  
  if (!missing(date)) {
    data <- filter(data, as.Date(ts) %in% date)
  }
  
  meanstd <- mean(data$cor1412he[data$Num == "S"])
  
  data <- data %>% 
    mutate(normFm = norm_gas(cor1412he, meanstd),
           dil_factor = case_when(str_ends(Sample.Name, "1") ~ 1,
                                 str_starts(Sample.Name, "Dil") ~ 3,
                                 TRUE ~ 0),
           pos_name = reorder(paste(Pos, Sample.Name, sep = " - "), as.numeric(Pos))) %>% 
    group_by(pos_name) %>% 
    mutate(cum_acqtime = cumsum(Cycles) / 10) %>% 
    ungroup()
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
    group_by(Pos, Sample.Name, dil_factor) %>%
    summarise(Cur = mean(he12C),
              Cur.sd = sd(he12C),
              mean = mean(normFm),
              sd = sd(normFm),
              exterr = normRunExtErr(normFm),
              interr = 1/sqrt(sum(CntTotGT)),
              acqtime = sum(Cycles)/10,
              N_acq = n()) 
}


#' Boxplot of samples in HGIS test, colored by dilution
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

#' Plot of normalized ratio vs time
#'
#' @param data A tibble in format of output from get_hgis_data().
#'
#' @return A ggplot object with norm ratio and current vs time,
#' faceted for each pos_name.
#' 
#' @export
#'
plot_hgis_time <- function(data) {
  ggplot(data, aes(cum_acqtime, normFm)) +
    geom_point() + 
    facet_wrap(vars(pos_name), scales = "free")
}


#' Calculate fraction of CO2 in a container vs. time when displaced by another gas
#' 
#' CO2 is displaced from a vial by introducing helium through one side 
#' of a double needle and allowing the mixture to flow out via the 
#' other side of the needle. Pressure is maintained at 1ATM as the 
#' outlet is open to air.
#'
#' TODO: Allow specification of a starting mixture of HE and CO2.
#' 
#' @param t time in same units as r
#' @param V volume of vessel in same units as r
#' @param r rate of displacement gas flow
#' @param flow rate of gas flow in capillary or outflow. Default of 1 returns fraction of CO2. Should be less than r.
#'
#' @return Fractional concentration of CO2 in vessel (and outflow).
#' @export
#'
#' @examples
concCO2 <- function(t, V = 7000, r = 100, flow = 1, initco2 = 1) {
  stopifnot(flow <= r)
  initco2 * exp(-(r/V)*t) * flow
}