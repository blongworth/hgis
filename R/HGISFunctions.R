# Functions for Hybrid GIS calculations


# Pouiselle's; convert delta pressure in Pa to flow in m^3/s
#' Pouiselle's equation for flow through small tubing.
#'
#' @param dp Pressure differential in Pa.
#' @param r Inner radius of capillary in m.
#' @param u Viscosity in kg/m.
#' @param l Length of capillary in m.
#'
#' @return Predicted capillary flow in m^3s-1.
#' @export
#'
#' @examples
prestoflow <- function(dp, r, u, l) {
  # r = radius of capillary
  # l = length of capillary
  # u = viscosity kg/m
  flow <- dp * pi * r^4 / 8 / u / l
  flow
}

#' Calcultate flow through a capillary or tube.
#'
#' @param pres Differential pressure in kPa.
#' @param ... 
#'
#' @return Flow in uL/min.
#' @export
#'
#' @examples
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


# calculate efficiency from ml/m CO2 and uA C12
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

# Normalize each run of a gas target using the mean of another target(s)
norm_gas <- function(sample, standard, stdrat = 1.0398) {
  sample/standard * stdrat
}

# get and process hgis data from local data folder
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
           pos_name = paste(Pos, Sample.Name, sep = " - "))
}
# Summary table of hgis data
sum_hgis <- function(data) {
  data %>% 
    group_by(Pos, Sample.Name, dil_factor) %>%
    summarise(Cur = mean(he12C),
              mean = mean(normFm),
              sd = sd(normFm),
              interr = 1/sqrt(sum(CntTotGT)),
              acqtime = sum(Cycles)/10,
              N_acq = n()) 
}


# boxplot of samples in HGIS test, colored by dilution
plot_hgis <- function(data) {
  
  data %>% 
    mutate(Fm = ifelse(normFm > .15, "modern", "dead"),
           Name = factor(pos_name, levels = unique(pos_name[order(Fm, dil_factor)])),
           Fm = ordered(Fm, levels = c("modern", "dead")),
           Dillution_Ratio = ordered(dil_factor, levels = c(0, 1, 3))) %>%  
    ggplot(aes(Name, normFm, fill = Dillution_Ratio)) +
      geom_boxplot() +
      geom_hline(data = data.frame(yint=1,Fm=ordered("modern", levels = c("modern", "dead"))),
                 aes(yintercept = yint)) + 
      geom_hline(data = data.frame(yint=0,Fm=ordered("dead", levels = c("modern", "dead"))), 
                 aes(yintercept = yint)) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      scale_fill_manual(values = c("0" = "slategray1", "1" = "slategray2", "3" = "slategray3")) +
      labs(x = NULL,
           y = "Fraction Modern") +
      facet_grid(Fm~., scales = "free_y")
}

# plot of normalized ratio vs time
plot_hgis_time <- function(data) {
  ggplot(data, aes(ts, normFm)) +
    geom_point() + 
    facet_wrap(vars(pos_name), scales = "free")
}