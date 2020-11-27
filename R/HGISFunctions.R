# Functions for Hybrid GIS calculations


# Pouiselle's; convert delta pressure in Pa to flow in m^3/s
prestoflow <- function(dp, r, u, l) {
  # r = radius of capillary
  # l = length of capillary
  # u = viscosity kg/m
  flow <- dp * pi * r^4 / 8 / u / l
  flow
}

# convert gauge kPa to uL/min
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
    mungeResfile() %>% 
filter(as.Date(ts) %in% date)
  
  meanstd <- mean(data$cor1412he[data$Num == "S"])
  
  data <- data %>% 
    mutate(normFm = norm_gas(cor1412he, meanstd),
           dil_factor = case_when(str_ends(Sample.Name, "1") ~ 1,
                                 str_starts(Sample.Name, "Dil") ~ 3,
                                 TRUE ~ 0))
}
# Summary table of hgis data
sum_hgis <- function(data) {
  data %>% 
    group_by(Pos, Sample.Name, dil_factor) %>%
    summarise(Cur = mean(he12C),
              mean = mean(normFm),
              sd = sd(normFm),
              interr = 1/sqrt(sum(CntTotGT)),
              N = n()) 
}

# boxplot of samples in HGIS test, colored by dilution
# TODO: put in lines for expected value of modern gas standard
plot_hgis <- function(data) {
  
  data <- data %>% 
    mutate(Fm = ifelse(normFm > .15, "modern", "dead"),
           Fm = ordered(Fm, levels = c("modern", "dead")),
           Dillution_Ratio = as.factor(dil_factor)) 
  
  ord <- data %>% 
    select(Sample.Name, Fm, dil_factor) %>% 
    arrange(desc(Fm), dil_factor) %>% 
    pull(Sample.Name) %>% 
    unique()
  
  ggplot(data, aes(factor(Sample.Name, level = ord), normFm, fill = Dillution_Ratio)) +
    geom_boxplot() +
    geom_hline(data = data.frame(yint=1,Fm=ordered("modern", levels = c("modern", "dead"))), aes(yintercept = yint)) + 
    geom_hline(data = data.frame(yint=0,Fm=ordered("dead", levels = c("modern", "dead"))), aes(yintercept = yint)) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
    labs(x = NULL,
         y = "Fraction Modern") +
    facet_grid(Fm~., scales = "free_y")
}

