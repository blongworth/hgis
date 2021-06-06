# Functions for modelling gas behavior


#' Poiseuille's equation for flow through small tubing
#'
#' Convert delta pressure in Pa to flow in m^3/s
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