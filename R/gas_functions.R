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
#' @param ... Capillary parameters passed to `prestoflow()`
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
#' @param initco2 initial fraction co2
#'
#' @return Fractional concentration of CO2 in vessel (and outflow).
#' @export
#'
concCO2 <- function(t, V = 7000, r = 100, flow = 1, initco2 = 1) {
  stopifnot(flow <= r)
initco2 * exp(-(r/V)*t) * flow
}