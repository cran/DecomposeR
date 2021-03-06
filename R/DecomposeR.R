#' @name DecomposeR
#' @title DecomposeR: Empirical Mode Decomposition for Cyclostratigraphy
#' @description This package provides tools to apply Ensemble Empirical Mode
#' Decomposition (EEMD) for cyclostratigraphy purposes. It proposes a new
#' algorithm, that performs EEMD in seconds, a linear interpolation algorithm
#' using the greatest rational common divisor of depth or time, different
#' algorithms to compute instantaneous amplitude, frequency and ratios of
#' frequencies, and functions to verify and visualise the outputs.
#'
#' @details Package: DecomposeR
#'
#' Type: R package
#'
#' Version: 1.0.4 (Fall 2020)
#'
#' License: GPL-3
#'
#' @note
#' If you want to use this package for publication or research
#' purposes, please cite Wouters, S., Da Silva, A.C. Crucifix,
#' M., Sinnesael, M., Zivanovic, M., Boulvain, F., Devleeschouwer, X., 2019,
#' Litholog generation with the StratigrapheR package and signal decomposition
#' for cyclostratigraphic purposes. Geophysical Research Abstracts Vol. 21,
#' EGU2019-5520, 2019, EGU General Assembly 2019.
#' <http://hdl.handle.net/2268/234402>
#'
#' @author Sebastien Wouters
#'
#' Maintainer: Sebastien Wouters \email{sebastien.wouters@@doct.uliege.be}
#' @importFrom utils globalVariables
#' @importFrom stats approx lm rnorm runif sd spline

utils::globalVariables(c("glo_var_sim", "glo_var_i", "glo_var_abso"))
