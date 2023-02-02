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
#' Version: 1.0.6 (begin of 2023)
#'
#' License: GPL-3
#'
#' @note
#' If you want to use this package for publication or research
#' purposes, please cite Wouters, S., Crucifix, M., Sinnesael, M., Da Silva,
#' A.C., Zeeden, C., Zivanovic, M., Boulvain, F., Devleeschouwer, X., 2022,
#' "A decomposition approach to cyclostratigraphic signal processing".
#' Earth-Science Reviews 225 (103894).
#' <doi:10.1016/j.earscirev.2021.103894>.
#'
#' @author Sebastien Wouters
#'
#' Maintainer: Sebastien Wouters \email{wouterseb@@gmail.com}
#' @importFrom utils globalVariables
#' @importFrom stats approx lm rnorm runif sd spline

utils::globalVariables(c("glo_var_sim", "glo_var_i", "glo_var_abso"))
