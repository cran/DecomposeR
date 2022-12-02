#' @title Derive instantaneous frequency
#'
#' @description Calculates instantaneous frequency from an analytic signal.
#'
#' @param asig Analytic signal produced by \code{\link{HilbertTransform}}
#' @param tt Sample times
#' @param method How the instantaneous frequency is calculated. "\code{arctan}"
#' uses the arctangent of the real and imaginary parts of the Hilbert transform,
#' taking the numerical derivative of phase for frequency. "\code{chain}" uses
#' the  analytical derivative of the arctangent function prior to performing the
#' numerical calculation.
#' @param lag Differentiation lag, see the \code{diff} function in the
#' \code{base} package.
#'
#' @return instfreq	Instantaneous frequency in 1/time
#'
#' @note The "\code{arctan}" method was adapted from the \code{hilbertspec}
#' function in the \code{EMD} package.
#'
#' !!IMPORTANT!! The numeric differentiation may be unstable for certain
#' signals. For example, high frequency sinusoids near the Nyquist frequency
#' can give inaccurate results when using the "\code{chain}" method. When in
#' doubt, use the \code{\link{PrecisionTester}} function to check your results!
#'
#' @author Daniel C. Bowman (in the hht package)
#'
#' @seealso \code{\link{PrecisionTester}}
#'
#' @export

InstantaneousFrequency <- function (asig, tt, method = "arctan", lag = 1)
{
  if (!method %in% c("arctan", "chain")) {
    stop(paste("Did not recognize frequency calculation method:",
               method, "Please use either arctan or chain.", sep = " "))
  }
  if (method == "arctan") {
    phase = atan2(Im(asig), Re(asig))
    d = c(0, -diff(phase))
    p = 2 * pi * ((d > pi) - (d < -pi))
    unphase = phase + cumsum(p)
    instfreq = abs(diff(unphase, lag)/diff(tt, lag))
    instfreq = abs(instfreq[-length(instfreq)] + instfreq[-1])/2
    instfreq = c(rep(instfreq[1], 1), instfreq, rep(instfreq[length(instfreq)],
                                                    lag))
  }
  if (method == "chain") {
    dsig = diff(Re(asig), lag)/diff(tt, lag)
    dsig = c(rep(dsig[1], 1), dsig, rep(dsig[length(dsig)],
                                        lag - 1))
    dhsig = diff(Im(asig), lag)/diff(tt, lag)
    dhsig = c(rep(dhsig[1]), dhsig, rep(dhsig[length(dhsig)],
                                        lag - 1))
    instfreq = (Re(asig) * dhsig - Im(asig) * dsig)/((Re(asig)^2) +
                                                       (Im(asig)^2))
  }
  invisible(instfreq/(2 * pi))
}
