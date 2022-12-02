#' @title Test numerically determined instantaneous frequency against exact
#' instantaneous frequency
#'
#' @description This function compares the performance of
#' \code{\link{InstantaneousFrequency}} against signals of known instantaneous
#' frequency. The known signal is of the form
#'
#' \deqn{ x(t) = a\sin(\omega_{1} + \varphi_{1}) + b\sin(\omega_{2} +
#' \varphi_{2}) + c}{a * \sin(omega_1 + phi_1) + b * sin(omega_2 + phi_2) + c}
#'
#' One can create quite complicated signals by choosing the various amplitude,
#' frequency, and phase constants.
#'
#' @param tt Sample times.
#' @param method How the numeric instantaneous frequency is calculated, see
#' \code{\link{InstantaneousFrequency}}
#' @param lag Differentiation lag, see the \code{diff} function in the
#' \code{base} package
#' @param a Amplitude coefficient for the first sinusoid.
#' @param b Amplitude coefficient for the second sinusoid.
#' @param c DC shift
#' @param omega.1 Frequency of the first sinusoid.
#' @param omega.2 Frequency of the second sinusoid.
#' @param phi.1 Phase shift of the first sinusoid.
#' @param phi.2 Phase shift of the second sinusoid.
#' @param plot.signal Whether to show the time series.
#' @param plot.instfreq Whether to show the instantaneous frequencies, comparing
#' the numerical and analytical result.
#' @param plot.error Whether to show the difference between the numerical and
#' analytical result.
#' @param new.device Whether to open each plot as a new plot window (defaults to
#' TRUE). However, Sweave doesn't like \code{dev.new()}. If you want to use
#' \code{PrecisionTester} in Sweave, be sure that new.device = FALSE
#' @param ... Plotting parameters
#'
#' @return
#' \item{instfreq$sig}{The time series}
#' \item{instfreq$analytic}{The exact instantaneous frequency}
#' \item{instfreq$numeric}{The numerically-derived instantaneous frequency from
#' \code{\link{InstantaneousFrequency}}}
#'
#' @author Daniel C. Bowman (in the hht package)
#'
#' @seealso \code{\link{InstantaneousFrequency}}
#'
#' @examples
#' #Simple signal
#'
#' tt <- seq(0, 10, by = 0.01)
#' a <- 1
#' b <- 0
#' c <- 0
#' omega.1 <- 30 * pi
#' omega.2 <- 0
#' phi.1 <- 0
#' phi.2 <- 0
#'
#' PrecisionTester(tt, method = "arctan", lag = 1, a, b, c,
#'                 omega.1, omega.2, phi.1, phi.2, new.device = FALSE)
#'
#' #That was nice - what happens if we use the "chain" method...?
#'
#' PrecisionTester(tt, method = "chain", lag = 1, a, b, c,
#'                 omega.1, omega.2, phi.1, phi.2, new.device = FALSE)
#'
#' #Big problems!  Let's increase the sample rate
#'
#' tt <- seq(0, 10, by = 0.0005)
#' PrecisionTester(tt, method = "chain", lag = 1, a, b, c,
#'                 omega.1, omega.2, phi.1, phi.2, new.device = FALSE)
#'
#' #That's better
#'
#' #Frequency modulations caused by signal that is not symmetric about 0
#'
#' tt <- seq(0, 10, by = 0.01)
#' a <- 1
#' b <- 0
#' c <- 0.25
#' omega.1 <- 2 * pi
#' omega.2 <- 0
#' phi.1 <- 0
#' phi.2 <- 0
#'
#' PrecisionTester(tt, method = "arctan", lag = 1, a, b, c,
#'                 omega.1, omega.2, phi.1, phi.2, new.device = FALSE)
#'
#' #Non-uniform sample rate
#' set.seed(628)
#' tt <- sort(runif(500, 0, 10))
#' a <- 1
#' b <- 0
#' c <- 0
#' omega.1 <- 2 * pi
#' omega.2 <- 0
#' phi.1 <- 0
#' phi.2 <- 0
#'
#' PrecisionTester(tt, method = "arctan", lag = 1, a, b, c,
#'                 omega.1, omega.2, phi.1, phi.2, new.device = FALSE)
#' @importFrom grDevices dev.new
#' @export


PrecisionTester <- function (tt = seq(0, 10, by = 0.01), method = "arctan", lag = 1,
                             a = 1, b = 1, c = 1, omega.1 = 2 * pi, omega.2 = 4 * pi,
                             phi.1 = 0, phi.2 = pi/6, plot.signal = TRUE, plot.instfreq = TRUE,
                             plot.error = TRUE, new.device = TRUE, ...)
{
  A = sin(omega.1 * tt + phi.1)
  B = sin(omega.2 * tt + phi.2)
  C = cos(omega.1 * tt + phi.1)
  D = cos(omega.2 * tt + phi.2)
  sig = a * A + b * B + c
  num = omega.1 * a^2 + omega.2 * b^2 + a * b * (omega.1 +
                                                   omega.2) * (A * B + C * D) + c * (omega.1 * a * A + omega.2 *
                                                                                       b * B)
  denom = a^2 + 2 * a * b * (A * B + C * D) + 2 * c * (a *
                                                         A + b * B) + b^2 + c^2
  analytic.instfreq = num/(denom * 2 * pi)
  asig = HilbertTransform(sig)
  numeric.instfreq = InstantaneousFrequency(asig, tt, method = method,
                                            lag = lag)
  if (plot.signal) {
    if (new.device) {
      dev.new()
    }
    plot(tt, sig, type = "l", xlab = "Time", ylab = "Amplitude",
         main = "Time series", ...)
  }
  if (plot.instfreq) {
    if (new.device) {
      dev.new()
    }
    ylow = min(c(min(analytic.instfreq), min(numeric.instfreq)))
    yhigh = max(c(max(analytic.instfreq), max(numeric.instfreq)))
    plot(tt, analytic.instfreq, type = "l", col = "red",
         ylim = c(ylow, yhigh), xlab = "Time", ylab = "Frequency",
         main = "Analytically and numerically derived values for instantaneous frequency",
         ...)
    lines(tt, numeric.instfreq, lty = 2, col = "black")
    legend("topright", lty = c(1, 2), legend = c("Analytic",
                                                 "Numeric"), col = c("red", "black"))
  }
  if (plot.error) {
    if (new.device) {
      dev.new()
    }
    plot(tt, analytic.instfreq - numeric.instfreq, type = "l",
         xlab = "Time", ylab = "Frequency Error", main = "Numerically derived instantaneous frequency subtracted from analytically derived instantaneous frequency",
         ...)
  }
  instfreq = list(sig = sig, analytic = analytic.instfreq,
                  numeric = numeric.instfreq)
  invisible(instfreq)
}
