#' @title The Hilbert transform
#'
#' @description Creates the analytic signal using the Hilbert transform.
#'
#' @param sig Signal to transform.
#'
#' @details Creates the real and imaginary parts of a signal.
#'
#' @return asig Analytic signal
#'
#' @author Daniel C. Bowman (in the hht package)
#'
#' @seealso \code{\link{HilbertEnvelope}}, \code{\link{InstantaneousFrequency}}
#'
#' @examples
#' tt   <- seq(1000) * 0.01
#' sig  <- sin(pi * tt)
#' asig <- HilbertTransform(sig)
#'
#' plot(tt, sig, xlim = c(0, 12))
#'
#' lines(tt, Re(asig), col = "green")
#' lines(tt, Im(asig), col = "red")
#' legend("topright", col = c("black", "green", "red"),
#'        lty = c(NA, 1, 1), pch = c(1, NA, NA),
#'        legend = c("Signal", "Real", "Imaginary"))
#'
#' @importFrom stats fft
#' @export

HilbertTransform <- function (sig)
{
  ndata = length(sig)
  h = rep(0, ndata)
  if (ndata%%2 == 0) {
    h[c(1, ndata/2 + 1)] = 1
    h[2:(ndata/2)] = 2
  }
  else {
    h[1] = 1
    h[2:((ndata + 1)/2)] = 2
  }
  asig = fft(h * fft(sig), inverse = TRUE)/ndata
  invisible(asig)
}

