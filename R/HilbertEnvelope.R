#' @title Instantaneous amplitude
#'
#' @description Generates the instantaneous amplitude of an analytic signal
#' given by \code{\link{HilbertTransform}}
#'
#' @param asig The analytic signal returned by \code{\link{HilbertTransform}}
#'
#' @return envelope	Instantaneous amplitude
#'
#' @author Daniel C. Bowman (in the hht package)
#'
#' @seealso \code{\link{HilbertTransform}}, \code{\link{InstantaneousFrequency}}
#'
#' @examples
#' tt <- seq(1000) * 0.01
#' sig <- sin(4 * pi * tt) + sin(3.4 * pi * tt)
#' asig <- HilbertTransform(sig)
#' env <- HilbertEnvelope(asig)
#' plot(tt, sig, type = "l")
#' lines(tt, env, col = "red")
#' lines(tt, -env, col = "red")
#'
#' @export

HilbertEnvelope <- function (asig)
{
  envelope = abs(asig)
  invisible(envelope)
}
