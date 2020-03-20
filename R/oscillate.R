#' @title Modify a signal using a Van der Pol oscillator
#'
#' @param xy initial signal (vector or matrix)
#' @param dt depth/time (same length than length/rows of xy)
#' @param period the period of the oscillator (length 1 or n)
#' @param delta the sampling interval for iteration (length 1 or n)
#' @param damp damping parameter
#' @param normalise whether to recenter the output signal on the initial signal
#' @param limit whether to warn when parameters are irrealistic (subjective)
#' @param f.noise a factor of the amount of noise (length 1 or n)
#' @param f.signal a factor of the amount of signal (length 1 or n)
#' @param dx,dy the differentials used in the oscillator. They should be
#' provided as functions needing x, y, beta (2*pi/period) and damp (damping)
#' parameters
#' @param xi the initial x value
#' @param yi the initial y value
#'
#' @examples
#' set.seed(42)
#'
#' n <- 800
#'
#' dt <- seq(0,n, 1)
#'
#' p1 <- 100
#' p2 <- 40
#'
#' xy <- (1 + 0.6 * sin(dt*2*pi/p1)) * sin(dt*2*pi/p2)  + 2 * sin(dt*2*pi/p1) + 1
#'
#' xyout <- oscillate(xy, dt, period = 30)
#'
#' opar <- par("mfrow")
#'
#' par(mfrow = c(1,1))
#'
#' plot(xy, dt, type = "l",
#'      main = "Initial signal (bold) & oscillated signal (dashed)",
#'      lwd = 2, xlim = c(-4, 6))
#' lines(xyout, dt, type = "l", col = "grey50", lwd = 2, lty = 5)
#'
#' par(mfrow = opar)
#'
#' @importFrom StratigrapheR homogenise
#' @export

oscillate <- function(xy, dt, period, delta = 0.05, damp = 0.00005,
                      f.noise = 5, f.signal = 0.95,
                      dx = function(x, y, beta, damp) beta * y - x * (x^2 + y^2 - 1) * damp,
                      dy = function(x, y, beta, damp) -beta * x  - y * (x^2 + y^2 - 1) * damp,
                      xi = if(length(xy) != 0) xy[1] else 0.5,
                      yi = if(length(xy) != 0) xy[1] else 0.5,
                      normalise = TRUE, limit = TRUE)
{

  if(length(xy) != length(dt)) stop("'xy' and 'dt' should be of same length")

  if(isTRUE(limit)){

    if(damp > 0.05) warning("'damp' should be less than 0.05")

  }

  para <- homogenise(i = xy, cycle = F,
                     l = list(period = period, damp = damp,
                              f.noise = f.noise, f.signal = f.signal))

  period   <- para$period
  damp     <- para$damp
  f.noise  <- para$f.noise
  f.signal <- para$f.signal

  beta <- 2*pi/period

  # Respace

  mp <- as.matrix(data.frame(xy = xy, beta = beta, damp = damp,
                             f.noise = f.noise, f.signal = f.signal))

  ns <- respace(dt, xy = mp, delta)

  nxy <- ns$xy[,1]
  ndt <- ns$dt

  nbeta     <- ns$xy[,2]
  ndamp     <- ns$xy[,3]
  nf.noise  <- ns$xy[,4]
  nf.signal <- ns$xy[,5]

  # Noise ----

  nl <- length(nxy)

  noise1 <- rnorm(nl, sd = 1) * nf.noise
  noise2 <- rnorm(nl, sd = 1) * nf.noise

  # Iterate ----

  # Initial conditions ----

  xa <- xi
  ya <- yi

  # Loop ----

  speaky <- round(seq(0, length(ndt) - 1, length.out = 21))

  for(i in seq_len(length(ndt) - 1))
  {
    if(i %in% speaky) print(paste(round(100*i/length(ndt)), "%"))

    xn <- xi + dx(x = xi, y = yi, beta = nbeta[i], damp = ndamp[i]) * delta  +
      noise1[i] * sqrt(delta) +
      nf.signal[i] * nxy[i] * delta

    yn <- yi + dy(x = xi, y = yi, beta = nbeta[i], damp = ndamp[i]) * delta +
      noise2[i] * sqrt(delta)

    xa <- c(xa, xn)
    ya <- c(ya, yn)

    xi <- xn
    yi <- yn

  }

  xyout <- xa[ns$initial]

  # Normalise ----

  if(normalise) {
    xyout1 <- (xyout - mean(xyout))/(max(xyout) - min(xyout))
    xyout2 <- (max(xy) - min(xy)) * xyout1
    xyout <- mean(xy) + xyout2
  }

  return(xyout)

}





