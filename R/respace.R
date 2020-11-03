#' @title Interpolate with even spacing
#'
#' @description Interpolate with even spacing. Can determine on its own the
#' most conservative sampling interval (using the Greatest Common Rational
#' Divisor)
#'
#' @param dt depth/time (same length than length/rows of xy)
#' @param xy signal (vector or matrix)
#' @param delta the new sampling interval. If NULL, uses the Greatest Common
#' Rational Divisor
#' @param tolerance,relative parameters for the \code{divisor} function
#' (\code{StratigrapheR} package), to compute the Greatest Common
#' Rational Divisor
#' @param n.warn the amount of interpolated points in between the largest
#' interval above which a warning is provided. This warning can be useful to
#' avoid needlessly long outputs, which might make any subsequent computation
#' take too much time.
#'
#' @return a list of interpolated xy and dt values ($xy and $dt), plus a vector
#' of logicals indicating whether each point was part of the initial input or
#' was added by interpolation
#'
#' @examples
#' set.seed(42)
#'
#' n <- 50
#' t <- seq_len(n)
#'
#' xy <- (1 + 0.6 * sin(t*0.025)) * sin(t*0.2)  + 2 * sin(t*0.025) +
#'         rnorm(n, sd = 0.5)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5), 1)
#'
#' dt <- cumsum(inter_dt)
#'
#' res <- respace(xy = xy, dt = dt)
#'
#' opar <- par("mfrow")
#'
#' par(mfrow = c(1,1))
#'
#' plot(res$xy, res$dt, type = "l")
#' points(res$xy[res$initial], res$dt[res$initial], pch = 19, col = "green")
#' points(res$xy[!res$initial], res$dt[!res$initial],
#'        pch = 19, col = "red", cex = 0.5)
#'
#' par(mfrow = opar)
#'
#' @importFrom StratigrapheR divisor is.divisor
#' @export

respace <- function(dt, xy = NULL, delta = NULL,
                    tolerance = 8, relative = TRUE, n.warn = 100)
{

  dec.dt <- min(dt)

  dt <- dt - dec.dt

  if(!is.null(xy)) xy <- as.matrix(xy)

  if(any(is.na(xy))) stop("'xy' should not have any NA values")
  if(any(is.na(dt))) stop("'dt' should not have any NA values")

  if(!is.null(delta)){

    if(!is.divisor(dt, delta, tolerance = tolerance, relative = relative)) {
      stop(paste("'delta' (", delta, ") is not a divisor of the Greatest ",
                 "Common Rational Divisor of the 'dt' values (GCRD = ",
                 divisor(dt, tolerance = tolerance,
                         relative = relative, speak = F), ")",
                 sep = ""))
    }

    interval <- delta

  } else {

    interval <- divisor(dt, tolerance = tolerance,
                        relative = relative, speak = F)

  }

  if(is.unsorted(dt)){
    ord <- order(dt)
    dt  <- dt[ord]
    if(!is.null(xy)) xy  <- xy[ord,,drop = F]
  }

  ran <- range(dt)

  ran.int <- ran/interval

  s.int.prime <- as.integer(ran.int[2] - ran.int[1] + 1.1)

  error <- (max(dt) - as.integer(ran.int[2] + 0.1) * interval)/s.int.prime

  if(!is.null(xy)) nc <- ncol(xy) else nc <- 1

  s <- (seq_len(s.int.prime * nc) * interval) +
    (seq_len(s.int.prime) * error)

  nap <- s - s[1] + ran[1]

  mdt <- matrix(rep(dt, nc), ncol = nc)

  increase <- cumsum(c(0, rep(max(dt) - min(dt) + interval -
                                (error * (s.int.prime - 1)),
                              nc - 1)))
  add      <- matrix(rep(increase, each = length(dt)), ncol = nc)

  ndt <- as.vector(mdt + add)
  if(!is.null(xy)) nxy <- as.vector(xy)

  in.int <- seq_len(s.int.prime) - 1 + as.integer(dt[1]/interval + 0.1)

  initial <- in.int %in% as.integer(0.1 + dt/interval)

  approx_dt <- (in.int * interval) + (in.int * error)

  approx_dt[initial] <- dt

  if(!is.null(xy)){

    approx_xy <- approx(x = ndt, y = nxy, xout = nap, rule = 2)$y

    minitial            <- rep(initial, nc)
    approx_xy[minitial] <- nxy

    if(nc != 1) approx_xy <- matrix(approx_xy, ncol = nc)

    res <- list(dt = approx_dt, xy = approx_xy, initial = initial)

  } else {

    res <- list(dt = approx_dt, initial = initial)

  }

  res$dt <- res$dt + dec.dt

  # Warning if the amount of points in between the largest interval is too
  # important

  if(is.null(delta) & !is.null(n.warn)){

    if(!(inherits(n.warn,"numeric") | inherits(n.warn, "integer"))){
      stop("The 'n.warn' parameter should be a numeric or an integer")
    }

    in.max <- as.integer((max(abs(lag(dt) - dt), na.rm = T)/interval) - 1 + 0.1)

    if(in.max > n.warn) {

      dl <- length(res$dt)

      warning("There are ", in.max,
              " interpolated points between the largest interval.",
              "\nThis exceeds the warning threshold of ", n.warn ,
              ".\nThis brings the amount of points to a total of ", dl,
              ".\nThis might make any subsequent computation",
              " go more slowly than needed.",
              "\nTo solve this problem you can:",
              "\n - Round the dt values at a reasonable scale.",
              "\nTo avoid this warning you can:",
              "\n - Set the 'n.warn' parameter at a higher value.",
              "\n - Set the 'n.warn' parameter as NULL.")

    }

  }

  return(res)
}



