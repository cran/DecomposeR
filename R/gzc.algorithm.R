#' @title Calculates instantaneous frequency of simplified IMF using the GZC
#' method
#'
#' @description Calculates instantaneous frequency of simplified IMF using the
#' Generalised Zero-Crossing method from Huang et al., 2009.
#'
#' @param xy a matrix of amplitude
#' @param dt a vector of depth or time values
#'
#' @references Huang, Norden E., Zhaohua Wu, Steven R. Long, Kenneth C. Arnold,
#' Xianyao Chen, and Karin Blank. 2009. ‘On Instantaneous Frequency’. Advances
#' in Adaptive Data Analysis 01 (02): 177–229.
#' https://doi.org/10.1142/S1793536909000096.
#'
#' @details the GZC method is precise to 1/4th of a period, so the results are
#' provided between left and right points, i.e. either an extrema or a
#' zero-crossing.
#'
#' @return a list of $ldt (left position), $rdt (right position), $f
#' (frequency) and $a (amplitude)
#'
#' @examples
#' xyi  <- c(0.5,0,-0.5,0,0.5,0,-0.5,0,0.5,0,-0.5,0,0.5,0,-0.5,0,0.5,0,-0.5,0,
#'          1,1,0,-1,-1,0,1,1,0,-1,-1,0,1,1,0,-1,-1)
#'
#' dti <- 1:length(xyi)
#'
#' d <- simp.emd(m = xyi, dt = dti)
#'
#' xy <- d$xy
#' dt <- d$dt
#'
#' res <- gzc.algorithm(xy, dt)
#'
#' opar <- par('mfrow')
#'
#' par(mfrow = c(2,1))
#'
#' plot(dti, xyi, pch = 19, type = "o", ylab = "xy", xlab = "dt")
#' points(dt, xy, pch = 19, col = "green")
#' points(res$ldt, res$a, pch = 19, col = "red")
#' points(res$rdt, res$a, pch = 19, col = "red")
#'
#' plot(dt, rep(max(res$f, na.rm = TRUE), length(dt)), type = "n",
#'      ylab = "Frequency", xlab = "dt",
#'      ylim = c(0, 2 * max(res$f, na.rm = TRUE)))
#' points(res$ldt, res$f, pch = 19)
#' points(res$rdt, res$f, pch = 19)
#'
#' par(mfrow = opar)
#'
#' @export

gzc.algorithm <- function(xy, dt)
{
  if(!is.simp.emd(xy)) stop("Input is not a simplified IMF")

  test <- apply(dt, MARGIN = 2, function(v) !is.unsorted(v[!is.na(v)]))

  if(!all(test)) stop("The 'dt' columns should be ordered")

  im3 <- abs(lag(dt,3)  - lag(dt,2))
  im2 <- abs(lag(dt,2)  - lag(dt,1))
  im1 <- abs(lag(dt,1)  - dt)
  i0  <- abs(dt         - lead(dt,1))
  ip1 <- abs(lead(dt,1) - lead(dt,2))
  ip2 <- abs(lead(dt,2) - lead(dt,3))
  ip3 <- abs(lead(dt,3) - lead(dt,4))

  freq <- (1/12) * (4/(4*i0) + 2/(2*(im1 + i0)) + 2/(2*(i0  + ip1)) +
                      1/(im3 + im2 + im1 + i0)  + 1/(im2 + im1 + i0  + ip1) +
                      1/(im1 + i0  + ip1 + ip2) + 1/(i0  + ip1 + ip2 + ip3))

  afourth <- (abs(xy) + abs(lead(xy,1)))

  ahalf <- (1/3) * (abs(lag(xy,1))  +     abs(lead(xy,2)) +
                      2*abs(xy)       +   2*abs(lead(xy,1)))

  afull <- (1/10) * (abs(lag(xy,3))   +   abs(lead(xy,4)) +
                       2*abs(lag(xy,2)) + 2*abs(lead(xy,3)) +
                       3*abs(lag(xy,1)) + 3*abs(lead(xy,2)) +
                       4*abs(xy)        + 4*abs(lead(xy,1)))

  ampl <- (1/7) * (4 * afourth + 2 * ahalf + afull)

  ldt <- dt
  rdt <- lead(dt)

  out <- is.na(ampl)

  ldt[out] <- NA
  rdt[out] <- NA

  res <- list(ldt = ldt, rdt = rdt, f= freq, a = ampl)

  return(res)
}
