#' @title Calculates instantaneous frequency of freqeuncy carriers using the DQ
#' method
#'
#' @description Calculates instantaneous frequency of frequency carriers using
#' the direct quadrature method from Huang et al., 2009.
#'
#' @param fc a matrix of amplitude between -1 and 1, making up the frequency
#' carrier
#' @param dt a vector of depth or time values
#'
#' @references Huang, Norden E., Zhaohua Wu, Steven R. Long, Kenneth C. Arnold,
#' Xianyao Chen, and Karin Blank. 2009. "On Instantaneous Frequency". Advances
#' in Adaptive Data Analysis 01 (02): 177–229.
#' https://doi.org/10.1142/S1793536909000096.
#'
#' @return a list of the depth/time (dt), frequency (f), and identity tuning
#' (idt), i.e. depths adapted to transform the frequency carrier into a cosine
#' of period 1.
#'
#' @examples
#' n <- 600
#'
#' t <- seq_len(n)
#'
#' p1 <- 30
#'
#' xy <- sin(t*2*pi/p1 + 50)
#'
#' int <- c(rep(1, 99 + 100), seq(1,3,2/100), seq(3,1,-2/100), rep(1,100 + 99))
#'
#' dt <- cumsum(int)
#'
#' cond <- dt < 75
#'
#' xy <- xy[!cond]
#' dt <- dt[!cond]/1.2 - 62.5
#'
#' res <- dq.algorithm(xy, dt)
#'
#' opar <- par("mfrow")
#'
#' par(mfrow = c(3,1))
#'
#' plot(dt, xy, type = "o", pch = 19, main = "Frequency carrier")
#'
#' plot(dt, 1/res$f, pch = 19, type = "l", log = "y", lwd = 2, ylim = c(25,80),
#'      main = "Period (Direct Quadrature method)", ylab = "Period")
#'
#' plot(res$idt[,1], xy, type = "o", pch = 19,
#'      main = "Identity tuning", axes = FALSE, ylab = "xy", xlab = "dt")
#'
#' ap <- approx(x = dt, y = res$idt[,1], xout = seq(0,600, by = 20))
#'
#' axis(1, at = ap$y, labels = ap$x)
#' axis(2)
#' box()
#'
#' par(mfrow = opar)
#'
#' @export
#' @importFrom dplyr lead lag

dq.algorithm <- function(fc, dt)
{
  if(any(fc > 1 | fc < - 1)) stop(paste("fc values should be between -1 and 1;",
                                        "fc needs to be a frequency carrier"))

  fc <- as.matrix(fc)

  nc <- ncol(fc)

  li <- nrow(fc)

  if(li != length(dt)) stop("dt should be as long as fc (in length or rows)")

  dtna <- rep(c(dt, NA), nc)

  fcna <- as.vector(rbind(fc, matrix(rep(NA, nc), ncol = nc)))

  extr <- extremist(fcna, local = F, zc = F)

  maxe <- extr$maxindex
  mine <- extr$minindex

  l <- length(fcna)

  ext <- seq_len(l) %in% c(maxe$l, maxe$r,mine$l, mine$r)

  fcna[ext] <- sign(fcna[ext])

  ext[is.na(fcna)] <- NA

  extc <- ext & lead(!ext) & lag(!ext)

  extc[1] <- F
  extc[l] <- F

  df <- data.frame(ext = ext, extc = extc, dt = dtna, fc = fcna)

  df$phase <- atan(df$fc/sqrt(1-df$fc^2))

  df$phaselag  <- lag(df$phase)
  df$phaselead <- lead(df$phase)

  df$dtlag  <- lag(df$dt)
  df$dtlead <- lead(df$dt)

  df$IF <- (abs(df$phase - df$phaselag) + abs(df$phase - df$phaselead))/
    abs(2*pi*(df$dtlag - df$dtlead))

  df$correct <- df$extc * abs(df$phase - df$phaselag)/(df$IF * 2 * pi)

  df$ndt <- df$dt

  df$ndt[df$extc] <- df$dtlag[df$extc] + df$correct[df$extc]

  df$ndtlag  <- lag(df$ndt)
  df$ndtlead <- lead(df$ndt)

  df$nIF <- (abs(df$phase - df$phaselag) + abs(df$phase - df$phaselead))/
    abs(2*pi*(df$ndtlag - df$ndtlead))

  first <- seq(1, nc * (li+1), li + 1)
  last  <- seq(li, nc * (li+1), li + 1)


  df$nIF[first] <- abs(df$phase[first] - df$phaselead[first])/
    abs(2*pi*(df$ndt[first] - df$ndtlead[first]))

  df$nIF[last] <- abs(df$phase[last] - df$phaselag[last])/
    abs(2*pi*(df$ndt[last] - df$ndtlag[last]))

  df$dp <- abs(df$phase - df$phaselag)

  # Extract direct quadrature IF ----

  dq <- matrix(df$nIF, ncol = nc)

  dq <- dq[-(li+1),,drop = F]

  # Extract identity tuning----

  dp <- matrix(df$dp, ncol = nc)

  dp <- dp[-(li+1),,drop = F]

  exm <- matrix(extc, ncol = nc)

  dp[1,] <- rep(0, nc)

  idt <- apply(dp, 2, cumsum)/(2*pi)

  fex <- apply(exm, 2, function(x) which(x)[1])

  sidt <- idt[fex, seq_len(nc)]
  sfc  <- sign(fc[fex, seq_len(nc)])

  add <- rep(0, nc)
  add[sfc == -1] <- 0.5

  dephase <- matrix(rep(sidt + add, each = li), ncol = nc)

  res <- list(dt = dt, f = dq, idt = idt - dephase + 1)

  return(res)
}



