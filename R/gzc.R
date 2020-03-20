#' @title Calculates instantaneous frequency using the GZC method
#'
#' @description Calculates instantaneous frequency using the Generalised
#' Zero-Crossing method from Huang et al., 2009. General wrapper for the
#' \code{\link{gzc.algorithm}} function that does all the actual work.
#'
#' @param emd emd-type object
#' @param ini an optional vector of length n of the eventual initial Intrinsic
#' Mode Function xy would be a demodulation of, if it is a demodulation. It will
#' be integrated to the results as mode 1.
#' @param m a matrix of the amplitude values (xy) of the components, each
#' column being a component. Each column should have the same number of non NA
#' values. Vectors, for 1 component, are accepted. Is overridden by emd.
#' @param dt the depth or time value. Is overridden by emd.
#' @param repl the amount of replicates in m. Is overridden by emd.
#' @param mode the mode sequence index to give to each replicated IMFs
#' @param dtout the dt values to sample the frequency and amplitude from if
#' \code{output = 2}.
#' @param output the style of the output, whether 0, 1 or 2. 0 provides the raw
#' output of \code{\link{gzc.algorithm}}, 1 and 2 provides a matrix with $dt
#' (depth/time), $f (frequency) and $a ()amplitude, but with \code{output = 1}
#' the matrix provides the dt only at the extremas and zero-crossings, whereas
#' with \code{output = 2} the dt values
#' are the ones provided with the \code{dtout} parameter. 1 is better for
#' plots, 2 allows easier calculations to be performed downstream.
#' @param warn whether to warn if the sampling interval defined by the
#' \code{dtout} parameter is to small (redirected from
#' \code{StratigrapheR::tie.lim})
#'
#' @return depending on the output parameter:
#'
#' \code{output = 0} provides the raw output of \code{\link{gzc.algorithm}},
#' with $ldt and $rdt (the left and right boundaries of the depth/time
#' intervals), $f (frequency) and $a (amplitude). To that are added $repl (the
#' replicate id) and $mode (the mode id)
#'
#' \code{output = 1} or \code{2} provides a matrix with $dt,
#' $f and $a, but with \code{output = 1} the matrix provides the dt only at
#' the extremas and zero-crossings, whereas with \code{output = 2} the dt values
#' are the ones provided with the \code{out} parameter. \code{1} is better for
#' plots, \code{2} allows easier calculations to be performed downstream.
#'
#'
#' @references Huang, Norden E., Zhaohua Wu, Steven R. Long, Kenneth C. Arnold,
#' Xianyao Chen, and Karin Blank. 2009. "On Instantaneous Frequency". Advances
#' in Adaptive Data Analysis 01 (02): 177–229.
#' https://doi.org/10.1142/S1793536909000096.
#'
#' @examples
#' set.seed(42)
#'
#' n <- 600
#' t <- seq_len(n)
#'
#' p1 <- 30
#' p2 <- 240
#'
#' xy <- (1 + 0.6 * sin(t*2*pi/p2)) * sin(t*2*pi/p1)  + 2 * sin(t*2*pi/p2) +
#'   rnorm(n, sd = 0.5) + t * 0.01
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#' dec <- extricate(xy, dt, nimf = 7, repl = 1, comb = 50,
#'                   factor_noise = 10, speak = TRUE)
#'
#' \donttest{
#' plot_emd(dec, dir = tempdir())}
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' res <- gzc(dec)
#'
#' numb <- 4
#'
#' opar <- par('mfrow')
#'
#' par(mfrow = c(1,2))
#'
#' plot(dec$m[,numb], dec$dt, type = "l",
#'      main = paste("Mode", numb, " + Amplitude"),
#'      xlab = "xy", ylab = "dt", ylim = c(0, 600))
#' lines(res$a[,numb], res$dt[,numb], col = "red", lwd = 2)
#'
#' plot(1/res$f[,numb], res$dt[,numb], ylim = c(0,600),
#'      xlab = "Period", ylab = "dt", log = "x",
#'      type = "l", col = "red", lwd = 2, main = "Period")
#'
#' par(mfrow = opar)
#'
#' @export
#' @importFrom StratigrapheR tie.lim

gzc <- function(emd = NULL, ini = NULL, m = NULL, dt = NULL,
                repl = 1, mode = NULL,
                dtout = NULL, output = 1, warn = TRUE)
{

  if(is.null(emd) & is.null(m) & is.null(dt)) {

    stop("Missing 'emd', 'm' or 'dt' argument")

  } else if(!is.null(emd)) {

    if(!is.emd(emd)) stop("Incorrect 'emd' object")

  } else if (is.null(m) | is.null(dt)){

    stop("If 'emd' is NULL, m & dt should be provided")

  } else {

    fxy <- as.matrix(m)

    if(ncol(fxy) >= 2) fxy <- rowSums(m)

    emd <- as.emd(xy = fxy, dt = dt, ini = ini, imf = m,
                  repl = repl, mode = mode)

  }

  if(isFALSE(last)) emd <- mode.out(emd, lose = max(emd$mode))

  if(!is.null(emd$ini)) emd <- mode.in(emd, emd$ini, mode = 1)

  if(is.null(dtout)) dtout <- emd$dt

  simp  <- simp.emd(emd = emd)

  raw <- gzc.algorithm(simp$xy, simp$dt)

  if(output == 2){

    raw$id <- matrix(rep(seq_len(ncol(raw$f)), each = nrow(raw$f)),
                     ncol = ncol(raw$f))

    freq <- tie.lim(l = raw$ldt, r = raw$rdt, y = raw$f,
                    id = raw$id, xout = dtout, warn = warn)

    ampl <- tie.lim(l = raw$ldt, r = raw$rdt, y = raw$a,
                    id = raw$id, xout = dtout, warn = F)

    nr <- matrix(rep(emd$repl[1,],nrow(freq$y)), nrow = nrow(freq$y), byrow = T)
    nm <- matrix(rep(emd$mode[1,],nrow(freq$y)), nrow = nrow(freq$y), byrow = T)

    res <- list(dt = dtout, f= freq$y, a = ampl$y, repl = nr, mode = nm)

  } else if(output == 0){

    nr <- matrix(rep(emd$repl[1,],nrow(raw$f)), nrow = nrow(raw$f), byrow = T)
    nm <- matrix(rep(emd$mode[1,],nrow(raw$f)), nrow = nrow(raw$f), byrow = T)

    res <- list(ldt = raw$ldt, rdt = raw$rdt, f = raw$f, a = raw$a,
                repl = nr, mode = nm)

  } else {

    ndt <- c(raw$ldt,raw$rdt)
    ndt <- ndt[seq_mult(length(ndt), length(ndt)/2)]
    ndt <- matrix(ndt, ncol = ncol(raw$ldt))

    nf <- c(raw$f,raw$f)
    nf <- nf[seq_mult(length(nf), length(nf)/2)]
    nf <- matrix(nf, ncol = ncol(raw$f))

    na <- c(raw$a,raw$a)
    na <- na[seq_mult(length(na), length(na)/2)]
    na <- matrix(na, ncol = ncol(raw$a))

    nr <- matrix(rep(emd$repl[1,],nrow(na)), nrow = nrow(na), byrow = T)
    nm <- matrix(rep(emd$mode[1,],nrow(na)), nrow = nrow(na), byrow = T)

    res <- list(dt = ndt, f = nf, a = na, repl = nr, mode = nm)

  }

  return(res)

}
