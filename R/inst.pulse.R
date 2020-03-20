#' @title Computes instantaneous frequency using the Hilbert transform
#'
#' @description Calculates instantaneous frequency using the Hilbert transform
#' (HT), normalised Hilbert transform (NHT) or the direct quadrature (DQ)
#' methods. Normalisation is done for NHT and DQ using Huang et al., 2009
#' algorithm, but the empirical normalisation scheme can fail due to overshoot
#' or undershoot of the spline. Additional research is necessary for that last
#' feature.
#'
#' @param emd an emd object
#' @param imf a matrix of same frequency modes to calculate the frequency from.
#' Is overridden by emd. This allows to calculate and visualise the results
#' for single IMFs more clearly than in a population plot.
#' @param m a matrix of the modes to calculate the frequency from. Is overridden
#' by emd and imf.
#' @param dt the depth or time. Is overridden by emd.
#' @param ini an optional vector of length n of the eventual initial Intrinsic
#' Mode Function xy would be a demodulation of, if it is a demodulation. It will
#' be integrated to the results as mode 1.
#' @param repl the amount of replicates in m. Is overridden by emd.
#' @param mode the mode sequence index to give to each replicated IMFs.
#' Is overridden by emd.
#' @param last whether to use the last mode (trend/residue).
#' @param plot whether to have a plot summary of the output.
#' @param method the IF calculation method: "HT" for Hilbert transform
#' (default), "NHT" for normalised Hilbert transform, and "DQ" for direct
#' quadrature. The two last require normalisation, which can sometimes fail.
#' @param delta,tolerance,relative parameters to feed to \code{\link{respace}}
#' for interpolation
#' @param breaks,bins,cut parameter for the plots: \code{breaks} is fed to
#' \code{\link{plot_hist}}, \code{bins} is fed to \code{\link{plot_hex}}, and
#' cut defines the number of color cuts for \code{\link{plot_hex}}. For better
#' control use \code{\link{plot_hist}} and \code{\link{plot_hex}} directly.
#' @param lines the period of lines to be added to the plots for better
#' visualisation
#'
#' @return a list made of $dt (depth/time), $f (instantaneous frequency), $a
#'(instantaneous amplitude),$repl (the replicate id of each point) and
#' $mode (the mode id of each point)
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
#' dec <- extricate(xy, dt, nimf = 7, repl = 10, comb = 10,
#'                   factor_noise = 10, speak = TRUE)
#'
#' \donttest{
#' plot_emd(dec, dir = tempdir())}
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' ht   <- inst.pulse(dec, lines = c(30, 240))
#' gzcr <- gzc(dec)
#'
#' imf <- dec$m[,4]
#'
#' inst.pulse(imf = imf, dt = dt, method = "DQ")
#'
#' @importFrom hht HilbertTransform InstantaneousFrequency HilbertEnvelope
#' @export

inst.pulse <- function(emd = NULL, imf = NULL, m = NULL, dt = NULL,
                       ini = NULL, repl = 1, mode = NULL,
                       last = FALSE, plot = TRUE, method = "HT",
                       delta = NULL, tolerance = 8, relative = TRUE,
                       breaks = 500, bins = 100, cut = 18, lines = NULL)
{


  if(is.null(emd) & is.null(imf) & is.null(m) & is.null(dt)) {

    stop("Missing 'emd', 'imf', m' or 'dt' argument")

  } else if(!is.null(emd)) {

    if(!is.emd(emd)) stop("Incorrect 'emd' object")

  } else if (is.null(dt) | (is.null(imf) & is.null(m))){

    stop("If 'emd' is NULL, 'imf' or 'm' should be provided, along with 'dt'")

  } else {

    if(!is.null(imf)){

      fxy <- as.matrix(imf)

      emd <- as.emd(xy = rep(0,nrow(fxy)), dt = dt, ini = ini, imf = imf,
                    repl = seq_len(ncol(fxy)), mode = mode)

      last <- T

    } else {

      fxy <- as.matrix(m)

      if(ncol(fxy) >= 2) fxy <- rowSums(m)

      emd <- as.emd(xy = fxy, dt = dt, ini = ini, imf = m,
                    repl = repl, mode = mode)

    }

  }

  if(isFALSE(last)) emd <- mode.out(emd, lose = max(emd$mode))

  if(!is.null(emd$ini)) emd <- mode.in(emd, emd$ini, mode = 1)

  m   <- emd$m
  dt  <- emd$dt

  if(ncol(m) == 0) {
    stop("No signal is treated: set parameter 'last' as TRUE or check 'xy'")
  }

  if(method == "DQ"){

    norm <- normalise(dt = dt, m = m, last = T, speak = F)

    dq   <- dq.algorithm(norm$fc, dt)

    freq <- dq$f
    idt  <- dq$idt
    amp  <- norm$a

    res <- as.pulse(dt = dt, m = m, f = freq, a = amp, idt = idt,
                    repl = emd$repl, mode = unique(as.vector(emd$mode)))

  } else {

    spaced <- respace(dt = dt, xy = m, delta = delta, tolerance = tolerance,
                      relative = relative)

    sm  <- as.matrix(spaced$xy)
    sdt <- spaced$dt

    HilbertFre <- function(ifm){
      ht  <- HilbertTransform(ifm)
      out <- InstantaneousFrequency(ht, sdt, method = "arctan", lag = 1)
    }

    if(method == "NHT") {

      norm <- normalise(dt = sdt, m = sm, last = T, speak = F)
      freq <- as.matrix(apply(norm$fc, 2, HilbertFre))
      amp  <- norm$a

    } else {

      HilbertAmp <- function(ifm){
        ht  <- HilbertTransform(ifm)
        out <- HilbertEnvelope(ht)
      }

      freq <- as.matrix(apply(sm, 2, HilbertFre))
      amp  <- as.matrix(apply(sm, 2, HilbertAmp))

    }

    freq <- freq[spaced$initial,,drop = F]
    amp  <- amp[spaced$initial,,drop = F]

    res <- as.pulse(dt = dt, m = m, f = freq, a = amp,
                    repl = emd$repl)

  }

  if(isTRUE(plot)){

    if(!is.null(imf)){

      plot_imf(res)

    } else {

      plot_pulse(res, breaks = breaks, bins = bins, cut = cut, lines = lines)

    }

  }

  return(invisible(res))

}







