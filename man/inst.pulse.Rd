% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inst.pulse.R
\name{inst.pulse}
\alias{inst.pulse}
\title{Computes instantaneous frequency using the Hilbert transform}
\usage{
inst.pulse(
  emd = NULL,
  imf = NULL,
  m = NULL,
  dt = NULL,
  ini = NULL,
  repl = 1,
  mode = NULL,
  last = FALSE,
  plot = TRUE,
  method = "HT",
  delta = NULL,
  tolerance = 8,
  relative = TRUE,
  breaks = 500,
  bins = 100,
  cut = 18,
  lines = NULL
)
}
\arguments{
\item{emd}{an emd object}

\item{imf}{a matrix of same frequency modes to calculate the frequency from.
Is overridden by emd. This allows to calculate and visualise the results
for single IMFs more clearly than in a population plot.}

\item{m}{a matrix of the modes to calculate the frequency from. Is overridden
by emd and imf.}

\item{dt}{the depth or time. Is overridden by emd.}

\item{ini}{an optional vector of length n of the eventual initial Intrinsic
Mode Function xy would be a demodulation of, if it is a demodulation. It will
be integrated to the results as mode 1.}

\item{repl}{the amount of replicates in m. Is overridden by emd.}

\item{mode}{the mode sequence index to give to each replicated IMFs.
Is overridden by emd.}

\item{last}{whether to use the last mode (trend/residue).}

\item{plot}{whether to have a plot summary of the output.}

\item{method}{the IF calculation method: "HT" for Hilbert transform
(default), "NHT" for normalised Hilbert transform, and "DQ" for direct
quadrature. The two last require normalisation, which can sometimes fail.}

\item{delta, tolerance, relative}{parameters to feed to \code{\link{respace}}
for interpolation}

\item{breaks, bins, cut}{parameter for the plots: \code{breaks} is fed to
\code{\link{plot_hist}}, \code{bins} is fed to \code{\link{plot_hex}}, and
cut defines the number of color cuts for \code{\link{plot_hex}}. For better
control use \code{\link{plot_hist}} and \code{\link{plot_hex}} directly.}

\item{lines}{the period of lines to be added to the plots for better
visualisation}
}
\value{
a list made of $dt (depth/time), $f (instantaneous frequency), $a
(instantaneous amplitude),$repl (the replicate id of each point) and
$mode (the mode id of each point)
}
\description{
Calculates instantaneous frequency using the Hilbert transform
(HT), normalised Hilbert transform (NHT) or the direct quadrature (DQ)
methods. Normalisation is done for NHT and DQ using Huang et al., 2009
algorithm, but the empirical normalisation scheme can fail due to overshoot
or undershoot of the spline. Additional research is necessary for that last
feature.
}
\examples{
set.seed(42)

n <- 600
t <- seq_len(n)

p1 <- 30
p2 <- 240

xy <- (1 + 0.6 * sin(t*2*pi/p2)) * sin(t*2*pi/p1)  + 2 * sin(t*2*pi/p2) +
  rnorm(n, sd = 0.5) + t * 0.01

inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)

dt <- cumsum(inter_dt)
dec <- extricate(xy, dt, nimf = 7, repl = 10, comb = 10,
                  factor_noise = 10, sifting = 10, speak = FALSE)
\dontrun{
plot_emd(dec, dir = tempdir())}

integrity(xy, dec)
parsimony(dec)

ht   <- inst.pulse(dec, lines = c(30, 240))
gzcr <- gzc(dec)

imf <- dec$m[,4]

inst.pulse(imf = imf, dt = dt, method = "DQ")

}
\references{
Huang, Norden E., Zhaohua Wu, Steven R. Long, Kenneth C. Arnold,
Xianyao Chen, and Karin Blank. 2009. "On Instantaneous Frequency". Advances
in Adaptive Data Analysis 01 (02): 177–229.
https://doi.org/10.1142/S1793536909000096.
}
