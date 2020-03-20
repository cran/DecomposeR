#' @title Visualise the instantaneous frequencies and amplitudes of a
#' decomposition
#'
#' @param pulse a pulse object (created by \code{\link{inst.pulse}} or
#' \code{\link{as.pulse}})
#' @param style whether to plot the distribution of frequency ('d'), the
#' spectral population ('p') or both ('b', is the default)
#' @param breaks,bins,cut parameter for the plots: \code{breaks} is fed to
#' \code{\link{plot_hist}}, \code{bins} is fed to \code{\link{plot_hex}}, and
#' \code{cut} defines the number of color cuts for \code{\link{plot_hex}}.
#' For better control use \code{\link{plot_hist}} and \code{\link{plot_hex}}
#' directly.
#' @param lines the period of lines to be added to the plots for better
#' visualisation
#' @param keep,lose which modes to plot or to not (keep overrides lose)
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
#'                  factor_noise = 10, speak = TRUE)
#'
#' \donttest{
#' plot_emd(dec, dir = tempdir())}
#'
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' ht   <- inst.pulse(dec, plot = FALSE)
#'
#' plot_pulse(ht, lines = c(30, 240))
#'
#' @export

plot_pulse <- function(pulse, style = "b", breaks = 500, bins = 100, cut = 18,
                       lines = NULL, keep = NULL, lose = NULL)
{
  if(!is.pulse(pulse)) stop("The 'pulse' parameter should be an object similar",
                            " to an output of 'inst.pulse' or 'as.pulse'")

  pulse <- mode.out(pulse, lose = lose, keep = keep)

  if(style == "d" | style == "b"){
    opar <- par('mfrow')

    par(mfrow = c(1,1))

    plot_hist(x = 1/pulse$f, breaks = breaks,
              xlog = T, xlab = "Period", main = "Period distribution")

    if(length(unique(as.vector(pulse$mode))) > 1){
      plot_hist(x = 1/pulse$f, breaks = breaks, id = pulse$mode,
                xlog = T, text = T, add = T, line = T, pile = F)
    }

    abline(v = lines, col = "red", lwd = 2, lty = 5)

    par(mfrow = opar)
  }

  if(style == "p" | style == "b"){

    if(length(lines) > 0){
      plot_hex(x = 1/pulse$f, y = pulse$a, bins = bins,
               ybnds = c(0, max(pulse$a)),
               log = "x", trans = log10, inv = function(x) 10^x,
               colorcut = seq(0, 1, length = cut + 1),
               main = "Spectral Population",
               xlab = "Period", ylab = "Amplitude",
               l = list(x = rep(lines, each = 2),
                        y = unit(rep(c(0,1),length(lines)), "npc"),
                        id = rep(seq_len(length(lines)), each = 2),
                        gp = gpar(col = "red", lwd = 2, lty = 5)))
    } else {
      plot_hex(x = 1/pulse$f, y = pulse$a, bins = bins,
               ybnds = c(0, max(pulse$a)),
               log = "x", trans = log10, inv = function(x) 10^x,
               colorcut = seq(0, 1, length = cut + 1),
               main = "Spectral Population",
               xlab = "Period", ylab = "Amplitude")
    }

  }
}
