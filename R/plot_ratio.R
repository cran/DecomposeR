#' @title Visualise the instantaneous frequencies ratios of a
#' decomposition
#'
#' @param ratio a ratio object (created by \code{\link{inst.ratio}}
#' @param sqrt.rpwr whether to use the squared ratio power (i.e. the squared
#' multiplication of the instantaneous amplitudes of the modes two by two)
#' rather than the ratio power itself.
#' @param style whether to plot a single plot in the graphics device ('s'), the
#' to plot an ensemble of all the ratios combinations in a pdf ('e'), or both
#' ('b', is the default)
#' @param select the groups of ratios combinations to plot in the single plot
#' (in the "1/2" form)
#' @param bins,cut parameter for the plots: \code{bins} is fed to
#' \code{\link{plot_hex}}, and \code{cut} defines the number of color cuts for
#' \code{\link{plot_hex}}. For better control use \code{\link{plot_hex}}
#' directly.
#' @param lines the ratio of lines to be added to the plots for better
#' visualisation
#' @param plot whether to plot. Otherwise output a grob of the single plot.
#' @param width,height the width  and height in inches of each separate plot
#' in the ensemble of all the ratios combinations
#' @param name,ext,dir,track,openfile parameters for the pdfDisplay function,
#' namely the name of the pdf file, its extension (if you want to make a .svg
#' file you can), the directory of the file, whether to track the changes
#' (if you use sumatrapdf as a default pdf reader you can set it to F and it
#' will avoid creating too many pdf files), and whether to directly open the
#' file
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
#' dec <- extricate(xy, dt, nimf = 7, sifting = 10,
#'                  repl = 10, comb = 10,
#'                  factor_noise = 10, speak = TRUE)
#' \dontrun{
#' plot_emd(dec, dir = tempdir())}
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' ht    <- inst.pulse(dec, plot = FALSE)
#' ratio <- inst.ratio(ht, plot = FALSE)
#'
#' plot_ratio(ratio, lines = c(8), style = "s")
#' plot_ratio(ratio, lines = c(8), style = "s", select = c("4/6"))
#' \dontrun{
#' plot_ratio(ratio, lines = c(8), style = "e", dir = tempdir())}
#'
#' @export


plot_ratio <- function(ratio, sqrt.rpwr = TRUE, style = "b", select = NA,
                       bins = 100, cut = 18, lines = NULL, plot = TRUE,
                       width = 10, height = 10, name = "Ratio", ext = ".pdf",
                       dir = tempdir(), track = TRUE, openfile = TRUE)
{

  if(!is.ratio(ratio)) {
    stop("The 'ratio' parameter should be a ratio object like the ones",
         " made by the inst.ratio() function, see ?inst.ratio and ?is.ratio",
         " for further information")
  }

  x <- ratio$ratio

  if(isTRUE(sqrt.rpwr)) {
    y <- sqrt(ratio$rpwr)
    ylab <- "Squared Ratio Power"
  } else {
    y <- ratio$rpwr
    ylab <- "Ratio Power"
  }

  if(!is.null(lines)){

    ll <- length(lines)

    l <- list(x = rep(lines, each = 2),
              y = unit(rep(c(0,1), ll), "npc"),
              id = rep(seq_len(ll), each = 2),
              gp =  gpar(col = "red", lwd = 2, lty = 5))

  } else {
    l <- list()
  }

  if(style == "b" | style == "s" | isFALSE(plot)){

    if(isTRUE(plot)){
      plot_hex(x = x, y = y, bins = bins, colorcut = seq(0, 1, length = cut + 1),
               log = "x", trans = log10, inv = function(x) 10^x,
               id = ratio$lr, select = select,
               main = "Ratio Population", xlab = "Ratio", ylab = ylab, l = l)
    } else {
      out <- plot_hex(x = x, y = y, bins = bins,
                      colorcut = seq(0, 1, length = cut + 1),
                      log = "x", trans = log10, inv = function(x) 10^x,
                      id = ratio$lr, select = select,
                      main = "Ratio Population", xlab = "Ratio", ylab = ylab,
                      l = l, plot = F)
      return(out)
    }

  }

  if((style == "b" | style == "e") & isTRUE(plot)){

    all <- unique(ratio$lr[1,])

    basis <- length(unique(ratio$l[1,]))

    matx <- matrix(rep(seq_len(basis), basis), nrow = basis, byrow = T)
    maty <- matrix(rep(seq_len(basis), basis), nrow = basis)

    keep <- as.vector(t(upper.tri(matx, diag = T)))

    px <- as.vector(t(matx))[keep]
    py <- as.vector(t(maty))[keep]

    ref <- seq_len(length(all))

    g <- function(){

      grid.newpage()

      lay <- grid.layout(ncol = basis, nrow = basis)

      vp1 <- viewport(layout = lay, name = "vp1")

      pushViewport(vp1)

      for(i in ref)
      {

        pushViewport(viewport(layout.pos.row = py[i],
                              layout.pos.col = px[i]))

        j <- plot_hex(x = x, y = y, bins = bins,
                      colorcut = seq(0, 1, length = cut + 1), id = ratio$lr,
                      log = "x", trans = log10, inv = function(x) 10^x,
                      main = paste("Ratio Population", all[i]),
                      select = all[i], xlab = "Ratio", ylab = ylab,
                      l = l, plot = F)

        grid.rect()

        grid.draw(j)

        seekViewport("vp1")

      }

    }

    pdfDisplay(g(), name = name, ext = ext,
               width = width * basis, height = height * basis,
               track = track, openfile = openfile, dir = dir)

  }

}




