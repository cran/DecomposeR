#' @title Plot a decomposition
#'
#' @description General plot for a complete decomposition (that can be
#' summed back to the original signal)
#'
#' @param emd an emd object
#' @param xy the original signal. Is overridden by emd.
#' @param ini an optional vector of length n of the eventual initial Intrinsic
#' Mode Function xy would be a demodulation of, if it is a demodulation.
#' @param dt the depth/time. Is overridden by emd.
#' @param m a matrix with columns of same length that xy, made of the
#' decomposition of the signal. Is overridden by emd.
#' @param repl the replication of decompositions in m. Is overridden by emd.
#' @param size.xy,size.dt the size i inches of each individual plot in pdf
#' @param style whether to not plot the original signal (style = 0), to plot it
#' as the first signal (style = 1), or to plot it before each individual mode
#' (style = 2, is the default)
#' @param xylim,dtlim,inilim the boundaries for the plots (inilim stands for the
#' xy boundaries of the plot of the initial IMF xy is a demodulation of, if
#' applicable)
#' @param vertical whether to have the depth/time [dt] axis vertically
#' (geologist convention) or horizontaly (climatologist convention)
#' @param adapt.axis whether to let the plot adapt the axis to see the
#' variability of the decompositions. The default os to have a comparable x axis
#' for each plots
#' @param adapt.last whether to adapt the last plot as a residue (if TRUE the
#' x axis will be identical to the one of the signal, not centered on 0)
#' @param select the components to plot
#' @param mode which modes/decompositions to plot
#' @param over which modes/decompositions will be cumulated and added to the
#' signal plotted at their left or above them (if style = 2)
#' @param s,o,i,e lists of parameters to feed lines, for the original signal,
#' the cumulated modes/decompositions overlapping it, the
#' modes/decompositions themselves, and the enveloppe of the initial signal
#' used for demodulation if it applies, respectively.
#' @param la,ls,li lists of parameters to provide the abline function (makes
#' personalised lines for you to have a better grasp of the data). la will plot
#' on all panels, ls on the signal ones, and li on the modes ones.
#' @param box whether to draw boxes around the plots
#' @param ax,ay lists of parameters to feed minorAxis, the function making the
#' axes, for the x and y axes
#' @param parg list of parameters to feed par
#' @param title whether to write titles
#' @param t1 the title for the signal
#' @param t2 the title for the modes
#' @param pdf whether to plot as a pdf
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
#'   rnorm(n, sd = 0.5) + 0.01 * t
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' dec <- extricate(xy, dt, nimf = 7,
#'                  repl = 10, comb = 10, factor_noise = 10,
#'                  speak = TRUE)
#'
#' plot_emd(dec, select = c(4,6), pdf = FALSE)
#'\donttest{
#' plot_emd(dec, dir = tempdir())}
#'
#'
#' @importFrom StratigrapheR merge_list minorAxis pdfDisplay multilines
#' @import graphics
#' @export

plot_emd <- function(emd = NULL, xy = NULL, ini = NULL, dt = NULL, m = NULL,
                     mode = NULL, repl = 1,
                     size.xy = 5, size.dt = 25, style = 2,
                     xylim =  NULL, dtlim = NULL, inilim = NULL, vertical = TRUE,
                     adapt.axis = FALSE, adapt.last = TRUE, select = NULL, over = NULL,
                     s = list(type = "o", pch = 19, cex = 0.5),
                     o = list(type = "l", col = "blue", lwd = 2),
                     i = list(type = "o", pch = 19, cex = 0.5),
                     e = list(type = "l", col = "red", lwd = 2),
                     la = list(h = c(), v = c(), col = "red", xpd = FALSE),
                     ls = list(), li = list(col = "grey", lty = 5),
                     box = TRUE, ax = list(), ay = list(), parg = list(),
                     title = TRUE, t1 = "Signal", t2 = "Mode",
                     pdf = TRUE, name = "EMD", ext = ".pdf", dir = tempdir(),
                     track = TRUE, openfile = TRUE)
{

  if(is.null(emd) & is.null(m) & is.null(xy) & is.null(dt)) {

    stop("Missing 'emd', 'xy', 'dt' or 'm' argument")

  } else if(!is.null(emd)) {

    if(!is.emd(emd)) stop("Incorrect 'emd' object")

  } else if (is.null(xy) | is.null(m) | is.null(dt)){

    stop("If 'emd' is NULL, xy, dt & m should be provided")

  } else {

    emd <- as.emd(xy = xy, dt = dt, imf = m, ini = ini, repl = repl)

  }

  xy    <- emd$xy
  ini   <- emd$ini
  dt    <- emd$dt
  m     <- emd$m
  repl  <- unique(emd$repl[1,])
  repln <- length(repl)
  mode  <- unique(emd$mode[1,])

  if(is.null(select)) select <- mode

  if(is.null(over)) over <- mode[-1]

  if(any(integrity(xy = xy, m = m, repl = repl)/sd(xy) > 10^-8) &
     length(over) != 0 & style == 2){
    warning(paste("Integrity of decomposition is not negligible,",
                  "cumulated modes overplotted to the signal may not be",
                  "representative. Set 'style' to 0 or 1 or 'over' to NULL to",
                  "avoid misrepresentation. Type ?additionality for further",
                  "information."))
  }

  if(!isTRUE(vertical) & !isFALSE(vertical)) {
    stop("'vertical' should be TRUE or FALSE")
  }

  if(!isTRUE(adapt.axis) & !isFALSE(adapt.axis)) {
    stop("'adapt.axis' should be TRUE or FALSE")
  }

  if(!isTRUE(title) & !isFALSE(title)) {
    stop("'title' should be TRUE or FALSE")
  }

  if(!isTRUE(box) & !isFALSE(box)) {
    stop("'box' should be TRUE or FALSE")
  }

  if(!isTRUE(pdf) & !isFALSE(pdf)) {
    stop("'pdf' should be TRUE or FALSE")
  }

  nc <- ncol(m)/repln

  id <- rep(seq_len(repln), each = nrow(m))

  opar <- par(unique(c("mfrow", "cex", names(parg))))

  on.exit(par(opar))

  style <- as.integer(style)

  if(style == 2) {
    plots <- 2*length(select)
    if(!is.null(ini)) plots <- plots + 1
  } else if (style == 0) {
    plots <- length(select)
  } else {
    plots <- 1 + length(select)
    if(!is.null(ini)) plots <- plots + 1
  }

  devxy    <- max(xy) - min(xy)
  rangexys <- c(min(xy) - 0.04 * devxy, max(xy) + 0.04 * devxy)

  if(vertical){
    width  <- plots * size.xy
    height <- size.dt
  } else {
    width  <- size.dt
    height <- plots * size.xy
  }

  if(is.null(xylim)){

    if((par("xaxs") == "r" & vertical) | (par("yaxs") == "r" & !vertical)) {
      xylims <- rangexys
    } else {
      xylims <- c(min(xy), max(xy))
    }

  } else {

    xylims <- xylim

  }

  if(!adapt.axis) {
    ranger <- (xylims[2] - xylims[1])/2
    xylimi  <- c(-ranger, ranger)
  }

  if(is.null(dtlim))  dtlim  <- range(dt)
  if(is.null(inilim)) inilim <- c(-max(xy), max(xy))

  si <- merge_list(s, list(type = "o", pch = 19, cex = 0.5))
  ii <- merge_list(i, list(type = "o", pch = 19, cex = 0.5))
  se <- merge_list(e, list(type = "l", col = "red", lwd = 2))

  if(length(over) != 0){
    oi <- merge_list(o, list(type = "l", col = "blue", lwd = 2))
  }

  ax <- merge_list(list(side = 1), ax, list())
  ay <- merge_list(list(side = 2), ay, list())

  vmatrix  <- unlist(m[,seq_mult(ncol(m), ncol(m)/repln)], use.names = F)

  rematrix <- matrix(vmatrix, ncol = ncol(m)/repln)

  m.depth <- rep(dt, repln)

  la <- merge_list(la, list(h = c(), v = c(), col = "red", xpd = F))
  li <- merge_list(li, list(col = "grey", lty = 5))

  graphical_function_plot_emd <- function()
  {

    if(vertical) par(mfrow = c(1, plots)) else par(mfrow = c(plots, 1))

    if(style == 2) {

      if(!is.null(ini)) {

        plot.new()

        if(vertical){


          plot.window(xlim = inilim, ylim = dtlim)

          sinia <- merge_list(list(x = ini, y = dt), si)
          seap  <- merge_list(list(x = xy, y = dt), se)
          seam  <- merge_list(list(x = -xy, y = dt), se)

        } else {

          plot.window(xlim = dtlim, ylim = inilim)

          sinia <- merge_list(list(x = dt, y = ini), si)
          seap  <- merge_list(list(x = dt, y = xy), se)
          seam  <- merge_list(list(x = dt, y = -xy), se)

        }

        if(box) box()

        do.call(graphics::lines, sinia)
        do.call(graphics::lines, seap)
        do.call(graphics::lines, seam)

        do.call(abline, ls)

        do.call(abline, la)

        do.call(minorAxis, ax)
        do.call(minorAxis, ay)

        if(title) title(main = paste("Initial", t1, sep = " "))

      }

      for(imf in select){

        plot.new()

        imff <- which(mode %in% imf)

        if(vertical){

          plot.window(xlim = xylims, ylim = dtlim)

          sj <- merge_list(list(x = xy, y = dt), si)

          if(imf %in% over) {
            os <- rowSums(rematrix[,seq(imff, nc), drop = F])
            oj <- merge_list(list(i = id, x = os, y = m.depth), oi)
          }

        } else {

          plot.window(xlim = dtlim, ylim = xylims)

          sj <- merge_list(list(x = dt, y = xy), si)

          if(imf %in% over) {
            os <- rowSums(rematrix[,seq(imff, nc), drop = F])
            oj <- merge_list(list(i = id, x = m.depth, y = os), oi)
          }

        }

        if(box) box()

        do.call(graphics::lines, sj)

        do.call(abline, ls)

        do.call(abline, la)

        do.call(minorAxis, ax)
        do.call(minorAxis, ay)

        if(imf %in% over) {

          if(title) title(main = paste(t1, "&", t2, imf))

          do.call(multilines, oj)

        } else {
          if(title) title(main = t1)
        }

        imf_main <- paste(t2, imf)

        plot.new()

        if(adapt.axis) {
          xylimi <- range(rematrix[,imff])
        } else if(imf == mode[nc] & isTRUE(adapt.last)){
          xylimi <- xylims
        }

        if(vertical){

          plot.window(xlim = xylimi, ylim = dtlim)

          ij <- merge_list(list(i = id, x = rematrix[,imff], y = m.depth), i, ii)

        } else {

          plot.window(xlim = dtlim, ylim = xylimi)

          ij <- merge_list(list(i = id, x = m.depth, y = rematrix[,imff]), i, ii)

        }

        if(box) box()

        do.call(multilines, ij)

        do.call(abline, li)

        do.call(abline, la)

        do.call(minorAxis, ax)
        do.call(minorAxis, ay)

        if(title) title(main = paste(t2, imf))

      }

    } else {

      if(style != 0) {

        if(!is.null(ini)) {

          plot.new()

          if(vertical){

            plot.window(xlim = inilim, ylim = dtlim)

            sinia <- merge_list(list(x = ini, y = dt), si)
            seap  <- merge_list(list(x = xy, y = dt), se)
            seam  <- merge_list(list(x = -xy, y = dt), se)

          } else {

            plot.window(xlim = dtlim, ylim = inilim)

            sinia <- merge_list(list(x = dt, y = ini), si)
            seap  <- merge_list(list(x = dt, y = xy), se)
            seam  <- merge_list(list(x = dt, y = -xy), se)

          }

          if(box) box()

          do.call(graphics::lines, sinia)
          do.call(graphics::lines, seap)
          do.call(graphics::lines, seam)

          do.call(abline, ls)

          do.call(abline, la)

          do.call(minorAxis, ax)
          do.call(minorAxis, ay)

          if(title) title(main = paste("Initial", t1, sep = " "))

        }

        plot.new()

        if(vertical){

          plot.window(xlim = xylims, ylim = dtlim)

          sj <- merge_list(list(x = xy, y = dt), si)

        } else {

          plot.window(xlim = dtlim, ylim = xylims)

          sj <- merge_list(list(x = dt, y = xy), si)

        }

        if(box) box()

        do.call(graphics::lines, sj)

        do.call(abline, ls)

        do.call(abline, la)

        do.call(minorAxis, ax)
        do.call(minorAxis, ay)

        if(title) title(main = t1)

      }

      for(imf in select){

        imff <- which(mode %in% imf)

        imf_main <- paste(t2, imf)

        plot.new()

        if(adapt.axis) {
          xylimi <- range(rematrix[,imff])
        } else if(imf == mode[nc] & isTRUE(adapt.last)){
          xylimi <- xylims
        }

        if(vertical){

          plot.window(xlim = xylimi, ylim = dtlim)

          ij <- merge_list(list(i = id, x = rematrix[,imff], y = m.depth), i, ii)

        } else {

          plot.window(xlim = dtlim, ylim = xylimi)

          ij <- merge_list(list(i = id, x = m.depth, y = rematrix[,imff]), i, ii)

        }

        if(box) box()

        do.call(multilines, ij)

        do.call(abline, li)

        do.call(abline, la)

        do.call(minorAxis, ax)
        do.call(minorAxis, ay)

        if(title) title(main = paste(t2, imf))

      }

    }

  }

  if(pdf){
    pdfDisplay(graphical_function_plot_emd(), name = name, ext = ext, dir = dir,
               width = width, height = height, parg = parg, track = track,
               openfile = openfile)
  } else {
    par(parg)
    graphical_function_plot_emd()
  }

  par(opar)

}



