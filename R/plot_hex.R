#' @title Group and/or log-scale hexagonal binning
#'
#' @description Group and/or log-scale hexagonal binning. Provides a legend
#' indicating the count representations. USES THE GRID GRAPHICAL SYSTEM, BASE
#' GRAPHICS NOT SUPPORTED. To add lines, polygons or text, use the l, g and t
#' arguments.
#'
#' @param x,y vectors giving the coordinates of the bivariate data points to be
#' binned.
#' @param id a vector of ids for each x value, to separate different groups of
#' data
#' @param select the groups of ids to plot
#' @param uniform whether to keep the creaks defined by the entire matrixes when
#' selecting only a part of it
#' @param bins the number of bins partitioning the range of xbnds.
#' @param xbnds,ybnds horizontal and vertical limits of the binning region in
#' x or y units respectively; must be numeric vector of length 2.
#' @param xlim,ylim the limits of the plot
#' @param log a character string which contains "x" if the x axis is to be
#' logarithmic, "y" if the y axis is to be logarithmic and "xy" or "yx" if both
#' axes are to be logarithmic.
#' @param shape the theoretical shape = yheight/xwidth of the plotting. This
#' adapts the form of the hexagons accordingly.
#' @param mincnt,maxcnt fraction of cell area for the lowest and largest count,
#' respectively
#' @param colorcut vector of values covering [0, 1] that determine hexagon color
#' class boundaries and hexagon legend size boundaries. Alternatively, an
#' integer (<= maxcnt) specifying the number of equispaced colorcut values in
#' [0,1].
#' @param colramp function accepting an integer n as an argument and returning n
#' colors.
#' @param trans a transformation function for the counts such as
#' \code{\link{log10}}
#' @param inv	the inverse transformation function (if
#' \code{trans = \link{log10}}, \code{inv} should for instance be
#' \code{function(x) 10^x}.
#' @param border the color of the border of the hexagons. By default it will be
#' the color of the filling
#' @param lwd the width of the border of the hexagons.
#' @param cex the magnification of text.
#' @param main main title.
#' @param xlab,ylab x and y axis labels respectively.
#' @param xaxis,yaxis whether to plot the x and y axes respectively.
#' @param xaxs,yaxs The style of axis interval calculation to be used for the
#' axes. By default the style "r" (regular) first extends the data range by 4
#' percent at each end and then finds an axis with pretty labels that fits
#' within the extended range. Style "i" (internal) just finds an axis with
#' pretty labels that fits within the original data range.
#' @param box whether to plot a box.
#' @param mar a numerical vector of the form c(bottom, left, top, right) which
#' gives the room the give to the margins in Normalised Parent Coordinates
#' (see \code{grid} package for more information)
#' @param legend whether to plot the legend.
#' @param leg_sep the distance between hexagons and text f the legend in
#' Normalised Parent Coordinates left on the right margin
#' @param xpd_hex factor to expand the legend hexagons
#' @param xpd_leg factor to expand the height of the legend
#' @param l a list of arguments to feed to \code{grid::grid.polyline}
#' ATTENTION the grid package has to be loaded
#' @param g a list of arguments to feed to \code{grid::grid.polygon}
#' ATTENTION the grid package has to be loaded
#' @param t a list of arguments to feed to \code{grid::grid.text}
#' ATTENTION the grid package has to be loaded
#' @param plot whether to plot. If FALSE, returns a grob.
#'
#' @examples
#' library(grid) # To use the gpar function
#'
#' set.seed(42)
#'
#' n <- 600
#' t <- seq_len(n)
#'
#' p1 <- 30
#' p2 <- 240
#'
#' xy <- (1 + 0.6 * sin(t*2*pi/p2)) * sin(t*2*pi/p1)  + 2 * sin(t*2*pi/p2) +
#'         rnorm(n, sd = 0.5)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' dec <- extricate(xy, dt, nimf = 7, sifting = 10,
#'                 repl = 10, comb = 10, factor_noise = 10,
#'                  speak = FALSE)
#'
#' \dontrun{
#' plot_emd(dec, dir = tempdir())}
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' ht  <- inst.pulse(dec, plot = FALSE)
#'
#' plot_hex(x = 1/ht$f, y = ht$a, bins = 100, ybnds = c(0,2),
#'          log = "x", trans = log10, inv = function(x) 10^x,
#'          main = "Spectral Population", xlab = "Period", ylab = "Amplitude")
#'
#' plot_hex(x = 1/ht$f, y = ht$a, bins = 100, ybnds = c(0,2),
#'          log = "x", trans = log10, inv = function(x) 10^x,
#'          main = "Spectral Population", xlab = "Period", ylab = "Amplitude",
#'          id = ht$mode, select = c(4,6,7),
#'          l = list(x = c(30, 30, 240, 240), y = unit(c(0,1,0,1), "npc"),
#'                  id = c(1,1,2,2), gp = gpar(col = c("red", "blue"), lwd = 2)),
#'         g = list(x = c(18, 50, 50, 18, 18, 50, 50, 18),
#'                  y = c(0, 0, 1.9, 1.9, 2.05, 2.05, 1.95, 1.95),
#'                  id = c(1,1,1,1,2,2,2,2),
#'                  gp = gpar(col = c("red", NA), fill = c(NA, "white"), lwd = 2)),
#'         t = list(label = "Mode 4", x = 30, y = 2, gp = gpar(col = "red")))
#'
#' @export
#' @import grid
#' @import hexbin
#' @importFrom StratigrapheR seq_log
#' @importFrom colorRamps matlab.like

plot_hex <- function(x, y, id = NA, select = NA, uniform = TRUE, bins = 60,
                     xbnds = range(x, na.rm = TRUE), ybnds = range(y, na.rm = TRUE),
                     xlim = xbnds, ylim = ybnds, log = "", shape = 1,
                     mincnt = 1, maxcnt = NA,
                     colorcut = seq(0,1, length = 17),
                     colramp = function(n) matlab.like(length(colorcut)-1),
                     trans = NULL, inv = NULL,
                     border = NULL, lwd = 0.1, cex = 1,
                     main = "", xlab = "x", ylab = "y",
                     xaxis = TRUE, yaxis = TRUE, xaxs = "r", yaxs = "r", box = TRUE,
                     mar = c(0.15, 0.125, 0.15, 0.2),
                     legend = TRUE, leg_sep = 0.1, xpd_hex = 0.75, xpd_leg = 1.5,
                     l = list(x = NULL, y = NULL, default.units = "native"),
                     g = list(x = NULL, y = NULL, default.units = "native"),
                     t = list(label = NULL, default.units = "native"),
                     plot = TRUE)
{

  x  <- as.vector(x)
  y  <- as.vector(y)
  id <- as.vector(id)

  xbnds <- xbnds
  ybnds <- ybnds

  xlim  <- xlim
  ylim  <- ylim

  trans <- trans

  if(max(xbnds) < max(x) | min(xbnds) > min(x)){
    losex <- which(x > max(xbnds) | x < min(xbnds))
    x  <- x[-losex]
    y  <- y[-losex]
    id <- id[-losex]
  }

  if(max(ybnds) < max(y) | min(ybnds) > min(y)){
    losey <- which(y > max(ybnds) | y < min(ybnds))
    x  <- x[-losey]
    y  <- y[-losey]
    id <- id[-losey]
  }

  # Logarithmic distribution if needed ----

  if(length(g$gp) != 0){
    gpl  <- merge_list(as.list(g$gp), list(fill = NA))
    g$gp <- do.call(gpar, gpl)
  }

  l <- merge_list(l, list(default.units = "native"))
  g <- merge_list(g, list(default.units = "native"))
  t <- merge_list(t, list(default.units = "native"))

  if(log == "x" | log == "xy" | log == "yx") {
    xplot <- log10(x)
    xlim  <- log10(xlim)
    xbnds <- log10(xbnds)

    if(length(l$x) != 0){
      if(inherits(l$x, "unit")) {
        if(attributes(l$x)$unit == "native")  l$x <- log10(l$x)
      } else if (!inherits(l$x, "unit") &  l$default.units == "native"){
        l$x <- log10(l$x)
      }
    }

    if(length(g$x) != 0){
      if(inherits(g$x, "unit")) {
        if(attributes(g$x)$unit == "native")  g$x <- log10(g$x)
      } else if (!inherits(g$x, "unit") &  g$default.units == "native"){
        g$x <- log10(g$x)
      }
    }

    if(length(t$x) != 0){
      if(inherits(t$x, "unit")) {
        if(attributes(t$x)$unit == "native")  t$x <- log10(t$x)
      } else if (!inherits(t$x, "unit") &  t$default.units == "native"){
        t$x <- log10(t$x)
      }
    }

  } else {
    xplot <- x
  }

  if(log == "y" | log == "xy" | log == "yx") {
    yplot <- log10(y)
    ylim  <- log10(ylim)
    ybnds <- log10(ybnds)

    if(length(l$y) != 0){
      if(inherits(l$y, "unit")) {
        if(attributes(l$y)$unit == "native")  l$y <- log10(l$y)
      } else if (!inherits(l$y, "unit") &  l$default.units == "native"){
        l$y <- log10(l$y)
      }
    }

    if(length(g$y) != 0){
      if(inherits(g$y, "unit")) {
        if(attributes(g$y)$unit == "native")  g$y <- log10(g$y)
      } else if (!inherits(g$y, "unit") &  g$default.units == "native"){
        g$y <- log10(g$y)
      }
    }

    if(length(t$y) != 0){
      if(inherits(t$y, "unit")) {
        if(attributes(t$y)$unit == "native")  t$y <- log10(t$y)
      } else if (!inherits(t$y, "unit") &  t$default.units == "native"){
        t$y <- log10(t$y)
      }
    }

  } else {
    yplot <- y
  }

  if(xaxs != "i"){
    xlim[which.max(xlim)] <- max(xlim) + (max(xlim) - min(xlim)) * 0.04
    xlim[which.min(xlim)] <- min(xlim) - (max(xlim) - min(xlim)) * 0.04
  }

  if(yaxs != "i"){
    ylim[which.max(ylim)] <- max(ylim) + (max(ylim) - min(ylim)) * 0.04
    ylim[which.min(ylim)] <- min(ylim) - (max(ylim) - min(ylim)) * 0.04
  }

  # Theoretical hexbin object,
  # and reselect the columns of the matrix if needed ----

  hexbin.theo <- hexbin(xplot, yplot, xbins = bins,
                        xbnds = xbnds, ybnds = ybnds,
                        shape = shape)

  if(is.na(select[1]) & length(select) == 1) {

    select <- seq_len(length(unique(id)))

  } else {

    if(inherits(select, "numeric") | inherits(select,"integer")){
      if(any(!(select %in% seq_len(length(unique(id)))))) {
        stop("If 'select' is numeric, it should stand for the index of the id,",
             " i.e. the position in unique(id). So it should be limited to",
             " 1:(amount of different ids)")
      }
    } else {
      reselect  <- match(select, unique(id))
      na.select <- is.na(reselect)
      if(any(na.select)){
        wrong <- select[which(na.select)]
        stop("The following identifier(s) in 'select' is(are)",
             " not found in 'id': ", wrong)
      }
      select <- reselect
    }

    keepid <- id %in% unique(id)[select]

    xplot <- xplot[keepid]
    yplot <- yplot[keepid]

  }

  hexbin.obj <- hexbin(xplot, yplot, xbins = bins,
                       xbnds = xbnds, ybnds = ybnds,
                       shape = shape)

  # Parameters determination ----

  cnt   <- hexbin.obj@count
  tmp   <- hcell2xy(hexbin.obj)

  if(is.na(maxcnt)){
    if(uniform) maxcnt <- max(hexbin.theo@count) else maxcnt <- max(cnt)
  }

  if(is.null(trans)){
    rcnt  <- (cnt - mincnt)/(maxcnt - mincnt)
  } else {
    rcnt <- (trans(cnt) - trans(mincnt))/(trans(maxcnt) -
                                            trans(mincnt))
    if (any(is.na(rcnt))) stop("bad count transformation")

    ucnt <- unique(cnt)

    if(any(as.integer(ucnt + 0.1) != as.integer(inv(trans(ucnt)) + 0.1))){
      stop("The 'inv' function is not the inverse of the 'trans' function'")
    }
  }

  good  <- mincnt <= cnt & cnt <= maxcnt
  xnew  <- tmp$x[good]
  ynew  <- tmp$y[good]

  xbins <- hexbin.theo@xbins
  shape <- hexbin.theo@shape
  sx    <- xbins/diff(hexbin.theo@xbnds)
  sy    <- (xbins * shape)/diff(hexbin.theo@ybnds)

  # Colorcut determination ----

  nc <- length(colorcut)
  if (colorcut[1] > colorcut[nc]) {
    colorcut[1] <- colorcut[1] + 1e-08
    colorcut[nc] <- colorcut[nc] - 1e-08
  } else {
    colorcut[1] <- colorcut[1] - 1e-08
    colorcut[nc] <- colorcut[nc] + 1e-08
  }

  colgrp <- cut(rcnt, colorcut, labels = FALSE)

  # if (any(is.na(colgrp))) colgrp <- ifelse(is.na(colgrp), 0, colgrp)

  clrs <- colramp(length(colorcut) - 1)
  pen  <- clrs[colgrp]

  if(is.null(border)) {
    border     <- pen
    borderline <- T
  }

  # Coordinates of hexagons ----

  inner <- 0.5
  outer <- (2 * inner)/sqrt(3)
  dx    <- inner/sx
  dy    <- outer/(2 * sy)
  hexC  <- hexcoords(dx, dy, sep = NULL)
  n     <- length(cnt)
  n6    <- rep.int(6, n)

  pltx  <- rep.int(hexC$x, n) + rep.int(xnew, n6)
  plty  <- rep.int(hexC$y, n)  + rep.int(ynew, n6)
  plti  <- rep(seq_len(length(pltx)/6), each = 6)

  # Create grob ----

  vp1 <- viewport(name = "vp1")

  vp2 <- viewport(x = mar[2], y = mar[1],
                  width = 1 - mar[4] - mar[2],
                  height = 1 - mar[3] - mar[1],
                  just = c("left", "bottom"),
                  xscale = xlim, yscale = ylim, name = "vp2")

  aGrob <- grobTree()

  if(xaxis){

    if(log == "x" | log == "xy" | log == "yx"){

      xscalet <- 10^xlim

      xticks <- seq_log(xscalet[1], xscalet[2], divide = T)

      if(length(xticks[[1]]) >= 2){

        xtmaj <- xaxisGrob(at = log10(xticks[[1]]), vp = vp2)
        xtmin <- xaxisGrob(at = log10(xticks[[2]]), vp = vp2)

        xtmin$children$labels$label <- ""
        xtmin$children$ticks$y1     <- unit(-0.5 * 0.5, "lines")

        xtmaj$children$labels$label <- xticks[[1]]

        aGrob <- grobTree(aGrob, xtmaj, xtmin)

      } else {

        nxticks <- pretty(xscalet)
        nxticks <- nxticks[!(nxticks < min(xscalet) |
                               nxticks > max(xscalet))]
        xtmaj   <- xaxisGrob(at = log10(nxticks), vp = vp2)
        xtmaj$children$labels$label <- nxticks

        aGrob <- grobTree(aGrob, xtmaj)

      }

    } else {

      xmarks <- pretty(xlim)
      xmarks <- xmarks[xmarks >= min(xlim) & xmarks <= max(xlim)]

      xtmaj <- xaxisGrob(at = xmarks, vp = vp2)

      aGrob <- grobTree(aGrob, xtmaj)

    }

  }

  if(yaxis){

    if(log == "y" | log == "xy" | log == "yx"){

      yscalet <- 10^ylim

      yticks <- seq_log(yscalet[1], yscalet[2], divide = T)

      if(length(yticks[[1]]) >= 2){

        ytmaj <- yaxisGrob(at = log10(yticks[[1]]), vp = vp2)
        ytmin <- yaxisGrob(at = log10(yticks[[2]]), vp = vp2)

        ytmin$children$labels$label <- ""
        ytmin$children$ticks$x1     <- unit(-0.5 * 0.5, "lines")

        ytmaj$children$labels$label <- yticks[[1]]

        aGrob <- grobTree(aGrob, ytmaj, ytmin)

      } else {

        nyticks <- pretty(yscalet)
        nyticks <- nyticks[!(nyticks < min(yscalet) |
                               nyticks > max(yscalet))]
        ytmaj   <- yaxisGrob(at = log10(nyticks), vp = vp2)
        ytmaj$children$labels$label <- nyticks

        aGrob <- grobTree(aGrob, ytmaj)

      }

    } else {

      ymarks <- pretty(ylim)
      ymarks <- ymarks[ymarks >= min(ylim) & ymarks <= max(ylim)]

      ytmaj <- yaxisGrob(at = ymarks, vp = vp2)

      aGrob <- grobTree(aGrob, ytmaj)

    }

  }

  t1 <- textGrob(xlab, y = unit(-3, "line"), gp = gpar(cex = cex), vp = vp2)
  t2 <- textGrob(ylab, x = unit(-3.5, "line"), rot = 90,
                 gp = gpar(cex = cex), vp = vp2)

  vp3 <- viewport(x = mar[2], y = mar[1],
                  width = 1 - mar[4] - mar[2],
                  height = 1 - mar[3] - mar[1],
                  just = c("left", "bottom"),
                  xscale = xlim, yscale = ylim, clip = T)

  p1 <- polygonGrob(x = pltx, y = plty, id=plti, default.units="native",
                    gp=gpar(fill = pen, col = border, lwd = lwd), vp = vp3)

  aGrob <- grobTree(aGrob, p1, t1, t2)


  if(!(length(l$x) == 0 & length(l$y) == 0)) {

    l     <- merge_list(list(vp = vp3), l)
    l1    <- do.call(polylineGrob, l)
    aGrob <- grobTree(aGrob, l1)

  }

  if(!(length(g$x) == 0 & length(g$y) == 0)) {

    g     <- merge_list(list(vp = vp3), g)
    g1    <- do.call(polygonGrob, g)
    aGrob <- grobTree(aGrob, g1)

  }


  if(!(length(g$x) == 0 & length(g$y) == 0 & length(labels) != 0)) {

    t     <- merge_list(list(vp = vp3), t)
    t3    <- do.call(textGrob, t)
    aGrob <- grobTree(aGrob, t3)

  }

  if(box){

    r1 <- rectGrob(vp = vp3,x = 0, y = 0, width = 1, height = 1,
                   just = c("left", "bottom"), gp = gpar(fill = NA))

    aGrob <- grobTree(aGrob, r1)

  }

  t4 <- textGrob(main, x = 0.5, y = 1 - mar[3]/2,
                 gp = gpar(font = 2, cex = 1.2 * cex), vp = vp1)

  aGrob <- grobTree(aGrob, t4)

  # Add a legend (former function) ----

  if(legend){

    legend.vp <- viewport(x = 1, y = 0, width = mar[4], height = 1,
                          just = c("right", "bottom"), name = "vpl")

    inter <- 0.4 * xpd_leg * (length(colorcut)-1)/30
    y1       <- 0.5 - inter
    y2       <- 0.5 + inter
    xmid     <- 0.5
    unit     <- "npc"
    xcor     <- shape * 1/(mar[4])

    if (is.null(trans)) {
      sc <- maxcnt - mincnt
      bnds <- round(mincnt + sc * colorcut)
    } else {
      if (!is.function(trans) && !is.function(inv))
        stop("'trans' and 'inv' must both be functions if 'trans' is not NULL")
      con  <- trans(mincnt)
      sc   <- trans(maxcnt) - con
      bnds <- round(inv(con + sc * colorcut))
    }

    n       <- length(colorcut)
    spacing <- abs(y2 - y1)/(n)
    inner   <- xpd_hex * (sqrt(3) * spacing)/2

    dx   <- inner/2
    dy   <- dx/sqrt(3)
    hexC <- hexcoords(dx, dy, n = 1, sep = NULL)

    textx <- xmid + leg_sep/2
    tx    <- (hexC$x*xcor) + xmid - max(hexC$x*xcor) - leg_sep/2

    for (i in seq(length = n - 1)) {

      if(borderline) border <- colramp()[i] else border <- border

      pl <- polygonGrob(x = tx,
                        y = min(c(y1, y2)) + hexC$y + i * spacing,
                        default.units = unit, id.lengths = 6,
                        gp = gpar(fill = colramp()[i], col = border, lwd = lwd),
                        vp = legend.vp)

      tl1 <- textGrob(as.character(bnds[i]),
                      x = textx,
                      y = min(c(y1, y2)) + (i - 0.5) * spacing,
                      default.units = unit, just = "left", gp = gpar(cex = cex),
                      vp = legend.vp)

      aGrob <- grobTree(aGrob, pl, tl1)
    }

    tl2 <- textGrob(as.character(bnds[n]), textx,
                    min(c(y1, y2)) +  (n - 0.5) * spacing,
                    default.units = unit, just = "left", gp = gpar(cex = cex),
                    vp = legend.vp)

    aGrob <- grobTree(aGrob, tl2)

  }

  if(plot){
    grid.newpage()
    grid.draw(aGrob)
  } else {
    return(aGrob)
  }

}

