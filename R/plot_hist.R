#' @title Group and/or log-scale histogram
#'
#' @description Specialised histogram: allows to work in log-scale (for x) and
#' to distinguish different groups of data
#'
#' @param x vector or matrix
#' @param breaks one of:
#' \itemize{
#' \item{a vector giving the breakpoints between histogram cells,}
#' \item{a function to compute the vector of breakpoints,}
#' \item{a single number giving the number of cells for the histogram,}
#' \item{a character string naming an algorithm to compute the number of cells
#' (see ‘Details’ in \code{\link{hist}}),}
#' \item{a function to compute the number of cells.}
#' }
#' In the last three cases the number is a suggestion only; as the breakpoints
#' will be set to pretty values, the number is limited to 1e6 (with a warning if
#' it was larger). If breaks is a function, the x vector is supplied to it as
#' the only argument (and the number of breaks is only limited by the amount of
#' available memory).
#' @param id a vector of ids for each x value, to separate different groups of
#' data
#' @param select a vector of id values idenifying the groups of data to plot and
#' their order
#' @param pile whether to cumulate the different one on the other
#' @param line whether to plot as lines or rectangles
#' @param mids if lines is TRUE, whether the nodes of the lines are the middle
#' positions or the upper corner of the rectangles.
#' @param xlim,ylim the boundaries for the plots. If ylim = NA the upper ylim
#' will be increased by 10\% to allow for text (see 'text' parameter)
#' @param xlog whether to set the x axis in log scale
#' @param axes whether to plot the axes
#' @param xa,ya list of arguments to feed minorAxis for the x and y axes
#' respectively
#' @param main,xlab,ylab the main title and the labels of the x and y axes
#' @param col a function or a character vector defining the colors of the
#' different modes
#' @param border the colour of the borders, by default identical to col
#' @param text if there are different groups, whether to add a number above
#' each of them to distinguish them
#' @param labels the labels to put on top of each group
#' @param t a list of parameters to feed text()
#' @param add whether to add the plot to a preexisting plot
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
#'
#' dec <- extricate(xy, dt, nimf = 7, sifting = 10,
#'                  repl = 10, comb = 10, factor_noise = 10,
#'                  speak = TRUE)
#' \dontrun{
#' plot_emd(dec, dir = tempdir())}
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' ht  <- inst.pulse(dec, plot = FALSE)
#'
#' opar <- par('mfrow')
#'
#' par(mfrow = c(2,1))
#'
#' plot_hist(x = 1/ht$f, breaks = 500,
#'           xlog = TRUE, xlab = "Period")
#'
#' plot_hist(x = 1/ht$f, breaks = 500, id = ht$mode,
#'           xlog = TRUE, text = TRUE, add = TRUE, line = TRUE, pile = FALSE)
#'
#' abline(v = c(p1, p2), col = "red", lwd = 2, lty = 5)
#'
#' plot_hist(x = 1/ht$f, breaks = 500, id = ht$mode,
#'           xlog = TRUE, text = TRUE, xlab = "Period")
#'
#' abline(v = c(p1, p2), col = "red", lwd = 2, lty = 5)
#'
#' par(mfrow = opar)
#'
#' @importFrom StratigrapheR merge_list minorAxis
#' @export

plot_hist <- function(x, breaks = 100, id = NA,
                      select = NA, pile = TRUE, line = FALSE, mids = FALSE,
                      xlim = NA, ylim = NA, xlog = FALSE, axes = TRUE,
                      xa = list(), ya = list(), main = "",
                      xlab = "X", ylab = "Counts", col = NA, border = NA,
                      text = FALSE, labels = NA,
                      t = list(adj = c(0.5, -2), font = 2),
                      add = FALSE)
{

  if(!inherits(breaks, "numeric") | !inherits(breaks, "integer") &
     length(breaks) != 1){
    stop("'breaks' should be a numeric or integer of length 1")
  }

  if(!(isFALSE(pile) | isTRUE(pile))) {
    stop("The 'pile' parameter should be TRUE or FALSE'")
  }

  if(!(isFALSE(line) | isTRUE(line))) {
    stop("The 'line' parameter should be TRUE or FALSE'")
  }

  if(!(isFALSE(mids) | isTRUE(mids))) {
    stop("The 'mids' parameter should be TRUE or FALSE'")
  }

  x  <- as.vector(x)
  id <- as.vector(id)

  if(!(length(id) == length(x) | (is.na(id[1]) & length(id) == 1))) {
    stop("'id' should be of same length than 'x' or simply 'NA'")
  }

  if(is.na(select[1]) & length(select) == 1) {

    select <- seq_len(length(unique(id)))

  } else {

    if(class(select) == "numeric" | class(select) == "integer"){
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

  }

  nselect  <- length(unique(as.vector(id)))

  if(class(breaks) == "function") breaks <- breaks(x)

  if(xlog & length(breaks) == 1){
    all      <- hist(log10(x), breaks = breaks, plot = F)
    all      <- hist(x, breaks = 10^all$breaks, plot = F)
  } else {
    all <- hist(x, breaks = breaks, plot = F)
  }

  if(xlog) log = "x" else log = ""

  if(length(id) == length(x)) {

    separate <- T

    ix <- split(x, id)

    rx <- lapply(ix[rev(select)],
                 function(x, plot, breaks) {
                   hist(x, plot = plot, breaks = breaks)$counts
                 },
                 plot = F, breaks = all$breaks)

    rx <- matrix(unlist(rx), ncol = length(select))

    rcounts <- rx

    if(pile) {
      rcounts <- rcounts[,rev(seq_len(ncol(rcounts))), drop = F]

      if(length(select) == 1) {
        rcounts <- matrix(apply(rcounts, 1, cumsum))
      } else {
        rcounts <- t(apply(rcounts, 1, cumsum))
      }

      rcounts <- rcounts[,rev(seq_len(ncol(rcounts))), drop = F]
      if(text) maxhist <- rcounts[,1]
    }

  } else {
    separate <- F
  }

  if(is.na(xlim[1]) & length(xlim) == 1) {
    xlim <-  c(min(all$breaks), max(all$breaks))
  }

  if(is.na(ylim[1]) & length(ylim) == 1) {
    ylim <-  c(0, max(all$counts) * 1.1)
  }

  if(is.na(col)){
    if(separate) col <- matlab.like(length(unique(id))) else col <- "black"
  }

  if(separate) {
    col <- homogenise(unique(as.vector(id)), l = list(rev(col)))[[1]]
    col <- rev(col)[rev(select)]
  }

  if(text & separate){

    if(is.na(labels[1]) & length(labels) == 1) {
      labels <- unique(id)
    } else {
      labels <- homogenise(unique(as.vector(id)),
                           l = list(labels = labels), cycle = F)[[1]]
    }

    labels <- labels[select]

    irx <- apply(rx, 2, which.max)

    if(pile){
      tt <- merge_list(list(x = all$mids[irx],
                            y = maxhist[irx]),
                       t, list(adj = c(0.5, -2), font = 2,
                               col = col, labels = rev(labels)))
    } else {
      tt <- merge_list(list(x = all$mids[irx],
                            y = rx[matrix(c(irx,seq_len(length(irx))),
                                          ncol = 2)]),
                       t, list(adj = c(0.5, -2), font = 2,
                               col = col, labels = rev(labels)))
    }

  }

  if(!add){

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, log = log)

    if(axes){
      if(xlog){

        xt <- seq_log(10^par("usr")[1], 10^par("usr")[2], divide = T)

        if(length(xt[[1]]) >= 2){
          xae <- merge_list(list(side = 1), xa,
                            list(at.maj = xt[[1]], at.min = xt[[2]]))
        } else {

          xscalet <- 10^par("usr")[c(1,2)]
          nxticks <- axTicks(side = 1, log = T)
          nxticks <- nxticks[!(nxticks < min(xscalet) |
                                 nxticks > max(xscalet))]

          xae <- merge_list(list(side = 1), xa,
                            list(at.maj = nxticks))

        }

        do.call(minorAxis, xae)

      } else {
        xae <- merge_list(list(side = 1), xa)
        do.call(minorAxis, xae)
      }
      yae <- merge_list(list(side = 2), ya)
      do.call(minorAxis, yae)
    }

    title(main = main, xlab = xlab, ylab = ylab)

  }

  if(separate){

    if(line){

      vcounts <- as.vector(rcounts)

      if(mids){

        lid <- as.vector(matrix(rep(select, nrow(rcounts)),
                                ncol = ncol(rcounts), byrow = T))

        multilines(lid, rep(all$mids, length(select)),vcounts, col = col)

      } else {

        lv <- length(vcounts)

        l <- as.vector(rep(all$breaks[-length(all$breaks)], length(select)))
        r <- as.vector(rep(all$breaks[-1], length(select)))

        py <- rep(vcounts, 2)[seq_mult(2*lv, lv)]
        px <- c(l, r)[seq_mult(2*lv, lv)]

        lid <- as.vector(matrix(rep(select, 2*nrow(rcounts)),
                                ncol = ncol(rcounts), byrow = T))

        multilines(lid, px, py, col = col)

      }

    } else {

      for(i in seq_len((length(select))))
      {

        cur.hist        <- all
        cur.hist$counts <- rcounts[,i]

        suppressWarnings(plot(cur.hist, add = T, freq = T,
                              col = col[i], border = border))
      }

    }

    if(text) do.call("text", tt)

  } else {

    if(line){

      val <- all$counts

      if(mids){

        lines(all$mids, val)

      } else {

        l <- all$breaks[-length(all$breaks)]
        r <- all$breaks[-1]

        py <- rep(val, 2)[seq_mult(2*length(val), length(val))]
        px <- c(l, r)[seq_mult(2*length(val), length(val))]

        lines(px, py)

      }

    } else {

      if(is.na(border) & length(border) == 1) border <- col

      suppressWarnings(plot(all, freq = T, add = T, col = col,
                            border = border))

    }

  }

}





