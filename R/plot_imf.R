#' @title Plot IMFs characteristics
#'
#' @description General plot for the envelope, instantaneous frequency (period)
#' and identity tuning of an intrinsic mode function (IMF)
#'
#' @param pulse a pulse object
#' @param dtlim,xylim,flim,fclim the boundaries for the plots, respectively for
#' the depth/time, amplitude, frequency and freqeuncy carrier
#' @param dtline,fline,fcline coordinates to add vertical/horizontal lines
#' @param vertical whether to have the depth/time [dt] axis vertically
#' @param n the the number of intervals defined by minor ticks
#' (geologist convention) or horizontaly (climatologist convention)
#' @param at.maj the positions at which major tick-marks are to be drawn.
#' @param ls,le1,le2,lid,lcos lists of parameters to feed lines, for the
#' original signal, the upper and lower envelope, the identity tuning, and the
#' cosine line in the identity tuning
#' @param ldt,lf,lfc lists of parameters to provide the abline function (makes
#' personalised lines for you to have a better grasp of the data).
#' @param box whether to draw boxes around the plots
#'
#' @details the line in the identity tuning plot is a genuine cosine,
#' independent from the signal. This is evident when riding waves generate
#' dephasing.
#'
#' @examples
#' n <- 600
#'
#' t <- seq_len(n)
#'
#' p1 <- 30
#' p2 <- 40 * 21
#'
#' am <- sin(t*2*pi/p2 + 50) + 0.03
#'
#' xy <- sin(t*2*pi/p1 + 50) * 3 * am
#'
#' int <- c(rep(1, 99 + 100), seq(1,3,2/100), seq(3,1,-2/100), rep(1,100 + 99))
#'
#' dt <- cumsum(int)
#'
#' samp <- approx(dt, xy, xout = seq(1,802, by = 2))
#'
#' xy <- samp$y
#' dt <- samp$x
#'
#' e <- normalise(m = xy, dt = dt)$a
#'
#' cond <- dt < 75
#'
#' xy <- xy[!cond]
#' dt <- (dt[!cond] - 75) / 1.2
#' e  <- e[!cond]
#'
#' dq   <- dq.algorithm(xy/e, dt)
#'
#' pulse <- as.pulse(dt = dt, m = xy, f = dq$f, a = e, idt = dq$idt,
#'                   repl = 1)
#'
#' plot_imf(pulse, fline = 25, dtline = c(222, 489))
#'
#' @importFrom StratigrapheR merge_list minorAxis minorAxisTicks
#' @import graphics
#' @export

plot_imf <- function(pulse,
                     dtlim = NULL, xylim = NULL, flim = NULL, fclim = NULL,
                     dtline = NULL, fline = NULL, fcline = NULL,
                     vertical = FALSE, n = 10, at.maj = NULL,
                     ls = list(type = "o", pch = 19),
                     le1 = list(lwd = 2), le2 = list(lty = 2),
                     lid = list(type = "p", pch = 19), lcos = list(),
                     ldt = list(lty = 5, lwd = 2),
                     lf = list(lty = 5), lfc = list(lty = 5),
                     box = TRUE)
{

  if(!is.pulse(pulse)) stop("The 'pulse' parameter should be a pulse object, ",
                            "see ?is.pulse for more details")

  plots <- c("s")

  if(!is.null(pulse$a))   plots <- c(plots, "a")
  if(!is.null(pulse$f))   plots <- c(plots, "f")
  if(!is.null(pulse$idt)) plots <- c(plots, "idt")

  nc <- length(unique(pulse$repl))

  if(is.null(dtlim)) dtlim <- c(min(pulse$dt), max(pulse$dt))

  if(is.null(xylim)) xylim <- c(min(pulse$m), max(pulse$m))
  if(is.null(flim))  flim  <- c(min(1/pulse$f), max(1/pulse$f))
  if(is.null(fclim))  fclim  <- c(-1,1)

  # Graphical part ----

  if(vertical){
    dtline1 <- merge_list(list(h = dtline), ldt)
    fcline1  <- merge_list(list(v = fcline), lfc)
  } else {
    dtline1 <- merge_list(list(v = dtline), ldt)
    fcline1  <- merge_list(list(h = fcline), lfc)
  }

  opar <- par('mfrow')

  lp <- length(plots)

  if(vertical) par(mfrow = c(1, lp)) else par(mfrow = c(lp, 1))

  # Plot signal and envelope ----

  plot.new()

  if(vertical) {
    plot.window(xlim = xylim, ylim = dtlim)
    title(xlab = "xy", ylab = "dt")
    minorAxis(2, n = n, at.maj = at.maj, las = 1)
    axis(1)
  } else {
    plot.window(xlim = dtlim, ylim = xylim)
    title(xlab = "dt", ylab = "xy")
    minorAxis(1, n = n, at.maj = at.maj)
    axis(2)
  }

  if(box) box()

  repl <- seq_len(nc)

  for(i in repl){

    if(vertical) {
      ps1 <- merge_list(list(y = pulse$dt, x = pulse$m[,i]), ls)
      pe1 <- merge_list(list(y = pulse$dt, x = pulse$a[,i]), le1)
      pe2 <- merge_list(list(y = pulse$dt, x = -pulse$a[,i]), le2)
    } else {
      ps1 <- merge_list(list(x = pulse$dt, y = pulse$m[,i]), ls)
      pe1 <- merge_list(list(x = pulse$dt, y = pulse$a[,i]), le1)
      pe2 <- merge_list(list(x = pulse$dt, y = -pulse$a[,i]), le2)
    }

    do.call('lines', ps1)
    do.call('lines', pe1)
    do.call('lines', pe2)

  }

  do.call(abline, dtline1)

  # Plot frequency carrier ----

  if("a" %in% plots){

    plot.new()

    if(vertical){
      plot.window(xlim = fclim, ylim = dtlim)
      title(ylab = "dt")
      minorAxis(2, n = n, at.maj = at.maj, las = 1)
      axis(1)
    } else {
      plot.window(xlim = dtlim, ylim = fclim)
      title(xlab = "dt")
      minorAxis(1, n = n, at.maj = at.maj)
      axis(2)
    }

    if(box) box()


    for(i in repl){

      if(vertical){
        ps2 <- merge_list(list(y = pulse$dt, x = pulse$m[,i]/pulse$a[,i]), ls)
      } else {
        ps2 <- merge_list(list(x = pulse$dt, y = pulse$m[,i]/pulse$a[,i]), ls)
      }

      do.call('lines', ps2)

    }

    do.call(abline, dtline1)
    do.call(abline, fcline1)

  }

  # Plot IF period ----

  if("f" %in% plots){

    plot.new()

    if(vertical){
      plot.window(xlim = flim, ylim = dtlim, log = "x")
      title(xlab = "Period", ylab = "dt")
      minorAxis(2, n = n, at.maj = at.maj, las = 1)
      axis(1)
      fline1 <- merge_list(list(v = fline), lf)
    } else {
      plot.window(xlim = dtlim, ylim = flim, log = "y")
      title(xlab = "dt", ylab = "Period")
      minorAxis(1, n = n, at.maj = at.maj)
      axis(2)
      fline1 <- merge_list(list(h = fline), lf)
    }

    if(box) box()

    for(i in repl){

      if(vertical){
        ps3 <- merge_list(list(y = pulse$dt, x = 1/pulse$f[,i]), ls)
      } else {
        ps3 <- merge_list(list(x = pulse$dt, y = 1/pulse$f[,i]), ls)
      }

      do.call('lines', ps3)

    }

    do.call(abline, fline1)
    do.call(abline, dtline1)

  }

  # Plot identity tuning ----

  if("idt" %in% plots){

    if(is.null(dtlim)){
      c(min(pulse$idt), max(pulse$idt))
    } else {
      idtlim <- approx(y = pulse$idt[,1], x = pulse$dt, xout = dtlim)$y

      if(is.na(idtlim[1])) idtlim[1] <- min(pulse$idt)
      if(is.na(idtlim[2])) idtlim[2] <- max(pulse$idt)
    }

    plot.new()

    if(vertical){
      plot.window(xlim = fclim, ylim = idtlim)
      title(ylab = "dt")
      axis(4, las = 1)
      axis(4, at = c(1), las = 1)
      axis(1)
    } else {
      plot.window(xlim = idtlim, ylim = fclim)
      title(xlab = "Cycle Number")
      axis(2)
      axis(1, at = c(1))
      axis(1)
    }

    if(box) box()

    ints <- seq(0, ceiling(max(pulse$idt)), by = 0.01)

    if(vertical){
      pcos <- merge_list(list(y = ints, x = cos(2*pi*ints)), lcos)
    } else {
      pcos <- merge_list(list(x = ints, y = cos(2*pi*ints)), lcos)
    }

    do.call('lines', pcos)

    # MAKE SAMPLE RATE AGAIN

    for(i in repl) {

      if(vertical){
        pid <- merge_list(list(y = pulse$idt[,i],
                               x = pulse$m[,i]/pulse$a[,i]), lid)
      } else {
        pid <- merge_list(list(x = pulse$idt[,i],
                               y = pulse$m[,i]/pulse$a[,i]), lid)
      }

      do.call('lines', pid)

    }

    if(nc == 1){

      mt <- minorAxisTicks(usr = c(dtlim[1] - 0.04*diff(dtlim),
                                   dtlim[2] + 0.04*diff(dtlim)),
                           n = n, at.maj = at.maj, extend = F)

      mt1 <- approx(y = pulse$idt[,1], x = pulse$dt, xout = mt[[1]])
      mt2 <- approx(y = pulse$idt[,1], x = pulse$dt, xout = mt[[2]])

      exc <- !is.na(mt1$y)

      if(vertical) {
        minorAxis(2, at.maj = mt1$y[exc], labels.maj = mt1$x[exc], at.min = mt2$y,
                  las = 1)
      } else {
        minorAxis(3, at.maj = mt1$y[exc], labels.maj = mt1$x[exc], at.min = mt2$y)
      }

      do.call(abline, fcline1)

      if(!is.null(dtlim)){

        mt3     <- approx(y = pulse$idt[,1], x = pulse$dt, xout = dtline)

        if(vertical) {
          dtline2 <- merge_list(list(h = mt3$y), ldt)
        } else {
          dtline2 <- merge_list(list(v = mt3$y), ldt)
        }

        do.call(abline, dtline2)

      }
    }
  }

  par(mfrow = opar)

}




