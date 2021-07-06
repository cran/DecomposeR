#' @title Correlation of time-series with different sampling rate
#'
#' @description Allows to correlate time-series having different sampling rate,
#' if they have a comparable depth or time scale
#'
#' @param xy1 intensity values for the first data set
#' @param dt1 depth or time scale for the first data set
#' @param xy2 intensity values for the second data set
#' @param dt2 depth or time scale for the second data set
#' @param plot whether to plot
#' @param output whether to output
#' @param type type of points in the plot (see help page of \code{lines()} for
#' details)
#' @param ... additional parameters to feed to the \code{lines()} function
#'
#' @return a list of correlation ($cor), slope ($slope), intercept ($intercept)
#' (two values for each: interpolation to fit dt1 and dt2 respectively), and of
#' the xy1 and xy2 values, interpolated for dt1 ($df1) and df2 ($df2)
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
#' xy.pure <- (1 + 0.6 * sin(t*2*pi/p2)) * sin(t*2*pi/p1)  + 2 * sin(t*2*pi/p2)
#'
#' xy <- xy.pure + rnorm(n, sd = 0.5)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt.pure <- cumsum(inter_dt)
#'
#' keep <- runif(length(dt.pure)) < 0.5
#'
#' xy <- xy[keep]
#' dt <- dt.pure[keep] + rnorm(sum(keep), -0.2, 0.2)
#'
#' par(mfrow = c(1,2))
#'
#' plot(xy, dt, type = "o", pch = 19)
#'
#' plot(xy.pure, dt.pure, type = "o", pch = 19)
#'
#' par(mfrow = c(1,1))
#'
#' out <- approx.cor(xy, dt, xy.pure, dt.pure)
#'
#' out$cor
#' out$slope
#' out$intercept
#'
#' @export

approx.cor <- function(xy1, dt1, xy2, dt2,
                       plot = T, output = T, type = "p", ...)
{

  rem1 <- is.na(xy1) | is.na(dt1)
  rem2 <- is.na(xy2) | is.na(dt2)

  xy1p <- xy1[!rem1]
  dt1p <- dt1[!rem1]
  xy2p <- xy2[!rem2]
  dt2p <- dt2[!rem2]

  if(is.unsorted(dt1p) & is.unsorted(rev(dt1p))) {
    stop("The 'dt1' parameter should be sorted")
  }

  if(is.unsorted(dt2p) & is.unsorted(rev(dt2p))) {
    stop("The 'dt2' parameter should be sorted")
  }

  if(!is.unsorted(rev(dt1p))){
    xy1p <- rev(xy1p)
    dt1p <- rev(dt1p)
  }

  if(!is.unsorted(rev(dt2p))){
    xy2p <- rev(xy2p)
    dt2p <- rev(dt2p)
  }

  ddt1p <- duplicated(dt1p)
  ddt2p <- duplicated(dt2p)

  if(any(ddt1p)){
    warning(paste("Duplicated dt1 values: ",
                  paste(signif(dt1p[which(ddt1p)], 7), collapse = ", ")))
  }

  if(any(ddt2p)){
    warning(paste("Duplicated dt2 values: ",
                  paste(signif(dt2p[which(ddt2p)], 7), collapse = ", ")))
  }

  xy1i <- approx(dt1p, xy1p, dt2p)$y
  xy2i <- approx(dt2p, xy2p, dt1p)$y

  df1 <- data.frame(dt = dt1p, xy1 = xy1p, xy2 = xy2i)
  df2 <- data.frame(dt = dt2p, xy2 = xy2p, xy1 = xy1i)


  cor.out <- c(cor(xy2i, xy1p, use = "complete.obs"),
               cor(xy1i, xy2p, use = "complete.obs"))

  coefa <- lm(xy2i ~ xy1p)$coef
  coefb <- lm(xy2p ~ xy1i)$coef

  if(isTRUE(plot)){

    plot.new()
    plot.window(xlim = range(xy1p, na.rm = T), ylim = range(xy2p, na.rm = T))

    axis(1)
    axis(2)

    title(xlab = "xy1", ylab = "xy2")

    box()

    lines(xy1i, xy2p, type = type, ...)
    lines(xy1p, xy2i, type = type, ...)

    intercept <- c(coefa[1], coefb[1])
    slope     <- c(coefa[2], coefb[2])

    names(intercept) <- NULL
    names(slope)     <- NULL

    abline(a = intercept[1], b = slope[1], col = "red")
    abline(a = intercept[2], b = slope[2], col = "red")

  }

  output.list <- list(cor = cor.out,
                      slope = slope, intercept = intercept,
                      d1 = df1, d2 = df2)

  if(isTRUE(output)) return(output.list)

}


