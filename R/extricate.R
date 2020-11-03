#' @title Extricate a signal: an EEMD algorithm
#'
#' @description Performes EEMD
#'
#' @param xy signal, detrended and demeaned, maybe linearly interpolated to have
#' regular sampling interval
#' @param dt depth/time
#' @param nimf number of modes/components/intrinsic mode functions to decompose
#' the signal into
#' @param ini an optional vector of length n of the eventual initial Intrinsic
#' Mode Function xy would be a demodulation of, if it is a demodulation. In that
#' case the mode indexes will start at 2.
#' @param repl the amount of decompositions to output
#' @param comb the amount of decompositions each output decomposition will be a
#' combination of. Has to be a multiple of 2 (even and odd extension stacks
#' have to be combined in any case)
#' @param mirror_noise whether to generate a mirrored noise signal (for even and
#' odd extension) that will cancel perfectly when combining the decompositions
#' @param factor_noise a factor for the amplitude of white noise (finite
#' amplitude obtained via \code{\link{runif}}). By default it will be
#' multiplied with the mean of the lagged-one difference to define the noise
#' amplitude
#' @param unit_noise whether to multiply factor_noise by the mean of the
#' lagged-one difference (unit_noise = "1stdiff") or not (unit_noise =
#' "native")
#' @param sifting amount of iterations of the sifting process
#' @param output_sifting whether to output each sifting
#' @param remove whether to remove the linear trend (remove = "trend") or the
#' mean (remove = "mean") prior to decomposition. The removed part will be
#' added back after the decomposition. If remove is anything else, nothing will
#' be removed, which can be problematic for the even and odd extension scheme
#' used.
#' @param bind whether to bind the removed trend or mean to the last
#' component (T), or to add it as another component (F)
#' @param speak whether to print a sentence at each sifting: it gives the stack
#' (even or odd), the mode number and sifting number
#' @param plot_process whether to have a plot of the entire sifting process.
#' This slows down the algorithm, use with low 'repl' and 'comb' values for
#' visualisation purposes
#' @param pdf whether the plot be directly set as a pdf file
#' @param name,ext,dir,width,height,track,openfile arguments to provide to
#' pdfDisplay if plot_process and pdf are TRUE
#'
#' @return a list made of $xy (original signal), $dt (depth/time), $m (a matrix
#' of the decomposition), $repl (the replicate id of each point) and
#' $mode (the mode id of each point). If output_sifting is TRUE, additional
#' $even_sifting and $odd_sifting data.tables are provided, giving the
#' condensed siftings for the even and odd extensions.
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
#' dec <- extricate(xy, dt, nimf = 7, repl = 1, comb = 40, factor_noise = 10,
#'                  sifting = 10, speak = TRUE, output_sifting = TRUE)
#'
#' integrity(xy, dec)
#'
#' parsimony(dec)
#'
#' plot_emd(dec, select = c(4, 6), pdf = FALSE)
#' \dontrun{
#' plot_emd(dec, li = list(v = 0), dir = tempdir())}
#'
#' @export
#' @importFrom dplyr lag
#' @importFrom StratigrapheR homogenise seq_mult multilines pdfDisplay

extricate <- function(xy, dt, nimf, ini = NULL, repl = 1, comb = 100,
                      mirror_noise = TRUE, factor_noise = 3,
                      unit_noise = "1stdiff", sifting  = 1, output_sifting = FALSE,
                      remove = "trend", bind = FALSE, speak = FALSE, plot_process = FALSE,
                      pdf = TRUE, name = "extricate", ext = ".pdf", dir = tempdir(),
                      width = 10, height = 20, track = TRUE, openfile = TRUE)
{
  # Conditions ----

  repl  <- as.integer(repl)
  comb  <- as.integer(comb)

  sifting <- homogenise(seq_len(nimf), l = list(sifting = sifting))$sifting

  if(comb/2 != round(comb/2)) stop("The comb parameter should be a multiple of 2")

  xy <- as.vector(xy)

  lx <- length(xy)

  rep <- repl * comb

  l_number <- 2 + rep/2

  iteration_size <- l_number * lx

  if(iteration_size > 10^7) {
    warning(paste("The signal you want to decompose is made of", lx,
                  "data points. It will be repeated", l_number,
                  "times (2 + (repl + comb)/2), making an object of",
                  iteration_size, "data points. This is a lot, you may",
                  "want to iterate this function with less repetitions and",
                  "different seeds."))

  }

  if(remove == "mean") {
    trend   <- mean(xy)
    xyo     <- xy
    xy      <- xy - trend
  } else if(remove == "trend"){
    lm.0  <- lm(xy ~ dt)
    trend <- lm.0$coeff[2] * dt + lm.0$coeff[1]
    xyo   <- xy
    xy    <- xy - trend
  }

  # Pile Up ----

  stack <- pile.up(xy, dt, rep, warn = F)

  sil <- length(xy)
  stl <- nrow(stack)

  # Prepare noise ----

  stack_length <- nrow(stack)

  if(unit_noise == "native"){
    nlevel <- factor_noise
  } else if(unit_noise == "1stdiff"){
    mean_diff <- mean(abs(xy - lag(xy))[-1])
    nlevel <- factor_noise * mean_diff/2
  } else {
    stop("'unit_noise' should be 'native' or '1stunit'.")
  }

  noise1 <- runif(stack_length, -nlevel, nlevel)

  if(mirror_noise){

    noise2 <- noise1 * -stack$invert

  } else {

    noise2 <- runif(stack_length, -mean_diff/2, mean_diff/2)

  }

  # Prepare constants ----

  depth   <- stack$ndt

  repl_seq <- seq_len(repl)

  # Graphical wrapper of the general process ----

  g <- function()
  {

    if(plot_process){

      opar <- par("mfrow")
      par(mfrow = c(1,2))
      on.exit(par(mfrow = opar))

      col.cy <- c("tan1" ,
                  homogenise(n = repl*comb/2,
                             l = list(col = c("grey10", "grey40")))$col,
                  "tan1")

      lty.cy <- c(1 ,
                  homogenise(n = repl*comb/2,
                             l = list(lty = c(5, 1)))$lty,
                  1)

      plot(stack$even, depth, type = "n")
      multilines(stack$id, stack$even, depth, col = col.cy, lty = lty.cy)

      plot(stack$odd, depth, type = "n")
      multilines(stack$id, stack$odd, depth, col = col.cy, lty = lty.cy)

    }

    # i <- 1

    for(i in 1:2)
    {

      # Work with even and odd repetitions

      if(i == 1) {
        subject <- stack$even
        noise   <- noise1
        even    <- T
      } else if(i == 2) {
        subject <- stack$odd
        noise   <- noise2
        even <- F
      }

      subject <- subject + noise

      accu <- data.frame(matrix(nrow = sil, ncol = 0))

      if(output_sifting) accu_s <- data.frame(matrix(nrow = sil, ncol = 0))

      # j <- 1

      for(j in seq_len(nimf))
      {

        # Work with each IMF

        proto_imf <- subject

        # s <- 1

        for(s in seq_len(sifting[j]))
        {

          # Work with each sifting

          if(i == 1) odden <- "even" else odden <- "odd"

          sid <- paste("Stack: ",odden,", Mode: ", j, ", Sifting: ", s,
                       sep = "")

          if(speak) print(sid)

          ext <- extremist(proto_imf, bound = T, zc = F)

          cond <- nrow(ext$maxindex) + nrow(ext$maxindex) - 4

          maxima_index <- c(ext$maxindex$l,ext$maxindex$r)
          l_max        <- length(maxima_index)/2
          maxima_index <- unique(maxima_index[kronecker(1:l_max, c(0,l_max ), "+")])

          minima_index <- c(ext$minindex$l,ext$minindex$r)
          l_min        <- length(minima_index)/2
          minima_index <- unique(minima_index[kronecker(1:l_min, c(0,l_min ), "+")])

          max_spline <- spline(y = proto_imf[maxima_index],
                               x = depth[maxima_index],
                               method = "fmm",
                               xout = depth)$y

          min_spline <- spline(y = proto_imf[minima_index],
                               x = depth[minima_index],
                               method = "fmm",
                               xout = depth)$y

          removed <- (min_spline + max_spline)/2

          if(plot_process){

            plot(proto_imf, depth, type = "n", main = "Sifting",
                 xlab = "xy", ylab = "dt")
            multilines(stack$id, proto_imf, depth, col = col.cy, lty = lty.cy)
            lines(max_spline, depth, col = "red")
            lines(min_spline, depth, col = "blue")
            lines(removed, depth, col = "green")

          }

          proto_imf <- proto_imf - removed

          if(plot_process) {

            plot(proto_imf, depth, type = "n", main = sid,
                 xlab = "xy", ylab = "dt")
            multilines(stack$id, proto_imf, depth, col =  col.cy, lty = lty.cy)

          }

          if(output_sifting) {

            condensed_s <- pile.down(x = proto_imf, stack = stack,
                                     even = even, n = comb/2)

            nfilled_s <- ncol(accu_s)
            cname_s   <- paste("M", j, "S", s, "R", repl_seq, sep = "")

            accu_s <- cbind(accu_s, condensed_s)

            colnames(accu_s)[nfilled_s + repl_seq] <- cname_s

          }

        }

        condensed <- pile.down(x = proto_imf, stack = stack,
                               even = even, n = comb/2)

        nfilled <- ncol(accu)
        cname   <- paste("M", j, "S", s, "R", repl_seq, sep = "")

        accu <- cbind(accu, condensed)

        colnames(accu)[nfilled + repl_seq] <- cname

        subject <- subject - proto_imf

      }

      residue <- pile.down(x = subject, stack = stack,
                           even = even, n = comb/2)

      nfilled <- ncol(accu)
      cname   <- paste("Residue_R", repl_seq, sep = "")

      accu <- cbind(accu, residue)

      colnames(accu)[nfilled + repl_seq] <-cname

      if(i == 1) even_emd <- accu else if(i == 2) odd_emd <- accu

      if(output_sifting){
        if(i == 1) even_emd_s <- accu_s else if(i == 2) odd_emd_s <- accu_s
      }

    }

    # Mean even and odd ----

    mean <- (even_emd + odd_emd)/2

    if(output_sifting) {
      res <- list(mean = mean, even_emd_s = even_emd_s, odd_emd_s = odd_emd_s)
    } else {
      res <- list(mean = mean)
    }

    if(plot_process) par(mfrow = opar)

    return(res)

  }

  # ----

  if(plot_process & pdf){

    res <- pdfDisplay(g(), name = name, ext = ext, dir = dir,
                      width = width, height = height, track = track,
                      openfile = openfile, output = T)

  } else {

    res <- g()

  }

  res$mean <- res$mean[seq_mult(ncol(res$mean), repl)]

  modeid <- rep(seq_len(nimf + 1), repl)
  if(!is.null(ini)) modeid <- modeid + 1
  modeid <- matrix(rep(modeid, nrow(res$mean)), nrow = nrow(res$mean),
                   byrow = T)

  reps <- matrix(rep(rep(seq_len(repl), each = nimf + 1), lx),
                 nrow = lx, byrow = T)

  if(is.null(ini)){
    output <- list(xy = xyo, dt = dt, m = res$mean, repl = reps, mode = modeid)
  } else {
    output <- list(xy = xyo, ini = ini, dt = dt, m = res$mean,
                   repl = reps, mode = modeid)
  }

  if(output_sifting){
    addput <- list(even_sifting = res$even_emd_s, odd_sifting = res$odd_emd_s)
    output <- merge_list(output, addput)
  }

  output$m    <- as.matrix(output$m)
  output$repl <- as.matrix(output$repl)
  output$mode <- as.matrix(output$mode)

  if(remove == "trend" | remove == "mean"){
    if(isTRUE(bind)){
      output <- mode.bind(output, xy = trend, adjust = F,
                          name = "Trend")
    } else{
      output <- mode.in(output, xy = trend, adjust = F,
                        name = "Trend")
    }
  }

  return(output)

}
