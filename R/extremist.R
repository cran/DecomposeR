#' @title Gives local extrema and zero crossings intervals
#'
#' @description Gives local minimas, maximas and zero crossings. Optimised for
#' large data sets; the sky is the limit (and by the sky I mean the ability of
#' R and your computer to memorise large data sets; but within this limit the
#' algorithm can handle millions of points quickly).
#'
#' @param xy the values where to find the local extremas
#' @param bound whether to consider the first and last points as both
#' minima and maxima, for special purposes. Default is F, has it should be.
#' @param local whether to consider the first and last points as local
#' minima and maxima, if TRUE by default, otherwise these first and last points
#' will be ignored
#' @param zc whether to return the zero crossings
#'
#' @return a list of the indexes of the left (l) and right (r) boundaries for
#' the minima (minindex), maxima (maxindex) and zero crossing (cross), along
#' with the number of extrema and zero crossings
#'
#' @examples
#' # Function script ----
#'
#' xy <- c(1,0,0,0,4,5,5,0.5,-0.5,0.5,0,2,2,1,-1,-1,1,1,0,0,-4,-2,2,1,0,0.5,0,
#'         NA, 0.5,0,-0.5,3,2,3,0,0.5,4,4,0)
#'
#' impressme <- 0 # Increase up to 5 or 6 to be impressed (bugs if your system
#'                # can't handle the size of the data).
#'                # If you increase it, do not run the plot script.
#'
#' xy <- rep(xy, round(10^impressme))
#'
#' print(paste("You are running ", length(xy), " points", sep = ""))
#'
#' res <- extremist(xy)
#'
#'
#'
#' # Plot script: do not run if you increase the impressme parameter ----
#'
#' mini <- unique(c(res$minindex[[1]], res$minindex[[2]]))
#' maxi <- unique(c(res$maxindex[[1]], res$maxindex[[2]]))
#' zeri <- unique(c(res$cross[[1]], res$cross[[2]]))
#'
#' l <- length(xy)
#'
#' opar <- par("mfrow")
#'
#' par(mfrow = c(3,1))
#'
#' plot(1:l, xy, type = "o",pch = 19)
#' points(mini, xy[mini], pch = 19, col = "blue")
#'
#' plot(1:l, xy, type = "o",pch = 19)
#' points(maxi, xy[maxi], pch = 19, col = "red")
#'
#' plot(1:l, xy, type = "o",pch = 19)
#' points(zeri, xy[zeri], pch = 19, col = "green")
#' abline(h = 0, col = "grey")
#'
#' par(mfrow = opar)
#'
#' @export
#' @importFrom dplyr lead lag last

extremist <- function(xy, bound = FALSE, local = TRUE, zc = TRUE)
{

  if(!(isTRUE(bound) | isFALSE(bound))) {
    stop("The 'bound' parameter should be T or F")
  }

  if(!(isTRUE(local) | isFALSE(local))) {
    stop("The 'local' parameter should be T or F")
  }

  if(!(isTRUE(zc) | isFALSE(zc))) {
    stop("The 'zc' parameter should be T or F")
  }

  if(is.na(last(xy))){
    rem <- seq_len(min(which((duplicated(cumsum(rev(is.na(xy))))))) - 1) - 1
    xy <- xy[-(length(xy) - rem)]
  }

  xy <- c(NA, xy, NA)

  l   <- length(xy)

  lxy         <- xy
  length(lxy) <- l-1
  rxy         <- xy[-1]

  li  <- 1:(l-1)
  ri  <- 2:l

  d <- data.frame(lxy = lxy, rxy = rxy, li = li, ri = ri)

  # Extrema calulation ----

  # Simplification to intervals of same slope sign

  d$diff  <- d$rxy - d$lxy
  d$slope <- sign(d$diff)
  d$off   <- d$slope - lag(d$slope)
  off     <- which(d$off == 0)

  sd <- data.frame(l = d$li[-off], r = d$ri[-off + 1], id = d$slope[-off])

  # Removal of intervals of zero slope surrounded by intervals of same slope

  sd$out <- abs(lead(sd$id) + lag(sd$id)) - 2*abs(sd$id)
  out    <- which(sd$out == 2)

  if(length(out) != 0)
  {
    outl <- unique(c(out, out+1))
    outr <- unique(c(out-1, out))
    ssd  <- data.frame(l = sd$l[-outl], r = sd$r[-outr], id = sd$id[-outr])
  } else {
    ssd <- sd
  }

  # Merging of consecutive NAs

  ssd$na  <- is.na(ssd$id)
  ssd$out <- ssd$na & lag(ssd$na)
  out2    <- which(ssd$out)

  if(length(out2) != 0)
  {
    out2l <- out2
    out2r <- out2-1
    ssdm  <- data.frame(l = ssd$l[-out2l], r = ssd$r[-out2r],
                        id = ssd$id[-out2r])
  } else {
    ssdm <- ssd
  }

  ssdm <- ssdm[ssdm$id != 0 | is.na(ssdm$id),]

  # Inversion of intervals to find peaks, separation in minima and maxima

  ext <- data.frame(l = c(1, ssdm$r), r =  c(ssdm$l, l),
                    id = c(-ssdm$id[1], ssdm$id))

  ext$id[lead(is.na(ext$id))] <- NA

  if(is.na(xy[1])) ext <- ext[-1,]
  if(is.na(last(xy))) ext <- ext[-nrow(ext),]

  # Gestion of boundaries

  if(bound){

    minloc <- unique(c(1, which(ext$id == -1 | is.na(ext$id)), nrow(ext)))
    maxloc <- unique(c(1, which(ext$id == 1 | is.na(ext$id)), nrow(ext)))

  } else {

    if(local)
    {

      low  <- which(is.na(ext$id) & is.na(lead(ext$id)))
      high <- which(is.na(ext$id) & is.na(lag(ext$id)))

      ext$id[low]  <- -ext$id[low - 1]
      ext$id[high] <- -ext$id[high + 1]

      minloc <- which(ext$id == -1)
      maxloc <- which(ext$id == 1)

    } else {

      ext <- ext[-c(1,nrow(ext)),]

      minloc <- which(ext$id == -1)
      maxloc <- which(ext$id == 1)

    }

  }

  minindex <- ext[minloc,c(1,2)]
  maxindex <- ext[maxloc,c(1,2)]

  if(zc)
  {

    # Zero crossing calculation ----

    f <- data.frame(xy = xy, i = 1:l)

    f$sign <- sign(f$xy)

    f$lead <- lead(f$sign)
    f$lag  <- lag(f$sign)

    f$iso  <- f$sign == 0 &
      (f$lead != 0 | is.na(f$lead)) &
      (f$lag != 0 | is.na(f$lag))

    f$rem <- abs(f$sign) + abs(f$lead) + abs(f$lag)

    f$sim  <- f$lead == f$lag

    # Removal of zero points surrounded by zeros

    keep1 <- f$rem != 0 | is.na(f$rem)

    sf <- data.frame(xy = f$xy[keep1], i = f$i[keep1], sign = f$sign[keep1],
                     lead = f$lead[keep1], lag = f$lag[keep1],
                     iso = f$iso[keep1], sim = f$sim[keep1])

    sf$diff1 <- sf$sign + sf$lead
    sf$diff2 <- sf$sign + sf$lag

    # Removal of points surrounded by points of identical sign, if not zero

    keep2 <- sf$diff1 == 0 | sf$diff2 == 0 | sf$sign == 0
    keep2[is.na(keep2)] <- F

    if(any(keep2)){

      ssf <- data.frame(xy = sf$xy[keep2], i = sf$i[keep2],
                        iso = sf$iso[keep2], sim = sf$sim[keep2])

      if(is.na(ssf$i[1])) ssf <- ssf[-1,]
      if(is.na(ssf$i[length(ssf$i)])) ssf <- ssf[-nrow(ssf),]

      # Rules for l and r indexes

      double <- (ssf$iso | ssf$sim)
      simple <- !double

      ssf$simple<- is.na(simple) | simple

      rs <- every_nth(which(ssf$simple), 2, empty = F)
      ls <- every_nth(which(ssf$simple), 2, empty = F, inverse = T)

      si <- 1:nrow(ssf)

      ssf$ls <- si %in% ls | !ssf$simple
      ssf$rs <- si %in% rs | !ssf$simple

      cross <- data.frame(l = ssf$i[ssf$ls], r = ssf$i[ssf$rs])

    } else {

      cross           <- data.frame(matrix(nrow = 0, ncol = 2))
      colnames(cross) <- c("l", "r")

      ncross <- 0

    }

    # Return ----

    res <- list(minindex = minindex - 1, maxindex = maxindex - 1,
                cross = cross - 1)

  } else{

    res <- list(minindex = minindex - 1, maxindex = maxindex - 1)

  }

  return(res)

}



