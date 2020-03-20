#' @title Simplifies the components of an EMD
#'
#' @description Simplifies the component of an EMD to only extremas and
#' zero-crossings, and outputs problematic extrema: multiple extrema (extrema
#' not separated by zero-crossings) and crossing extrema (extrema at zero).
#'
#' @param emd emd-type object
#' @param m a matrix of the amplitude values (xy) of the components, each
#' column being a component. Each column should have the same number of non NA
#' values. Vectors, for 1 component, are accepted. Is overridden by emd.
#' @param dt the depth or time value. Is overridden by emd.
#' @param repl the amount of replicates in m. Is overridden by emd.
#' @param use.names whether to use the column names to identify problematic
#' extrema
#'
#' @return a list of the depth or time values ($dt) of the simplified IMF
#' (Intrinsic Mode Function), of their amplitude ($xy), and of the position
#' and component of problematic multiple extrema ($multiple_extrema) and
#' crossing extrema ($crossing_extrema)
#'
#' @examples
#' xytest <- c(0.5, 1,-1,-0.85,-0.5,-1,-0.5,-1,1,0.5,0,0,
#'             1,-1,0,1,2,-2,1,2,1,3,0,-1,-1,3,0)
#'
#' repeatafterme <- 2
#'
#' m  <- matrix(rep(xytest,repeatafterme), ncol = repeatafterme)
#' dt <- 1:length(xytest)
#'
#' res <- simp.emd(m = m, dt = dt, repl = repeatafterme)
#'
#' opar <- par("mfrow")
#'
#' par(mfrow = c(1,1))
#'
#' plot(dt, xytest, type = "o", pch = 19)
#' abline(h = 0, col = "grey")
#'
#' me <- res$multiple_extrema$dt[res$multiple_extrema$repl == 1]
#' ce <- res$crossing_extrema$dt[res$multiple_extrema$repl == 1]
#'
#' abline(v = me, col = "orange")
#' abline(v = ce, col = "darkred")
#'
#' points(res$dt[,1], res$xy[,1], col = "red", pch = 19)
#'
#' par(mfrow = opar)
#'
#' @export
#' @importFrom dplyr arrange left_join desc
#' @importFrom StratigrapheR every_nth enlarge

simp.emd <- function(emd = NULL, m = NULL, dt = NULL, repl = 1, use.names = FALSE)
{
  # Preping data ----

  if(is.null(emd) & is.null(m) & is.null(dt)) {

    stop("Missing 'emd', 'm' or 'dt' argument")

  } else if(!is.null(emd)) {

    if(!is.emd(emd)) stop("Incorrect 'emd' object")

  } else if (is.null(m) | is.null(dt)){

    stop("If 'emd' is NULL, m & dt should be provided")

  } else {

    emd <- as.emd(xy = rowSums(as.matrix(m)), dt = dt, imf = m, repl = repl)

  }

  m  <- emd$m
  dt <- emd$dt

  m1 <- as.matrix(m)

  coln <- seq_len(ncol(m1))

  n <- ncol(m1)
  l <- nrow(m1)

  if(is.vector(dt)){

    if(length(dt) != l) {
      stop(paste("If the 'dt' parameter is a vector it should have the ",
                 "length of the number of rows of 'm'", sep = ""))
    }

    dt1 <- matrix(rep(dt, n), ncol = n)

  } else if(is.matrix(dt)){

    if(ncol(dt) != n | nrow(dt) != l) {
      stop(paste("If the 'dt' parameter is a matrix it should have the ",
                 "same number of rows and columns than 'm'", sep = ""))
    }

    dt1 <- dt

  }

  xyf <- as.vector(rbind(m1,rep(NA, ncol = ncol(m1))))
  dtf <- as.vector(rbind(dt1,rep(NA, ncol = ncol(dt1))))

  # Identification of each point (in extension of data) ----

  sepid1 <- matrix(rep(seq_len(n), l), ncol = n, byrow = T)
  sepid  <- as.vector(rbind(sepid1,rep(0, ncol = ncol(sepid1))))

  # Get extremas and zero-crossings ----

  rex <- extremist(xyf, local = F)

  resmin <- rex$minindex
  resmax <- rex$maxindex
  rescro <- rex$cross

  # Averaging extremas on two or more points ----

  minxy <- xyf[resmin$l]
  mindt <- (dtf[resmin$l] + dtf[resmin$r])/2

  maxxy <- xyf[resmax$l]
  maxdt <- (dtf[resmax$l] + dtf[resmax$r])/2

  # Averaging zero-crossings of several points (e.g. 2 succesive zeros) ----

  lxy <- xyf[rescro$l]
  rxy <- xyf[rescro$r]

  ldt <- dtf[rescro$l]
  rdt <- dtf[rescro$r]

  slope   <- (ldt - rdt)/(lxy - rxy)
  crossdt <- ldt - slope * lxy

  replace <- is.nan(crossdt)

  crossraw <- (ldt + rdt)/2

  crossdt[replace] <- crossraw[replace]

  crossxy <- rep(0,length(crossdt))

  # Set d to Identify problems owith the -min cross max cross- pattern ----

  d <- data.frame(dt = c(mindt,maxdt,crossdt), xy = c(minxy, maxxy, crossxy))

  d$id <- c(rep("min", length(mindt)),
            rep("max", length(maxdt)),
            rep("cross", length(crossdt)))

  d$coln <- c(sepid[resmin$l], sepid[resmax$l], sepid[rescro$l])

  d <- arrange(d, coln, dt)

  alt <- data.frame(id = c("min", "max", "cross"), cond = c(1,1,-1),
                    stringsAsFactors =  F)

  d <- left_join(d,alt, by = "id")

  # Identify crossing extrema (extrema at zero) ---

  keepce         <- which(d$dt == lead(d$dt))

  if(length(keepce) != 0){
    d$cond[keepce] <- -1
    d$id[keepce]   <- "cross ext"

    d <- d[-(keepce + 1),]
  }

  # Identify multiple extrema (not separated by zero crossing) ----

  d$ct  <- d$cond + lead(d$cond)

  d$trans <- d$coln != lead(d$coln)

  d$ct[d$trans]   <- 0

  d$int <- cumsum(as.numeric(d$ct == 2 & lag(d$ct) == 0))

  me <- which(d$ct == 2)
  me <- sort(unique(c(me, me + 1)))

  d$me <- seq_len(nrow(d)) %in% me

  dout   <- d

  d      <- d[,c(1:3,9,4:8)]

  # Finding multiple extrema to remove ----

  if(length(me) != 0){

    w <- data.frame(i = me, dt = d$dt[me], xy = d$xy[me],
                    coln = d$coln[me], int = d$int[me])

    w$id    <- paste(w$coln, "C", w$int, sep = "")
    w$glo_var_abso  <- abs(w$xy)
    w$glo_var_sim   <- paste(w$id, "V", w$glo_var_abso, sep = "")

    wss <- w[xor(duplicated(w$glo_var_sim),
                 duplicated(w$glo_var_sim, fromLast=TRUE)), ]

    wss$glo_var_i <- wss$i

    wss <- arrange(wss, glo_var_sim, glo_var_i)

    even <- every_nth(seq_len(nrow(wss)), 2, empty = F)
    odd  <- every_nth(seq_len(nrow(wss)), 2, empty = F, inverse = T)

    swsw      <- data.frame(i = wss$i[odd],
                            ldt = wss$dt[odd],
                            rdt = wss$dt[even])

    swsw$mean <- (swsw$ldt + swsw$rdt)/2

    ws <- arrange(w, desc(glo_var_abso))

    removal <- duplicated(ws$id)

    # Multiple extrema of same xy value ----

    iso <- ws[which(!removal),]

    both <- match(iso$i, swsw$i)

    swsw <- swsw[both[!is.na(both)], ]

    dout$dt[swsw$i] <- swsw$mean

    removed <- d[sort(ws$i[removal]),]

    rr <- removed$id == "max" & sign(removed$xy) < 0 |
      removed$id == "min" & sign(removed$xy) > 0 | removed$xy == 0

    removed <- removed[!rr,]

    if(length(ws$i[removal]) != 0) {
      dout <- dout[-ws$i[removal],]
    }

    # Signal multiple extrema ----

    multiple_extrema <- d[d$me & d$xy != 0,c(1,5)]

    colnames(multiple_extrema) <- c("dt", "name")

  } else {

    multiple_extrema <- data.frame(matrix(ncol = 2, nrow = 0))

    colnames(multiple_extrema) <- c("dt", "name")

  }

  if(length(keepce) != 0){

    crossing_extrema <- d[d$id == "cross ext",c(1,5)]

    colnames(crossing_extrema) <- c("dt", "name")

    if(nrow(crossing_extrema) != 0){
      warning(paste("Extrema at zero detected, see $crossing_extrema",
                    " in the results for position", sep = ""))
    }

  } else {

    crossing_extrema <- data.frame(matrix(ncol = 2, nrow = 0))

    colnames(crossing_extrema) <- c("dt", "name")

  }

  # Get signal to distinguished form ----

  df <- dout[,-c(3,5,6,7,8)]

  nr <- max(tabulate(df$coln))

  spdt <- split(df$dt, df$coln)
  spxy <- split(df$xy, df$coln)

  add.NA <- function(v, n) c(v, rep(NA, n-length(v)))

  listdt <- lapply(spdt, add.NA, n = nr)
  resdt  <- do.call(cbind, listdt)

  listxy <- lapply(spxy, add.NA, n = nr)
  resxy  <- do.call(cbind, listxy)

  remain <- unique(d$coln)

  resdtf <- data.frame(matrix(nrow = nrow(resdt), ncol = length(coln)))
  resdtf[,remain] <- resdt

  resxyf <- data.frame(matrix(nrow = nrow(resxy), ncol = length(coln)))
  resxyf[,remain] <- resxy

  cn <- colnames(m)

  if(isTRUE(use.names) & !is.null(cn)){

    colnames(resdtf) <- cn
    colnames(resxyf) <- cn

    multiple_extrema$name <- cn[multiple_extrema$name]
    crossing_extrema$name <- cn[crossing_extrema$name]

  } else {

    dref <- data.frame(name = coln, mode = emd$mode[1,], repl = emd$repl[1,])

    multiple_extrema <- left_join(multiple_extrema, dref, by = "name")[,-2]
    crossing_extrema <- left_join(crossing_extrema, dref, by = "name")[,-2]

  }

  nmode <- matrix(rep(emd$mode[1,], nrow(resdtf)),
                  nrow = nrow(resdtf), byrow = T)

  nrepl <- matrix(rep(emd$repl[1,], nrow(resdtf)),
                  nrow = nrow(resdtf), byrow = T)

  names(resdtf) <- colnames(m)
  names(resxyf) <- colnames(m)

  resdtf <- as.matrix(resdtf)
  resxyf <- as.matrix(resxyf)

  res <- list(dt = resdtf, xy = resxyf, mode = nmode, repl = nrepl,
              multiple_extrema = multiple_extrema,
              crossing_extrema = crossing_extrema)

  return(res)

}
