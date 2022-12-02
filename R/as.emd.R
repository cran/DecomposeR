#' @title Create / Check emd objects
#'
#' @description Allows to convert the result of a decomposition into a standard
#' list. The warnings of the is.emd checking function allow to identify the
#' problems.
#'
#' @param xy a vector of length n for the original signal at each dt
#' @param dt a vector of length n for the depth or time reference
#' @param imf a data.frame or matrix of n rows of the IMFs
#' @param residue a vector of length n for the residue of the decomposition
#' @param ini an optional vector of length n of the eventual initial Intrinsic
#' Mode Function xy would be a demodulation of, if it is a demodulation.
#' @param mode the mode sequence index to give to each replicated IMFs
#' @param repl the id of each replicates. The length of unique(repl) defines the
#' amount of replicates.
#' @param order the order of the imf, typically from higher frequency to lower
#' frequency
#' @param emd an emd object to check
#'
#' @return a list made of $xy (original signal), $dt (depth/time), $m (a matrix
#' of the decomposition), $repl (the replicate id of each point) and
#' $mode (the mode id of each point).
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
#' s30  <- (1 + 0.6 * sin(t*2*pi/p2)) * sin(t*2*pi/p1)
#' s240 <- 2 * sin(t*2*pi/p2)
#' sn   <- rnorm(n, sd = 0.5)
#'
#' xy <- s30 + s240 + sn
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt  <- cumsum(inter_dt)
#'
#' dec <- as.emd(xy = xy, dt = dt, imf = matrix(c(sn, s30, s240), ncol = 3))
#'
#' plot_emd(dec, pdf = FALSE)
#'
#' is.emd(dec)
#'
#' \dontrun{
#' dec$xy <- 1
#' is.emd(dec)}
#'
#' @export

as.emd <- function(xy, dt, imf, residue = NULL, ini = NULL,
                   mode = NULL, repl = 1, order = NA)
{

  ldt <- length(dt)

  xy  <- as.vector(xy)
  imf <- as.matrix(imf)

  if(ldt != length(xy)) {
    stop("The 'dt' and 'xy' parameters should be of same length")
  }

  if(ldt != nrow(imf)) {
    stop("The 'imf' parameter should be convertible into a matrix of n rows")
  }

  if(!is.null(residue)) {

    if(ldt != length(residue)){
      stop("If not NULL the 'residue' parameter should be of length n")
    }

    if(!(inherits(residue, "numeric") | inherits(residue, "integer"))){
      stop("If not NULL 'residue' should be of class numeric or integer")
    }

  }

  if(!is.null(ini)) {

    if(ldt != length(ini)){
      stop("If not NULL the 'ini' parameter should be of length n")
    }

    if(!(inherits(ini, "numeric") | inherits(ini,"integer"))){
      stop("If not NULL 'ini' should be of class numeric or integer")
    }

  }

  nimf <- ncol(imf)

  repl  <- unique(as.vector(repl))
  repln <- length(repl)

  if(!(is.na(order[[1]]) & length(order) == 1) &
     length(order) == nimf & is.numeric(order)) imf <- imf[,order]

  if(length(residue) == 0) {

    m  <- as.matrix(imf)
    nm <- nimf

  } else {

    m  <- as.matrix(imf)
    m  <- m[,seq_mult(nimf, repln, inv = T)]
    m  <- cbind(m, matrix(rep(matrix(residue), repln), ncol = repln))
    nm <- ncol(m)
    m   <- m[,seq_mult(nm, repln)]

  }

  if(!is.null(mode)){

    if(length(mode) != nm/repln) {
      stop("'mode' should have ", nm/repln, " elements")
    }

    if(!(inherits(mode, "numeric") | inherits(mode, "integer"))){
      stop("'mode' should be of class numeric or integer")
    }

    mode <- matrix(rep(mode, repln * nrow(m)),
                   ncol = ncol(m), byrow = T)

  } else {
    mode <- matrix(rep(seq_len(nm/repln), repln * nrow(m)),
                   ncol = nm, byrow = T)
  }

  if(!(inherits(repl,"numeric") | inherits(repl, "integer"))){
    stop("'repl' should be of class numeric or integer")
  }

  reps <- matrix(rep(rep(repl, each = nm/repln), ldt),
                 nrow = ldt, byrow = T)

  if(is.null(ini)){
    res <- list(xy = xy, dt = dt, m = m, repl = reps, mode = mode)
  } else {
    res <- list(xy = xy, ini = ini, dt = dt, m = m, repl = reps, mode = mode)
  }

  return(res)

}


#' @rdname as.emd
#' @export

is.emd <- function(emd)
{

  name <- deparse(substitute(emd))

  if(!all(c("xy", "dt", "m", "repl", "mode") %in% names(emd))){
    warning("The emd object should have $xy, $dt, $m, $repl and $mode elements")
    return(F)
  }

  res <- T

  tc1 <- inherits(emd$xy, "numeric") | inherits(emd$xy, "integer")
  tc2 <- inherits(emd$dt, "numeric") | inherits(emd$dt, "integer")

  if(!(tc1 & tc2)) {
    warning(name, "$xy & ", name, "$dt should be of class numeric or integer")
    res <- F
  }

  tc3 <- inherits(emd$m, "matrix")
  tc4 <- inherits(emd$repl, "matrix")
  tc5 <- inherits(emd$mode, "matrix")

  if(!(tc3 & tc4 & tc5)) {
    warning(name, "$m, ", name, "$repl & ", name,
            "$mode should be of class matrix")
    res <- F
  }

  lxy <- length(emd$xy)
  ldt <- length(emd$dt)

  nr <- length(unique(as.vector(emd$repl)))
  nm <- length(unique(as.vector(emd$mode)))

  dm    <- dim(emd$m)
  drepl <- dim(emd$repl)
  dmode <- dim(emd$mode)

  tl1 <- ldt == lxy
  tl2 <- dm[1] == lxy
  tl3 <- drepl[1] == lxy
  tl4 <- dmode[1] == lxy

  if(!tl1) {
    warning(name, "$dt should be of same length than ", name, "$xy")
    res <- F
  }

  if(!(tl2 & tl3 & tl4)) {
    warning(name, "$m,  ", name, "$repl & ", name, "$mode should have as many",
            " rows as ", name, "$xy has elements")
    res <- F
  }

  tw1 <- dm[2] == nr*nm
  tw2 <- drepl[2] == nr*nm
  tw3 <- dmode[2] == nr*nm

  if(!(tw1 & tw2 & tw3)) {
    warning(name, "$m,  ", name, "$repl & ", name, "$mode should have ", nr*nm,
            " columns, which is the amount of replicates multiplied by the",
            " amount of modes")
    res <- F
  }

  tu1 <- length(unlist(apply(emd$repl, 2, unique))) == nr*nm
  tu2 <- length(unlist(apply(emd$mode, 2, unique))) == nr*nm

  if(!(tu1 & tu2)) {
    warning("Each row in ", name, "$repl & ", name, "$mode should be identical")
    res <- F
  }

  sr <- rep(unique(emd$repl[1,]), each = nm)

  ts1 <- all(emd$repl[1,] == sr)

  if(!ts1) {
    warning("Each row in ", name, "$repl should be ", paste(sr, collapse = " "))
    res <- F
  }

  sm <- rep(emd$mode[1,seq_len(nm)], nr)

  ts2 <- all(emd$mode[1,] == sm)

  if(!ts2) {
    warning("Each row in ", name, "$repl should be a repetition of identical ",
            "mode sequence such as ", sm)
    res <- F
  }

  if(!is.null(emd$ini)){

    tc6 <- inherits(emd$ini, "numeric") | inherits(emd$ini, "integer")

    if(!(tc6)) {
      warning(name, "ini should be NULL or of class numeric or integer")
      res <- F
    }

    lini <- length(emd$ini)

    tl5 <- lini == lxy

    if(!tl5) {
      warning(name, "$ini should be NULL or of same length than ", name, "$xy")
      res <- F
    }

  }

  return(res)

}









