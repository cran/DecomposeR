#' @title Add / Remove / Bind modes in emd objects
#'
#' @param emd emd-type object
#' @param obj emd or pulse type object
#' @param xy an Instrinsic Mode Function to add
#' @param mode,keep,lose [mode.in] the position where to add the mode /
#' [mode.out] the modes to keep or lose / [mode.bind] the modes to merge
#' @param adjust whether to adapt the initial signal of an emd object ($xy in
#' the emd object) when adding or removing a mode
#' @param name the name of the new mode
#' @param reorder whether to reinitialise the index of modes when suppressing
#' one
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
#'   rnorm(n, sd = 0.5)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' dec <- extricate(xy, dt, nimf = 7, sifting = 10,
#'                  repl = 10, comb = 10, factor_noise = 10,
#'                  speak = TRUE)
#'
#' opar <- par('mfrow')
#'
#' par(mfrow = c(2,1))
#'
#' integrity(xy, dec)
#'
#' ht  <- inst.pulse(dec, plot = FALSE)
#'
#' plot_hist(x = 1/ht$f, breaks = 500, id = ht$mode,
#'           xlog = TRUE, text = TRUE, xlab = "Period",
#'           main = "Initial Decomposition")
#'
#' bound <- mode.bind(dec, mode = c(6,7))
#'
#' ht2  <- inst.pulse(bound, plot = FALSE)
#'
#' plot_hist(x = 1/ht2$f, breaks = 500, id = ht2$mode,
#'           xlog = TRUE, text = TRUE, xlab = "Period",
#'           main = "Binding of modes 6 and 7")
#'
#' par(mfrow = opar)
#'
#' \dontrun{
#' plot_emd(bound, dir = tempdir(), adapt.axis = TRUE)}
#'
#' @importFrom StratigrapheR seq_mult
#' @export

mode.in <- function(emd, xy, mode = NA, adjust = TRUE, name = "Added")
{

  if(!is.emd(emd)) stop("Incorrect 'emd' object")

  repl  <- unique(emd$repl[1,])
  repln <- length(repl)

  nimf <- ncol(emd$m)/repln
  nc   <- nrow(emd$m)

  omode <- unique(emd$mode[1,])

  if(any(is.na(mode))) mode <- max(omode) + 1

  xy <- as.vector(xy)

  lxy <- length(xy)

  if(lxy == nc){

    if(adjust) rexy <- emd$xy + xy else rexy <- emd$xy

    addxy <- rep(xy, repln)

  } else if(lxy != repln * nc){

    stop("'xy' should have the same number of elements than the original ",
         "signal [length(emd$dt)], or that length multiplied by the number of",
         " replicates [max(emd$repl) * length(emd$dt)]")

  } else {

    if(adjust) {
      stop("to add 'xy' to the original signal, it should have ",
           "the same number of elements than the original ",
           "signal [length(emd$dt)]")
    } else {
      rexy  <- emd$xy
      addxy <- xy
    }

  }

  addxy   <- matrix(addxy, ncol = repln)

  colnames(addxy) <- paste(name ,"M", mode, "R",
                           repl, sep = "")

  # change m and mode ----

  seq1 <- seq_mult(nimf * repln, nimf)
  seq2 <- seq_mult((nimf + 1) * repln, repln)

  o <- emd$mode[,seq1]

  low  <- o[1,] < mode
  high <- o[1,] >= mode

  nmode <- emd$mode[1, seq_len(nimf)]

  hmode <- nmode[nmode >= mode]
  if(any(omode == mode))  hmode <- hmode +1
  nmode <- c(nmode[nmode < mode], mode,  hmode)

  om <- emd$m[,seq1]
  nm <- cbind(om[,low], addxy,om[,high])[,seq2]

  output <- as.emd(rexy, emd$dt, imf = nm, ini = emd$ini,
                   repl = repl, mode = nmode)

  return(output)

}

#' @rdname mode.in
#' @export

mode.out <- function(obj, keep = NULL, lose = NULL, adjust = F, reorder = F)
{

  if(is.null(keep) & is.null(lose)) return(obj)

  isemd <- suppressWarnings(is.emd(obj))
  ispul <- suppressWarnings(is.pulse(obj))

  if(!isemd & !ispul) {
    stop("The 'obj' parameter should be an emd or pulse object. Use is.emd()",
         "or is.pulse on the object for more precision")
  } else if(isemd & ispul){
    stop("The 'obj' parameter should be either an emd object or a pulse",
         " object. Not both. But nice try anyway !")
  }

  repl  <- unique(obj$repl[1,])
  repln <- length(repl)

  if(isemd & adjust & repln != 1) {
    stop("To adjust, the number of replicates in the 'emd' object should be 1")
  }

  mode <- keep

  if(is.null(keep)){
    is   <- unique(obj$mode[1,])
    mode <- is[!(is %in% lose)]
  }

  out    <- obj$mode[1,] %in% mode

  keep_mode <- rep(mode, repln)
  keep_mode <- matrix(rep(keep_mode, each = nrow(obj$mode)),
                      nrow = nrow(obj$mode))

  if(isTRUE(reorder)) modeout <- NULL else modeout <- mode

  if(isemd){

    keep_m <- obj$m[,out,drop = F]

    output <- as.emd(obj$xy, obj$dt, imf = keep_m, ini = obj$ini,
                     mode = modeout, repl = repl)

    if(isemd & adjust){

      if(sum(!out) == 1) {
        remove_m <- obj$m[,!out]
      } else {
        remove_m  <- rowSums(obj$m[,!out])
      }

      rexy      <- obj$xy - remove_m
      output$xy <- rexy

    }

  } else if(ispul){

    keep_f <- obj$f[,out,drop = F]
    keep_a <- obj$a[,out,drop = F]

    output <- as.pulse(obj$dt, f = keep_f, a = keep_a,
                       mode = modeout, repl = repl)

  }

  return(output)

}

#' @rdname mode.in
#' @export

mode.bind <- function(emd, mode = NA, xy = NULL, adjust = T,
                      name = "bound")
{

  if(!is.emd(emd)) stop("Incorrect 'emd' object")

  repl  <- unique(emd$repl[,1])
  repln <- length(repl)

  if(all(is.na(mode))) mode <- max(unique(emd$mode[1,]))

  if(length(xy) == 1) xy <- rep(xy, nrow(emd$m))

  if(length(xy) != 0) {
    emd  <- mode.in(emd, xy, mode = max(mode) + 1, adjust = adjust)
    mode <- c(mode, max(mode) + 1)
  }

  bound   <- condense(emd$m[,emd$mode[1,] %in% mode],
                      length(mode), fun = "sum")
  cleaned <- mode.out(emd, lose = mode, adjust = F, reorder = T)
  output  <- mode.in(cleaned, xy = bound, mode = min(mode), adjust = F,
                     name = name)

  return(output)

}




