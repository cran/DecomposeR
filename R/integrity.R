#' @title Integrity of a decomposition
#'
#' @description The function additions each component of a decomposition by
#' depth/time, subtract it with the original signal, and provides the absolute
#' of this subtraction. This is allows to verify if the
#' decomposition is computed correctly.
#'
#' The bulk value is the cumulated value of this proxy. If the decomposition
#' is done right the value should be very small, but non-zero due to the
#' floating-point arithmetics used by computers that generate tiny errors. Its
#' actually interesting: the first computations of the orbital solutions were
#' strongly affected by this error, as the chaotic behaviour of the equations
#' enhanced the effect of these tiny tiny errors.
#'
#' @param xy the signal
#' @param emd an emd object to test. The emd$xy original signal is not used,
#' to avoid confusion: you always have to provide the xy signal yourself.
#' @param m a matrix with columns of same length that xy, made of the
#' decomposition of the signal. Is overridden by emd.
#' @param repl the replication of decompositions in m. Is overridden by emd.
#' @param bulk whether to have a bulk value each decomposition replication, or
#' for each dt of each replication
#'
#' @return a matrix with each column being a replication, or a list of bulk
#' values for each replication
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
#' dec <- extricate(xy, dt, nimf = 7, repl = 10, comb = 10, factor_noise = 10,
#'                  sifting = 10, speak = TRUE, output_sifting = TRUE)
#'
#' integrity(xy, dec)
#'
#' @export

integrity <- function(xy, emd = NULL, m = NULL, repl = 1, bulk = TRUE)
{

  if(is.null(emd) & is.null(m)) {

    stop("Missing 'emd' or 'm' argument")

  } else if(!is.null(emd)) {

    if(!is.emd(emd)) stop("Incorrect 'emd' object")

    m    <- emd$m
    repl <- unique(emd$repl[1,])

  }

  repln <- length(repl)

  a <- condense(m, ncol(m)/repln, fun = "sum")

  b <- matrix(rep(xy, repln), ncol = repln)

  d <- abs(a - b)

  if(bulk) d <- colSums(d)/length(xy)

  return(d)

}


