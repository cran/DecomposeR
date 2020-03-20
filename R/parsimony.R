#' @title Parsimony of a decomposition
#'
#' @description The function additions the absolute values of each component of
#' a decomposition by depth/time, and computes the ratio of that with the
#' absolute values of the signal. This is done either by depth/time or on the
#' time/depth-cumulated signal (i.e. the bulk signal).
#'
#' This is a proxy for parsimony: it is the factor of amplitude added by the
#' decomposition. A perfect decomposition, that does not 'invent' wiggles,
#' should approach 1, but will logically always be higher. However it is
#' influenced by the absolute value of the initial signal: if the original
#' signal is not centered around 0, the parsimony is not significative (it
#' will artificially be closer to 1). To correct for that, the residue (part of
#' the decomposition that is not centered around zero) has to be removed from
#' the original signal.
#'
#' @param emd an emd object
#' @param xy the signal
#' @param m a matrix with columns of same length that xy, made of the
#' decomposition of the signal
#' @param mode the mode sequence index to give to each replicated IMFs
#' @param repl the replication of decompositions in m
#' @param bulk whether to have a bulk value each decomposition replication, or
#' for each dt of each replication
#' @param correct the modes to remove from the original signal and decomposition
#' for a significative parsimony calculation. If NA,
#' it removes the last mode, considered as the residue. Can be a vector of
#' several integers, standing for the columns of m. If NULL, no mode is removed
#'
#' @return a matrix with each column being a replication, or a list of bulk
#' values for each replication
#' @examples
#' set.seed(42)
#'
#' n <- 500
#'
#' dt <- seq_len(n)
#' xy <- rnorm(n, mean = 0, sd = 1) + 10
#'
#' dec <- extricate(xy, dt, nimf = 7, comb = 10,
#'                  factor_noise = 1, speak = TRUE)
#'
#' \donttest{
#' plot_emd(dec, dir = tempdir())}
#'
#' parsimony(dec, correct = NULL)
#'
#' parsimony(dec)
#'
#' @importFrom stats approx lm rnorm runif sd spline
#' @export

parsimony <- function(emd = NULL, xy = NULL, m = NULL, mode = NULL, repl = 1,
                        bulk = TRUE, correct = NA)
{

  if(is.null(emd) & is.null(m) & is.null(xy)) {

    stop("Missing 'emd', 'm' or 'xy' argument")

  } else if(!is.null(emd)) {

    if(!is.emd(emd)) stop("Incorrect 'emd' object")

  } else if (is.null(xy) | is.null(m)){

    stop("If 'emd' is NULL, xy & m should be provided")

  } else {

    emd <- as.emd(xy = xy, dt = seq_len(length(xy)), imf = m,
                  mode = mode, repl = repl)

  }

  repln <- length(unique(emd$repl[1,]))

  nimf <- ncol(emd$m)/repln

  mode <- unique(emd$mode[1,])

  if(any(is.na(correct))) correct <- mode[nimf]

  b <- matrix(rep(emd$xy, repln), ncol = repln)

  if(length(correct) != 0) {
    reseq   <- mode[!(mode %in% correct)]
    residue <- mode.out(emd, lose = reseq, adjust = F)$m
    residue <- condense(residue, length(correct), fun = "sum")
    b       <- b - residue

    emd     <- mode.out(emd, lose = correct, adjust = F)
  }

  a <- condense(abs(emd$m), ncol(emd$m)/repln, fun = "sum")

  if(mean(b) > 0.2 * sd(b) | mean(b) < - 0.2 * sd(b)){
    warning("parsimony might not be significative, ",
            "try correcting by a trend or a residue so that xy is centered ",
            "around zero")
  }

  b <- abs(b)

  if(!bulk) d <- a/b else d <- colSums(a)/colSums(b)

  return(d)

}


