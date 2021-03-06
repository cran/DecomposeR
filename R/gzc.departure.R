#' @title departure of instantaneous frequency to generalized zero-crossing
#'
#' @description departure of instantaneous frequency to generalized
#' zero-crossing of instantaneous freqeuncy. The departure is calculated as the
#' exponential of the absolute difference of logarithms
#' of frequencies obtained using a robust generalized zero-crossing method
#' through the \code{\link{gzc}} function (where the components are simplified
#' into extrema separated by zero-crossings) and instantaneous frequency
#' computed from another method
#'
#' @param pulse a pulse object object
#' @param dt the depth or time. Is overridden by pulse.
#' @param m a matrix of the modes to calculate the gzc frequency from. Is
#' overridden by pulse.
#' @param f a matrix of the frequencies to compare to gzc.
#' @param repl the amount of replicates in m. Is overridden by emd.
#' @param mode the mode sequence index to give to each replicated IMFs.
#' Is overridden by emd.
#' @param simplify whether to average the value for each component of each
#' replicate
#'
#' @return If simplify is TRUE, the function returns the average gzc departure
#' as a data frame where the columns stand for the modes and the rows for the
#' replicates. If simplify if FALSE, the function returns the functions returns
#' local gzc departure.
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
#' dec1 <- extricate(xy, dt, nimf = 5, repl = 1, comb = 10, sifting = 1,
#'                  factor_noise = 10, bind = TRUE, speak = TRUE)
#'
#' dec2 <- extricate(xy, dt, nimf = 6, repl = 1, comb = 100, sifting = 5,
#'                   factor_noise = 50, bind = TRUE, speak = TRUE)
#'
#' \dontrun{
#' plot_emd(dec1, name = "EMD 1", dir = tempdir())
#' plot_emd(dec2, name = "EMD 2", dir = tempdir())}
#'
#' parsimony(dec1)
#' parsimony(dec2)
#'
#' f1 <- inst.pulse(dec1, plot = FALSE)
#' f2 <- inst.pulse(dec2, plot = FALSE)
#'
#' gzc.departure(f1)
#' gzc.departure(f2)
#'
#' @export


gzc.departure <- function(pulse = NULL, dt = NULL, m = NULL, f = NULL,
                          repl = 1, mode = NULL, simplify = TRUE)
{

  if(is.null(pulse) & is.null(dt) & is.null(m) & is.null(f)) {

    stop("Missing 'emd', 'dt', 'm' or 'f' argument")

  } else if(!is.null(pulse)) {

    if(!is.pulse(pulse)) stop("Incorrect 'pulse' object")

  } else if (is.null(dt) | is.null(m) | is.null(f)){

    stop("If 'pulse' is NULL, 'dt', 'm' and 'f' should be provided")

  } else {

    pulse <- as.pulse(dt = dt, f = f, m = m, repl = repl, mode = mode)

  }

  gzc <- gzc(m = pulse$m, dt = pulse$dt, repl = pulse$repl,
             mode = unique(pulse$mode[1,]),
             output = 2, warn = F)

  gf <- gzc$f

  base <- 2

  comp <- abs(log(gf, base = base) - log(pulse$f, base = base))

  if(isTRUE(simplify)){

    res <- matrix(apply(comp, 2, mean, na.rm = T),
                  nrow = max(unique(pulse$repl[1,])), byrow = T)

    colnames(res) <- paste("C", unique(gzc$mode[1,]), sep = "")
    rownames(res) <- paste("R", unique(gzc$repl[1,]), sep = "")

  } else {

    res <- comp

  }

  res <- base^res

  return(res)

}







