#' @title Repeat and stack a signal in central and line symmetry
#'
#' @description Repeats and stacks a signal duplicated in central (even) and
#' line (odd) symmetry to apply Ensemble Empirical Mode Decomposition (EEMD) on
#' one single vector following the simple boundary rule of Zeng and He (2004).
#' This allows to avoid the iterations that are typical of EEMD. A complete
#' set of signal is added by default at the upper and lower part of the stack,
#' to be removed in the end process.
#'
#' @param xy the signal
#' @param dt the depth/time positions of each xy
#' @param n the number of replicates you want. It has to be a multiple of two,
#' as you will generate two stacks: the even and the odd one.
#' @param warn whether you want to be annoyed
#'
#' @return a dataframe of the original dt (odt), the stack-modified dt (ndt),
#' the inversion factor to change the even stack into the odd one and
#' vice-versa (invert), the even xy stack (even) and the odd one (odd)
#'
#' @examples
#' set.seed(42)
#'
#' n <- 200
#' t <- seq_len(n)
#'
#' p1 <- 25
#' p2 <- 75
#'
#' xy <- (1 + 0.6 * sin(t*2*pi/p2)) * sin(t*2*pi/p1)  + 2 * sin(t*2*pi/p2) +
#'   rnorm(n, sd = 0.5)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#' inter_dt[20] <- 20
#'
#' dt <- cumsum(inter_dt)
#'
#' opar <- par()$mfrow
#' par(mfrow = c(1,1))
#'
#' res <- pile.up(xy, dt, 4)
#'
#' par(mfrow = c(2,1))
#' plot(res$ndt, res$even, type = "l", col = "blue")
#' plot(res$ndt, res$odd,  type = "l", col = "red")
#'
#' par(mfrow = c(opar))
#'
#' @export
#' @importFrom dplyr lag

pile.up <- function(xy, dt, n, warn = TRUE)
{
  # Conditions ----

  n  <- as.integer(n)

  if(n/2 != round(n/2)) stop("The n parameter should be a multiple of 2")

  if(is.unsorted(dt) & is.unsorted(-dt)){
    stop("The 'dt' values should be sorted")
  }

  # Determine smallest interval ----

  intervals <- dt - lag(dt)

  sep <- min(abs(intervals), na.rm = T)

  # Determine the number of repetitions ----
  # By default each iteration will deal with at max 10 million points, any bigger
  # and the data should be divided in chunks of less than 10 million points

  lx <- length(xy)

  l_number <- 2 + n/2

  iteration_size <- l_number * lx

  if(warn == T){
    if(iteration_size > 10^7) {
      warning(paste("The signal you want to decompose is made of", lx,
                    "data points. It will be repeated", l_number,
                    "times (2 + n/2), making an object of",
                    iteration_size, "data points. This is a lot, you may",
                    "want to iterate this function with less repetitions and",
                    "different seeds."))
    }
  }

  # Is the number of total repetitions odd or even ? ----

  if(l_number/2 != round(l_number/2)) odd_n <- T else odd_n <- F

  ln2 <- floor(l_number/2)

  # Y reversals

  y_nor <- c(sep, intervals[-1])
  y_rev <- c(sep, rev(intervals[-1]))

  eoy <- rep(c(y_nor, y_rev), ln2)

  # Y position ----

  ory <- rep(c(dt, rev(dt)), ln2)

  # X reversals and inversions ----

  rev_x <- rev(xy)

  ex <- rep(c(xy, rev_x), ln2)

  invert <- rep(c(1,-1), ln2, each = lx)

  # If odd number of repetitions ----

  if(odd_n){
    ex     <- c(ex, xy)
    invert <- c(invert, rep(1, lx))
    eoy    <- c(eoy, y_nor)
    ory    <- c(ory, dt)
  }

  ox  <- ex * invert
  eoy <- cumsum(eoy)

  # Markers ----

  em <- seq_len(2 * ln2) - 1

  if(odd_n) em <- c(em, "out_2") else em[length(em)] <- "out_2"

  em[1] <- "out_1"

  em <- rep(em, each = lx)

  res <- data.frame(odt = ory, ndt = eoy, id = em, invert = invert,
                    even = ex, odd = ox, stringsAsFactors = F)

  return(res)

}


