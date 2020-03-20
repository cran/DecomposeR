#' @title Destacks a pile.up() signal
#'
#' @description Destacks a signal stacked by \code{\link{pile.up}} by averaging
#' each repetition back to n multiples.
#'
#' @param x Treated signal
#' @param stack Initial stack from which the x signal is from
#' @param even Whether the x signal comes from even extension part of the
#' initial stack (if FALSE, it would come from the odd extension part)
#' @param n The multiple of destacking (has to be a multiple of n/2 (n being the
#' parameter used in \code{\link{pile.up}}), in other words a multiple of
#' length(unique(stack$id)) - 2 (minus 2 as the upper an lower extension are to
#' be removed)
#'
#' @return a matrix or a vector of the destacked signal
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
#' # Small number of repetitions ----
#'
#' opar <- par("mfrow")
#' par(mfrow = c(1,2))
#'
#' stack <- pile.up(xy, dt, 10)
#'
#' signal <- stack$even + runif(length(stack$even), -3, 3)
#'
#' res <- pile.down(signal, stack, even = TRUE, n = 5)
#'
#' plot(xy, dt, type = "l", lwd = 2, main = "Low number of repetitions")
#' lines(res, dt,  type = "l", lty = 5, col = "red")
#'
#' # High number of repetitions ----
#'
#' stack <- pile.up(xy, dt, 1000)
#'
#' signal <- stack$even + runif(length(stack$even), -3, 3)
#'
#' res <- pile.down(signal, stack, even = TRUE, n = 500)
#'
#' plot(xy, dt, type = "l", lwd = 2, main = "High number of repetitions")
#' lines(res, dt,  type = "l", lty = 5, col = "red")
#'
#' par(mfrow = c(opar))
#'
#' @export
#' @importFrom StratigrapheR seq_mult

pile.down <- function(x, stack, even, n = length(unique(stack$id)) - 2)
{
  if(class(stack) != "data.frame"){
    stop("'stack' should be a data.frame made by pile.up()")
  }

  if(paste(colnames(stack), collapse = ", ") !=
     paste(c("odt", "ndt", "id", "invert","even", "odd"), collapse = ", ")){
    stop("'stack' should be a data.frame made by pile.up()")
  }

  if(length(x) != nrow(stack)){
    stop(paste("'x' should be of same length than",
               "the number of rows of 'stack'"))
  }

  if(!isTRUE(even) & !isFALSE(even)) {
    stop("'even' should be TRUE or FALSE")
  }

  colnumb <- length(unique(stack$id))

  k.ord <- seq_mult(2 * ceiling(colnumb/2), 2)

  if(colnumb/2 != round(colnumb/2)) k.ord <- k.ord[-length(k.ord)]

  inverse_columns <- seq(2,colnumb, by = 2)
  normal_columns  <- seq(1,colnumb, by = 2)

  rev_order <- rev(seq_len(length(inverse_columns)))

  if(!even) redress <- x * stack$invert else redress <- x

  mat_red <- matrix(redress, ncol = colnumb)


  mat_red1 <- matrix(rev(mat_red[,inverse_columns]),
                     ncol = length(inverse_columns))

  mat_red2 <- mat_red1[,rev_order]

  redress1 <- cbind(mat_red[,normal_columns], mat_red2)

  redress2 <- redress1[, k.ord]

  redress3 <- redress2[,-c(1, colnumb)]

  condensed <- condense(redress3, n)

  return(condensed)

}
