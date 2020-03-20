#' @title Symmetry of components
#'
#' @description The function returns the highest factor of amplitude either in
#' negative or positive values. This quantifies the symmetry of components.
#'
#' @param xy signal (vector or matrix)
#' @param names the names to use for the resulting vector. If NULL no names are
#' provided, if NA its the names of the columns of the xy matrix, if "num" it
#' the column index of the matrix xy
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
#'                  speak = TRUE, output_sifting = TRUE)
#'
#' symmetry(dec$m)
#'
#' plot_emd(dec, select = c(6,8,9), pdf = FALSE, adapt.axis = TRUE)
#'
#' @export

symmetry <- function(xy, names = "num")
{
  m <- as.matrix(xy)

  s <- sign(m)

  pos <- m * (s + 1)/2
  neg <- m * (1 - s)/2

  cp <- colSums(pos)
  cn <- -colSums(neg)

  pp <- cp/(cp+cn)
  pn <- cn/(cp+cn)

  maxp <- pp >= 0.5

  res <- pn
  res[maxp] <- pp[maxp]

  if(is.null(names)) {
    names(res) <- NULL
  } else if(names == "num"){
    names(res) <- seq_len(ncol(m))
  }

  return(res)

}

