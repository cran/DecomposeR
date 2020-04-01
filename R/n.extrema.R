#' @title Number of extrema/zero-crossings
#'
#' @description Computes the number of extrema and zero-crossings for different
#' groups of data, by their id or separated by NA values
#'
#' @param xy signal or decomposed signal
#' @param id the id for different groups. If any NA value is in xy, it will
#' also separate two groups of data
#' @param use.names whether to use the names in id
#' @param bound,local,zc parameters to feed to \code{\link{extremist}}
#'
#' @return a list of the number of minima ($n.min), maxima ($n.max), and, if
#' zc = TRUE, zero-crossings ($n.cross)
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
#' xy <- xy - mean(xy)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' dec <- extricate(xy, dt, nimf = 7, repl = 1, comb = 40, factor_noise = 10,
#'                 speak = TRUE)
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' n.extrema(dec$m, dec$mode)
#'
#' plot_emd(dec, select = c(6,8,9), pdf = FALSE, adapt.axis = TRUE)
#' \donttest{
#' plot_emd(dec, li = list(v = 0), adapt.axis = TRUE, dir = tempdir())}
#'
#' @export


n.extrema <- function(xy, id = NULL, use.names = TRUE,
                      bound = FALSE, local = FALSE, zc = TRUE)
{

  if(!(isTRUE(use.names) | isFALSE(use.names))) {
    stop("The 'use.names' parameter should be T or F")
  }

  xy <- as.vector(as.matrix(xy))

  if(!is.null(id)) {
    id <- as.vector(as.matrix(id))
    if(length(id) != length(xy)) stop("'xy' and 'id' should be of same length")

    if(any(duplicated(rle(id)$values))){
      warning("All identical values in 'id' should normally be consecutive ",
              "(like 1,1,2,2,3,3 not 1,1,2,2,1,1)")
    }

  }

  lid <- length(id)

  has_id <- F

  if(lid == length(xy)){
    idfactor <- factor(id, levels=unique(id))
    s        <- split(xy, idfactor)
    xy       <- unlist((sapply(s,function(x) append(x, NA))), use.names = F)
    has_id   <- T
  } else if(lid != 0) {
    stop("'id' should be NULL or of same length than 'xy' (", length(xy), ")")
  }

  ind <- seq_len(length(xy))

  idna     <- is.na(xy)
  idn       <- cumsum(idna)
  idn[idna] <- NA

  resex <- extremist(xy, bound = bound, local = local, zc = zc)

  mini <- resex$minindex$l
  maxi <- resex$maxindex$l

  minc <- sapply(split(as.integer(ind %in% mini), idn), sum)
  maxc <- sapply(split(as.integer(ind %in% maxi), idn), sum)

  if(has_id & use.names){
    names(minc) <- names(s)
    names(maxc) <- names(s)
  } else {
    names(minc) <- NULL
    names(maxc) <- NULL
  }

  if(zc){

    croi <- resex$cross$l
    croc   <- sapply(split(as.integer(ind %in% croi), idn), sum)

    if(has_id & use.names){
      names(croc) <- names(s)
    } else {
      names(croc) <- NULL
    }

    res <- list(n.min = minc, n.max = maxc, n.cross  = croc)

  } else {

    res <- list(n.min = minc, n.max = maxc)

  }

  return(res)

}
