#' @title Remove / Bind replicates in emd objects
#'
#' @param emd emd-type object
#' @param keep,lose the modes to keep or lose
#' @param reorder whether to reinitialise the index of replicates when
#' suppressing one
#' @param comb the number of replicates that have to be bound together
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
#' dec <- extricate(xy, dt, nimf = 7, repl = 20, comb = 2, factor_noise = 10,
#'                  speak = TRUE, output_sifting = TRUE)
#'
#' reduced  <- repl.out(dec, keep = c(3,4))
#'
#' parsimony(reduced)
#'
#' plot_emd(reduced, pdf = FALSE, select = c(4,6))
#'
#' combined <- repl.bind(dec, 10)
#'
#' parsimony(combined)
#'
#' plot_emd(combined, pdf = FALSE, select = c(4,6))
#'
#' @export

repl.out <- function(emd, keep = NULL, lose = NULL, reorder = FALSE)
{

  if(!is.emd(emd)) stop("Incorrect 'emd' object")

  if(is.null(keep) & is.null(lose)) return(emd)

  replin <- keep

  if(is.null(keep)){
    is     <- unique(emd$repl[1,])
    replin <- is[!(is %in% lose)]
  }

  out    <- emd$repl[1,] %in% replin
  keep_m <- emd$m[,out,drop = F]

  repln <- length(unique(emd$repl[1,]))

  keep_repl <- rep(seq_len(ncol(keep_m)/repln), repln)
  keep_repl <- matrix(rep(keep_repl, each = nrow(keep_m)), nrow = nrow(keep_m))

  if(isTRUE(reorder)) replout <- seq_len(repln) else replout <- replin

  nm <- ncol(emd$m)/repln

  output <- as.emd(emd$xy, emd$dt, imf = keep_m, ini = emd$ini,
                   mode = unique(emd$mode[1,]), repl = replout)

  return(output)

}

#' @rdname repl.out
#' @export

repl.bind <- function(emd, comb)
{

  if(!is.emd(emd)) stop("Incorrect 'emd' object")

  repln <- length(unique(emd$repl[1,]))

  ratio <- repln/comb

  if(ratio != round(ratio) | ratio < 1 | comb > repln) {
    stop(paste("The 'comb' parameter should be an integer and ",
               "divisor of the number of columns of 'm'"))
  }

  no  <- seq_mult(ncol(emd$m), repln, inv = T)

  con <- condense(emd$m[,no], comb)

  bo  <- seq_mult(ncol(con), repln/comb)

  nm  <- con[,bo]

  output <- as.emd(xy = emd$xy, dt = emd$dt, imf = nm, ini = emd$ini,
                   mode = unique(emd$mode[1,]), repl = seq_len(ratio))

}

