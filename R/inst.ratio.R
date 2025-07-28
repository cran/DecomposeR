#' @title Computes instantaneous ratio of frequency
#'
#' @param pulse a pulse object (created by inst.pulse for instance)
#' @param dt depth/time. Is overridden by pulse.
#' @param f instantaneous frequency. Is overridden by pulse.
#' @param a instantaneous amplitude. Is overridden by pulse.
#' @param repl number of replicates in f
#' @param plot whether to plot an output
#' @param sqrt.rpwr,style,select,bins,cut,lines,width,height parameters to feed
#' to \code{\link{plot_ratio}} for the plots
#' @param name,ext,dir,track,openfile parameters to feed to
#' \code{\link[StratigrapheR]{pdfDisplay}} in \code{\link{plot_ratio}} for pdf plot.
#'
#' @return a list of depth/time ($dt), frequency ($f), ratio of frequency
#' ($ratio), if a is provided; the ratio power ($rpwr) i.e. the multiplication
#' of the instantaneous amplitudes of the modes two by two, the replicates id
#' ($repl)and id for the first and second frequency modes used for the ratio
#' ($l for the first, $r for the second, or $lr for the two combined)
#'
#' @examples
#' set.seed(42)
#'
#' n    <- 600
#' time <- seq_len(n)
#'
#' p1 <- 30
#' p2 <- 240
#'
#' xy <- (1 + 0.6 * sin(time *2*pi/p2)) * sin(time *2*pi/p1)  +
#'   2 * sin(time *2*pi/p2) + rnorm(n, sd = 0.5)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' dec <- extricate(xy, dt, nimf = 7, sifting = 10,
#'                 repl = 10, comb = 10, factor_noise = 10,
#'                 speak = TRUE)
#'
#' \dontrun{
#' plot_emd(dec, dir = tempdir())}
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' ht    <- inst.pulse(dec, lines = c(30, 240))
#' ratio <- inst.ratio(ht, style = "s", lines = 8)
#'
#' @importFrom StratigrapheR seq_mult
#' @export


inst.ratio <-  function(pulse = NULL, dt = NULL, f = NULL, a = NULL, repl = 1,
                        plot = TRUE, sqrt.rpwr = TRUE, style = "b", select = NA,
                        bins = 100, cut = 18, lines = NULL,
                        width = 10, height = 10, name = "Ratio", ext = ".pdf",
                        dir = tempdir(), track = TRUE, openfile = TRUE){

  if(is.null(pulse) & is.null(f) & is.null(dt)) {

    stop("Missing 'pulse', 'f' or 'dt' argument")

  } else if(!is.null(pulse)) {

    if(!is.pulse(pulse)) stop("Incorrect 'pulse' object")

  } else if (is.null(f) | is.null(dt)){

    stop("If 'pulse' is NULL, 'f' & 'dt' should be provided")

  } else {

    pulse <- as.pulse(dt = dt, f = f, a = a, repl = repl)

  }

  dt   <- pulse$dt
  f    <- pulse$f
  a    <- pulse$a
  repl <- max(pulse$repl)
  mode <- unique(pulse$mode)

  x <- ncol(f)

  if(repl > 1){
    f <- f[,seq_mult(x, x/repl)]
    f <- matrix(f, ncol = x/repl)

    if(!is.null(a)){
      a <- a[,seq_mult(x, x/repl)]
      a <- matrix(a, ncol = x/repl)
    }

    x <- ncol(f)
  }

  tri.mat <- matrix(rep(seq_len(x - 1) + 1, x-1), nrow = x - 1)

  r <- tri.mat[lower.tri(tri.mat, diag = T)]
  l <- rep(seq_len(x - 1), times = rev(seq_len(x - 1)))

  lm <- matrix(rep(mode[l], nrow(f)), nrow = nrow(f), byrow = T)
  rm <- matrix(rep(mode[r], nrow(f)), nrow = nrow(f), byrow = T)

  lr <- matrix(paste(lm, rm,sep = "/"), nrow = nrow(f))

  rpl <- rep(seq_len(repl), each = nrow(f)/repl)

  fl   <- f[,l]
  fr   <- f[,r]
  accf <- fl/fr
  colnames(accf) <- NULL

  if(!isTRUE(is.na(a))){
    al   <- a[,l]
    ar   <- a[,r]
    acca <- al*ar
    colnames(acca) <- NULL
  }

  if(is.null(a)){
    res <- list(dt = dt, f = f,ratio = accf,  repl = rpl,
                l = lm, r = rm, lr = lr)
  } else {
    res <- list(dt = dt, f = f,ratio = accf,  rpwr = acca, repl = rpl,
                l = lm, r = rm, lr = lr)
  }

  if(isTRUE(plot)){
    plot_ratio(res, sqrt.rpwr = sqrt.rpwr, style = style, select = select,
               bins = bins, cut = cut, lines = lines,
               width = width, height = height, name = name, ext = ext,
               dir = dir, track = track, openfile = openfile)
  }

  return(invisible(res))

}



