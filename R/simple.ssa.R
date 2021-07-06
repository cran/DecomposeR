#' @title Simple SSA decomposition
#'
#' @description Simple wrapper for Singular Spectrum Analysis, using the
#' functions of the Rssa package (which is not installed by default by the
#' DecomposeR package, you should install it independently). This function
#' allows unevenly sampled data.
#'
#' @param xy signal to be decomposed
#' @param dt depth/time
#' @param n maximum amount of components
#' @param remove whether to remove a linear trend ("trend", is the default),
#' a mean value ("mean"), or to decompose as is (any other value)
#' @param groups which components to regroup (list of the indices of elementary
#' components to be regrouped, the entries of the list can be named, see
#' the reconstruct() function in the Rssa package for more information)
#' @param plot whether to show a visualisation of the importance of each
#' component
#' @param ... any arguments to by given to the ssa() function (see Rssa package
#' for more information)
#'
#' @return a list made of $xy (original signal), $dt (depth/time), $m (a matrix
#' of the decomposition), $repl (the replicate id of each point) and $mode (the
#' mode id of each point).
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
#'   rnorm(n, sd = 0.5) + 0.01 * t
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' res <- simple.ssa(xy, dt, groups = list(c(1,2), c= 3:10))
#'
#' parsimony(res)
#'
#' integrity(xy, res)
#'
#' \dontrun{
#' plot_emd(res, style = 1)}
#'
#' @importFrom stats residuals
#' @export

simple.ssa <- function(xy, dt, n = 10,
                       remove = "trend", groups = list(), plot = T,
                       ...)
{

   error <- requireNamespace("Rssa", quietly = T)

   if(!error) {
     stop("For this function to run you should install the Rssa package")
  }

  if(remove == "mean") {
    trend <- mean(xy)
    nxy   <- xy - trend
  } else if(remove == "trend"){
    lm.0  <- lm(xy ~ dt)
    trend <- lm.0$coeff[2] * dt + lm.0$coeff[1]
    nxy    <- xy - trend
  } else {
    nxy <- xy
  }

  rsp <- respace(dt, nxy)

  # ssa.res <- ssa(rsp$x, neig = n)
  ssa.res <- Rssa::ssa(rsp$x, neig = n, ...)

  if(isTRUE(plot)){
    plot(ssa.res$sigma, ylab = "Norms",
         main = "Component Norms", type = "o", pch = 19)
  }

  li <- as.list(seq_len(n))

  names(li) <- paste("C_complex_name", seq_len(n), sep = "")

  if(length(groups) != 0){

    grpd  <- unlist(groups)
    ugrpd <- unique(grpd)

    if(length(grpd) != length(ugrpd)){
      stop("The values in the 'groups' parameter should only appear once")
    }

    lin <- c(groups, li[-grpd])

  } else {

    lin <- li

  }

  rec <- Rssa::reconstruct(ssa.res, groups = lin)

  comp <- data.frame(rec[seq_len(length(lin))])
  resi <- residuals(rec)

  comp.rsp <- comp[rsp$initial,,drop = F]
  resi.rsp <- data.frame(NC = resi[rsp$initial])

  comps <- cbind(comp.rsp, resi.rsp)

  if(remove == "mean" | remove == "trend"){
    res <- as.emd(xy, dt, imf = comps, residue = trend)
  } else {
    res <- as.emd(xy, dt, imf = comps)
  }

  return(res)

}

