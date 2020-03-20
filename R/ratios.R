#' @title Computes ratios of numerical values
#'
#' @description Computes ratios of numerical values
#'
#' @param x values to compute the ratio from
#'
#' @return a dataframe of $ratio, $x1 and $x2
#'
#' @examples
#' ratios(c(20,40,100,400))
#'
#' @export

ratios <- function(x)
{

  t <- sort(x)

  lx <- length(x)

  tri.mat <- matrix(rep(seq_len(lx - 1) + 1, lx-1), nrow = lx - 1)

  r <- tri.mat[lower.tri(tri.mat, diag = T)]
  l <- rep(seq_len(lx - 1), times = rev(seq_len(lx - 1)))

  ratio <- t[r]/t[l]

  xl <- t[l]
  xr <- t[r]

  res <- list(ratio = ratio, x1 = xl, x2 = xr)

  res <- as.data.frame(res)

  return(res)

}
