#' @title Condenses columns of matrix
#'
#' @description Condenses columns of a matrix by averaging or summing
#' them. The condensing can be done partially: a multiple of the repetitions can
#' be averaged or summed to keep some repetitions.
#'
#' @param m matrix of repeated signal, each column being a repetition
#' @param n the number of repetitions that will be averaged/summed
#' @param fun the function to apply to each repetition: "mean" or "sum".
#'
#' @return a matrix with n times less columns
#' @examples
#' m <- matrix(rep(seq(100, 800, 100), each = 10) + rep(1:10, 8), ncol = 8)
#'
#' m
#'
#' condense(m, 4)
#'
#' @export

condense <- function(m, n, fun = "mean")
{

  m <- as.matrix(m)

  ratio <- ncol(m)/n

  if(ratio != round(ratio) | ratio < 1) {
    stop(paste("The 'n' parameter should be a",
               "integer and divisor of the number of columns of 'm'"))
  }

  m1 <- matrix(as.vector(t(m)), ncol = n, byrow = T)

  if(fun == "mean") {
    m2 <- rowMeans(m1)
  } else if (fun == "sum") {
    m2 <- rowSums(m1)
  } else {
    stop("The 'fun' parameter should be 'mean' or' sum'")
  }

  res <- matrix(as.vector(t(m2)), ncol = ncol(m)/n, byrow = T)

  return(res)

}
