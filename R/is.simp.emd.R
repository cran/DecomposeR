#' @title Tests for simplified EMD
#'
#' @description Tests whether each column of a matrix is an alternation of
#' -minima zero-crossing maxima zero-crossing-
#'
#' @param xy a vector or matrix of values to test
#'
#' @examples
#' xytest1 <- c(0.5, 1,-1,-0.85,-0.5,-1,-0.5,-1,1,0.5,0,-1,0,
#'              1,-1,0,1,2,-2,1,2,1,3,0,-1,-1,3,0)
#'
#' xytest2 <- c(0, 1,-1,-0.85,-0.5,-1,-0.5,-1,1,0.5,0,0,
#'              1,1,1,1,2,-2,1,2,1,3,0,-1,-1,3,0)
#'
#' dat1 <- simp.emd(m = xytest1, dt = 1:length(xytest1))
#'
#' dat2 <- simp.emd(m = xytest2, dt = 1:length(xytest2))
#'
#' is.simp.emd(dat1$xy)
#'
#' is.simp.emd(dat2$xy)
#'
#' # There is a problem when two maxima or minima are separeted by a point at 0
#' # that does not cross any further, creating a false simplified IMF. THis is
#' # not considered as a simplified IMF by this function. However this scenario
#' # should be very rare in EMDs, but you never really know.
#'
#' @importFrom StratigrapheR mat.lag
#' @export

is.simp.emd <- function(xy)
{
  patrix <- sign(xy)
  fatrix <- abs(patrix)

  EMD       <- mat.lag(patrix, 2) + patrix
  zerocross <- mat.lag(fatrix, 1) + fatrix

  if(all(zerocross == 1 | is.na(zerocross)) & all(EMD == 0 | is.na(EMD))){
    return(TRUE)
  } else {
    return(FALSE)
  }
}





