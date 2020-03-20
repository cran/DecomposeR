#' @title Check ratio objects
#'
#' @param ratio a ratio object to check
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
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' dec <- extricate(xy, dt, nimf = 7,
#'                  repl = 10, comb = 10, factor_noise = 10,
#'                  speak = TRUE)
#'
#' ht    <- inst.pulse(dec)
#' ratio <- inst.ratio(ht, plot = FALSE)
#'
#' is.ratio(ratio)
#'
#' @export

is.ratio <- function(ratio)
{

  name <- deparse(substitute(ratio))

  ra <- ratio

  if(!all(c("dt", "f", "ratio", "repl", "l", "r", "lr") %in% names(ra))){
    warning("The ratio object should have $dt, $f, $ratio, $repl,",
            " $l, $r and $lr elements")
    return(F)
  }

  res <- T

  tc1 <- inherits(ra$dt, "numeric") | inherits(ra$dt, "integer")
  tc2 <- inherits(ra$repl, "numeric") | inherits(ra$repl, "integer")

  if(!(tc1 & tc2)) {
    warning(name, "$dt and ", name, "$repl should be of class numeric or integer")
    res <- F
  }

  tc3 <- inherits(ra$f, "matrix")
  tc4 <- inherits(ra$ratio, "matrix")
  tc5 <- inherits(ra$l, "matrix")
  tc6 <- inherits(ra$r, "matrix")
  tc7 <- inherits(ra$lr, "matrix")

  if(!(tc3 & tc4 & tc5 & tc6 & tc7)) {
    warning(name, "$f, ", name, "$ratio, ", name,
            "$repl, ", name, "$l, ", name, "$r and ", name,
            "$lr should be of class matrix")
    res <- F
  }

  ldt <- length(ra$dt)
  nr  <- length(unique(as.vector(ra$repl)))
  df  <- dim(ra$f)

  lrat <- sum(seq_len(df[2])[-df[2]])

  drepl  <- length(ra$repl)
  dratio <- dim(ra$ratio)
  drpwr  <- dim(ra$rpwr)
  dl     <- dim(ra$l)
  dr     <- dim(ra$r)
  dlr    <- dim(ra$lr)

  tl1 <- drepl == ldt * nr
  tl2 <- dratio[1] == ldt * nr
  tl3 <- drpwr[1] == ldt * nr
  tl4 <- dl[1] == ldt * nr
  tl5 <- dr[1] == ldt * nr
  tl6 <- dlr[1] == ldt * nr

  if(!tl1) {
    warning(name, "$repl should be of same length than ", name,
            "$dt multiplied by the number of replicates")
    res <- F
  }

  if(!(tl2 & tl3 & tl4 & tl5 & tl6)) {
    warning(name, "$ratio, ", name, "$rpwr, ", name, "$l, ", name, "$r, & ",
            name ,"$lr should have as many rows as the amount of elements in ",
            name, "$dt multiplied by the number of replicates")
    res <- F
  }

  tw1 <- dratio[2] == lrat
  tw2 <- drpwr[2] == lrat
  tw3 <- dl[2] == lrat
  tw4 <- dr[2] == lrat
  tw5 <- dlr[2] == lrat

  if(!(tw1 & tw2 & tw3 & tw4 & tw5)) {
    warning(name, "$ratio, ", name, "$rpwr, ", name, "$l, ", name, "$r, & ",
            name ,"$lr should should have ", lrat,
            " columns, for each combination of the columns of ", name, "$f")
    res <- F
  }

  tu1 <- length(unlist(apply(ratio$l, 2, unique))) == lrat
  tu2 <- length(unlist(apply(ratio$r, 2, unique))) == lrat
  tu3 <- length(unlist(apply(ratio$lr, 2, unique))) == lrat

  if(!(tu1 & tu2 & tu3)) {
    warning("Each row in ", name, "$l, ", name, "$r & ",name,
            "$lr should be identical")
    res <- F
  }

  replm <- matrix(ra$repl, ncol = nr)

  tu4 <- length(unlist(apply(replm, 2, unique))) == nr

  if(!tu4) {
    warning("If repl is converted into a matrix having the same amount of ",
            "columns that ", name, "$repl has unique indexes, these ",
            "indexes should be separated into one column for each")
    res <- F
  }

  tlr1 <- length(unique(ratio$l[1,])) == df[2] - 1
  tlr2 <- length(unique(ratio$r[1,])) == df[2] - 1

  if(!(tlr1 & tlr2)) {
    warning(name, "$l & ", name, "$r should be made of ",df[2] - 1,
            " unique elements")
    res <- F
  }

  tlr3 <- length(unique(ratio$lr[1,])) == lrat

  if(!(tlr3)) {
    warning(name, "$lr should be made of ", lrat," unique elements")
    res <- F
  }

  tlri <- all(paste(ratio$l[1,], ratio$r[1,], sep = "/") == ratio$lr[1,])

  if(!(tlri)) {
    warning(name, "$lr should be made of ", name, "$l and ", name,
            "$r pasted and separated by '/'")
    res <- F
  }

  return(res)

}
