#' @title Check an EMD object
#'
#' @description Provides an ensemble of check on the quality of a decomposition
#' presented as an emd object (see \code{\link{as.emd}} for more information)
#'
#' @param emd an amd object to test
#' @param xy the original signal that was decomposed: this parameter is simply
#' to insure that you are indeed comparing the decomposition to the original
#' signal, and not cheating by providing the sum of your decomposition
#' @param timelimit a time limit for the computation of the greatest common
#' rational divisor. A too long time may be indicative of a problem, typically
#' depth/time values that are not rounded adequately.
#'
#' @examples
#' set.seed(50)
#'
#' h <- rnorm(n = 1000)
#'
#' dt <- seq_len(length(h))
#'
#' alpha <- 0.95
#'
#' for(i in dt[-1]) h[i] <- alpha *  h[i-1] + h[i]
#'
#' set.seed(42)
#'
#' em <- extricate(h, dt, nimf = 7, repl = 1, comb = 100, sifting = 4,
#'                 factor_noise = 20, unit_noise = "native", speak = TRUE)
#'
#'\dontrun{
#' plot_emd(em, adapt.axis = TRUE)}
#'
#' check.emd(em, h)
#'
#' @importFrom StratigrapheR divisor seq_mult
#' @importFrom tictoc tic toc
#' @importFrom dplyr left_join
#' @importFrom stats cor
#' @import usethis
#' @export

check.emd <- function(emd, xy = NULL, timelimit = 15){

  if(!is.null(xy)){

    xy <- as.vector(xy)

    if(length(emd$dt) != length(xy)){
      stop("Incorrect xy length")
    }

    if(!(inherits(xy, "numeric") | inherits(xy, "integer"))){
      stop("Incorrect xy class")
    }

  }

  myTryCatch <- function(expr) {
    warn <- err <- NULL
    value <- withCallingHandlers(
      tryCatch(expr, error=function(e) {
        err <<- e
        NULL
      }), warning=function(w) {
        warn <<- w
        invokeRestart("muffleWarning")
      })
    list(value=value, warning=warn, error=err)
  }

  error.emd <- myTryCatch(is.emd(emd))

  error.emd.text <- "Testing if this is an 'emd' object"

  if(isTRUE(error.emd$value)) {

    ui_done(error.emd.text)

    m  <- emd$m
    id <- paste(emd$mode, emd$repl, sep = "-%uniquesep*-")

    nr <- length(unique(as.vector(emd$repl)))
    nm <- length(unique(as.vector(emd$mode)))

    ui_line(paste0("  Modes: ", "{ui_value(nm)}"))
    ui_line(paste0("  Replicates: ", "{ui_value(nr)}"))

    range.round <- signif(range(emd$dt), 3)

    ui_line("  Time/Depth (dt) range: {ui_value(range.round)}")

    ui_line("  dt GCRD (Greatest Common Rational Divisor) computing;")
    ui_line(paste("  if error, consider increasing 'timelimit'",
                  "parameter or rounding dt values;"))
    ui_line("  computing...............")

    after.hours <- function()
    {

      ui_line("  {ui_value(length(emd$dt))} initial points")

      if(is.null(xy)) {
        xyi <- emd$xy
      } else {
        xyi <- xy
      }

      integrity.value <- integrity(xyi, emd)

      intergrity.sep <- trunc(log10(max(xyi))) - trunc(log10(integrity.value))

      integrity.value.round <- signif(integrity.value, digits = 3)

      range.round <- signif(range(xyi), 3)

      ui_line("  Intensity (xy) range: {ui_value(range.round)}")

      if(all(intergrity.sep > 13)){
        ui_done(paste("Range >>",
                      "Integrity ({ui_value(integrity.value.round)})"))
      } else {
        ui_oops(paste("Integrity ({ui_value(integrity.value.round)})",
                      "is not low",
                      "enough compared to the range"))
      }

      if(is.null(xy)) {
        ui_oops(paste("Please provide xy manually to safely compute",
                      "integrity based on original signal"))
      } else {
        ui_done("xy provided independently")
      }

      error.parsimony <- myTryCatch(parsimony(emd))

      parsimony.value <- error.parsimony$value

      parsimony.value.round <- signif(parsimony.value, 3)

      if(is.null(error.parsimony$warning)){

        if(all(parsimony.value < 2.25)){

          ui_done("Parsimony ({ui_value(parsimony.value.round)}) < 2.25")

        } else if(all(parsimony.value < 3)){

          ui_todo("Parsimony ({ui_value(parsimony.value.round)}) < 3")

        } else {

          ui_oops("Parsimony ({ui_value(parsimony.value.round)}) >= 3")

        }

      } else {

        ui_oops("Warning when computing parsimony")

        warning(error.parsimony$warning)

        ui_line("Parsimony values: {ui_value(parsimony.value.round)}")

      }

      nex <- n.extrema(m, id)

      ro <- seq_mult(l = nr * nm, mult = nm)

      nexmi <- matrix(nex$n.min, ncol = nr)
      nexma <- matrix(nex$n.max, ncol = nr)

      nexsu <- nexmi + nexma

      nexsu0 <- nexsu == 0
      nexsu1 <- nexsu == 1

      all0 <- apply(nexsu0, 1,  all)
      all1 <- apply(nexsu1, 1,  all)

      any0 <- apply(nexsu0, 1,  any)
      any1 <- apply(nexsu1, 1,  any)

      w1 <- which(all0)
      w2 <- which(all1)
      w3 <- which(any0)
      w4 <- which(any1)

      l1 <- length(w1)
      l2 <- length(w2)
      l3 <- length(w3)
      l4 <- length(w4)

      trend.i   <- NULL
      is.linear <- F

      if(l1 == 1 & l3 == 1){

        ui_done(paste0("Unique trend found in mode ", "{ui_value(w3)}"))

        trend.i <- w3

        trends     <- m[ , rep(trend.i, nr) + nm * (seq(nr) - 1), drop = F]
        cors       <- apply(trends, 2, function(x) cor(x, emd$dt))
        cors.round <- round(cors, 3)

        is.linear <- all(abs(cors) > 0.999)

        if(is.linear){
          ui_line(paste0("  Trend is linear: correlation with ",
                         "dt ({ui_value(cors.round)}) > 0.999 or < -0.999"))
        } else {
          ui_line(paste0("  Trend is nonlinear: correlation with ",
                         "dt ({ui_value(cors.round)}) between 0.999 and -0.999"))
        }

      } else if(l3 == 1){

        ui_todo(paste0("Unique trend found in mode ", "{ui_value(w3)}",
                       ", but not for all replicates"))

        trend.i <- w3

      } else if(l3 > 1){

        ui_todo(paste0("Multiple trends found at modes ", "{ui_value(w3)}",
                       ", consider summing them into 1,",
                       " or reducing the amount of IMFs"))

        trend.i <- w3

      } else if(l3 == 0){

        ui_todo("No trend found, consider increasing the amount of IMFs")

      }

      if(l2 == 1 & l4 == 1){

        ui_done(paste0("Unique residue found in mode ", "{ui_value(w4)}"))

      } else if(l4 == 1){

        ui_todo(paste0("Unique residue found in mode ", "{ui_value(w4)}",
                       ", but not for all replicates"))

      } else if(l4 > 1){

        ui_todo(paste0("Multiple residues found at modes ", "{ui_value(w4)}",
                       ", consider summing them into 1,",
                       " or reducing the amount of IMFs"))

      } else if(l4 == 0){

        if(is.linear){
          ui_todo(paste0("No residue found: as the trend is linear, consider ",
                         "increasing the amount of IMFs"))
        } else {
          ui_line(paste0("  No residue found"))

        }
      }

      sym <- matrix(symmetry(m), byrow = T, nrow = nr)

      sym.round <- signif(sym, 3)

      any.sym <- apply(sym, 2, function(x) all(x < 0.6))

      if(!is.null(trend.i)){
        any.sym.t <- any.sym[-trend.i]
      } else {
        any.sym.t <- any.sym
      }

      if(any(!any.sym.t)){

        sym.cor <- which(!any.sym.t)

        for(i in sym.cor){

          ui_oops(paste("Symmetry ({ui_value(sym.round[,i])})",
                        "> 0.6 for mode {ui_value(i)}"))
        }

      } else{

        ui_done("Symmetry < 0.6")

      }

      se <- simp.emd(emd)

      if(nrow(se$multiple_extrema) != 0){

        seid <- matrix(paste(se$mode, se$repl, sep = "-%uniquesep*-"), ncol = nm)

        sp.simp <- split(se$xy, seid, lex.order = T)

        simp.ex <- function(x) length(which(x != 0 & !is.na(x)))

        real.ex <- lapply(sp.simp, simp.ex)

        d1 <- data.frame(sep = names(real.ex),
                         real.ex =  unlist(real.ex, use.names = F))

        d2 <- data.frame(sep = names(nex$n.min),
                         min =  unname(nex$n.min),
                         max = unname(nex$n.max))

        d2$n.tot <- d2$min + d2$max

        d3 <- left_join(d2, d1, by = "sep")

        d3$mode <- sub("-%uniquesep\\*-.+", "", d3$sep)
        d3$repl <- sub(".+-%uniquesep\\*-", "", d3$sep)

        d3$exc <- d3$n.tot - d3$real.ex

        umode <- unique(d3$mode)

        excm <- umode[umode %in% d3$mode[d3$exc != 0]]

        d4 <- d3[d3$mode %in% excm,]

        if(inherits(as.vector(emd$mode), "numeric") |
           inherits(as.vector(emd$mode), "integer")){
          excm <- as.numeric(excm)
        }

        for(i in excm){

          d4find <- which(d4$mode == i)

          ui_todo(paste("Excessive extrema ({ui_value(d4$exc[d4find])}) in",
                        "mode {ui_value(i)} ; total extrema:",
                        "{ui_value(d4$n.tot[d4find])}"))

        }

      } else {

        ui_done("No riding waves detected")

      }

      if(nrow(se$crossing_extrema) != 0){
        ui_oops(paste("Extrema detected at 0 intensity value; you are either",
                      "very unlucky, or deliberately doing shit.",
                      "It will probably be a problem. Deal with it."))
      } else {
        ui_done("No extrema at 0 intensity value")
      }

    }

    on.exit({

      setTimeLimit()

      after.hours()

    })

    setTimeLimit(elapsed = timelimit, transient = T)

    # Sys.sleep(14)

    tic()

    dt.div <- divisor(emd$dt)

    calc.time <- toc(quiet = T)

    setTimeLimit()

    inter.num <- as.integer(abs(((max(emd$dt) - min(emd$dt))/dt.div)) + 1.1)

    elapsed <- signif(calc.time$toc - calc.time$tic,3)

    pr.div <- signif(dt.div, 5)

    ui_done("dt GCRD: {ui_value(pr.div)} (computed in {ui_value(elapsed)} sec)")

    if(inter.num > 1e5){

      ui_todo(paste("{ui_value(inter.num)} interpolated points for regular",
                    "sampling: consider rounding the dt values"))

    } else {

      ui_done(paste("{ui_value(inter.num)} interpolated",
                    "points for regular sampling"))

    }

  } else {
    ui_oops(error.emd.text)

    if(!is.null(error.emd$warning)){
      warning(gsub(".+: ", "", error.emd$warning))
    }

    if(!is.null(error.emd$error)){
      stop(gsub(".+: ", "", error.emd$error))
    }
  }
}




