#' @title Empirical AM and FM decomposition
#'
#' @description Applies the normalisation scheme of Huang et al., 2009 to
#' decompose any Intrinsic Mode Functions obtained (usually via Empirical Mode
#' Decomposition) into an Frequency Modulated component of amplitude 1, also
#' called carrier, and its Amplitude Modulated enveloppe. The carrier can then
#' be used to compute the instantaneous frequency via the Normalised Hilbert
#' Transform (NHT) or by calculating its Direct Quadrature (DQ) (Huang et al.,
#' 2009). HOWEVER THIS FUNCTION CAN FAIL due to overshoot or undershoot of the
#' spline fitting. Additional research is necessary.
#'
#' @param emd an emd object
#' @param m a matrix of the modes to calculate the amplitude and the frequency
#' carrier from. Is overridden by emd.
#' @param dt the depth or time. Is overridden by emd.
#' @param repl the amount of replicates in m. Is overridden by emd.
#' @param last whether to use the last mode (trend/residue).
#' @param speak whether to print a sentence at each iteration
#'
#' @references Huang, Norden E., Zhaohua Wu, Steven R. Long, Kenneth C. Arnold,
#' Xianyao Chen, and Karin Blank. 2009. ‘On Instantaneous Frequency’. Advances
#' in Adaptive Data Analysis 01 (02): 177–229.
#' https://doi.org/10.1142/S1793536909000096.
#'
#' @return a list of two matrices: $fc (frequency carrier) and $a (instantaneous
#' amplitude)
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
#'         rnorm(n, sd = 0.5)
#'
#' inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)
#'
#' dt <- cumsum(inter_dt)
#'
#' dec <- extricate(xy, dt, nimf = 7,
#'                repl = 1, comb = 100, factor_noise = 10,
#'                speak = TRUE)
#'
#' plot_emd(dec, pdf = FALSE, select = 4)
#'
#' integrity(xy, dec)
#' parsimony(dec)
#'
#' m  <- dec$m
#'
#' res <- normalise(dt = dt, m = m, last = FALSE)
#'
#' numb <- 4
#'
#' opar <- par('mfrow')
#'
#' par(mfrow = c(1,2))
#'
#' plot(m[,numb], dt, type = "l", xlab = "xy",
#'      main = paste("Mode", numb, "and AM enveloppe"))
#' lines(res$a[,numb], dt, col = "red", lty = 5, lwd = 2)
#'
#' plot(res$fc[,numb], dt, type = "l", xlab = "xy",
#'      main = "FM carrier")
#'
#' par(mfrow = opar)
#'
#' @export
#' @importFrom StratigrapheR in.lim as.lim

normalise <- function(emd = NULL, m = NULL, dt = NULL, repl = 1,
                      last = TRUE, speak = TRUE)
{

  if(is.null(emd) & is.null(m) & is.null(dt)) {

    stop("Missing 'emd', 'm' or 'dt' argument")

  } else if(!is.null(emd)) {

    if(!is.emd(emd)) stop("Incorrect 'emd' object")

  } else if (is.null(m) | is.null(dt)){

    stop("If 'emd' is NULL, m & dt should be provided")

  } else {

    fxy <- as.matrix(m)

    if(ncol(fxy) >= 2) fxy <- rowSums(m)

    emd <- as.emd(xy = fxy, dt = dt, imf = m, repl = repl)

  }

  if(isFALSE(last)){
    emd <- mode.out(emd, lose = max(emd$mode))
  }

  m  <- emd$m
  dt <- emd$dt

  sig <- as.matrix(abs(m))

  nid <- matrix(paste(emd$mode,":",emd$repl, sep = ""), ncol = ncol(emd$m))

  nextr <- n.extrema(m, nid)

  n.ext <- nextr$n.min + nextr$n.max
  n.cro <- nextr$n.cross

  if(any(n.ext == 0)){
    pos <- which(n.ext == 0)
    stop("The mode(s) ", paste(emd$mode[1,pos], collapse = ", "),
         " respectively in the replicate(s) ",
         paste(emd$repl[1,pos], collapse = ", "), " has (have) no extrema.")
  }

  if(any(n.cro == 0)){
    pos <- which(n.cro == 0)
    stop("The mode(s) ", paste(emd$mode[1,pos], collapse = ", "),
         " respectively in the replicate(s) ",
         paste(emd$repl[1,pos], collapse = ", "),
         " has (have) no zero crossing.")
  }

  nm  <- ncol(sig)
  ne  <- nrow(sig)

  depth <- matrix(rep(dt, nm), ncol = nm)

  mna           <- matrix(rep(NA, nm), ncol = nm)
  colnames(sig) <- NULL
  signa         <- as.vector(rbind(sig, mna))

  orisig    <- matrix(seq_len(nm*ne), ncol = nm)
  signature <- as.vector(rbind(orisig, mna))

  orid <- seq_len(ncol(orisig))

  # Repetitions ----

  i <- 1
  repsigna   <- signa
  normalised <- sig

  general_spline <- function(xys, dts)
  {
    spline(y = xys, x = dts, method = "natural",xout = dt)$y
  }

  repeat{

    if(speak) print(paste("Iteration", i))

    ext <- extremist(repsigna , local = F, zc = F)

    mlen <- nrow(ext$maxindex)
    mind <- (c(ext$maxindex[,1], ext$maxindex[,2])[seq_mult(2*mlen,mlen)])
    mind <- signature[unique(mind)]

    maxin <- in.lim(mind, as.lim(l = orisig[1,],
                                 r = orisig[nrow(orisig),],
                                 id = orid))

    maxxy <- normalised[maxin$x]
    maxdt <- depth[maxin$x]

    xys <- split(maxxy, maxin$id)
    dts <- split(maxdt, maxin$id)

    splined <- mapply(general_spline, xys, dts)

    # Include the boundaries that exceed the spline, interpolated maxima
    # where the minima exceed the spline or where it is negative ----

    over <- which(normalised > splined)

    add <- over[over %in% c(orisig[1,], orisig[nrow(orisig),])]

    minlen <- nrow(ext$minindex)
    minind <- (c(ext$minindex[,1], ext$minindex[,2])[seq_mult(2*minlen,minlen)])
    minpos <- signature[unique(minind)]

    excmin <- minpos[which(normalised[minpos] > splined[minpos])]

    neg_first  <- which(splined < 0)

    if(length(c(add, excmin, neg_first)) != 0){

      mind2 <- unique(sort(c(add, mind, excmin)))

      maxin2 <- in.lim(mind2, as.lim(l = orisig[1,],
                                     r = orisig[nrow(orisig),],
                                     id = seq_len(ncol(orisig))))

      maxdt2 <- as.vector(depth[maxin2$x])
      maxxy2 <- as.vector(normalised)[maxin2$x]

      d <- data.frame(id = maxin2$x,
                      xy = maxxy2, lagxy = lag(maxxy2), leadxy = lead(maxxy2),
                      dt = maxdt2, lagdt = lag(maxdt2), leaddt = lead(maxdt2),
                      ismin = maxin2$x %in% c(excmin))

      d <- d[d$ismin,]

      d$a <- (d$leadxy - d$lagxy)/(d$leaddt - d$lagdt)
      d$b <- d$lagxy - d$a * d$lagdt

      d$inter <- d$a*d$dt + d$b

      maxxy2[maxin2$x %in% excmin] <- d$inter

      dts2 <- split(maxdt2, maxin2$id)
      xys2 <- split(maxxy2, maxin2$id)

      splined <- mapply(general_spline, xys2, dts2)

      neg  <- which(splined < 0)

      if(length(neg) != 0){

        negint <- in.lim(neg, as.lim(l = maxin$x[-length(maxin$x)],
                                     r = maxin$x[-1], b = "]["))

        d2 <- data.frame(id = negint$x, lagid = negint$l, leadid = negint$r)

        d2$xy     <- normalised[d2$id]
        d2$lagxy  <- normalised[d2$lagid]
        d2$leadxy <- normalised[d2$leadid]

        d2$dt     <- depth[d2$id]
        d2$lagdt  <- depth[d2$lagid]
        d2$leaddt <- depth[d2$leadid]

        d2$a <- (d2$leadxy - d2$lagxy)/(d2$leaddt - d2$lagdt)
        d2$b <- d2$lagxy - d2$a * d2$lagdt

        d2$inter <- d2$a*d2$dt + d2$b

        mind3 <- sort(c(d2$id, maxin2$x))

        maxin3 <- in.lim(mind3, as.lim(l = orisig[1,],
                                       r = orisig[nrow(orisig),],
                                       id = seq_len(ncol(orisig))))

        maxdt3 <- as.vector(depth[maxin3$x])
        maxxy3 <- as.vector(normalised)[maxin3$x]

        maxxy3[maxin3$x %in% excmin] <- d$inter
        maxxy3[maxin3$x %in% neg]    <- d2$inter

        dts3 <- split(maxdt3, maxin3$id)
        xys3 <- split(maxxy3, maxin3$id)

        splined <- mapply(general_spline, xys3, dts3)

      }

    }

    normalised <- normalised/splined

    if(!any(normalised > 1 | normalised < 0)) break

    if(i == 15) stop("Normalisation failed (15 iterations)")

    i <- i + 1

    repsigna <- as.vector(rbind(normalised, mna))

  }

  fm <- sign(m) * normalised

  am <- sig/normalised

  colnames(am) <- colnames(fm)

  res <- list(fc = fm, a = am)

  return(res)

}




