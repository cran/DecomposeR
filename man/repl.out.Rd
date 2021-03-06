% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repl.out.R
\name{repl.out}
\alias{repl.out}
\alias{repl.bind}
\title{Remove / Bind replicates in emd objects}
\usage{
repl.out(emd, keep = NULL, lose = NULL, reorder = FALSE)

repl.bind(emd, comb)
}
\arguments{
\item{emd}{emd-type object}

\item{keep, lose}{the modes to keep or lose}

\item{reorder}{whether to reinitialise the index of replicates when
suppressing one}

\item{comb}{the number of replicates that have to be bound together}
}
\description{
Remove / Bind replicates in emd objects
}
\examples{
set.seed(42)

n <- 600
t <- seq_len(n)

p1 <- 30
p2 <- 240

xy <- (1 + 0.6 * sin(t*2*pi/p2)) * sin(t*2*pi/p1)  + 2 * sin(t*2*pi/p2) +
  rnorm(n, sd = 0.5) + t * 0.01

inter_dt <- round(runif(length(xy), min = 0.5, max = 1.5),1)

dt <- cumsum(inter_dt)

dec <- extricate(xy, dt, nimf = 7, sifting = 10,
                 repl = 20, comb = 2, factor_noise = 10,
                 speak = TRUE, output_sifting = TRUE)

reduced  <- repl.out(dec, keep = c(3,4))

parsimony(reduced)

plot_emd(reduced, pdf = FALSE, select = c(4,6))

combined <- repl.bind(dec, 10)

parsimony(combined)

plot_emd(combined, pdf = FALSE, select = c(4,6))

}
