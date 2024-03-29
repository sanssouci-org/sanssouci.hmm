% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{viterbi_log}
\alias{viterbi_log}
\title{Perform the forward backward algorithm (see Rabiner 89)}
\usage{
viterbi_log(m, A_log, f0x_log, f1x_log, Pi_log)
}
\arguments{
\item{m}{the number of positions (hypothesis)}

\item{A}{a matrix 2 * 2 the transition probabilities}

\item{f0x}{a vector of the values of the density under the null hypothesis on the observations}

\item{f1x}{a vector of the values of the density under the alternative hypothesis on the observations}

\item{Pi}{a vector of the initial state probabilities}
}
\value{
alpha the forward variables, the lines corespond to the position, the first column is for state 0 and the second for state one

beta the backward variables (same as alpha)

gamma  matrix such that gamma\[i , 0\] = \eqn{P(\theta_i = 0 | X)}, gamma\[i,1\]  = \eqn{P(\theta_i = 1 | X)}

ksi  matrix such that ksi\[i , 0\] = \eqn{P(\theta_i = 0 | X)}, ksi\[i,1\]  = \eqn{P(\theta_i = 1 | X)}
}
\description{
Perform the forward backward algorithm (see Rabiner 89)
}
\examples{
 m <-  10
 A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
 f0 <- c(0, 1)
 f1 <- c(2, 1)
 Pi <- c( 0.9, 0.1)
 rdata <- simulate.data.hmm.2states(m, Pi, A, f0, f1)
 x <- rdata$x
 theta <- rdata$theta
 mod <- for_back(m, A, f0x, f1x, Pi)
}
