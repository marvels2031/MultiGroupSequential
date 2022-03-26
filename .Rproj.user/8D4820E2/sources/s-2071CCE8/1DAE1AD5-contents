\name{xcseqgxgs}
\alias{xcseqgxgs}
\title{Sequential graphical procedure for multiple endpoints based on group-sequential p-values}
\description{Sequential graphical procedure for multiple endpoints based on group-sequential p-values}
\usage{
xcseqgxgs(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),
          alpha=0.025,graphin=BonferroniHolm(2))
}
\arguments{
  \item{pm}{a matrix of group-sequential p-values for different endpoints (in row) at different times (in column)}
  \item{alpha}{overall familywise error rate}
  \item{graphin}{graph to be used, this is graph object defined by the gMCP package}
}
\details{Let \eqn{\tau_{min}} and \eqn{\tau_{max}} correspond to the minimum follow-up time \code{taumin} and the maximum follow-up time \code{taumax}. Let \eqn{T_f}, \eqn{C}, \eqn{E} and \eqn{R} be the follow-up time, the drop-out time, the study entry time and the total recruitment period(\eqn{R} is the last element of \code{ut}). For type 2 follow-up \eqn{T_f=min(C,\tau_{max})}. For type 3 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max})}. For type 4 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max},C)}. Let \eqn{f} be the density of \eqn{T_f}. We calculate
\deqn{\int_0^{\infty} \{at/(1+at)\}f(t)dt.}
}
\value{\item{rejected}{the index set of rejected hypotheses}
       \item{decisionsm}{rejection decision for each endpoint (row) at each timepoint (column)}
       \item{cumdecisionsm}{cumulative rejection decision for each endpoint (row) at each timepoint (column)}
}
\author{Xiaodong Luo}
\note{This provides the calculation for the variance. }
\examples{
Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-17.0.2/")
library(gMCP)
pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)
graphin=gMCP::BonferroniHolm(2)
xcseqgxgs(pm=pm,alpha=0.025,graphin=graphin)
}
\keyword{group-sequential}

