\name{xcmaurerbretz}
\alias{xcmaurerbretz}
\title{Maurer-Bretz sequential graphical approach}
\description{This will conduct group-sequential testing for multiple endpoints based on Maurer-Bretz approach}
\usage{
xcmaurerbretz(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4)),
                         informationm=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4),
                         spending=rep("OBF",4),param.spending=rep(1,4),
                         alpha=0.025,direction=-1,graphin=BonferroniHolm(nrow(xm)),
                         alpha.low=1e-10,retrospective=0)
}
\arguments{
  \item{xm}{a matrix of test statistics for each endpoint (in row) and each interim (in column)}
  \item{informationm}{a matrix of information fractions for the statistics \code{xm}}
  \item{spending}{spending functions for each endpoint}
  \item{param.spending}{parameters in the spending functions}
  \item{alpha}{overall familywise error rate}
  \item{direction}{-1: (one-sided)reject if test stat is smaller than the critical value;
                    1: (one-sided)reject if test stat is smaller than the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
  \item{graphin}{a graph object generated from gMCP}
  \item{alpha.low}{default is 1e-10, if allocated alpha is smaller than this number, the corrsponding critical value will be set to infinity}
  \item{retrospective}{retrospective: 0 (default) only compares the current test statistic with the updated critical value, 1 compares all the test statistics up to the current one with the updated critical values. Even though retrospectively looking at the values are statsitical valid in terms of control of the type-1 error rate, not retrospectively looking at the past comparisons avoids the dilemma of retrospectively inflating the alpha level}
}
\details{Let \eqn{\tau_{min}} and \eqn{\tau_{max}} correspond to the minimum follow-up time \code{taumin} and the maximum follow-up time \code{taumax}. Let \eqn{T_f}, \eqn{C}, \eqn{E} and \eqn{R} be the follow-up time, the drop-out time, the study entry time and the total recruitment period(\eqn{R} is the last element of \code{ut}). For type 2 follow-up \eqn{T_f=min(C,\tau_{max})}. For type 3 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max})}. For type 4 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max},C)}. Let \eqn{f} be the density of \eqn{T_f}. We calculate
\deqn{\int_0^{\infty} \{at/(1+at)\}f(t)dt.}
}
\value{\item{Hrej}{rejcted hypotheses}
       \item{rejected}{rejcted hypotheses}
       \item{decisionsm}{rejcted hypotheses}
       \item{cumdecisionsm}{rejcted hypotheses at each interim}
}
\author{Xiaodong Luo}
\note{This provides the calculation for the variance. }
\examples{
%library(rJava)
Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-17.0.2/")
library(gMCP)
xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4))
im=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4)
spending=rep("OBF",4)
param.spending=rep(1,4)
graphin=gMCP::BonferroniHolm(nrow(xm))
xcmaurerbretz(xm=xm,
                informationm=im,
                spending=spending,
                param.spending=param.spending,
                graphin=graphin)
}
\keyword{graphical procedure}
\keyword{group-sequential}
\keyword{Maurer-Bretz}

