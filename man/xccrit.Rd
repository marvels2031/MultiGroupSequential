\name{xccrit}
\alias{xccrit}
\title{Calculate critical values}
\description{This untility function calculates the critical values}
\usage{
xccrit(direction=-1,alpha=0.025,informationRates=c(0.4,0.7,1),
                 userAlphaSpending=c(0.01,0.015,0.025),alpha.low=1e-10)
  
}
\arguments{
  \item{direction}{-1: (one-sided)reject if test stat is smaller than the critical value;
                    1: (one-sided)reject if test stat is smaller than the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
  \item{alpha}{overall familywise error rate}
  \item{informationRates}{information fractions}
  \item{userAlphaSpending}{alpha spent at each interim}
  \item{alpha.low}{default is 1e-10, if allocated alpha is smaller than this number, the corrsponding critical value will be set to infinity}
}
\details{Let \eqn{\tau_{min}} and \eqn{\tau_{max}} correspond to the minimum follow-up time \code{taumin} and the maximum follow-up time \code{taumax}. Let \eqn{T_f}, \eqn{C}, \eqn{E} and \eqn{R} be the follow-up time, the drop-out time, the study entry time and the total recruitment period(\eqn{R} is the last element of \code{ut}). For type 2 follow-up \eqn{T_f=min(C,\tau_{max})}. For type 3 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max})}. For type 4 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max},C)}. Let \eqn{f} be the density of \eqn{T_f}. We calculate
\deqn{\int_0^{\infty} \{at/(1+at)\}f(t)dt.}
}
\value{\item{crit}{critical values}
}
\author{Xiaodong Luo}
\note{This provides the calculation for the variance. }
\examples{
xccrit(direction=-1,alpha=0.025,informationRates=c(0.4,0.7,1),
                 userAlphaSpending=c(0.01,0.015,0.025),alpha.low=1e-10)
}
\keyword{group-sequential}
\keyword{critical values}
\keyword{efficacy boundary}

