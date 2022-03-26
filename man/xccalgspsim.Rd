\name{xccalgspsim}
\alias{xccalgspsim}
\title{Calculate group-sequential p-values via simulation}
\description{This will calculate group-sequential p-values via simulation}
\usage{
xccalgspsim(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
                   alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                   informationm=matrix(rep(c(0.4,0.8,1),each=2),ncol=3,nrow=2),
                   r.seed=rep(17,2),nsample=1e+6,direction=-1)
}
\arguments{
  \item{xm}{a matrix of test statistics for each endpoint (in row) and each interim (in column)}
  \item{alpham}{a matrix of alpha spending for the statistics \code{xm}}
  \item{informationm}{a matrix of information fractions for the statistics \code{xm}}
  \item{r.seed}{random seeds for each endpoints}
  \item{nsample}{number of random samples}
  \item{direction}{-1: (one-sided)reject if test stat is smaller than the critical value;
                    1: (one-sided)reject if test stat is smaller than the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
}
\details{Let \eqn{\tau_{min}} and \eqn{\tau_{max}} correspond to the minimum follow-up time \code{taumin} and the maximum follow-up time \code{taumax}. Let \eqn{T_f}, \eqn{C}, \eqn{E} and \eqn{R} be the follow-up time, the drop-out time, the study entry time and the total recruitment period(\eqn{R} is the last element of \code{ut}). For type 2 follow-up \eqn{T_f=min(C,\tau_{max})}. For type 3 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max})}. For type 4 follow-up, \eqn{T_f=min(R+\tau_{min}-E,\tau_{max},C)}. Let \eqn{f} be the density of \eqn{T_f}. We calculate
\deqn{\int_0^{\infty} \{at/(1+at)\}f(t)dt.}
}
\value{\item{pm}{group-sequential p-values}
       \item{critm}{critical values}
}
\author{Xiaodong Luo}
\note{This provides the calculation for the group-sequential p-values in case there is an issue in using rpact package.}
\examples{
xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4))
im=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4)
alpham=matrix(0,nrow=4,ncol=3)
for (i in 1:4){
   alpham[i,]=xcspending(alpha=0.025,fractions=im[i,],family="OBF",rho=(i+1)/2)$aseq
}
xccalgspsim(xm=xm,alpham=alpham,informationm=im,r.seed=rep(17,4),direction=-1)
}
\keyword{group-sequential p-values}

