\name{xccalgspcor}
\alias{xccalgspcor}
\title{Calculate group-sequential p-values given the correlation matrix of the test statistics}
\description{This will calculate group-sequential p-values given the correlation matrix of the test statistics}
\usage{
xccalgspcor(xm=qnorm(c(0.03,0.04,0.01)),
                   alpham=c(0.02,0.03,0.05),
                   corrm=diag(length(xm)),direction=-1,tol=1e-10)
                
}
\arguments{
  \item{xm}{a vector of test statistics at each analysis}
  \item{alpham}{a vector of alpha spending for the statistics \code{xm}}
  \item{corrm}{correlation matrix of the statistics \code{xm}}
  \item{direction}{-1: (one-sided)reject if test stat is smaller than or equal to the critical value;
                    1: (one-sided)reject if test stat is greater than or equal to the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
  \item{tol}{accuracy tolerance when calculating the quantiles}
}
\details{}
\value{\item{pm}{group-sequential p-values}
       \item{critm}{critical values}
}
\author{Xiaodong Luo}
\note{}
\examples{
xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2))
ir=c(0.4,0.8,1)
corrm=diag(length(ir))
for (i in 1:length(ir))for(j in 1:length(ir))corrm[i,j]=sqrt(ir[pmin(i,j)]/ir[pmax(i,j)])
xccalgsp(xm=xm)$critm[1,]
xccalgspsim(xm=xm)$critm[1,]
xccalgspcor(xm=xm[1,],corrm=corrm)$critm
}
\keyword{group-sequential p-values}

