\name{xccalgsp}
\alias{xccalgsp}
\title{Calculate group-sequential p-values}
\description{This will calculate group-sequential p-values}
\usage{
xccalgsp(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
                alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                informationm=matrix(rep(c(0.4,0.8,1),each=2),ncol=3,nrow=2),direction=-1)
                
}
\arguments{
  \item{xm}{a matrix of test statistics for each endpoint (in row) and each interim (in column)}
  \item{alpham}{a matrix of alpha spending for the statistics \code{xm}}
  \item{informationm}{a matrix of information fractions for the statistics \code{xm}}
  \item{direction}{-1: (one-sided)reject if test stat is smaller than the critical value;
                    1: (one-sided)reject if test stat is smaller than the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
}
\details{}
\value{\item{pm}{group-sequential p-values}
       \item{critm}{critical values}
}
\author{Xiaodong Luo}
\note{}
\examples{
xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4))
im=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4)
alpham=matrix(0,nrow=4,ncol=3)
for (i in 1:4){
   alpham[i,]=xcspending(alpha=0.025,fractions=im[i,],family="OBF",rho=(i+1)/2)$aseq
}
xccalgsp(xm=xm,alpham=alpham,informationm=im,direction=-1)
}
\keyword{group-sequential p-values}

