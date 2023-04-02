\name{xccalgspsim1}
\alias{xccalgspsim1}
\title{Calculate group-sequential p-values via simulation}
\description{This utility function will be called by function "xccalgspsim" to calculate group-sequential p-values via simulation for single endpoint}
\usage{
xccalgspsim1(xm=qnorm(c(0.03,0.04,0.01)),alpham=c(0.02,0.03,0.05),
                      informationm=c(0.4,0.8,1),
                      r.seed=17,nsample=1e+6,direction=0)
}
\arguments{
  \item{xm}{test statistics}
  \item{alpham}{alpha spending}
  \item{informationm}{information fractions}
  \item{r.seed}{random seed}
  \item{nsample}{number of random samples}
  \item{direction}{-1: (one-sided)reject if test stat is smaller than or equal to the critical value;
                    1: (one-sided)reject if test stat is greater than or equal to the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
}
\details{}
\value{\item{crit.value}{critical values}
       \item{p.value.gs}{group-sequential p-values}
       \item{xm}{test statistics}
       \item{alpham}{alpha spending}
       \item{informationm}{information fractions}
}
\author{Xiaodong Luo}
\note{This provides the calculation for the group-sequential p-values in case there is an issue in using rpact package.}
\examples{
xm=qnorm(c(0.03,0.04,0.01))
im=c(0.4,0.8,1)
alpham=xcspending(alpha=0.025,fractions=im,family="OBF",rho=2)$aseq
xccalgspsim1(xm=xm,alpham=alpham,informationm=im,direction=-1)
}
\keyword{group-sequential p-values}

