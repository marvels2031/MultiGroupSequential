\name{xccrit}
\alias{xccrit}
\title{Calculate critical values}
\description{This utility function calculates the critical values}
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
  \item{alpha.low}{default is 1e-10, if allocated alpha is smaller than this number, the corresponding critical value will be set to infinity}
}
\details{}
\value{\item{crit}{critical values}
}
\author{Xiaodong Luo}
\note{}
\examples{
xccrit(direction=-1,alpha=0.025,informationRates=c(0.4,0.7,1),
                 userAlphaSpending=c(0.01,0.015,0.025),alpha.low=1e-10)
}
\keyword{group-sequential}
\keyword{critical values}
\keyword{efficacy boundary}

