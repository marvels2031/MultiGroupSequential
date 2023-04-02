\name{xcspending}
\alias{xcspending}
\title{Calculate alpha spending}
\description{This utility function calculates alpha spending. Note that the OBF and Pocock spending functions are not the originally proposed ones, they are the modified ones that are closely resemble the original versions. That being said, you might still see some differences}
\usage{
xcspending(alpha,fractions=seq(0.2,1,by=0.2),family="OBF",rho=1)
}
\arguments{
  \item{alpha}{overall familywise error rate}
  \item{fractions}{information fractions}
  \item{family}{family of spending functions, one of "OBF", "pocock", "power"}
  \item{rho}{parameter of the spending function}
}
\details{
OBF: \eqn{2\{1-\Phi(\Phi^{-1}(1-\alpha/2)/t^{\rho/2})\}};
pocock: \eqn{\alpha \log\{1+(e-1)*t\}};
power: \eqn{\alpha*t^{\rho}}
}
\value{\item{aseq}{alpha spending}
}
\author{Xiaodong Luo}
\note{}
\examples{
xcspending(alpha=0.025,fractions=seq(0.2,1,by=0.2),family="OBF",rho=1)
}
\keyword{group-sequential}
\keyword{critical values}
\keyword{efficacy boundary}
\keyword{alpha spending}


