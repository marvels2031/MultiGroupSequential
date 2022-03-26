\name{xcseqhh}
\alias{xcseqhh}
\title{Sequential Hochberg/Hommel procedure for multiple endpoints based on q-values}
\description{Sequential Hochberg/Hommel procedure for multiple endpoints based on q-values}
\usage{
xcseqhh(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),
        alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
        epsilon=1.0e-10,precision=10,method='Hochberg')          
          
}
\arguments{
  \item{pm}{a matrix of group-sequential p-values for different endpoints (in row) at different times (in column)}
  \item{alpham}{a matrix of alpha spending corresponding to the p-values \code{pm}}
  \item{epsilon}{lower bound for the alpha}
  \item{precision}{precision of the values}
  \item{method}{"Hochberg" or "Hommel"}
}
\details{}
\value{\item{rejected}{the index set of rejected hypotheses}
       \item{decisionsm}{rejection decision for each endpoint (row) at each timepoint (column)}
       \item{cumdecisionsm}{cumulative rejection decision for each endpoint (row) at each timepoint (column)}
       \item{alphaused}{alpha levels actually used for each endpoint (row) at each timepoint (column)}
       
}
\author{Xiaodong Luo}
\note{This provides the calculation for the variance. }
\examples{
pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)
alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2)
xcseqhh(pm=pm,alpham=alpham,method="Hochberg")
xcseqhh(pm=pm,alpham=alpham,method="Hommel")
}
\keyword{group-sequential}
\keyword{Hochberg}
\keyword{Hommel}
\keyword{q-values}


