## calgspsim: a function to calculate the group-sequential p-values for each endpoint via simulation
xccalgspsim=function(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
                   alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                   informationm=matrix(rep(c(0.4,0.8,1),each=2),ncol=3,nrow=2),
                   r.seed=rep(17,2),nsample=1e+6,direction=-1){
  #xm: matrix of test statistics, assumed to be multivariate normal
  #alpham: matrix of alpha-levels, each row is for one hypothesis at different time points, for each row, alpha levels must be non-decreasing
  n=nrow(xm)
  s=ncol(xm)
  mseq=seq(1,n,by=1)
  critm=pm=xm
  
  for (j in 1:n){
    abc=xccalgspsim1(xm=xm[j,],alpham=alpham[j,],informationm=informationm[j,],
                   r.seed=r.seed[j]+j,nsample=nsample,direction=direction)
    pm[j,]=abc$p.value.gs
    critm[j,]=abc$crit.value
  }
  list(pm=pm,critm=critm)
}


