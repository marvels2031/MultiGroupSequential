xccalgspsim1=function(xm=qnorm(c(0.03,0.04,0.01)),alpham=c(0.02,0.03,0.05),
                      informationm=c(0.4,0.8,1),
                      r.seed=17,nsample=1e+6,direction=0){
  # direction=0 two-sided, equal tails
  # direction=-1 one-sided, less than or equal to critical value
  # direction=1 one-sided, greater than or equal to critical value
  s=length(xm) #number of analyses
  crit.value=rep(0,s)
  p.value=p.value.gs=rep(0,s)
  
  if (direction==-1){
    crit.value[1]=qnorm(alpham[1])
    p.value[1]=pnorm(xm[1])
  }
  else if (direction==0){
    crit.value[1]=qnorm(1-alpham[1]/2)
    p.value[1]=2*pnorm(-abs(xm[1]))
  }
  else if (direction==1){
    crit.value[1]=qnorm(1-alpham[1])
    p.value[1]=pnorm(xm[1],lower.tail = FALSE)
  }
  p.value.gs[1]=ifelse(p.value[1]<=alpham[1],p.value[1],1)
  
  if (s>1){
    set.seed(r.seed)
    rmatrix=matrix(0,nrow=nsample,ncol=s)
    rindex=rep(1,nsample)
    rmatrix[,1]=sqrt(informationm[1])*rnorm(nsample)
    for (j in 2:s)rmatrix[,j]=rmatrix[,j-1]+sqrt(informationm[j]-informationm[j-1])*rnorm(nsample)
    rmatrix=rmatrix%*%diag(1/sqrt(informationm))
    ai=1
    for (j in 2:s){
      #j=3
      alp=(alpham[j]-alpham[j-1])/(1-alpham[j-1])
      if (direction==-1){
        rindex=rindex*(rmatrix[,j-1]>crit.value[j-1])
        crit.value[j]=quantile(rmatrix[rindex==1,j],prob=alp)
        ai=ai*(xm[j-1]>crit.value[j-1]) 
        ptemp=sum(rmatrix[rindex==1,j]<=xm[j])/nsample+alpham[j-1]
        if (j<s){
          p.value.gs[j]=p.value.gs[j-1]-ai*(xm[j]<=crit.value[j])*(1-ptemp)
        }
        else {
          p.value.gs[j]=p.value.gs[j-1]-ai*(1-ptemp)
        }
      }
      else if (direction==0){
        rindex=rindex*(abs(rmatrix[,j-1])<crit.value[j-1])
        crit.value[j]=quantile(abs(rmatrix[rindex==1,j]),prob=(1-alp/2))
        ai=ai*(abs(xm[j-1])<crit.value[j-1])
        ptemp=sum(abs(rmatrix[rindex==1,j])>=abs(xm[j]))/nsample+alpham[j-1]
        if (j<s){
          p.value.gs[j]=p.value.gs[j-1]-ai*(abs(xm[j])>crit.value[j])*(1-ptemp)
        }
        else {
          p.value.gs[j]=p.value.gs[j-1]-ai*(1-ptemp)
        }
      }
      else if (direction==1){
        rindex=rindex*(rmatrix[,j-1]<crit.value[j-1])
        crit.value[j]=quantile(rmatrix[rindex==1,j],prob=alp,lower.tail=FALSE)
        ai=ai*(xm[j-1]<crit.value[j-1]) 
        ptemp=sum(rmatrix[rindex==1,j]>=xm[j])/nsample+alpham[j-1]
        if (j<s){
          p.value.gs[j]=p.value.gs[j-1]-ai*(xm[j]>=crit.value[j])*(1-ptemp)
        }
        else {
          p.value.gs[j]=p.value.gs[j-1]-ai*(1-ptemp)
        }
      }
    }
  }
  list(crit.value=crit.value,p.value.gs=p.value.gs,xm=xm,alpham=alpham,informationm=informationm)
}
