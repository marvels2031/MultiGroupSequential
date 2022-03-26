## xccalgsp: a function to calculate the group-sequential p-values for each endpoint
xccalgsp=function(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
                alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                informationm=matrix(rep(c(0.4,0.8,1),each=2),ncol=3,nrow=2),direction=-1){
  #xm: matrix of test statistics, assumed to be multivariate normal
  #alpham: matrix of alpha-levels, each row is for one hypothesis at different time points, for each row, alpha levels must be non-decreasing
  # direction=0 two-sided, reject if the absolute value of test stat >= the critical value
  # direction=-1 one-sided, reject if the test stat <= the critical value
  # direction=1 one-sided, reject if the test stat >= the critical value
  n=nrow(xm)
  s=ncol(xm)
  mseq=seq(1,n,by=1)
  critm=pm=xm
  
  sided=ifelse(direction==0,2,1) #direction=0 is two-sided test
  
  for (j in 1:n){
    design =rpact::getDesignGroupSequential(sided = sided, alpha=alpham[j,s], 
                                            informationRates = informationm[j,],
                                            typeOfDesign = "asUser",
                                            userAlphaSpending = alpham[j,])
    crit=design$criticalValues
    if (direction==-1)crit=-design$criticalValues
    critm[j,]=crit
    ir=informationm[j,]
    nir=length(ir)
    dcorr=diag(nir)
    
    for (u in 1:nir){
      for (v in u:nir){
        dcorr[u,v]=dcorr[v,u]=sqrt(ir[u]/ir[v])
      }
    }
    
    if (direction==-1){
      aa=rep(0,s)
      aa[1]=pnorm(xm[j,1])
      pm[j,1]=(xm[j,1]<=crit[1])*aa[1]+(xm[j,1]>crit[1])*1
      if (s>=3){
        for (i in 2:(s-1)){
          aa[i]=OpenMx::omxMnor(covariance=dcorr[1:i,1:i], means=rep(0,i), lbound=c(crit[1:(i-1)],-Inf), ubound=c(rep(Inf,i-1),xm[j,i]))+alpham[j,i-1]
          #aa[i]=mvtnorm::pmvnorm(lower=c(crit[1:(i-1)],-Inf),upper=c(rep(Inf,i-1),xm[j,i]),sigma=dcorr[1:i,1:i])+alpham[j,i-1]
          pm[j,i]=pm[j,i-1]-prod(xm[j,1:(i-1)]>crit[1:(i-1)])*(xm[j,i]<=crit[i])*(1-aa[i])
        }
      }
      aa[s]=OpenMx::omxMnor(covariance=dcorr[1:s,1:s], means=rep(0,s), lbound=c(crit[1:(s-1)],-Inf), ubound=c(rep(Inf,s-1),xm[j,s]))+alpham[j,s-1]
      #aa[s]=mvtnorm::pmvnorm(lower=c(crit[1:(s-1)],-Inf),upper=c(rep(Inf,s-1),xm[j,s]),sigma=dcorr[1:s,1:s])+alpham[j,s-1]
      pm[j,s]=pm[j,s-1]-prod(xm[j,1:(s-1)]>crit[1:(s-1)])*(1-aa[s])
    }
    else if (direction==1){
      aa=rep(0,s)
      aa[1]=1-pnorm(xm[j,1])
      pm[j,1]=(xm[j,1]>=crit[1])*aa[1]+(xm[j,1]<crit[1])*1
      if (s>=3){
        for (i in 2:(s-1)){
          aa[i]=OpenMx::omxMnor(covariance=dcorr[1:i,1:i], means=rep(0,i), lbound=c(rep(-Inf,i-1),xm[j,i]), ubound=c(crit[1:(i-1)],Inf))+alpham[j,i-1]
          #aa[i]=mvtnorm::pmvnorm(lower=c(rep(-Inf,i-1),xm[j,i]),upper=c(crit[1:(i-1)],Inf),sigma=dcorr[1:i,1:i])+alpham[j,i-1]
          pm[j,i]=pm[j,i-1]-prod(xm[j,1:(i-1)]<crit[1:(i-1)])*(xm[j,i]>=crit[i])*(1-aa[i])
        }
      }
      aa[s]=OpenMx::omxMnor(covariance=dcorr[1:s,1:s], means=rep(0,s), lbound=c(rep(-Inf,s-1),xm[j,s]),ubound=c(crit[1:(s-1)],Inf))+alpham[j,s-1]
      #aa[s]=mvtnorm::pmvnorm(lower=c(rep(-Inf,s-1),xm[j,s]),upper=c(crit[1:(s-1)],Inf),sigma=dcorr[1:s,1:s])+alpham[j,s-1]
      pm[j,s]=pm[j,s-1]-prod(xm[j,1:(s-1)]<crit[1:(s-1)])*(1-aa[s])
    }
    else if (direction==0){
      aa=rep(0,s)
      aa[1]=2*(1-pnorm(abs(xm[j,1])))
      pm[j,1]=(abs(xm[j,1])>=crit[1])*aa[1]+(abs(xm[j,1])<crit[1])*1
      if (s>=3){
        for (i in 2:(s-1)){
          aa[i]=OpenMx::omxMnor(covariance=dcorr[1:i,1:i], means=rep(0,i), lbound=c(-crit[1:(i-1)],abs(xm[j,i])), ubound=c(crit[1:(i-1)],Inf))
          aa[i]=aa[i]*2+alpham[j,i-1]
          #aa[i]=mvtnorm::pmvnorm(lower=c(-crit[1:(i-1)],abs(xm[j,i])),upper=c(crit[1:(i-1)],Inf),sigma=dcorr[1:i,1:i])
          #aa[i]=aa[i]*2+alpham[j,i-1]
          pm[j,i]=pm[j,i-1]-prod(abs(xm[j,1:(i-1)])<crit[1:(i-1)])*(abs(xm[j,i])>=crit[i])*(1-aa[i])
        }
      }
      aa[s]=OpenMx::omxMnor(covariance=dcorr[1:s,1:s], means=rep(0,s), lbound=c(-crit[1:(s-1)],abs(xm[j,s])),ubound=c(crit[1:(s-1)],Inf))
      aa[s]=aa[s]*2+alpham[j,s-1]
      #aa[s]=mvtnorm::pmvnorm(lower=c(-crit[1:(s-1)],abs(xm[j,s])),upper=c(crit[1:(s-1)],Inf),sigma=dcorr[1:s,1:s])
      #aa[s]=aa[s]*2+alpham[j,s-1]
      pm[j,s]=pm[j,s-1]-prod(abs(xm[j,1:(s-1)])<crit[1:(s-1)])*(1-aa[s])
    }
  }
  list(pm=pm,critm=critm)
}


