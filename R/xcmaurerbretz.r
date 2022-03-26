xcmaurerbretz=function(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4)),
                         informationm=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4),
                         spending=rep("OBF",4),param.spending=rep(1,4),
                         alpha=0.025,direction=-1,graphin=BonferroniHolm(nrow(xm)),alpha.low=1e-10,retrospective=0){
  #xm: matrix of test statistics, assumed to be multivariate normal, each row is for one hypothesis at different time points, for each row, alpha levels must be non-decreasing
  #informationm: information for each endpoints
  #spending: type of spending function for each endpoint
  #param.spending: parameter in the spending function
  #retrospective: 0 (default) only compares the current test statistic with the updated critical value, 1 compares all the test statistics up to the current one with the updated critical values. Even though retrospectively looking at the values are statsitical valid in terms of control the type-1 error, not retrospectively looking at the past comparisons avoids the dilemma of retrospectively inflating the alpha level.  
  
  n=nrow(xm) #number of endpoints
  s=ncol(xm) #number of analyses (interims+final)
  mseq=seq(1,n,by=1)
  pm=xm
  crit=matrix(0,nrow=n,ncol=s)
  Hrej=rep(0,n)
  aset=NULL
  g=graphin
  
  galpha=as.numeric(alpha*g@weights)
  for (n1 in 1:n){
    seq.alpha=xcspending(alpha=galpha[n1],fractions=informationm[n1,],family=spending[n1],rho=param.spending[n1])$aseq
    seq.alpha=pmin(seq.alpha,rep(galpha[n1],s))
    crit[n1,]=xccrit(direction=direction,alpha=galpha[n1],informationRates=informationm[n1,],
                      userAlphaSpending = seq.alpha,alpha.low=alpha.low)$crit
  }
  
  if (direction==-1){
    for (s0 in 1:s){
      aset0=aset1=mseq[Hrej>0]
      ec=1;iter=0
      while (ec==1&iter<nrow(xm)){
        iter=iter+1
        if (retrospective==1) {ax=rowSums(as.matrix(xm[,1:s0])<=as.matrix(crit[,1:s0]))}
        else if (retrospective==0) {ax=as.numeric(xm[,s0]<=crit[,s0])}
        if (sum(ax)==0){ec=0}
        else{
          aset1=union(aset0,mseq[ax>0])
          adiff=setdiff(aset1,aset0)
          an=length(adiff)
          aset0=aset1
          if (an==0){ec=0}
          else{Hrej[adiff]=s0;for (j in 1:an)g=gMCP::rejectNode(g, paste0('H',adiff[j]))
          galpha=as.numeric(alpha*g@weights)
          for (n1 in 1:n){
            if (!n1 %in%aset0){
              seq.alpha=xcspending(alpha=galpha[n1],fractions=informationm[n1,],
                                    family=spending[n1],rho=param.spending[n1])$aseq
              seq.alpha=pmin(seq.alpha,rep(galpha[n1],s))
              crit[n1,]=xccrit(direction=direction,alpha=galpha[n1],
                                informationRates=informationm[n1,],
                                userAlphaSpending = seq.alpha,alpha.low=alpha.low)$crit
            }
          }
          }
        }
      }
    }
  }
  else if (direction==1){
    for (s0 in 1:s){
      aset0=aset1=mseq[Hrej>0]
      ec=1;iter=0
      while (ec==1&iter<nrow(xm)){
        iter=iter+1
        if (retrospective==1) {ax=rowSums(as.matrix(xm[,1:s0])>=as.matrix(crit[,1:s0]))}
        else if (retrospective==0) {ax=as.numeric(xm[,s0]>=crit[,s0])}
        if (sum(ax)==0){ec=0}
        else{
          aset1=union(aset0,mseq[ax>0])
          adiff=setdiff(aset1,aset0)
          an=length(adiff)
          aset0=aset1
          if (an==0){ec=0}
          else{Hrej[adiff]=s0;for (j in 1:an)g=gMCP::rejectNode(g, paste0('H',adiff[j]))
          galpha=as.numeric(alpha*g@weights)
          for (n1 in 1:n){
            if (!n1 %in%aset0){
              seq.alpha=xcspending(alpha=galpha[n1],fractions=informationm[n1,],
                                    family=spending[n1],rho=param.spending[n1])$aseq
              seq.alpha=pmin(seq.alpha,rep(galpha[n1],s))
              crit[n1,]=xccrit(direction=direction,alpha=galpha[n1],
                                informationRates=informationm[n1,],
                                userAlphaSpending = seq.alpha,alpha.low=alpha.low)$crit
            }
          }
          }
        }
      }
    }
  }
  else if (direction==0){
    for (s0 in 1:s){
      aset0=aset1=mseq[Hrej>0]
      ec=1;iter=0
      while (ec==1&iter<nrow(xm)){
        iter=iter+1
        if (retrospective==1) {ax=rowSums(as.matrix(abs(xm[,1:s0]))>=as.matrix(crit[,1:s0]))}
        else if (retrospective==0) {ax=as.numeric(abs(xm[,s0])>=crit[,s0])}
        if (sum(ax)==0){ec=0}
        else{
          aset1=union(aset0,mseq[ax>0])
          adiff=setdiff(aset1,aset0)
          an=length(adiff)
          aset0=aset1
          if (an==0){ec=0}
          else{Hrej[adiff]=s0;for (j in 1:an)g=gMCP::rejectNode(g, paste0('H',adiff[j]))
          galpha=as.numeric(alpha*g@weights)
          for (n1 in 1:n){
            if (!n1 %in%aset0){
              seq.alpha=xcspending(alpha=galpha[n1],fractions=informationm[n1,],
                                    family=spending[n1],rho=param.spending[n1])$aseq
              seq.alpha=pmin(seq.alpha,rep(galpha[n1],s))
              crit[n1,]=xccrit(direction=direction,alpha=galpha[n1],
                                informationRates=informationm[n1,],
                                userAlphaSpending = seq.alpha,alpha.low=alpha.low)$crit
            }
          }
          }
        }
      }
    }
  }
  rejected=NULL
  if (sum(Hrej)>0) rejected=sort(mseq[Hrej>0])
  cumdecisionsm=matrix(0,nrow=n,ncol=s)
  for (s0 in 1:s){
    cumdecisionsm[Hrej<=s0&Hrej>=1,s0]=1
  }
  decisionsm=cumdecisionsm[,s]
  
  list(Hrej=Hrej,rejected=rejected,decisionsm=decisionsm,cumdecisionsm=cumdecisionsm)
}
