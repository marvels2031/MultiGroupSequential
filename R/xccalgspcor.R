xccalgspcor=function(xm=qnorm(c(0.03,0.04,0.01)),
                   alpham=c(0.02,0.03,0.05),
                   corrm=diag(length(xm)),direction=-1,tol=1e-10){
  #xm: A vector of test statistics at different (ascending) timepoints, assumed to be multivariate normal
  #alpham: A vector of cumulative alpha levels at different (ascending) timepoints alpha levels must be non-decreasing
  #corrm: correlation matrix of the test statistics
  # direction=0 two-sided, reject if the absolute value of test stat >= the critical value
  # direction=-1 one-sided, reject if the test stat <= the critical value
  # direction=1 one-sided, reject if the test stat >= the critical value
  s=length(xm)
  critm=pm=xm
  
  if (direction==-1){
    aa=rep(0,s)
    aa[1]=pnorm(xm[1])
    critm[1]=qnorm(alpham[1])
    pm[1]=(aa[1]<=alpham[1])*aa[1]+(aa[1]>alpham[1])*1
    if (s>=3){
      for (i in 2:(s-1)){
        aa[i]=OpenMx::omxMnor(covariance=corrm[1:i,1:i], means=rep(0,i), lbound=c(critm[1:(i-1)],-Inf), ubound=c(rep(Inf,i-1),xm[i]))+alpham[i-1]
        pm[i]=pm[i-1]-prod(xm[1:(i-1)]>critm[1:(i-1)])*(aa[i]<=alpham[i])*(1-aa[i])
        xf=function(x){ 
          ax=OpenMx::omxMnor(covariance=corrm[1:i,1:i], means=rep(0,i), lbound=c(critm[1:(i-1)],-Inf), ubound=c(rep(Inf,i-1),x))+alpham[i-1]
          (ax-alpham[i])^2
        }
        xmin <- optimize(xf, c(-20, 0), tol = tol)
        critm[i]=xmin$minimum
      }
    }
    aa[s]=OpenMx::omxMnor(covariance=corrm[1:s,1:s], means=rep(0,s), lbound=c(critm[1:(s-1)],-Inf), ubound=c(rep(Inf,s-1),xm[s]))+alpham[s-1]
    pm[s]=pm[s-1]-prod(xm[1:(s-1)]>critm[1:(s-1)])*(1-aa[s])
    xf=function(x){ 
      ax=OpenMx::omxMnor(covariance=corrm[1:s,1:s], means=rep(0,s), lbound=c(critm[1:(s-1)],-Inf), ubound=c(rep(Inf,s-1),x))+alpham[s-1]
      (ax-alpham[s])^2
    }
    xmin <- optimize(xf, c(-20, 0), tol = tol)
    critm[s]=xmin$minimum
  }
  else if (direction==1){
    aa=rep(0,s)
    aa[1]=1-pnorm(xm[1])
    critm[1]=qnorm(1-alpham[1])
    pm[1]=(aa[1]<=alpham[1])*aa[1]+(aa[1]>alpham[1])*1
    if (s>=3){
      for (i in 2:(s-1)){
        aa[i]=OpenMx::omxMnor(covariance=corrm[1:i,1:i], means=rep(0,i), lbound=c(rep(-Inf,i-1),xm[i]), ubound=c(critm[1:(i-1)],Inf))+alpham[i-1]
        pm[i]=pm[i-1]-prod(xm[1:(i-1)]<critm[1:(i-1)])*(aa[i]<=alpham[i])*(1-aa[i])
        xf=function(x){ 
          ax=OpenMx::omxMnor(covariance=corrm[1:i,1:i], means=rep(0,i), lbound=c(rep(-Inf,i-1),x), ubound=c(critm[1:(i-1)],Inf))+alpham[i-1]
          (ax-alpham[i])^2
        }
        xmin <- optimize(xf, c(0, 20), tol = tol)
        critm[i]=xmin$minimum
      }
    }
    aa[s]=OpenMx::omxMnor(covariance=corrm[1:s,1:s], means=rep(0,s), lbound=c(rep(-Inf,s-1),xm[s]), ubound=c(critm[1:(s-1)],Inf))+alpham[s-1]
    pm[s]=pm[s-1]-prod(xm[1:(s-1)]<critm[1:(s-1)])*(1-aa[s])
    xf=function(x){ 
      ax=OpenMx::omxMnor(covariance=corrm[1:s,1:s], means=rep(0,s), lbound=c(rep(-Inf,s-1),x), ubound=c(critm[1:(s-1)],Inf))+alpham[s-1]
      (ax-alpham[s])^2
    }
    xmin <- optimize(xf, c(0, 20), tol = tol)
    critm[s]=xmin$minimum
  }
  else if (direction==0){
    aa=rep(0,s)
    aa[1]=2*(1-pnorm(abs(xm[1])))
    critm[1]=qnorm(1-alpham[1]/2)
    pm[1]=(aa[1]<=alpham[1])*aa[1]+(aa[1]>alpham[1])*1
    if (s>=3){
      for (i in 2:(s-1)){
        aa[i]=OpenMx::omxMnor(covariance=corrm[1:i,1:i], means=rep(0,i), lbound=c(-critm[1:(i-1)],abs(xm[i])), ubound=c(critm[1:(i-1)],Inf))*2+alpham[i-1]
        pm[i]=pm[i-1]-prod(abs(xm[1:(i-1)])<critm[1:(i-1)])*(aa[i]<=alpham[i])*(1-aa[i])
        xf=function(x){ 
          ax=OpenMx::omxMnor(covariance=corrm[1:i,1:i], means=rep(0,i), lbound=c(-critm[1:(i-1)],x), ubound=c(critm[1:(i-1)],Inf))*2+alpham[i-1]
          (ax-alpham[i])^2
        }
        xmin <- optimize(xf, c(0, 20), tol = tol)
        critm[i]=xmin$minimum
      }
    }
    aa[s]=OpenMx::omxMnor(covariance=corrm[1:s,1:s], means=rep(0,s), lbound=c(-critm[1:(s-1)],abs(xm[s])), ubound=c(critm[1:(s-1)],Inf))*2+alpham[s-1]
    pm[s]=pm[s-1]-prod(abs(xm[1:(s-1)])<critm[1:(s-1)])*(1-aa[s])
    xf=function(x){ 
      ax=OpenMx::omxMnor(covariance=corrm[1:s,1:s], means=rep(0,s), lbound=c(-critm[1:(s-1)],x), ubound=c(critm[1:(s-1)],Inf))*2+alpham[s-1]
      (ax-alpham[s])^2
    }
    xmin <- optimize(xf, c(0, 20), tol = tol)
    critm[s]=xmin$minimum
  }
  list(pm=pm,critm=critm)
}

