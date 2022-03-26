## sequential Graphical procedure based on group-sequential p-values
xcseqgxgs=function(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),alpha=0.025,graphin=BonferroniHolm(2)){
  #pm: group-sequential p-values
  #alpha: FWER
  n=nrow(pm)
  s=ncol(pm)
  mseq=seq(1,n,by=1)
  decisionsm=alphaused=cumdecisionsm=matrix(0,nrow=n,ncol=s)
  
  t=1
  abc=gMCP(graphin, pvalues=pm[,t], alpha=alpha)
  g=as.matrix(abc@rejected)[,1]
  rejected=mseq[g]
  decisionsm[,t]=cumdecisionsm[,t]=g
  if (s>1){
     for (t in 2:s){
       #graph2 <- setWeights(graphin, c(1/6, 1/2, 1/3))
       #abc=gMCP(graph2, pvalues=pm[,t], alpha=alpha) # Having different graphs for different interims will not protect the FWER
       abc=gMCP(graphin, pvalues=pm[,t], alpha=alpha)
       g=as.matrix(abc@rejected)[,1]
       rejected=union(rejected,mseq[g])
       decisionsm[,t]=g
       cumdecisionsm[,t]=(g|cumdecisionsm[,t-1])
     }
  }
  if (length(rejected)>0) rejected=sort(rejected)
  list(rejected=rejected,decisionsm=decisionsm,cumdecisionsm=cumdecisionsm)
}