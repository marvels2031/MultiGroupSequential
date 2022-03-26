xccrit=function(direction=-1,alpha=0.025,informationRates=c(0.4,0.7,1),
                 userAlphaSpending=c(0.01,0.015,0.025),alpha.low=1e-10){
  
  crit=rep(1,length(userAlphaSpending))
  if (alpha<1e-06&alpha>alpha.low)alpha=1e-06
  else if (alpha<=alpha.low){alpha=0;crit=rep(Inf,length(userAlphaSpending))*direction}
  
  if (alpha>alpha.low){
    if (direction!=0){
      design=rpact::getDesignGroupSequential(sided = 1, alpha=alpha, 
                                             informationRates = informationRates,
                                             typeOfDesign = "asUser",
                                             userAlphaSpending = userAlphaSpending)
      crit=direction*design$criticalValues
    }
    else if (direction==0){
      design=rpact::getDesignGroupSequential(sided = 2, alpha=alpha, 
                                             informationRates = informationRates,
                                             typeOfDesign = "asUser",
                                             userAlphaSpending = userAlphaSpending)
      crit=design$criticalValues
    }
  }
  list(crit=crit)
}
