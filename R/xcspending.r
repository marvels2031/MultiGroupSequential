## Note that the OBF and Pocock spending functions are not the originally proposed ones, they are the modified ones that are closely resemble the original versions. That being said, you might still see some differences 
xcspending=function(alpha,fractions=seq(0.2,1,by=0.2),family="OBF",rho=1){
  if (family=='OBF'){
    qa=qnorm(1-alpha/2)/fractions^(rho/2)
    aseq=2*(1-pnorm(qa))
  }
  else if (family=='pocock'){
    qa=1+(exp(1)-1)*fractions
    aseq=alpha*log(qa)
  }
  else if (family=='power'){
    aseq=alpha*fractions^rho
  }
  list(aseq=aseq)
}
