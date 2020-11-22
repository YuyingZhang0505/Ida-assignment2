load("dataex2.Rdata")
require(maxLik)

log_like=function(mu,data){
  sum=0
  x=data[ ,1]
  R=data[ ,2]
  for (i in 1:length(data[ ,1])) {
    sum=sum+R[i]*dnorm(x[i],mean=mu,sd=1.5,log=TRUE)+(1-R[i])*pnorm(x[i],mean=mu,sd=1.5,log.p=TRUE)
  }
  sum
}

# the maximum likelihood estimate of mu found by maxLik function is mu=5.5328
mle = maxLik(logLik=log_like,data=dataex2,start=c(6))
summary(mle)
mle$estimate