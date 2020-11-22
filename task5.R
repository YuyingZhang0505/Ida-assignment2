load("dataex5.Rdata")

mixture_density=function(data,theta0,eps) {
  n=length(data)
  theta=theta0
  p=theta[1];mu=theta[2];sig=sqrt(theta[3]);lamda=theta[4]
  diff=1
  while(diff>eps){
    theta.old=theta
    
    # E-step
    ptilde1=p*dlnorm(data,mean=mu,sd=sig)
    ptilde2=(1-p)*dexp(data,rate=lamda)
    ptilde=ptilde1/(ptilde1+ptilde2)
    
    # M-step
    p=mean(ptilde)
    mu=sum(ptilde*log(data))/sum(ptilde)
    sig=sqrt(sum(ptilde*(log(data)-mu)^2)/sum(ptilde))
    lamda=sum(1-ptilde)/sum((1-ptilde)*data)
    theta=c(p,mu,sig,lamda)
    diff=sum(abs(theta-theta.old))
  }
  return(theta)
}

res=mixture_density(data=dataex5,theta0=c(0.1,1,0.5^2,2),eps=0.00001)
pest=res[1];muest=res[2];sigest=res[3];lamdaest=res[4]
#the maximum likelihood estimate of p is 0.4796;mu is 2.0131;sigma is 0.9294;lamda is 1.0331
pest;muest;sigest;lamdaest

# plot
hist(dataex5, main = "The mixture distribution for random sample Y",
     xlab = "Y", 
     ylab = "Density", 
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
     breaks=35,
     freq = F, ylim = c(0,0.15))
curve(pest*dlnorm(x, mean = muest, sd = sigest) + (1 - pest)*dexp(x, rate=lamdaest),
      add = TRUE, xlim=c(0,120), lwd = 2, col = "red")