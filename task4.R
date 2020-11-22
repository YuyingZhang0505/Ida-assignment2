load("dataex4.Rdata")
require(maxLik)

EM_function = function(data,beta_0,eps){
  #the starting point for beta
  beta=beta_0
  diff=1
  while(diff>eps){
    beta.old=beta
    
    # E-step, the expression of Q-function
    Q_function=function(beta,data){
      x=data[ ,1];y=data[ ,2]
      beta0=beta[1];beta1=beta[2]
      Q=0
      for (i in 1:length(x)) {
        Q=Q-log(1+exp(beta0+x[i]*beta1))
        # the value of y is missing
        if (is.na(y[i])==TRUE) {
          Q=Q+(beta0+beta1*x[i])*exp(beta.old[1]+x[i]*beta.old[2])/(1+exp(beta.old[1]+x[i]*beta.old[2]))
          # the value of y is not missing
        } else {
          Q=Q+y[i]*(beta0+beta1*x[i])
        }
      }
      Q
    }
    
    # M-step
    mle = maxLik(logLik=Q_function,data=dataex4,start=beta)
    beta0=mle$estimate[1];beta1=mle$estimate[2]
    beta=c(beta0,beta1)
    diff=sum(abs(beta-beta.old))
    
  }
  return(beta)
}

res=EM_function(data=dataex4,beta_0=c(0,0),eps=0.00001)
beta0est=res[1];beta1est=res[2]
#the maximum likelihood estimate of beta0 is 0.9755; beta1 is -2.4804
beta0est;beta1est