rm(list=ls())
library(DEoptim)
load("C:/Users/User/Documents/CASdatasets/data/swautoins.rda")
swautoins
y=swautoins$Payment
y=y[y>0]
z=swautoins$Claims
z=z[z>0]
x=y/z
x=x/100

#x=swaggclaim

# x=swautoins$Payment / swautoins$Insured
# x=x[x>0]
# x=x/100

min(x)
max(x)
mean(x)
std(x)
#this data is not the same as the used data in the paper.
#The mean and the range are different.
n=length(x)
g=function(x,beta,theta){
t=0
for (i in 1:10){
t1=((beta*x)^(i-1))*(theta^i)/(factorial(i)*gamma(i))
t=t+t1
}
return(t)
}
g=Vectorize(g,"x")
#g(2,2,3)
#g(x,2,3)



L=function(alpha){
beta=alpha[1]
theta=alpha[2]
l=sum(log(beta*exp(-beta*x)*exp(-theta)*g(x,beta,theta)/(1-exp(-theta))))
return(-l)
}
h=DEoptim(L,lower=c(0.0001,0.0001),upper=c(100,100),control = DEoptim.control(iter=500,trace = T))

	par.hat=h$optim$bestmem
names(par.hat)=NULL
	betahat=par.hat[1]
      thetahat=par.hat[2]
      
	value=-h$optim$bestval
      AIC=4-2*value

      par.hat
value
AIC

+-----------------------------------


hist(x,prob=T,ylim=c(0,0.02),xlab=expression(x))
f=function(x,beta,theta){
beta*exp(-beta*x)*exp(-theta)*g(x,beta,theta)/(1-exp(-theta))
}
#curve(f(x, 0.06979151, 3.52653845))
curve(f(x,betahat,thetahat),add=T,lty=1,col=2)


