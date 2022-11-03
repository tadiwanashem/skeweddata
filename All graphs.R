g=function(x,beta,theta){
  t=0
  for (i in 1:10){
    t1=((beta*x)^(i-1))*(theta^i)/(factorial(i)*gamma(i))
    t=t+t1
  }
  return(t)
}
g=Vectorize(g,"x")


f=function(x,beta,theta){
  
  beta*exp(-beta*x)*exp(-theta)*g(x,beta,theta)/(1-exp(-theta))
}




dbga <- function(x,beta,theta) {
  f <- beta*exp(-beta*x)*exp(-theta)*g(x,beta,theta)/(1-exp(-theta))
  return(f)
}
library(latex2exp)

xrange=seq(0, 20, 0.001)

plot(xrange, f(xrange, 1, 0.0001), type="l", col="blue", lwd=4, ylab="PDF", xlab = "x")
lines(xrange, f(xrange, 1, 0.005), type="l", col="green", lwd=3)
lines(xrange, f(xrange, 1, 1), type="l", col="red", lwd=3)
lines(xrange, f(xrange, 1, 2), type="l", col="yellow", lwd=3)
lines(xrange, f(xrange, 1, 5), type="l", col="pink", lwd=3)

legend("topright", legend = c("ZTP-EXP(1,0.0001)", "ZTP-EXP(1,0.005)", "ZTP-EXP(1,1)", "ZTP-EXP(1,2)", "ZTP-EXP(1,5)"), 
       fill =c("blue", "green", "red", "yellow", "pink"), cex=0.9)


