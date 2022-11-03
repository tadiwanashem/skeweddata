rm(list = ls())
library(plyr) 
library(mcmc) 
library(gibbs.met) 
library(LearnBayes) 
library(coda) 
library(ggmcmc) 
library(MCMCvis) 
library(ggplot2) 
library(pracma) 
library(DEoptim)

#-------------------------------------
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

dbga2 <- function(x,beta,theta) {
  
  f <-beta*exp(-beta*x)*exp(-theta)*g(x,beta,theta)/(1-exp(-theta))
  return(f)
}

logdbga <- function(x,beta,theta) {
  
  f <- beta*exp(-beta*x)*exp(-theta)*g(x,beta,theta)/(1-exp(-theta))
  f=log(f)
  id <- (x < 0) | (x > 20)
  f[id] <- -Inf
  return(f)
}

logdbga2 <- function(x,beta,theta) {
  f <- beta*exp(-beta*x)*exp(-theta)*g(x,beta,theta)/(1-exp(-theta))
  logl=log(f)
  return(logl)
}


mloglbga0 <- function(param, x) {
  beta <- param[1]
  theta <- param[2]
  m <- -sum(logdbga2(x,beta,theta))
  return(m)
}



fitbga <- function(param, x) {
  
  h <- DEoptim(mloglbga0,
               x = x, lower = c(0, 0),
               upper = c(20,20), control = DEoptim.control(iter = 500, trace = F)
  )
  par.hat <- c(h$optim$bestmem)
  bias <- par.hat - param
  Mse <- (par.hat - param)^2
  names(bias) <- paste("bias", names(bias), sep = "_")
  names(Mse) <- paste("Mse", names(Mse), sep = "_")
  out <- c(par.hat, bias, Mse)
  return(out)
}



sim_mle_bga <- function(n, Bp, x0, beta,theta) {
  p00 <- c(beta=beta,theta=theta)
  x=rbga_rwmetrop(n = n, Bp = Bp, x0 = x0, beta=beta,theta=theta)
  #x=rbga_gibbs_met(n = n, Bp = Bp, x0 = x0, beta=beta,theta=theta)
  fit <- fitbga(p00, x)
  return(fit)
}

#----------------------

rbga_rwmetrop <- function(n, Bp, x0, beta,theta) {
  proposal <- list(var = diag(1), scale = 1)
  xy <- rwmetrop(
    logpost = logdbga, proposal = proposal, start = x0, m = n + Bp,
    beta=beta,theta=theta
  )
  v <- var(xy$par)
  proposal <- list(var = v, scale = 1)
  xy <- rwmetrop(
    logpost = logdbga, proposal = proposal, start = x0, m = n + Bp,
    beta=beta,theta=theta
  )
  xy <- xy$par[-(1:Bp), ]
  
  return(xy)
}


rbga_gibbs_met <- function(n, Bp, x0, beta,theta) {
  xy_mcmc <- gibbs_met(
    log_f = logdbga, no_var = 1, ini_value = x0, iters = n +
      Bp, iters_met = 10, stepsizes_met = c(0.5, 0.5), beta=beta,theta=theta
  )
  xy <- xy_mcmc[-(1:(Bp + 1)), ]
  
  return(xy)
}


# simulation --------------------------------------------------------------
n <- 1000

Bp <- 10000  #convergence 10000

x0 <- 1


beta <- 1
theta <- 0.5

xmax <- 100


p00 <- c(beta=beta,theta=theta)

dbga(x = x0, beta=beta,theta=theta)
dbga2(x = x0, beta=beta,theta=theta)
logdbga(x = x0, beta=beta,theta=theta)

#approach 1

x=rbga_rwmetrop(n = n, Bp = Bp, x0 = x0, beta=beta,theta=theta)

# approach 2
#x=rbga_gibbs_met(n = n, Bp = Bp, x0 = x0, beta=beta,theta=theta)

hist(x,prob=T,xlab=expression(paste(theta)),main="",col="white",ylim=c(0,0.6))
curve(dbga(x, beta,theta),add=T,ylim=c(0,1))



#MLE, bias and mse------------------------------------

sim1 <- sim_mle_bga(n, Bp, x0, beta,theta)
sim1

repl <- 10
sim2 <- raply(repl, sim_mle_bga(n, Bp, x0, beta,theta), .progress = "text")
sim2
colMeans(sim2)

# --------------------------------




