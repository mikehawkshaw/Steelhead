#Model to estimate 50% run timing for steelhead over time series

setwd("C:/DFO-MPO/github/Steelhead/Model")
source("directories.R")
library("xtable")
library(R2OpenBUGS)
setwd(data_dir)

#Read in and set up data
albion_annual<-read.table("steelhead_albion.csv",header=T, sep=",")

years<-seq(1995,1995+n_col-2)

n_col<-dim(albion_annual)[2]
N<-n_col-1

m<-rep(0,N)
s<-rep(0,N)

for(i in 2:n_col) {
  #normalize y to be a probability (sum = 1)
  p = albion_annual[,i]/sum(albion_annual[,i],na.rm=T)
  #compute weighted mean and standard deviation
  m[i-1]=sum(albion_annual[,1]*p,na.rm=T)
  s[i-1]=sqrt(sum((albion_annual[,1] - m[i-1])^ 2*p,na.rm=T))
}

data <- list("N","m","s")
inits<- function()
    list (theta=rnorm(N,0,100),mu.theta=rnorm(1,0,100),sigma.theta=runif(1,0,100))
parameters<-c("theta","mu.theta", "sigma.theta")

runtiming<- bugs(data,inits,parameters,model.file="runtimingHBM_OpenBUGS.txt",n.chains=3,n.iter=10000)

runtiming$summary

# Make results available so can call individual parameter sims (sims.list)
attach.bugs(runtiming)
pth<-dnorm(albion_annual[,1],mean(mu.theta),sd(sigma.theta))
par(mar=c(4,4,2,1))
plot(albion_annual[,1],p,type="b", col="blue",main="50% estimate of steelhead run timing",xlab="Julian Date", ylab="Corrected Albion Catch",bty="l")
lines(albion_annual[,1],pth, col="dark red",lwd=2)

#Get run timing estimates
mean(mu.theta)
sd(sigma.theta)