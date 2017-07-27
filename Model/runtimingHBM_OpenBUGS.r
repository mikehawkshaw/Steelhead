#Model to estimate 50% run timing for steelhead over time series
setwd("C:/DFO-MPO/github/Steelhead/Model")
source("directories.R")
library(xtable)
library(reshape)
library(stringr)
library(R2OpenBUGS)
library(mcmcplots)
library(coda)

setwd(data_dir)

#Read in and set up data
albion_annual<-read.table("steelhead_albion.csv",header=T, sep=",")

years<-seq(1,1+n_col-2)

n_days<-dim(albion_annual)[1]
n_col<-dim(albion_annual)[2]
n_years<-n_col-1

annual_sums<-colSums(albion_annual, na.rm=T)

catch<-as.matrix(albion_annual[,2:n_col])

for(i in 1:dim(catch)[2])
{
  catch[,i]<-catch[,i]/annual_sums[1+i]
}

colSums(catch, na.rm=T)

#year_index<-matrix(1,n_days,n_years)
#for(i in 1:n_years){
#  year_index[,i]<-years[i]
#}

colnames(catch)<-seq(1,n_years,by=1)

mod_catch <- melt(data.frame(catch))

mod_catch[,1]<-as.character(mod_catch[,1])

for(i in 1:dim(mod_catch)[1]){
  
  if(nchar(mod_catch[i,1])==2){
  mod_catch[i,1]<-str_sub(mod_catch[i,1],-1)
  }
  mod_catch[i,1]<-str_sub(mod_catch[i,1],-2)
}

mod_catch[,1]<-as.numeric(mod_catch[,1])

colnames(mod_catch)<-c("year","catch")
mod_catch$day<- rep(seq(1,n_days,by=1),n_years)

#Change back to model directory
setwd(model_dir)

#Prep inputs to model
dat<-list(n_days=n_days,n_years=n_years,obs_abund=catch)
inits<-list(list(mu_rt_m=290,tau_rt_m=1,mu_rt_sd=15,tau_rt_sd=1,tau.c=1),list(mu_rt_m=280,tau_rt_m=1,mu_rt_sd=13,tau_rt_sd=1,tau.c=1),list(mu_rt_m=300,tau_rt_m=1,mu_rt_sd=10,tau_rt_sd=1,tau.c=1))
parameters<-c("rt_m","rt_sd","mu_rt_m","sigma_rt_m","mu_rt_sd","sigma_rt_sd","sigma")

RT_fit<-bugs(data=dat,inits=inits,parameters.to.save=parameters,model.file="runtimingHBM_OpenBUGS.txt",n.chains=3,n.burnin=0,n.iter=80000,debug=F)

# Make results available so can call individual parameter sims (sims.list)
attach.bugs(RT_fit)
print(RT_fit)
caterplot(RT_fit,parms=c("rt_m","mu_rt_m"),reorder=F)
caterplot(RT_fit,parms=c("rt_sd","mu_rt_sd"),reorder=F)

#----------------Alternate option for getting trace plots into R
RT_fit<-bugs(data=dat,inits=inits,parameters.to.save=parameters,model.file="runtimingHBM_OpenBUGS.txt",n.chains=3,n.iter=80000,debug=F,codaPkg=T)

## Produce a CODA object from the lineout output
line.coda <- read.bugs(RT_fit)
## Summary statistics from the CODA output (uses coda package)
summary(line.coda)
## to change line.coda to pretend that only 1 chain
line.coda2 <- as.mcmc(as.matrix(line.coda))
#Summarize across chains
HPDinterval(line.coda2)
## Use CODA to produce diagnostic plots for the MCMC
plot(line.coda)
traceplot(line.coda)
#---------------------------------------

p = albion_annual[,i]/sum(albion_annual[,i],na.rm=T)

pth<-dnorm(albion_annual[,1],mean(mu_rt_m),mean(mu_rt_sd))
par(mar=c(4,4,2,1))
plot(albion_annual[,1],p,type="b", col="blue",main="50% estimate of steelhead run timing",xlab="Julian Date", ylab="Corrected Albion Catch",bty="l")
lines(albion_annual[,1],pth, col="dark red",lwd=2)

#Get run timing estimates
mean(mu_rt_m)
mean(mu_rt_sd)