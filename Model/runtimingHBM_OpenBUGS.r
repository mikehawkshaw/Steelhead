#Model to estimate 50% run timing for steelhead over time series
setwd("C:/DFO-MPO/github/Steelhead/Model")
source("directories.R")
library("xtable")
library("reshape")
library("stringr")
library(R2OpenBUGS)
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
  
  if(nchar(mod_catch[i,1]==2)){
  mod_catch[i,1]<-str_sub(mod_catch[i,1],-1)
  }
  mod_catch[i,1]<-str_sub(mod_catch[i,1],-2)
}

mod_catch[,1]<-as.numeric(mod_catch[,1])

colnames(mod_catch)<-c("year","catch")

#Change back to model directory
setwd(model_dir)

#Prep inputs to model
dat<-list(n_days=n_days, n_years=n_years, obs_abund=catch)
inits<-list(list(mu_rt_m=290,sig_rt_m=20,mu_rt_sd=15,sig_rt_sd=8),list(mu_rt_m=280,sig_rt_m=22,mu_rt_sd=13,sig_rt_sd=7),list(mu_rt_m=300,sig_rt_m=10,mu_rt_sd=10,sig_rt_sd=9))
parameters<-c("mu_rt_m","sig_rt_m","mu_rt_sd","sig_rt_sd","mu.year","sigma.year")

RT_fit<- bugs(data=dat,inits=inits,parameters.to.save=parameters,model.file="runtimingHBM_OpenBUGSv4.txt",n.chains=3,n.iter=5000,debug=T)

RT_fit$summary

# Make results available so can call individual parameter sims (sims.list)
attach.bugs(RT_fit)
pth<-dnorm(albion_annual[,1],mean(mu),sd(sigma.year))
par(mar=c(4,4,2,1))
plot(albion_annual[,1],p,type="b", col="blue",main="50% estimate of steelhead run timing",xlab="Julian Date", ylab="Corrected Albion Catch",bty="l")
lines(albion_annual[,1],pth, col="dark red",lwd=2)

#Get run timing estimates
mean(mu.theta)
sd(sigma.theta)