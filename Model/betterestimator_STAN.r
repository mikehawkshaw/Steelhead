#Model to estimate 50% run timing for steelhead, both over time series and year-specific.

source("directories.R")
library("xtable")
library(rstan)
library(ggplot2)
library(ggmcmc)
library(extrafont)

#extrafont::font_import() Getting a weird error message to do with default fonts - this fixed it but might be unique to my computer
extrafont::loadfonts() # need some fonts tloaded to stop errors w/ default font settings in Rmd

setwd(data_dir)

#model

stanmodelcode <- "
data {
int<lower=1> n_days;
int<lower=0> n_years;
int index[n_days,n_years];
real obs_catch[n_days,n_years];
}

parameters {
//hyper_parameters
real<lower=1,upper=360> mu_rt_m;
real<lower=1,upper=100> tau_rt_m;
real<lower=1,upper=50> mu_rt_sd;
real<lower=1,upper=50> tau_rt_sd;
real<lower=1e-9,upper=5> sigma;
  
//annual parameters
real<lower=1,upper=360> rt_m[n_years];
real<lower=1,upper=100> rt_sd[n_years];
} 
transformed parameters{
  real sig_rt_m;
  real sig_rt_sd;

  sig_rt_m = 1.0 / sqrt(tau_rt_m);
  sig_rt_sd = 1.0 / sqrt(tau_rt_sd);


}
model {

real pred_abundance[n_days,n_years];

//hyper_parameters

mu_rt_m~uniform(240,330);
tau_rt_m~gamma(0.001,0.001);
mu_rt_sd~uniform(1,50);
tau_rt_sd~gamma(1,50);


//annual parameters
for(y in 1:n_years){
rt_m[y]~normal(mu_rt_m,sig_rt_m);
rt_sd[y]~normal(mu_rt_sd,sig_rt_sd);

}

//model

for (y in 1:n_years){
for(d in 1:n_days)
{
//pred_abundance[d,y]=(1/(1+exp(-1.7*((d+244)-rt_m[y])/(rt_sd[y]))))-(1/(1+exp(-1.7*((d+243)-rt_m[y])/(rt_sd[y]))));

pred_abundance[d,y]=(1/(rt_sd[y]*pow((22.0/7.0)*2,0.5)))*exp(-pow(((d+243)-rt_m[y]),2)/(2*pow(rt_sd[y],2)));

if(index[d,y]==1){obs_catch[d,y]~normal(pred_abundance[d,y], sigma);}
}
}

}
generated quantities {
}
"



#Read in and set up data
setwd("C:/DFO-MPO/GitHub/Steelhead/Data")
albion_annual<-read.table("steelhead_albion.csv",header=T, sep=",")
annual_sums<-colSums(albion_annual, na.rm=T)
n_days<-dim(albion_annual)[1]
n_years<-dim(albion_annual)[2]-1
catch<-as.matrix(albion_annual[,2:20])
for(i in 1:dim(catch)[2])
{
  catch[,i]<-catch[,i]/annual_sums[1+i]
}
colSums(catch, na.rm=T)
fisheries_index<-matrix(1,103,19)
fisheries_index[is.na(catch)]<-0

for(i in 1:dim(catch)[2])
{
  catch[,i]<-catch[,i]*mean(fisheries_index[,i])
}


catch[is.na(catch)]<-999  #STAN hates NAs #Winbugs/OpenBugs it would be better to leave NA's
#Change back to model directory
setwd("C:/DFO-MPO/GitHub/Steelhead/Model")

dat<-list(n_days=n_days, n_years=n_years, index=fisheries_index, obs_catch=catch)
inits<-list(list(mu_rt_m=290,mu_rt_sd=15),list(mu_rt_m=280,mu_rt_sd=13),list(mu_rt_m=300,mu_rt_sd=10))

fit <- stan(model_code = stanmodelcode, data = dat, iter = 10000, chains = 3,  verbose = TRUE,control = list(adapt_delta = 0.99)) 
fit_ggobj<-ggs(fit)

mu_sd_summary <- summary(fit, pars = c("mu_rt_m", "sig_rt_m","mu_rt_sd", "sig_rt_sd"), probs = c(0.1,0.5, 0.9))$summary
print(mu_sd_summary) #use these values for the Steelhead Run Timing inputs to the steelhead.r file  
par(mfcol=c(2,2))
plot(fit, pars=c("mu_rt_m"))
plot(fit, pars=c("sig_rt_m"))
plot(fit, pars=c("mu_rt_sd"))
plot(fit, pars=c("sig_rt_sd"))

