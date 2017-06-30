#Model to estimate 50% run timing for steelhead, both over time series and year-specific.

source("directories.R")
library("xtable")
library(rstan)
setwd(data_dir)

#model

stanmodelcode <- "
#Heirarch_RT.stan

data {
int<lower=1> n_days;
int<lower=0> n_years;
int index[n_days,n_years];
int obs_catch[n_days,n_years];
}
transformed data{

}
parameters {
#hyper_parameters
real<lower=1,upper=365> mu_rt_m;
real<lower=1,upper=100> sig_rt_m;
real<lower=1,upper=50> mu_rt_sd;
real<lower=1,upper=50> sig_rt_sd;
real<lower=1,upper=100> annual[n_years];
real sigma;  
#annual parameters
real<lower=1,upper=365> rt_m[n_years];
real<lower=1,upper=100> rt_sd[n_years];
} 

model {

real pred_abundance[n_days,n_years];

#hyper_parameters

mu_rt_m~normal(290,30);
sig_rt_m~uniform(1,100);
mu_rt_sd~uniform(1,50);
sig_rt_sd~uniform(1,50);
annual~uniform(1e-6,1e6);

#annual parameters
for(y in 1:n_years){
rt_m[y]~normal(mu_rt_m,sig_rt_m);
rt_sd[y]~normal(mu_rt_sd,sig_rt_sd);

}

#model

for (y in 1:n_years){
for(d in 1:n_days)
{
pred_abundance[d,y]=1/(1+exp(-1.7*((d+244)-rt_m[y])/(rt_sd[y])));
if(index[d,y]==1){obs_catch[d,y]~poisson(pred_abundance[d,y]*annual[y]);}
}
}

}
generated quantities {
}
"



#Read in and set up data
albion_annual<-read.table("steelhead_albion.csv",header=T, sep=",")

n_days<-dim(albion_annual)[1]
n_years<-dim(albion_annual)[2]-1
catch<-as.matrix(albion_annual[,2:20])
fisheries_index<-matrix(1,103,19)
fisheries_index[is.na(catch)]<-0
catch[is.na(catch)]<-99
#Change back to model directory
setwd(model_dir)

dat<-list(n_days=n_days, n_years=n_years, index=fisheries_index, obs_catch=catch)
inits<-list(list(mu_rt_m=290,sig_rt_m=20,mu_rt_sd=15,sig_rt_sd=8),list(mu_rt_m=280,sig_rt_m=22,mu_rt_sd=13,sig_rt_sd=7),list(mu_rt_m=300,sig_rt_m=10,mu_rt_sd=10,sig_rt_sd=9))

fit <- stan(model_code = stanmodelcode, data = dat, iter = 5000, chains = 3,  verbose = TRUE) 







