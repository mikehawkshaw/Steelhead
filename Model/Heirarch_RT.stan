#Heirarch_RT.stan

data {
  int<lower=1> n_days;
  int<lower=0> n_years;
  real index[n_days,n_years];
  real catch[n_days,n_years];
}
transformed data{

}
parameters {
#hyper_parameters
  real<lower=0.5,upper=4> mu_rt_m;
  real<lower=0.1,upper=5> sig_rt_m;
  real<lower=0.5,upper=4> mu_rt_sd;
  real<lower=0.1,upper=5> sig_rt_sd;
  real annual[n_years];
  real  [n_years];  
#annual parameters
  real<lower=0.5,upper=4> rt_m[n_years];
  real<lower=1e-9,upper=1e-2> rt_sd[n_years];
} 

model {

  real pred_abundance[n_days,n_years];

#hyper_parameters

  mu_rt_m~uniform(1,365);
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
    pred_abundance[d,y]=1/(1+exp(-1.7*(d-rt_m[y])/(rt_sd[y])));
    if(fishery_index[d,y]==1){catch~normal(pred_abundance[d,y]*annual[y],sigma)}
}
}

}
generated quantities {
}
