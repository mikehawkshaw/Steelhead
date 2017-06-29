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
#annual parameters
  real<lower=0.5,upper=4> rt_m[n_years];
  real<lower=1e-9,upper=1e-2> rt_sd[n_years];
} 

model {

#hyper_parameters

  mu_rt_m~uniform();
  sig_rt_m~uniform();
  mu_rt_sd~uniform();
  sig_rt_sd~uniform();


#annual parameters
  for(y in 1:n_years){
    rt_m[y]~normal(mu_rt_m,sig_rt_m)
    rt_sd[y]~normal(mu_rt_sd,sig_rt_sd)

}

#model



}
generated quantities {
}
