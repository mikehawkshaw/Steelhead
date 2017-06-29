#Model to estimate 50% run timing for steelhead, both over time series and year-specific.

source("directories.R")
library("xtable")
library(rstan)
setwd(data_dir)

#Read in and set up data
albion_annual<-read.table("steelhead_albion.csv",header=T, sep=",")

#Change back to model directory
setwd(data_dir)


