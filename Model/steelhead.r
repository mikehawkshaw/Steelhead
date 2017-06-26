#steelhead.r
#exposure model for steelhead

#read in data
source("directories.R")
library("xtable")
setwd(data_dir)

fishery_mat<-as.matrix(read.csv("2013Area E_openings.csv"))
fishery_mat<-as.matrix(read.csv("2013Area B_openings.csv"))

colnames(fishery_mat)<-NULL
fishery_mat<-fishery_mat[,2:3337]

sh_runtiming<-as.data.frame(read.csv("steelhead_runtiming.csv", header=T))

n_km<-521
n_hours<-3336

#set up a fake steelhead population
#IBM like model

n_fish<-1000
fish<-seq(1,n_fish,by=1)

#each fish has characteristics and they are in these vectors
exposure<-rep(NA,n_fish)
speeds<-rep(0,n_fish)
passage_date<-rep(0,n_fish)

#Get day of year for July 15 of year of interest (season start)
yr="2013" #put into variable so it can be made dynamic later when looping
seasonstart_doy <- as.numeric(strftime(paste(yr,"-07-15",sep=""), format = "%j"))

#population characteristics (these are the hypothesis about the population that will be tested)
rt_mean<-subset(sh_runtiming$mean,sh_runtiming$year=="2013")-seasonstart_doy
rt_sd<-subset(sh_runtiming$sd,sh_runtiming$year=="2013")
#passage date = the date that the fish passes Albion
passage_date<-(pmax(30,pmin(140,rnorm(fish,rt_mean,rt_sd))))
passage_hour<-passage_date*24

speed_mean<-20
speed_sd<-3
speeds<-(pmax(9,pmin(55,rnorm(fish,speed_mean,speed_sd))))/24	#speed in km/h

#move fish though fisheries 

for(ind in 1:n_fish)
{

exposure[ind]<-0
for(loc in 1:n_km){

start_time<-passage_hour[ind]
time_at_loc<-start_time+loc*speeds[ind]

#check exposure against fishery matrix

exposure[ind]<-exposure[ind]+fishery_mat[loc,round(time_at_loc)]

}
}


setwd(plots_dir)

#save(exposure, file="2014AreaE.Rdata")

pdf()

plot(density(exposure))
dev.off()
