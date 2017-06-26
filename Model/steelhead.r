#steelhead.r
#exposure model for steelhead

#read in data
source("directories.R")
library("xtable")
setwd(data_dir)

###########################################
#Read in fishery openings data and format
###########################################

fishery_mat<-as.matrix(read.csv("2013Area B_openings.csv"))
fishery_mat<-as.matrix(read.csv("2013Area D_openings.csv"))
fishery_mat<-as.matrix(read.csv("2013Area E_openings.csv"))
fishery_mat<-as.matrix(read.csv("2013Area G_openings.csv"))
fishery_mat<-as.matrix(read.csv("2013Area H_openings.csv"))

colnames(fishery_mat)<-NULL
fishery_mat<-fishery_mat[,2:3337]

sh_runtiming<-as.data.frame(read.csv("steelhead_runtiming.csv", header=T))

n_km<-521
n_hours<-3336

###################################
#set up a fake steelhead population
#IBM like model
###################################

n_fish<-1000
fish<-seq(1,n_fish,by=1)

#each fish has characteristics and they are in these vectors
exposure<-rep(NA,n_fish)
speeds<-rep(0,n_fish)
passage_date<-rep(0,n_fish)

#Get day of year for July 15 of year of interest (season start)
yr="2013" #put into variable so it can be made dynamic later when looping
seasonstart_doy <- as.numeric(strftime(paste(yr,"-07-15",sep=""), format = "%j"))

################################################################################################
#Population characteristics (these are the hypothesis about the population that will be tested)
################################################################################################

#Run-timing of the population. Based on mean and SD calculated in file "badestimator.r"
rt_mean<-subset(sh_runtiming$mean,sh_runtiming$year=="2013")-seasonstart_doy #subtract season start day to put in correct position in matrix
rt_sd<-subset(sh_runtiming$sd,sh_runtiming$year=="2013")

#passage_date = the date that the fish passes Albion
passage_date<-(pmax(30,pmin(140,rnorm(fish,rt_mean,rt_sd))))
passage_hour<-passage_date*24 #convert to hours. Hour 0 = midnight July 15

#Speed that the fish travel. Assumptions based on speed of other salmonids.
speed_mean<-20 #km/day
speed_sd<-3
speeds<-(pmax(9,pmin(55,rnorm(fish,speed_mean,speed_sd))))/24 #km/hr

#################################################
#Move fish BACKWARD from Albion through fisheries
#################################################

#Loop through each fishery

#Loop through each fish

for(ind in 1:n_fish)
{
  
  exposure[ind]<-0
  for(loc in 1:494){ #494 is km where Albion located
    
    passage_time<-passage_hour[ind]
    time_at_loc<-passage_time-(494-loc)/speeds[ind]
    
    #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
    
    exposure[ind]<-exposure[ind]+fishery_mat[loc,round(time_at_loc)]
    
  }
}

#################################################
#Move fish FORWARD from Albion through fisheries
#################################################

#Loop through each fishery

#Loop through each fish

for(ind in 1:n_fish)
{

exposure[ind]<-0
for(loc in 495:n_km){ #From Albion upstream, not including Albion start

passage_time<-passage_hour[ind]
time_at_loc<-passage_time+(loc-5494)/speeds[ind]

#check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery

exposure[ind]<-exposure[ind]+fishery_mat[loc,round(time_at_loc)]

}
}

##############################
#Print plots to pdf file
##############################

setwd(plots_dir)

#save(exposure, file="2014AreaE.Rdata")

pdf()

plot(density(exposure))
dev.off()
