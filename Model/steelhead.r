#steelhead.r
#exposure model for steelhead

#read in data
source("directories.R")
library("xtable")
setwd(data_dir)

###########################################
#Read in fishery openings data and format
###########################################

fishery_array<- array(as.numeric(NA), dim = c(521,3336,5,10)) #row, column, fishery, year

yr=2004

for(i in 1:10){ #Only goes to 2013 right now!! Need run timing for later years.
  
AreaB<-as.matrix(read.csv(paste(yr,"Area B_openings.csv",sep="")))
AreaD<-as.matrix(read.csv(paste(yr,"Area D_openings.csv",sep="")))
AreaE<-as.matrix(read.csv(paste(yr,"Area E_openings.csv",sep="")))
AreaG<-as.matrix(read.csv(paste(yr,"Area G_openings.csv",sep="")))
AreaH<-as.matrix(read.csv(paste(yr,"Area H_openings.csv",sep="")))
  
fishery_array[,,1,i]<-AreaB[,2:3337]
fishery_array[,,2,i]<-AreaD[,2:3337]
fishery_array[,,3,i]<-AreaE[,2:3337]
fishery_array[,,4,i]<-AreaG[,2:3337]
fishery_array[,,5,i]<-AreaH[,2:3337]

yr=yr+1
}

colnames(fishery_array)<-NULL

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
exposure<-array(as.numeric(NA),dim=c(n_fish,5,13))
speeds<-rep(0,n_fish)
passage_date<-rep(0,n_fish)

yr=2004 #re-initialize year variable

#Loop through each year
for(y in 1:10){
#Get day of year for July 15 of year of interest (season start)
seasonstart_doy <- as.numeric(strftime(paste(yr,"-07-15",sep=""), format = "%j"))

################################################################################################
#Population characteristics (these are the hypothesis about the population that will be tested)
################################################################################################

#Run-timing of the population. Based on mean and SD calculated in file "badestimator.r"
rt_mean<-subset(sh_runtiming$mean,sh_runtiming$year==yr)-seasonstart_doy #subtract season start day to put in correct position in matrix
rt_sd<-subset(sh_runtiming$sd,sh_runtiming$year==yr)

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

for(f in 1:5) {
  
   #Loop through each fish

   for(ind in 1:n_fish){
  
     exposure[ind,f,y]<-0
     for(loc in 1:494){ #494 is km where Albion located
    
       passage_time<-passage_hour[ind]
       time_at_loc<-passage_time-(494-loc)/speeds[ind]
      
       #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
       #If the fish is in the area before the time we care about, then obviously it's not exposed to any fisheries. May want to
       #edit later so that we can run it longer. Will need to make the opening matrices larger.
       if (time_at_loc<0 | time_at_loc>3335){
          exposure[ind,f,y]<-exposure[ind,f,y]
       } else{
         exposure[ind,f,y]<-exposure[ind,f,y]+fishery_array[loc+1,round(time_at_loc)+1,f,y]
       }
     }
   }
}

#################################################
#Move fish FORWARD from Albion through fisheries
#################################################

#Loop through each fishery

for(f in 1:5){

   #Loop through each fish

   for(ind in 1:n_fish)
   {

     for(loc in 495:515){ #From Albion to Mission, not including Albion start

     passage_time<-passage_hour[ind]
     time_at_loc<-passage_time+(loc-494)/speeds[ind]

     #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
     #If the fish is in the area after the time we care about, then obviously it's not exposed to any fisheries.
     if (time_at_loc>3335){
       exposure[ind,f,y]<-exposure[ind,f,y]
      } else {
       exposure[ind,f,y]<-exposure[ind,f,y]+fishery_array[loc+1,round(time_at_loc)+1,f,y]
      }
     }
   }
}
yr=yr+1
}

##############################
#Print plots to pdf file
##############################

setwd(plots_dir)

#save(exposure, file="2014AreaE.Rdata")

pdf(file = "Exposure by Fishery.pdf")

#plot(density(exposure))
par(mfcol=c(5,5))
par(mar=c(2,2,1,1))
yr=2004
for(y in 1:10){
hist(exposure[,1,y], breaks=50, main=paste("Area B ",yr),xlim=range(0,200))
hist(exposure[,2,y], breaks=50, main=paste("Area D ",yr), xlim=range(0,200))
hist(exposure[,3,y], breaks=50, main=paste("Area E ",yr), xlim=range(0,200))
hist(exposure[,4,y], breaks=50, main=paste("Area G ",yr), xlim=range(0,200))
hist(exposure[,5,y], breaks=50, main=paste("Area H ",yr), xlim=range(0,200))
yr=yr+1
}
dev.off()

pdf(file = "Exposure vs Passage Time.pdf")

par(mfcol=c(5,5))
par(mar=c(2,2,1,1))
yr=2004
for(y in 1:10){
  #Exposure to fisheries compared to run timing
  #xlab="Passage hour at Albion"
plot(passage_hour,exposure[,1,y], main=paste("Area B ",yr))
plot(passage_hour,exposure[,2,y], main=paste("Area D ",yr))
plot(passage_hour,exposure[,3,y], main=paste("Area E ",yr))
plot(passage_hour,exposure[,4,y], main=paste("Area G ",yr))
plot(passage_hour,exposure[,5,y], main=paste("Area H ",yr))
yr=yr+1
}
dev.off()
