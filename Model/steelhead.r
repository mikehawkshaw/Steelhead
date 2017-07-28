#steelhead.r
#Steelhead Predictors of In-river Timing and Exposure models

#read in data
source("directories.R")
library("xtable")
library("svMisc")
setwd(data_dir)

###########################################
#Read in fishery openings data and format
###########################################

data_source<-"FN" #Options: "Commercial", "FN", "REC"

if(data_source=="Commercial"){
  
fishery_array<- array(as.numeric(NA), dim = c(521,3336,5,13)) #row, column, fishery, year

yr=2004

for(i in 1:13){ #Updated to include up to 2016 since now using same run timing every year
  
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
n_km<-521
km_end<-515
n_fisheries<-5

}else if(data_source=="FN"){
  fishery_array<- array(as.numeric(NA), dim = c(625,3336,4,13)) #row, column, fishery, year
  
  yr=2004
  
  for(i in 1:13){
    
    APMBS<-as.matrix(read.csv(paste(yr,"APM BSEO_openings.csv",sep="")))
    APMDN<-as.matrix(read.csv(paste(yr,"APM DNEO_openings.csv",sep="")))
    APMSN<-as.matrix(read.csv(paste(yr,"APM SNEO_openings.csv",sep="")))
    BPMDN<-as.matrix(read.csv(paste(yr,"BPM DNEO_openings.csv",sep="")))
    
    fishery_array[,,1,i]<-APMBS[,2:3337]
    fishery_array[,,2,i]<-APMDN[,2:3337]
    fishery_array[,,3,i]<-APMSN[,2:3337]
    fishery_array[,,4,i]<-BPMDN[,2:3337]
    
    yr=yr+1
  }
  n_km<-624
  km_end<-624
  n_fisheries<-4
  
}else{ #data_source=="REC"
  
} 

colnames(fishery_array)<-NULL

sh_runtiming<-as.data.frame(read.csv("steelhead_runtiming.csv", header=T))

n_hours<-3336

###################################
#set up a fake steelhead population
#IBM like model
###################################
n_fish<-1000
n_reps<-100
fish<-seq(1,n_fish,by=1)

#each fish has characteristics and they are in these vectors
exposure<-array(0,dim=c(n_fish,n_fisheries,13,n_reps))
speeds<-rep(0,n_fish)
passage_date<-rep(0,n_fish)

################################################################################################
#Population characteristics (these are the hypotheses about the population that will be tested)
################################################################################################

#Start the clock
ptm <- proc.time()

for(i in 1:(n_reps)){
set.seed(i)
  
yr=2004 #re-initialize year variable

#Loop through each year
for(y in 1:13){
#Get day of year for July 15 of year of interest (season start)
seasonstart_doy <- as.numeric(strftime(paste(yr,"-07-15",sep=""), format = "%j"))

#Run timing based on Bayesian estimator (grand mean)
#rt_mean<-282.131681-seasonstart_doy #subtract season start day to put in correct position in matrix
#rt_mean_sd<-12.416564
#rt_sd<-16.510714
#rt_sd_sd<-5.970674

#Run-timing of the population. Based on mean and SD calculated with Bayesian estimator for each year.
rt_mean<-sh_runtiming$rt_mean[sh_runtiming$year==yr]-seasonstart_doy #subtract season start day to put in correct position in matrix
rt_mean_sd<-sh_runtiming$rt_mean_sd[sh_runtiming$year==yr]
rt_sd<-sh_runtiming$rt_sd[sh_runtiming$year==yr]
rt_sd_sd<-sh_runtiming$rt_sd_sd[sh_runtiming$year==yr]

#cumulative and daily proportions of the run vulnerable to each fishery
m_vec<-rnorm(n_reps,rt_mean,rt_mean_sd) 
s_vec<-rnorm(n_reps,rt_sd,rt_sd_sd)

#passage_date = the date that the fish passes Albion
passage_date<-(pmax(30,pmin(140,rnorm(fish,m_vec[i],s_vec[i]))))
passage_hour<-passage_date*24 #convert to hours. Hour 0 = midnight July 15

#Speed that the fish travel. Assumptions based on speed of other salmonids.
speed_mean<-20 #km/day
speed_sd<-3
speeds<-(pmax(9,pmin(55,rnorm(fish,speed_mean,speed_sd))))/24 #km/hr

#################################################
#Move fish BACKWARD from Albion through fisheries
#################################################

#Loop through each fishery

for(f in 1:n_fisheries) {
  
   #Loop through each fish

   for(ind in 1:n_fish){
     
     for(loc in 1:494){ #494 is km where Albion located
       
       time_at_loc<-passage_hour[ind]-(494-loc)/speeds[ind]
      
       #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
       #If the fish is in the area before the time we care about, then obviously it's not exposed to any fisheries. May want to
       #edit later so that we can run it longer. Will need to make the opening matrices larger.
       if (time_at_loc>0 & time_at_loc<3335){
         exposure[ind,f,y,i]<-exposure[ind,f,y,i]+fishery_array[loc+1,round(time_at_loc)+1,f,y]
       }
     }
   }
}

#################################################
#Move fish FORWARD from Albion through fisheries
#################################################

#Loop through each fishery

for(f in 1:n_fisheries){

   #Loop through each fish
   for(ind in 1:n_fish)
   {
     for(loc in 495:km_end){ #From Albion, not including Albion start
       
     time_at_loc<-passage_hour[ind]+(loc-494)/speeds[ind]

     #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
     #If the fish is in the area after the time we care about, then obviously it's not exposed to any fisheries.
     if (time_at_loc<3335){
       exposure[ind,f,y,i]<-exposure[ind,f,y,i]+fishery_array[loc+1,round(time_at_loc)+1,f,y]
      }
     }
   }
}
yr=yr+1
}
  progress(i,max.value=n_reps,init=1)
  Sys.sleep(0.01)
  if (i==n_reps) cat("Done!\n")
}
#Stop the clock
proc.time() - ptm

#Save iterations - change file name as appropriate
saveRDS(exposure,file="EO_exposure_1-100.RData")

##############################
#Manipulate exposure data
##############################

#Add multiple exposure runs together:
#This is not very dynamic but it is fine for now...

total_reps<-200

temp_exposure1<-readRDS("EO_exposure_1-100.RData")
temp_exposure2<-readRDS("EO_exposure_101-200.RData")

exposure<-array(as.numeric(NA),dim=c(n_fish,4,13,total_reps))

#There was a mistake made when re-opening the file, where exposure didn't get reset to 4 fisheries, so it had a 5th fishery
#still saved. The code below gets rid of that fifth set of data

for(i in 1:100){
 exposure[,,,i]<-temp_exposure1[,-5,,i]
}

for(i in 101:200){
  exposure[,,,i]<-temp_exposure2[,,,i-100]
}

#Get total exposure time by fishery
#Not sure this is a useful metric, but it's here in case we decide to use it
yr=2004

total_exposure<-array(as.numeric(NA),dim=c(n_fisheries,13))

for(y in 1:13){
  for(f in 1:n_fisheries){
    total_exposure[f,y]<-sum(exposure[,f,y]) #Incorrect num dimensions, old code
  }
  yr=yr+1
}

#-------------Get # of fish exposed by fishery

total_exposed<-array(as.numeric(NA),dim=c(n_fisheries,13,total_reps))
for(i in 1:total_reps){
 for(y in 1:13){
  for(f in 1:n_fisheries){
   total_exposed[f,y,i]<-sum(exposure[,f,y,i]>0)
   }
  }
}

#Convert to average #/% fish exposed by fishery each year
if(data_source=="Commercial"){
mean_exposed<-array(as.numeric(NA),dim=c(n_fisheries,13))
sd_exposed<-array(as.numeric(NA),dim=c(n_fisheries,13))

mean_perc_exposed<-array(as.numeric(NA),dim=c(n_fisheries,13))
sd_perc_exposed<-array(as.numeric(NA),dim=c(n_fisheries,13))

for(y in 1:13){
  for(f in 1:n_fisheries){
    mean_exposed[f,y]<-mean(total_exposed[f,y,])
    mean_perc_exposed[f,y]<-mean_exposed[f,y]/1000*100
    sd_exposed[f,y]<-sd(total_exposed[f,y,])
    sd_perc_exposed[f,y]<-sd_exposed[f,y]/1000*100
  }
}

#Change calculation to sum DN/SN exposure:
}else if(data_source=="FN"){
  total_exposed<-array(as.numeric(NA),dim=c(n_fisheries-1,13,n_reps))
  for(i in 1:n_reps){
    for(y in 1:13){
        total_exposed[1,y,i]<-sum(exposure[,1,y,i]>0)
        total_exposed[2,y,i]<-sum((exposure[,2,y,i]+exposure[,3,y,i])>0)
        total_exposed[3,y,i]<-sum(exposure[,4,y,i]>0)
    }
  }
  
  mean_exposed<-array(as.numeric(NA),dim=c(n_fisheries-1,13))
  sd_exposed<-array(as.numeric(NA),dim=c(n_fisheries-1,13))
  
  mean_perc_exposed<-array(as.numeric(NA),dim=c(n_fisheries-1,13))
  sd_perc_exposed<-array(as.numeric(NA),dim=c(n_fisheries-1,13))
  
  for(y in 1:13){
    for(f in 1:n_fisheries-1){
      mean_exposed[f,y]<-mean(total_exposed[f,y,])
      mean_perc_exposed[f,y]<-mean_exposed[f,y]/1000*100
      sd_exposed[f,y]<-sd(total_exposed[f,y,])
      sd_perc_exposed[f,y]<-sd_exposed[f,y]/1000*100
    }
  }  
}

#------------Get cumulative exposure to fisheries

if(data_source=="Commercial"){
  
  cml_exposure<-array(as.numeric(NA),dim=c(n_fish,13,n_reps))
  
  for(i in 1:n_reps){
    for(y in 1:13){
      for(n in 1:n_fish){
        cml_exposure[n,y,i]<-sum(exposure[n,,y,i]>0)
      }
    }
  }
}else if(data_source=="FN"){
  exposure_temp<-array(as.numeric(NA),dim=c(n_fish,n_fisheries-1,13,n_reps))
  
  exposure_temp[,1,,]<-exposure[,1,,]
  exposure_temp[,2,,]<-exposure[,2,,]+exposure[,3,,]
  exposure_temp[,3,,]<-exposure[,4,,]
  
  cml_exposure<-array(as.numeric(NA),dim=c(n_fish,13,n_reps))
  
  for(i in 1:n_reps){
    for(y in 1:13){
      for(n in 1:n_fish){
        cml_exposure[n,y,i]<-sum(exposure_temp[n,,y,i]>0)
      }
    }
  }
}else{ #data_source=="REC"
  
}

rm(exposure_temp) #takes up a lot of space, might as well remove

#------------Get iterative exposure to fisheries (Area B then D then H then E/BPM then APM…)


##############################
#Print plots to pdf file
##############################

setwd(plots_dir)

#Percentage of run exposed by year
###Need to edit
if(data_source=="Commercial"){
pdf(file="Population Percent Exposure by Year - Barplots.pdf")
par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(2,2,2,2), xpd=FALSE)

yr=2004
for(y in 1:13){
barplot(total_exposed[,y]/1000*100, main=yr, xlab="",ylab="", 
        names.arg=c("B", "D", "E", "G", "H"), axis.lty=1, ylim=range(0,100))
  yr=yr+1
 if(y==6 | y==12 | y==13){  
  mtext(text="Fishery Area",side=1,line=1,outer=TRUE)
  mtext(text="% Exposure",side=2,line=1,outer=TRUE)
  mtext(text="Population Percent Exposure by Year",side=3,line=1,outer=TRUE)
 }
}
dev.off()
}
#Percentage of run exposed to all fisheries by year

#Get # of fish exposed to any fishery

is_exposed<-0
total_exposed_all_fisheries<-rep(0,13)
for(y in 1:13){
  for(ind in 1:n_fish){
     for(f in 1:5){
        is_exposed <- is_exposed + exposure[ind,f,y]
     }
    if(is_exposed>0){
      total_exposed_all_fisheries[y]<- total_exposed_all_fisheries[y] + 1
     }
  }
}
###This results in 100% of the population being exposed to at least one fishery each year, so not bothering with a plot.


#--------Plots of probability of exposure with error bars

#Percentage of run exposed by fishery

x<-1:13
pdf(file=paste0("Population Percent Exposure by ",data_source," Fishery - Line plots w Error bars.pdf"))
#par(mfrow=c(1,1),mar=c(3,3,1,1), oma=c(5,5,3,1))
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1), oma=c(1,1,1,1))

if(data_source=="Commercial"){
  plot(mean_perc_exposed[1,], main="Area B", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(mean_perc_exposed[1,])
  arrows(x, mean_perc_exposed[1,]-sd_perc_exposed[1,], x, mean_perc_exposed[1,]+sd_perc_exposed[1,], length=0.05, angle=90, code=3, col="red")
plot(mean_perc_exposed[2,], main="Area D", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(mean_perc_exposed[2,])
  arrows(x, mean_perc_exposed[2,]-sd_perc_exposed[2,], x, mean_perc_exposed[2,]+sd_perc_exposed[2,], length=0.05, angle=90, code=3, col="red")
plot(mean_perc_exposed[3,], main="Area E", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(mean_perc_exposed[3,])
  arrows(x, mean_perc_exposed[3,]-sd_perc_exposed[3,], x, mean_perc_exposed[3,]+sd_perc_exposed[3,], length=0.05, angle=90, code=3, col="red")
plot(mean_perc_exposed[4,], main="Area G", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(mean_perc_exposed[4,])
  arrows(x, mean_perc_exposed[4,]-sd_perc_exposed[4,], x, mean_perc_exposed[4,]+sd_perc_exposed[4,], length=0.05, angle=90, code=3, col="red")
plot(mean_perc_exposed[5,], main="Area H", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(mean_perc_exposed[5,])
  arrows(x, mean_perc_exposed[5,]-sd_perc_exposed[5,], x, mean_perc_exposed[5,]+sd_perc_exposed[5,], length=0.05, angle=90, code=3, col="red")

}else if(data_source=="FN"){
  plot(mean_perc_exposed[1,], main="APM BSn", xlab="Year",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    points(mean_perc_exposed[1,])
    arrows(x, mean_perc_exposed[1,]-sd_perc_exposed[1,], x, mean_perc_exposed[1,]+sd_perc_exposed[1,], length=0.05, angle=90, code=3, col="red")

    y1<-array(as.numeric(NA),dim=c(1,13))
    for(i in 1:13){
      y1[i]<-max(0,mean_perc_exposed[2,i]-sd_perc_exposed[2,i])
    }
  plot(mean_perc_exposed[2,], main="APM GN", xlab="Year",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    points(mean_perc_exposed[2,])
    arrows(x, y1, x, mean_perc_exposed[2,]+sd_perc_exposed[2,], length=0.05, angle=90, code=3, col="red")
  plot(mean_perc_exposed[3,], main="BPM DN", xlab="Year",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    points(mean_perc_exposed[3,])
    arrows(x, mean_perc_exposed[3,]-sd_perc_exposed[3,], x, mean_perc_exposed[3,]+sd_perc_exposed[3,], length=0.05, angle=90, code=3, col="red")
  #Turn on if splitting APM DN and SN exposure (remember to rename mains for plot 2 and 3):
    #plot(mean_perc_exposed[4,], main="BPM DN", xlab="Year",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    #axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    #points(mean_perc_exposed[4,])
    #arrows(x, mean_perc_exposed[4,]-sd_perc_exposed[4,], x, mean_perc_exposed[4,]+sd_perc_exposed[4,], length=0.05, angle=90, code=3)
}else{ #data_source=="REC"
  
}
#Turning off for now, can turn on if you want to put more than one graph per page
#mtext(text="Year",side=1,line=1,outer=TRUE)
#mtext(text="% Exposed",side=2,line=1,outer=TRUE)
#mtext(text="Population Percent Exposure by Fishery",side=3,line=1,outer=TRUE)
dev.off()