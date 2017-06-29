#steelhead.r
#exposure model for steelhead

#read in data
source("directories.R")
library("xtable")
setwd(data_dir)

###########################################
#Read in fishery openings data and format
###########################################

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
for(y in 1:13){
#Get day of year for July 15 of year of interest (season start)
seasonstart_doy <- as.numeric(strftime(paste(yr,"-07-15",sep=""), format = "%j"))

################################################################################################
#Population characteristics (these are the hypothesis about the population that will be tested)
################################################################################################

#Run-timing of the population. Based on mean and SD calculated in file "badestimator.r"
rt_mean<-mean(sh_runtiming$mean)-seasonstart_doy #subtract season start day to put in correct position in matrix
rt_sd<-mean(sh_runtiming$sd)

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
for(y in 1:13){
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
for(y in 1:13){
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

#Get total exposure time by fishery
#Not sure this is a useful metric, but it's here in case we decide to use it
yr=2004

total_exposure<-array(as.numeric(NA),dim=c(5,13))

for(y in 1:13){
  for(f in 1:5){
    total_exposure[f,y]<-sum(exposure[,f,y])
  }
  yr=yr+1
}


#Get # of fish exposed by fishery
yr=2004

total_exposed<-array(as.numeric(NA),dim=c(5,13))

for(y in 1:13){
  for(f in 1:5){
   total_exposed[f,y]<-sum(exposure[,f,y]>0)
   }
yr=yr+1
}


#--------------Generate barplots---------------------

pdf(file="Population Exposure by Fishery - Barplots.pdf")
par(mfrow=c(2,3))
par(mar=c(2,2,2,2))

barplot(total_exposed[1,], main="Area B", xlab="Year",ylab="# Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,1000))
barplot(total_exposed[2,], main="Area D", xlab="Year",ylab="# Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,1000))
barplot(total_exposed[3,], main="Area E", xlab="Year",ylab="# Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,1000))
barplot(total_exposed[4,], main="Area G", xlab="Year",ylab="# Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,1000))
barplot(total_exposed[5,], main="Area H", xlab="Year",ylab="# Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,1000))

dev.off()

#Percentage of run exposed by fishery

pdf(file="Population Percent Exposure by Fishery - Barplots.pdf")
par(mfrow=c(2,3), oma=c(5,5,3,0),mar=c(2,2,2,2), xpd=FALSE)

barplot(total_exposed[1,]/1000*100, main="Area B", xlab="",ylab="", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,100))
barplot(total_exposed[2,]/1000*100, main="Area D", xlab="Year",ylab="% Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,100))
barplot(total_exposed[3,]/1000*100, main="Area E", xlab="Year",ylab="% Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,100))
barplot(total_exposed[4,]/1000*100, main="Area G", xlab="Year",ylab="% Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,100))
barplot(total_exposed[5,]/1000*100, main="Area H", xlab="Year",ylab="% Exposed", names.arg=seq(from=2004,to=2016, by=1), axis.lty=1,ylim=range(0,100))
  
mtext(text="Year",side=1,line=1,outer=TRUE)
mtext(text="% Exposed",side=2,line=1,outer=TRUE)
mtext(text="Population Percent Exposure by Fishery - Barplots",side=3,line=1,outer=TRUE)
dev.off()

#Percentage of run exposed by year
pdf(file="Population Percent Exposure by Year - Barplots.pdf")
par(mfrow=c(3,2), oma=c(5,5,3,0), mar=c(2,2,2,2), xpd=FALSE)

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

#-----------------------Generate line plots for Mike -------------------------------

pdf(file="Population Exposure by Fishery - Line plots.pdf")
par(mfrow=c(2,3),mar=c(3,3,1,1), oma=c(5,5,3,1))

plot(total_exposed[1,], main="Area B", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[1,])
plot(total_exposed[2,], main="Area D", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[2,])
plot(total_exposed[3,], main="Area E", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[3,])
plot(total_exposed[4,], main="Area G", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[4,])
plot(total_exposed[5,], main="Area H", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[5,])
  
mtext(text="Year",side=1,line=1,outer=TRUE)
mtext(text="# Exposed",side=2,line=1,outer=TRUE)
mtext(text="Population Exposure by Fishery",side=3,line=1,outer=TRUE)

dev.off()

#Percentage of run exposed by fishery

pdf(file="Population Percent Exposure by Fishery - Line plots.pdf")
par(mfrow=c(2,3),mar=c(3,3,1,1), oma=c(5,5,3,1))

plot(total_exposed[1,]/1000*100, main="Area B", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(total_exposed[1,]/1000*100)
plot(total_exposed[2,]/1000*100, main="Area D", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(total_exposed[2,]/1000*100)
plot(total_exposed[3,]/1000*100, main="Area E", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(total_exposed[3,]/1000*100)
plot(total_exposed[4,]/1000*100, main="Area G", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(total_exposed[4,]/1000*100)
plot(total_exposed[5,]/1000*100, main="Area H", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
  points(total_exposed[5,]/1000*100)

mtext(text="Year",side=1,line=1,outer=TRUE)
mtext(text="% Exposed",side=2,line=1,outer=TRUE)
mtext(text="Population Percent Exposure by Fishery",side=3,line=1,outer=TRUE)
dev.off()

pdf(file="Population Percent Exposure by Year (All Gear) - Line plots.pdf")
par(las=2)
par(mar=c(4,4,1,1), cex=0.75)
plot(total_exposed[1,]/1000*100, type="b",pch=21,bg=1,main="",lwd=1.1, xlab="",ylab="", xaxt="n", bty="n", lty=1, ylim=range(0,100))
lines(total_exposed[2,]/1000*100, type="b",pch=21,bg=2,lwd=1.1)
lines(total_exposed[3,]/1000*100, type="b",pch=21,bg=3,lwd=1.1)
lines(total_exposed[4,]/1000*100, type="b",pch=21,bg=4,lwd=1.1)
lines(total_exposed[5,]/1000*100, type="b",pch=21,bg=5,lwd=1.1)
axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
leg<-c("Area B","Area D","Area E","Area G","Area H")
legend(11.5,104,legend=leg, fill=c(1,2,3,4,5), cex=0.85)
dev.off()

#Percentage of run exposed by year
pdf(file="Population Percent Exposure by Year - Line plots.pdf")
par(mfrow=c(3,2), oma=c(5,5,3,1), mar=c(2,2,2,2), xpd=FALSE)

yr=2004
for(y in 1:13){
  plot(total_exposed[,y]/1000*100, main=yr, xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:5,labels=c("B", "D", "E", "G", "H"))
    points(total_exposed[,y]/1000*100)
  yr=yr+1
  if(y==6 | y==12 | y==13){  
    mtext(text="Fishery Area",side=1,line=1,outer=TRUE)
    mtext(text="% Exposure",side=2,line=1,outer=TRUE)
    mtext(text="Population Percent Exposure by Year",side=3,line=1,outer=TRUE)
  }
}
dev.off()