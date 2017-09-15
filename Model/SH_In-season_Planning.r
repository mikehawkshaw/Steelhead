###################################################################
#Steelhead Predictors of In-River Timing and Exposure model
#Modified for in-season estimation of exposure given fishing plans
#Authors: Brittany Jenewein, Mike Hawkshaw
###################################################################

library("xtable")
library("svMisc")

#Input fishery plan data
#Fishery Opening Times
#Can input any day between Sept. 15 to Nov. 15 as fishery opening. Must specify opening times (24 hr format)
#If you don't want to include all 10 openings, plug in default: "2017-09-15 00:00:00"

yr="2017"

OP1_S<-as.POSIXlt("2017-10-23 09:00:00")
OP1_E<-as.POSIXlt("2017-10-23 16:00:00")

OP2_S<-as.POSIXlt("2017-10-24 09:00:00")
OP2_E<-as.POSIXlt("2017-10-24 16:00:00")

OP3_S<-as.POSIXlt("2017-10-25 09:00:00")
OP3_E<-as.POSIXlt("2017-10-25 16:00:00")

OP4_S<-as.POSIXlt("2017-10-26 00:01:00")
OP4_E<-as.POSIXlt("2017-10-26 23:59:00")

OP5_S<-as.POSIXlt("2017-10-27 09:00:00")
OP5_E<-as.POSIXlt("2017-10-27 16:00:00")

OP6_S<-as.POSIXlt("2017-09-15 00:00:00")
OP6_E<-as.POSIXlt("2017-09-15 00:00:00")

OP7_S<-as.POSIXlt("2017-09-15 00:00:00")
OP7_E<-as.POSIXlt("2017-09-15 00:00:00")

OP8_S<-as.POSIXlt("2017-09-15 00:00:00")
OP8_E<-as.POSIXlt("2017-09-15 00:00:00")

OP9_S<-as.POSIXlt("2017-09-15 00:00:00")
OP9_E<-as.POSIXlt("2017-09-15 00:00:00")

OP10_S<-as.POSIXlt("2017-09-15 00:00:00")
OP10_E<-as.POSIXlt("2017-09-15 00:00:00")

#Fishing areas
#Default = "closed"
#BPM = Sandheads to Port Mann Bridge ; Use for MFN/TFN fisheries
#Area E = Sandheads to Mission; Use for Area E fisheries
#APM = Kanaka Creek to Hope (if Katzie/Kwantlen sign CFA, would need to adjust to include Port Mann to Kanaka Creek)

OP1_Area<-"BPM"
OP2_Area<-"Area E"
OP3_Area<-"BPM"
OP4_Area<-"APM"
OP5_Area<-"APM"
OP6_Area<-"closed"
OP7_Area<-"closed"
OP8_Area<-"closed"
OP9_Area<-"closed"
OP10_Area<-"closed"

opening_dates<-c(OP1_S,OP1_E,OP2_S,OP2_E,OP3_S,OP3_E,OP4_S,OP4_E,OP5_S,OP5_E,OP6_S,OP6_E,OP7_S,OP7_E,OP8_S,OP8_E,OP9_S,OP9_E,OP10_S,OP10_E)
#Contains fillers so the loop in function will work with opening dates
opening_areas<-c(OP1_Area,"blank",OP2_Area,"blank",OP3_Area,"blank",OP4_Area,"blank",OP5_Area,"blank",OP6_Area,"blank",OP7_Area,"blank",OP8_Area,"blank",OP9_Area,"blank",OP10_Area,"blank")

#This only needs the year value changed annually
season_start<-as.POSIXlt("2017-09-15 00:00:00")

###########################################
#Format fishery openings data
###########################################

fishery_array<-array(0, dim = c(158,1488)) #row (km), column (hr)
open<-array(0,dim=c(1,10))
close<-array(0,dim=c(1,10))
j<-1

for(i in seq(1,19,2)){
  open[j]<-opening_dates[i]-season_start
  close[j]<-opening_dates[i+1]-season_start
  
  if(open[j]!=0){ #If not including the opening, then skip
  
  #Convert date/time to hours
  open[j]<-as.numeric(floor(open[j]*24))
  close[j]<-as.numeric(ceiling(close[j]*24)-1) #Round up for cases like closing at 8:59, then subtract 1 so the close hour (e.g.9:00-10:00) is not included in matrix
  
  #Turn "on" fisheries in the fishery_array
  #Note: km start is never 0 because it messes up the arrays; added 1km to all km markers (see Lookup Tables)
    if(opening_areas[i]=="BPM"){
      km_start=1 #Sandheads
      km_end=37 #Port Mann Bridge
    }else if(opening_areas[i]=="APM"){
      km_start=56 #If need to include PM-Kanaka, =37; if start at Kanaka, = 56
      km_end=158 #Hope
    }else if(opening_areas[i]=="Area E"){
      km_start=1 #Sandheads
      km_end=80 #Mission
    }
  
    fishery_array[seq(km_start,km_end,1),seq(open[j],close[j],1)]<-1
  }
  j<-j+1
  }

###################################
#set up a fake steelhead population
#IBM like model
###################################
n_fish<-1000
n_reps<-10000
fish<-seq(1,n_fish,by=1)

#each fish has characteristics and they are in these vectors
exposure<-array(0,dim=c(n_fish,n_reps))
speeds<-rep(0,n_fish)
passage_date<-rep(0,n_fish)

################################################################################################
#Population characteristics (these are the hypotheses about the population that will be tested)
################################################################################################

#Start the clock
ptm <- proc.time()

for(i in 1:(n_reps)){
  set.seed(i)
  
    #Get day of year for Oct 15 of year of interest (season start)
    seasonstart_doy <- as.numeric(strftime(paste(yr,"-09-15",sep=""), format = "%j"))
    
    #Run timing based on Bayesian estimator (grand mean)
    rt_mean<-282.131681 
    rt_mean_sd<-12.416564
    rt_sd<-16.510714
    rt_sd_sd<-5.970674
    
    #cumulative and daily proportions of the run vulnerable to each fishery
    m_vec<-rnorm(n_reps,rt_mean,rt_mean_sd) 
    s_vec<-abs(rnorm(n_reps,rt_sd,rt_sd_sd))
    
    #passage_date = the date that the fish passes Albion
    #passage_date<-(pmax(0,pmin(46,rnorm(fish,m_vec[i],s_vec[i])-seasonstart_doy))) #subtract season start day to put in correct position in matrix
    passage_date<-rnorm(fish,m_vec[i],s_vec[i])-seasonstart_doy
    passage_hour<-passage_date*24 #convert to hours. Hour 0 = midnight July 15
    
    #Speed that the fish travel. Assumptions based on speed of chum in river.
    speed_mean<-20 #km/day
    speed_sd<-3
    speeds<-(pmax(9,pmin(55,rnorm(fish,speed_mean,speed_sd))))/24 #km/hr
    
    #################################################
    #Move fish BACKWARD from Albion through fisheries
    #################################################
      
      #Loop through each fish
        
        for(loc in 1:59){ #1 = Sandheads; 59 = Albion
          
          time_at_loc<-round(passage_hour-(59-loc)/speeds)
          
          #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
          #If the fish is in the area before the time we care about, then obviously it's not exposed to any fisheries. May want to
          #edit later so that we can run it longer. Will need to make the opening matrices larger.
          #if (time_at_loc>0 & time_at_loc<1488){
           # exposure[ind,i]<-exposure[ind,i]+fishery_array[loc,time_at_loc]
         # }
          exposure[,i]<-ifelse(time_at_loc>0 & time_at_loc<1488,exposure[,i]+fishery_array[loc,ifelse(time_at_loc>0 & time_at_loc<1488,time_at_loc,1)],exposure[,i])
        
      }
    
    #################################################
    #Move fish FORWARD from Albion through fisheries
    #################################################
    
      #Loop through each fish

        for(loc in 60:158){ #From Albion, not including Albion start, up to Hope
          
          time_at_loc<-round(passage_hour+(loc-59)/speeds)
          
          #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
          #If the fish is in the area after the time we care about, then obviously it's not exposed to any fisheries.
          #if (time_at_loc>0 & time_at_loc<1488){
           # exposure[ind,i]<-exposure[ind,i]+fishery_array[loc,time_at_loc]
          #}
          exposure[,i]<-ifelse(time_at_loc>0 & time_at_loc<1488,exposure[,i]+fishery_array[loc,ifelse(time_at_loc>0 & time_at_loc<1488,time_at_loc,1)],exposure[,i])
          
        }
      
  progress(i,max.value=n_reps,init=1)
  Sys.sleep(0.01)
  if (i==n_reps) cat("Done!\n")
}
#Stop the clock
proc.time() - ptm

#Save iterations - change file name as appropriate
saveRDS(exposure,file=paste0("com_plan_exposure_",Sys.Date(),".RData"))

#-------------Get #/% of fish exposed by fishing plan
total_exposed<-array(as.numeric(NA),dim=c(1,n_reps))
perc_exposed<-array(as.numeric(NA),dim=c(1,n_reps))

for(i in 1:n_reps){
      total_exposed[i]<-sum(exposure[,i]>0)
      perc_exposed[i]<-total_exposed[i]/1000*100
}

#-------------Probability of protecting 80% of the run:

print("Mean % exposed (over all iterations)")
print(mean(perc_exposed))
print("95%")
print(mean(perc_exposed)+1.96*sd(perc_exposed))
print("5%")
print(mean(perc_exposed)-1.96*sd(perc_exposed))
print("Probabilty of protecting 80% of the run")
prob_20_perc<-sum(perc_exposed<=20)/n_reps
print(prob_20_perc)
high_degree_of_confidence<-pnorm(0.2, mean(1-perc_exposed/100),sd(1-perc_exposed/100))
print(high_degree_of_confidence)

#generate a "predicted RT curve" based on the variability we see

#day of the year
x_pred<-seq(225,350)
y_pred_matrix<-matrix(NA,length(x_pred),n_reps)
cord_x<-array(0,dim=c(1,10))
cord_y<-array(0,dim=c(1,10))

for(i in 1:n_reps){
  y_pred_matrix[,i]<-dnorm(x_pred,m_vec[i],s_vec[i])
}

  par(mfcol=c(1,1))
  
  matplot(x_pred,y_pred_matrix[,1:50],type="l",col="light grey",  ylab="Proportion of run",
          xlab="Day of Year (250=Sept 7; p50=282/Oct 9)",xlim=c(225,350),lwd=0.2,bty="l")
  abline(v=rt_mean,col = "grey",lwd=2)
  lines(x_pred,dnorm(x_pred,rt_mean,rt_sd),lwd=2,lty=2, col="black")

  for(j in 1:10){
    if(open[j]!=0){
      cord_x <- c(floor(open[j]/24)+seasonstart_doy,seq(floor(open[j]/24)+seasonstart_doy,ceiling(close[j]/24)+seasonstart_doy,length=20),ceiling(close[j]/24)+seasonstart_doy)
      cord_y <- c(0,dnorm(seq(floor(open[j]/24)+seasonstart_doy,ceiling(close[j]/24)+seasonstart_doy,length=20),rt_mean,rt_sd),0)
      polygon(cord_x,cord_y,col=rgb(0, 0, 1,0.25))
      }
   }

plot(density(perc_exposed/100), main=paste("Predictive density of Steelhead exposed to fisheries \n(assume average 20km/day movement speed)"), bty="l", cex.main=0.85)
