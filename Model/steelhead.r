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

data_source<-"Commercial" #Options: "Commercial", "FN", "REC"

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

for(i in 101:(n_reps+100)){
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
passage_date<-(pmax(30,pmin(140,rnorm(fish,m_vec[i-100],s_vec[i-100]))))
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
         exposure[ind,f,y,i-100]<-exposure[ind,f,y,i-100]+fishery_array[loc+1,round(time_at_loc)+1,f,y]
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
       exposure[ind,f,y,i-100]<-exposure[ind,f,y,i-100]+fishery_array[loc+1,round(time_at_loc)+1,f,y]
      }
     }
   }
}
yr=yr+1
}
  progress(i,max.value=n_reps+100,init=101)
  Sys.sleep(0.01)
  if (i==n_reps+100) cat("Done!\n")
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

####NOTE: The code below the line only needs to be run if you re-run the model. The combined iterations are saved, 
####      and can be extracted with:

total_reps<-200

if(data_source=="Commercial"){
  exposure<-array(as.numeric(NA),dim=c(n_fish,5,13,total_reps))
  exposure<-readRDS("Com_exposure_1-200.RData")
}else if(data_source=="FN"){
  exposure<-array(as.numeric(NA),dim=c(n_fish,4,13,total_reps))
  exposure<-readRDS("EO_exposure_1-200.RData")
}

#This section adds the commercial and FN EO together:

total_reps<-200

  exposure_com<-array(as.numeric(NA),dim=c(n_fish,5,13,total_reps))
  exposure_com<-readRDS("Com_exposure_1-200.RData")

  exposure_FN<-array(as.numeric(NA),dim=c(n_fish,4,13,total_reps))
  exposure_FN<-readRDS("EO_exposure_1-200.RData")

  exposure<-array(as.numeric(NA),dim=c(n_fish,9,13,total_reps))
  
  for(i in 1:5){
    exposure[,i,,]<-exposure_com[,i,,]
  }
  
  for(i in 6:9){
    exposure[,i,,]<-exposure_FN[,i-5,,]
  }

#__________________________________________

if(data_source=="Commercial"){
  
temp_exposure1<-readRDS("Com_exposure_1-100.RData")
temp_exposure2<-readRDS("Com_exposure_101-200.RData")

exposure<-array(as.numeric(NA),dim=c(n_fish,5,13,total_reps))

#There was a mistake made when re-opening the file, where exposure didn't get reset to 4 fisheries, so it had a 5th fishery
#still saved. The code below gets rid of that fifth set of data

for(i in 1:100){
 exposure[,,,i]<-temp_exposure1[,,,i]
}

for(i in 101:200){
  exposure[,,,i]<-temp_exposure2[,,,i-100]
}

}else if(data_source=="FN"){
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

#-----Get total and mean exposure for all 8 fisheries together

n_fisheries<-8
n_reps<-200

total_exposed<-array(as.numeric(NA),dim=c(n_fisheries,13,n_reps))

for(i in 1:n_reps){
  for(y in 1:13){
    total_exposed[1,y,i]<-sum(exposure[,1,y,i]>0)
    total_exposed[2,y,i]<-sum(exposure[,2,y,i]>0)
    total_exposed[3,y,i]<-sum(exposure[,3,y,i]>0)
    total_exposed[4,y,i]<-sum(exposure[,4,y,i]>0)
    total_exposed[5,y,i]<-sum(exposure[,5,y,i]>0)
    total_exposed[6,y,i]<-sum(exposure[,6,y,i]>0)
    total_exposed[7,y,i]<-sum((exposure[,7,y,i]+exposure[,8,y,i])>0)
    total_exposed[8,y,i]<-sum(exposure[,9,y,i]>0)
  }
}

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
}else if(data_source=="AllCom"){  #For adding Com and FN together
  
  exposure_temp<-array(as.numeric(NA),dim=c(n_fish,n_fisheries,13,n_reps))
  
  exposure_temp[,1,,]<-exposure[,1,,]
  exposure_temp[,2,,]<-exposure[,2,,]
  exposure_temp[,3,,]<-exposure[,3,,]
  exposure_temp[,4,,]<-exposure[,4,,]
  exposure_temp[,5,,]<-exposure[,5,,]
  exposure_temp[,6,,]<-exposure[,6,,]
  exposure_temp[,7,,]<-exposure[,7,,]+exposure[,8,,]
  exposure_temp[,8,,]<-exposure[,9,,]
  
  cml_exposure<-array(as.numeric(NA),dim=c(n_fish,13,n_reps))
  
  for(i in 1:n_reps){
    for(y in 1:13){
      for(n in 1:n_fish){
        cml_exposure[n,y,i]<-sum(exposure_temp[n,,y,i]>0)
      }
    }
  }
  rm(exposure_temp)
  
}else{ #data_source=="REC"
  
}

if(data_source=="Commercial"){ #Or if "AllCom"
cml_exp_iters<-array(as.numeric(NA),dim=c(n_fisheries+1,13,total_reps))
mean_cml_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,13))
sd_cml_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,13))
mean_cml_perc_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,13))
sd_cml_perc_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,13))

for(f in 0:n_fisheries){
  for(y in 1:13){
   for(i in 1:total_reps){
    cml_exp_iters[f+1,y,i]<-sum(cml_exposure[,y,i]==f)
   }
  }
}

for(y in 1:13){
   for(f in 1:(n_fisheries+1)){
    mean_cml_exp[f,y]<-mean(cml_exp_iters[f,y,])
    mean_cml_perc_exp[f,y]<-mean_cml_exp[f,y]/1000*100
    sd_cml_exp[f,y]<-sd(cml_exp_iters[f,y,])
    sd_cml_perc_exp[f,y]<-sd_cml_exp[f,y]/1000*100
  }
}
}else if(data_source=="FN"){
  cml_exp_iters<-array(as.numeric(NA),dim=c(n_fisheries,13,total_reps))
  mean_cml_exp<-array(as.numeric(NA),dim=c(n_fisheries,13))
  sd_cml_exp<-array(as.numeric(NA),dim=c(n_fisheries,13))
  mean_cml_perc_exp<-array(as.numeric(NA),dim=c(n_fisheries,13))
  sd_cml_perc_exp<-array(as.numeric(NA),dim=c(n_fisheries,13))
  
  for(f in 0:(n_fisheries-1)){
    for(y in 1:13){
      for(i in 1:total_reps){
        cml_exp_iters[f+1,y,i]<-sum(cml_exposure[,y,i]==f)
      }
    }
  }
  
  for(y in 1:13){
    for(f in 1:n_fisheries){
      mean_cml_exp[f,y]<-mean(cml_exp_iters[f,y,])
      mean_cml_perc_exp[f,y]<-mean_cml_exp[f,y]/1000*100
      sd_cml_exp[f,y]<-sd(cml_exp_iters[f,y,])
      sd_cml_perc_exp[f,y]<-sd_cml_exp[f,y]/1000*100
    }
  }
}
rm(exposure_temp) #takes up a lot of space, might as well remove

#------------Get incremental exposure to fisheries

if(data_source=="Commercial"){

#Order of commercial fisheries when data imported: B D E G H
  
incr_exp<-array(0,dim=c(n_fisheries,13,total_reps))
mean_incr_exp<-array(0,dim=c(n_fisheries,13))
mean_incr_perc_exp<-array(0,dim=c(n_fisheries,13))
sd_incr_exp<-array(0,dim=c(n_fisheries,13))
sd_incr_perc_exp<-array(0,dim=c(n_fisheries,13))

####NOTE NEW ORDER OF FISHERIES!!! 1=G, 2=B, 3=D, 4=H, 5=E

#How many fish exposed to Area G fishery?
for(i in 1:total_reps){
  for(y in 1:13){
    incr_exp[1,y,i]<-sum(exposure[,4,y,i]>0)
  }
}

#How many fish exposed to Area B fishery that were not also exposed to Area G?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure[n,1,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[2,y,i]<-incr_exp[2,y,i]+1
      }else{
        incr_exp[2,y,i]<-incr_exp[2,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[2,y,i]<-incr_exp[2,y,i]+incr_exp[1,y,i]  remove for now, don't want the cumulative incremental exposure
  }
}

#How many fish exposed to Area D fishery that were not also exposed to Area G and B?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
          if((sum(exposure[n,2,y,i]>0)>sum(exposure[n,1,y,i]>0)) && (sum(exposure[n,2,y,i]>0)>sum(exposure[n,4,y,i]>0))){
              incr_exp[3,y,i]<-incr_exp[3,y,i]+1
          }else{
            incr_exp[3,y,i]<-incr_exp[3,y,i]
          }
    }
    #Add new fish to old total
    #incr_exp[3,y,i]<-incr_exp[3,y,i]+incr_exp[2,y,i]
  }
}


#How many fish exposed to Area H fishery that were not also exposed to Area G, B, and D?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
        if(sum(exposure[n,5,y,i]>0)>sum(exposure[n,2,y,i]>0) && sum(exposure[n,5,y,i]>0)>sum(exposure[n,1,y,i]>0) && sum(exposure[n,5,y,i]>0)>sum(exposure[n,4,y,i]>0)){
          incr_exp[4,y,i]<-incr_exp[4,y,i]+1
        }else{
          incr_exp[4,y,i]<-incr_exp[4,y,i]
        }
    }
    #Add new fish to old total
    #incr_exp[4,y,i]<-incr_exp[4,y,i]+incr_exp[3,y,i]
  }
}

#How many fish exposed to Area E fishery that were not also exposed to Area G, B, D, and H?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure[n,3,y,i]>0)>sum(exposure[n,5,y,i]>0) && sum(exposure[n,3,y,i]>0)>sum(exposure[n,2,y,i]>0) && sum(exposure[n,3,y,i]>0)>sum(exposure[n,1,y,i]>0) && sum(exposure[n,3,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[5,y,i]<-incr_exp[5,y,i]+1
      }else{
        incr_exp[5,y,i]<-incr_exp[5,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[5,y,i]<-incr_exp[5,y,i]+incr_exp[4,y,i]
  }
}

for(y in 1:13){
  for(f in 1:n_fisheries){
    mean_incr_exp[f,y]<-mean(incr_exp[f,y,])
    mean_incr_perc_exp[f,y]<-mean_incr_exp[f,y]/1000*100
    sd_incr_exp[f,y]<-sd(incr_exp[f,y,])
    sd_incr_perc_exp[f,y]<-sd_incr_exp[f,y]/1000*100
  }
}

}else if(data_source=="FN"){
  #Order of FN fisheries: BPM DN, APM GN, APM BSn  (needs refining)
  
  exposure_temp<-array(as.numeric(NA),dim=c(n_fish,n_fisheries-1,13,total_reps))
  
  exposure_temp[,1,,]<-exposure[,1,,]
  exposure_temp[,2,,]<-exposure[,2,,]+exposure[,3,,]
  exposure_temp[,3,,]<-exposure[,4,,]
  
  incr_exp<-array(0,dim=c(n_fisheries-1,13,total_reps))
  mean_incr_exp<-array(0,dim=c(n_fisheries-1,13))
  mean_incr_perc_exp<-array(0,dim=c(n_fisheries-1,13))
  sd_incr_exp<-array(0,dim=c(n_fisheries-1,13))
  sd_incr_perc_exp<-array(0,dim=c(n_fisheries-1,13))
  
  ####NOTE NEW ORDER OF FISHERIES!!! 1=BPM GN, 2=APM GN, 3=APM BSn
  
  #How many fish exposed to BPM GN fishery?
  for(i in 1:total_reps){
    for(y in 1:13){
      incr_exp[1,y,i]<-sum(exposure_temp[,3,y,i]>0)
    }
  }
  
  #How many fish exposed to APM GN fishery that were not also exposed to BPM GN?
  for(i in 1:total_reps){
    for(y in 1:13){
      for(n in 1:n_fish){
        if(sum(exposure_temp[n,2,y,i]>0)>sum(exposure_temp[n,3,y,i]>0)){
          incr_exp[2,y,i]<-incr_exp[2,y,i]+1
        }else{
          incr_exp[2,y,i]<-incr_exp[2,y,i]
        }
      }
      #Add new fish to old total
      #incr_exp[2,y,i]<-incr_exp[2,y,i]+incr_exp[1,y,i]
    }
  }
  
  #How many fish exposed to APM BSn fishery that were not also exposed to BPM GN and APM GN?
  for(i in 1:total_reps){
    for(y in 1:13){
      for(n in 1:n_fish){
        if((sum(exposure_temp[n,1,y,i]>0)>sum(exposure_temp[n,2,y,i]>0)) && (sum(exposure_temp[n,1,y,i]>0)>sum(exposure_temp[n,3,y,i]>0))){
          incr_exp[3,y,i]<-incr_exp[3,y,i]+1
        }else{
          incr_exp[3,y,i]<-incr_exp[3,y,i]
        }
      }
      #Add new fish to old total
      #incr_exp[3,y,i]<-incr_exp[3,y,i]+incr_exp[2,y,i]
    }
  }

  for(y in 1:13){
    for(f in 1:n_fisheries-1){
      mean_incr_exp[f,y]<-mean(incr_exp[f,y,])
      mean_incr_perc_exp[f,y]<-mean_incr_exp[f,y]/1000*100
      sd_incr_exp[f,y]<-sd(incr_exp[f,y,])
      sd_incr_perc_exp[f,y]<-sd_incr_exp[f,y]/1000*100
    }
  }
  
}else{ #data_source=="REC"
  
}



#-----------Incremental exposure of all 8 fisheries together--------------------------------------------------


#Order of commercial fisheries when data imported: B D E G H

incr_exp<-array(0,dim=c(n_fisheries,13,total_reps))
mean_incr_exp<-array(0,dim=c(n_fisheries,13))
mean_incr_perc_exp<-array(0,dim=c(n_fisheries,13))
sd_incr_exp<-array(0,dim=c(n_fisheries,13))
sd_incr_perc_exp<-array(0,dim=c(n_fisheries,13))

####NOTE NEW ORDER OF FISHERIES!!! 1=G, 2=B, 3=D, 4=H, 5=E

#How many fish exposed to Area G fishery?
for(i in 1:total_reps){
  for(y in 1:13){
    incr_exp[1,y,i]<-sum(exposure[,4,y,i]>0)
  }
}

#How many fish exposed to Area B fishery that were not also exposed to Area G?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure[n,1,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[2,y,i]<-incr_exp[2,y,i]+1
      }else{
        incr_exp[2,y,i]<-incr_exp[2,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[2,y,i]<-incr_exp[2,y,i]+incr_exp[1,y,i]  remove for now, don't want the cumulative incremental exposure
  }
}

#How many fish exposed to Area D fishery that were not also exposed to Area G and B?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if((sum(exposure[n,2,y,i]>0)>sum(exposure[n,1,y,i]>0)) && (sum(exposure[n,2,y,i]>0)>sum(exposure[n,4,y,i]>0))){
        incr_exp[3,y,i]<-incr_exp[3,y,i]+1
      }else{
        incr_exp[3,y,i]<-incr_exp[3,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[3,y,i]<-incr_exp[3,y,i]+incr_exp[2,y,i]
  }
}


#How many fish exposed to Area H fishery that were not also exposed to Area G, B, and D?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure[n,5,y,i]>0)>sum(exposure[n,2,y,i]>0) && sum(exposure[n,5,y,i]>0)>sum(exposure[n,1,y,i]>0) && sum(exposure[n,5,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[4,y,i]<-incr_exp[4,y,i]+1
      }else{
        incr_exp[4,y,i]<-incr_exp[4,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[4,y,i]<-incr_exp[4,y,i]+incr_exp[3,y,i]
  }
}

#How many fish exposed to Area E fishery that were not also exposed to Area G, B, D, and H?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure[n,3,y,i]>0)>sum(exposure[n,5,y,i]>0) && sum(exposure[n,3,y,i]>0)>sum(exposure[n,2,y,i]>0) && sum(exposure[n,3,y,i]>0)>sum(exposure[n,1,y,i]>0) && sum(exposure[n,3,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[5,y,i]<-incr_exp[5,y,i]+1
      }else{
        incr_exp[5,y,i]<-incr_exp[5,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[5,y,i]<-incr_exp[5,y,i]+incr_exp[4,y,i]
  }
}

####NOTE NEW ORDER OF FISHERIES!!! incr_exp[i,,] where i= 6=BPM GN, 7=APM GN, 8=APM BSn

#How many fish exposed to BPM GN fishery that were not also exposed to 5 commercial fisheries?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure_temp[n,8,y,i]>0)>sum(exposure[n,3,y,i]>0) && sum(exposure_temp[n,8,y,i]>0)>sum(exposure[n,5,y,i]>0) && sum(exposure_temp[n,8,y,i]>0)>sum(exposure[n,2,y,i]>0) && sum(exposure_temp[n,8,y,i]>0)>sum(exposure[n,1,y,i]>0) && sum(exposure_temp[n,8,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[6,y,i]<-incr_exp[6,y,i]+1
      }else{
        incr_exp[6,y,i]<-incr_exp[6,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[6,y,i]<-incr_exp[6,y,i]+incr_exp[5,y,i]
  }
}

#How many fish exposed to APM GN fishery that were not also exposed to BPM GN and 5 commercial fisheries?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure_temp[n,7,y,i]>0)>sum(exposure_temp[n,8,y,i]>0) && sum(exposure_temp[n,7,y,i]>0)>sum(exposure[n,3,y,i]>0) && sum(exposure_temp[n,7,y,i]>0)>sum(exposure[n,5,y,i]>0) && sum(exposure_temp[n,7,y,i]>0)>sum(exposure[n,2,y,i]>0) && sum(exposure_temp[n,7,y,i]>0)>sum(exposure[n,1,y,i]>0) && sum(exposure_temp[n,7,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[7,y,i]<-incr_exp[7,y,i]+1
      }else{
        incr_exp[7,y,i]<-incr_exp[7,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[7,y,i]<-incr_exp[7,y,i]+incr_exp[6,y,i]
  }
}

#How many fish exposed to APM BSn fishery that were not also exposed to BPM GN, APM GN, and 5 commercial fisheries?
for(i in 1:total_reps){
  for(y in 1:13){
    for(n in 1:n_fish){
      if(sum(exposure_temp[n,6,y,i]>0)>sum(exposure_temp[n,7,y,i]>0) && sum(exposure_temp[n,6,y,i]>0)>sum(exposure_temp[n,8,y,i]>0) && sum(exposure_temp[n,6,y,i]>0)>sum(exposure[n,3,y,i]>0) && sum(exposure_temp[n,6,y,i]>0)>sum(exposure[n,5,y,i]>0) && sum(exposure_temp[n,6,y,i]>0)>sum(exposure[n,2,y,i]>0) && sum(exposure_temp[n,6,y,i]>0)>sum(exposure[n,1,y,i]>0) && sum(exposure_temp[n,6,y,i]>0)>sum(exposure[n,4,y,i]>0)){
        incr_exp[8,y,i]<-incr_exp[8,y,i]+1
      }else{
        incr_exp[8,y,i]<-incr_exp[8,y,i]
      }
    }
    #Add new fish to old total
    #incr_exp[8,y,i]<-incr_exp[8,y,i]+incr_exp[7,y,i]
  }
}

for(y in 1:13){
  for(f in 1:n_fisheries){
    mean_incr_exp[f,y]<-mean(incr_exp[f,y,])
    mean_incr_perc_exp[f,y]<-mean_incr_exp[f,y]/1000*100
    sd_incr_exp[f,y]<-sd(incr_exp[f,y,])
    sd_incr_perc_exp[f,y]<-sd_incr_exp[f,y]/1000*100
  }
}


##############################
#Print plots to pdf file
##############################

setwd(plots_dir)

#Run timing plots
png(file="annual_runtiming.png")

yr=1995
y1<-array(as.numeric(NA),dim=c(1,22))
for(y in 1:22){
  y1[y]<-sh_runtiming$rt_mean[sh_runtiming$year==yr]-sh_runtiming$rt_sd[sh_runtiming$year==yr]
  yr=yr+1
}

yr=1995
y2<-array(as.numeric(NA),dim=c(1,22))
for(y in 1:22){
  y2[y]<-sh_runtiming$rt_mean[sh_runtiming$year==yr]+sh_runtiming$rt_sd[sh_runtiming$year==yr]
  yr=yr+1
}

#Plot of mean with sd around mean (not 50% date plus spread of run)
yr=1995
y1<-array(as.numeric(NA),dim=c(1,22))
for(y in 1:22){
  y1[y]<-sh_runtiming$rt_mean[sh_runtiming$year==yr]-sh_runtiming$rt_mean_sd[sh_runtiming$year==yr]
  yr=yr+1
}

yr=1995
y2<-array(as.numeric(NA),dim=c(1,22))
for(y in 1:22){
  y2[y]<-sh_runtiming$rt_mean[sh_runtiming$year==yr]+sh_runtiming$rt_mean_sd[sh_runtiming$year==yr]
  yr=yr+1
}


#Estimated mean annual run timing at Albion

par(mar=c(6,5,1,1))
x<-1:22
plot(sh_runtiming$rt_mean, main="", xlab="Year",ylab="Day of Year", xaxt="n", yaxt="n", type="l", bty="n", lty=1, ylim=range(240,340),cex.axis=1.5, cex.lab=1.5, mgp=c(4,1,0))
  axis(1, at=1:22,labels=c("1995","","1997","","1999","","2001","","2003","","2005","","2007","","2009","","2011","","2013","","2015",""), las=2, cex.axis=1.5)
  axis(2, at=c(240,250,260,270,280,290,300,310,320,330,340), labels=c("240","","260","","280","","300","","320","","340"), las=2, cex.axis=1.5)
  arrows(x, y1, x, y2, length=0.05, angle=90, code=3, col="red", lwd=2)
  points(sh_runtiming$rt_mean,cex=2, pch=16)

dev.off()

png(file="meanvssd_runtiming.png")
#mean vs sd to show no pattern between 50% date and spread of return
plot(sh_runtiming$rt_mean,sh_runtiming$rt_sd, xlab="Mean",ylab="Standard deviation")

dev.off()

#Percentage of run exposed by year

if(data_source=="Commercial"){
pdf(file="Population Percent Exposure to Commercial Fisheries by Year - Barplots.pdf")
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,3.1), oma=c(1,1,1,1), xpd=FALSE)

yr=2004
for(y in 1:13){
  
  y1<-array(as.numeric(NA),dim=c(1,5))
  for(f in 1:5){
    y1[f]<-max(0,mean_perc_exposed[f,y]-sd_perc_exposed[f,y])
  }
  
  y2<-array(as.numeric(NA),dim=c(1,5))
  for(f in 1:5){
    y2[f]<-min(100,mean_perc_exposed[f,y]+sd_perc_exposed[f,y])
  }
  
barCenters<-barplot(mean_perc_exposed[,y], main=paste0("Population Percent Exposure in ",yr), xlab="Fishery Area",ylab="% Exposure", 
        names.arg=c("B", "D", "E", "G", "H"), axis.lty=1, ylim=range(0,100))
  arrows(barCenters, y1, barCenters, y2, length=0.05, angle=90, code=3)

  yr=yr+1
 
#if(y==6 | y==12 | y==13){  
 # mtext(text="Fishery Area",side=1,line=1,outer=TRUE)
 #mtext(text="% Exposure",side=2,line=1,outer=TRUE)
 #mtext(text="Population Percent Exposure by Year",side=3,line=1,outer=TRUE)
 }

}else if(data_source=="FN"){
  pdf(file="Population Percent Exposure to FN Fisheries by Year - Barplots.pdf")
  par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,3.1), oma=c(1,1,1,1), xpd=FALSE)
  
  yr=2004
  for(y in 1:13){
    
    y1<-array(as.numeric(NA),dim=c(1,3))
    for(f in 1:3){
      y1[f]<-max(0,mean_perc_exposed[f,y]-sd_perc_exposed[f,y])
    }
    
    y2<-array(as.numeric(NA),dim=c(1,3))
    for(f in 1:3){
      y2[f]<-min(100,mean_perc_exposed[f,y]+sd_perc_exposed[f,y])
    }
    
    barCenters<-barplot(mean_perc_exposed[,y], main=paste0("Population Percent Exposure in ",yr), xlab="Fishery",ylab="% Exposure", 
                        names.arg=c("APM BSn","APM GN","BPM DN"), axis.lty=1, ylim=range(0,100))
    arrows(barCenters, y1, barCenters, y2, length=0.05, angle=90, code=3)
    
    yr=yr+1
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


#--------Plots of probability of exposure with error bars

#Percentage of run exposed by fishery

x<-1:13
pdf(file=paste0("Population Percent Exposure by ",data_source," Fishery - Line plots w Error bars.pdf"))
#par(mfrow=c(1,1),mar=c(3,3,1,1), oma=c(5,5,3,1))
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1), oma=c(1,1,1,1))

if(data_source=="Commercial"){
  plot(mean_perc_exposed[1,], main="Area B", xlab="",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    points(mean_perc_exposed[1,])
    arrows(x, mean_perc_exposed[1,]-sd_perc_exposed[1,], x, mean_perc_exposed[1,]+sd_perc_exposed[1,], length=0.05, angle=90, code=3, col="red")
  plot(mean_perc_exposed[2,], main="Area D", xlab="",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    points(mean_perc_exposed[2,])
    arrows(x, mean_perc_exposed[2,]-sd_perc_exposed[2,], x, mean_perc_exposed[2,]+sd_perc_exposed[2,], length=0.05, angle=90, code=3, col="red")
  plot(mean_perc_exposed[3,], main="Area E", xlab="",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    points(mean_perc_exposed[3,])
    arrows(x, mean_perc_exposed[3,]-sd_perc_exposed[3,], x, mean_perc_exposed[3,]+sd_perc_exposed[3,], length=0.05, angle=90, code=3, col="red")
  plot(mean_perc_exposed[4,], main="Area G", xlab="",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    points(mean_perc_exposed[4,])
    arrows(x, mean_perc_exposed[4,]-sd_perc_exposed[4,], x, mean_perc_exposed[4,]+sd_perc_exposed[4,], length=0.05, angle=90, code=3, col="red")
  plot(mean_perc_exposed[5,], main="Area H", xlab="",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
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

#-----------Plots of cumulative exposure--------------------------

pdf(file=paste0("Population Cumulative Perc Exposure to ",data_source," Fisheries - Bar Plots.pdf"))
#par(mfrow=c(1,1),mar=c(3,3,1,1), oma=c(5,5,3,1))
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,4.1), oma=c(1,1,1,1))

if(data_source=="Commercial"){
#Plot of all years together
barplot(mean_cml_perc_exp, main="Average Cumulative Exposure \nto Commercial Fisheries by Year",
          xlab="Year", col=c("lightblue","red","green4","darkgoldenrod1","mediumorchid4","darkgrey"),
          ylab="% Exposure",
          legend = c(0,1,2,3,4,5),names.arg=c("04","05","06","07","08","09","10","11","12","13","14","15","16"),args.legend=list(
            x=18,y=100,bty="n"))

#Plots for each year separately

yr=2004
for(y in 1:13){
  
  y1<-array(as.numeric(NA),dim=c(1,6))
  for(f in 1:6){
    y1[f]<-max(0,mean_cml_perc_exp[f,y]-sd_cml_perc_exp[f,y])
  }
  
  y2<-array(as.numeric(NA),dim=c(1,6))
  for(f in 1:6){
    y2[f]<-min(100,mean_cml_perc_exp[f,y]+sd_cml_perc_exp[f,y])
  }
  
barCenters<-barplot(mean_cml_perc_exp[,y], main=paste0("Average Cumulative Exposure \nto Commercial Fisheries in ",yr),
        xlab="# of fisheries each fish exposed to", col=c("lightblue","red","green4","darkgoldenrod1","mediumorchid4","darkgrey"),
        ylab="% Exposure", ylim=range(0,100),names.arg=c(0,1,2,3,4,5))
  arrows(barCenters, y1, barCenters, y2, length=0.05, angle=90, code=3)

  yr=yr+1
}

}else if(data_source=="FN"){
  #Plot of all years together
  barplot(mean_cml_perc_exp, main="Average Cumulative Exposure \nto FN Fisheries by Year",
          xlab="Year", col=c("lightblue","red","green4","darkgoldenrod1"),
          ylab="% Exposure",
          legend = c(0,1,2,3),names.arg=c("04","05","06","07","08","09","10","11","12","13","14","15","16"),args.legend=list(
            x=18,y=100,bty="n"))
  
  #Plots for each year separately
  
  yr=2004
  for(y in 1:13){
    
    y1<-array(as.numeric(NA),dim=c(1,4))
      for(f in 1:4){
      y1[f]<-max(0,mean_cml_perc_exp[f,y]-sd_cml_perc_exp[f,y])
      }
    
    y2<-array(as.numeric(NA),dim=c(1,4))
      for(f in 1:4){
        y2[f]<-min(100,mean_cml_perc_exp[f,y]+sd_cml_perc_exp[f,y])
      }
    
    barCenters<-barplot(mean_cml_perc_exp[,y], main=paste0("Average Cumulative Exposure \nto FN Fisheries in ",yr),
                        xlab="# of fisheries each fish exposed to", col=c("lightblue","red","green4","darkgoldenrod1"),
                        ylab="% Exposure", ylim=range(0,100),names.arg=c(0,1,2,3))
    arrows(barCenters, y1, barCenters, y2, length=0.05, angle=90, code=3)
    
    yr=yr+1
  }
}else if(data_source=="AllCom"){
  
#col=c("lightblue","red","green4","darkgoldenrod1","mediumorchid4","darkgrey","lightpink","chartreuse","royalblue","black")  
#col=c("gray90","gray80","gray70","gray60","gray50","gray40","gray30","gray20","gray10","gray0")
  
  #Plot of all years together
  barplot(mean_cml_perc_exp, main="Average Cumulative Exposure \nto Commercial Fisheries by Year",
          xlab="Year", col=c("red","orange","yellow","yellowgreen","green","darkcyan","blue","purple","black"),
          ylab="% Exposure",
          legend = c(0,1,2,3,4,5,6,7,8),names.arg=c("04","05","06","07","08","09","10","11","12","13","14","15","16"),args.legend=list(
            x=18,y=100,bty="n"))
  
  #Plots for each year separately
  #####Needs to be edited####
  
  yr=2004
  for(y in 1:13){
    
    y1<-array(as.numeric(NA),dim=c(1,6))
    for(f in 1:6){
      y1[f]<-max(0,mean_cml_perc_exp[f,y]-sd_cml_perc_exp[f,y])
    }
    
    y2<-array(as.numeric(NA),dim=c(1,6))
    for(f in 1:6){
      y2[f]<-min(100,mean_cml_perc_exp[f,y]+sd_cml_perc_exp[f,y])
    }
    
    barCenters<-barplot(mean_cml_perc_exp[,y], main=paste0("Average Cumulative Exposure \nto Commercial Fisheries in ",yr),
                        xlab="# of fisheries each fish exposed to", col=c("lightblue","red","green4","darkgoldenrod1","mediumorchid4","darkgrey"),
                        ylab="% Exposure", ylim=range(0,100),names.arg=c(0,1,2,3,4,5))
    arrows(barCenters, y1, barCenters, y2, length=0.05, angle=90, code=3)
    
    yr=yr+1
  }
  
}
dev.off()

#-------------Plots of incremental exposure------------------------

pdf(file=paste0("Population Incremental Exposure by ",data_source," Fishery - Line plots w Error bars.pdf"))
#par(mfrow=c(1,1),mar=c(3,3,1,1), oma=c(5,5,3,1))
#Default mar=c(5.1,4.1,4.1,2.1)
par(mfrow=c(2,1),mar=c(5.1,4.1,1,0.5), oma=c(1,1,1,1))

if(data_source=="Commercial"){
  
  y1<-array(as.numeric(NA),dim=c(5,13))
  for(y in 1:13){
    for(f in 1:n_fisheries){
      y1[f,y]<-max(0,mean_incr_perc_exp[f,y]-sd_incr_perc_exp[f,y])
    }
  }
  
  y2<-array(as.numeric(NA),dim=c(5,13))
  for(y in 1:13){
    for(f in 1:n_fisheries){
      y2[f,y]<-min(100,mean_incr_perc_exp[f,y]+sd_incr_perc_exp[f,y])
    }
  }
  
  yr=2004
  x<-1:5
  
  for(y in 1:13){ 
    plot(mean_incr_perc_exp[,y], main=paste0("Incremental Exposure to \nCommercial Fisheries in ",yr), xlab="Fishery",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:5,labels=c("Area G","Area B","Area D", "Area H", "Area E"), las=1)
    points(mean_incr_perc_exp[,y])
    arrows(x, y1[,y], x, y2[,y], length=0.05, angle=90, code=3, col="red")
    yr=yr+1
  }
  
}else if(data_source=="FN"){
  
  y1<-array(as.numeric(NA),dim=c(3,13))
  for(y in 1:13){
    for(f in 1:n_fisheries-1){
      y1[f,y]<-max(0,mean_incr_perc_exp[f,y]-sd_incr_perc_exp[f,y])
    }
  }
  
  y2<-array(as.numeric(NA),dim=c(3,13))
  for(y in 1:13){
    for(f in 1:n_fisheries-1){
      y2[f,y]<-min(100,mean_incr_perc_exp[f,y]+sd_incr_perc_exp[f,y])
    }
  }

  yr=2004
  x<-1:3
 for(y in 1:13){ 
  plot(mean_incr_perc_exp[,y], main=paste0("Incremental Exposure to FN Fisheries in ",yr), xlab="Fishery",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:3,labels=c("BPM DN","APM GN","APM BSn"), las=1)
    points(mean_incr_perc_exp[,y])
    arrows(x, y1[,y], x, y2[,y], length=0.05, angle=90, code=3, col="red")
    yr=yr+1
 }
    
}else if(data_source=="AllCom"){
  
  y1<-array(as.numeric(NA),dim=c(8,13))
  for(y in 1:13){
    for(f in 1:n_fisheries){
      y1[f,y]<-max(0,mean_incr_perc_exp[f,y]-sd_incr_perc_exp[f,y])
    }
  }
  
  y2<-array(as.numeric(NA),dim=c(8,13))
  for(y in 1:13){
    for(f in 1:n_fisheries){
      y2[f,y]<-min(100,mean_incr_perc_exp[f,y]+sd_incr_perc_exp[f,y])
    }
  }
  
  yr=2004
  x<-1:8
  
  for(y in 1:13){ 
    plot(mean_incr_perc_exp[,y], main=paste0("Incremental Exposure to Commercial Fisheries in ",yr), xlab="Fishery",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:8, cex.axis=0.8,labels=c("Area G","Area B","Area D", "Area H", "Area E", "BPM GN", "APM GN", "APM BSn"), las=1)
    arrows(x, y1[,y], x, y2[,y], length=0.05, angle=90, code=3, col="red")
    points(mean_incr_perc_exp[,y],pch=16)
    yr=yr+1
  }
  
}else{ #data_source=="REC"
  
}
dev.off()
