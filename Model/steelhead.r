#steelhead.r
#Steelhead Predictors of In-river Timing and Exposure models

#read in data
source("directories.R")
library("xtable")
library("svMisc")
setwd(data_dir)


##########################
#Set inputs
##########################
km_end<-624
n_fisheries<-24
n_hours<-3336
n_years<-13
n_fish<-1000
n_reps<-1000 #1000 reps with the full 24 fisheries array takes about 8 hrs
total_reps<-5000 #Total reps after running model multiple times

fishery_names<-c("Area B CM","Area B PK","Area B SK",
                 "Area D CM","Area D PK","Area D SK",
                 "Area E CM","Area E PK","Area E SK",
                 "Area G CM","Area G PK","Area G SK",
                 "Area H CM","Area H PK","Area H SK",
                 "BPM DN CM","BPM DN PK","BPM DN SK",
                 "APM GN CM","APM GN PK","APM GN SK",
                 "APM BSn CM","APM BSn PK","APM BSn SK")

###########################################
#Read in fishery openings data and format
###########################################

fishery_array<- array(as.numeric(NA), dim = c(625,3336,n_fisheries,n_years)) #row, column, fishery, year
openings<-array(NA,dim=c(n_fisheries,n_years))

yr<-2004

for(i in 1:n_years){
  
path = paste0("C:/DFO-MPO/github/Steelhead/Data/BySpecies/",yr,"/")

file.names <- dir(path, pattern =".csv")

for(k in 1:n_fisheries){
  fishery_array[,,k,i] <- as.matrix(read.csv(paste0(path,file.names[k]),colClasses=c("NULL",rep(NA,3336))))
  #Does this fishery have at least one opening?
    openings[k,i]<-sum(fishery_array[,,k,i])>0
}
yr=yr+1
}

colnames(fishery_array)<-NULL

colnames(openings)<-c("Fishery","Year","Openings","File_Name")
openings<-as.data.frame(openings)

saveRDS(openings,"Fisheries_with_openings.RData")

sh_runtiming<-as.data.frame(read.csv("steelhead_runtiming.csv", header=T))

###################################
#set up a fake steelhead population
#IBM like model
###################################

fish<-seq(1,n_fish,by=1)

#each fish has characteristics and they are in these vectors
exposure<-array(0,dim=c(n_fish,n_fisheries,n_years,n_reps))
speeds<-rep(0,n_fish)
passage_date<-rep(0,n_fish)

################################################################################################
#Population characteristics (these are the hypotheses about the population that will be tested)
################################################################################################

#Start the clock
ptm <- proc.time()

for(i in 4001:(n_reps+4000)){
set.seed(i)
  
yr=2004 #re-initialize year variable

#Loop through each year
for(y in 1:n_years){

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
passage_date<-(pmax(30,pmin(140,rnorm(fish,m_vec[i-4000],s_vec[i-4000]))))
passage_hour<-passage_date*24 #convert to hours. Hour 0 = midnight July 15

#Speed that the fish travel. Assumptions based on speed of other salmonids.
speed_SW_mean<-34.8 #km/day based on Pieter Van Will's chum tagging study 2000-2002 Area 12/13
speed_SW_sd<-4.1
speed_FW_mean<-20 #km/day applied to Fraser chum in river (need to find reference)
speed_FW_sd<-3
speeds_SW<-(pmax(9,pmin(55,rnorm(fish,speed_SW_mean,speed_SW_sd))))/24 #km/hr
speeds_FW<-(pmax(9,pmin(55,rnorm(fish,speed_FW_mean,speed_FW_sd))))/24 #km/hr

#################################################
#Move fish BACKWARD from Albion through fisheries
#################################################

for(loc in 1:494){ #494 is km where Albion located
       
  time_at_loc<-round(passage_hour-(494-loc)/speeds_SW)
  
  for(f in 1:n_fisheries){   
  #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
  #If the fish is in the area before the time we care about, then obviously it's not exposed to any fisheries. May want to
  #edit later so that we can run it longer. Will need to make the opening matrices larger.
  exposure[,f,y,i-4000]<-ifelse(time_at_loc>0 & time_at_loc<3335,
                          exposure[,f,y,i-4000]+fishery_array[loc+1,
                                            ifelse(time_at_loc>0 & time_at_loc<3335,time_at_loc+1,1),f,y],exposure[,f,y,i-4000])
  }
}

#################################################
#Move fish FORWARD from Albion through fisheries
#################################################

#Assuming switch from SW to FW speeds happens at Albion, otherwise model gets too complicated. 
#Probably need to explore impacts.

for(loc in 495:km_end){ #From Albion, not including Albion start
       
  time_at_loc<-round(passage_hour+(loc-494)/speeds_FW)
  
  for(f in 1:n_fisheries){
  #check exposure against fishery matrix - sum the number of times each fish passes through an area during an open fishery
  #If the fish is in the area after the time we care about, then obviously it's not exposed to any fisheries.
  exposure[,f,y,i-4000]<-ifelse(time_at_loc>0 & time_at_loc<3335,
                          exposure[,f,y,i-4000]+fishery_array[loc+1,
                                            ifelse(time_at_loc>0 & time_at_loc<3335,time_at_loc+1,1),f,y],exposure[,f,y,i-4000])
     }
  }
yr=yr+1
}
  progress(i,max.value=n_reps+4000,init=4001)
  Sys.sleep(0.01)
  if (i==n_reps+4000) cat("Done!\n")
}
#Stop the clock
proc.time() - ptm

#Save iterations - change file name as appropriate
saveRDS(exposure,file="ComEO_exposure_4001-5000.RData")

rm(list=ls()) #Reset global environment
gc() #garbage collector - releases memory back to computer

##############################
#Manipulate exposure data
##############################

#Add multiple exposure runs together:
#This is not very dynamic but it is fine for now...
####NOTE: The code below only needs to be run if you re-run the model. 

# Initialize array for # of fish exposed by fishery
total_exposed<-array(as.numeric(NA),dim=c(n_fisheries,n_years,total_reps))

#Initialize array for cumulative exposure to fisheries
cml_exposure<-array(as.numeric(NA),dim=c(n_fish,n_years,total_reps))

exposure_temp<-readRDS("ComEO_exposure_1-1000.RData")

for(i in 1:1000){
 for(y in 1:n_years){
  for(f in 1:n_fisheries){
   total_exposed[f,y,i]<-sum(exposure_temp[,f,y,i]>0)
   }
  }
}

for(i in 1:1000){
  for(y in 1:n_years){
    for(n in 1:n_fish){
      cml_exposure[n,y,i]<-sum(exposure_temp[n,,y,i]>0)
    }
  }
}

exposure_temp<-readRDS("ComEO_exposure_1001-2000.RData")

for(i in 1001:2000){
  for(y in 1:n_years){
    for(f in 1:n_fisheries){
      total_exposed[f,y,i]<-sum(exposure_temp[,f,y,i-1000]>0)
    }
  }
}

for(i in 1001:2000){
  for(y in 1:n_years){
    for(n in 1:n_fish){
      cml_exposure[n,y,i]<-sum(exposure_temp[n,,y,i-1000]>0)
    }
  }
}

exposure_temp<-readRDS("ComEO_exposure_2001-3000.RData")

for(i in 2001:3000){
  for(y in 1:n_years){
    for(f in 1:n_fisheries){
      total_exposed[f,y,i]<-sum(exposure_temp[,f,y,i-2000]>0)
    }
  }
}

for(i in 2001:3000){
  for(y in 1:n_years){
    for(n in 1:n_fish){
      cml_exposure[n,y,i]<-sum(exposure_temp[n,,y,i-2000]>0)
    }
  }
}

exposure_temp<-readRDS("ComEO_exposure_3001-4000.RData")

for(i in 3001:4000){
  for(y in 1:n_years){
    for(f in 1:n_fisheries){
      total_exposed[f,y,i]<-sum(exposure_temp[,f,y,i-3000]>0)
    }
  }
}

for(i in 3001:4000){
  for(y in 1:n_years){
    for(n in 1:n_fish){
      cml_exposure[n,y,i]<-sum(exposure_temp[n,,y,i-3000]>0)
    }
  }
}

exposure_temp<-readRDS("ComEO_exposure_4001-5000.RData")

for(i in 4001:5000){
  for(y in 1:n_years){
    for(f in 1:n_fisheries){
      total_exposed[f,y,i]<-sum(exposure_temp[,f,y,i-4000]>0)
    }
  }
}

for(i in 4001:5000){
  for(y in 1:n_years){
    for(n in 1:n_fish){
      cml_exposure[n,y,i]<-sum(exposure_temp[n,,y,i-4000]>0)
    }
  }
}

rm(exposure_temp) #Takes up a lot of memory, best to remove it when done
gc()

saveRDS(total_exposed,"ComEO_tot_exposure_1-5000.RData")
saveRDS(cml_exposure,"ComEO_cml_exposure_1-5000.RData")

#The combined iterations are saved, and can be extracted with:

total_exposed<-array(as.numeric(NA),dim=c(n_fisheries,n_years,total_reps))

total_exposed<-readRDS("ComEO_tot_exposure_1-5000.RData")

cml_exposure<-array(as.numeric(NA),dim=c(n_fish,n_years,total_reps))

cml_exposure<-readRDS("ComEO_cml_exposure_1-5000.RData")


#Convert to average #/% fish exposed by fishery each year

#-----Get total and mean exposure for all fisheries together

mean_exposed<-array(as.numeric(NA),dim=c(n_fisheries,n_years))
sd_exposed<-array(as.numeric(NA),dim=c(n_fisheries,n_years))

mean_perc_exposed<-array(as.numeric(NA),dim=c(n_fisheries,n_years))
sd_perc_exposed<-array(as.numeric(NA),dim=c(n_fisheries,n_years))

for(y in 1:n_years){
  for(f in 1:n_fisheries){
    mean_exposed[f,y]<-mean(total_exposed[f,y,],na.rm=T)
    mean_perc_exposed[f,y]<-mean_exposed[f,y]/n_fish*100
    sd_exposed[f,y]<-sd(total_exposed[f,y,],na.rm=T)
    sd_perc_exposed[f,y]<-sd_exposed[f,y]/n_fish*100
  }
}  

#------------Get cumulative exposure to fisheries
  
cml_exp_iters<-array(as.numeric(NA),dim=c(n_fisheries+1,n_years,total_reps))
mean_cml_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,n_years))
sd_cml_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,n_years))
mean_cml_perc_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,n_years))
sd_cml_perc_exp<-array(as.numeric(NA),dim=c(n_fisheries+1,n_years))

for(f in 0:n_fisheries){
  for(y in 1:n_years){
   for(i in 1:total_reps){
    cml_exp_iters[f+1,y,i]<-sum(cml_exposure[,y,i]==f)
   }
  }
}

for(y in 1:n_years){
   for(f in 1:(n_fisheries+1)){
    mean_cml_exp[f,y]<-mean(cml_exp_iters[f,y,],na.rm=T)
    mean_cml_perc_exp[f,y]<-mean_cml_exp[f,y]/n_fish*100
    sd_cml_exp[f,y]<-sd(cml_exp_iters[f,y,],na.rm=T)
    sd_cml_perc_exp[f,y]<-sd_cml_exp[f,y]/n_fish*100
  }
}

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
    incr_exp[2,y,i]<-incr_exp[2,y,i]+incr_exp[1,y,i]  #remove for now, don't want the cumulative incremental exposure
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
    incr_exp[3,y,i]<-incr_exp[3,y,i]+incr_exp[2,y,i]
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
    incr_exp[4,y,i]<-incr_exp[4,y,i]+incr_exp[3,y,i]
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
    incr_exp[5,y,i]<-incr_exp[5,y,i]+incr_exp[4,y,i]
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
    incr_exp[6,y,i]<-incr_exp[6,y,i]+incr_exp[5,y,i]
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
    incr_exp[7,y,i]<-incr_exp[7,y,i]+incr_exp[6,y,i]
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
    incr_exp[8,y,i]<-incr_exp[8,y,i]+incr_exp[7,y,i]
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
##############################
#Print plots to pdf file
##############################
##############################

setwd(plots_dir)

###############################
#Run timing plots
###############################
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


#####################################
#Percentage of run exposed by year
#####################################

pdf(file="Population Percent Exposure to Commercial Fisheries by Year - Barplots.pdf")
par(mfrow=c(1,1),mar=c(5.5,4.1,4.1,3.1), oma=c(1,1,1,1), xpd=FALSE)

yr=2004
for(y in 1:n_years){
  
  y1<-array(as.numeric(NA),dim=c(1,n_fisheries))
  for(f in 1:n_fisheries){
    y1[f]<-max(0,mean_perc_exposed[f,y]-sd_perc_exposed[f,y])
  }
  
  y2<-array(as.numeric(NA),dim=c(1,n_fisheries))
  for(f in 1:n_fisheries){
    y2[f]<-min(100,mean_perc_exposed[f,y]+sd_perc_exposed[f,y])
  }


barCenters<-barplot(mean_perc_exposed[,y], main=paste0("Population Percent Exposure in ",yr), xlab="",ylab="% Exposure", 
        names.arg=fishery_names, axis.lty=1, ylim=range(0,100),las=3)
  arrows(barCenters, y1, barCenters, y2, length=0.05, angle=90, code=3)
  
  for(f in 1:n_fisheries){
    if(openings[f,y]=="FALSE"){
       text(barCenters[f],3,"*")
    }
  }
  
  text(20,90,adj=c(0,0),"* No fishery openings")
  

  yr=yr+1
 
#if(y==6 | y==12 | y==13){  
 # mtext(text="Fishery Area",side=1,line=1,outer=TRUE)
 #mtext(text="% Exposure",side=2,line=1,outer=TRUE)
 #mtext(text="Population Percent Exposure by Year",side=3,line=1,outer=TRUE)
 }

dev.off()

#--------Plots of probability of exposure with error bars

#Percentage of run exposed by fishery

x<-1:13
pdf(file=paste0("Population Percent Exposure by Fishery - Line plots w Error bars.pdf"))
#par(mfrow=c(1,1),mar=c(3,3,1,1), oma=c(5,5,3,1))
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1), oma=c(1,1,1,1))

for(f in 1:n_fisheries){
  plot(mean_perc_exposed[f,], main=fishery_names[f], xlab="",ylab="% Exposed", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1), las=2)
    arrows(x, mean_perc_exposed[f,]-sd_perc_exposed[f,], x, mean_perc_exposed[f,]+sd_perc_exposed[f,], length=0.05, angle=90, code=3, col="red")
    points(mean_perc_exposed[f,],pch=ifelse(openings[f,x]=="TRUE",16,1))
    
  legend("topright", legend="No openings",pch=1,bty="n")
}

   
#Turning off for now, can turn on if you want to put more than one graph per page
#mtext(text="Year",side=1,line=1,outer=TRUE)
#mtext(text="% Exposed",side=2,line=1,outer=TRUE)
#mtext(text="Population Percent Exposure by Fishery",side=3,line=1,outer=TRUE)
dev.off()

#-----------Plots of cumulative exposure--------------------------

pdf(file=paste0("Population Cumulative Perc Exposure to Commercial Fisheries - Bar Plots.pdf"))
#par(mfrow=c(1,1),mar=c(3,3,1,1), oma=c(5,5,3,1))
par(mfrow=c(1,1),mar=c(5.1,3.1,3.1,4.1), oma=c(1,1,1,1))

palette(c(grDevices::rainbow(17)))

#Plot of all years together
barplot(mean_cml_perc_exp, main="Average Cumulative Exposure \nto Commercial Fisheries by Year",
          xlab="Year", col=1:17,
          ylab="% Exposure",
          legend = seq(0,16,by=1),names.arg=c("04","05","06","07","08","09","10","11","12","13","14","15","16"),
          args.legend=list(x=18,y=100,bty="n"))

#Plots for each year separately
par(mfrow=c(1,1),mar=c(5.1,4.1,3.1,0.1), oma=c(1,0.5,1,0.5))

yr=2004
for(y in 1:n_years){
  
  y1<-array(as.numeric(NA),dim=c(1,n_fisheries+1))
  for(f in 1:n_fisheries){
    y1[f]<-max(0,mean_cml_perc_exp[f,y]-sd_cml_perc_exp[f,y])
  }
  
  y2<-array(as.numeric(NA),dim=c(1,n_fisheries+1))
  for(f in 1:n_fisheries){
    y2[f]<-min(100,mean_cml_perc_exp[f,y]+sd_cml_perc_exp[f,y])
  }
  
barCenters<-barplot(mean_cml_perc_exp[1:17,y], main=paste0("Average Cumulative Exposure \nto Commercial Fisheries in ",yr),
        xlab="# of fisheries each fish exposed to", col=1:17,
        ylab="% Exposure", ylim=range(0,100),names.arg=seq(0,16,by=1))
  arrows(barCenters, y1, barCenters, y2, length=0.05, angle=90, code=3)

  yr=yr+1
}

dev.off()

#-------------Plots of incremental exposure------------------------

#---Currently set up to calculate CUMULATIVE incr. exposure in the AllCom section and PDF title

pdf(file=paste0("Population Cumulative Incremental Exposure by ",data_source," Fishery - Line plots w Error bars.pdf"))
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
    plot(mean_incr_perc_exp[,y], main=paste0("Cumulative Incremental Exposure to Commercial Fisheries in ",yr), xlab="Fishery",ylab="% Exposed", xaxt="n", type="l", cex.main=0.8, bty="n", lty=1, ylim=range(0,100))
    axis(1, at=1:8, cex.axis=0.8,labels=c("Area G","Area B","Area D", "Area H", "Area E", "BPM GN", "APM GN", "APM BSn"), las=1)
    arrows(x, y1[,y], x, y2[,y], length=0.05, angle=90, code=3, col="red")
    points(mean_incr_perc_exp[,y],pch=16)
    yr=yr+1
  }
  
}else{ #data_source=="REC"
  
}
dev.off()
