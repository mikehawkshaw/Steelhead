#Model to estimate 50% run timing for steelhead, both over time series and year-specific.

source("directories.R")
library("xtable")
setwd(data_dir)

#Read in and set up data
albion_annual<-read.table("steelhead_albion.csv",header=T, sep=",")

n_col<-dim(albion_annual)[2]

m<-rep(0,n_col)
s<-rep(0,n_col)

years<-seq(1995,1995+n_col-2)

#ignore first column because it has the day of year index and we don't want to fit a normal curve to that...

for(i in 2:n_col)
{
#normalize y to be a probability (sum = 1)
p = albion_annual[,i]/sum(albion_annual[,i],na.rm=T)
#compute weighted mean and standard deviation
m[i]=sum(albion_annual[,1]*p,na.rm=T)
s[i]=sqrt(sum((albion_annual[,1] - m[i])^ 2*p,na.rm=T))

#compute theoretical probabilities
xs = albion_annual[,1]
pth = dnorm(albion_annual[,1],m[i],s[i])

deviates[,i]<-p-pth

#plot fitted normal distribution to observations
#it's important to check these - because of the nature of the truncated data set there appears to be a bias towards a later run timing prediction.
#this is probably an artifact of fixing the IFR steelhead run timing with the decay curve - i don't knwo the BEST way to correct for this - though i have several ideas


plot(albion_annual[,1],p,type="b", col="blue",main=years[i-1],xlab="Julian Date", ylab="Corrected Albion Catch")
lines(albion_annual[,1],pth, col="dark red",lwd=2)

}

par(mfcol=c(1,3))
#If we're going to use these plots in publications/presentations we'll need to fix the x-axis title widths and/or add line breaks
m_mean = mean(m[2:n_col])
m_sd = sd(m[2:n_col])
m_aprox<-dnorm(min(albion_annual$jday):max(albion_annual$jday),m_mean,1.5*m_sd)
#plot(density(m[2:n_col]))
hist(m[2:n_col],prob=TRUE, breaks=20, main="",xlab="Distribution of means (N=20)")
lines(min(albion_annual$jday):max(albion_annual$jday),m_aprox,col="red")
mtext(" (a) ",cex=0.75,line=-1.5,adj=1)

s_mean = mean(s[2:n_col])
s_sd = sd(s[2:n_col])
s_aprox<-dnorm(5:50,s_mean,1.5*s_sd)
#plot(density(s[2:53]))
hist(s[2:n_col],prob=TRUE, breaks=20, main="",xlab="Distribution of standard deviations (N=20)")
lines(5:50,s_aprox,col="red")
mtext(" (b) ",cex=0.75,line=-1.5,adj=1)

plot(m[2:n_col],s[2:n_col],xlab="mean",ylab="standard deviations", bty="n")
mtext(" (c) ",cex=0.75,line=-1.5,adj=1)

par(mfcol=c(1,1))

plot(years,m[2:n_col], xlab="Year", ylab="50% date", main="No Time Trend in Date of 50% return")



#calculate mean p50 and SD from run timing data 
mean=m_mean
#SD*2.5 provides a better fit of the normal curve to the average daily catch (by eye! so don't use it but think about why?)
sd=s_mean

#Write each year's 50% date and SD to file
timing_df <- data.frame(years,m[2:n_col],s[2:n_col])
names(timing_df) <- c("year","mean","sd")
write.csv(file=paste(data_dir,"/steelhead_runtiming.csv", sep=""), x=timing_df,row.names=F)