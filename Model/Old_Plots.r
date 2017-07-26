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

#-----------------------Generate line plots for Mike -------------------------------

pdf(file=paste0("Population Exposure by Fishery ",data_source,"- Line plots.pdf"))
par(mfrow=c(1,1),mar=c(3,3,1,1), oma=c(5,5,3,1))

if(data_source=="Commercial"){
  plot(total_exposed[1,], main="Area B", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[1,,])
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
}else if(data_source=="FN"){
  plot(total_exposed[1,], main="APM BSn", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[1,,])
  plot(total_exposed[2,], main="APM DN", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[2,])
  plot(total_exposed[3,], main="APM SN", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[3,])
  plot(total_exposed[4,], main="BPM DN", xlab="",ylab="", xaxt="n", type="l", bty="n", lty=1, ylim=range(0,1000))
  axis(1, at=1:13,labels=seq(from=2004,to=2016, by=1))
  points(total_exposed[4,])
}
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
