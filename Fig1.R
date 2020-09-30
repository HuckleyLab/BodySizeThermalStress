#Analysis for a NCC news and views

library(zoo)
library("accelerometry")
library(viridisLite)
library(rphylopic) #for beetle icon

cols= viridis(10)[c(2,5,9)]

#load climate data
#Darrington, WA surface temp, 5 miute interval
#ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/subhourly01/2019/
#data from 2013-2019 available, use 2019

setwd("./data/")

# 4    LST_DATE                       YYYYMMDD
# 5    LST_TIME                       HHmm
# 9    AIR_TEMPERATURE                Celsius
# 13   SURFACE_TEMPERATURE            Celsius
# 14   ST_TYPE                        X
# 15   ST_FLAG                        X

clim= read.table("CRNS0101-05-2019-WA_Darrington_21_NNE.txt", na.strings = "-9999.0")
clim=clim[,c(4,5,9,13)]
names(clim)<- c("date","time","Tair","Tsurf")

#load beetle data
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12789
#Scaphinotus angusticollis

#body masses to model
#average dry mass= 31.3
# -0.12% change per year
#1915 to 2015
#/0.4 to convert to wet mass as in paper
#/1000 to convert mg to g
ml= (31.3 +0.0012*50*31.3)/0.4/1000
ms= (31.3 -0.0012*50*31.3)/0.4/1000

#estimate 0.5% quantile of mass
set.seed(1)
size 	<- rnorm(10000000,31.3,9.6)/0.4/1000
mf <- quantile(size,0.005)

#------------
#FUNCTIONS
#estimate Tcrit
tcrit= function(m,t){
  ctmax= 41.92 -1.65*log10(m)
  z= 2.85-0.45*log10(m)
  tc= ctmax-z*log10(t)
  return(tc)
}

#number heat stress events
#y= TSM, x=time
#find high resolution weather data
#use rolling average to calculate TSM for each time period

#Average temperature to different time periods then estimate maximum 
#n is with of period for rolling average
temp=clim$Tsurf

rollmax= function(n) {
  max(movingaves(temp, n), na.rm=TRUE)
}

ns=1:2016
#ns= c(1:12,(2:12)*12, seq(13,168,12)*12,168*12)

rt= sapply(ns, FUN=rollmax)

#----
#FIGURE

#Fig 1A. Tcrit and Ta over time
hrs=ns

#Estimate Tcrits for 3 sizes
tcrits.l= tcrit(ml,hrs*60)
tcrits.s= tcrit(ms,hrs*60)
tcrits.f= tcrit(mf,hrs*60)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/NCCnews/figures/")
pdf("Fig1.pdf", height = 5, width = 10)
par(mfrow=c(1,3), mgp=c(2.5,1,0))

#labels
#labs= c("1 hour","6 hours","1 day", "1 week")
#xs= log10( c(1,6,24, 168))
labs= c("5 minutes","1 hour","6 hours","1 day", "1 week")
xs= log10( c(1,1*12,6*12,24*12, 168*12))

plot(log10(hrs), tcrits.l, col=cols[3], type="l", 
     ylab="temperature (°C)", xlab="",cex.lab=1.8, xaxt="n", ylim=c(21,39), lwd=2, 
     cex.axis=1.4)
points(log10(hrs), tcrits.s, col=cols[2], type="l", lwd=2)
points(log10(hrs), tcrits.f, col=cols[1], type="l", lwd=2)
#plot temps
points(log10(hrs), rt, type="l", lty="dashed", lwd=2)
text(0.4,39.2,"a", cex=1.8)
#update axis
axis(1, at=xs, labels=labs, cex.axis=1.4)
  
# Add a legend
legend("bottomleft", legend=c("1915: 82.9mg", "2015: 73.5mg","smallest: 16.4mg"),
       col=cols[3:1], lty=c("solid","solid","solid"), cex=1.6, bty="n")
#"environmental temperature" ,"dashed"

#add beetle
#img <- image_data("dd92c3e8-ad3d-485d-b551-12961476e8d8", size = "1")[[1]]
#add_phylopic_base(img, 1, 27, 4)

#-----------------
#Fig 1B. Plot change in maximum metabolism at Tcrit as a function of exposure time
a <- 2
b <- 0.75
Q10 <- 2.5
M <- 10^c(-6:3)
c <- log10(Q10)/10
t <- c(1,60,24*60,7*24*60)
BM <- - 1.630
BtM <- - 0.45
bprim <- b + c*(BM - BtM*log10(t))

Bo=41.92
Bt= 2.85

log10MR= function(m,t) {log10(a)+c*(Bo-Bt*log10(t))+(b+c*(BM-BtM*log10(t)))*log10(m)}

plot(log10(hrs), log10MR(ml,hrs*60), type="l", ylab="log10 metabolic rate (W)", xlab="log10 time", 
     col=cols[3],cex.lab=1.8, xaxt="n", ylim= range(0.1,1), lwd=2, cex.axis=1.4) 
points(log10(hrs), log10MR(ms,hrs*60), col=cols[2], type="l", lwd=2)
points(log10(hrs), log10MR(mf,hrs*60), col=cols[1], type="l", lwd=2)
text(0,1.02,"b", cex=1.8)
#update axis
axis(1, at=xs, labels=labs, cex.axis=1.4)

#slopes
mc= c(ml, ms, mf)
lm(log10MR(mf,hrs*60)~ log10(hrs))$coefficients[2]

#text
text(1,0.9,"slope= -0.133",col=cols[3], srt=-35, cex=1.6)
text(1,0.78,"slope= -0.134",col=cols[2], srt=-35, cex=1.6)
text(1,0.4,"slope= -0.145",col=cols[1], srt=-35, cex=1.6)

#-----------------
#Fig 1C. Plot TSM as a function of exposure time

#Estimate TSM
tsm.l= tcrits.l -rt 
tsm.s= tcrits.s -rt 
tsm.f= tcrits.f -rt

plot(log10(hrs), tsm.l, type="l", ylab="thermal safety margin (TSM, °C)", xlab="", 
     col=cols[3],cex.lab=1.8, xaxt="n", lwd=2, cex.axis=1.4)
points(log10(hrs), tsm.s, col=cols[2], type="l", lwd=2)
points(log10(hrs), tsm.f, col=cols[1], type="l", lwd=2)
text(0,7.4,"c", cex=1.8)
#update axis
axis(1, at=xs, labels=labs, cex.axis=1.4)

dev.off()

#------
#metrics
#min TSM
mean(c(which.min(tsm.l)*5/60, which.min(tsm.s)*5/60, which.min(tsm.f)*5/60))
#TSM<0
mean(c(which.min(tsm.l<0)*5/60, which.min(tsm.l<0)*5/60, which.min(tsm.l<0)*5/60))

#difference in Tcrit at 5 minutes
tcrits.f[1]-tcrits.l[1]


