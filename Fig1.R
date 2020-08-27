#Analysis for a NCC news and views

library(zoo)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/NCCnews/data/")

#load climate data
clim= read.csv("BoulderCr.csv")
temp=clim$TempC
#potential sources to update climate data:
#https://climate.weather.gc.ca/historical_data/search_historic_data_e.html
#https://data.pacificclimate.org/portal/pcds/map/
#ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/subhourly01/2019/

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
#also assume 75% reduction
mf=ml*0.25

#FUNCTIONS
#estimate Tcrit
tcrit= function(m,t){
  ctmax= 41.93 -1.63*log10(m)
  z= 2.86-0.45*log10(m)
  tc= ctmax-z*log10(t)
  return(tc)
}

#number heat stress events
#y= TSM, x=time
#find high resolution weather data
#use rolling average to calculate TSM for each time period

#Average temperature to different time periods then estimate maximum 
#n is with of period for rolling average
rollmax= function(n) {
  max(rollmean(temp, n))
}

#----
#FIGURE

#Fig 1A. Tcrit and Ta over time
hrs=1:168
rt= sapply( hrs, FUN=rollmax)

#Estimate Tcrits for 3 sizes
tcrits.l= tcrit(ml,hrs*60)
tcrits.s= tcrit(ms,hrs*60)
tcrits.f= tcrit(mf,hrs*60)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/NCCnews/figures/")
pdf("Fig0.pdf", height = 6, width = 10)
par(mfrow=c(1,3))

#labels
labs= c("1 hour","6 hours","1 day", "1 week")
xs= log10( c(1,6,24, 168))

plot(log10(hrs), tcrits.l, col="blue", type="l", 
     ylab="temperature (°C)", xlab="log10 time",cex.lab=1.2, xaxt="n")
points(log10(hrs), tcrits.s, col="orange", type="l")
points(log10(hrs), tcrits.f, col="green", type="l")
#plot temps
points(log10(hrs), rt, type="l", lty="dashed")
text(0.4,37.5,"a")
#update axis
axis(1, at=xs, labels=labs)

# Add a legend
legend("topright", legend=c("1915: 82.9mg", "2015: 73.5mg","20.7mg", "Ta"),
       col=c("blue", "orange", "green","black"), lty=c("solid","solid","solid","dashed"), cex=1.2)

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
Bt= 2.85  #CHECK SIGN

log10MR= function(m,t) {log10(a)+c*(Bo-Bt*log10(t))+(b+c*(BM-BtM*log10(t)))*log10(m)}

plot(log10(hrs), log10MR(ml,hrs*60), type="l", ylab="log10 metabolic rate (W)", xlab="log10 time", 
     col="blue",cex.lab=1.2, xaxt="n", ylim= range(0.25,1) ) 
points(log10(hrs), log10MR(ms,hrs*60), col="orange", type="l")
points(log10(hrs), log10MR(mf,hrs*60), col="green", type="l")
text(0,1,"b")
#update axis
axis(1, at=xs, labels=labs)

#slopes
mc= c(ml, ms, mf)
lm(log10MR(mf,hrs*60)~ log10(hrs))$coefficients[2]

#text
text(1,0.9,"slope= -0.133",col="blue", srt=-35, cex=1.2)
text(1,0.78,"slope= -0.134",col="orange", srt=-35, cex=1.2)
text(1,0.47,"slope= -0.144",col="green", srt=-35, cex=1.2)

#-----------------
#Fig 1C. Plot TSM as a function of exposure time

#Estimate TSM
tsm.l= tcrits.l -rt 
tsm.s= tcrits.s -rt 
tsm.f= tcrits.f -rt

plot(log10(hrs), tsm.l, type="l", ylab="thermal safety margin (TSM, °C)", xlab="log10 time", 
     col="blue",cex.lab=1.2, xaxt="n")
points(log10(hrs), tsm.s, col="orange", type="l")
points(log10(hrs), tsm.f, col="green", type="l")
text(0,8,"c")
#update axis
axis(1, at=xs, labels=labs)

dev.off()
