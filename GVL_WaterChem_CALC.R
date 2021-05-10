# Project: Spatiotemporal Variation in Internal P Loading in a Hypereutrophic Reservoir (GVL) - Ecosystems
# Last Modified 29 April 2021
# Contributers: Ellen Albright
# Description: The following script details calculations used to determine total suspended solids and QAQC visualizations of water column nutrients.

# The code can be run with the following datasets:
# "GVL_WaterQuality_RAW.csv" - raw water quality data (water column TP, SRP, TN, NOx, TSS, ISS, VSS, Microcystin) 2014-2020

# REQUIRED PACKAGES FOR CODE - install as needed
# install.packages("tidyverse")
library(tidyverse)

# Clear environment, set working directory - UPDATE AS NEEDED --------------------------------------------------------------------------------------------------------------
rm(list=ls())
getwd()
# Set working directory to "WQ" folder within GVL 2020 > Data folders. This contains the current, compiled WQ csv file
setwd("C:/Users/Ellen/Desktop/Box Sync/Albright DISSERTATION/Albright_Chapter_2/DATA/DATA_PUBLISH")

# Working with complete nutrient dataset - includes 2020, 2019 projects and 319 sampling from 2014-5
wq<-read.csv("GVL_WaterQuality_RAW.csv")
str(wq)
names(wq)

# QAQC by project
# 1 "sediment" - 2020 sampling, surface and bottom water, sites 1, 4, and 10
sediment<-subset(wq,project=="sediment")
P_diff<-sediment$TP_ugL-sediment$SRP_ugL #check that TP is always greater than SRP
min(P_diff,na.rm=T) #minimum difference is 28.15 ug/L (no negative values)
N_diff<-sediment$TN_mgL-sediment$NOx_mgL
min(N_diff,na.rm=T) #minimum difference is 0.01 mg/L (no negative values)

TSScol<-"#1a1a1a"
VSScol<-"#80cdc1"
ISScol<-"#bf812d"
TPcol<-"#023858"
SRPcol<-"#a6bddb"
TNcol<-"#67000d"
NOXcol<-"#fb6a4a"

# SURFACE WATER - solids, phosphorus, and nitrogen by site (confirm dissolved nutrients are less than total and ISS+VSS is less than TSS)
par(mfrow=c(3,3),cex=1.1,oma=c(1,1,0,0), mar = c(1.5,2,1,1.5))
plot((sediment[sediment$depthID=="1"&sediment$siteID=="1","TSS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,bg=TSScol,lwd=3,col=TSScol,
     xlab="DOY",ylab="Suspended Solids (mg/L)",ylim=c(0,60),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","ISS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=21,bg='white',lwd=3,col=ISScol)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","VSS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=21,bg='white',lwd=3,col=VSScol)
legend(10,60,legend=c('TSS','VSS','ISS'),col=c(TSScol,VSScol,ISScol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,60,"Shallow Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="1"&sediment$siteID=="10","TSS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="10","doy"]),type='o',pch=17,cex=1.2,bg=TSScol,lwd=3,col=TSScol,
     xlab="DOY",ylab="Suspended Solids (mg/L)",ylim=c(0,60),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="10","ISS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="10","doy"]),type='o',pch=24,bg='white',lwd=3,col=ISScol)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="10","VSS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="10","doy"]),type='o',pch=24,bg='white',lwd=3,col=VSScol)
legend(10,60,legend=c('TSS','VSS','ISS'),col=c(TSScol,VSScol,ISScol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,60,"Middle Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="1"&sediment$siteID=="4","TSS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,bg=TSScol,lwd=3,col=TSScol,
     xlab="DOY",ylab="Suspended Solids (mg/L)",ylim=c(0,60),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","ISS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg='white',lwd=3,col=ISScol)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","VSS_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg='white',lwd=3,col=VSScol)
legend(10,60,legend=c('TSS','VSS','ISS'),col=c(TSScol,VSScol,ISScol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,60,"Deep Site",cex=0.8,font=2)

plot((sediment[sediment$depthID=="1"&sediment$siteID=="1","TP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,bg=TPcol,lwd=3,col=TPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","SRP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=21,bg='white',lwd=3,col=SRPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,400,"Shallow Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="1"&sediment$siteID=="10","TP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="10","doy"]),type='o',pch=17,cex=1.2,bg=TPcol,lwd=3,col=TPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="10","SRP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="10","doy"]),type='o',pch=24,bg='white',lwd=3,col=SRPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,400,"Middle Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="1"&sediment$siteID=="4","TP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,bg=TPcol,lwd=3,col=TPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","SRP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg='white',lwd=3,col=SRPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,400,"Deep Site",cex=0.8,font=2)

plot((sediment[sediment$depthID=="1"&sediment$siteID=="1","TN_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,bg=TNcol,lwd=3,col=TNcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,2),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","NOx_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=21,bg='white',lwd=3,col=NOXcol)
legend(20,2,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,2,"Shallow Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="1"&sediment$siteID=="10","TN_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="10","doy"]),type='o',pch=17,cex=1.2,bg=TNcol,lwd=3,col=TNcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,2),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="10","NOx_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="10","doy"]),type='o',pch=24,bg='white',lwd=3,col=NOXcol)
legend(20,2,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,2,"Middle Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="1"&sediment$siteID=="4","TN_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,bg=TNcol,lwd=3,col=TNcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,2),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","NOx_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg='white',lwd=3,col=NOXcol)
legend(20,2,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,2,"Deep Site",cex=0.8,font=2)

# BOTTOM WATER - solids, phosphorus, and nitrogen by site (confirm dissolved nutrients are less than total and ISS+VSS is less than TSS)
par(mfrow=c(3,3),cex=1.1,oma=c(1,1,0,0), mar = c(1.5,2,1,1.5))
plot((sediment[sediment$depthID=="4"&sediment$siteID=="1","TSS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,bg=TSScol,lwd=3,col=TSScol,
     xlab="DOY",ylab="Suspended Solids (mg/L)",ylim=c(0,60),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","ISS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=21,bg='white',lwd=3,col=ISScol)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","VSS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=21,bg='white',lwd=3,col=VSScol)
legend(10,60,legend=c('TSS','VSS','ISS'),col=c(TSScol,VSScol,ISScol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,60,"Shallow Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="4"&sediment$siteID=="10","TSS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="10","doy"]),type='o',pch=16,cex=1.2,bg=TSScol,lwd=3,col=TSScol,
     xlab="DOY",ylab="Suspended Solids (mg/L)",ylim=c(0,60),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="10","ISS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="10","doy"]),type='o',pch=21,bg='white',lwd=3,col=ISScol)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="10","VSS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="10","doy"]),type='o',pch=21,bg='white',lwd=3,col=VSScol)
legend(10,60,legend=c('TSS','VSS','ISS'),col=c(TSScol,VSScol,ISScol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,60,"Middle Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="4"&sediment$siteID=="4","TSS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=16,cex=1.2,bg=TSScol,lwd=3,col=TSScol,
     xlab="DOY",ylab="Suspended Solids (mg/L)",ylim=c(0,60),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","ISS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=21,bg='white',lwd=3,col=ISScol)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","VSS_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=21,bg='white',lwd=3,col=VSScol)
legend(10,60,legend=c('TSS','VSS','ISS'),col=c(TSScol,VSScol,ISScol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,60,"Deep Site",cex=0.8,font=2)

plot((sediment[sediment$depthID=="4"&sediment$siteID=="1","TP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=17,cex=1.2,bg=TPcol,lwd=3,col=TPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","SRP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=24,bg='white',lwd=3,col=SRPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,400,"Shallow Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="4"&sediment$siteID=="10","TP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="10","doy"]),type='o',pch=17,cex=1.2,bg=TPcol,lwd=3,col=TPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="10","SRP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="10","doy"]),type='o',pch=24,bg='white',lwd=3,col=SRPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,400,"Middle Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="4"&sediment$siteID=="4","TP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=17,cex=1.2,bg=TPcol,lwd=3,col=TPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","SRP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=24,bg='white',lwd=3,col=SRPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,400,"Deep Site",cex=0.8,font=2)

plot((sediment[sediment$depthID=="4"&sediment$siteID=="1","TN_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=15,cex=1.2,bg=TNcol,lwd=3,col=TNcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,2),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","NOx_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=22,bg='white',lwd=3,col=NOXcol)
legend(20,2,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,2,"Shallow Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="4"&sediment$siteID=="10","TN_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="10","doy"]),type='o',pch=15,cex=1.2,bg=TNcol,lwd=3,col=TNcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,2),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="10","NOx_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="10","doy"]),type='o',pch=22,bg='white',lwd=3,col=NOXcol)
legend(20,2,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,2,"Middle Site",cex=0.8,font=2)
plot((sediment[sediment$depthID=="4"&sediment$siteID=="4","TN_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,bg=TNcol,lwd=3,col=TNcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,2),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","NOx_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=22,bg='white',lwd=3,col=NOXcol)
legend(20,2,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=0.8, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,2,"Deep Site",cex=0.8,font=2)

# Closer look at bottom water nutrients
sediment$SRP_mgL<-sediment$SRP_ugL/1000
sediment$TP_mgL<-sediment$TP_ugL/1000
# Deep Site (bottom water depth id=4, site=4)
par(mfrow=c(1,1),cex=1.1,mar=c(5.1,4.1,1.5,1.5))
plot((sediment[sediment$depthID=="4"&sediment$siteID=="4","SRP_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col=SRPcol,
     xlab="DOY",ylab="... (mg/L)",ylim=c(0,1),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","TP_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col=TPcol)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","NOx_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col=NOXcol)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","TN_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col=TNcol)

plot((sediment[sediment$depthID=="4"&sediment$siteID=="4","TP_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","TSS_mgL"]))
plot((sediment[sediment$depthID=="4"&sediment$siteID=="4","SRP_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","TSS_mgL"]))


# Surface water (depth id = 1) nutrients at the deep site (4) to compare to beach monitoring results
par(mfrow=c(1,1),cex=1.1,mar=c(5.1,4.1,1.5,1.5))
plot((sediment[sediment$depthID=="1"&sediment$siteID=="4","SRP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col=SRPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","TP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col=TPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=1.2, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(50,400,"Deep Site",cex=1.2,font=2)

plot((sediment[sediment$depthID=="1"&sediment$siteID=="4","NOx_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col=NOXcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,1.75),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","TN_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col=TNcol)
legend(20,1.75,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=1.2, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(50,1.75,"Deep Site",cex=1.2,font=2)


# 2 "beachMonitoring" - IDNR beach samples site 9 (2019-2020) surface. 2019 just MC
beach20<-subset(wq,project=="beachMonitoring"&year=="2020")
P_diff_beach20<-beach20$TP_ugL-beach20$SRP_ugL #check that TP is always greater than SRP
min(P_diff_beach20,na.rm=T) #minimum difference is 76.15 ug/L (no negative values)
N_diff_beach20<-beach20$TN_mgL-beach20$NOx_mgL
min(N_diff_beach20,na.rm=T) #minimum difference is -0.04 mg/L (likely a below detection NOx that needs to be excluded)

par(mfrow=c(1,1),cex=1.1,mar=c(5.1,4.1,1.5,1.5))
plot((beach20$SRP_ugL~beach20$doy),type='o',pch=22,bg="white",cex=1.2,lwd=3,col=SRPcol,
     xlab="DOY",ylab="Phosphorus (ug/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((beach20$TP_ugL~beach20$doy),type='o',pch=15,cex=1.2,lwd=3,col=TPcol)
legend(20,400,legend=c('TP','SRP'),col=c(TPcol,SRPcol),lwd=3,cex=1.2, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,400,"Beach Site (IDNR)",cex=1.2,font=2)

plot((beach20$NOx_mgL~beach20$doy),type='o',pch=22,bg="white",cex=1.2,lwd=3,col=NOXcol,
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,1.75),xlim=c(35,300),las=1)
lines((beach20$TN_mgL~beach20$doy),type='o',pch=15,cex=1.2,lwd=3,col=TNcol)
legend(20,1.75,legend=c('TN','NOx'),col=c(TNcol,NOXcol),lwd=3,cex=1.2, y.intersp=0.7,x.intersp=0.3,bty='n',seg.len=0.5)
text(75,1.75,"Beach Site (IDNR)",cex=1.2,font=2)

# 3 "dnr319" - just TP at the deep site - will compare site 4 surface (1) water across all available years (2014, 2015, 2019, 2020)
wq2020<-filter(wq, year=="2020"&siteID=="4"&depthID=="1")
wq2019<-filter(wq, year=="2019"&siteID=="4"&depthID=="1")
wq2015<-filter(wq, year=="2015"&siteID=="4"&depthID=="1")
wq2014<-filter(wq, year=="2014"&siteID=="4"&depthID=="1")

wq2020b<-filter(wq, year=="2020"&siteID=="4"&depthID=="4")
wq2015b<-filter(wq, year=="2015"&siteID=="4"&depthID=="4")
wq2014b<-filter(wq, year=="2014"&siteID=="4"&depthID=="4")

plot((wq2020$TP_ugL~wq2020$doy),type="o",lwd="3",pch=15,
     xlab="DOY",ylab="Total Phosphorus (ug/L)",ylim=c(0,670),xlim=c(35,300),las=1)
lines((wq2019$TP_ugL~wq2019$doy),type="o",lwd="3",pch=15,col="#e7298a")
lines((wq2015$TP_ugL~wq2015$doy),type="o",lwd="3",pch=15,col="#d95f02")
lines((wq2014$TP_ugL~wq2014$doy),type="o",lwd="3",pch=15,col="#1b9e77")
legend(40,650,legend=c('2020 Sediment Project ','2019 GVL Monitoring','2015 DNR 319','2014 DNR 319'),col=c("black","#e7298a","#d95f02","#1b9e77"),lwd=3,cex=1, y.intersp=0.8,bty='n')
text(60,650,"Deep Site",cex=1.2,font=2)

par(mfrow=c(1,1),cex=1.1,mar=c(5.1,4.1,1.5,1.5))
plot((wq2020b$TP_ugL~wq2020b$doy),type="o",lwd="3",pch=15,
     xlab="DOY",ylab="Total Phosphorus (ug/L)",ylim=c(0,670),xlim=c(35,300),las=1)
lines((wq2015b$TP_ugL~wq2015b$doy),type="o",lwd="3",pch=15,col="#d95f02")
lines((wq2014b$TP_ugL~wq2014b$doy),type="o",lwd="3",pch=15,col="#1b9e77")
legend(40,650,legend=c('2020 Sediment Project ','2015 DNR 319','2014 DNR 319'),col=c("black","#d95f02","#1b9e77"),lwd=3,cex=1, y.intersp=0.8,bty='n')
text(60,650,"Deep Site",cex=1.2,font=2)

# 4 "GVLmonitoring"
GVL2019<-subset(wq,project=="GVLmonitoring")
P_diff_GVL2019<-GVL2019$TP_ugL-GVL2019$SRP_ugL #check that TP is always greater than SRP
min(P_diff_GVL2019,na.rm=T) #minimum difference is 23.65 ug/L (no negative values)
N_diff_GVL2019<-GVL2019$TN_mgL-GVL2019$NOx_mgL
min(N_diff_GVL2019,na.rm=T) #minimum difference is -0.32 mg/L (.....site 1, surface, doy 164. )

# Comparing 2019 and 2020 - site 4, surface and bottom water
par(mfrow=c(2,2),cex=1.1,mar=c(1.5,4.1,1.5,1))
plot((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","TP_ugL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col="gray35",
     xlab="DOY",ylab="Phosphorus (µg/L)",ylim=c(0,400),xlim=c(35,300),las=1)
lines((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","SRP_ugL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="gray60")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","TP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col="#2171b5")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","SRP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="#6baed6")
legend(40,400,legend=c('2019 TP','2019 SRP','2020 TP','2020 SRP'),col=c("gray35","gray60","#2171b5","#6baed6"),pch=c(15,22,15,22),pt.bg=c("gray35","white","#2171b5","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,400,"Deep Site SURFACE",cex=1,font=2)

plot((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","TN_mgL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col="gray35",
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,8),xlim=c(35,300),las=1)
lines((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","NOx_mgL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="gray60")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","TN_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col="#ef3b2c")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="4","NOx_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="#fc9272")
legend(170,8.5,legend=c('2019 TN','2019 Nox','2020 TN','2020 NOx'),col=c("gray35","gray60","#ef3b2c","#fc9272"),pch=c(15,22,15,22),pt.bg=c("gray35","white","#ef3b2c","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,8,"Deep Site SURFACE",cex=1,font=2)

plot((GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="4","SRP_ugL"])~(GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="gray60",
     xlab="DOY",ylab="Phosphorus (µg/L)",ylim=c(0,1000),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","TP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col="#2171b5")
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","SRP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="#6baed6")
legend(40,1000,legend=c('2019 SRP','2020 TP','2020 SRP'),col=c("gray60","#2171b5","#6baed6"),pch=c(22,15,22),pt.bg=c("white","#2171b5","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,1000,"Deep Site BOTTOM",cex=1,font=2)

plot((GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="4","NOx_mgL"])~(GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="gray60",
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,4),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","TN_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=15,cex=1.2,lwd=3,col="#ef3b2c")
lines((sediment[sediment$depthID=="4"&sediment$siteID=="4","NOx_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="4","doy"]),type='o',pch=22,bg="white",cex=1.2,lwd=3,col="#fc9272")
legend(40,4,legend=c('2019 Nox','2020 TN','2020 NOx'),col=c("gray60","#ef3b2c","#fc9272"),pch=c(22,15,22),pt.bg=c("white","#ef3b2c","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,4,"Deep Site BOTTOM",cex=1,font=2)

# Comparing 2019 and 2020 - site 1, surface and bottom water
par(mfrow=c(2,2),cex=1.1,mar=c(1.5,4.1,1.5,1))
plot((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","TP_ugL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","doy"]),type='o',pch=16,cex=1.2,lwd=3,col="gray35",
     xlab="DOY",ylab="Phosphorus (µg/L)",ylim=c(0,500),xlim=c(35,300),las=1)
lines((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","SRP_ugL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="gray60")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","TP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,lwd=3,col="#2171b5")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","SRP_ugL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="#6baed6")
legend(40,500,legend=c('2019 TP','2019 SRP','2020 TP','2020 SRP'),col=c("gray35","gray60","#2171b5","#6baed6"),pch=c(16,21,16,21),pt.bg=c("gray35","white","#2171b5","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,500,"Deep Site SURFACE",cex=1,font=2)

plot((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","TN_mgL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","doy"]),type='o',pch=16,cex=1.2,lwd=3,col="gray35",
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,8),xlim=c(35,300),las=1)
lines((GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","NOx_mgL"])~(GVL2019[GVL2019$depthID=="1"&GVL2019$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="gray60")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","TN_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,lwd=3,col="#ef3b2c")
lines((sediment[sediment$depthID=="1"&sediment$siteID=="1","NOx_mgL"])~(sediment[sediment$depthID=="1"&sediment$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="#fc9272")
legend(170,8.5,legend=c('2019 TN','2019 Nox','2020 TN','2020 NOx'),col=c("gray35","gray60","#ef3b2c","#fc9272"),pch=c(16,21,16,21),pt.bg=c("gray35","white","#ef3b2c","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,8,"Deep Site SURFACE",cex=1,font=2)

plot((GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="1","SRP_ugL"])~(GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="gray60",
     xlab="DOY",ylab="Phosphorus (µg/L)",ylim=c(0,500),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","TP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,lwd=3,col="#2171b5")
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","SRP_ugL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="#6baed6")
legend(40,500,legend=c('2019 SRP','2020 TP','2020 SRP'),col=c("gray60","#2171b5","#6baed6"),pch=c(21,16,21),pt.bg=c("white","#2171b5","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,500,"Deep Site BOTTOM",cex=1,font=2)

plot((GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="1","NOx_mgL"])~(GVL2019[GVL2019$depthID=="4"&GVL2019$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="gray60",
     xlab="DOY",ylab="Nitrogen (mg/L)",ylim=c(0,4),xlim=c(35,300),las=1)
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","TN_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=16,cex=1.2,lwd=3,col="#ef3b2c")
lines((sediment[sediment$depthID=="4"&sediment$siteID=="1","NOx_mgL"])~(sediment[sediment$depthID=="4"&sediment$siteID=="1","doy"]),type='o',pch=21,bg="white",cex=1.2,lwd=3,col="#fc9272")
legend(40,4,legend=c('2019 Nox','2020 TN','2020 NOx'),col=c("gray60","#ef3b2c","#fc9272"),pch=c(21,16,21),pt.bg=c("white","#ef3b2c","white"),lwd=3,cex=1, y.intersp=0.8,x.intersp=0.5,bty='n')
text(95,4,"Deep Site BOTTOM",cex=1,font=2)


# Background: TSS, VSS, and ISS calculation approach
# Clear environment, set working directory - UPDATE AS NEEDED --------------------------------------------------------------------------------------------------------------
rm(list=ls())
getwd()
# Set working directory to "WQ" folder within GVL 2020 > Data folders. This contains the current, compiled WQ csv file
setwd("C:/Users/Ellen/Desktop/Box Sync/Albright DISSERTATION/Albright_Chapter_2/DATA/WQ")

# Read in raw TSS data
tss_raw<-read.csv("GVL_TSS_2020.csv")
tss<-tss_raw %>% 
   filter(!is.na(year)) #remove NA rows... 

# Calculate TSS (mg/L) = (A-B)/V 
# A = mass of filter and dried residue - mg
# B = mass of filter (tare weight) - mg
# V = volume of sample filtered - L

tss$TSS_mgL<-((tss$Dried_Weight_g*1000)-(tss$Filter_Weight_g*1000))/(tss$Volume_Filtered_mL/1000)
tss$VSS_mgL<-((tss$Dried_Weight_g*1000)-(tss$Combusted_Weight_g*1000))/(tss$Volume_Filtered_mL/1000)
tss$ISS_mgL<-tss$TSS_mgL-tss$VSS_mgL

# QAQC - blanks and standards
blanks<-tss %>% 
   filter(sample_id=="Blank")
# Goal - 0 mg/L. DOY 39-223: blanks were within 1 mg/L of 0. DOY 298: blank was low (-6.4 mg/L) - DATA FLAG
standards<-tss %>% 
   filter(sample_id=="Standard")
# Goal - I can't remember what the accepted range is.... messaged Elena to find out. Note from Julia that 31.2 passed. So likely DOY 39, 156, 181 pass. 223 high and 117 low.

# Now remove blanks and standards 
tss<-tss %>% 
   filter(!sample_id %in% c("Blank","Standard"))
# Summarize average of lab duplicates for each site and event
tss_sum<-tss %>% 
   mutate(doy=factor(doy,levels=c("39","117","156","181","223","298"))) %>% 
   mutate(siteID=factor(siteID,levels=c("1","10","4"))) %>% 
   mutate(depthID=factor(depthID,levels=c("1","4"))) %>% 
   group_by(doy, siteID,depthID) %>% 
   summarize(mean_TSS = mean(TSS_mgL,na.rm=T),mean_VSS = mean(VSS_mgL,na.rm=T),mean_ISS = mean(ISS_mgL,na.rm=T))

