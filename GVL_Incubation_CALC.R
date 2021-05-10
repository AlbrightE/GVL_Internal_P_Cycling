# Project: Spatiotemporal Variation in Internal P Loading in a Hypereutrophic Reservoir (GVL) - Ecosystems
# Last Modified 5 May 2021
# Contributers: Ellen Albright
# Description: The following script details calculations used to determine a daily sediment phosphorus (P) release rate from sediment core incubation data
#              It also details the data QAQC procedure used to produce the final dataset used for analysis

# The code can be run with the following datasets:
   # "GVL_Incubation_RAW.csv" - raw, hand-entered data from 5 sediment core incbuation experiments from Feb-Oct, 2020 as well as July 2019
# The code produces the following datasets:
   # "GVL_Incubation_Tidy.csv" - final dataset of calculated sediment P release rates
   # "GVL_Incubation_Sum.csv" - final dataset of calculated sediment P release rates, averaged by sediment core over the course of each incubation (used for analyses)

# Clear environment, set working directory - UPDATE AS NEEDED --------------------------------------------------------------------------------------------------------------
rm(list=ls())
getwd()
setwd('C:/Users/Ellen/Desktop/Box Sync/Albright DISSERTATION/Albright_Chapter_2/DATA/DATA_PUBLISH')

# REQUIRED PACKAGES FOR CODE - install as needed
# install.packages("tidyverse")
# install.packages("RColorBrewer")
library(tidyverse)
library(RColorBrewer)

### PART 1 - CALCULATIONS AND DATA CLEANING -------------------------------------------------------------------------------------------------------------------------------
### 2020, 2019 INCUBATIONS ------------------------------------------------------------------------------------------------------------------------------------------------
# Read in raw datafile
raw<-read.csv('GVL_Incubation_RAW.csv')
View(raw)
names(raw)

#1. Calculate the average of lab TP duplicates for the samples and the replacement water
TP_mgL<-(raw$TP_mgL_1+raw$TP_mgL_2)/2
replace_TP_mgL<-rowMeans(raw[,c(22,23)],na.rm=TRUE)

#2. Calculate the average volume of the water column for each core.
#   Average the pre and post incubation water column height values
water_column_m<-rowMeans(raw[,c(12,13)],na.rm=TRUE)
view(water_column_m) #make sure the NAs behaved
#   Multiply the average water column height (m) by the surface area of the sediment (m2; equals the interior cross-sectional area of the core sleeve, and therefore water column)
#   m3 --> L multiply by 1000
water_vol_L<-water_column_m*raw$sediment_sa_m2*1000

#3. Calculate the mass of TP in the water column of each sediment core immediately after a daily sample has been taken
#   Multiply the TP (mg/L) concentration by the volume of water (L) minus the replacement volume (L)
TP_mg<-TP_mgL*(water_vol_L-raw$vol_replace_L)

#4. Calculate the mass of TP in the hypolimnetic water used to replace the volume of daily sampling.
replaceTP_mg<-replace_TP_mgL*raw$vol_replace_L

#5. Calculate the new concentration of TP in the core following the addition of replacement water
#   Add the mass of TP (mg) in the water column (calculated in 3.) and the mass TP (mg) in the replacement water (calculated in 4.) and divide by the total volume of water (L, calc in 2.).
newTP_mgL<-(TP_mg+replaceTP_mg)/water_vol_L

#6. Bind it all together!
incubation<-cbind(raw,TP_mgL,replace_TP_mgL,water_column_m,water_vol_L,TP_mg,replaceTP_mg,newTP_mgL)

#7. Calculate daily difference in TP
# Take the TP concentration (mg/L) and subtract the "new" concentration for the previous day
incubation2<-incubation %>%
  arrange(core_unique) %>%
  group_by(core_unique) %>%
  mutate(TPDiff = TP_mgL - lag(newTP_mgL))

#8. Calculate flux (release rate) based on the change in water column TP using the following equation: 
#      Prr = (Ct - C0) × V / A / d
# Prr is the net P release (positive values) or retention (negative values) rate per unit surface area of sediment (mg P/m2/d),
# Ct is the TP concentration in the water column at time t (mg/L)
# C0 is the 'new' TP concentration the day before (mg/L)
# V is the volume of water in the water column of the core tube (L) 
# A is the planar surface area of the sediment cores (0.001735 m2)
# d is the number of days between samples Ct and C0

#NOTE - use difference in TP calcualted in #5 as Ct-C0.  mg/L.

incubation2$PrrTP<-((incubation2$TPDiff)*incubation2$water_vol_L)/0.001735/1
names(incubation2)

#9. select the relevent columns to include in the complete data table for QAQC visualizations
incubation_final<-incubation2 %>%
  select(lake,year,month,sample_id,site,hypo_temp_c,hypo_DO_mgL,hypo_DO_sat,core_id,core_unique,doy,
                         pH,temp_c,DO_mgL,DO_sat,PrrTP)

tapply(incubation_final$PrrTP,incubation_final$core_unique,mean,na.rm=TRUE)
dat<-data.frame(incubation_final)

### PART 2 - QAQC VISUALIZATIONS ---------------------------------------------------------------------------------------------------------------------------------------------
# work with "dat" dataframe
# 1. SEDIMENT P RELEASE RATES (Prr) OVER TIME (doy) --------------------------------------------------------------------------------------------------------------------------
#   a. FEBRUARY 2020 (Cores 1-3 = deep site (4). Cores 4-6 = middle site (10))
plot(NA, NA, xlim=c(40,44), ylim=c(-15,40), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
abline(h=0,col="grey50",lwd=3)
mtext(side=2, line=2.5, expression("Total P Release Rate (mg P " *m^-2*day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
legend("topleft",legend=c("Deep Site", "Middle Site"),pch=c(15,17),col=c("#fc8d62","#66c2a5"),cex=1.1,y.intersp=0.8,bty='n')
#       DEEP SITE
lines((dat[dat$core_unique=="1","PrrTP"]~(dat[dat$core_unique=="1","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="2","PrrTP"]~(dat[dat$core_unique=="2","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="3","PrrTP"]~(dat[dat$core_unique=="3","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
#       MIDDLE SITE
lines((dat[dat$core_unique=="4","PrrTP"]~(dat[dat$core_unique=="4","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="5","PrrTP"]~(dat[dat$core_unique=="5","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="6","PrrTP"]~(dat[dat$core_unique=="6","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)

#   b. APRIL 2020 (Cores 7-9 = deep site (4). Cores 10-12 = middle site (10). Cores 13-15 = shallow site (1))
plot(NA, NA, xlim=c(119,121), ylim=c(-15,40), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
abline(h=0,col="grey50",lwd=3)
mtext(side=2, line=2.5, expression("Total P Release Rate (mg P " *m^-2*day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
legend("topleft",legend=c("Deep Site", "Middle Site","Shallow Site"),pch=c(15,17,16),col=c("#fc8d62","#66c2a5","#8da0cb"),cex=1.1,y.intersp=0.8,bty='n')
#       DEEP SITE
lines((dat[dat$core_unique=="7","PrrTP"]~(dat[dat$core_unique=="7","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="8","PrrTP"]~(dat[dat$core_unique=="8","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="9","PrrTP"]~(dat[dat$core_unique=="9","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
#       MIDDLE SITE
lines((dat[dat$core_unique=="10","PrrTP"]~(dat[dat$core_unique=="10","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="11","PrrTP"]~(dat[dat$core_unique=="11","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="12","PrrTP"]~(dat[dat$core_unique=="12","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
#       SHALLOW SITE
lines((dat[dat$core_unique=="13","PrrTP"]~(dat[dat$core_unique=="13","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="14","PrrTP"]~(dat[dat$core_unique=="14","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="15","PrrTP"]~(dat[dat$core_unique=="15","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)

#   c. JUNE 2020 (Cores 16-18 = deep site (4). Cores 19-21 = middle site (10). Cores 22-24 = shallow site (1))
plot(NA, NA, xlim=c(183,185), ylim=c(-25,45), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
abline(h=0,col="grey50",lwd=3)
mtext(side=2, line=2.5, expression("Total P Release Rate (mg P " *m^-2*day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
legend("topright",legend=c("Deep Site", "Middle Site","Shallow Site"),pch=c(15,17,16),col=c("#fc8d62","#66c2a5","#8da0cb"),cex=1.1,y.intersp=0.8,bty='n')
#       DEEP SITE
lines((dat[dat$core_unique=="16","PrrTP"]~(dat[dat$core_unique=="16","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="17","PrrTP"]~(dat[dat$core_unique=="17","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="18","PrrTP"]~(dat[dat$core_unique=="18","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
#       MIDDLE SITE
lines((dat[dat$core_unique=="19","PrrTP"]~(dat[dat$core_unique=="19","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="20","PrrTP"]~(dat[dat$core_unique=="20","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="21","PrrTP"]~(dat[dat$core_unique=="21","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
#       SHALLOW SITE
lines((dat[dat$core_unique=="22","PrrTP"]~(dat[dat$core_unique=="22","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="23","PrrTP"]~(dat[dat$core_unique=="23","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="24","PrrTP"]~(dat[dat$core_unique=="24","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)

#   d. AUGUST 2020 (Cores 25-27 = deep site (4). Cores 28-30 = middle site (10). Cores 31-33 = shallow site (1))
plot(NA, NA, xlim=c(225,227), ylim=c(-100,60), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
abline(h=0,col="grey50",lwd=3)
mtext(side=2, line=2.5, expression("Total P Release Rate (mg P " *m^-2*day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
legend("bottomleft",legend=c("Deep Site", "Middle Site","Shallow Site"),pch=c(15,17,16),col=c("#fc8d62","#66c2a5","#8da0cb"),cex=1.1,y.intersp=0.8,bty='n')
#       DEEP SITE
lines((dat[dat$core_unique=="25","PrrTP"]~(dat[dat$core_unique=="25","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="26","PrrTP"]~(dat[dat$core_unique=="26","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="27","PrrTP"]~(dat[dat$core_unique=="27","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
#       MIDDLE SITE
lines((dat[dat$core_unique=="28","PrrTP"]~(dat[dat$core_unique=="28","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="29","PrrTP"]~(dat[dat$core_unique=="29","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="30","PrrTP"]~(dat[dat$core_unique=="30","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
#       SHALLOW SITE
lines((dat[dat$core_unique=="31","PrrTP"]~(dat[dat$core_unique=="31","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="32","PrrTP"]~(dat[dat$core_unique=="32","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="33","PrrTP"]~(dat[dat$core_unique=="33","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)

#   e. OCTOBER 2020 (Cores 34-36 = deep site (4). Cores 37-39 = middle site (10). Cores 40-42 = shallow site (1))
plot(NA, NA, xlim=c(300,302), ylim=c(-10,40), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
abline(h=0,col="grey50",lwd=3)
mtext(side=2, line=2.5, expression("Total P Release Rate (mg P " *m^-2*day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
legend("topleft",legend=c("Deep Site", "Middle Site","Shallow Site"),pch=c(15,17,16),col=c("#fc8d62","#66c2a5","#8da0cb"),cex=1.1,y.intersp=0.8,bty='n')
#       DEEP SITE
lines((dat[dat$core_unique=="34","PrrTP"]~(dat[dat$core_unique=="34","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="35","PrrTP"]~(dat[dat$core_unique=="35","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="36","PrrTP"]~(dat[dat$core_unique=="36","doy"])),col="#fc8d62",type='o',pch=15,cex=2,lwd=4)
#       MIDDLE SITE
lines((dat[dat$core_unique=="37","PrrTP"]~(dat[dat$core_unique=="37","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="38","PrrTP"]~(dat[dat$core_unique=="38","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="39","PrrTP"]~(dat[dat$core_unique=="39","doy"])),col="#66c2a5",type='o',pch=17,cex=2,lwd=4)
#       SHALLOW SITE
lines((dat[dat$core_unique=="40","PrrTP"]~(dat[dat$core_unique=="40","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="41","PrrTP"]~(dat[dat$core_unique=="41","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)
lines((dat[dat$core_unique=="42","PrrTP"]~(dat[dat$core_unique=="42","doy"])),col="#8da0cb",type='o',pch=17,cex=2,lwd=4)

#   f. JULY 2019 (Site 1 = Cores 43-44. Site 2 = 45-46. Site 3 = 47,48. Site 4 = 49-50. Site 5 = 51-52. Site 6 = 53-54)
plot(NA, NA, xlim=c(204,208), ylim=c(-50,20), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
abline(h=0,col="grey50",lwd=3)
mtext(side=2, line=2.5, expression("Total P Release Rate (mg P " *m^-2*day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
#legend("topleft",legend=c("Deep Site", "Middle Site","Shallow Site"),pch=c(15,17,16),col=c("#fc8d62","#66c2a5","#8da0cb"),cex=1.1,y.intersp=0.8,bty='n')
#       DEEP SITE (4)
lines((dat[dat$core_unique=="49","PrrTP"]~(dat[dat$core_unique=="49","doy"])),col="#F34226",type='o',pch=15,cex=2,lwd=4)
lines((dat[dat$core_unique=="50","PrrTP"]~(dat[dat$core_unique=="50","doy"])),col="#F34226",type='o',pch=15,cex=2,lwd=4)
#       SHALLOW SITE (1)
lines((dat[dat$core_unique=="43","PrrTP"]~(dat[dat$core_unique=="43","doy"])),col="#005867",type='o',pch=16,cex=2,lwd=4)
lines((dat[dat$core_unique=="44","PrrTP"]~(dat[dat$core_unique=="44","doy"])),col="#005867",type='o',pch=16,cex=2,lwd=4)
#       SITE (2)
lines((dat[dat$core_unique=="45","PrrTP"]~(dat[dat$core_unique=="45","doy"])),col="#F34226",type='o',pch=1,cex=2,lwd=4)
lines((dat[dat$core_unique=="46","PrrTP"]~(dat[dat$core_unique=="46","doy"])),col="#F34226",type='o',pch=1,cex=2,lwd=4)
#       SITE (3)
lines((dat[dat$core_unique=="47","PrrTP"]~(dat[dat$core_unique=="47","doy"])),col="#F34226",type='o',pch=1,cex=2,lwd=4)
lines((dat[dat$core_unique=="48","PrrTP"]~(dat[dat$core_unique=="48","doy"])),col="#F34226",type='o',pch=1,cex=2,lwd=4)
#       SITE (5)
lines((dat[dat$core_unique=="51","PrrTP"]~(dat[dat$core_unique=="51","doy"])),col="#F34226",type='o',pch=1,cex=2,lwd=4)
lines((dat[dat$core_unique=="52","PrrTP"]~(dat[dat$core_unique=="52","doy"])),col="#F34226",type='o',pch=1,cex=2,lwd=4)
#       SITE (6)
lines((dat[dat$core_unique=="53","PrrTP"]~(dat[dat$core_unique=="53","doy"])),col="#005867",type='o',pch=1,cex=2,lwd=4)
lines((dat[dat$core_unique=="54","PrrTP"]~(dat[dat$core_unique=="54","doy"])),col="#005867",type='o',pch=1,cex=2,lwd=4)


gvl<-subset(dat,month=="July")
gvl$site.f<-as.factor(gvl$site)
ggplot(gvl,aes(x=site.f, y=PrrTP,fill=hypo_DO_mgL))+
  geom_boxplot(width=0.3)+
  geom_jitter(shape=21,position=position_jitter(0.05),size=4)+
  ylim(-60,25)+theme_classic(base_size = 18)+
  scale_fill_gradient2(low="#a50026",mid="#ffffbf",high="#4575b4",midpoint=2)+
  theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+
  xlab(NULL)+ylab("TP Flux (mgP/m2/d)")+ggtitle("GVL 2019")+geom_hline(yintercept=0,linetype="dashed",color="black")
gvl_match<-subset(gvl,site.f=="1"|site.f=="4")
ggplot(gvl_match,aes(x=site.f, y=PrrTP,fill=hypo_DO_mgL))+
  geom_boxplot(width=0.3)+
  geom_jitter(shape=21,position=position_jitter(0.05),size=4)+
  ylim(-60,25)+theme_classic(base_size = 18)+
  scale_fill_gradient2(low="#a50026",mid="#ffffbf",high="#4575b4",midpoint=2)+
  theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+
  xlab(NULL)+ylab("TP Flux (mgP/m2/d)")+ggtitle("GVL 2019")+geom_hline(yintercept=0,linetype="dashed",color="black")

tapply(gvl_match$PrrTP,gvl_match$site.f,mean,na.rm=TRUE)
tapply(gvl_match$PrrTP,gvl_match$core_unique,mean,na.rm=TRUE)


# 2. COMPARE INCUBATION CONDITIONS (TEMP, DO) TO HYPOLIMNETIC CONDITIONS AT TIME OF SAMPLING ----------------------------------------------------------------------------
#   a. 1:1 PLOTS - Incubation DO ~ Hypo DO 
#          i. All data plotted together
plot(NA, NA, xlim=c(0,13), ylim=c(0,13),xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
abline(a=0,b=1,lwd=3)
mtext(side=2, line=2.5, expression("Incubation Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Hypolimnetic Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
points(dat$hypo_DO_mgL,dat$DO_mgL,cex=1.5) 
#         ii. Visualize by site
plot(NA, NA, xlim=c(0,13), ylim=c(0,13),xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
abline(a=0,b=1,lwd=3)
mtext(side=2, line=2.5, expression("Incubation Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Hypolimnetic Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
points(dat[dat$site=="4", "hypo_DO_mgL"], dat[dat$site=="4", "DO_mgL"],pch=15, cex=1.5, col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="10", "hypo_DO_mgL"], dat[dat$site=="10", "DO_mgL"],pch=17, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="1", "hypo_DO_mgL"], dat[dat$site=="1", "DO_mgL"],pch=16, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
legend("bottomright",legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),cex=1.1,y.intersp=0.8,bty='n')
#         iii. Visualize by site AND samplng event
# Preparing color palette
pttrans=0.5 #transparency value for points
ptpal=brewer.pal(7,"Dark2")
ptrgb<-col2rgb(ptpal)
ptpal<-rgb(red=ptrgb[1,],green=ptrgb[2,],blue=ptrgb[3,],alpha=pttrans*255, maxColorValue = 255) #color palette for transparent points
lnpal=brewer.pal(7,"Dark2") #color palette for solid lines

plot(NA, NA, xlim=c(0,13), ylim=c(0,13),xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
abline(a=0,b=1,lwd=3)
mtext(side=2, line=2.5, expression("Incubation Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Hypolimnetic Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
legend(7,7,legend=c("February","April","June","August","October"),pch=18,col=c(ptpal[4],ptpal[2],ptpal[5],ptpal[1],ptpal[3]),cex=1.1, y.intersp=0.8,bty='n')
legend(7,3,legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),cex=1.1,y.intersp=0.8,bty='n')
# DEEP SITE
points(dat[dat$site=="4"&dat$month=="February", "hypo_DO_mgL"], dat[dat$site=="4"&dat$month=="February", "DO_mgL"],pch=15, col=ptpal[4], cex=1.5)
points(dat[dat$site=="4"&dat$month=="April", "hypo_DO_mgL"], dat[dat$site=="4"&dat$month=="April", "DO_mgL"],pch=15, col=ptpal[2], cex=1.5)
points(dat[dat$site=="4"&dat$month=="June", "hypo_DO_mgL"], dat[dat$site=="4"&dat$month=="June", "DO_mgL"],pch=15, col=ptpal[5], cex=1.5)
points(dat[dat$site=="4"&dat$month=="August", "hypo_DO_mgL"], dat[dat$site=="4"&dat$month=="August", "DO_mgL"],pch=15, col=ptpal[1], cex=1.55)
points(dat[dat$site=="4"&dat$month=="October", "hypo_DO_mgL"], dat[dat$site=="4"&dat$month=="October", "DO_mgL"],pch=15, col=ptpal[3], cex=1.55)
# MIDDLE SITE
points(dat[dat$site=="10"&dat$month=="February", "hypo_DO_mgL"], dat[dat$site=="10"&dat$month=="February", "DO_mgL"],pch=17, col=ptpal[4], cex=1.5)
points(dat[dat$site=="10"&dat$month=="April", "hypo_DO_mgL"], dat[dat$site=="10"&dat$month=="April", "DO_mgL"],pch=17, col=ptpal[2], cex=1.5)
points(dat[dat$site=="10"&dat$month=="June", "hypo_DO_mgL"], dat[dat$site=="10"&dat$month=="June", "DO_mgL"],pch=17, col=ptpal[5], cex=1.5)
points(dat[dat$site=="10"&dat$month=="August", "hypo_DO_mgL"], dat[dat$site=="10"&dat$month=="August", "DO_mgL"],pch=17, col=ptpal[1], cex=1.55)
points(dat[dat$site=="10"&dat$month=="October", "hypo_DO_mgL"], dat[dat$site=="10"&dat$month=="October", "DO_mgL"],pch=17, col=ptpal[3], cex=1.55)
# SHALLOW SITE
points(dat[dat$site=="1"&dat$month=="April", "hypo_DO_mgL"], dat[dat$site=="1"&dat$month=="April", "DO_mgL"],pch=16, col=ptpal[2], cex=1.5)
points(dat[dat$site=="1"&dat$month=="June", "hypo_DO_mgL"], dat[dat$site=="1"&dat$month=="June", "DO_mgL"],pch=16, col=ptpal[5], cex=1.5)
points(dat[dat$site=="1"&dat$month=="August", "hypo_DO_mgL"], dat[dat$site=="1"&dat$month=="August", "DO_mgL"],pch=16, col=ptpal[1], cex=1.55)
points(dat[dat$site=="1"&dat$month=="October", "hypo_DO_mgL"], dat[dat$site=="1"&dat$month=="October", "DO_mgL"],pch=16, col=ptpal[3], cex=1.55)

#   b.  1:1 PLOTS - Incubation Temp ~ Hypo Temp
plot(NA, NA, xlim=c(0,30), ylim=c(0,30),xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
abline(a=0,b=1,lwd=3)
mtext(side=2, line=2.5, expression("Incubation Temperature (C)"), cex=1.2)
mtext(side=1, line=2.5, expression("Hypolimnetic Temperature (C)"), cex=1.2)
points(dat$hypo_temp_c,dat$temp_c,cex=1.5) 

#   c.Prr as a function of core DO, temperature, and pH
#          i. All data plotted together
#   DISSOLVED OXYGEN
plot(NA, NA, xlim=c(0,15), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
points(dat$DO_mgL,dat$PrrTP,cex=1.5) 
#   TEMPERATURE
plot(NA, NA, xlim=c(0,30), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Temperature (C)"), cex=1.2)
points(dat[dat$year=="2020", "temp_c"], dat[dat$year=="2020", "PrrTP"],cex=1.5) 
#   pH
plot(NA, NA, xlim=c(6.5,9.5), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("pH"), cex=1.2)
points(dat$pH,dat$PrrTP,cex=1.5)

#         ii. Visualize by site
#   DISSOLVED OXYGEN
plot(NA, NA, xlim=c(0,15), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
legend("bottomright",legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),cex=1.1,y.intersp=0.8,bty='n')
points(dat[dat$site=="4", "DO_mgL"], dat[dat$site=="4", "PrrTP"],pch=15, cex=1.5, col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="10", "DO_mgL"], dat[dat$site=="10", "PrrTP"],pch=17, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="1", "DO_mgL"], dat[dat$site=="1", "PrrTP"],pch=16, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
#   TEMPERATURE
plot(NA, NA, xlim=c(0,30), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Temperature (C)"), cex=1.2)
legend("bottomleft",legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),cex=1.1,y.intersp=0.8,bty='n')
points(dat[dat$site=="4"&dat$year=="2020", "temp_c"], dat[dat$site=="4"&dat$year=="2020", "PrrTP"],pch=15, cex=1.5, col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="10"&dat$year=="2020", "temp_c"], dat[dat$site=="10"&dat$year=="2020", "PrrTP"],pch=17, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="1"&dat$year=="2020", "temp_c"], dat[dat$site=="1"&dat$year=="2020", "PrrTP"],pch=16, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
#   pH
plot(NA, NA, xlim=c(6.5,9.5), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("pH"), cex=1.2)
legend("bottomright",legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),cex=1.1,y.intersp=0.8,bty='n')
points(dat[dat$site=="4", "pH"], dat[dat$site=="4", "PrrTP"],pch=15, cex=1.5, col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="10", "pH"], dat[dat$site=="10", "PrrTP"],pch=17, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))
points(dat[dat$site=="1", "pH"], dat[dat$site=="1", "PrrTP"],pch=16, cex=1.5,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5))

#         iii. Visualize by site AND samplng event
# Preparing color palette
pttrans=0.5 #transparency value for points
ptpal=brewer.pal(7,"Dark2")
ptrgb<-col2rgb(ptpal)
ptpal<-rgb(red=ptrgb[1,],green=ptrgb[2,],blue=ptrgb[3,],alpha=pttrans*255, maxColorValue = 255) #color palette for transparent points
lnpal=brewer.pal(7,"Dark2") #color palette for solid lines
#   DISSOLVED OXYGEN
plot(NA, NA, xlim=c(0,15), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Dissolved Oxygen (mg " *L^-1*")"), cex=1.2)
legend(8,-20,legend=c("February","April","June","August","October"),pch=18,col=c(ptpal[4],ptpal[2],ptpal[5],ptpal[1],ptpal[3]),cex=1.1, y.intersp=0.8,bty='n')
legend(8,-70,legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),cex=1.1,y.intersp=0.8,bty='n')

points(dat[dat$site=="4"&dat$month=="February", "DO_mgL"], dat[dat$site=="4"&dat$month=="February", "PrrTP"],pch=15, col=ptpal[4], cex=1.5)
points(dat[dat$site=="4"&dat$month=="April", "DO_mgL"], dat[dat$site=="4"&dat$month=="April", "PrrTP"],pch=15, col=ptpal[2], cex=1.5)
points(dat[dat$site=="4"&dat$month=="June", "DO_mgL"], dat[dat$site=="4"&dat$month=="June", "PrrTP"],pch=15, col=ptpal[5], cex=1.5)
points(dat[dat$site=="4"&dat$month=="August", "DO_mgL"], dat[dat$site=="4"&dat$month=="August", "PrrTP"],pch=15, col=ptpal[1], cex=1.55)
points(dat[dat$site=="4"&dat$month=="October", "DO_mgL"], dat[dat$site=="4"&dat$month=="October", "PrrTP"],pch=15, col=ptpal[3], cex=1.55)

points(dat[dat$site=="10"&dat$month=="February", "DO_mgL"], dat[dat$site=="10"&dat$month=="February", "PrrTP"],pch=17, col=ptpal[4], cex=1.5)
points(dat[dat$site=="10"&dat$month=="April", "DO_mgL"], dat[dat$site=="10"&dat$month=="April", "PrrTP"],pch=17, col=ptpal[2], cex=1.5)
points(dat[dat$site=="10"&dat$month=="June", "DO_mgL"], dat[dat$site=="10"&dat$month=="June", "PrrTP"],pch=17, col=ptpal[5], cex=1.5)
points(dat[dat$site=="10"&dat$month=="August", "DO_mgL"], dat[dat$site=="10"&dat$month=="August", "PrrTP"],pch=17, col=ptpal[1], cex=1.55)
points(dat[dat$site=="10"&dat$month=="October", "DO_mgL"], dat[dat$site=="10"&dat$month=="October", "PrrTP"],pch=17, col=ptpal[3], cex=1.55)

points(dat[dat$site=="1"&dat$month=="February", "DO_mgL"], dat[dat$site=="1"&dat$month=="February", "PrrTP"],pch=16, col=ptpal[4], cex=1.5)
points(dat[dat$site=="1"&dat$month=="April", "DO_mgL"], dat[dat$site=="1"&dat$month=="April", "PrrTP"],pch=16, col=ptpal[2], cex=1.5)
points(dat[dat$site=="1"&dat$month=="June", "DO_mgL"], dat[dat$site=="1"&dat$month=="June", "PrrTP"],pch=16, col=ptpal[5], cex=1.5)
points(dat[dat$site=="1"&dat$month=="August", "DO_mgL"], dat[dat$site=="1"&dat$month=="August", "PrrTP"],pch=16, col=ptpal[1], cex=1.55)
points(dat[dat$site=="1"&dat$month=="October", "DO_mgL"], dat[dat$site=="1"&dat$month=="October", "PrrTP"],pch=16, col=ptpal[3], cex=1.55)

#   TEMPERATURE
plot(NA, NA, xlim=c(0,30), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Temperature (C)"), cex=1.2)
legend(4,-20,legend=c("February","April","June","August","October"),pch=18,col=c(ptpal[4],ptpal[2],ptpal[5],ptpal[1],ptpal[3]),cex=1.1, y.intersp=0.8,bty='n')
legend(4,-70,legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),cex=1.1,y.intersp=0.8,bty='n')

points(dat[dat$site=="4"&dat$month=="February", "temp_c"], dat[dat$site=="4"&dat$month=="February", "PrrTP"],pch=15, col=ptpal[4], cex=1.5)
points(dat[dat$site=="4"&dat$month=="April", "temp_c"], dat[dat$site=="4"&dat$month=="April", "PrrTP"],pch=15, col=ptpal[2], cex=1.5)
points(dat[dat$site=="4"&dat$month=="June", "temp_c"], dat[dat$site=="4"&dat$month=="June", "PrrTP"],pch=15, col=ptpal[5], cex=1.5)
points(dat[dat$site=="4"&dat$month=="August", "temp_c"], dat[dat$site=="4"&dat$month=="August", "PrrTP"],pch=15, col=ptpal[1], cex=1.5)
points(dat[dat$site=="4"&dat$month=="October", "temp_c"], dat[dat$site=="4"&dat$month=="October", "PrrTP"],pch=15, col=ptpal[3], cex=1.5)

points(dat[dat$site=="10"&dat$month=="February", "temp_c"], dat[dat$site=="10"&dat$month=="February", "PrrTP"],pch=17, col=ptpal[4], cex=1.5)
points(dat[dat$site=="10"&dat$month=="April", "temp_c"], dat[dat$site=="10"&dat$month=="April", "PrrTP"],pch=17, col=ptpal[2], cex=1.5)
points(dat[dat$site=="10"&dat$month=="June", "temp_c"], dat[dat$site=="10"&dat$month=="June", "PrrTP"],pch=17, col=ptpal[5], cex=1.5)
points(dat[dat$site=="10"&dat$month=="August", "temp_c"], dat[dat$site=="10"&dat$month=="August", "PrrTP"],pch=17, col=ptpal[1], cex=1.5)
points(dat[dat$site=="10"&dat$month=="October", "temp_c"], dat[dat$site=="10"&dat$month=="October", "PrrTP"],pch=17, col=ptpal[3], cex=1.5)

points(dat[dat$site=="1"&dat$month=="April", "temp_c"], dat[dat$site=="1"&dat$month=="April", "PrrTP"],pch=16, col=ptpal[2], cex=1.5)
points(dat[dat$site=="1"&dat$month=="June", "temp_c"], dat[dat$site=="1"&dat$month=="June", "PrrTP"],pch=16, col=ptpal[5], cex=1.5)
points(dat[dat$site=="1"&dat$month=="August", "temp_c"], dat[dat$site=="1"&dat$month=="August", "PrrTP"],pch=16, col=ptpal[1], cex=1.5)
points(dat[dat$site=="1"&dat$month=="October", "temp_c"], dat[dat$site=="1"&dat$month=="October", "PrrTP"],pch=16, col=ptpal[3], cex=1.5)

#   pH
plot(NA, NA, xlim=c(6.5,9.5), ylim=c(-100,60), 
     xlab="", ylab="", cex.lab=1.2, cex.axis=1.1, las=1)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey95",border=NA)
mtext(side=2, line=2.5, expression(P~Release~Rate~"("*mg~P~m^-2~day^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("pH"), cex=1.2)
legend(8.2,-20,legend=c("February","April","June","August","October"),pch=18,col=c(ptpal[4],ptpal[2],ptpal[5],ptpal[1],ptpal[3]),cex=1.1, y.intersp=0.8,bty='n')
legend(8.2,-65,legend=c("Deep Site", "Middle Site", "Shallow Site"),pch=c(15,17,16),cex=1.1,y.intersp=0.8,bty='n')
# 
points(dat[dat$site=="4"&dat$month=="February", "pH"], dat[dat$site=="4"&dat$month=="February", "PrrTP"],pch=15, col=ptpal[4], cex=1.5)
points(dat[dat$site=="4"&dat$month=="April", "pH"], dat[dat$site=="4"&dat$month=="April", "PrrTP"],pch=15, col=ptpal[2], cex=1.5)
points(dat[dat$site=="4"&dat$month=="June", "pH"], dat[dat$site=="4"&dat$month=="June", "PrrTP"],pch=15, col=ptpal[5], cex=1.5)
points(dat[dat$site=="4"&dat$month=="August", "pH"], dat[dat$site=="4"&dat$month=="August", "PrrTP"],pch=15, col=ptpal[1], cex=1.5)
points(dat[dat$site=="4"&dat$month=="October", "pH"], dat[dat$site=="4"&dat$month=="October", "PrrTP"],pch=15, col=ptpal[3], cex=1.5)

points(dat[dat$site=="10"&dat$month=="February", "pH"], dat[dat$site=="10"&dat$month=="February", "PrrTP"],pch=17, col=ptpal[4], cex=1.5)
points(dat[dat$site=="10"&dat$month=="April", "pH"], dat[dat$site=="10"&dat$month=="April", "PrrTP"],pch=17, col=ptpal[2], cex=1.5)
points(dat[dat$site=="10"&dat$month=="June", "pH"], dat[dat$site=="10"&dat$month=="June", "PrrTP"],pch=17, col=ptpal[5], cex=1.5)
points(dat[dat$site=="10"&dat$month=="August", "pH"], dat[dat$site=="10"&dat$month=="August", "PrrTP"],pch=17, col=ptpal[1], cex=1.5)
points(dat[dat$site=="10"&dat$month=="October", "pH"], dat[dat$site=="10"&dat$month=="October", "PrrTP"],pch=17, col=ptpal[3], cex=1.5)

points(dat[dat$site=="1"&dat$month=="April", "pH"], dat[dat$site=="1"&dat$month=="April", "PrrTP"],pch=16, col=ptpal[2], cex=1.5)
points(dat[dat$site=="1"&dat$month=="June", "pH"], dat[dat$site=="1"&dat$month=="June", "PrrTP"],pch=16, col=ptpal[5], cex=1.5)
points(dat[dat$site=="1"&dat$month=="August", "pH"], dat[dat$site=="1"&dat$month=="August", "PrrTP"],pch=16, col=ptpal[1], cex=1.5)
points(dat[dat$site=="1"&dat$month=="October", "pH"], dat[dat$site=="1"&dat$month=="October", "PrrTP"],pch=16, col=ptpal[3], cex=1.5)

# 3. Investigate low Prr for the shallow site, August sampling event 
# Note: low outliers on last day of incubation for cores 32, 33
# Temperature
plot(NA, NA, xlim=c(224,227), ylim=c(23,26), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
mtext(side=2, line=2.5, expression("Temperature (C)"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
lines((dat[dat$core_unique=="32","temp_c"]~(dat[dat$core_unique=="32","doy"])),col="red",type='o',pch=16,cex=2,lwd=4)
lines((dat[dat$core_unique=="33","temp_c"]~(dat[dat$core_unique=="33","doy"])),col="orange",type='o',pch=16,cex=2,lwd=4)
lines((dat[dat$core_unique=="31","temp_c"]~(dat[dat$core_unique=="31","doy"])),type='o',pch=16,cex=2,lwd=4)
# Core without a low outlier (31) follows the same pattern as the cores with low outliers on the last day of the incubation
# DO
plot(NA, NA, xlim=c(224,227), ylim=c(3,8), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
mtext(side=2, line=2.5, expression("Dissolved Oxygen (mg" *L^-1*")"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
lines((dat[dat$core_unique=="32","DO_mgL"]~(dat[dat$core_unique=="32","doy"])),col="red",type='o',pch=16,cex=2,lwd=4)
lines((dat[dat$core_unique=="33","DO_mgL"]~(dat[dat$core_unique=="33","doy"])),col="orange",type='o',pch=16,cex=2,lwd=4)
lines((dat[dat$core_unique=="31","DO_mgL"]~(dat[dat$core_unique=="31","doy"])),type='o',pch=16,cex=2,lwd=4)
# Core without a low outlier (31) follows the same pattern as core 33... no clear pattern
# pH
plot(NA, NA, xlim=c(224,227), ylim=c(7,8.5), xlab="", ylab="", cex.lab=1.2, cex.axis=1.1,las=1)
mtext(side=2, line=2.5, expression("pH"), cex=1.2)
mtext(side=1, line=2.5, expression("Incubation DOY"), cex=1.2)
lines((dat[dat$core_unique=="32","pH"]~(dat[dat$core_unique=="32","doy"])),col="red",type='o',pch=16,cex=2,lwd=4)
lines((dat[dat$core_unique=="33","pH"]~(dat[dat$core_unique=="33","doy"])),col="orange",type='o',pch=16,cex=2,lwd=4)
lines((dat[dat$core_unique=="31","pH"]~(dat[dat$core_unique=="31","doy"])),type='o',pch=16,cex=2,lwd=4)
# Core without a low outlier (31) follows the same pattern as the cores with low outliers on the last day of the incubation

### PART 3 - FINAL DATA TABLE -----------------------------------------------------------------------------------------------------------------------------------------------
# Exclude incubation doy 40, 41 for cores 1, 2, 3, (February incubation, site 4) - core DO too high above ambient conditions at time of sampling
# Exclude incubation doy 227 for cores 32, 33 (August incubation, site 1) - low values likely due to laboratory error
# Exclude incubation doy 205 for core 43 (July 2019 incubation, site 1) - low value likely due to laboratory error
# Select only sample site that were studied in 2020 or 2020 and 2019 (so exclude 2, 3, 5, 6)

dat2<-dat %>%
  filter(!sample_id %in% c("C20041040401","C20041041401","C20041040402","C20041041402","C20041040403","C20041041403")) %>%
  filter(!sample_id %in% c("C20041227108","C20041227109")) %>% 
  filter(!sample_id %in% c("C1904120501o1")) %>% 
  filter(!site %in% c("2","3","5","6"))

# Write final data table
write.csv(dat2,file="GVL_Incubation_Tidy.csv")

# Write a SUMMARY data table - takes the average of temporal replicates (i.e. over the course of the incubation) for each core
sum<-dat2 %>%
  group_by(core_unique) %>%
  mutate(mean_hypo_temp = mean(hypo_temp_c),mean_hypo_DO = mean(hypo_DO_mgL),mean_hypo_DOsat = mean(hypo_DO_sat),
         mean_temp = mean(temp_c), mean_DO = mean(DO_mgL), mean_DOsat = mean(DO_sat), mean_pH = mean(pH),
         mean_Prr = mean(PrrTP,na.rm=T)) %>%
  select(lake,year,month,site,core_unique,mean_hypo_temp,mean_hypo_DO,mean_hypo_DOsat,mean_temp,mean_DO,mean_DOsat,mean_pH,mean_Prr) %>%
  slice(n=1)

write.csv(sum,file="GVL_Incubation_Sum.csv")

# Visualize how excluding outliers and summarizing temporal replicates affects distribution of Prr
plot(density(dat$PrrTP,na.rm = TRUE))
plot(density(dat2$PrrTP,na.rm= TRUE))
plot(density(test$mean_Prr))
plot(density(test2$mean_Prr))



