# Project: Spatiotemporal Variation in Internal P Loading in a Hypereutrophic Reservoir (GVL) - Ecosystems
# Last Modified 12 October 2021
# Contributers: Ellen Albright, Grace Wilkinson
# Description: The following script details all analyses and figures for our manuscript: 
#     "Sediment phosphorus speciation controls hot spots and hot moments of internal loading in a temperate reservoir"

# Data citation: 

# The code can be run with the following datasets:
# "GVL_Incubation_Sum.csv" - Sediment Core Incubation, Summarized Data, Green Valley Lake (2019-2020)
#       Final dataset of calculated sediment P release rates, averaged by sediment core over the course of each incubation
# "GVL_WaterChem_RAW.csv" - Water Chemistry, Raw Data, Green Valley Lake (2014-2020)
#       Water column chemical parameters (TP, SRP, NOx, TN, TSS/ISS/VSS, Microcystin). Spatial and temporal resolution varies considerably across variables. 
# "GVL_SedimentP_Sum.csv" - Sediment P Speciation, Summarized Data, Green Valley Lake (2020)
#       Final dataset of calculated concentrations of various sediment P species, averaged over three laboratory replicates
# "GVL_SedTP_Sum.csv" - Sediment Total Phosphorus, Summarized Data, Green Valley Lake (2020)
#       Final dataset of calculated concentrations of sediment total P, averaged over three laboratory replicates
# "GVL_bathymetry.shp" - Green Valley Lake Bathymetry, Spatial Data
#       Shapefile of lakebed bathymetry every 2 feet, surveyed by Iowa DNR Fisheries Bureau
# "GVL_sites.csv" - Coordinates of sampling sites
# "GVL_SedimentAreaCALC.csv" - Necessary morphometric information to calculate areas of lakebed associated with each sampling site


# Clear environment, set working directory - UPDATE AS NEEDED --------------------------------------------------------------------------------------------------------------
rm(list=ls())
setwd('C:/Users/Ellen/Desktop/Box Sync/Albright DISSERTATION/Albright_Chapter_2/DATA/DATA_PUBLISH')

# REQUIRED PACKAGES FOR CODE - install as needed
# install.packages("tidyverse")
# install.packages("RColorBrewer")
# install.packages("gridExtra")
# install.packages('fGarch')
# install.pacakges("robCompositions")
# install.packages("vegan")
# install.packages("sf")
# install.packages("ggspatial")
# install.packages("spData")
# install.packages("cowplot")

library(tidyverse) #majority of figures made in ggplot2, will need this package to recreate figures
library(RColorBrewer) #used for some figures
library(gridExtra) #used for some figures
library(fGarch) # used for analysis of distribution skew
library(robCompositions) #used for compositional data analysis (CoDA) of sediment P chemistry
library(vegan) #used for CoDA of sediment P chemistry (PCA on covariance matrix)
library(sf) #needed for reading in shapefiles, other spatial analyses
library(ggspatial) #needed to make site map 
library(spData) #used US map available in this package to make inset map
library(cowplot) #needed to make site map with inset map


### OVERVIEW OF SCRIPT ORGANIZATION ### ------------------------------------------------------------------------------------------------------------------------------------
# PART 1: Spatiotemporal variation in water column chemistry and sediment P flux rates
#      1a. Summarize mean flux rate and standard error of the mean across replicate cores for each site and event (Table S2)
#      1b. Plot timeseries of sediment P flux rates and hypolimnetic TP and SRP (Figure 2)
#      1c. Explore epi- and hypolimnetic P patterns over time (2014, 2015, 2019, 2020) (Figure S1)
# PART 2: Explore hot spots and hot moments of sediment P release
#      2a. Visualize the P flux rate distribution and identify outliers (Figure S2)
#      2b. Visualize the P flux distribution by site and sampling event (Figure S3)
#      2c. Quantify the distribution skew, leave-one-out analysis (Table S3)
# PART 3: Sediment P flux rates and dissolved oxygen status at the sediment-water interface (Figure 3)
# PART 4: Sediment P Pools
#      4a. Mean sediment P composition by sampling site (Figure S4)
#      4b. Timeseries of mobile sediment P fractions by sampling site (Figure 4)
#      4c. CoDA of sediment P composition to understand variation over space and time (Figure 4)
# PART 5: Scaling P fluxes across the lakebed to estimate load; Site map (Figure 5, Table 2, Figure 1)
#      5a. Processing shapefile of bathymetric data
#      5b. Site map (Figure 1)
#      5c. P load estimation (Figure 5, Table 2)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### PART 1 - SPATIOTEMPORAL VARIATION IN WATER COLUMN CHEMISTRY AND SEDIMENT P FLUX RATES (TABLE S2; FIGURES 2, S1) --------------------------------------------------------
# READ IN DATA: Summarized sediment core incubation data (summarizes mean daily, aerial total P release rate for each sediment core over the course of the incubation)
sum<-read.csv("GVL_Incubation_Sum.csv")
# READ IN DATA: water chemistry data to plot hypolimnetic TP and SRP 
waterchem <- read.csv("GVL_WaterChem_RAW.csv")

### Part 1A - Summarize mean flux rate and standard error of the mean across replicate cores for each site and event (Table S2) ------------------

# Investigate "sum" dataframe and add necessary columns for visualizations
str(sum) # site number is read as integer, need to make a new column with sampling site as a factor
sum$site.f<-as.factor(sum$site) # make a new column with sampling site as a factor, not an integer
sum$n_cores[sum$year=="2020"] <- 3 # add a column that indicates the number of replicate cores used in the experiment (3 for 2020)...
sum$n_cores[sum$year=="2019"] <- 2 #...and 2 for 2019. Will be used to calculate the standard error of the mean
sum$doy[sum$month=="February"] <- 39 # add day of year when samples were collected to make timeseries
sum$doy[sum$month=="April"] <- 117
sum$doy[sum$month=="June"] <- 181
sum$doy[sum$month=="August"] <- 223
sum$doy[sum$month=="October"] <- 298
sum$doy[sum$month=="July"] <- 203

# Further summarize the mean and standard error of the mean across replicate cores for each site and event (Values for TABLE S2)
sum2<-sum %>% 
  mutate(month=factor(month,levels=c("February","April","June","July","August","October"))) %>% 
  mutate(site.f=factor(site.f,levels=c("1","10","4"))) %>% 
  group_by(month, site.f) %>% 
  mutate(sem_Prr=(sd(mean_Prr))/sqrt(n_cores),average_Prr=mean(mean_Prr),hypo_DO=mean(mean_hypo_DO)) %>% 
  select(year,month,doy,site.f,hypo_DO,average_Prr,sem_Prr) %>% 
  slice(n=1)

### Part 1B - Plot timeseries of sediment P flux rates and hypolimnetic TP and SRP (Figure 2) ----------------------------------------------------
# Subset both sum2 and waterchem dataframes to work with years/sites individually
deep20 <- dplyr::filter(sum2, year=="2020" & site.f=="4") # P flux rates, deep site, 2020
deep20wq <-dplyr::filter(waterchem, year=="2020"&siteID=="4"&depthID=="4") # water chem, deep site, hypolimnetic, 2020
mid20 <- dplyr::filter(sum2, year=="2020" & site.f=="10") # P flux rates, intermediate depth site, 2020
mid20wq <-dplyr::filter(waterchem, year=="2020"&siteID=="10"&depthID=="4") # water chem, intermediate depth site, hypolimnetic, 2020
shallow20 <- dplyr::filter(sum2, year=="2020" & site.f=="1") # P flux rates, shallow site, 2020
shallow20wq <-dplyr::filter(waterchem, year=="2020"&siteID=="1"&depthID=="4") # water chem, shallow site, hypolimnetic, 2020


# FIGURE 2 - Sediment P Flux Rates and Hypolimnetic P (2020) 
colors<-c("SRP_ugL"="gray30","TP_ugL"="gray15")

p1<-ggplot(data=shallow20) +
  geom_hline(yintercept=0,color="gray40",linetype = 'dotted',size=1) +
  geom_line(aes(x=doy,y=average_Prr),size=1.2,color="#1B9E77")+
  geom_point(aes(x=doy,y=average_Prr),color="#1B9E77",size=4,shape=19,stroke=1.4)+
  geom_errorbar(aes(x=doy,ymin=average_Prr-sem_Prr,ymax=average_Prr+sem_Prr),color="#1B9E77",width=0, position=position_dodge(0.05))+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlim(30,305) + ylim(-15,25) + theme(legend.position = "none") +
  xlab(" ") + ylab(bquote('P Flux Rate (mg ' *~m^-2~day^-1*')'))+theme(axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  ggtitle("Shallow")+theme(plot.title=element_text(face="bold",size=10))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,0),"cm"))
p2 <-ggplot(data=shallow20wq) +
  geom_line(aes(x=doy,y=SRP_ugL,color="SRP_ugL"),size=1.2) +
  geom_point(aes(x=doy,y=SRP_ugL),color="gray30",size=4,shape=19)+
  geom_line(aes(x=doy,y=TP_ugL,color="TP_ugL"),size=1.2) +
  geom_point(aes(x=doy,y=TP_ugL),color="gray15",size=4,shape=19)+
  scale_color_manual(values=colors,
                     name=" ",
                     breaks=c("SRP_ugL","TP_ugL"),
                     labels=c("SRP","TP"))+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  theme(legend.title=element_blank(), legend.position = c(0.25,0.8),legend.text = element_text(size=8),legend.key.size=unit(0.3,"cm"))+
  xlim(30,305) + ylim(0,400) + 
  xlab("Day of Year") + ylab(bquote('Hypolimnetic P (µg ' *L^-1*')'))+
  theme(axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,0),"cm"))
p3<-ggplot(data=mid20) +
  geom_hline(yintercept=0,color="gray40",linetype = 'dotted',size=1) +
  geom_line(aes(x=doy,y=average_Prr),size=1.2,color="#D95F02")+
  geom_point(aes(x=doy,y=average_Prr),color="#D95F02",size=4,shape=17,stroke=1.4)+
  geom_errorbar(aes(x=doy,ymin=average_Prr-sem_Prr,ymax=average_Prr+sem_Prr),color="#D95F02",width=0, position=position_dodge(0.05))+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlim(30,305) + ylim(-15,25) + theme(legend.position = "none") + 
  xlab(" ") + ylab(bquote(' '))+theme(axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  ggtitle("Intermediate Depth")+theme(plot.title=element_text(face="bold",size=10))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,-0.1),"cm"))
p4 <-ggplot(data=mid20wq) +
  geom_line(aes(x=doy,y=SRP_ugL,color="SRP_ugL"),size=1.2) +
  geom_point(aes(x=doy,y=SRP_ugL),color="gray30",size=4,shape=17)+
  geom_line(aes(x=doy,y=TP_ugL,color="TP_ugL"),size=1.2) +
  geom_point(aes(x=doy,y=TP_ugL),color="gray15",size=4,shape=17)+
  scale_color_manual(values=colors,
                     name=" ",
                     breaks=c("SRP_ugL","TP_ugL"),
                     labels=c("SRP","TP"))+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  theme(legend.title=element_blank(), legend.position = c(0.25,0.8),legend.text = element_text(size=8),legend.key.size=unit(0.3,"cm"))+
  xlim(30,305) + ylim(0,400) + 
  xlab("Day of Year") + ylab(bquote(' '))+theme(axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,-0.1),"cm"))
p5<-ggplot(data=deep20) +
  geom_hline(yintercept=0,color="gray40",linetype = 'dotted',size=1) +
  geom_line(aes(x=doy,y=average_Prr),size=1.2,color="#E7298A")+
  geom_point(aes(x=doy,y=average_Prr),color="#E7298A",size=4,shape=15,stroke=1.4)+
  geom_errorbar(aes(x=doy,ymin=average_Prr-sem_Prr,ymax=average_Prr+sem_Prr),color="#E7298A",width=0, position=position_dodge(0.05))+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlim(30,305) + ylim(-15,25) + theme(legend.position = "none") +
  xlab("") + ylab(bquote(' '))+theme(axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  ggtitle("Deep")+theme(plot.title=element_text(face="bold",size=10))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,-0.1),"cm"))
p6 <-ggplot(data=deep20wq) +
  geom_line(aes(x=doy,y=SRP_ugL,color="SRP_ugL"),size=1.2) +
  geom_point(aes(x=doy,y=SRP_ugL),color="gray30",size=4,shape=15)+
  geom_line(aes(x=doy,y=TP_ugL,color="TP_ugL"),size=1.2) +
  geom_point(aes(x=doy,y=TP_ugL),color="gray15",size=4,shape=15)+
  scale_color_manual(values=colors,
                     name=" ",
                     breaks=c("SRP_ugL","TP_ugL"),
                     labels=c("SRP","TP"))+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  theme(legend.title=element_blank(), legend.position = c(0.25,0.8),legend.text = element_text(size=8),legend.key.size=unit(0.3,"cm"))+
  xlim(30,305) + ylim(0,400) + xlab("Day of Year") + ylab(bquote(' '))+theme(axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,-0.1),"cm"))

windows(height=3.3,width=6.5)  
plotfull<-grid.arrange(p1,p3,p5,p2,p4,p6,nrow=2) #need to name the grid.arrange to use in ggsave
ggsave("GVL_2020_Prr_HypoP.png", plotfull, width=6.5, height=3.3, units="in", dpi=300) #FIGURE 2

### Part 1C - Explore epi- and hypolimnetic P patterns over time (2014, 2015, 2019, 2020) (Figure S1) --------------------------------------------
# Epilimnetic TP & SRP
wq_surface_deepsite<-dplyr::filter(waterchem, siteID=="4",depthID=="1")
wq_surface_deepsite$year.f<-as.factor(wq_surface_deepsite$year)

wq1<-ggplot(data=wq_surface_deepsite, aes(x=doy, y=TP_ugL, group=year.f))+
  geom_line(aes(color=year.f),size=1.5)+geom_point(aes(color=year.f),shape=15,size=2)+
  scale_color_manual(values=c("#9ecae1","#4292c6","#08519c","#e7298a"))+theme_classic()+
  xlab("Day of Year") + ylab(bquote('Epilimnetic Total P (µg ' *L^-1*')'))+theme(axis.text=element_text(color="black",size=12),axis.title=element_text(size=14))+
  theme(legend.title=element_blank(), legend.position = c(0.3,0.8),legend.text = element_text(size=14))
wq2<-ggplot(data=wq_surface_deepsite, aes(x=doy, y=SRP_ugL, group=year.f))+
  geom_line(aes(color=year.f),size=1.5)+geom_point(aes(color=year.f),shape=15,size=2)+
  scale_color_manual(values=c("#9ecae1","#4292c6","#08519c","#e7298a"))+theme_classic()+
  xlab("Day of Year") + ylab(bquote('Epilimnetic Soluable Reactive P (µg ' *L^-1*')'))+theme(axis.text=element_text(color="black",size=12),axis.title=element_text(size=14))+
  theme(legend.title=element_blank(), legend.position = c(0.3,0.8),legend.text = element_text(size=14))

# Hypolimnetic TP & SRP
wq_bottom_deepsite<-dplyr::filter(waterchem, siteID=="4",depthID=="4")
wq_bottom_deepsite$year.f<-as.factor(wq_bottom_deepsite$year)

wq3<-ggplot(data=wq_bottom_deepsite, aes(x=doy, y=TP_ugL, group=year.f))+
  geom_line(aes(color=year.f),size=1.5)+geom_point(aes(color=year.f),shape=15,size=2)+
  scale_color_manual(values=c("#9ecae1","#4292c6","#08519c","#e7298a"))+theme_classic()+
  xlab("Day of Year") + ylab(bquote('Hypolimnetic Total P (µg ' *L^-1*')'))+theme(axis.text=element_text(color="black",size=12),axis.title=element_text(size=14))+
  theme(legend.title=element_blank(), legend.position = c(0.3,0.8),legend.text = element_text(size=14))
wq4<-ggplot(data=wq_bottom_deepsite, aes(x=doy, y=SRP_ugL, group=year.f))+
  geom_line(aes(color=year.f),size=1.5)+geom_point(aes(color=year.f),shape=15,size=2)+
  scale_color_manual(values=c("#9ecae1","#4292c6","#08519c","#e7298a"))+theme_classic()+
  xlab("Day of Year") + ylab(bquote('Hypolimnetic Soluable Reactive P (µg ' *L^-1*')'))+theme(axis.text=element_text(color="black",size=12),axis.title=element_text(size=14))+
  theme(legend.title=element_blank(), legend.position = c(0.3,0.8),legend.text = element_text(size=14))

plotfull_wq<-grid.arrange(wq1,wq2,wq3,wq4,nrow=2) #need to name the grid.arrange to use in ggsave
ggsave("GVL_FigS1.png",plotfull_wq,width=8,height=8,units="in",dpi=300)

### PART 2 - EXPLORE HOT SPOTS AND HOT MOMENTS OF SEDIMENT P RELEASE (Figures S2, S3; Table S3) ---------------------------------------------------------------------------
### Part 2A - Visualize the P flux rate distribution and identify outliers (Figure S2) -----------------------------------------------------------

# Subset "sum" data frame for 2020 sampling study year (data=sum, subset for 2020. all sites, all cores, whole year)
sum20 <- dplyr::filter(sum,year=="2020")

# Visualize boxplot and histogram of flux rate distribution, identify outliers (Figure S2)
png("GVL_Incubation_Distribution.png",width=6,height=8,units="in",res=300)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8)) # Layout to split the screen
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(sum20$mean_Prr, horizontal=TRUE , ylim=c(-15,25), xaxt="n" , col="#E8D19D", pars=list(boxwex=0.8,cex=2,outcol="#BA8F2D",pch=16), frame=F)
par(mar=c(4, 3.1, 1.1, 2.1))
hist(sum20$mean_Prr , breaks=20, prob=TRUE, col="#7F7F7F", las=1, border=F, main="", ylab="", xlab='P Flux Rate (mg P' *m^-2~day^-1*')', xlim=c(-15,25),cex.lab=1.2,cex.axis=1.2)
lines(density(sum$mean_Prr),lwd=3,col="#3A5750")
dev.off()

### Part 2B - Visualize the P flux distribution by site and sampling event (Figure S3) -----------------------------------------------------------
pttrans=0.5 #transparency value for points
ptpal=brewer.pal(7,"Dark2")
ptrgb<-col2rgb(ptpal)
ptpal<-rgb(red=ptrgb[1,],green=ptrgb[2,],blue=ptrgb[3,],alpha=pttrans*255, maxColorValue = 255) #color palette for transparent points
lnpal=brewer.pal(7,"Dark2") #color palette for solid lines

sum_ordered<-sum %>% 
  mutate(month=factor(month,levels=c("February","April","June","August","October"))) %>% 
  mutate(site.f=factor(site.f,levels=c("4","10","1"))) 

dotplot<-ggplot(sum_ordered, aes(y=mean_Prr,x=site.f,shape=site.f,col=month,alpha=2))+
  geom_jitter(position=position_jitter(0.1),cex=5)+coord_flip()+
  theme_classic(base_size=18)+
  scale_color_manual(values=c(ptpal[4],ptpal[2],ptpal[5],ptpal[1],ptpal[3]))+
  scale_shape_manual(values=c(15,17,16))+theme(axis.text=element_text(color="black",size=12),axis.title=element_text(size=14))+
  labs(x=NULL,y='P Flux Rate (mg P' *m^-2~day^-1*')',cex=3)+theme(legend.position = "none")
ggsave("GVL_FigS2.png",dotplot,width=10,height=8,units="in",dpi=300)

### Part 2C - Quantify the distribution skew, leave-one-out analysis (Table S3) ------------------------------------------------------------------
# Code adapted from: https://helleng.github.io/Data_Mgt_Analysis_and_Graphics_R/Data_Analysis/chap2.html 

# DATA: Work with full incubation dataset, subset for 2020 ("sum" subset for year==2020) - already subset in Part 2a ("sum20")

# visualize distribution of whole incubation dataset (histogram)
right <- rgb(16,78,139, max = 255, alpha = 120, names = "poscol")
hist(sum20$mean_Prr, col = right, pch=19, breaks = 20, lwd=2, lty=1, border = "white", xlim=c(-20,20))
#First statistical moment (i.e., mean) - all data
mean(sum20$mean_Prr) #3.38 mg/m2/d

# SKEW: Calculating skewness of P flux rates (data are moderately positively skewed, m3 = 0.82)
# Functions to compute third statistical moment (m3) and standardized third statistical moment 
m3 <- function(x) {
  n <- length(x)
  dif <- (x - mean(x, na.rm=T))^3
  sum(dif)/n
}

# Standardized third statistical moment
m3_std <- function(x) {
  s <- sd(x, na.rm=T)
  m3(x)/(s^3)
}

# Standardized m3 of distribution of P flux rates
skew_PrrTP <- m3_std(sum20$mean_Prr)
skew_PrrTP #standardized m3=0.08182019

# How to interpret m3 values:
skewness_interpreter <- function(x) {
  if(x == 0) {
    return("symmetric (not skewed)")
  } else if (x > -0.5 & x < 0.5) {
    return("approximately symmetric")
  } else if (x <= -0.5 & x >= -1) {
    return("moderately (negatively) skewed")
  } else if (x >= 0.5 & x <= 1) {
    return("moderately (positively) skewed")
  } else if (x < -1 | x > 1) {
    if (x < -1) {
      return("highly negatively skewed")
    } else {
      return("highly positively skewed")
    }
  } else {
    return("Can't interpret that, I need one numerical value.")
  }
}

skewness_interpreter(skew_PrrTP)

# LEAVE-ONE-OUT ANALYSIS - What HSHM values need to be left out to make the distribution not skewed? (Table S3)
dat <- sum20$mean_Prr
dat<-sort(dat,decreasing=FALSE)

# 1. Skewness - leave highest out
dif <- (dat[1:41] - mean(dat[1:41], na.rm=T))^3
m3 <- sum(dif)/41
# To compute standardized third statistical moment
s <- sd(dat[1:41], na.rm=T)
m3_std <-m3/(s^3)
skewness_interpreter(m3_std) # "moderately (positively) skewed"

# 2. Skewness - leave two highest out
dif <- (dat[1:40] - mean(dat[1:40], na.rm=T))^3
m3 <- sum(dif)/40
# To compute standardized third statistical moment
s <- sd(dat[1:40], na.rm=T)
m3_std <-m3/(s^3)
skewness_interpreter(m3_std) # "moderately (positively) skewed"

# 3. Skewness - leave three highest out
dif <- (dat[1:39] - mean(dat[1:39], na.rm=T))^3
m3 <- sum(dif)/39
# To compute standardized third statistical moment
s <- sd(dat[1:39], na.rm=T)
m3_std <-m3/(s^3)
skewness_interpreter(m3_std) # "moderately (positively) skewed"

# 4. Skewness - leave four highest out
dif <- (dat[1:38] - mean(dat[1:38], na.rm=T))^3
m3 <- sum(dif)/38
# To compute standardized third statistical moment
s <- sd(dat[1:38], na.rm=T)
m3_std <-m3/(s^3)
skewness_interpreter(m3_std) # "moderately (positively) skewed"


#5. Skewness - leave five highest out
dif <- (dat[1:37] - mean(dat[1:37], na.rm=T))^3
m3 <- sum(dif)/37
# To compute standardized third statistical moment
s <- sd(dat[1:37], na.rm=T)
m3_std <-m3/(s^3)
skewness_interpreter(m3_std) #"approximately symmetric"

#First statistical moment if we leave out the 5 high values that cause skewness
mean(dat[1:37]) #1.25 - VERY! different mean value if we hadn't captured HSHM

### PART 3 - SEDIMENT P FLUX RATES AND DISSOLVED OXYGEN STATUS AT THE SEDIMENT-WATER INTERFACE (Figure 3) ----------------------------------------------------------------
# DATA: Subset the summarized "sum2" dataframe generated in Part 1A by sampling year (year==2020)
sum_DO<-dplyr::filter(sum2,year=="2020")

DO_fig<-ggplot(data=sum_DO)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  ylim(-15,20)+xlim(0,13)+
  xlab(bquote('Hypolimnetic DO (mg '*L^-1*')')) + ylab(bquote('P Release Rate (mg P' *m^-2~day^-1*')'))+
  theme(axis.text=element_text(color="black",size=9),axis.title=element_text(size=10))+
  geom_hline(yintercept=0,color="gray40",linetype="solid",size=1.1)+
  geom_errorbar(aes(x=hypo_DO,ymin=average_Prr-sem_Prr,ymax=average_Prr+sem_Prr),color="gray20",width=0,size=1,position=position_dodge(0.05))+
  geom_point(aes(x=hypo_DO,y=average_Prr,shape=month,fill=site.f,color=site.f),size=5)+
  scale_fill_manual(values=c("#1B9E77","#D95F02","#E7298A"),
                    name="Sampling Site",
                    breaks=c(1,10,4),
                    labels=c("Shallow","Intermediate Depth","Deep"))+
  scale_color_manual(values=c("#1B9E77","#D95F02","#E7298A"),
                     name="Sampling Site",
                     breaks=c(1,10,4),
                     labels=c("Shallow","Intermediate Depth","Deep"))+
  scale_shape_manual(values=c(21,24,22,23,25),
                     name="Month",
                     breaks=c("February","April","June","August","October"),
                     labels=c("Winter","Spring","Summer","Late Summer","Autumn"))+
  theme(legend.position=c(0.5,0.08),legend.direction="horizontal",legend.box="vertical",legend.spacing.y=unit(0.12,"cm"),legend.spacing.x=unit(0.05,"cm"),legend.margin=margin(0,0,0,0),
        legend.text=element_text(size=8),legend.title=element_text(size=8,face="bold"))

DO_fig
ggsave("GVL_Figure4.png",DO_fig,width=4.75,height=4.5,units="in", dpi=300)

### PART 4 - SEDIMENT P POOLS (Figure 4, S4) -----------------------------------------------------------------------------------------------------------------------------
# READ IN DATA: Sediment P Speciation, Summarized Data, Green Valley Lake (2020)
sed_sum<-read.csv("GVL_SedimentP_Sum.csv")
# READ IN DATA: Sediment Total P, Summarized Data, Green Valley Lake (2020)
sed_tp<-read.csv("GVL_SedTP_Sum.csv")

### Part 4A -  Mean sediment P composition by sampling site (Figure S4), visualized as a stacked barplot of sediment P speciation --------------

# calculate the average, annual sediment TP concentration by site
sed_tp_sum <- sed_tp %>% 
  mutate(site=factor(site,levels=c("4","10","1"))) %>% 
  group_by(site) %>% 
  mutate(annualmean_sedTP_ug = mean(mean_TP_ug)) %>% 
  select(site,annualmean_sedTP_ug) %>%
  slice(n=1)
# calculate the average, annual concentration of each sediment P species by site
sed_sum2 <- sed_sum %>%
  mutate(site=factor(site,levels=c("4","10","1"))) %>%
  group_by(site) %>%
  mutate(annual_looseP_ug=mean(mean_looseP_ug), annual_redoxP_ug = mean(mean_redoxP_ug), annual_AlP_ug = mean(mean_AlP_ug), annual_orgP_ug = mean(mean_orgP_ug)) %>% 
  select(site, annual_looseP_ug,annual_redoxP_ug,annual_AlP_ug,annual_orgP_ug) %>% 
  slice(n=1)
# join tables and calculate the residual P (total P minus the sum of the measured fractions)
sed_bar<-full_join(sed_sum2,sed_tp_sum,by="site")  
sed_bar$residP_ug <- (sed_bar$annualmean_sedTP_ug-(sed_bar$annual_looseP_ug+sed_bar$annual_redoxP_ug+sed_bar$annual_AlP_ug+sed_bar$annual_orgP_ug))  
# pivot the data table into the long format for making the bar plot
sed_bar2 <- sed_bar %>% 
  pivot_longer(!site, names_to = "P_form", values_to = "concentration_ug") #it worked! I learned pivot!
sed_bar2 <- dplyr::filter(sed_bar2, !P_form=="annualmean_sedTP_ug") #remove total P to make stacked bar plot of other fractions

sed_bar_3 <- sed_bar2 %>% 
  mutate(P_form=factor(P_form,levels=c("annual_looseP_ug","residP_ug","annual_orgP_ug","annual_AlP_ug","annual_redoxP_ug")))

# calculate abundance of each fraction as a percent of the total sediment P
sed_bar_per <- sed_bar %>% 
  mutate(per_loose=(annual_looseP_ug/annualmean_sedTP_ug*100),per_redox=(annual_redoxP_ug/annualmean_sedTP_ug*100),
         per_al=(annual_AlP_ug/annualmean_sedTP_ug*100),per_org=(annual_orgP_ug/annualmean_sedTP_ug*100),per_ref=(residP_ug/annualmean_sedTP_ug*100))

# FIGURE S4
ggplot(data = sed_bar_3, aes(x = site, y = concentration_ug),color="white") +
  geom_col(aes(fill = P_form), color="white",width = 0.7)+ coord_flip()+
  theme_classic(base_size=18)+
  xlab(NULL) + ylab(bquote('P Concentration (µg P ' *g^-1*' Dry Sediment)'))+
  theme(legend.position = c(0.8,0.8),
        axis.text=element_text(color="black",size=14),axis.title=element_text(size=16))+
  scale_fill_manual(values=c("#d53e4f","gray80","#51AAAE","#3288bd","#f46d43"),labels = c("Loosely-Bound", "Refractory", "Labile Organic", "Aluminum-Bound","Redox-Sensitive"))+
  labs(fill="P Species")
ggsave("GVL_barplot.png",width=8,height=8,units="in",dpi=300)

### Part 4B - Timeseries of mobile sediment P fractions by sampling site (Figure 4, top panels) ------------------------------------------------
# work with sed_sum
sed_sum<-sed_sum %>% 
  mutate(site=factor(site,levels=c("1","10","4")))
#Redox P
redox_p<-
  ggplot(data=sed_sum) +
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  ylim(0,1700)+
  geom_line(aes(x=doy,y=mean_redoxP_ug,color=site),size=1) +
  geom_point(aes(x=doy,y=mean_redoxP_ug,size=0.7,color=site,shape=site),stroke=0.7,fill="white") +
  geom_errorbar(aes(x=doy,ymin=mean_redoxP_ug-sem_redoxP_ug,ymax=mean_redoxP_ug+sem_redoxP_ug,color=site),width=0, position=position_dodge(0.05))+
  scale_shape_manual(values=c(19,17,15)) + scale_color_manual(values=c("#1B9E77","#D95F02","#E7298A"))+
  geom_text(x=35,y=1640,label="Redox-Sensitive P",size=3,hjust=0)+
  xlab(" ") + ylab(bquote('P (µg ' *g^-1*' Dry Sediment)'))+
  theme(legend.position = "none",axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,0),"cm"))

#Organic P
org_p<-
  ggplot(data=sed_sum) +
  ylim(0,450)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_line(aes(x=doy,y=mean_orgP_ug,color=site),size=1) +
  geom_point(aes(x=doy,y=mean_orgP_ug,size=0.7,color=site,shape=site),stroke=0.7,fill="white") +
  geom_errorbar(aes(x=doy,ymin=mean_orgP_ug-sem_orgP_ug,ymax=mean_orgP_ug+sem_orgP_ug,color=site),width=0, position=position_dodge(0.05))+
  scale_shape_manual(values=c(19,17,15)) + scale_color_manual(values=c("#1B9E77","#D95F02","#E7298A"))+
  geom_text(x=35,y=430,label="Labile Organic P",size=3,hjust=0)+
  xlab(" ") + ylab(bquote(' '))+
  theme(legend.position = "none",axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,-0.1),"cm"))

#Al-bound P
al_p<-
  ggplot(data=sed_sum) +
  ylim(0,900)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_line(aes(x=doy,y=mean_AlP_ug,color=site),size=1) +
  geom_point(aes(x=doy,y=mean_AlP_ug,size=0.7,color=site,shape=site),stroke=0.7,fill="white") +
  geom_errorbar(aes(x=doy,ymin=mean_AlP_ug-sem_AlP_ug,ymax=mean_AlP_ug+sem_AlP_ug,color=site),width=0, position=position_dodge(0.05))+
  scale_shape_manual(values=c(19,17,15)) + scale_color_manual(values=c("#1B9E77","#D95F02","#E7298A"))+
  geom_text(x=35,y=850,label="Aluminum-Bound P",size=3,hjust=0)+
  xlab("Day of Year") + ylab(bquote('P (µg ' *g^-1*' Dry Sediment)'))+
  theme(legend.position = "none",axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,0),"cm"))

#Loosely-bound P
loose_p<-
  ggplot(data=sed_sum) +
  ylim(0,45)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_line(aes(x=doy,y=mean_looseP_ug,color=site),size=1) +
  geom_point(aes(x=doy,y=mean_looseP_ug,size=0.7,color=site,shape=site),stroke=0.7,fill="white") +
  geom_errorbar(aes(x=doy,ymin=mean_looseP_ug-sem_looseP_ug,ymax=mean_looseP_ug+sem_looseP_ug,color=site),width=0, position=position_dodge(0.05))+
  geom_point(aes(x=doy,y=mean_looseP_ug,size=0.7,color=site,shape=site),stroke=0.7,fill="white") +
  scale_shape_manual(values=c(19,17,15)) + scale_color_manual(values=c("#1B9E77","#D95F02","#E7298A"))+
  geom_text(x=35,y=42,label="Loosely-Bound P",size=3,hjust=0)+
  xlab("Day of Year") + ylab(bquote(''))+
  theme(legend.position = "none",axis.text=element_text(color="black",size=8),axis.title=element_text(size=8))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,-0.1),"cm"))

windows(height=3.5,width=5)  
plotfull_p<-grid.arrange(redox_p, org_p, al_p, loose_p,nrow=2) #need to name the grid.arrange to use in ggsave
ggsave("GVL_2020_SedimentP.png", plotfull_p, width=5, height=3.5, units="in", dpi=300)

### Part 4C - CoDA of sediment P composition to understand variation over space and time (Figure 4, bottom panels) -----------------------------
# for compositional data analysis - will need a dataframe with the mean (over replicates) concentration of P fractions
# loosely-bound, redox-sensitive, al-bound, labile organic P ---- use sed_sum

CoDA_df<-sed_sum %>% 
  mutate(month=factor(month,levels=c("February","April","June","July","August","October"))) %>% 
  mutate(site=factor(site,levels=c("1","10","4"))) %>% 
  select(month,doy,site,mean_looseP_ug,mean_redoxP_ug,mean_AlP_ug,mean_orgP_ug)
# compositional data in columns 4:7

# Take the centered log raio (CLR) of the compositional data
CLR<-cenLR(CoDA_df[, 4:7])

# the cenLR function returns two lists: x.clr gives the CLR of each concentration, gm gives the geometric means for each row (sample)
# for now, we want the CLR values in a data frame so that we can add these values to the main dataset and then perform a PCA
CLRdf<-data.frame(CLR$x.clr)

# Join CLR dataframe to "CoDA_df" dataframe - need to change CLRdf column names first
CLRdf2<-rename(CLRdf,CLR_loose=mean_looseP_ug,CLR_redox=mean_redoxP_ug,CLR_Al=mean_AlP_ug,CLR_Organic=mean_orgP_ug)
core_CLR<-cbind(CoDA_df,CLRdf2)
View(core_CLR)# CLR values in columns 8-11

# PCA on covariance matrix using rda() function, need vegan package
core.rda<-rda(core_CLR[,8:11],scale=FALSE) #selects all rows but only columns 8-11, scale=FALSE uses covariance matrix
summary(core.rda) # PC1 explains 76.77% of variation, PC2 explains 14.35% (sum=91.11%)

# PCA Biplot (Figure 2)
with(CoDA_df,levels(site))
brewer.pal(3,"Dark2")
col_10<-"#D95F02" #orange
col_4<-"#E7298A" #pink
col_1<- "#1B9E77" #teal
colvec<-c(col_1,col_10,col_4)

png("GVL_CoDA.png",width=5,height=4,units="in",res=300)
#windows(height=4,width=5)
par(mfrow=c(1,1),mar=c(3.2,3.3,0.8,0.8),cex.lab=0.7,cex.axis=0.7) #bottom, left, top, right
biplot(core.rda,choices=c(1,2),display="sites",type="points",scaling=3,xlab="",ylab="",ylim=c(-0.6,0.6),xlim=c(-0.8,0.8))#plot PCA, symmetrical scaling, "sites" only
title(xlab="PC1 (76.77)%",ylab="PC2 (14.35%)",line=2) #show % variation explained by each principal component
#ordihull(core.rda,group=CoDA_df$site,scaling=3,label=FALSE,draw="polygon",border=FALSE,alpha=0.2,col=c(col_1,col_10,col_4))
#with(CoDA_df,points(core.rda,display="sites",col="black",pch=c(22,24,21),scaling=3,bg=colvec[site],cex=1.5)) #color points by site
ordihull(core.rda,group=CoDA_df$month,scaling=3,label=FALSE,draw="polygon",lwd=2,col="gray75",border="gray75",alpha=0.2)
with(CoDA_df,points(core.rda,display="sites",col="black",pch=c(22,24,21),scaling=3,bg=colvec[site],cex=1.8)) #color points by site

text(-0.15,0.31,labels="Winter",col="gray30",cex=0.8,font=3)
text(-0.38,-0.15,labels="Spring",col="gray30",cex=0.8,font=3)
text(0.2,-0.5,labels="Late Summer",col="gray30",cex=0.8,font=3)
text(0.47,0.2,labels="Summer",col="gray30",cex=0.8,font=3)
text(0.6,-0.05,labels="Autumn",col="gray30",cex=0.8,font=3)


# envfit() function to add environmental data 
x_clr<-envfit(core.rda,core_CLR[,4:7],choices=c(1,2))
plot(x_clr,choices=c(1,2),at=c(0,0),axis=FALSE,p.max=NULL,add=TRUE,col="black",labels=c(" "," "," "," "))
text(-0.75,0.12,labels="Redox-
     Sensitive",cex=0.8)
text(-0.6,-0.48,labels="Al-Bound",cex=0.8)
text(0.33,-0.13,labels="Labile Organic",cex=0.8)
text(-0.65,-0.05,labels="Loosely-
     Bound",cex=0.8)
dev.off()

### PART 5 - SCALING P FLUXES ACROSS THE LAKEBED TO ESTIMATE LOAD (Figure 5, Table 2) AND SITE MAP (Figure 1) --------------------------------------------------------
### Part 5a. Processing shapefile of bathymetric data ------------------------------------------------------------------------------------------

# Step 1 - Determine the area of various depth contours across the lakebed
# READ IN DATA: shapefile of bathymetric map (make sure all spatial files are in appropriate directory)
bathy<-st_read("GVL_bathymetry.shp") 
st_crs(bathy) # provides metadata on shapefile
# Geometry type = LINESTRING. That means that the various depth contours are represented by line features... so we cannot calculate an area
bathy_mls<-st_cast(bathy, "MULTILINESTRING") #First, cast the linestring shapefile into a multilinestring geometry type
bathy_mp<-st_cast(bathy_mls, "MULTIPOLYGON") #Next, cast that multilinestring into a multipolygon geometry type. This is what we want!

bathy_mp$area<-st_area(bathy_mp) #add a column with the calculated area of each depth contour polygon.
# NOTE. This will be the total, horizontal area at that depth contour

bathy_area<-as.data.frame(bathy_mp) #make into a dataframe
bathy_area<-subset(bathy_area,select = -geometry) #remove the geometry column because it makes saving the csv a pain

write.csv(bathy_area,file="GVL_bathy_area.csv") #write csv to work with

# Step 2 - refine the above csv to focus on main body of reservoir (exclude areas north of sediment retention dikes) and make hyposographic curve
GVL_bath<-read.csv("GVL_bathy_area.csv")
GVL_bath<-GVL_bath %>% 
  select(CONTOUR,area) %>% 
  mutate(depth_m = CONTOUR/3.281) %>% 
  rename(area_m2 = area) %>% 
  mutate(CONTOUR = factor(CONTOUR,levels=c("2","4","6","8","10","12","14","16","18","20","22","24"))) %>% 
  arrange(CONTOUR)

# We want to exclude the areas in the little bays north of the two sediment retention dikes. 
# We will do this by removing the smallest area(s) from the 2, 4, 6, 8, 10, and 12 contours (informed from stuying labeled IDNR bathy map)
GVL_bathy<-GVL_bath[-c(1,2,4,5,7,8,10,11,12,14,15,18,20),]
GVL_bathy2<-GVL_bathy %>% 
  mutate(area_km2 = area_m2/1000000) %>% 
  mutate(CONTOUR = factor(CONTOUR,levels=c("24","22","20","18","16","14","12","10","8","6","4","2"))) %>% 
  arrange(CONTOUR) %>% 
  group_by(CONTOUR) %>% 
  summarise(sum_area_km2=sum(area_km2),sum_area_m2=sum(area_m2))

# Tangent figures - hypsographic curve
ggplot(GVL_bathy2, aes(x=CONTOUR, y=sum_area_km2))+
  geom_bar(stat="identity",fill="cadetblue4",alpha=0.5)+
  geom_point(stat="identity",size=1.5)+
  coord_flip()+theme_classic()+xlab("Depth Contour (feet)")+ylab(bquote('Planar Area ('*km^2*')'))

### Part 5B - Site Map (Figure 1) --------------------------------------------------------------------------------------------------------------
# READ IN DATA: csv with coordinates of sampling sites
sites<-read.csv("GVL_sites.csv")

# Convert sites csv to spatial files
sites_sf<-st_as_sf(sites,coords=c("LONG","LAT"),crs=4893,agr="constant") #convert dataframe with sampling sites to a spatial file (4893 is the CRS code for NAD83 - what shp projected in)

# SELECT lake outline contour and countours for depth intervals associated with 3 sampling sites
outline_shp <- bathy_mp %>% 
  dplyr::filter(CONTOUR =="2")
deep_shp <- bathy_mp %>% 
  dplyr::filter(CONTOUR=="18")
mid_shp <- bathy_mp %>% 
  dplyr::filter(CONTOUR=="12")
shallow_shp <- bathy_mp %>% 
  dplyr::filter(CONTOUR == "4")

# Preparing inset map to help folks find Iowa (tutorial = https://geocompr.github.io/post/2019/ggplot2-inset-maps/)
# Use us_states data from the spData package
data("us_states",package="spData")
us_states_2163 = st_transform(us_states, crs = 2163)
iowa <- us_states_2163 %>% 
  dplyr::filter(NAME=="Iowa")

us <- ggplot() + 
  geom_sf(data = us_states_2163, fill = "white") +
  geom_sf(data = iowa, fill="#e6ab02") +
  theme_void()

# Site Map (Figure 1)
cols<-c("Deep"="#E7298A","Middle"="#D95F02","Shallow"="#1B9E77")
gvl_map <- ggplot() + 
  geom_sf(data = outline_shp, size = 0.7, color = "black", fill = "white")+
  geom_sf(data = shallow_shp, size = 0.7, color = "black", fill = "gray90")+
  geom_sf(data = mid_shp, size=0.7, color = "black", fill = "gray70") +
  geom_sf(data = deep_shp, size=0.7, color="black", fill = "gray50") + 
  geom_sf(data=sites_sf,aes(shape=Site,fill=Site),size=5,stroke=1.1) + 
  scale_shape_manual(values=c(22,24,21),
                     name="Sampling Site",
                     breaks=c("Deep","Middle","Shallow"),
                     labels=c("Deep","Intermediate Depth","Shallow")) + 
  scale_fill_manual(values=c("#E7298A","#D95F02","#1B9E77"),
                    name="Sampling Site",
                    breaks=c("Deep","Middle","Shallow"),
                    labels=c("Deep","Intermediate Depth","Shallow"))+
  annotation_scale(bar_cols=c("gray50","white"),text_cex=0.6) + 
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position=c(0.68,0.66),legend.text=element_text(size=9),legend.title=element_text(size=10,face="bold"),
        legend.key.size=unit(0.6,"cm"),legend.background = element_blank(),
        axis.text=element_text(size=8))

windows(height=5,width=3.5)  
GVL_SiteMap<-ggdraw() +
  draw_plot(gvl_map) +
  draw_plot(us, x=0.5, y = 0.64, width=0.46, height=0.46)

ggsave(filename = "GVL_SiteMap.png", plot = GVL_SiteMap,
       width = 3.5, height = 5, dpi = 300)

### Part 5C - P load estimation (Figure 5, Table 2) -------------------------------------------------------------------------------------------
# READ IN DATA: perimeter and areas of representative depth contours (calculated in Part 5A)
sed_area<-read.csv("GVL_SedimentAreaCALC.csv")
# Calculate the area of sediment for each area of interest. Imagine the area is a ribbon of sediment wrapping around the lakebed. It is a trapezoid.
# ((sed_area$top_perimeter_m+sed_area$bottom_perimeter_m)/2) takes the average of the bases of the trapezoid - the perimeters of the top and bottom contours
# the second term gives the height of the trapezoid via the pythagorean theorem h=sqrt(a^2 + b^2)
#     a=the difference in depth (sed_area$bottom_depth_m-sed_area$top_depth_m)
#     b=the horizontal distance between each depth area (estimate as difference between sqare root of areas) 
sed_area$sediment_area_m2<-((sed_area$top_perimeter_m+sed_area$bottom_perimeter_m)/2)*(sqrt(((sed_area$bottom_depth_m-sed_area$top_depth_m)^2)+(((sqrt(sed_area$top_area_m2))-(sqrt(sed_area$bottom_area_m2)))^2)))
# results: area of representative contours for each sampling site: (Table 2)
#          shallow site (1.2-3.6 m contour) = 3,014,983.6 m2
#          middle site (3.6-5.6 m contour) = 	1,593,676.1 m2
#          deep site (5.6-7.3 m contour) = 	801,158.7 m2 ***Then need to add area of >7.3m contour (836.757 m2). TOTAL = 801,995.5 m2

# Work with incubation summary table ("GVL_Incubation_Sum.csv", df=sum), add areas for each sediment contour, calculate loads
# Read in summarized sediment core incubation data (summarizes mean daily, aerial total P release rate for each sediment core over the course of the incubation experiment)
sum<-read.csv("GVL_Incubation_Sum.csv")
str(sum) # site number is read as integer, need to make a new column with sampling site as a factor
sum$site.f<-as.factor(sum$site) # make a new column with sampling site as a factor, not an integer
sum$n_cores[sum$year=="2020"] <- 3 # add a column that indicates the number of replicate cores used in the experiment (3 for 2020)...
sum$n_cores[sum$year=="2019"] <- 2 #...and 2 for 2019. Will be used to calculate the standard error of the mean
sum$doy[sum$month=="February"] <- 39 # add day of year when samples were collected to make timeseries
sum$doy[sum$month=="April"] <- 117
sum$doy[sum$month=="June"] <- 181
sum$doy[sum$month=="August"] <- 223
sum$doy[sum$month=="October"] <- 298
sum$doy[sum$month=="July"] <- 203
load<-sum %>% 
  dplyr::filter(year=="2020") %>% 
  select(month,doy,site.f,core_unique,n_cores,mean_Prr)
# add the area of each depth contour to the cooresponding sampling site
load$area_m2[load$site.f=="4"] <- 801995.5
load$area_m2[load$site.f=="10"] <- 1593676.1
load$area_m2[load$site.f=="1"] <- 3014983.6
# calculate the daily P load (mgP/day) as the product of the flux rate (mg P/m2/day) and the sediment area (m2)
load$Pload_mgday<-load$mean_Prr*load$area_m2
load$Pload_kgday<-load$Pload_mgday/1000000

# sumamrize the mean and sem P load for each site and event. Calculate the net load each month and propagate the uncertainty (uncertainties add in quadrature)
load_sum<-load %>% 
  mutate(site.f=factor(site.f,levels=c("1","10","4"))) %>% 
  group_by(month, site.f) %>% 
  mutate(mean_load_kgday=mean(Pload_kgday),sem_load_kgday=(sd(Pload_kgday))/sqrt(n_cores)) %>% 
  select(month,doy,site.f,mean_load_kgday,sem_load_kgday) %>% 
  slice(n=1)
load_sum2<-load_sum %>% 
  group_by(month) %>% 
  mutate(net_load_kgday=sum(mean_load_kgday),net_load_error=sqrt((sem_load_kgday)^2)) %>% 
  select(month,doy,net_load_kgday,net_load_error) %>% 
  slice(n=1)

# Daily P Load Estimation (Figure 5)
load_fig<-ggplot()+
  geom_bar(data=load_sum,aes(x=doy, y=mean_load_kgday, fill=site.f), 
           stat="identity", width=30, position=position_dodge2(preserve="single",padding=0.1))+
  scale_fill_manual(values=c("#1B9E77","#D95F02","#E7298A"),
                    name="Mean P Load by Site",
                    breaks=c("1","10","4"),
                    labels=c("Shallow","Intermediate Depth","Deep"))+
  geom_point(data=load_sum2,aes(x=doy,y=net_load_kgday),size=4,col="gray10")+
  geom_errorbar(data=load_sum2,aes(x=doy,ymin=net_load_kgday-net_load_error,ymax=net_load_kgday+net_load_error),col="gray10",size=1,width=0, position=position_dodge(0.05))+
  geom_hline(yintercept=0,col="black")+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlab("Day of Year") + ylab(bquote('Sediment P Load (kg ' * day^-1*')'))+
  theme(axis.text=element_text(color="black",size=9),axis.title=element_text(size=9))+
  theme(legend.position=c(0.3,0.83),legend.text=element_text(size=9),legend.title=element_text(size=9,face="bold"),legend.key.size=unit(0.3,"cm"))

ggsave("GVL_LoadEstimation.png",load_fig,width=3.5,height=3.25,units="in",dpi=300)
