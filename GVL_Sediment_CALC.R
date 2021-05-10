# Project: Spatiotemporal Variation in Internal P Loading in a Hypereutrophic Reservoir (GVL) - Ecosystems
# Last Modified 19 April 2021
# Contributers: Ellen Albright
# Description: The following script details calculations used to determine the following sediment characteristics
#                 1. Physical: LOI organic matter content, moisture content, bulk density
#                 2. Phosphorus speciation: concentrations of loosely-bound, redox-sensitive, labile organic, and aluminum-bound P
#                 3. Sediment total phosphorus content
#              It also details the data QAQC procedure used to produce the final dataset used for analysis

# The code can be run with the following datasets:
#     "GVL_SedimentPhys_RAW.csv" - Sediment Physical Data - RAW, hand-entered data from 5 sediment sampling events from Feb-Oct, 2020 as well as 1 from July 2019
#     "GVL_SedimentP_RAW.csv" - Sediment P Speciation Data - RAW, hand-entered data from 5 sediment sampling events from Feb-Oct, 2020
# The code produces the following datasets:
#     "GVL_SedimentPhys_Tidy.csv" - calculated LOI OM, moisture content, and sediment bulk density
#     "GVL_SedimentP_Tidy.csv" - calculated concentrations of loosely-bound, redox sensitive, labile organic, and aluminum-bound P
#     "GVL_SedimentP_Sum.csv" - summary of the mean and standard error of three laboratory replicates for each sediment P fraction at each site and sampling event

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
### 1a. SEDIMENT PHYSICAL DATA --------------------------------------------------------------------------------------------------------------------------------------------
# Read in raw sediment physical data
phys.raw<-read.csv("GVL_SedimentPhys_RAW.csv")
names(phys.raw)

# Calculations of sediment physical characteristics
#   1. Dry Mass (g) = the mass of the tin and the dry sediment sample - the tin mass
dry_g<-phys.raw$dry_tin_g-phys.raw$tin_g
#   2. Moisture Content (fraction) = (wet mass - dry mass)/wet mass 
MCfrac<-(phys.raw$wet_g-dry_g)/phys.raw$wet_g
#   3. Combusted Mass (g) = mass of tin and ashed sediment - the tin mass
ash_g<-(phys.raw$ash_tin_g-phys.raw$tin_g)
#   4. LOI Organic Matter Content (%) = (dry-ashed mass)/dry mass *100
OM_per<-((dry_g-ash_g)/dry_g)*100
#   5. Bulk Density (g/cm3) - estimated using MC and LOI OM, Principles of Lake Sedimentology (Håkanson and Jansson 2002)
bd<-260/(100+(1.6*((MCfrac*100)+(OM_per/100)*(100-(MCfrac*100)))))

# Bind all columns to write complete physical dataset ('phys')
phys<-cbind(phys.raw,dry_g,MCfrac,ash_g,OM_per,bd)

# Select only the rows and columns that you want for the final physical data table
phys<-phys %>% 
  filter(!is.na(year)) %>% 
  select(sample_id,year,lake,month,doy,site,lab_rep,MCfrac,OM_per,bd)
  
# Write a lil csv with the final (ie everything calcualted) dataset
write.csv(phys,file="GVL_SedimentPhys_Tidy.csv")

# Determine the average sediment moisture content for each sampling site, per season (for use in sediment P extraction calculations)
tapply(phys$MCfrac[phys$month=="February"],phys$site[phys$month=="February"],mean)   #4 = 0.85 #10 = 0.69 
tapply(phys$MCfrac[phys$month=="April"],phys$site[phys$month=="April"],mean)   #4 = 0.86 #10 = 0.78 #1 = 0.72
tapply(phys$MCfrac[phys$month=="June"],phys$site[phys$month=="June"],mean)   #4 = 0.87 #10 = 0.77 #1 = 0.71
tapply(phys$MCfrac[phys$month=="August"],phys$site[phys$month=="August"],mean)   #4 = 0.87 #10 = 0.79 #1 = 0.73
tapply(phys$MCfrac[phys$month=="October"],phys$site[phys$month=="October"],mean) #4 = 0.87 #10 = 0.83 #1 = 0.71
tapply(phys$MCfrac[phys$month=="July"],phys$site[phys$month=="July"],mean) #1=0.76 #2=0.73 #3=0.85 #4=0.87 #5=0.84 #6=0.70

# Quick visualization
# Sediment organic matter content (OM) and dry bulk density (bd)
par(mfrow=c(3,2),cex=1.2,omi=c(0.3,0.2,0.2,0.2), mai = c(0.2,1,0.1,0.1))
boxplot((phys[phys$Site.ID=="1","OM"])~(phys[phys$Site.ID=="1","DOY"]),boxwex=0.5,las=1,col=site1,ylim=c(6,18),ylab=NULL,xaxt='n')
text(3.5,17,labels="West Inlet",cex=1)
boxplot((phys[phys$Site.ID=="1","bd"])~(phys[phys$Site.ID=="1","DOY"]),boxwex=0.5,las=1,col=site1,ylim=c(1,1.25),ylab=NULL,xaxt='n')
text(3.5,1.23,labels="West Inlet",cex=1)
boxplot((phys[phys$Site.ID=="3","OM"])~(phys[phys$Site.ID=="3","DOY"]),boxwex=0.5,las=1,col=site3,ylim=c(6,18),ylab="Organic Matter (%)",xaxt='n')
text(3.5,17,labels="West Arm",cex=1)
boxplot((phys[phys$Site.ID=="3","bd"])~(phys[phys$Site.ID=="3","DOY"]),boxwex=0.5,las=1,col=site3,ylim=c(1,1.25),ylab="Bulk Density (g/cm3)",xaxt='n')
text(3.5,1.23,labels="West Arm",cex=1)
boxplot((phys[phys$Site.ID=="4","OM"])~(phys[phys$Site.ID=="4","DOY"]),boxwex=0.5,las=1,col=site4,ylim=c(6,18),ylab=NULL,xlab="DOY")
text(3.5,17,labels="Deep Site",cex=1)
boxplot((phys[phys$Site.ID=="4","bd"])~(phys[phys$Site.ID=="4","DOY"]),boxwex=0.5,las=1,col=site4,ylim=c(1,1.25),ylab=NULL,xlab="DOY")
text(3.5,1.23,labels="Deep Site",cex=1)

### 1b. SEDIMENT PHOSPHORUS SPECIATION DATA --------------------------------------------------------------------------------------------------------------------------------------------
# Read in raw sediment P extraction data
raw.sed<-read.csv("GVL_SedimentP_RAW.csv")
names(raw.sed)

# Calcuations of sediment P extraction results
### NOTES: Units will be in mgP/g dry sediment for first calculation. Can multiply by 1000 if you want to work in ug/L
#   1. Determine the dry mass equivalent of the fresh (wet) sediment subsample used based on average moisture content for that site (see physical data calculations above)
DryMassEq<-raw.sed$wet_g*(1-raw.sed$MCfrac)

#   2. NaCl Extraction for pore water and surface-sorbed P
looseP_mg<-((((raw.sed$NaCl_TP1_mgL+raw.sed$NaCl_TP2_mgL)/2)*raw.sed$NaCl_dil)*raw.sed$NaCl_L)/DryMassEq
looseP_ug<-looseP_mg*1000
### Calculation steps: average TP concentrations (lab dups on AQ2), multiply by dilution factor, multiply this corrected TP concentration (mg/L) by volume NaCl used (0.1 L), divide by dry mass sediment used (g)

#   3. BD Extraction for redox-sensitive P
redoxP_mg<-((((raw.sed$BD_TP1_mgL+raw.sed$BD_TP2_mgL)/2)*raw.sed$BD_dil)*raw.sed$BD_L)/DryMassEq
redoxP_ug<-redoxP_mg*1000
### Calculation steps: average TP concentrations (lab dups on AQ2), multiply by dilution factor, multiply this corrected TP concentration (mg/L) by volume BD+NaCl used (0.2 L), divide by dry mass sediment used (g)

#4. NaOH Extraction for Al-bound P and labile organic P
#4a. Al-bound P (also includes P associated with non-reducible Fe and Mn-oxides). Use SRP results
AlP_mg<-((((raw.sed$NaOH_SRP1_mgL+raw.sed$NaOH_SRP2_mgL)/2)*raw.sed$NaOH_SRP_dil)*raw.sed$NaOH_L)/DryMassEq
AlP_ug<-AlP_mg*1000
#4b. Labile organic P. Use TP results (will give concentration of Al-P AND labile oragnic P), subtract Al-P concentration from 4a.
orgP_mg<-(((((raw.sed$NaOH_TP1_mgL+raw.sed$NaOH_TP2_mgL)/2)*raw.sed$NaOH_TP_dil)*raw.sed$NaOH_L)/DryMassEq)-AlP_mg
orgP_ug<-orgP_mg*1000

# Bind all columns to write complete sediment P dataset ('sed'). Select relevant columns for tidy data table. Write csv.
sed<-cbind(raw.sed,DryMassEq,looseP_mg,looseP_ug,redoxP_mg,redoxP_ug,AlP_mg,AlP_ug,orgP_mg,orgP_ug)
sed<-sed %>% 
  select(sample_id,year,lake,month,doy,site,lab_rep,looseP_mg,looseP_ug,redoxP_mg,redoxP_ug,AlP_mg,AlP_ug,orgP_mg,orgP_ug)
write.csv(sed,file="GVL_SedimentP_Tidy.csv")

# Write a summary data table that includes the mean and sem for each fraction, site, event. Summarizing the three laboratory replicates makes later visualizations easier!
tidyP<-read.csv("GVL_SedimentP_Tidy.csv")
names(tidyP)

sumP<-tidyP %>%
  group_by(doy,site) %>% 
  mutate(mean_looseP_ug = mean(looseP_ug), sem_looseP_ug = (sd(looseP_ug))/sqrt(3),
            mean_redoxP_ug = mean(redoxP_ug), sem_redoxP_ug = (sd(redoxP_ug))/sqrt(3),
            mean_AlP_ug = mean(AlP_ug), sem_AlP_ug = (sd(AlP_ug))/sqrt(3),
            mean_orgP_ug = mean(orgP_ug), sem_orgP_ug = (sd(orgP_ug))/sqrt(3)) %>% 
  select(year,lake,month,doy,site,mean_looseP_ug,sem_looseP_ug,mean_redoxP_ug,sem_redoxP_ug,mean_AlP_ug,sem_AlP_ug,mean_orgP_ug,sem_orgP_ug) %>% 
  slice(n=1)
write.csv(sumP,file="GVL_SedimentP_Sum.csv")

### 1c. SEDIMENT TOTAL PHOSPHORUS DATA --------------------------------------------------------------------------------------------------------------------------------------------
rawTP<-read.csv("GVL_SedTP_RAW.csv",as.is=T)
View(rawTP)

#Calculate average TP concentration (average dups), corrected for dilution factor (multiply by) - mg/L
meanTPdil<-((rawTP$TP_mgL+rawTP$TP_dup_mgL)/2)*rawTP$dil_fac

#Calculate pH-corrected TP concentration (mg/L). Multiply the corrected TP concentration from above by the ratio of the change in mass due to the pH-adjustment
corrected_TP<-(meanTPdil*((rawTP$post_g-rawTP$jar_g)/(rawTP$pre_g-rawTP$jar_g)))

#Calculate sediment TP concentration (mg P/g dry sediment). Multiply pH-corrected TP concentration by dilution/extraction volume (L) and divide by dry mass sediment used (g)
TP<-(corrected_TP*rawTP$dilution_L)/rawTP$sed_g

#convert to ug P/g dry sediment
TP_ug<-TP*1000

#Bind it all together
tp<-cbind(rawTP,meanTPdil,corrected_TP,TP,TP_ug)
tp<-tp %>% 
  select(sample_id,year,lake,month,doy,site,lab_rep,TP_ug)

#Write CSV
write.csv(tp,file="GVL_SedTP_Tidy.csv")

# Write a summary data table that includes the mean sediment TP for each site, event. Summarizing the three laboratory replicates makes later visualizations easier!

sumTP<-tp %>%
  group_by(doy,site) %>% 
  mutate(mean_TP_ug = mean(TP_ug)) %>% 
  select(year,lake,month,doy,site,mean_TP_ug) %>% 
  slice(n=1)
write.csv(sumTP,file="GVL_SedTP_Sum.csv")

