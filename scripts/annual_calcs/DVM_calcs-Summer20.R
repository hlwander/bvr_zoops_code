# DVM calcs for 2020 MSN
# Created 17Nov2020

### Description of data --> full water column tows and epi tows from FCR + BVR summer 2020 (outside of 12-13Aug 20 MSN)
    #includes samples collected at macrophyes (BVR_l) and pelagic site (BVR_50_p for epi tows collected during MSN ONLY; BVR_50 for full water column tows and tows outside of 24-hour campaigns)
    #samples collected from noon (x1), midnight (x1), sunset (x4), and sunrise (x4)

#read in libraries
pacman::p_load(plyr,plotrix,lubridate,dplyr,ggplot2,scales,tidyr,viridis)

#read in zoop summary csv
zoop<- read.csv('output/FCR_ZooplanktonSummary2020.csv',header = TRUE)

#create function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#make sure sample_ID is class character
zoop$sample_ID<- as.character(zoop$sample_ID)

#merge collect_date and hour in a new column
zoop$date<- paste(zoop$collect_date,zoop$Hour,sep=" ")
#get times into date format (character here)
zoop$date<- format(as.POSIXct(zoop$date,format="%Y-%m-%d %H:%M"), format="%Y-%m-%d %H:%M:%S")
#convert to posixct date format
zoop$date<- as.POSIXct(zoop$date, format="%Y-%m-%d %H:%M")

#order by site, then hour
zoop<- arrange(zoop,zoop$site_no,zoop$date)

#pull rep # off as new column
zoop$rep <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2" | substrEnd(zoop$sample_ID,4)=="rep3"| substrEnd(zoop$sample_ID,4)=="rep4",
   substrEnd(zoop$sample_ID,1),NA)

#drop rep# from sample ID
zoop$sample_ID <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2" | substrEnd(zoop$sample_ID,4)=="rep3" | substrEnd(zoop$sample_ID,4)=="rep4",
                  substr(zoop$sample_ID,1,nchar(zoop$sample_ID)-5),zoop$sample_ID)

#get hour into character format for grouping
zoop$Hour <- format(round(strptime(paste0(zoop$collect_date, zoop$Hour), format="%Y-%m-%d %H:%M"),units="hours"),format="%H:%M")
#manually change hour of some samples (rounding problems)
zoop$Hour[zoop$sample_ID=="B_pel_12Aug20_sunset_epi_h1"] <- "18:00"
zoop$Hour[zoop$sample_ID=="B_pel_13Aug20_midnight"] <- "00:00"
zoop$Hour[zoop$sample_ID=="B_pel_13Aug20_noon"] <- "12:00"
zoop$Hour[zoop$sample_ID=="B_pel_12Aug20_noon"] <- "12:00"
zoop$Hour[zoop$sample_ID=="B_pel_12Aug20_sunset_epi_h3"] <- "20:00"
zoop$Hour[zoop$sample_ID=="B_pel_12Aug20_sunset_epi_h4"] <- "21:00"
            
#drop 20 um samples
zoop <- zoop[substrEnd(zoop$sample_ID,2)!="20",]

#drop schindler samples
zoop <- zoop[substr(zoop$sample_ID,1,20)!="B_pel_12Aug20_schind",]

#drop the 06-29 samples (these were test samples)
zoop <- zoop[substrEnd(zoop$sample_ID,4)!="filt",]

##### Create new df to combine reps over 24 hours
zoop.repmeans <- zoop %>% select(sample_ID,site_no,collect_date,Hour, Volume_L, Volume_unadj, proportional_vol, ZoopDensity_No.pL, OverallCount_n, TotalBiomass_ug,
                                 BiomassConcentration_ugpL,Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL, CladoceraCount_n, Cladocera_totalbiomass_ug, Cladocera_PercentOfTotal,
                                 Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, CyclopoidaCount_n, Cyclopoida_totalbiomass_ug, Cyclopoida_PercentOfTotal,
                                 Rotifera_density_NopL,Rotifera_BiomassConcentration_ugpL, RotiferaCount_n, Rotifera_totalbiomass_ug, Rotifera_PercentOfTotal, 
                                 Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL, CalanoidaCount_n, Calanoida_PercentOfTotal, Calanoida_totalbiomass_ug,
                                 Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL, CopepodaCount_n, Copepoda_PercentOfTotal, Copepoda_totalbiomass_ug,
                                 nauplius_density_NopL, nauplius_BiomassConcentration_ugpL, naupliusCount_n, nauplius_PercentOfTotal, nauplius_totalbiomass_ug) %>%
  group_by(sample_ID, site_no, collect_date, Hour) %>%
  summarise_at(vars(Volume_L:nauplius_totalbiomass_ug,), funs(rep.mean=mean, rep.SE=stderr)) #CHECK WARNING HERE!!
#'`funs()` was deprecated in dplyr 0.8.0.
#'Please use a list of either functions or lambdas:

# Simple named list: list(mean = mean, median = median)

# Auto named with `tibble::lst()`: tibble::lst(mean, median)

# Using lambdas list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))

#get hour into posixct for graphing
zoop.repmeans$Hour <- strptime(paste0(as.character(zoop.repmeans$collect_date), zoop.repmeans$Hour),format="%Y-%m-%d %H:%M")
zoop.repmeans$Hour <- as.POSIXct(zoop.repmeans$Hour)

#make sure zoop.repmeans is a dataframe
zoop.repmeans <- data.frame(zoop.repmeans)

#--------------------------------------#
#            BVR hypo calcs            #
#Note: oxycline = epi sample depth (4m)#
#--------------------------------------#

#sum up counts by sample/site/day for DVM analyses + figs
BVR_counts <- zoop %>% select(sample_ID,site_no,collect_date,Hour, OverallCount_n,
  CladoceraCount_n, CyclopoidaCount_n, RotiferaCount_n, CalanoidaCount_n, CopepodaCount_n, naupliusCount_n) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(OverallCount_n:naupliusCount_n), list(rep.mean=mean, rep.SE=stderr))

#add unadjusted volume
BVR_counts$Volume_unadj<- (zoop %>% select(sample_ID,site_no,collect_date,Hour, Volume_unadj) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(Volume_unadj),list(rep.mean=mean)))$rep.mean

#add proportional volume for numerator hypo calc
BVR_counts$proportional_vol<- (zoop %>% select(sample_ID,site_no,collect_date,Hour, proportional_vol) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(proportional_vol),list(rep.mean=mean)))$rep.mean

#get BVR_counts df in same order as zoop.repmeans df
BVR_counts<- BVR_counts[order(match(paste0(BVR_counts$sample_ID,BVR_counts$site_no,BVR_counts$collect_date), 
                                    paste0(zoop.repmeans$sample_ID,zoop.repmeans$site_no,zoop.repmeans$collect_date))),]

#add counts and vol to zoop.repmeans
zoop.repmeans[,paste0(colnames(BVR_counts[5:20]))]<- BVR_counts[5:20]

#new dfs for DVM data (BVR_pelagic_DVM_raw is just raw #/ug; BVR_pelagic_DVM_vol_calculated is #/L and ug/L)
BVR_pelagic_DVM<- zoop.repmeans[(zoop.repmeans$site_no=="BVR_50" |zoop.repmeans$site_no=="BVR_50_p") &
                                 (substrEnd(zoop.repmeans$sample_ID,5)=="night" | substrEnd(zoop.repmeans$sample_ID,4)=="noon" |
                                 substrEnd(zoop.repmeans$sample_ID,9)=="night_epi" | substrEnd(zoop.repmeans$sample_ID,8)=="noon_epi" |
                                 substrEnd(zoop.repmeans$sample_ID,5)=="oxy"),] 

#only select volume, count, and ug cols (plus SE cols)
BVR_pelagic_DVM_raw<- BVR_pelagic_DVM[,c(1:4,79,80,which(substrEnd(colnames(BVR_pelagic_DVM),10)=="n_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),11)=="ug_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),8)=="n_rep.SE"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),9)=="ug_rep.SE"))]

BVR_pelagic_DVM_vol_calculated <- BVR_pelagic_DVM[,c(1:4,which(substrEnd(colnames(BVR_pelagic_DVM),16)=="y_No.pL_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),15)=="y_NopL_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),15)=="n_ugpL_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),13)=="n_ugpL_rep.SE"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),14)=="y_No.pL_rep.SE"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),13)=="y_NopL_rep.SE"))]

#another one for percent calcs
BVR_pelagic_DVM_percent<- BVR_pelagic_DVM[,c(1:4,79,80,which(substrEnd(colnames(BVR_pelagic_DVM),14)=="Total_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),12)=="Total_rep.SE"))]
  
#initialize df
BVR.DVM.calcs<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))
  
#for loop to fill out epi vs hypo calcs 
#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. 
column.names<- colnames(BVR_pelagic_DVM_vol_calculated[,c(5:18)])
variables<- colnames(BVR_pelagic_DVM_raw[,c(7:20)])
percent<- colnames(BVR_pelagic_DVM_percent[,c(7:12)])
for(i in 1:length(variables)){
  BVR.DVM.calcs[,paste0(column.names,"_epi")[i]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  BVR.DVM.calcs[,paste0(column.names,"_hypo")[i]] <- (((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,"proportional_vol"]) * BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]] * (1/0.035)) - #combined 2020 and 2021 neteff
                                                     ((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "proportional_vol"]) *  BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]]))/ #no neteff needed for epi tows bc pretty much 100% efficient
                                                     (BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj"] - BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "Volume_unadj"])  

}

density.percent<- colnames(BVR_pelagic_DVM_vol_calculated[,c(6:11)])
for(i in 1:length(density.percent)){
for(j in 1:length(unique(BVR.DVM.calcs$Hour))){
    BVR.DVM.calcs[j,paste0(density.percent,"_epi_percent_density")[i]]<- (BVR.DVM.calcs[3,paste0(density.percent,"_epi")][4]/ sum(BVR.DVM.calcs[3,paste0(density.percent,"_epi")[4]],BVR.DVM.calcs[3,paste0(density.percent,"_hypo")[4]])) *100
    BVR.DVM.calcs[j,paste0(density.percent,"_hypo_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_hypo")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) * 100

}       
}

#initialize df
BVR.DVM.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))

#not sure if this is right, but calculating SE of difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}

#pull only noon/midnight samples
DVM_samples_raw <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | substrEnd(zoop$sample_ID,5)=="night" | substrEnd(zoop$sample_ID,8)=="noon_epi" | substrEnd(zoop$sample_ID,9)=="night_epi") & zoop$site_no!="BVR_l",]
matchingcols <- match(substr(colnames(BVR_pelagic_DVM_raw[1:20]),1,14),substr(colnames(DVM_samples_raw),1,14))
DVM_samples_raw<- DVM_samples_raw[,unique(matchingcols)]
                        
DVM_samples_dens <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | substrEnd(zoop$sample_ID,5)=="night" | substrEnd(zoop$sample_ID,8)=="noon_epi" | substrEnd(zoop$sample_ID,9)=="night_epi") & zoop$site_no!="BVR_l",]
matchingcols <- match(substr(colnames(BVR_pelagic_DVM_vol_calculated [,c(1:4,6:11)]),1,14),substr(colnames(DVM_samples_dens),1,14))
DVM_samples_dens<- DVM_samples_dens[,unique(matchingcols)] 

#separate full vs epi samples
FullSamples <- DVM_samples_raw$sample_ID[c(3,7,11)]
EpiSamples<- DVM_samples_raw$sample_ID[c(1,5,9)]

#calculate hypo SE 
SEonly<- colnames(DVM_samples_raw)[7:20]
Percentdens <- colnames(DVM_samples_dens)[5:10]

for(i in 1:length(SEonly)){
  BVR.DVM.calcs.SE[,paste0(column.names,"_epi_SE")[i]] <- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),2)=="SE"][i]
  BVR.DVM.calcs.SE[1,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[1],paste0(SEonly)[i]],
                                                                      DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[1],paste0(SEonly)[i]])
  BVR.DVM.calcs.SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[2],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[2],paste0(SEonly)[i]])
  BVR.DVM.calcs.SE[3,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[3],paste0(SEonly)[i]],
                                                                       DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[3],paste0(SEonly)[i]])
  } #Note: just repeating hypo SE calc for each date because I can't figure out the for loop across rows and columns...


  for (j in 1:length(Percentdens)){
  BVR.DVM.calcs.SE[,paste0(Percentdens,"_epi_percent_density_SE")[j]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),11)=="NopL_rep.SE"][j]
  BVR.DVM.calcs.SE[1,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[1],paste0(Percentdens)[j]],
                                                                                      DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[1],paste0(Percentdens)[j]])
  BVR.DVM.calcs.SE[2,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[2],paste0(Percentdens)[j]],
                                                                                       DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[2],paste0(Percentdens)[j]])
  BVR.DVM.calcs.SE[3,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[3],paste0(Percentdens)[j]],
                                                                                          DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[3],paste0(Percentdens)[j]])
  }

#wide to long for both dfs separately
BVR.DVM.calcs.long <-  BVR.DVM.calcs %>%
    gather(metric,value, ZoopDensity_No.pL_rep.mean_epi:nauplius_density_NopL_rep.mean_hypo_percent_density) %>%
    mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))
BVR.DVM.calcs.SE.long <-  BVR.DVM.calcs.SE %>%
gather(taxa.metric,SE,ZoopDensity_No.pL_rep.mean_epi_SE:nauplius_density_NopL_hypo_percent_density_SE)

#add the SE column from df2 to df1 for combined df
BVR.DVM.calcs.long$SE <- BVR.DVM.calcs.SE.long$SE

#add watercolumn and hour columns
BVR.DVM.calcs.long$WaterColumn <- ifelse((substrEnd(BVR.DVM.calcs.long$metric,3)=="epi" |substrEnd(BVR.DVM.calcs.long$metric,19)=="epi_percent_density") ,"epilimnion","hypolimnion")
BVR.DVM.calcs.long$Hour <- ifelse(substrEnd(BVR.DVM.calcs.long$DateTime,5)=="12:00","noon","midnight")
BVR.DVM.calcs.long$Taxa <- substr(BVR.DVM.calcs.long$metric,1,9)

#shorten date time to just date
BVR.DVM.calcs.long$DateTime <- substr(BVR.DVM.calcs.long$DateTime,1,nchar(BVR.DVM.calcs.long$DateTime)-6)

#replace NAN with 0
BVR.DVM.calcs.long$value[is.nan(BVR.DVM.calcs.long$value)] <- 0

#export 2020 dvm stats
write.csv(BVR.DVM.calcs.long, "output/DVM_2020_zoops.csv")

