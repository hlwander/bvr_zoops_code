# DVM calcs for 2019 MSNs
# Created 10Oct19

### Description of data --> full water column tows and epi tows from BVR summer 2019 (Jun 10-11)
    #includes samples collected at macrophyes (BVR_l), dam (BVR_d) and pelagic site (BVR_50_p for epi tows collected during MSN ONLY; BVR_50 for full water column tows and tows outside of 24-hour campaigns)
    #samples collected from at noon (x1), midnight (x1), sunset (x4), and sunrise (x4)

#libarries
pacman::p_load(plyr,plotrix,lubridate,dplyr,ggplot2,scales,tidyr,viridis)

#read in zoop summary csv
zoop<- read.csv('output/FCR_ZooplanktonSummary2019.csv',header = TRUE)

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
zoop$rep <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2" | 
                     substrEnd(zoop$sample_ID,4)=="rep3"| substrEnd(zoop$sample_ID,4)=="rep4",
   substrEnd(zoop$sample_ID,1),NA)

#drop rep# from sample ID
zoop$sample_ID <- ifelse(substrEnd(zoop$sample_ID,4)=="rep1" |substrEnd(zoop$sample_ID,4)=="rep2" | 
                           substrEnd(zoop$sample_ID,4)=="rep3" | substrEnd(zoop$sample_ID,4)=="rep4",
                  substr(zoop$sample_ID,1,nchar(zoop$sample_ID)-5),zoop$sample_ID)

#get hour into character format for grouping
zoop$Hour <- format(round(strptime(paste0(zoop$collect_date, zoop$Hour), 
                                   format="%Y-%m-%d %H:%M"),units="hours"),format="%H:%M")

##### Create new df to combine reps over 24 hours
zoop.repmeans <- zoop %>% select(sample_ID,site_no,collect_date,Hour, ZoopDensity_No.pL, BiomassConcentration_ugpL,
                                 TotalBiomass_ug,Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL, Cladocera_totalbiomass_ug, Cladocera_PercentOfTotal,
                                 Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, Cyclopoida_totalbiomass_ug, Cyclopoida_PercentOfTotal,
                                 Rotifera_density_NopL, Rotifera_totalbiomass_ug, Rotifera_BiomassConcentration_ugpL,Rotifera_PercentOfTotal,
                                 Calanoida_PercentOfTotal, Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL, Calanoida_totalbiomass_ug,
                                 Copepoda_PercentOfTotal, Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL, Copepoda_totalbiomass_ug,
                                 nauplius_PercentOfTotal, nauplius_density_NopL, nauplius_BiomassConcentration_ugpL, nauplius_totalbiomass_ug) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(ZoopDensity_No.pL:nauplius_totalbiomass_ug,), list(rep.mean=mean, rep.SE=stderr))

#get hour into posixct for graphing
zoop.repmeans$Hour <- strptime(paste0(as.character(zoop.repmeans$collect_date), zoop.repmeans$Hour),format="%Y-%m-%d %H:%M")
zoop.repmeans$Hour <- as.POSIXct(zoop.repmeans$Hour)

#------------------------------------------------------------------------------#

  #### BVR epi vs hypo figures!! 
  #there should be 6 samples per MSN (epi and full tows from 2 noons and 1 midnight)
  #MSN1 oxycline is ~6m; MSN2 oxycline is 5.2(ish)m

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
                                 substrEnd(zoop.repmeans$sample_ID,3)=="oxy"),] 

#change B_pel_11Jul19_midnight_oxy hour to 00 (not 01)
BVR_pelagic_DVM$Hour[BVR_pelagic_DVM$sample_ID=="B_pel_11Jul19_midnight_oxy"] <-
  "2019-07-11 00:00:00 EDT"

#only select volume, count, and ug cols (plus SE cols)
BVR_pelagic_DVM_raw<- BVR_pelagic_DVM[,c(1:4,73,74,which(substrEnd(colnames(BVR_pelagic_DVM),10)=="n_rep.mean"),
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
BVR_pelagic_DVM_percent<- BVR_pelagic_DVM[,c(1:4,73,74,which(substrEnd(colnames(BVR_pelagic_DVM),14)=="Total_rep.mean"),
                                     which(substrEnd(colnames(BVR_pelagic_DVM),12)=="Total_rep.SE"))]

#drop oxy samples so only have epi or full - NOTE: need new df if looking at meta data
BVR_pelagic_DVM_vol_calculated <- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)!="oxy",]
BVR_pelagic_DVM_raw <- BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="oxy",]
BVR_pelagic_DVM_percent <- BVR_pelagic_DVM_percent[substrEnd(BVR_pelagic_DVM_percent$sample_ID,3)!="oxy",]

#calculate hypo SE for 11Jul19 midnight 
SE.hypo.calcs.raw<- zoop[zoop$sample_ID=="B_pel_11Jul19_midnight" |zoop$sample_ID=="B_pel_11Jul19_midnight_epi",]
matchingcols<- match(substr(colnames(BVR_pelagic_DVM_raw[1:20]),1,14),substr(colnames(SE.hypo.calcs.raw),1,14))
SE.hypo.calcs.raw<- SE.hypo.calcs.raw[,unique(matchingcols)]

SE.hypo.calcs.dens <- zoop[zoop$sample_ID=="B_pel_11Jul19_midnight" |zoop$sample_ID=="B_pel_11Jul19_midnight_epi",c(1:5,11,12,which(substrEnd(colnames(zoop),4)=="NopL"))]
matchingcols.dens<- match(substr(colnames(BVR_pelagic_DVM_vol_calculated[,c(1:4,6:11)]),1,14),substr(colnames(SE.hypo.calcs.dens),1,14))
SE.hypo.calcs.dens<- SE.hypo.calcs.dens[,matchingcols.dens]

#not sure if this is right, but calculating SE of difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}
  
#initialize df
BVR.DVM.calcs<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))
  
#averaged from bvr 2020+2021 data
NetEfficiency2019<- c(0.03664618)

#for loop to fill out epi vs hypo calcs 
#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. 
column.names<- colnames(BVR_pelagic_DVM_vol_calculated[,c(5:18)])
variables<- colnames(BVR_pelagic_DVM_raw[,c(7:20)])
percent<- colnames(BVR_pelagic_DVM_percent[,c(7:12)])
for(i in 1:length(variables)){
  BVR.DVM.calcs[,paste0(column.names,"_epi")[i]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[i]]
  BVR.DVM.calcs[,paste0(column.names,"_hypo")[i]] <- (((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi"  ,"proportional_vol"]) * BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,paste0(variables)[i]] * (1/NetEfficiency2019)) - 
                                                     ((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "proportional_vol"]) *  BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",paste0(variables)[i]]))/ #no neteff needed for epi tows bc pretty much 100% efficient
                                                     (BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,"Volume_unadj"] - BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "Volume_unadj"])  

}
density.percent<- colnames(BVR_pelagic_DVM_vol_calculated[,c(6:11)])
for(i in 1:length(density.percent)){
for(j in 1:length(unique(BVR.DVM.calcs$Hour))){
    BVR.DVM.calcs[j,paste0(density.percent,"_epi_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_epi")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) *100
    BVR.DVM.calcs[j,paste0(density.percent,"_hypo_percent_density")[i]]<- (BVR.DVM.calcs[j,paste0(density.percent,"_hypo")][i]/ sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[i]],BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[i]])) * 100

}       
}


#initialize df
BVR.DVM.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))

#can only calculate hypo SE for 11JUl19 midnight because only time that had reps of epi and full samples
SEonly<- colnames(SE.hypo.calcs.raw)[7:20]
Percentdens <- colnames(SE.hypo.calcs.dens)[5:10]
for(i in 1:length(SEonly)){
  BVR.DVM.calcs.SE[,paste0(column.names,"_epi_SE")[i]] <- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),2)=="SE"][i]
  BVR.DVM.calcs.SE[2,paste0(column.names,"_hypo_SE")[i]]<- SE.diffMean(SE.hypo.calcs.raw[SE.hypo.calcs.raw$sample_ID=="B_pel_11Jul19_midnight",paste0(SEonly)[i]],
                                                           SE.hypo.calcs.raw[SE.hypo.calcs.raw$sample_ID=="B_pel_11Jul19_midnight_epi",paste0(SEonly)[i]])
}
  for (j in 1:length(Percentdens)){
  BVR.DVM.calcs.SE[,paste0(Percentdens,"_epi_percent_density_SE")[j]]<- BVR_pelagic_DVM_vol_calculated[substrEnd(BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),11)=="NopL_rep.SE"][j]
  BVR.DVM.calcs.SE[2,paste0(Percentdens,"_hypo_percent_density_SE")[j]] <- SE.diffMean(SE.hypo.calcs.dens[SE.hypo.calcs.dens$sample_ID=="B_pel_11Jul19_midnight",paste0(Percentdens)[j]],
                                                                                       SE.hypo.calcs.dens[SE.hypo.calcs.dens$sample_ID=="B_pel_11Jul19_midnight_epi",paste0(Percentdens)[j]])
  }
    
#reorder BVR.DVM.calcs so order matches BVR.DVM.calcs.SE
#BVR.DVM.calcs<- BVR.DVM.calcs[,c(1:22,26,23,27,24,28,25,29)]

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

#ridiculous way to work with these stupid facet labels
facet_labeller_bot <- function(variable, value) {
  rep("",8)
}

facet_labeller_top <- function(variable, value) {
  c("","","","Midnight","","","","Noon")
}

#export 2019 dvm stats
write.csv(BVR.DVM.calcs.long,"output/DVM_2019_zoops.csv",row.names = FALSE)

