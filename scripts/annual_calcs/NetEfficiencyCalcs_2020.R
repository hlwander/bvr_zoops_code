#script for calculating zooplankton net efficiency in Beaverdam Reservoir
#Created: 20 November 2020
#for some reason, did not collect full water column tow at the same time as schindlers (??!!)
#options are to either average noon and midnight full densities or just use one or the other

#load in packages
pacman::p_load(dplyr,ggplot2)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data 
zoop<- read.csv('output/FCR_ZooplanktonSummary2020.csv',header = TRUE)

#select cols to keep
keep <- c("sample_ID","site_no","collect_date","Hour","DepthOfTow_m","Zooplankton_No.","Volume_L","Volume_unadj","mesh_size_μm")

zoop_totalcount <- zoop[!is.na(zoop$mesh_size_μm) & zoop$site_no!="BVR_l" & zoop$site_no!="BVR_trap",keep]
                   
#convert hour to character
zoop_totalcount$Hour <- as.character(zoop_totalcount$Hour)

#round to nearest hour
zoop_totalcount$Hour <- format(strptime(zoop_totalcount$Hour,format="%H"),"%H")
zoop_totalcount$Hour <- ifelse(zoop_totalcount$Hour=="23"| zoop_totalcount$Hour=="01", "00", 
                               ifelse(zoop_totalcount$Hour=="11" |zoop_totalcount$Hour=="13"| zoop_totalcount$Hour=="14", "12", 
                                      ifelse(zoop_totalcount$Hour=="03","04",zoop_totalcount$Hour)))

#manually change hour for 3 schindlers
zoop_totalcount$Hour[zoop_totalcount$sample_ID=="B_pel_12Aug20_schind_10.0_rep2"]<- 19
zoop_totalcount$Hour[zoop_totalcount$sample_ID=="B_pel_12Aug20_schind_8.0_rep2"]<- 19
zoop_totalcount$Hour[zoop_totalcount$sample_ID=="B_pel_12Aug20_schind_9.0_rep2"]<- 19


#replace vol_unadj column with NA for schindler traps
zoop_totalcount$Volume_unadj[zoop_totalcount$mesh_size_μm ==61] <- NA

#------------------------------------------------------------------------------#
#### BVR and FCR 2020 net efficiency ####

#separate tow data from schindler data
Schindler_totalCount <- zoop_totalcount[zoop_totalcount$mesh_size_μm==61,]
Tow_totalCount_final <- zoop_totalcount[zoop_totalcount$mesh_size_μm==80 & zoop_totalcount$DepthOfTow_m>=9,]
Epi_tow_totalCount_final <- zoop_totalcount[zoop_totalcount$mesh_size_μm==80 & zoop_totalcount$DepthOfTow_m<9,]

#add a column for rep
Schindler_totalCount$Rep <- substrEnd(Schindler_totalCount$sample_ID,1)

#sum total count for rep1 vs rep2 for each reservoir
Schindler_totalCount_final <- Schindler_totalCount %>% select(collect_date,site_no,Hour,DepthOfTow_m,Zooplankton_No.,Rep) %>%
  group_by(site_no,collect_date,Rep)  %>% summarise(DepthOfTow_m=max(DepthOfTow_m),TotalCount_n=sum(Zooplankton_No.))
  
#new df based on epi tows
Schindler_totalCount_epi_final <- Epi_tow_totalCount_final[,c(1:5,7,8)]

#add rep column
Schindler_totalCount_epi_final$Rep <- substrEnd(Schindler_totalCount_epi_final$sample_ID,1)

#if no rep, then change to 1
Schindler_totalCount_epi_final$Rep <- ifelse(Schindler_totalCount_epi_final$Rep=="i" | Schindler_totalCount_epi_final$Rep=="y" , 1,Schindler_totalCount_epi_final$Rep)

#convert date column to date format
Schindler_totalCount_epi_final$collect_date <- as.Date(Schindler_totalCount_epi_final$collect_date)

#loop to sum all depths less than the rounded depth
for(i in 1:length(Schindler_totalCount_epi_final$DepthOfTow_m)){
  Schindler_totalCount_epi_final$Zooplankton_No.[i] <- sum(Schindler_totalCount$Zooplankton_No.[which(Schindler_totalCount$DepthOfTow_m <= Schindler_totalCount_epi_final$DepthOfTow_m[i] & 
                                                                                                        as.Date(Schindler_totalCount$collect_date) %in% Schindler_totalCount_epi_final$collect_date[i] &
                                                                                                        Schindler_totalCount$Hour %in% Schindler_totalCount_epi_final$Hour[i] &
                                                                                                        Schindler_totalCount$Rep %in% Schindler_totalCount_epi_final$Rep[i] &
                                                                                                        substr(Schindler_totalCount$site_no,1,1) %in% substr(Schindler_totalCount_epi_final$site_no,1,1)[i])])
}

#remove rows with NA so we only include tows that match up with schindlers
Schindler_totalCount_epi_final <- Schindler_totalCount_epi_final[Schindler_totalCount_epi_final$Zooplankton_No.>0,]

#initialize df
Density.neteff <- data.frame("sample_ID"=unique(Tow_totalCount_final$sample_ID))
Density.neteff_epi <- data.frame("sample_ID"=unique(Schindler_totalCount_epi_final$sample_ID))

#### Calculate APPARENT density from vertical tows ####
#count / volume_unadj (because we want the whole total volume not volume counted)
for(i in 1:length(Density.neteff$sample_ID)){
  Density.neteff$Apparent_dens[i] <-  Tow_totalCount_final$Zooplankton_No.[i] / Tow_totalCount_final$Volume_unadj[i]
}

for(i in 1:length(Density.neteff_epi$sample_ID)){
  Density.neteff_epi$Apparent_dens[i] <-  Schindler_totalCount_epi_final$Zooplankton_No.[i] / Schindler_totalCount_epi_final$Volume_unadj[i]
}

#initialize actual dens column
Density.neteff$Actual_dens <- NA
Density.neteff_epi$Actual_dens <- NA

#### Calculate ACTUAL density from schindler traps ####
#count / 30L (vol of schindler trap) * # of samples (pooled from all depths)
#repeating these for noon and midnight samples to see if there is a difference

  Density.neteff$Actual_dens[c(1,3,5)] <-  Schindler_totalCount_final$TotalCount_n[1] / (30*11)
  Density.neteff$Actual_dens[c(2,4,6)] <-  Schindler_totalCount_final$TotalCount_n[2] / (30*11)
  Density.neteff$Actual_dens[c(8,9)] <-  Schindler_totalCount_final$TotalCount_n[3] / (30*10)
  Density.neteff$Actual_dens[c(10,11)] <-  Schindler_totalCount_final$TotalCount_n[4] / (30*10)
  #NOTE - check that these still match up each time this script is run!


for(i in 1:length(Schindler_totalCount_epi_final$DepthOfTow_m[1:2])){
  Density.neteff_epi$Actual_dens[i] <-  Schindler_totalCount_epi_final$Zooplankton_No.[i]/ (30*11)
}

#manually calculate actual density for the 2 fcr samples becuase there are 10 (not 11) schindler depths
for(i in 3:length(Schindler_totalCount_epi_final$DepthOfTow_m)){
  Density.neteff_epi$Actual_dens[i] <-  Schindler_totalCount_epi_final$Zooplankton_No.[i]/ (30*10)
}

#### Net Efficiency = APPARENT density / ACTUAL density ####
for(i in 1:length(Density.neteff$sample_ID)){
  Density.neteff$NetEff[i] <- Density.neteff$Apparent_dens[i]/Density.neteff$Actual_dens[i]
}

for(i in 1:length(Density.neteff_epi$sample_ID)){
  Density.neteff_epi$NetEff[i] <- Density.neteff_epi$Apparent_dens[i]/Density.neteff_epi$Actual_dens[i]
}

#add depths for epi net efficiency
Density.neteff_epi$Depth_m <- Schindler_totalCount_epi_final$DepthOfTow_m


#Going to take the average net efficiency across both reps (and average noon/midnight) because values are pretty close to each other
NetEfficiency2020 <- c(mean(Density.neteff$NetEff[1:6]),mean(Density.neteff$NetEff[8:11])) #BVR, FCR
  
#---------------------------#
#visualize density by depth #
#---------------------------#
#order depth by decreasing number
Schindler_totalCount <- Schindler_totalCount[with(Schindler_totalCount,order(DepthOfTow_m)),]

#jpeg("Figures/Schindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_totalCount,aes(x=Zooplankton_No./30, y=DepthOfTow_m,color=Hour)) + geom_point() +
  scale_y_reverse() + geom_path() + facet_grid(~site_no+collect_date)
#dev.off()

#summarize schindler_totalCount so one # per depth
Schindler_avgCount <- Schindler_totalCount %>% group_by(site_no, DepthOfTow_m, collect_date) %>% summarise(mean_num = mean(Zooplankton_No.))

#jpeg("Figures/AvgSchindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_avgCount,aes(x=mean_num, y=DepthOfTow_m)) + geom_point() +
  scale_y_reverse() + geom_path() + xlab("Zooplankton (#)") + ylab("Depth (m)")  + facet_grid(~site_no+collect_date)
#dev.off()

