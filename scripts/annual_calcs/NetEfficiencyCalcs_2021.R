#script for calculating zooplankton net efficiency in Beaverdam Reservoir
#Created: 20 November 2020
#for some reason, did not collect full water column tow at the same time as schindlers (??!!)
#options are to either average noon and midnight full densities or just use one or the other
#updated for 2021 samples 22Aug22

#load in packages
pacman::p_load(dplyr,ggplot2)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data 
zoop<- read.csv('output/FCR_ZooplanktonSummary2021.csv',header = TRUE)

#select cols to keep
keep <- c("sample_ID","site_no","collect_date","Hour","DepthOfTow_m","Zooplankton_No.","Volume_L","Volume_unadj","mesh_size_μm")

zoop_totalcount <- zoop[!is.na(zoop$mesh_size_μm) & zoop$site_no!="BVR_l" & zoop$site_no!="BVR_trap",keep]

#convert hour to character
zoop_totalcount$Hour <- as.character(zoop_totalcount$Hour)

#round to nearest hour
zoop_totalcount$Hour <- format(strptime(zoop_totalcount$Hour,format="%H"),"%H")
zoop_totalcount$Hour <- ifelse(zoop_totalcount$Hour=="23"| zoop_totalcount$Hour=="01", "00", 
                               ifelse(zoop_totalcount$Hour=="11" |zoop_totalcount$Hour=="13"| zoop_totalcount$Hour=="14", "12", 
                                      ifelse(zoop_totalcount$Hour=="03","04" ,zoop_totalcount$Hour)))

#replace vol_unadj column with NA for schindler traps
zoop_totalcount$Volume_unadj[zoop_totalcount$mesh_size_μm ==61] <- NA

#------------------------------------------------------------------------------#
#### BVR and FCR 2021 net efficiency ####

#separate tow data from schindler data and epi tows in separate dfs
Schindler_totalCount <- zoop_totalcount[zoop_totalcount$mesh_size_μm==61,]
Tow_totalCount_final <- zoop_totalcount[zoop_totalcount$mesh_size_μm==80 & zoop_totalcount$DepthOfTow_m>=9,]
Epi_tow_totalCount_final <- zoop_totalcount[zoop_totalcount$mesh_size_μm==80 & zoop_totalcount$DepthOfTow_m<9,]

#add a column for rep
Schindler_totalCount$Rep <- substrEnd(Schindler_totalCount$sample_ID,1)

#sum total count for rep1 vs rep2 for each reservoir ot each of the different sampling times
Schindler_totalCount_final <- Schindler_totalCount %>% select(collect_date, site_no,Hour,DepthOfTow_m,Zooplankton_No.,Rep, Volume_L,Volume_unadj) %>%
  group_by(site_no,collect_date,Hour,Rep)  %>% summarise(DepthOfTow_m=max(DepthOfTow_m),Volume_L=mean(Volume_L),TotalCount_n=sum(Zooplankton_No.))


#new df based on epi tows
Schindler_totalCount_epi_final <- Epi_tow_totalCount_final[,c(1:5,7,8)]

#add rep column
Schindler_totalCount_epi_final$Rep <- substrEnd(Schindler_totalCount_epi_final$sample_ID,1)

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
#count / volume_unadj (because we want to whole total volume not volume counted)
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
  Density.neteff$Actual_dens[1] <-  Schindler_totalCount_final$TotalCount_n[4] / (30*11)
  Density.neteff$Actual_dens[2] <-  Schindler_totalCount_final$TotalCount_n[5] / (30*11)
  Density.neteff$Actual_dens[3:4] <- Schindler_totalCount_final$TotalCount_n[6] / (30*11) #only one schindler rep was collected on this date at midnight
  Density.neteff$Actual_dens[5] <- Schindler_totalCount_final$TotalCount_n[4] / (30*11) #using noon schindlers from previous day because not collected this day
  Density.neteff$Actual_dens[6] <- Schindler_totalCount_final$TotalCount_n[5] / (30*11) #using noon schindlers from previous day because not collected this day
  Density.neteff$Actual_dens[7] <- Schindler_totalCount_final$TotalCount_n[1] / (30*11)
  Density.neteff$Actual_dens[8] <- Schindler_totalCount_final$TotalCount_n[2] / (30*11)
  Density.neteff$Actual_dens[9] <- Schindler_totalCount_final$TotalCount_n[3] / (30*11)#only one midnight schindler taken so using for both reps     
  Density.neteff$Actual_dens[10] <- Schindler_totalCount_final$TotalCount_n[3] / (30*11) #only one midnight schindler taken so using for both reps     
  Density.neteff$Actual_dens[11] <-  Schindler_totalCount_final$TotalCount_n[1] / (30*11)#using noon schindlers from previous day because not collected this day
  Density.neteff$Actual_dens[12] <-  Schindler_totalCount_final$TotalCount_n[2] / (30*11)#using noon schindlers from previous day because not collected this day
  
  Density.neteff$Actual_dens[13] <-  Schindler_totalCount_final$TotalCount_n[11] / (30*10)
  Density.neteff$Actual_dens[14] <-  Schindler_totalCount_final$TotalCount_n[12] / (30*10)
  Density.neteff$Actual_dens[15] <-  Schindler_totalCount_final$TotalCount_n[13] / (30*10)
  Density.neteff$Actual_dens[16] <-  Schindler_totalCount_final$TotalCount_n[14] / (30*10)
  
  Density.neteff$Actual_dens[17] <-  Schindler_totalCount_final$TotalCount_n[9] / (30*10)
  Density.neteff$Actual_dens[18] <-  Schindler_totalCount_final$TotalCount_n[10] / (30*10)
  Density.neteff$Actual_dens[19] <-  Schindler_totalCount_final$TotalCount_n[7] / (30*10)
  Density.neteff$Actual_dens[20] <-  Schindler_totalCount_final$TotalCount_n[8] / (30*10)
  
  #NOTE - need to check this every time because I'm annoying and didn't make this robust/able to deal with changes to data or rows switching positions
  
  for(i in 1:length(Schindler_totalCount_epi_final$DepthOfTow_m[1:10])){
  Density.neteff_epi$Actual_dens[i] <-  Schindler_totalCount_epi_final$Zooplankton_No.[i]/ (30*11)
  }
  
  #manually calculate actual density for the 4 fcr samples becuase there are 10 (not 11) schindlers collected
  for(i in 11:length(Schindler_totalCount_epi_final$DepthOfTow_m)){
  Density.neteff_epi$Actual_dens[i] <-  Schindler_totalCount_epi_final$Zooplankton_No.[i]/ (30*10)
  }
  
#### Net Efficiency = APPARENT density / ACTUAL density ####
for(i in 1:length(Density.neteff$sample_ID)){
  Density.neteff$NetEff[i] <- Density.neteff$Apparent_dens[i]/Density.neteff$Actual_dens[i]
}
  
  for(i in 1:length(Density.neteff_epi$sample_ID)){
    Density.neteff_epi$NetEff[i] <- Density.neteff_epi$Apparent_dens[i]/Density.neteff_epi$Actual_dens[i]
  }

  #add depths for epi net efficiencies
  Density.neteff_epi$Depth_m <- Schindler_totalCount_epi_final$DepthOfTow_m
  
#Going to take the average net efficiency across both reps because values are super close to each other
NetEfficiency2021 <- c(mean(Density.neteff$NetEff[1:12]),mean(Density.neteff$NetEff[c(13:20)])) #BVR, FCR

epi_neteff_2021 <- mean(Density.neteff_epi$NetEff[15:18]) #just for FCR

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
Schindler_avgCount <- Schindler_totalCount %>% group_by(site_no, collect_date, DepthOfTow_m) %>% summarise(mean_num = mean(Zooplankton_No.))

#jpeg("Figures/AvgSchindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_avgCount,aes(x=mean_num/30, y=DepthOfTow_m)) + geom_point() +
  scale_y_reverse() + geom_path() + xlab("Zooplankton (#/L)") + ylab("Depth (m)")  + facet_grid(~site_no+collect_date)
#dev.off()

