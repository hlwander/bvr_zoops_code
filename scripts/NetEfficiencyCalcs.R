#script for calculating zooplankton net efficiency in Beaverdam Reservoir

#Load packages and functions
source("scripts/install.R")

#only want 2020-2022 bc those are only years with paired schindler and net tows
zoop20<- read.csv('output/FCR_ZooplanktonSummary2020.csv',
                header = TRUE)

zoop21<- read.csv('output/FCR_ZooplanktonSummary2021.csv',
                  header = TRUE)

#combine dfs
zoop <- bind_rows(zoop20, zoop21)

#drop oxycline samples
zoop <- zoop[!grepl("oxy", zoop$sample_ID),]

#select cols to keep
keep <- c("sample_ID","site_no","collect_date","Hour","DepthOfTow_m",
          "Zooplankton_No.","Volume_L","Volume_unadj","mesh_size_μm")

zoop_totalcount <- zoop[!is.na(zoop$mesh_size_μm) & zoop$site_no!="BVR_l" & 
                          zoop$site_no!="BVR_trap",keep]

#convert hour to character
zoop_totalcount$Hour <- as.character(zoop_totalcount$Hour)

#round to nearest hour
zoop_totalcount$Hour <- format(strptime(zoop_totalcount$Hour,format="%H"),"%H")
zoop_totalcount$Hour <- ifelse(zoop_totalcount$Hour=="23"| 
                                 zoop_totalcount$Hour=="01", "00", 
                               ifelse(zoop_totalcount$Hour=="11" | 
                                        zoop_totalcount$Hour=="13"| 
                                        zoop_totalcount$Hour=="14", "12", 
                                      ifelse(zoop_totalcount$Hour=="03","04" ,
                                             zoop_totalcount$Hour)))


#manually change hour for 3 schindlers	
zoop_totalcount$Hour[zoop_totalcount$sample_ID=="B_pel_12Aug20_schind_10.0_rep2"]<- 19	
zoop_totalcount$Hour[zoop_totalcount$sample_ID=="B_pel_12Aug20_schind_8.0_rep2"]<- 19	
zoop_totalcount$Hour[zoop_totalcount$sample_ID=="B_pel_12Aug20_schind_9.0_rep2"]<- 19	

#replace vol_unadj column with NA for schindler traps
zoop_totalcount$Volume_unadj[zoop_totalcount$mesh_size_μm ==61] <- NA

#order by date
zoop_totalcount <-  zoop_totalcount[order(as.Date(zoop_totalcount$collect_date,
                                                 format = "%Y-%m-%d")),]

#------------------------------------------------------------------------------#
#### BVR and FCR net efficiency ####

#separate tow data from schindler data and epi tows in separate dfs
Schindler_totalCount <- zoop_totalcount[zoop_totalcount$mesh_size_μm==61,]
Tow_totalCount_final <- zoop_totalcount[zoop_totalcount$mesh_size_μm==80 & 
                                          zoop_totalcount$DepthOfTow_m>=9,]
Epi_tow_totalCount_final <- zoop_totalcount[zoop_totalcount$mesh_size_μm==80 & 
                                              zoop_totalcount$DepthOfTow_m<9,]

#only keep tows collected on schindler days
Tow_totalCount_final <- Tow_totalCount_final[which(
  Tow_totalCount_final$collect_date %in% Schindler_totalCount$collect_date),]

Epi_tow_totalCount_final <- Epi_tow_totalCount_final[which(
  Epi_tow_totalCount_final$collect_date %in% Schindler_totalCount$collect_date),]

#multiply # of zoops by 2 bc each schindler trap is only ~0.5m
Schindler_totalCount$Zooplankton_No. <-
  Schindler_totalCount$Zooplankton_No. * 2

#add a column for rep
Schindler_totalCount$Rep <- substrEnd(Schindler_totalCount$sample_ID,1)

#sum total count for rep1 vs rep2 for each reservoir ot each of the different sampling times
Schindler_totalCount_final <- Schindler_totalCount %>% select(collect_date, 
                                 site_no,Hour,DepthOfTow_m,Zooplankton_No.,Rep, 
                                 Volume_L,Volume_unadj) %>%
  group_by(collect_date,Hour,Rep,site_no)  %>% 
  summarise(DepthOfTow_m=max(DepthOfTow_m), Volume_L=mean(Volume_L),
            TotalCount_n=sum(Zooplankton_No.))

#new df based on epi tows
Schindler_totalCount_epi_final <- Epi_tow_totalCount_final |> select(sample_ID,
                                      site_no, collect_date, Hour, DepthOfTow_m,
                                      Volume_L, Volume_unadj)

#add rep column
Schindler_totalCount_epi_final$Rep <- 
  substrEnd(Schindler_totalCount_epi_final$sample_ID,1)

#if no rep, then change to 1
Schindler_totalCount_epi_final$Rep <- 
  ifelse(Schindler_totalCount_epi_final$Rep=="i" | 
           Schindler_totalCount_epi_final$Rep=="y" , 1,
         Schindler_totalCount_epi_final$Rep)

#convert date column to date format
Schindler_totalCount_epi_final$collect_date <- 
  as.Date(Schindler_totalCount_epi_final$collect_date)

#loop to sum all depths less than the rounded depth
for(i in 1:length(Schindler_totalCount_epi_final$DepthOfTow_m)){
  Schindler_totalCount_epi_final$Zooplankton_No.[i] <- 
    sum(Schindler_totalCount$Zooplankton_No.[
      which(Schindler_totalCount$DepthOfTow_m <= 
              Schindler_totalCount_epi_final$DepthOfTow_m[i] & 
              as.Date(Schindler_totalCount$collect_date) %in%  
              Schindler_totalCount_epi_final$collect_date[i] &
              Schindler_totalCount$Hour %in% Schindler_totalCount_epi_final$Hour[i] &
              Schindler_totalCount$Rep %in% Schindler_totalCount_epi_final$Rep[i] &
              substr(Schindler_totalCount$site_no,1,1) %in% 
              substr(Schindler_totalCount_epi_final$site_no,1,1)[i])])
}

#remove rows with NA so we only include tows that match up with schindlers
Schindler_totalCount_epi_final <- Schindler_totalCount_epi_final[
    Schindler_totalCount_epi_final$Zooplankton_No. > 0,]

#initialize df
Density.neteff <- data.frame("sample_ID"=Tow_totalCount_final$sample_ID,
                             "collect_date" = Tow_totalCount_final$collect_date,
                             "Hour" = Tow_totalCount_final$Hour,
                             "Rep" = substrEnd(Tow_totalCount_final$sample_ID,1),
                             "site_no" = Tow_totalCount_final$site_no)
Density.neteff_epi <- data.frame("sample_ID"=
                                   Schindler_totalCount_epi_final$sample_ID,
                             "collect_date" = Schindler_totalCount_epi_final$collect_date,
                             "Hour" = Schindler_totalCount_epi_final$Hour,
                             "Rep" = substrEnd(Schindler_totalCount_epi_final$sample_ID,1),
                             "site_no" = Schindler_totalCount_epi_final$site_no)

#if rep is i or n, replace with 1
Density.neteff$Rep <- ifelse(Density.neteff$Rep=="i" | 
                             Density.neteff$Rep=="n" , 1, 
                             Density.neteff$Rep)

Density.neteff_epi$Rep <- ifelse(Density.neteff_epi$Rep=="i" | 
                                 Density.neteff_epi$Rep=="n" , 1, 
                                 Density.neteff_epi$Rep)

#### Calculate APPARENT density from vertical tows ####
#count / volume_unadj (because we want the whole total volume not volume counted)
for(i in 1:length(Density.neteff$sample_ID)){
  Density.neteff$Apparent_dens[i] <-  Tow_totalCount_final$Zooplankton_No.[i] / 
    Tow_totalCount_final$Volume_unadj[i]
}

for(i in 1:length(Density.neteff_epi$sample_ID)){
  Density.neteff_epi$Apparent_dens[i] <-  
    Schindler_totalCount_epi_final$Zooplankton_No.[i] / 
    Schindler_totalCount_epi_final$Volume_unadj[i]
}

#initialize actual dens column
Density.neteff$Actual_dens <- NA
Density.neteff_epi$Actual_dens <- NA

#shorten site name
Density.neteff$site_no <- substr(Density.neteff$site_no,1,3)
Schindler_totalCount_final$site_no <- substr(Schindler_totalCount_final$site_no,1,3)

#set two times to 12 so we have more 2020 neteff calcs (schindlers taken at 1800 and 1900)
Schindler_totalCount_final$Hour <- ifelse(Schindler_totalCount_final$Hour==18 |
                                     Schindler_totalCount_final$Hour==19, 12,
                                     Schindler_totalCount_final$Hour)

#only select rows with paired schindler data (date, hour, site, and rep match)
Density.neteff_paired <- plyr::match_df(Density.neteff, 
                        Schindler_totalCount_final)
  
#### Calculate ACTUAL density (AD) from schindler traps ####
#count / 30L (vol of schindler trap) * # of samples (pooled from all depths)
Density.neteff_paired$Actual_dens <- 
  ifelse(Density.neteff_paired$site_no=="BVR", Density.neteff_paired$Actual_dens <-
         Schindler_totalCount_final$TotalCount_n / (30 * 11), #11 schindlers at bvr
         Schindler_totalCount_final$TotalCount_n / (30 * 10)) #10 schindlers at fcr

#now calculate AD for epi df
  Density.neteff_epi$Actual_dens <-  
    ifelse(Schindler_totalCount_epi_final$site_no=="BVR", 
           Density.neteff_epi$Actual_dens <-
      Schindler_totalCount_epi_final$Zooplankton_No. / (30 * 11), #11 schindlers at bvr
      Schindler_totalCount_epi_final$Zooplankton_No. / (30 * 10)) #10 schindlers at fcr


#### Net Efficiency = APPARENT density / ACTUAL density ####
  Density.neteff_paired$NetEff <- 
    Density.neteff_paired$Apparent_dens / 
    Density.neteff_paired$Actual_dens

    Density.neteff_epi$NetEff <- 
      Density.neteff_epi$Apparent_dens / 
      Density.neteff_epi$Actual_dens

  #add depths for epi net efficiency
  Density.neteff_epi$Depth_m <- Schindler_totalCount_epi_final$DepthOfTow_m
  
#Going to take the average net efficiency across both reps because values are super close to each other
  NetEfficiency <- c(mean(Density.neteff_paired$NetEff[Density.neteff_paired$site_no=="BVR"]),
                       mean(Density.neteff_paired$NetEff[Density.neteff_paired$site_no=="FCR"])) 

  epi_neteff <- c(mean(Density.neteff_epi$NetEff[Density.neteff_epi$site_no=="BVR_50"]),
                mean(Density.neteff_epi$NetEff[Density.neteff_epi$site_no=="FCR_50"])) 

#---------------------------#
#visualize density by depth #
#---------------------------#
#order depth by decreasing number
Schindler_totalCount <- 
  Schindler_totalCount[with(Schindler_totalCount,order(DepthOfTow_m)),]

#jpeg("Figures/Schindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_totalCount,aes(x=Zooplankton_No./30, 
                                     y=DepthOfTow_m,color=Hour)) + geom_point() +
  scale_y_reverse() + geom_path() + facet_grid(~site_no+collect_date)
#dev.off()

#summarize schindler_totalCount so one # per depth
Schindler_avgCount <- Schindler_totalCount %>% 
  group_by(site_no, collect_date, DepthOfTow_m) %>% 
  summarise(mean_num = mean(Zooplankton_No.))

#jpeg("Figures/AvgSchindler_density_vs_depth.jpg", width = 6, height = 5, units = "in",res = 300)
ggplot(data=Schindler_avgCount,aes(x=mean_num/30, y=DepthOfTow_m)) + geom_point() +
  scale_y_reverse() + geom_path() + xlab("Zooplankton (#/L)") + 
  ylab("Depth (m)")  + facet_grid(~site_no+collect_date)
#dev.off()
