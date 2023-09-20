# DVM calcs for all five 24-hr sampling campaigns 

### Description of data --> full water column tows and epi tows from BVR 
    #summer 2019 (Jun 10-11); includes samples collected at macrophyes (BVR_l), 
    #dam (BVR_d) and pelagic site (BVR_50_p for epi tows collected during MSN ONLY;
    #BVR_50 for full water column tows and tows outside of 24-hour campaigns)
    #samples collected from at noon (x1), midnight (x1), sunset (x4), and sunrise (x4)

#Load packages and functions
source("scripts/install.R")

#list of years to loop through
year <- c(2019, 2020, 2021)

for(i in 1:length(year)){
  
#averaged from 2020+2021 data
bvr_neteff<- 0.02056477

#for loop to fill out epi vs hypo calcs 
#hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
#NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. 
column.names<- colnames(BVR_pelagic_DVM_vol_calculated[,grepl("pL_rep.mean",
                        colnames(BVR_pelagic_DVM_vol_calculated)) |
                                              grepl("ugpL_rep.mean",
                        colnames(BVR_pelagic_DVM_vol_calculated))])
  
variables<- colnames(BVR_pelagic_DVM_raw[,grepl("Count_n_rep.mean",
                                      colnames(BVR_pelagic_DVM_raw)) |
                                        grepl("ug_rep.mean",
                                              colnames(BVR_pelagic_DVM_raw))])

percent<- colnames(BVR_pelagic_DVM_percent[,grepl("PercentOfTotal_rep.mean",
                                        colnames(BVR_pelagic_DVM_percent))])

for(r in 1:length(variables)){
  BVR.DVM.calcs[,paste0(column.names,"_epi")[r]]<- 
    BVR_pelagic_DVM_vol_calculated[substrEnd(
    BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[r]]
  
  BVR.DVM.calcs[,paste0(column.names,"_hypo")[r]] <- 
    (((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi"  ,
          "proportional_vol"]) * BVR_pelagic_DVM_raw[
      substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,
      paste0(variables)[r]] * (1/bvr_neteff)) - 
      ((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",
          "proportional_vol"]) *  BVR_pelagic_DVM_raw[
      substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",
      paste0(variables)[r]]))/ #no neteff needed for epi tows bc pretty much 100% efficient
      (BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,
           "Volume_unadj"] - BVR_pelagic_DVM_raw[
             substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "Volume_unadj"])  
}

density.percent<- colnames(BVR_pelagic_DVM_vol_calculated[,grepl("NopL_rep.mean",
                            colnames(BVR_pelagic_DVM_vol_calculated))])

for(p in 1:length(density.percent)){
for(j in 1:length(unique(BVR.DVM.calcs$Hour))){
    BVR.DVM.calcs[j,paste0(density.percent,"_epi_percent_density")[p]]<- 
      (BVR.DVM.calcs[j,paste0(density.percent,"_epi")][p] / 
         sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[p]],
             BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[p]])) *100
   
     BVR.DVM.calcs[j,paste0(density.percent,"_hypo_percent_density")[p]]<- 
       (BVR.DVM.calcs[j,paste0(density.percent,"_hypo")][p] / 
          sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[p]],
              BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[p]])) * 100
  }       
 }

#pull only noon/midnight samples
DVM_samples_raw <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | 
                           substrEnd(zoop$sample_ID,5)=="night" | 
                           substrEnd(zoop$sample_ID,8)=="noon_epi" | 
                           substrEnd(zoop$sample_ID,9)=="night_epi") & 
                          zoop$site_no!="BVR_l",]

matchingcols <- match(substr(c("sample_ID","site_no","Hour","collect_date",
                               variables),1,14),
                      substr(colnames(DVM_samples_raw),1,14))

DVM_samples_raw<- DVM_samples_raw[,unique(matchingcols)]

DVM_samples_dens <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | 
                            substrEnd(zoop$sample_ID,5)=="night" | 
                            substrEnd(zoop$sample_ID,8)=="noon_epi" | 
                            substrEnd(zoop$sample_ID,9)=="night_epi") & 
                           zoop$site_no!="BVR_l",]

matchingcols <- match(substr(c("sample_ID","site_no","Hour","collect_date",
                               density.percent) ,1,14),
                      substr(colnames(DVM_samples_dens),1,14))

DVM_samples_dens<- DVM_samples_dens[,unique(matchingcols)] 

#select full water column tows with at least one other rep
if(year[i] == 2019){
       FullSamples<- c("B_pel_11Jul19_midnight","B_pel_10Jul19_noon") 
       
       EpiSamples <- c("B_pel_11Jul19_midnight_epi","B_pel_10Jul19_noon_epi")
       
  }else if(year[i] == 2020) {
         
       FullSamples<- c("B_pel_12Aug20_noon", "B_pel_13Aug20_midnight",
                       "B_pel_13Aug20_noon") 
       
       EpiSamples <- c("B_pel_12Aug20_noon_epi", "B_pel_13Aug20_midnight_epi",
                       "B_pel_13Aug20_noon_epi")
 
  }else {
       
       FullSamples<- c("B_pel_16Jun21_midnight", "B_pel_16Jun21_noon",
                       "B_pel_07Jul21_noon", "B_pel_08Jul21_midnight",
                       "B_pel_08Jul21_noon") 
         
         EpiSamples <- c("B_pel_16Jun21_midnight_epi", "B_pel_16Jun21_noon_epi",
                         "B_pel_07Jul21_noon_epi", "B_pel_08Jul21_midnight_epi",
                         "B_pel_08Jul21_noon_epi")
       
}


#initialize df
BVR.DVM.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))

#can only calculate hypo SE when there are reps of epi and full samples
SEonly<- colnames(DVM_samples_raw)[grepl("t_n", colnames(DVM_samples_raw)) |
                                 grepl("_ug", colnames(DVM_samples_raw))]

Percentdens <- colnames(DVM_samples_dens)[grepl("NopL", 
                                                colnames(DVM_samples_dens))]

#loop to calculate epi and hypo SE
for(k in 1:length(FullSamples)){
  for(l in 1:length(SEonly)){
    BVR.DVM.calcs.SE[,paste0(column.names,"_epi_SE")[l]] <- 
      BVR_pelagic_DVM_vol_calculated[substrEnd(
        BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",
        substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),2)=="SE"][l]
    
  #averaged from 2020+2021 data
  bvr_neteff<- 0.04112955
  
  #for loop to fill out epi vs hypo calcs 
  #hypo density and biomass calculated by subtracting epi raw zoop # from full zoop # and then dividing by the (full volume - epi volume) 
  #NOTE: using epi density/L and biomass/L but calculating hypo using raw # and ug values. 
  column.names<- colnames(BVR_pelagic_DVM_vol_calculated[,grepl("pL_rep.mean",
                          colnames(BVR_pelagic_DVM_vol_calculated)) |
                                                grepl("ugpL_rep.mean",
                          colnames(BVR_pelagic_DVM_vol_calculated))])
    
  variables<- colnames(BVR_pelagic_DVM_raw[,grepl("Count_n_rep.mean",
                                        colnames(BVR_pelagic_DVM_raw)) |
                                          grepl("ug_rep.mean",
                                                colnames(BVR_pelagic_DVM_raw))])
  
  percent<- colnames(BVR_pelagic_DVM_percent[,grepl("PercentOfTotal_rep.mean",
                                          colnames(BVR_pelagic_DVM_percent))])
  
  for(r in 1:length(variables)){
    BVR.DVM.calcs[,paste0(column.names,"_epi")[r]]<- 
      BVR_pelagic_DVM_vol_calculated[substrEnd(
      BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",paste0(column.names)[r]]
    
    BVR.DVM.calcs[,paste0(column.names,"_hypo")[r]] <- 
      (((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi"  ,
            "proportional_vol"]) * BVR_pelagic_DVM_raw[
        substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,
        paste0(variables)[r]] * (1/bvr_neteff)) - 
        ((1/BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",
            "proportional_vol"]) *  BVR_pelagic_DVM_raw[
        substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi",
        paste0(variables)[r]]))/ #no neteff needed for epi tows bc pretty much 100% efficient
        (BVR_pelagic_DVM_raw[substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)!="epi" ,
             "Volume_unadj"] - BVR_pelagic_DVM_raw[
               substrEnd(BVR_pelagic_DVM_raw$sample_ID,3)=="epi", "Volume_unadj"])  
  }
  
  density.percent<- colnames(BVR_pelagic_DVM_vol_calculated[,grepl("NopL_rep.mean",
                              colnames(BVR_pelagic_DVM_vol_calculated))])
  
  for(p in 1:length(density.percent)){
  for(j in 1:length(unique(BVR.DVM.calcs$Hour))){
      BVR.DVM.calcs[j,paste0(density.percent,"_epi_percent_density")[p]]<- 
        (BVR.DVM.calcs[j,paste0(density.percent,"_epi")][p] / 
           sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[p]],
               BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[p]])) *100
     
       BVR.DVM.calcs[j,paste0(density.percent,"_hypo_percent_density")[p]]<- 
         (BVR.DVM.calcs[j,paste0(density.percent,"_hypo")][p] / 
            sum(BVR.DVM.calcs[j,paste0(density.percent,"_epi")[p]],
                BVR.DVM.calcs[j,paste0(density.percent,"_hypo")[p]])) * 100
    }       
   }
  
  #pull only noon/midnight samples
  DVM_samples_raw <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | 
                             substrEnd(zoop$sample_ID,5)=="night" | 
                             substrEnd(zoop$sample_ID,8)=="noon_epi" | 
                             substrEnd(zoop$sample_ID,9)=="night_epi") & 
                            zoop$site_no!="BVR_l",]
  
  matchingcols <- match(substr(c("sample_ID","site_no","Hour","collect_date",
                                 variables),1,14),
                        substr(colnames(DVM_samples_raw),1,14))
  
  DVM_samples_raw<- DVM_samples_raw[,unique(matchingcols)]
  
  DVM_samples_dens <- zoop[(substrEnd(zoop$sample_ID,4)=="noon" | 
                              substrEnd(zoop$sample_ID,5)=="night" | 
                              substrEnd(zoop$sample_ID,8)=="noon_epi" | 
                              substrEnd(zoop$sample_ID,9)=="night_epi") & 
                             zoop$site_no!="BVR_l",]
  
  matchingcols <- match(substr(c("sample_ID","site_no","Hour","collect_date",
                                 density.percent) ,1,14),
                        substr(colnames(DVM_samples_dens),1,14))
  
  DVM_samples_dens<- DVM_samples_dens[,unique(matchingcols)] 
  
  #select full water column tows with at least one other rep
  if(year[i] == 2019){
         FullSamples<- c("B_pel_11Jul19_midnight","B_pel_10Jul19_noon") 
         
         EpiSamples <- c("B_pel_11Jul19_midnight_epi","B_pel_10Jul19_noon_epi")
         
    }else if(year[i] == 2020) {
           
         FullSamples<- c("B_pel_12Aug20_noon", "B_pel_13Aug20_midnight",
                         "B_pel_13Aug20_noon") 
         
         EpiSamples <- c("B_pel_12Aug20_noon_epi", "B_pel_13Aug20_midnight_epi",
                         "B_pel_13Aug20_noon_epi")
   
    }else {
         
         FullSamples<- c("B_pel_16Jun21_midnight", "B_pel_16Jun21_noon",
                         "B_pel_07Jul21_noon", "B_pel_08Jul21_midnight",
                         "B_pel_08Jul21_noon") 
           
           EpiSamples <- c("B_pel_16Jun21_midnight_epi", "B_pel_16Jun21_noon_epi",
                           "B_pel_07Jul21_noon_epi", "B_pel_08Jul21_midnight_epi",
                           "B_pel_08Jul21_noon_epi")
         
  }
  
  
  #initialize df
  BVR.DVM.calcs.SE<- data.frame("Hour"=unique(BVR_pelagic_DVM_raw$Hour))
  
  #can only calculate hypo SE when there are reps of epi and full samples
  SEonly<- colnames(DVM_samples_raw)[grepl("t_n", colnames(DVM_samples_raw)) |
                                   grepl("_ug", colnames(DVM_samples_raw))]
  
  Percentdens <- colnames(DVM_samples_dens)[grepl("NopL", 
                                                  colnames(DVM_samples_dens))]
  
  #loop to calculate epi and hypo SE
  for(k in 1:length(FullSamples)){
    for(l in 1:length(SEonly)){
      BVR.DVM.calcs.SE[,paste0(column.names,"_epi_SE")[l]] <- 
        BVR_pelagic_DVM_vol_calculated[substrEnd(
          BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",
          substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),2)=="SE"][l]
      
      BVR.DVM.calcs.SE[k,paste0(column.names,"_hypo_SE")[l]]<- 
        SE.diffMean(DVM_samples_raw[DVM_samples_raw$sample_ID==FullSamples[k],
                                    paste0(SEonly)[l]],
                    DVM_samples_raw[DVM_samples_raw$sample_ID==EpiSamples[k],
                                    paste0(SEonly)[l]])
      
    }}
  
  #loop to calculate epi and hypo percent density SE
  for(k in 1:length(FullSamples)){
    for (l in 1:length(Percentdens)){
      BVR.DVM.calcs.SE[,paste0(Percentdens,"_epi_percent_density_SE")[l]] <- 
        BVR_pelagic_DVM_vol_calculated[substrEnd(
          BVR_pelagic_DVM_vol_calculated$sample_ID,3)=="epi",
          substrEnd(colnames(BVR_pelagic_DVM_vol_calculated),11)=="NopL_rep.SE"][l]
      
      BVR.DVM.calcs.SE[k,paste0(Percentdens,"_hypo_percent_density_SE")[l]] <- 
        SE.diffMean(DVM_samples_dens[DVM_samples_dens$sample_ID==FullSamples[k],
                                     paste0(Percentdens)[l]],
                    DVM_samples_dens[DVM_samples_dens$sample_ID==EpiSamples[k],
                                     paste0(Percentdens)[l]])
    }
  }
  
  
  #wide to long for both dfs separately
  BVR.DVM.calcs.long <-  BVR.DVM.calcs %>% gather(metric,value, 
    ZoopDensity_No.pL_rep.mean_epi:nauplius_density_NopL_rep.mean_hypo_percent_density) %>%
      mutate(DateTime = strftime(Hour, "%m-%d-%Y %H:%M"))
  
  BVR.DVM.calcs.SE.long <-  BVR.DVM.calcs.SE %>% gather(taxa.metric,SE,
    ZoopDensity_No.pL_rep.mean_epi_SE:nauplius_density_NopL_hypo_percent_density_SE)
  
  #add the SE column from df2 to df1 for combined df
  BVR.DVM.calcs.long$SE <- BVR.DVM.calcs.SE.long$SE
  
  #add watercolumn and hour columns
  BVR.DVM.calcs.long$WaterColumn <- ifelse(
    substrEnd(BVR.DVM.calcs.long$metric,3)=="epi" |
      substrEnd(BVR.DVM.calcs.long$metric,19)=="epi_percent_density" ,
    "epilimnion","hypolimnion")
  
  BVR.DVM.calcs.long$Hour <- ifelse(substrEnd(
    BVR.DVM.calcs.long$DateTime,5)=="12:00","noon","midnight")
  
  BVR.DVM.calcs.long$Taxa <- substr(BVR.DVM.calcs.long$metric,1,9)
  
  #shorten date time to just date
  BVR.DVM.calcs.long$DateTime <- 
    substr(BVR.DVM.calcs.long$DateTime,1,nchar(BVR.DVM.calcs.long$DateTime)-6)
  
  #export 2019 dvm stats
  write.csv(BVR.DVM.calcs.long,paste0("output/DVM_",year[i], "_zoops.csv"),
            row.names = FALSE)

}
