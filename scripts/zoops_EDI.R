#script for getting all summary files into EDI format

source("scripts/01_install.R")

#read in zoop data from all 3 years
zoops2019<- read.csv('output/FCR_ZooplanktonSummary2019.csv',header = TRUE)
zoops2020<- read.csv('output/FCR_ZooplanktonSummary2020.csv',header = TRUE)
zoops2021<- read.csv('output/FCR_ZooplanktonSummary2021.csv',header = TRUE)
#all_zoops <- read.csv('output/FCR_ZooplanktonSummary2015_biomassONLY.csv',header=TRUE) |> 
#  select(!contains("Calanoida"))


#merge all files
all_zoops <- bind_rows(zoops2019, zoops2020, zoops2021)

#remove horizontal trap data
all_zoops <- all_zoops[!c(grepl("top", all_zoops$sample_ID) | 
                           grepl("bot", all_zoops$sample_ID)),]

#remove 20um samples
all_zoops <- all_zoops[all_zoops$mesh_size_μm!=20,]

#drop percent of total cols, count cols, and se cols
all_zoops <- all_zoops[,!c(grepl("Percent", colnames(all_zoops)) |
                             grepl("Count", colnames(all_zoops)) |
                             grepl("SE", colnames(all_zoops)))]

#drop the total zoop dens/biomass cols
all_zoops <- all_zoops |> select(-c(ZoopDensity_No.pL,OverallMeanSize_mm,
                                    TotalBiomass_ug, BiomassConcentration_ugpL))

#add rep column
all_zoops$Rep <- ifelse(substrEnd(all_zoops$sample_ID,4)=="rep1" |
                          substrEnd(all_zoops$sample_ID,4)=="rep2" | 
                          substrEnd(all_zoops$sample_ID,4)=="rep3"| 
                          substrEnd(all_zoops$sample_ID,4)=="rep4",
                        substrEnd(all_zoops$sample_ID,1),1)
  

#now convert from wide to long for dens
all_zoops_dens <- all_zoops |> 
  select(sample_ID:DepthOfTow_m, Rep,
         colnames(all_zoops[grepl("density",colnames(all_zoops))])) |> 
  pivot_longer(cols = colnames(all_zoops[grepl("density", colnames(all_zoops))]), 
               names_to = "Taxon", values_to = "Density_IndPerL") |> 
  arrange(site_no,collect_date,Hour,DepthOfTow_m,Rep,Taxon)

#convert from wide to long for size
all_zoops_size <- all_zoops |> 
  select(sample_ID:DepthOfTow_m, Rep,
         colnames(all_zoops[grepl("MeanSize",colnames(all_zoops))])) |> 
  pivot_longer(cols = colnames(all_zoops[grepl("MeanSize", colnames(all_zoops))]), 
               names_to = "Taxon", values_to = "MeanLength_mm") |> 
  arrange(site_no,collect_date,Hour,DepthOfTow_m,Rep,Taxon)

#convert from wide to long for weight
all_zoops_weight <- all_zoops |> 
  select(sample_ID:DepthOfTow_m, Rep,
         colnames(all_zoops[grepl("totalbiomass_ug",colnames(all_zoops))])) |> 
  pivot_longer(cols = colnames(all_zoops[grepl("totalbiomass_ug", colnames(all_zoops))]), 
               names_to = "Taxon", values_to = "MeanWeight_ug") |> 
  arrange(site_no,collect_date,Hour,DepthOfTow_m,Rep,Taxon)

#convert from wide to long for biomass
all_zoops_biom <- all_zoops |> 
  select(sample_ID:DepthOfTow_m, Rep,
         colnames(all_zoops[grepl("BiomassConcentration",colnames(all_zoops))])) |> 
  pivot_longer(cols = colnames(all_zoops[grepl("BiomassConcentration", 
               colnames(all_zoops))]), names_to = "Taxon",
               values_to = "Biomass_ugL") |> 
  arrange(site_no,collect_date,Hour,DepthOfTow_m,Rep,Taxon)

#combine dfs
all_zoops_final <- all_zoops_dens
all_zoops_final$MeanLength_mm <- all_zoops_size$MeanLength_mm
all_zoops_final$MeanWeight_ug <- all_zoops_weight$MeanWeight_ug
all_zoops_final$Biomass_ugL <- all_zoops_biom$Biomass_ugL

#add reservoir, site, and datetime cols
all_zoops_final$Reservoir <- ifelse(substr(all_zoops_final$sample_ID,1,1) == 
                                      "B", "BVR", "FCR")
all_zoops_final$Site <- ifelse(grepl("dam", all_zoops_final$sample_ID), 49,
                               ifelse(grepl("mac", all_zoops_final$sample_ID), 
                                      51, 50))
all_zoops_final$DateTime <- as.POSIXct(paste(all_zoops_final$collect_date, 
                                       all_zoops_final$Hour), 
                                       format = "%Y-%m-%d %H:%M")

#drop density_nopl from end of taxon
all_zoops_final$Taxon <- substr(all_zoops_final$Taxon,1,
                                nchar(all_zoops_final$Taxon)-13)

#add collection method col
all_zoops_final$CollectionMethod <- ifelse(grepl("schind", 
                                all_zoops_final$sample_ID), "Schindler", "Tow")

#rename depth of tow to StartDepth_m
all_zoops_final <- all_zoops_final |> rename(StartDepth_m = DepthOfTow_m)

#end depth is 0 for tows and start depth for schindlers 
all_zoops_final$EndDepth_m <- ifelse(all_zoops_final$CollectionMethod=="Tow", 0,
                                     all_zoops_final$StartDepth_m)

#drop unnecessary cols
all_zoops_final <- all_zoops_final |> select(-c(sample_ID, Project, site_no, 
                                                collect_date, Hour))

#change order of cols
all_zoops_final <- all_zoops_final |> select(Reservoir, Site, DateTime, 
                                             StartDepth_m, EndDepth_m, Rep,
                                             CollectionMethod, Taxon, 
                                             Density_IndPerL, MeanLength_mm,
                                             MeanWeight_ug, Biomass_ugL)

#if density is not 0, but weight/biomass is 0, set to NA
all_zoops_final$MeanWeight_ug <- ifelse(all_zoops_final$Density_IndPerL!=0 & 
                                          all_zoops_final$MeanWeight_ug==0, 
                                        NA, all_zoops_final$MeanWeight_ug)

all_zoops_final$Biomass_ugL <- ifelse(all_zoops_final$Density_IndPerL!=0 & 
                                        all_zoops_final$Biomass_ugL==0, 
                                        NA, all_zoops_final$Biomass_ugL)

#add flag cols
all_zoops_final$Flag_Length <- ifelse(is.na(all_zoops_final$MeanLength_mm), 1, 0)
all_zoops_final$Flag_Weight <- ifelse(is.na(all_zoops_final$MeanWeight_ug), 1, 0)
all_zoops_final$Flag_Biomass <- ifelse(is.na(all_zoops_final$Biomass_ugL), 1, 0)

#order by reservoir, site, date, then collection method
all_zoops_final <- all_zoops_final |> 
  arrange(Reservoir, DateTime, Site, CollectionMethod)

#convert datetime to character while keeping the 00 hour for midnight samples
all_zoops_final$DateTime <- 
  ifelse(hour(all_zoops_final$DateTime)==0,
    paste0(as.Date(all_zoops_final$DateTime), " 00:00:01"),
    as.character(all_zoops_final$DateTime))

#convert back to posixct
all_zoops_final$DateTime <- as.POSIXct(all_zoops_final$DateTime,
                                       format = "%Y-%m-%d %H:%M:%S", tz = "EST")

#export file 
write.csv(all_zoops_final, 'output/EDI_zoop_taxa_2019-2022.csv', row.names = FALSE)

#-------------------------------------------------------------------------------#
#raw density data sheet for EDI

#read in density csv from all 3 years
zoopdens19 <- read.csv('zoop_data/FCR2019_ZooplanktonCounting_Density_DataEntry.csv',
                       header = TRUE)
zoopdens20 <- read.csv('zoop_data/FCR2020_ZooplanktonCounting_Density_DataEntry.csv',
                       header = TRUE)
zoopdens21 <- read.csv('zoop_data/FCR2021_ZooplanktonCounting_Density_DataEntry.csv',
                       header = TRUE)

#combine files
zoopdens <- rbind(zoopdens19, zoopdens20, zoopdens21)

#drop 20um samples and horizontal trap samples
zoopdens <- zoopdens[!c(is.na(zoopdens$mesh_size_μm) | zoopdens$mesh_size_μm==20),]

#drop other test sample that doesn't have an hour
zoopdens <- zoopdens[zoopdens$Hour!="",]

#add reservoir
zoopdens$Reservoir <- ifelse(substr(zoopdens$sample_ID,1,1) == "B", "BVR", "FCR")

#add site
zoopdens$Site <- zoopdens$Site <- ifelse(grepl("dam", zoopdens$sample_ID), 49,
                                                ifelse(grepl("mac", 
                                                zoopdens$sample_ID), 51, 50))

#add datetime
zoopdens$DateTime <- as.POSIXct(paste(zoopdens$collect_date, zoopdens$Hour), 
                                       format = "%d-%b-%y %H:%M")

#add collection method
zoopdens$CollectionMethod <- ifelse(grepl("schind", zoopdens$sample_ID), 
                                    "Schindler", "Tow")

#start depth
zoopdens <- zoopdens |> rename(StartDepth_m = DepthOfTow_m)

#end depth is 0.1 for tows and start depth for schindlers 
zoopdens$EndDepth_m <- ifelse(zoopdens$CollectionMethod=="Tow", 0.1,
                              zoopdens$StartDepth_m)

#add column for rep # - need this bc not all reps have different collection times recorded...
zoopdens$Rep <- ifelse(grepl("rep",zoopdens$sample_ID), 
                       substrEnd(zoopdens$sample_ID,1), 1)

#drop unnecessary cols
zoopdens <- zoopdens |> select(-c(sample_ID, Project, site_no, collect_date, 
                                  Hour, date_processed, mesh_size_μm, INT, Notes,
                                  PhantomMidgePresence_YesNo))

#change order of cols
zoopdens <- zoopdens |> select(Reservoir, Site, DateTime, StartDepth_m, 
                               EndDepth_m, Rep, CollectionMethod, 
                               InitialSampleVolume_mL, Subsample,
                               SubsampleVolume_mL, Zooplankton_No.)

#export density file
write.csv(zoopdens, 'output/EDI_zoop_raw_dens_2019-2022.csv', row.names = FALSE)

#-------------------------------------------------------------------------------#
#raw biomass data sheet for EDI

#read in density csv from all 3 years
zoopbiom19 <- read.csv('zoop_data/FCR2019_ZooplanktonCounting_SizeID_DataEntry.csv',
                       header = TRUE)
zoopbiom20 <- read.csv('zoop_data/FCR2020_ZooplanktonCounting_SizeID_DataEntry.csv',
                       header = TRUE)
zoopbiom21 <- read.csv('zoop_data/FCR2021_ZooplanktonCounting_SizeID_DataEntry.csv',
                       header = TRUE)

#drop 17jun19 bc was a practice sample
zoopbiom19 <- zoopbiom19[zoopbiom19$collect_date!="17-Jun-19",]

#combine files
zoopbiom <- rbind(zoopbiom19, zoopbiom20, zoopbiom21)

#add column for rep #
zoopbiom$Rep <- ifelse(grepl("rep",zoopbiom$sample_ID), 
                       substrEnd(zoopbiom$sample_ID,1), 1)

#drop 20um samples and horizontal trap samples
zoopbiom <- zoopbiom[!c(grepl("top", zoopbiom$sample_ID) |
                          grepl("bot", zoopbiom$sample_ID) |
                          grepl("_20_", zoopbiom$sample_ID)),]

#add reservoir
zoopbiom$Reservoir <- ifelse(substr(zoopbiom$sample_ID,1,1) == "B", "BVR", "FCR")

#add site
zoopbiom$Site <- zoopbiom$Site <- ifelse(grepl("dam", zoopbiom$sample_ID), 49,
                                  ifelse(grepl("mac", zoopbiom$sample_ID), 51, 50))

#add datetime
zoopbiom$DateTime <- as.POSIXct(paste(zoopbiom$collect_date, zoopbiom$Hour), 
                                format = "%d-%b-%y %H:%M")

#add collection method
zoopbiom$CollectionMethod <- ifelse(grepl("schind", zoopbiom$sample_ID), 
                                    "Schindler", "Tow")

#add start and end depths
zoopbiom$StartDepth_m <- ifelse((grepl("sunset", zoopbiom$sample_ID) |
                                  grepl("sunrise", zoopbiom$sample_ID) |
                                   grepl("epi", zoopbiom$sample_ID)) &
                                  grepl("B_pel", zoopbiom$sample_ID) &
                                  !grepl("F", zoopbiom$sample_ID), 4, 999)

zoopbiom$StartDepth_m <- ifelse(grepl("B_mac", zoopbiom$sample_ID) |
                                grepl("B_dam", zoopbiom$sample_ID), 2, 
                                zoopbiom$StartDepth_m)

zoopbiom$StartDepth_m <- ifelse(grepl("B_pel", zoopbiom$sample_ID) &
                                  !(grepl("sunset", zoopbiom$sample_ID) |
                                   grepl("sunrise", zoopbiom$sample_ID) |
                                   grepl("oxy", zoopbiom$sample_ID) |
                                   grepl("epi", zoopbiom$sample_ID) |
                                   grepl("schind", zoopbiom$sample_ID)), 10, 
                                zoopbiom$StartDepth_m)
  
zoopbiom$StartDepth_m <- ifelse(grepl("schind", zoopbiom$sample_ID), 
                                substr(zoopbiom$sample_ID,
                                       nchar(zoopbiom$sample_ID)-7,
                                       nchar(zoopbiom$sample_ID)-5), 
                                zoopbiom$StartDepth_m)

#now replace the 0.0 with 10
zoopbiom$StartDepth_m <- ifelse(zoopbiom$StartDepth_m=="0.0",10,
                                zoopbiom$StartDepth_m)

#if sample name ends in m, pull that depth
zoopbiom$StartDepth_m <- ifelse(substrEnd(zoopbiom$sample_ID,1)=="m",
                                substrEnd(zoopbiom$sample_ID,5),
                                zoopbiom$StartDepth_m)

#now fix the weird subsetting issues
zoopbiom$StartDepth_m <- ifelse(zoopbiom$StartDepth_m=="_3.5m", 3.5,
                         ifelse(zoopbiom$StartDepth_m=="_9.0m" |
                                  zoopbiom$StartDepth_m=="21_9m" |
                                  zoopbiom$StartDepth_m=="22_9m", 9,
                         ifelse(zoopbiom$StartDepth_m=="_1.5m", 1.5,
                         ifelse(zoopbiom$StartDepth_m=="_3.0m", 3,
                         ifelse(zoopbiom$StartDepth_m=="_8.0m" | 
                                  zoopbiom$StartDepth_m=="19_8m", 8,
                         ifelse(zoopbiom$StartDepth_m=="_2.0m", 2,
                         ifelse(zoopbiom$StartDepth_m=="_2.5m", 2.5,
                         ifelse(zoopbiom$StartDepth_m=="9_10m" |
                                  zoopbiom$StartDepth_m=="0_10m" |
                                  zoopbiom$StartDepth_m=="1_10m", 10,
                         ifelse(zoopbiom$StartDepth_m=="_4.0m", 4,
                         ifelse(zoopbiom$StartDepth_m=="_6.0m", 6,
                         ifelse(zoopbiom$StartDepth_m=="_4.5m", 4.5,
                         ifelse(zoopbiom$StartDepth_m=="0_11m", 11,
                         ifelse(zoopbiom$StartDepth_m=="11.5m", 11.5,
                                zoopbiom$StartDepth_m)))))))))))))

#make all depths numeric
zoopbiom$StartDepth_m <- as.numeric(zoopbiom$StartDepth_m)

#and now do the last ~50 samples manually
zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F01Jul22_noon_epi_rep1",
                                                "F01Jul22_noon_epi_rep2",
                                                "F11Sep20_noon_epi")] <- 2.2

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F01Jul22_midnight_epi_rep1",
                                                "F01Jul22_midnight_epi_rep2",
                                                "F15Sep20_noon_epi")] <- 2.4

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F10Sep20_midnight_epi")] <- 2.5

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F28Jun20_midnight_epi",
                                                "F29Jun20_noon_epi")] <- 2.6

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F14Sep20_midnight_epi")] <- 2.7

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F30Sep20_noon_epi")] <- 4

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F29Jun20_noon_oxy")] <- 3.4

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F27Jul20_noon_epi")] <- 3.6

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F20Jul20_noon_epi")] <- 3.7

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F08Jun20_noon_epi")] <- 4.2

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F10Jun21_noon_epi_rep1",
                                                "F10Jun21_noon_epi_rep2")] <- 4.3

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F10Jun21_midnight_epi_rep1",
                                                "F10Jun21_midnight_epi_rep2")] <- 4.4

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("B_pel_08Jul21_midnight_oxy_rep1",
                                                "B_pel_08Jul21_midnight_oxy_rep2")] <- 5.2

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("B_pel_07Jul21_noon_oxy_rep1",
                                                "B_pel_07Jul21_noon_oxy_rep2",
                                                "B_pel_08Jul21_noon_oxy_rep1",
                                                "B_pel_08Jul21_noon_oxy_rep2")] <- 5.3

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("B_pel_10Jul19_noon_oxy", 
                                                "B_pel_11Jul19_midnight_oxy",
                                                "B_pel_11Jul19_noon_oxy")] <- 5.5

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("B_pel_16Jun21_midnight_oxy_rep1", 
                                                "B_pel_16Jun21_midnight_oxy_rep2")] <- 6.3

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("B_pel_15Jun21_noon_oxy_rep1", 
                                                "B_pel_15Jun21_noon_oxy_rep2")] <- 6.6

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("B_pel_16Jun21_noon_oxy_rep1", 
                                                "B_pel_16Jun21_noon_oxy_rep2")] <- 6.8

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("F01Jul22_noon_rep1", 
                                                "F01Jul22_noon_rep2",
                                                "F01Jul22_midnight_rep1",
                                                "F01Jul22_midnight_rep2",
                                                "F10Jun21_noon_rep1",
                                                "F10Jun21_noon_rep2",
                                                "F10Jun21_midnight_rep1",
                                                "F10Jun21_midnight_rep2",
                                                "F11Sep20_noon",
                                                "F10Sep20_midnight",
                                                "F29Jun20_noon_rep1",
                                                "F29Jun20_noon_rep2",
                                                "F28Jun20_midnight",
                                                "F14Sep20_midnight",
                                                "F15Sep20_noon",
                                                "F08Jun20_noon",
                                                "F30Sep20_noon",
                                                "F27Jul20_noon",
                                                "F20Jul20_noon")] <- 9

zoopbiom$StartDepth_m[zoopbiom$sample_ID %in% c("B_pel_08Jul21_midnight_rep1", 
                                                "B_pel_08Jul21_midnight_rep2",
                                                "B_pel_08Jul21_noon_rep1",
                                                "B_pel_08Jul21_noon_rep2",
                                                "B_pel_07Jul21_noon_rep1",
                                                "B_pel_07Jul21_noon_rep2",
                                                "B_pel_15Jun21_noon_rep1",
                                                "B_pel_15Jun21_noon_rep2",
                                                "B_pel_16Jun21_midnight_rep1",
                                                "B_pel_16Jun21_midnight_rep2",
                                                "B_pel_16Jun21_noon_rep1",
                                                "B_pel_16Jun21_noon_rep2")] <- 10

#end depth
zoopbiom$EndDepth_m <- ifelse(zoopbiom$CollectionMethod=="Tow",0, zoopbiom$StartDepth_m)

#drop unnecessary cols
zoopbiom <- zoopbiom |> select(-c(sample_ID, Project, site_no, collect_date, 
                                  Hour, MarksInOcularMicrometer_Width_No., 
                                  MarksInOcularMicrometer_Height_No., Initials, Notes))
                                  

#add ocular mag column
zoopbiom$OcularMagnification <- "10x"

#change order of cols
zoopbiom <- zoopbiom |> select(Reservoir, Site, DateTime, StartDepth_m, EndDepth_m, Rep,
                               CollectionMethod, Subsample, LowestTaxonomicLevelOfID,
                               TaxaID, Nauplius, ObjectiveMagnification,
                               OcularMagnification, MarksInOcularMicrometer_No.)

#export density file
write.csv(zoopbiom, 'output/EDI_zoop_raw_biom_2019-2022.csv', row.names = FALSE)

