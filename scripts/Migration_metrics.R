#Script to calculate DVM and DHM metrics
#Created 5Dec2022

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,lubridate, 
               scales, colorblindcheck, RColorBrewer)

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#read in DHM csv
all_DHM <- read.csv("output/All_MSN_DHM.csv")

#convert DHM hour to character
all_DHM$Hour <- ifelse(grepl("noon",all_DHM$sample_ID),"noon",
                       ifelse(grepl("midnight",all_DHM$sample_ID),"midnight",
                              ifelse(grepl("sunrise",all_DHM$sample_ID),"sunrise","sunset")))

#add taxa column
all_DHM$Taxa <- substr(all_DHM$metric,1,9)

#rename collect_date to DateTime
names(all_DHM)[names(all_DHM)=="collect_date"] <- "DateTime"

#add sunrise and sunset hours
all_DHM$Hour <- ifelse(grepl("sunrise_h1",all_DHM$sample_ID) | grepl("sunrise_epi_h1",all_DHM$sample_ID),"sunrise1",
                       ifelse(grepl("sunrise_h2",all_DHM$sample_ID) | grepl("sunrise_epi_h2",all_DHM$sample_ID),"sunrise2",
                       ifelse(grepl("sunrise_h3",all_DHM$sample_ID) | grepl("sunrise_epi_h3",all_DHM$sample_ID),"sunrise3",
                       ifelse(grepl("sunrise_h4",all_DHM$sample_ID) | grepl("sunrise_epi_h4",all_DHM$sample_ID),"sunrise4",
                       ifelse(grepl("sunset_h1",all_DHM$sample_ID) | grepl("sunset_epi_h1",all_DHM$sample_ID),"sunset1",
                       ifelse(grepl("sunset_h2",all_DHM$sample_ID) | grepl("sunset_epi_h2",all_DHM$sample_ID),"sunset2",
                       ifelse(grepl("sunset_h3",all_DHM$sample_ID) | grepl("sunset_epi_h3",all_DHM$sample_ID),"sunset3",
                       ifelse(grepl("sunset_h4",all_DHM$sample_ID) | grepl("sunset_epi_h4",all_DHM$sample_ID),"sunset4",all_DHM$Hour))))))))

#drop sample_id
all_DHM <- all_DHM %>% select(!sample_ID)

#remove unnecessary letters at the end of each metric name
all_DHM$metric <- substr(all_DHM$metric,1,nchar(all_DHM$metric)-9)

#work with density and biomass to calculate metrics
all_DHM <- all_DHM[!grepl("percent",all_DHM$metric, ignore.case = TRUE),]

#add MSN# column
all_DHM$MSN <- ifelse(all_DHM$DateTime=="2019-07-10" | all_DHM$DateTime=="2019-07-11",1,
                      ifelse(all_DHM$DateTime=="2019-07-24" | all_DHM$DateTime=="2019-07-25",2,
                             ifelse(all_DHM$DateTime=="2020-08-12" | all_DHM$DateTime=="2020-08-13",3,
                                    ifelse(all_DHM$DateTime=="2021-06-15" | all_DHM$DateTime=="2021-06-16",4,5))))

#change second noon for each MSN to noon2
all_DHM$Hour <- ifelse(all_DHM$Hour=="noon" & (all_DHM$DateTime=="2019-07-11" | 
                                                 all_DHM$DateTime=="2019-07-25" |
                                                 all_DHM$DateTime=="2020-08-13" |
                                                 all_DHM$DateTime=="2021-06-16" |
                                                 all_DHM$DateTime=="2021-07-08"), "noon2", all_DHM$Hour)


#split up dfs for noon1 vs noon2 data
all_DHM_noon <- all_DHM[all_DHM$Hour=="noon" | all_DHM$Hour=="midnight",]
all_DHM_noon2 <- all_DHM[all_DHM$Hour=="noon2" | all_DHM$Hour=="midnight",]

#-------------------------------------------------------------------------------
#read in DVM annual csvs
DVM_2021 <- read.csv("output/DVM_2021_zoops.csv") |> select(-X)
DVM_2020 <- read.csv("output/DVM_2020_zoops.csv") |> select(-X)
DVM_2019 <- read.csv("output/DVM_2019_zoops.csv") 

#combine annual DVM dfs
all_DVM <- rbind(DVM_2019,DVM_2020,DVM_2021)

#change date format
all_DVM$DateTime <- format(as.Date(all_DVM$DateTime, "%m-%d-%Y"), "%Y-%m-%d")

#work with raw density and biomass to calculate metrics
all_DVM <- all_DVM[!grepl("percent",all_DVM$metric),]

#add MSN# column
all_DVM$MSN <- ifelse(all_DVM$DateTime=="2019-07-10" | all_DVM$DateTime=="2019-07-11",1,
                            ifelse(all_DVM$DateTime=="2019-07-24" | all_DVM$DateTime=="2019-07-25",2,
                                   ifelse(all_DVM$DateTime=="2020-08-12" | all_DVM$DateTime=="2020-08-13",3,
                                          ifelse(all_DVM$DateTime=="2021-06-15" | all_DVM$DateTime=="2021-06-16",4,5))))

#change second noon for each MSN to noon2
all_DVM$Hour <- ifelse(all_DVM$Hour=="noon" & (all_DVM$DateTime=="2019-07-11" | 
                                                 all_DVM$DateTime=="2019-07-25" |
                                                 all_DVM$DateTime=="2020-08-13" |
                                                 all_DVM$DateTime=="2021-06-16" |
                                                 all_DVM$DateTime=="2021-07-08"), "noon2", all_DVM$Hour)

#not sure if this is the right call, but going to set all negative values to 0 bc they are not real and due to additive error across all steps
all_DVM$value[all_DVM$value < 0] <- 0 #n=20

#shorten metric name
all_DVM$metric <- ifelse(all_DVM$WaterColumn=="epilimnion", substr(all_DVM$metric,1,nchar(all_DVM$metric)-13),substr(all_DVM$metric,1,nchar(all_DVM$metric)-14))

#split up dfs for noon1 vs noon2 data
all_DVM_noon1 <- all_DVM[all_DVM$Hour!="noon2",]
all_DVM_noon2 <- all_DVM[all_DVM$Hour!="noon",]

#NOTE - should Io bring in the proportional volume for epi, hypo, and lit samples???
# BVR total volume = 1357140.624 m3
# epi vol = 806200.056 m3
# hypo vol = 550940.568 m3
# not sure what to do with littoral though - maybe I say the littoral is ~30% of the pelagic epi volumen and then split that one up??

#-------------------------------------------------------------------------------
####                         DVM METRICS                                    ####
#-------------------------------------------------------------------------------
#Calculate proportion of zoops in epi at noon and midnight
DVM_proportion_noon1 <-  plyr::ddply(all_DVM_noon1, c("metric", "MSN", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon1 = x$value[x$WaterColumn=="epilimnion"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

DVM_proportion_noon2 <-  plyr::ddply(all_DVM_noon2, c("metric", "MSN", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon2 = x$value[x$WaterColumn=="epilimnion"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#initialize df for DVM metrics for each day
#migration_df <- data.frame("MSN"= unique(DVM_proportion$MSN))

#DVM metrics for a single day --> DVM = (Depi / Depi + Dhypo)Night - (Depi / Depi + Dhypo)Day
DVM_noon1 <-  plyr::ddply(DVM_proportion_noon1, c("metric", "MSN"), function(x) {
  data.frame(
    DVM_metric_noon1 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


DVM_noon2 <-  plyr::ddply(DVM_proportion_noon2, c("metric", "MSN"), function(x) {
  data.frame(
DVM_metric_noon2 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon2"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#merge both dfs 
DVM_metrics <- right_join(DVM_noon1,DVM_noon2, by=c("metric","MSN"))

#-------------------------------------------------------------------------------
####                         DHM METRICS                                    ####
#-------------------------------------------------------------------------------
#Calculate proportion of zoops in epi at noon and midnight (or sunrise/sunset)
DHM_proportion_noon <-  plyr::ddply(all_DHM_noon, c("metric", "MSN", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon = x$value[x$site_no=="BVR_50"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

DHM_proportion_noon2 <-  plyr::ddply(all_DHM_noon2, c("metric", "MSN", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon2 = x$value[x$site_no=="BVR_50"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#DHM metrics for a single day --> DHM = (Dpelepi / Dpelepi + Dlit)Night - (Dpelepi / Dpelepi + Dlit)Day
DHM_metrics_noon1 <-  plyr::ddply(DHM_proportion_noon, c("metric", "MSN"), function(x) {
  data.frame(
    DHM_metric_noon1 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon"]
    )
}, .progress = plyr::progress_text(), .parallel = FALSE)

DHM_metrics_noon2 <-  plyr::ddply(DHM_proportion_noon2, c("metric", "MSN"), function(x) {
  data.frame(
    DHM_metric_noon2 = x$proportion_epi[x$Hour=="midnight"] - x$proportion_epi[x$Hour=="noon2"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#combine both DHM dfs
DHM_metrics <- DHM_metrics_noon1%>%
  full_join(DHM_metrics_noon2)

#initialize final migration df
migration_df <- data.frame(metric = DHM_metrics$metric, MSN = DHM_metrics$MSN)

#get DVM df into same order as DHM df
DVM_metrics <- DVM_metrics[order(DHM_metrics$metric, DHM_metrics$MSN),]

#add average and SE of noon DVM and DHM calcs - sunrise for DHM is still separate
migration_df$DVM_avg <- rowMeans(DVM_metrics[,c(3,4)], na.rm=TRUE)
migration_df$DHM_avg <- rowMeans(DHM_metrics[,c(3,4)], na.rm=TRUE)

#calculate stderr across for migration metrics
migration_df$DVM_SE <- apply(DVM_metrics[,c(3,4)], 1, stderr)
migration_df$DHM_SE <- apply(DHM_metrics[,c(3,4)], 1, stderr)

#replace NAN with 0 bc none of that taxa were found
migration_df$DHM_avg[is.nan(migration_df$DHM_avg)] <- 0

#now convert from wide to long
metrics_avg <- migration_df %>% gather(migration, value, DVM_avg:DHM_avg)
metrics_se <- migration_df %>% gather(migration, value, DVM_SE:DHM_SE)

migration_long <- metrics_avg[,c(1,2,5,6)]
migration_long$SE <- metrics_se$value
  
#export migration metrics
#write.csv(migration_long,"./output/migration_metrics.csv",row.names = FALSE)

#-------------------------------------------------------------------------------#

metric_taxa <-c("Total","Calanoida","Calanoida","Cladocera","Cladocera",
                "Copepoda","Copepoda","Cyclopoida", "Cyclopoida",
                 "nauplius","nauplius","Rotifera", "Rotifera", "Total")
names(metric_taxa) <- c(unique(migration_long$metric))

#reorder taxa
migration_long$metric <- factor(migration_long$metric, levels = c(
  "ZoopDensity_No.pL",unique(migration_long$metric)[1:13]))

#change facet labels
metric_taxa <-c("1", "2", "3", "a: cladocerans", "a: cladocerans",
                "b: copepods","b: copepods","8","9","10","11",
                "c: rotifers","c: rotifers","14")
names(metric_taxa) <- c(unique(migration_long$metric))

#plot migration metrics
ggplot(subset(migration_long, grepl("density",metric, ignore.case=T) & 
                metric %in% c("Cladocera_density_NopL","Copepoda_density_NopL","Rotifera_density_NopL")), 
              aes(x=as.factor(MSN), y=value, fill=as.factor(MSN), shape=migration)) + 
  scale_shape_manual("",values = c(24,21), labels = c("DHM","DVM")) + xlab("") +
  scale_x_discrete(breaks=c("1","2","3","4","5"),
                   labels=c("10-11 Jul 2019", "24-25 Jul 2019", "12-13 Aug 2020",
                            "15-16 Jun 2021", "7-8 Jul 2021")) +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE, color=as.factor(MSN)), width=.4,position=position_dodge(.9), linewidth = 0.8) +
  geom_point(position=position_dodge(.9), size=2) + theme_bw() + geom_hline(yintercept = 0, linetype="dotted")+
  scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
  scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0,3,0,0), 'lines'),
        legend.position = c(0.92,0.94), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line")) + guides(fill="none", color='none') +
  facet_wrap(~metric, labeller = labeller(metric=metric_taxa)) + ylab("Density migration metric") +
  geom_text(aes(x=5.9, y=c(rep(0,28),0.4,-0.4), label=c(rep(NA,28),"Typical \nMigration", "Reverse \nMigration")), 
            hjust = 0, size = 3, color="black") + coord_cartesian(xlim = c(1, 5), clip = 'off')
#ggsave("figures/BVR_MSNs_migration_metrics_dens_3taxa.jpg", width=5, height=4) 

ggplot(subset(migration_long, grepl("biomass",metric, ignore.case = TRUE) &
                metric %in% c("Cladocera_BiomassConcentration_ugpL","Copepoda_BiomassConcentration_ugpL","Rotifera_BiomassConcentration_ugpL")), 
       aes(x=MSN, y=value, color=metric, shape=migration)) + 
  geom_point(position=position_dodge(.9)) + theme_bw() + geom_hline(yintercept = 0, linetype="dotted")+
  scale_shape_manual("",values = c(1, 19), labels = c("DHM","DVM")) +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2,position=position_dodge(.9)) +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        axis.text.x = element_text(vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.92,0.94), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line")) + guides(color="none") +
  scale_color_manual("",values=c("#088f84","#ccae8f","#a61737"))+
  facet_wrap(~metric, labeller = labeller(metric=metric_taxa)) + ylab("Biomass migration metric")
#ggsave("figures/BVR_MSNs_migration_metrics_biom_3taxa.jpg", width=5, height=4) 

#-------------------------------------------------------------------------------
#create df with proportion of total zoops (both density and biomass) over time

#now calculate the proportion in each habitat
Hourly_prop <- plyr::ddply(all_DHM, c("metric", "MSN", "Hour","DateTime"), function(x) {
  data.frame(
    proportion_lit = x$value[x$site_no=="BVR_l"] / sum(x$value),
    proportion_pel = x$value[x$site_no=="BVR_50"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#export proportion df
#write.csv(Hourly_prop,"output/Hourly_proportions_pelvslit.csv",row.names = FALSE)

#-------------------------------------------------------------------------------
#plot migration metric vs PAR to see if sunny/cloudy days affect migration metric

#first, read in ctd
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1))

# MEL comment: is this data actually used in manuscript? if not, would remove
# if it IS used, we need to figure out a reliable way of downloading it, at a minimum
# the timeout needs to be extended; here are the errors I got below:
# 
# sh: 1: curl: not found
# Error in download.file(inUrl1, infile1, method = "curl") : 
#   'curl' call had nonzero exit status
# In addition: Warning message:
#   In system(paste("curl", paste(extra, collapse = " "), shQuote(url),  :
#                     error in running command
#                   Error in download.file(inUrl1, infile1) : 
#                     download from 'https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf' failed
#                   In addition: Warning messages:
#                     1: In download.file(inUrl1, infile1) :
#                     downloaded length 0 != reported length 0
#                   2: In download.file(inUrl1, infile1) :
#                     URL 'https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf': Timeout of 60 seconds was reached
# if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "DateTime",     
                 "Depth_m",     
                 "Temp_C",     
                 "DO_mgL",     
                 "DOsat_percent",     
                 "Cond_uScm",     
                 "SpCond_uScm",     
                 "Chla_ugL",     
                 "Turbidity_NTU",     
                 "pH",     
                 "ORP_mV",     
                 "PAR_umolm2s",     
                 "DescRate_ms",     
                 "Flag_DateTime",     
                 "Flag_Temp_C",     
                 "Flag_DO_mgL",     
                 "Flag_DOsat_percent",     
                 "Flag_Cond_uScm",     
                 "Flag_SpCond_uScm",     
                 "Flag_Chla_ugL",     
                 "Flag_Turbidity_NTU",     
                 "Flag_pH",     
                 "Flag_ORP_mV",     
                 "Flag_PAR_umolm2s",     
                 "Flag_DescRate_ms"    ), check.names=TRUE)

unlink(infile1)

#next filter just BVR PAR data from 0.1m
par_bvr = dt1%>%
  mutate(Date = as.Date(DateTime))%>%
  filter(Depth_m >= 0.1 & Depth_m <= 0.11,
         Site == 50,
         Reservoir=="BVR",
         Date %in% c(as.Date("2019-07-10"), as.Date("2019-07-11"),  as.Date("2019-07-24"),
                          as.Date("2019-07-25"), as.Date("2020-08-12"), as.Date("2020-08-13"),
                          as.Date("2021-06-15"), as.Date("2021-06-16"), as.Date("2021-07-07"),
                          as.Date("2021-07-08")))

#average surface PARs by date
avg_par_bvr <- par_bvr %>% group_by(as.Date(DateTime)) %>% 
  summarise(avg_par = mean(PAR_umolm2s))

#add msn col and do a more hacky avg
avg_par_bvr$msn <- c(1,1,2,2,3,4)

avg_par_bvr <- avg_par_bvr %>% group_by(msn) %>%
  summarise(avg_par = mean(avg_par))

#plot migration metric vs par
plot(migration_df$DVM_avg[migration_df$metric=="ZoopDensity_No.pL"] ~
       c(avg_par_bvr$avg_par,NA), pch=19, 
     xlab = "PAR (ummolm2s)", 
     ylab = "DVM metric")

plot(migration_df$DHM_avg[migration_df$metric=="ZoopDensity_No.pL"] ~
       c(avg_par_bvr$avg_par,NA), pch=19,
     xlab = "PAR (ummolm2s)", 
     ylab = "DHM metric")

#hmm potential inverse relationship between migration metric and PAR??

#------------------------------------------------------------------------------#
#read in secchi 
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "DateTime",     
                 "Secchi_m",     
                 "Flag_DateTime",     
                 "Flag_Secchi_m"    ), check.names=TRUE)

unlink(infile1)

#create df for secchi data
secchi_bvr = dt1 %>% 
  mutate(Date = as.Date(DateTime))%>%
  filter(Site == 50,
         Reservoir=="BVR",
         Date %in% c(as.Date("2019-07-10"), as.Date("2019-07-11"),  as.Date("2019-07-24"),
                     as.Date("2019-07-25"), as.Date("2020-08-12"), as.Date("2020-08-13"),
                     as.Date("2021-06-15"), as.Date("2021-06-16"), as.Date("2021-07-07"),
                     as.Date("2021-07-08"))) %>%
  select(Date,Secchi_m)

#add the day 5 secchi (calculated using zoop rope at pelagic site)
secchi_bvr <- data.frame("msn" = c(1,2,3,4,5),
                         "secchi_m" = c(secchi_bvr$Secchi_m, 1.8))

#now plot
plot(migration_df$DVM_avg[migration_df$metric=="Rotifera_density_NopL"] ~
       secchi_bvr$secchi_m, pch=19, 
     xlab = "Secchi (m)", 
     ylab = "DVM metric")

plot(migration_df$DHM_avg[migration_df$metric=="Rotifera_density_NopL"] ~
       secchi_bvr$secchi_m, pch=19,
     xlab = "Secchi (m)", 
     ylab = "DHM metric")
#cool! positive relationship between migration metric and secchi 
#normal migration with high water clarity and reverse for lower water clarity
