#Script to calculate DVM and DHM metrics
#Created 5Dec2022

source("scripts/install.R")

#read in DHM csv
all_DHM <- read.csv("output/All_MSN_tows_DHM.csv")

#new col for time
all_DHM$time <- ifelse(substr(all_DHM$Hour,12,13)=="12", "noon",
                       ifelse(substr(all_DHM$Hour,12,13)=="", "midnight",
                              ifelse(substr(all_DHM$Hour,12,13) %in% c(
                                "04", "05", "06", "07"), "sunrise", "sunset")))
  

#rename Hour to DateTime
names(all_DHM)[names(all_DHM)=="Hour"] <- "DateTime"

#add sunrise and sunset hours
all_DHM$time <- ifelse(substr(all_DHM$DateTime,12,13) == "04", "sunrise1",
                ifelse(substr(all_DHM$DateTime,12,13) == "05", "sunrise2",
                ifelse(substr(all_DHM$DateTime,12,13) == "06", "sunrise3",
                ifelse(substr(all_DHM$DateTime,12,13) == "07", "sunrise4",
                ifelse(substr(all_DHM$DateTime,12,13) == "18", "sunset1",
                ifelse(substr(all_DHM$DateTime,12,13) == "19", "sunset2",
                ifelse(substr(all_DHM$DateTime,12,13) == "20", "sunset3",
                ifelse(substr(all_DHM$DateTime,12,13) == "21", "sunset4",
                       all_DHM$time))))))))
  

#add MSN# column
all_DHM$MSN <- ifelse(all_DHM$date=="2019-07-10" | 
                        all_DHM$date=="2019-07-11",1,
               ifelse(all_DHM$date=="2019-07-24" | 
                        all_DHM$date=="2019-07-25",2,
               ifelse(all_DHM$date=="2020-08-12" | 
                        all_DHM$date=="2020-08-13",3,
               ifelse(all_DHM$date=="2021-06-15" | 
                        all_DHM$date=="2021-06-16",4,5))))

#change second noon for each MSN to noon2
all_DHM$time <- ifelse(all_DHM$time=="noon" & 
                         (all_DHM$date=="2019-07-11" | 
                         all_DHM$date=="2019-07-25" |
                         all_DHM$date=="2020-08-13" |
                         all_DHM$date=="2021-06-16" |
                         all_DHM$date=="2021-07-08"), "noon2", all_DHM$time)


#split up dfs for noon1 vs noon2 data
all_DHM_noon <- all_DHM[all_DHM$time=="noon" | all_DHM$time=="midnight",]
all_DHM_noon2 <- all_DHM[all_DHM$time=="noon2" | all_DHM$time=="midnight",]

#-------------------------------------------------------------------------------
#read in DVM annual csvs
DVM_2021 <- read.csv("output/DVM_2021_zoops.csv") 
DVM_2020 <- read.csv("output/DVM_2020_zoops.csv") 
DVM_2019 <- read.csv("output/DVM_2019_zoops.csv") 

#combine annual DVM dfs
all_DVM <- rbind(DVM_2019,DVM_2020,DVM_2021)

#change date format
all_DVM$DateTime <- format(as.Date(all_DVM$DateTime, "%m-%d-%Y"), "%Y-%m-%d")

#work with raw density and biomass to calculate metrics
all_DVM <- all_DVM[!grepl("percent",all_DVM$metric),]

#add MSN# column
all_DVM$MSN <- ifelse(all_DVM$DateTime=="2019-07-10" | 
                        all_DVM$DateTime=="2019-07-11",1,
               ifelse(all_DVM$DateTime=="2019-07-24" | 
                        all_DVM$DateTime=="2019-07-25",2,
               ifelse(all_DVM$DateTime=="2020-08-12" | 
                        all_DVM$DateTime=="2020-08-13",3,
               ifelse(all_DVM$DateTime=="2021-06-15" | 
                        all_DVM$DateTime=="2021-06-16",4,5))))

#change second noon for each MSN to noon2
all_DVM$Hour <- ifelse(all_DVM$Hour=="noon" & 
                         (all_DVM$DateTime=="2019-07-11" | 
                          all_DVM$DateTime=="2019-07-25" |
                          all_DVM$DateTime=="2020-08-13" |
                          all_DVM$DateTime=="2021-06-16" |
                          all_DVM$DateTime=="2021-07-08"), "noon2", all_DVM$Hour)

#rename Hour to time
names(all_DVM)[names(all_DVM)=="Hour"] <- "time"

#not sure if this is the right call, but going to set all negative values to 0 bc they are not real and due to additive error across all steps
all_DVM$value[all_DVM$value < 0] <- 0 #n=6

#shorten metric name
all_DVM$metric <- ifelse(all_DVM$WaterColumn=="epilimnion", substr(all_DVM$metric,1,nchar(all_DVM$metric)-13),substr(all_DVM$metric,1,nchar(all_DVM$metric)-14))

#split up dfs for noon1 vs noon2 data
all_DVM_noon1 <- all_DVM[all_DVM$time!="noon2",]
all_DVM_noon2 <- all_DVM[all_DVM$time!="noon",]

#-------------------------------------------------------------------------------
####                         DVM METRICS                                    ####
#-------------------------------------------------------------------------------
#Calculate proportion of zoops in epi at noon and midnight
DVM_proportion_noon1 <-  plyr::ddply(all_DVM_noon1, c("metric", "MSN", "time","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon1 = x$value[x$WaterColumn=="epilimnion"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

DVM_proportion_noon2 <-  plyr::ddply(all_DVM_noon2, c("metric", "MSN", "time","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon2 = x$value[x$WaterColumn=="epilimnion"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#initialize df for DVM metrics for each day
#migration_df <- data.frame("MSN"= unique(DVM_proportion$MSN))

#DVM metrics for a single day --> DVM = (Depi / Depi + Dhypo)Night - (Depi / Depi + Dhypo)Day
DVM_noon1 <-  plyr::ddply(DVM_proportion_noon1, c("metric", "MSN"), function(x) {
  data.frame(
    DVM_metric_noon1 = x$proportion_epi[x$time=="midnight"] - 
      x$proportion_epi[x$time=="noon"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


DVM_noon2 <-  plyr::ddply(DVM_proportion_noon2, c("metric", "MSN"), function(x) {
  data.frame(
DVM_metric_noon2 = x$proportion_epi[x$time=="midnight"] - 
  x$proportion_epi[x$time=="noon2"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#merge both dfs 
DVM_metrics <- right_join(DVM_noon1,DVM_noon2, by=c("metric","MSN"))

#-------------------------------------------------------------------------------
####                         DHM METRICS                                    ####
#-------------------------------------------------------------------------------
#Calculate proportion of zoops in epi at noon and midnight (or sunrise/sunset)
DHM_proportion_noon <-  plyr::ddply(all_DHM_noon, c("Taxon", "MSN", "time","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon_dens = x$Density_IndPerL_rep.mean[x$Site=="50"] / 
      sum(x$Density_IndPerL_rep.mean),
    proportion_epi_noon_biom = x$Biomass_ugL_rep.mean[x$Site=="50"] / 
      sum(x$Biomass_ugL_rep.mean)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

DHM_proportion_noon2 <-  plyr::ddply(all_DHM_noon2, c("Taxon", "MSN", "time","DateTime"), function(x) {
  data.frame(
    proportion_epi_noon2_dens = x$Density_IndPerL_rep.mean[x$Site=="50"] / 
      sum(x$Density_IndPerL_rep.mean),
    proportion_epi_noon2_biom = x$Biomass_ugL_rep.mean[x$Site=="50"] / 
      sum(x$Biomass_ugL_rep.mean)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#DHM metrics for a single day --> DHM = (Dpelepi / Dpelepi + Dlit)Night - (Dpelepi / Dpelepi + Dlit)Day
DHM_metrics_noon1 <-  plyr::ddply(DHM_proportion_noon, c("Taxon", "MSN"), function(x) {
  data.frame(
    DHM_metric_noon1_dens = x$proportion_epi_noon_dens[x$time=="midnight"] - 
      x$proportion_epi_noon_dens[x$time=="noon"],
    DHM_metric_noon1_biom = x$proportion_epi_noon_biom[x$time=="midnight"] - 
      x$proportion_epi_noon_biom[x$time=="noon"]
    )
}, .progress = plyr::progress_text(), .parallel = FALSE)

DHM_metrics_noon2 <-  plyr::ddply(DHM_proportion_noon2, c("Taxon", "MSN"), function(x) {
  data.frame(
    DHM_metric_noon2_dens = x$proportion_epi_noon2_dens[x$time=="midnight"] - 
      x$proportion_epi_noon2_dens[x$time=="noon2"],
    DHM_metric_noon2_biom = x$proportion_epi_noon2_biom[x$time=="midnight"] - 
      x$proportion_epi_noon2_biom[x$time=="noon2"]
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#combine both DHM dfs
DHM_metrics <- DHM_metrics_noon1 |> 
  full_join(DHM_metrics_noon2)

#only select cladocerans, copepods, and rotifers
DHM_metrics <- DHM_metrics |> 
  filter(Taxon %in% c("Cladocera", "Copepoda", "Rotifera"))

DVM_metrics <- DVM_metrics |> 
  filter(grepl("Cladocera",DVM_metrics$metric) |
           grepl("Copepoda",DVM_metrics$metric) |
           grepl("Rotifera", DVM_metrics$metric))

#average noons for dvm and dhm dfs
DHM_metrics$density_NopL <- rowMeans(DHM_metrics[,c(3,5)], na.rm=TRUE)
DHM_metrics$BiomassConcentration_ugpL <- rowMeans(DHM_metrics[,c(4,6)], na.rm=TRUE)

DVM_metrics$dvm_avg <- rowMeans(DVM_metrics[,c(3,4)], na.rm=TRUE)

#calculate stderr across for migration metrics
DVM_metrics$DVM_SE <- apply(DVM_metrics[,c(3,4)], 1, stderr)

DHM_metrics$DHM_SE_dens <- apply(DHM_metrics[,c(3,5)], 1, stderr)
DHM_metrics$DHM_SE_biom <- apply(DHM_metrics[,c(4,6)], 1, stderr)

#now convert dhm df from wide to long
DHM_metrics_final <- DHM_metrics |> group_by(Taxon, MSN) |> 
  pivot_longer(cols=density_NopL:BiomassConcentration_ugpL, names_to="group")

DHM_metrics_se_final <- DHM_metrics |> group_by(Taxon, MSN) |> 
  pivot_longer(cols=DHM_SE_dens:DHM_SE_biom, names_to="group")

#add metric col (combined taxon + group)
DHM_metrics_final$metric <- paste0(DHM_metrics_final$Taxon, "_", 
                                   DHM_metrics_final$group)

DHM_metrics_se_final$metric <- paste0(DHM_metrics_final$Taxon, "_", 
                                      DHM_metrics_final$group)

#drop unnecessary cols
DHM_metrics_final <- DHM_metrics_final |> select(metric, MSN, value)

#get DVM df into same order as DHM df
DVM_metrics <- DVM_metrics |> arrange(metric,MSN)
DHM_metrics_final <- DHM_metrics_final |> arrange(metric, MSN)
DHM_metrics_se_final <- DHM_metrics_se_final |> arrange(metric, MSN)

#initialize final migration df
migration_df <- data.frame(metric = DVM_metrics$metric, MSN = DHM_metrics$MSN)

#add average for noon DVM and DHM calcs 
migration_df$DVM_avg <- DVM_metrics$dvm_avg
migration_df$DHM_avg <- DHM_metrics_final$value

#add SE
migration_df$DVM_SE <- DVM_metrics$DVM_SE
migration_df$DHM_SE <- DHM_metrics_se_final$value

#replace NAN with 0 bc none of that taxa were found
migration_df$DHM_avg[is.nan(migration_df$DHM_avg)] <- 0

#now convert from wide to long
metrics_avg <- migration_df %>% gather(migration, value, DVM_avg:DHM_avg)
metrics_se <- migration_df %>% gather(migration, value, DVM_SE:DHM_SE)

migration_long <- metrics_avg[,c(1,2,5,6)]
migration_long$SE <- metrics_se$value
  
#export migration metrics
#write.csv(migration_long,"output/migration_metrics.csv",row.names = FALSE)

#-------------------------------------------------------------------------------#
#change facet labels
metric_taxa <-c("cladocera","cladocera","copepoda",
                "copepoda","rotifera","rotifera")
names(metric_taxa) <- c(unique(migration_long$metric))

#plot migration metrics --> Figure 9
ggplot(subset(migration_long, grepl("density",metric, ignore.case=T)), 
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

#Figure S10
ggplot(subset(migration_long, grepl("biomass",metric, ignore.case = TRUE)), 
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
  facet_wrap(~metric, labeller = labeller(metric=metric_taxa)) + ylab("Biomass migration metric") +
geom_text(aes(x=5.9, y=c(rep(0,28),0.4,-0.4), label=c(rep(NA,28),"Typical \nMigration", "Reverse \nMigration")), 
          hjust = 0, size = 3, color="black") + coord_cartesian(xlim = c(1, 5), clip = 'off')
#ggsave("figures/BVR_MSNs_migration_metrics_biom_3taxa.jpg", width=5, height=4) 

#-------------------------------------------------------------------------------
#create df with proportion of total zoops (both density and biomass) over time

library(plyr)

#now calculate the proportion in each habitat
Hourly_prop <- plyr::ddply(all_DHM, c("Taxon", "MSN", "time","DateTime"), function(x) {
  data.frame(
    proportion_lit = x$value[x$site_no=="BVR_l"] / sum(x$value),
    proportion_pel = x$value[x$site_no=="BVR_50"] / sum(x$value)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#export proportion df
#write.csv(Hourly_prop,"output/Hourly_proportions_pelvslit.csv",row.names = FALSE)

detach('package:plyr')
