#Zoop density change over 24hrs for all 5 MSNs
#created 16 Oct 2022

source("scripts/install.R")

#get date into posixct format
zoop_summary$DateTime <- as.POSIXct(zoop_summary$DateTime, 
                                    format="%Y-%m-%d %H:%M:%S", tz="UTC")

#create column with just date for subsetting
zoop_summary$date <- as.Date(zoop_summary$DateTime)

#only select 2019-2022 MSN data
zoop_summary <- zoop_summary |> 
  filter(date %in% c("2019-07-10", "2019-07-11", "2019-07-24",
                     "2019-07-25", "2020-08-12", "2020-08-13",
                     "2021-06-15", "2021-06-16", "2021-07-07",
                     "2021-07-08"))

#select cols to keep (dens, biom, and size)
zoop_summary <- zoop_summary |>  
  select("Reservoir","Site","DateTime","StartDepth_m","EndDepth_m", "Rep",
         "CollectionMethod","Taxon", "Density_IndPerL", "Biomass_ugL",
         "MeanLength_mm")

#drop schindler data and dam site (49)
zoop_summary <- zoop_summary |> filter(CollectionMethod!="Schindler",
                                       Site!=49)

#add date column
zoop_summary$date <- as.Date(zoop_summary$DateTime)

#add hour column
zoop_summary$Hour <- substr(zoop_summary$DateTime,12,13)

#change a couple hours so they group properly
zoop_summary$Hour[zoop_summary$Hour=="11"] <- "12"
zoop_summary$Hour[zoop_summary$Hour=="23"] <- "0"
zoop_summary$Hour[zoop_summary$Hour=="03"] <- "04"
zoop_summary$Hour[as.character(zoop_summary$DateTime)=="2019-07-10 19:04:00"] <- "18"

#change 00 to 0 for grouping purposes
zoop_summary$Hour[zoop_summary$Hour=="00"] <- "0"

#order by site, then hour
zoop_summary<- arrange(zoop_summary,zoop_summary$Site,zoop_summary$date)

##### Create new df to combine reps over 24 hours
zoop.repmeans <- zoop_summary %>% select(Site, date, Hour, StartDepth_m,
                                         EndDepth_m, Taxon, Density_IndPerL, 
                                         Biomass_ugL, MeanLength_mm) |> 
  group_by(Site, date, Hour, Taxon, StartDepth_m) %>%
  summarise_at(vars(Density_IndPerL:MeanLength_mm), 
               list(rep.mean=mean, rep.SE=stderr))

#standardize datetime for all MSNs so we can visualize on a single plot
zoop.repmeans$date_temp <- ifelse(zoop.repmeans$date=="2019-07-10" | 
                                zoop.repmeans$date=="2019-07-24"| 
                                zoop.repmeans$date=="2020-08-12" | 
                                zoop.repmeans$date=="2021-06-15" |
                                zoop.repmeans$date=="2021-07-07",
                              "2022-10-15","2022-10-16")

#only select hour and then add arbitrary dates for plotting
zoop.repmeans$Hour <- format(paste0(zoop.repmeans$Hour,":00"), format='%H:%M')

#combine hour and date 
zoop.repmeans$Hour <- strptime(paste0(as.character(zoop.repmeans$date_temp), 
                                      zoop.repmeans$Hour),format="%Y-%m-%d %H:%M")
zoop.repmeans$Hour <- as.POSIXct(zoop.repmeans$Hour)

#make sure zoop.repmeans is a dataframe
zoop.repmeans <- data.frame(zoop.repmeans)

#order by hour for plotting
zoop.repmeans <- zoop.repmeans[order(zoop.repmeans$Hour),]

#Export all zoop data df
write.csv(zoop.repmeans,"output/All_MSN_zoop_tows.csv",row.names = FALSE)

#only select epi samples for DHM plots
zoop_epi <- zoop.repmeans |> filter(StartDepth_m <=4)

#convert new dfs from tibble to dataframe 
zoop_DHM <- data.frame(zoop_epi)

#Export DHM csv
write.csv(zoop_DHM,"output/All_MSN_tows_DHM.csv",row.names = FALSE)

#add column for MSN #
zoop_DHM$MSN <- ifelse(zoop_DHM$date=="2019-07-10" | 
                         zoop_DHM$date=="2019-07-11",1,
                     ifelse(zoop_DHM$date=="2019-07-24" | 
                              zoop_DHM$date=="2019-07-25",2,
                     ifelse(zoop_DHM$date=="2020-08-12" | 
                              zoop_DHM$date=="2020-08-13",3,
                     ifelse(zoop_DHM$date=="2021-06-15" | 
                              zoop_DHM$date=="2021-06-16",4,5))))

sites <- c("Pelagic","Littoral")
names(sites) <- c(50, 51)

#------------------------------------------------------------------------------#
#Figure for zoop density for each MSN 24-hours --> Figure S4
ggplot(subset(zoop_DHM, Taxon %in% c("Cladocera",
              "Copepoda","Rotifera")),
                aes(Hour,Density_IndPerL_rep.mean, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),
                xmax=as.POSIXct("2022-10-15 20:41:00"), ymin=-Inf, ymax= Inf, 
                fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),
                xmax=as.POSIXct("2022-10-16 06:10:00"), 
                ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),
                xmax=as.POSIXct("2022-10-16 12:30:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(Site~Taxon,scales="free_y",
                labeller = labeller(Site=sites)) + 
  xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.1,0.92), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line"))+ 
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), 
                     labels=c("10-11 Jul 2019","24-25 Jul 2019",
                              "12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), 
                     guide=guide_legend(order=1)) + 
  geom_line()+ ylab("Density (Individuals/L)") + 
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+
  geom_errorbar(aes(ymin=Density_IndPerL_rep.mean - Density_IndPerL_rep.SE, 
                    ymax=Density_IndPerL_rep.mean + Density_IndPerL_rep.SE), 
                width=.2,position=position_dodge(.9))
#ggsave("figures/BVR_MSNs_taxa_density.jpg", width=5, height=4) 

#-------------------------------------------------------------------------------#
#new df to standardize density among taxa for ALL days
zoop_dens_stand <- data.frame(subset(zoop_DHM, 
                          Taxon %in% c("Cladocera", "Copepoda", "Rotifera")))

#calculate % density of max within a day for each taxa as x / max
zoop_dens_stand <- zoop_dens_stand %>% group_by(Taxon) %>%
  mutate(value_max_std = Density_IndPerL_rep.mean / 
           max(Density_IndPerL_rep.mean))

#Manuscript Figure 4
ggplot(zoop_dens_stand, aes(Hour,value_max_std, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),
                xmax=as.POSIXct("2022-10-15 20:41:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),
                xmax=as.POSIXct("2022-10-16 06:10:00"), 
                ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),
                xmax=as.POSIXct("2022-10-16 12:30:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(Site~Taxon,scales="free_y",
                      labeller = labeller(Site=sites)) +
  xlab("Hour")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.1,0.44), 
        legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line"))+ 
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT"))+
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), 
                     labels=c("10-11 Jul 2019","24-25 Jul 2019",
                              "12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), 
                     guide=guide_legend(order=1)) + 
  geom_line()+ ylab("Standardized density") + 
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")
#ggsave("figures/BVR_MSNs_taxa_percent_density_over_max_std.jpg", width=5, height=4) 

#numbers in results text
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$Taxon=="Cladocera"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$Taxon=="Copepoda"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$Taxon=="Rotifera"])

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==1])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==2])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==3])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==4])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==5])

mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Cladocera" & 
    zoop_dens_stand$Site=="50"])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Copepoda" & 
    zoop_dens_stand$Site=="50"])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Rotifera" & 
    zoop_dens_stand$Site=="50"])

mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Cladocera" & 
    zoop_dens_stand$Site=="51"])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Copepoda" & 
    zoop_dens_stand$Site=="51"])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Rotifera" & 
    zoop_dens_stand$Site=="51"])

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$Site=="51"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$Site=="50"])

#day vs night differences
#pull out hour + make separate column (20:00 to 6:00 is night)
zoop_dens_stand$hr <- hour(zoop_dens_stand$Hour)

#list of night hours
night <- c(21,22,23,0,1,2,3,4,5)

mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Cladocera" & 
    zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Copepoda" & 
    zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Rotifera" & 
    zoop_dens_stand$hr %in% night])

mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Cladocera" & 
    !zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Copepoda" & 
    !zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[
  zoop_dens_stand$Taxon=="Rotifera" & 
    !zoop_dens_stand$hr %in% night])

#-------------------------------------------------------------------------------------#
# plot df for avg size

#replace 0 sizes with NA (bc we do not have the data and zoops are most definitely not 0 mm long)
zoop_DHM$MeanLength_mm_rep.mean[zoop_DHM$MeanLength_mm_rep.mean==0] <- NA

#Figure S4
ggplot(subset(zoop_DHM, Taxon %in% c("Cladocera", "Copepoda", "Rotifera")),
       aes(Hour,MeanLength_mm_rep.mean, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),
                xmax=as.POSIXct("2022-10-15 20:41:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),
                xmax=as.POSIXct("2022-10-16 06:10:00"), 
                ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),
                xmax=as.POSIXct("2022-10-16 12:30:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + 
  facet_grid(Site~Taxon,scales="free_y",labeller = labeller(Site=sites)) + 
  xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.9,0.95), 
        legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line")) + 
  scale_x_datetime(expand = c(0,0),
                   labels = date_format("%H-%M",tz="EST5EDT")) +
  scale_color_manual("",
                     values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), 
                     labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020",
                              "15-16 Jun 2021","7-8 Jul 2021"), 
                     guide=guide_legend(order=1)) + 
  geom_line()+ ylab("Size (mm)") + 
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+
  geom_errorbar(aes(ymin=MeanLength_mm_rep.mean-MeanLength_mm_rep.SE, 
                    ymax=MeanLength_mm_rep.mean+MeanLength_mm_rep.SE), 
                width=.2,position=position_dodge(.9))
#ggsave("figures/BVR_MSNs_taxa_size.jpg", width=5, height=4) 

