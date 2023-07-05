#Zoop density change over 24hrs for all 5 MSNs
#created 16 Oct 2022

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, 
               ggplot2,tidyr,lubridate, scales, colorblindcheck, viridis)

#create function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#read in zoop data from all 3 years
zoops2019<- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2019.csv'),header = TRUE)
zoops2020<- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2020.csv'),header = TRUE)
zoops2021<- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/FCR_ZooplanktonSummary2021.csv'),header = TRUE)

#drop holopedium becuase only in 2019
zoops2019 <- zoops2019[,!c(grepl("Holopedium",colnames(zoops2019)))]

#combine all zoop datasets
zoop <- rbind(zoops2019,zoops2020,zoops2021)

#make sure sample_ID is class character
zoop$sample_ID<- as.character(zoop$sample_ID)

#change BVR_50_p --> BVR_50
zoop$site_no[zoop$site_no=="BVR_50_p"] <- "BVR_50"

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

#drop 20um samples
zoop <- zoop[!c(substrEnd(zoop$sample_ID,4)=="filt" | substrEnd(zoop$sample_ID,2)=="20"),]

#get hour into character format for grouping
zoop$Hour <- format(round(strptime(paste0(zoop$collect_date, zoop$Hour), format="%Y-%m-%d %H:%M"),units="hours"),format="%H:%M")
#manually change hour of some samples (rounding problems)
zoop$Hour[grepl("midnight",zoop$sample_ID,ignore.case = TRUE)] <- "00:00"
zoop$Hour[grepl("noon",zoop$sample_ID,ignore.case = TRUE)] <- "12:00"
zoop$Hour[grepl("sunrise_h1",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h1",zoop$sample_ID,ignore.case = TRUE)] <- "04:00"
zoop$Hour[grepl("sunrise_h2",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h2",zoop$sample_ID,ignore.case = TRUE)] <- "05:00"
zoop$Hour[grepl("sunrise_h3",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h3",zoop$sample_ID,ignore.case = TRUE)] <- "06:00"
zoop$Hour[grepl("sunrise_h4",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunrise_epi_h4",zoop$sample_ID,ignore.case = TRUE)] <- "07:00"
zoop$Hour[grepl("sunset_h1",zoop$sample_ID,ignore.case = TRUE) | 
            grepl("sunset_epi_h1",zoop$sample_ID,ignore.case = TRUE)] <- "18:00"
zoop$Hour[grepl("sunset_h2",zoop$sample_ID,ignore.case = TRUE) | 
            grepl("sunset_epi_h2",zoop$sample_ID,ignore.case = TRUE)] <- "19:00"
zoop$Hour[grepl("sunset_h3",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunset_epi_h3",zoop$sample_ID,ignore.case = TRUE)] <- "20:00"
zoop$Hour[grepl("sunset_h4",zoop$sample_ID,ignore.case = TRUE) |
            grepl("sunset_epi_h4",zoop$sample_ID,ignore.case = TRUE)] <- "21:00"

#drop the 06-29 samples (these were test samples)
zoop <- zoop[substrEnd(zoop$sample_ID,4)!="filt",]

#drop schindler samples and only select BVR_50 epi samples
zoop <- zoop[substrEnd(zoop$site_no,6)!="schind" & zoop$site_no!="BVR_d" & zoop$site_no!="BVR_dam" & zoop$site_no!="BVR_trap" & zoop$site_no!="FCR_50",]

##### Create new df to combine reps over 24 hours
zoop.repmeans <- zoop %>% select(sample_ID,site_no,collect_date,Hour, ZoopDensity_No.pL, BiomassConcentration_ugpL,
                                 Cladocera_density_NopL, Cladocera_BiomassConcentration_ugpL, Cladocera_PercentOfTotal,
                                 Cyclopoida_density_NopL, Cyclopoida_BiomassConcentration_ugpL, Cyclopoida_PercentOfTotal,
                                 Rotifera_density_NopL, Rotifera_BiomassConcentration_ugpL,Rotifera_PercentOfTotal,
                                 Calanoida_PercentOfTotal, Calanoida_density_NopL, Calanoida_BiomassConcentration_ugpL,
                                 Copepoda_PercentOfTotal, Copepoda_density_NopL, Copepoda_BiomassConcentration_ugpL,
                                 nauplius_PercentOfTotal, nauplius_density_NopL, nauplius_BiomassConcentration_ugpL) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(ZoopDensity_No.pL:nauplius_BiomassConcentration_ugpL), list(rep.mean=mean, rep.SE=stderr))

#merge collect_date and hour in a new column
zoop.repmeans$datetime<- paste(zoop.repmeans$collect_date,zoop.repmeans$Hour,sep=" ")
#get times into date format (character here)
zoop.repmeans$datetime<- format(as.POSIXct(zoop.repmeans$datetime,format="%Y-%m-%d %H:%M"), format="%Y-%m-%d %H:%M:%S")
#convert to posixct date format
zoop.repmeans$datetime<- as.POSIXct(zoop.repmeans$datetime, format="%Y-%m-%d %H:%M")

zoop.repmeans$dates <- ifelse(zoop.repmeans$collect_date=="2019-07-10" | zoop.repmeans$collect_date=="2019-07-24"| 
                                zoop.repmeans$collect_date=="2020-08-12" | zoop.repmeans$collect_date=="2021-06-15" |
                                zoop.repmeans$collect_date=="2021-07-07","2022-10-15","2022-10-16")

#only select hour and then add arbitrary dates for plotting
zoop.repmeans$Hour <- format(zoop.repmeans$Hour, format='%H:%M')

#combine hour and date
zoop.repmeans$Hour <- strptime(paste0(as.character(zoop.repmeans$dates), zoop.repmeans$Hour),format="%Y-%m-%d %H:%M")
zoop.repmeans$Hour <- as.POSIXct(zoop.repmeans$Hour)

#make sure zoop.repmeans is a dataframe
zoop.repmeans <- data.frame(zoop.repmeans)

#order by hour for plotting
zoop.repmeans <- zoop.repmeans[order(zoop.repmeans$Hour),]

#Export all zoop data df
#write.csv(zoop.repmeans,"./Summer2021-DataAnalysis/SummaryStats/All_MSN_zoops.csv",row.names = FALSE)

#only select epi samples for DHM plots
zoop_epi <- zoop.repmeans[grepl("epi",zoop.repmeans$sample_ID) |grepl("sunrise",zoop.repmeans$sample_ID) | grepl("sunset",zoop.repmeans$sample_ID) | zoop.repmeans$site_no=="BVR_l",]

#convert new dfs from tibble to dataframe 
zoop_DHM <- data.frame(zoop_epi)

#convert df from wide to long (kinda hacky way bc having problems doing this)
df1 <- zoop_DHM %>% gather(metric,value,ZoopDensity_No.pL_rep.mean:nauplius_BiomassConcentration_ugpL_rep.mean)
df2 <- zoop_DHM %>% gather(metric.SE,value.SE, ZoopDensity_No.pL_rep.SE:nauplius_BiomassConcentration_ugpL_rep.SE)

##cut and paste to merge df
zoop_DHM_long <- df1[,c(1:4,27,28)]
#zoop_DHM_long$metric.SE <- df2$metric.SE #use this as a check to make sure rows match up
zoop_DHM_long$value.SE <- df2$value.SE

#Export DHM csv
#write.csv(zoop_DHM_long,"./Summer2021-DataAnalysis/SummaryStats/All_MSN_DHM.csv",row.names = FALSE)

#reset DHM df so can only look at density
zoop_DHM_long <- NA

variables <- c("ZoopDensity_No.pL_rep.mean","Cladocera_density_NopL_rep.mean", "Cladocera_PercentOfTotal_rep.mean",
               "Cyclopoida_density_NopL_rep.mean","Cyclopoida_PercentOfTotal_rep.mean","Rotifera_density_NopL_rep.mean",
               "Rotifera_PercentOfTotal_rep.mean","Calanoida_density_NopL_rep.mean","Calanoida_PercentOfTotal_rep.mean",
               "Copepoda_density_NopL_rep.mean","Copepoda_PercentOfTotal_rep.mean", 
               "nauplius_density_NopL_rep.mean","nauplius_PercentOfTotal_rep.mean")
SE <- c("ZoopDensity_No.pL_rep.SE","Cladocera_density_NopL_rep.SE", "Cladocera_PercentOfTotal_rep.SE",
        "Cyclopoida_density_NopL_rep.SE","Cyclopoida_PercentOfTotal_rep.SE","Rotifera_density_NopL_rep.SE",
        "Rotifera_PercentOfTotal_rep.SE","Calanoida_density_NopL_rep.SE","Calanoida_PercentOfTotal_rep.SE",
        "Copepoda_density_NopL_rep.SE","Copepoda_PercentOfTotal_rep.SE",
        "nauplius_density_NopL_rep.SE","nauplius_PercentOfTotal_rep.SE")

#remove biomass too
zoop_DHM <- zoop_DHM[,-c(which(grepl("ugpL",colnames(zoop_DHM))))]

#convert df from wide to long (kinda hacky way bc having problems doing this)
df1 <- zoop_DHM %>% gather(metric,value,all_of(variables))
df2 <- zoop_DHM %>% gather(metric.SE,value.SE, all_of(SE))

#cut and paste to merge df
zoop_DHM_long <- df1[,c(1:4,18,20:21)]
#zoop_DHM_long$metric.SE <- df2$metric.SE #use this as a check to make sure rows match up
zoop_DHM_long$value.SE <- df2$value.SE

#drop _rep.mean from all metric names
zoop_DHM_long$metric <- substr(zoop_DHM_long$metric,1,nchar(zoop_DHM_long$metric)-9)

#add column for MSN #
zoop_DHM_long$MSN <- ifelse(zoop_DHM_long$collect_date=="2019-07-10" | zoop_DHM_long$collect_date=="2019-07-11",1,
                            ifelse(zoop_DHM_long$collect_date=="2019-07-24" | zoop_DHM_long$collect_date=="2019-07-25",2,
                                   ifelse(zoop_DHM_long$collect_date=="2020-08-12" | zoop_DHM_long$collect_date=="2020-08-13",3,
                                          ifelse(zoop_DHM_long$collect_date=="2021-06-15" | zoop_DHM_long$collect_date=="2021-06-16",4,5))))

metric_taxa <-c("ZoopDensity","Cladocera","Cladocera","Cyclopoida",
                "Cyclopoida", "Rotifera", "Rotifera", "Calanoida",
                "Calanoida","Copepoda","Copepoda","Nauplius","Nauplius")
names(metric_taxa) <- c(unique(zoop_DHM_long$metric))

sites <- c("Pelagic","Littoral")
names(sites) <- c("BVR_50","BVR_l")

#cb_friendly_2 <- c("#8C510A", "#BF812D","#DFC27D", "#C7EAE5", "#35978F")

#------------------------------------------------------------------------------#
#Figure for zoop density for each MSN 24-hours 
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Copepoda_density_NopL","Rotifera_density_NopL")),
                aes(Hour,value, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),xmax=as.POSIXct("2022-10-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),xmax=as.POSIXct("2022-10-16 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),xmax=as.POSIXct("2022-10-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(metric=metric_taxa, site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.92), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line"))+ scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  #scale_colour_viridis_d("",option="viridis", labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021")) + 
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) + 
  geom_line()+ ylab("Density (Individuals/L)") + scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_MSNs_taxa_density.jpg"), width=5, height=4) 

ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_PercentOfTotal","Copepoda_PercentOfTotal","Rotifera_PercentOfTotal")), #removing calanoids because they are super low %
                aes(Hour,value, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),xmax=as.POSIXct("2022-10-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),xmax=as.POSIXct("2022-10-16 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),xmax=as.POSIXct("2022-10-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",labeller = labeller(metric=metric_taxa, site_no=sites)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.44), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line"))+ scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT"))+
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) + 
  geom_line()+ ylab("% Density") + scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+ 
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_MSNs_taxa_percent_density.jpg"), width=5, height=4) 

range(zoop_DHM_long$value[zoop_DHM_long$metric=="Copepoda_PercentOfTotal"]) # 77.1%, 68.9%, 85.9%

#-------------------------------------------------------------------------------#
#new df to standardize density among taxa and days
zoop_dens_stand <- data.frame(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Copepoda_density_NopL","Rotifera_density_NopL")))

#calculate % density of max within a day for each taxa as x / max
zoop_dens_stand <- zoop_dens_stand %>% group_by(metric, MSN) %>%
  mutate(value_max_std = value / max(value))

#change facet labels
#metric_taxa <-c("cladocerans","rotifers","copepods")
#names(metric_taxa) <- c(unique(zoop_dens_stand$metric))

ggplot(zoop_dens_stand, aes(Hour,value_max_std, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),xmax=as.POSIXct("2022-10-15 20:41:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),xmax=as.POSIXct("2022-10-16 06:10:00"), 
                ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),xmax=as.POSIXct("2022-10-16 12:30:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(site_no~metric,scales="free_y",
                                               labeller = labeller(metric=metric_taxa, site_no=sites)) +
  xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.1,0.44), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line"))+ 
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT"))+
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), 
                     labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) + 
  geom_line()+ ylab("density / max dens") + 
  scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_MSNs_taxa_percent_density_over_max_std.jpg"), width=5, height=4) 


mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Cladocera_density_NopL"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Copepoda_density_NopL"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Rotifera_density_NopL"])

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==1])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==2])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==3])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==4])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$MSN==5])

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Cladocera_density_NopL" & zoop_dens_stand$site_no=="BVR_50"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Copepoda_density_NopL" & zoop_dens_stand$site_no=="BVR_50"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Rotifera_density_NopL" & zoop_dens_stand$site_no=="BVR_50"])

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Cladocera_density_NopL" & zoop_dens_stand$site_no=="BVR_l"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Copepoda_density_NopL" & zoop_dens_stand$site_no=="BVR_l"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Rotifera_density_NopL" & zoop_dens_stand$site_no=="BVR_l"])

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$site_no=="BVR_l"])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$site_no=="BVR_50"])

#day vs night differences
#pull out hour + make separate column (20:00 to 6:00 is night)
zoop_dens_stand$hr <- hour(zoop_dens_stand$Hour)

#list of night hours
night <- c(21,22,23,0,1,2,3,4,5)

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Cladocera_density_NopL" & zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Copepoda_density_NopL" & zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Rotifera_density_NopL" & zoop_dens_stand$hr %in% night])

mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Cladocera_density_NopL" & !zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Copepoda_density_NopL" & !zoop_dens_stand$hr %in% night])
mean(zoop_dens_stand$value_max_std[zoop_dens_stand$metric=="Rotifera_density_NopL" & !zoop_dens_stand$hr %in% night])

#-------------------------------------------------------------------------------------#
#looking at each MSN separately to look for evidence of DHM
ggplot(subset(zoop_DHM_long, metric %in% c("Cladocera_density_NopL","Copepoda_density_NopL","Rotifera_density_NopL") &
                MSN==5 & site_no=="BVR_l"), aes(Hour,value)) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),xmax=as.POSIXct("2022-10-15 20:41:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),xmax=as.POSIXct("2022-10-16 06:10:00"), ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),xmax=as.POSIXct("2022-10-16 12:30:00"), ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + facet_grid(~metric,scales="free_y",labeller = labeller(metric=metric_taxa)) + xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), legend.background = element_blank(), legend.key = element_blank(), legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(fill = "transparent"), legend.position = c(0.1,0.92), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line"))+ scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  geom_line()+ ylab("Density (Individuals/L)") + scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
#take-homes: DHM for rotifers on MSN #1,3, cyclopoids on MSN #2, cladocerans do not
#side note that MSN #5 is weird becuase there are basically no clads/copes and rotifer patterns are weird??

#-------------------------------------------------------------------------------------#
# DVM calcs for epi and hypo density across all 5 campaigns

#read in all 3 DVM tables from 2019-2021
DVM_2019 <- read.csv(file.path(getwd(),'Summer2019-DataAnalysis/SummaryStats/DVM_2019_zoops.csv'),header = TRUE)
DVM_2020 <- read.csv(file.path(getwd(),'Summer2020-DataAnalysis/SummaryStats/DVM_2020_zoops.csv'),header = TRUE) %>% select(-X)
DVM_2021 <- read.csv(file.path(getwd(),'Summer2021-DataAnalysis/SummaryStats/DVM_2021_zoops.csv'),header = TRUE) %>% select(-X)

#merge all zoop files
DVM_all <- rbind(DVM_2019,DVM_2020,DVM_2021)

#get datetime column into right format
DVM_all$date <- format(as.POSIXct(DVM_all$DateTime,format="%m-%d-%Y"), format="%Y-%m-%d %H:%M")
#convert to posixct date format
DVM_all$date<- as.Date(DVM_all$date, format="%Y-%m-%d %H:%M")

#change second noon for each MSN to noon2
DVM_all$Hour <- ifelse(DVM_all$Hour=="noon" & (DVM_all$date=="2019-07-11" | 
                                                 DVM_all$date=="2019-07-25" |
                                                 DVM_all$date=="2020-08-13" |
                                                 DVM_all$date=="2021-06-16" |
                                                 DVM_all$date=="2021-07-08"), "noon2", DVM_all$Hour)


#split up into percents df and dens/biom df
DVM_percent <- DVM_all[grepl("percent", DVM_all$metric),]
  
DVM <- DVM_all[!grepl("percent", DVM_all$metric),]

#drop _rep.mean_epi or _rep_mean_hypo from all metric names
DVM$metric <- ifelse(DVM$WaterColumn=="epilimnion",
                         substr(DVM$metric,1,nchar(DVM$metric)-13), 
                         substr(DVM$metric,1,nchar(DVM$metric)-14))

#same for percent df #FIX sampleid and watercolumn discrepancy
DVM_percent$metric <- ifelse(as.factor(DVM_percent$WaterColumn)=="epilimnion",
                     substr(DVM_percent$metric,1,nchar(DVM_percent$metric)-29), 
                     substr(DVM_percent$metric,1,nchar(DVM_percent$metric)-30))
                         
#set negative values to 0
DVM_percent$value[DVM_percent$value < 0] <- 0


MSN1 <- c("2019-07-10","2019-07-11")
MSN2 <- c("2019-07-24","2019-07-25")
MSN3 <- c("2020-08-12","2020-08-13")
MSN4 <- c("2021-06-15","2021-06-16")
MSN5 <- c("2021-07-07","2021-07-08")


ggplot(subset(DVM, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Calanoida_density_NopL","Rotifera_density_NopL") &
                date %in% as.Date(MSN3)), aes(x=WaterColumn, y=value, fill=Taxa)) +
  geom_bar(aes(fill=Taxa, group=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_rect(data=subset(DVM, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Calanoida_density_NopL","Rotifera_density_NopL") &
                     date %in% as.Date(MSN3) & Hour %in% c("midnight")), 
            aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  facet_grid(~Taxa+Hour, scales= 'free') + #labeller=as_labeller(function(variable,value){c("","","","","Noon")})
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),16)) + labs(title=MSN1[1]) + 
  scale_fill_manual(values=viridis(6),labels=c("Calanoida", "Cladocera", "Cyclopoida","Rotifera")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=7,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_epimetahypo_density_MSN1.jpg"), width=5, height=4) 


ggplot(subset(DVM_percent, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Calanoida_density_NopL","Rotifera_density_NopL") &
                date %in% as.Date(MSN5)), aes(x=WaterColumn, y=value, fill=Taxa)) +
  geom_bar(aes(fill=Taxa, group=WaterColumn), stat="identity", position=position_dodge(),show.legend = TRUE) + theme_bw() +
  geom_rect(data=subset(DVM, metric %in% c("Cladocera_density_NopL","Cyclopoida_density_NopL","Calanoida_density_NopL","Rotifera_density_NopL") &
                          date %in% as.Date(MSN5) & Hour %in% c("midnight")), 
            aes(fill=Hour),xmin=-Inf ,xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'black', alpha = 0.053,inherit.aes = FALSE) +
  facet_grid(~Taxa+Hour, scales= 'free') + #labeller=as_labeller(function(variable,value){c("","","","","Noon")})
  theme(strip.text.y = element_text(size = 10,margin = margin(0, -0.01, 0,0.1, "cm")),strip.background = element_blank(), plot.margin = margin(2.8,0,2,0,unit = "cm")) + 
  scale_x_discrete(name="", labels=rep(c("epi","hypo"),16)) + labs(title=MSN5[1]) + 
  scale_fill_manual(values=viridis(6),labels=c("Calanoida", "Cladocera", "Cyclopoida","Rotifera")) +
  geom_hline(yintercept=0) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                                  axis.text.x=element_text(size=7,family="Times"), plot.title = element_text(hjust = 0.5)) + ylab("Density (individual/L)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 1, 2, 1), legend.title = element_text(size=10),legend.key.size = unit(0.3,"cm"),
        legend.box="vertical", legend.box.spacing = unit(0.1,"cm"),legend.text = element_text(size=8),axis.text.y = element_text(size=10, family="Times"))
ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_epimetahypo_percent_density_MSN5.jpg"), width=5, height=4) 
# DVM for MSN 1 - cyclopoids and calanoids, MSN 2 - cladocera and cyclopoids, 
      #MSN 3 - cladocerans and cyclopoids (reverse DVM) and rotifers,
      #BUT, nothing for MSN 4/5 

#-------------------------------------------------------------------------------------#
# create new df for avg size
zoop_size <- zoop %>% select(sample_ID,site_no,collect_date,Hour, which(grepl("MeanSize",colnames(zoop)))) %>%
  group_by(sample_ID, site_no, Hour, collect_date) %>%
  summarise_at(vars(OverallMeanSize_mm:LecaneMeanSize_mm), list(rep.mean=mean, rep.SE=stderr))

#merge collect_date and hour in a new column
zoop_size$datetime<- paste(zoop_size$collect_date,zoop_size$Hour,sep=" ")
#get times into date format (character here)
zoop_size$datetime<- format(as.POSIXct(zoop_size$datetime,format="%Y-%m-%d %H:%M"), format="%Y-%m-%d %H:%M:%S")
#convert to posixct date format
zoop_size$datetime<- as.POSIXct(zoop_size$datetime, format="%Y-%m-%d %H:%M")

zoop_size$dates <- ifelse(zoop_size$collect_date=="2019-07-10" | zoop_size$collect_date=="2019-07-24"| 
                            zoop_size$collect_date=="2020-08-12" | zoop_size$collect_date=="2021-06-15" |
                            zoop_size$collect_date=="2021-07-07","2022-10-15","2022-10-16")

#only select hour and then add arbitrary dates for plotting
zoop_size$Hour <- format(zoop_size$Hour, format='%H:%M')

#combine hour and date
zoop_size$Hour <- strptime(paste0(as.character(zoop_size$dates), zoop_size$Hour),format="%Y-%m-%d %H:%M")
zoop_size$Hour <- as.POSIXct(zoop_size$Hour)

#make sure zoop_size is a dataframe
zoop_size <- data.frame(zoop_size)

#order by hour for plotting
zoop_size <- zoop_size[order(zoop_size$Hour),]

#only select epi samples for DHM plots
zoop_size_epi <- zoop_size[grepl("epi",zoop_size$sample_ID) |grepl("sunrise",zoop_size$sample_ID) | grepl("sunset",zoop_size$sample_ID) | zoop_size$site_no=="BVR_l",]

#convert new dfs from tibble to dataframe 
zoop_size_DHM <- data.frame(zoop_size_epi)

#convert df from wide to long (kinda hacky way bc having problems doing this)
df1_size <- zoop_size_DHM %>% gather(metric,value,OverallMeanSize_mm_rep.mean:LecaneMeanSize_mm_rep.mean)
df2_size <- zoop_size_DHM %>% gather(metric.SE,value.SE, OverallMeanSize_mm_rep.SE:LecaneMeanSize_mm_rep.SE)

##cut and paste to merge df
zoop_size_DHM_long <- df1_size[,c(1:4,31,32)]
zoop_size_DHM_long$value.SE <- df2_size$value.SE

#drop _rep.mean from all metric names
zoop_size_DHM_long$metric <- substr(zoop_size_DHM_long$metric,1,nchar(zoop_size_DHM_long$metric)-9)

#add column for MSN #
zoop_size_DHM_long$MSN <- ifelse(zoop_size_DHM_long$collect_date=="2019-07-10" | zoop_size_DHM_long$collect_date=="2019-07-11",1,
                            ifelse(zoop_size_DHM_long$collect_date=="2019-07-24" | zoop_size_DHM_long$collect_date=="2019-07-25",2,
                                   ifelse(zoop_size_DHM_long$collect_date=="2020-08-12" | zoop_size_DHM_long$collect_date=="2020-08-13",3,
                                          ifelse(zoop_size_DHM_long$collect_date=="2021-06-15" | zoop_size_DHM_long$collect_date=="2021-06-16",4,5))))

#change facet labels
size_taxa <-c("Total","Daphniidae","Copepoda","Calanoida",
                "Cladocera", "Cyclopoida", "Rotifera", "Keratella",
                "Kellicottia","Crustacea","Bosminidae","Nauplius","Ceriodaphnia",
                "Daphnia", "Bosmina", "Ploima", "Gastropidae", "Collothecidae",
                "Conochilidae", "Synchaetidae", "Trichocercidae", "Lepadella",
                "Monostyla", "Lecane")
names(size_taxa) <- c(unique(zoop_size_DHM_long$metric))

#plot
ggplot(subset(zoop_size_DHM_long, metric %in% c("CladoceraMeanSize_mm","CopepodaMeanSize_mm","RotiferaMeanSize_mm")),
       aes(Hour,value, color=as.factor(MSN))) + 
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 11:30:00"),xmax=as.POSIXct("2022-10-15 20:41:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-15 20:42:00"),xmax=as.POSIXct("2022-10-16 06:10:00"), 
                ymin=-Inf, ymax= Inf, fill= "Midnight"),color=NA) +
  geom_rect(aes(xmin=as.POSIXct("2022-10-16 06:11:00"),xmax=as.POSIXct("2022-10-16 12:30:00"), 
                ymin=-Inf, ymax= Inf, fill= "Noon"),color=NA) +
  geom_point(size=2) + theme_bw() + 
  facet_grid(site_no~metric,scales="free_y",labeller = labeller(metric=size_taxa, site_no=sites)) + 
  xlab("")+ coord_cartesian(clip = 'off') +
  theme(text = element_text(size=8), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.9,0.95), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.7,"line")) + 
  scale_x_datetime(expand = c(0,0),labels = date_format("%H-%M",tz="EST5EDT")) +
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), 
                     labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020",
                              "15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) + 
  geom_line()+ ylab("Size (mm)") + scale_fill_manual("",values=c("#CCCCCC","white"), guide = "none")+
  geom_errorbar(aes(ymin=value-value.SE, ymax=value+value.SE), width=.2,position=position_dodge(.9))
#ggsave(file.path(getwd(),"Summer2021-DataAnalysis/Figures/BVR_MSNs_taxa_size.jpg"), width=5, height=4) 

