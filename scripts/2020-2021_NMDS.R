#new ordination w/o 2019 MSNs so we can look at both hourly + annual/MSN drivers of variability
#12Jun2023

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, ggplot2,tidyr,
               viridis, egg, ggordiplots, splancs, ggpubr, FSA, rcompanion, ggrepel)

#function to count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#read in zoop data from all 3 years
zoops2020<- read.csv('output/FCR_ZooplanktonSummary2020.csv',header = TRUE)
zoops2021<- read.csv('output/FCR_ZooplanktonSummary2021.csv',header = TRUE)

#select density cols to keep
zoops2020 <- zoops2020 %>% select("sample_ID","site_no","collect_date","DepthOfTow_m","Hour","mesh_size_μm","ZoopDensity_No.pL","Calanoida_density_NopL",
                                  "Cyclopoida_density_NopL", "Keratella_density_NopL","Kellicottia_density_NopL", "Bosmina_density_NopL",
                                  "Daphnia_density_NopL", "Ceriodaphnia_density_NopL","nauplius_density_NopL", "Collothecidae_density_NopL",
                                  "Synchaetidae_density_NopL", "Conochilidae_density_NopL")

zoops2021 <- zoops2021 %>% select("sample_ID","site_no","collect_date","DepthOfTow_m","Hour","mesh_size_μm","ZoopDensity_No.pL","Calanoida_density_NopL",
                                  "Cyclopoida_density_NopL", "Keratella_density_NopL","Kellicottia_density_NopL", "Bosmina_density_NopL",
                                  "Daphnia_density_NopL", "Ceriodaphnia_density_NopL","nauplius_density_NopL", "Collothecidae_density_NopL",
                                  "Synchaetidae_density_NopL", "Conochilidae_density_NopL")


#combine all zoop datasets
zoops <- rbind(zoops2020,zoops2021)

#ignore 20 um and horizontal trap samples
zoops <- zoops %>% filter(mesh_size_μm >20, na.rm=TRUE)

#manually change a noon and midnight hour so  summarizing below actually works
zoops$Hour[zoops$sample_ID=="B_pel_10Jul19_noon_epi_rep1"] <- "12:00"
zoops$Hour[zoops$sample_ID=="B_pel_08Jul21_midnight_epi_rep1"] <- "0:00"

#create df for temporal epi tows
zoop_epi_tows <- zoops[zoops$site_no!="FCR_50"& zoops$site_no!="BVR_d" & zoops$site_no!="BVR_dam" & 
                         (grepl("epi",zoops$sample_ID) |grepl("sunrise",zoops$sample_ID) | 
                          grepl("sunset",zoops$sample_ID) | zoops$site_no=="BVR_l"), ] %>%
                    mutate(Hour=substr(Hour,1,2)) %>% 
                    select(!c(DepthOfTow_m, sample_ID)) %>% group_by(site_no,collect_date,Hour) %>%
                    summarise(across(everything(),list(mean)))

zoop_epi_tows$time <-ifelse(zoop_epi_tows$Hour=="12" | zoop_epi_tows$Hour=="11", "noon", ifelse(
  zoop_epi_tows$Hour =="0:" | zoop_epi_tows$Hour =="23", "midnight",ifelse(zoop_epi_tows$Hour=="18"|
                                                                             zoop_epi_tows$Hour=="19" | zoop_epi_tows$Hour=="20" | zoop_epi_tows$Hour=="21", "sunset", "sunrise")))

zoop_epi_tows$site <- ifelse(substrEnd(zoop_epi_tows$site_no,1)=="l","lit","pel")

#also add columns to group times by sampling event and time and sites together
zoop_epi_tows$groups <- ifelse(zoop_epi_tows$collect_date=="2019-07-10" | zoop_epi_tows$collect_date=="2019-07-11","1",
                               ifelse(zoop_epi_tows$collect_date=="2019-07-24" | zoop_epi_tows$collect_date=="2019-07-25","2",
                                      ifelse(zoop_epi_tows$collect_date=="2020-08-12" | zoop_epi_tows$collect_date=="2020-08-13","3",
                                             ifelse(zoop_epi_tows$collect_date=="2021-06-15" | zoop_epi_tows$collect_date=="2021-06-16","4","5"))))

zoop_epi_tows$timesite <- paste0(zoop_epi_tows$time,zoop_epi_tows$site)

zoop_epi_tows$timegroup <- paste0(zoop_epi_tows$time,zoop_epi_tows$groups)

#------------------------------------------------------------------------------#
# set up data for NMDS
#relies on rank orders for ordination, no assumptions of linear relationship
zoop_temporal_dens <- zoop_epi_tows[,c(grepl("density_NopL",colnames(zoop_epi_tows)))]  

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis (jaccard might be another one to try later on to look at absolute)
zoop_temporal_dens_trans <- hellinger(zoop_temporal_dens)


#-------------------------------------------------------------------------------#
#           Averaging zoops by time and campaign/day for new NMDS               #
#-------------------------------------------------------------------------------#
#first average times for each 24-hour campaign so there are 11 points per day (basically just averaging noon and midnight)
zoop_epi_tows$order <- ifelse(zoop_epi_tows$Hour=="11" | zoop_epi_tows$Hour=="12",1, 
                              ifelse(zoop_epi_tows$Hour=="18",2, ifelse(zoop_epi_tows$Hour=="19",3,
                                                                        ifelse(zoop_epi_tows$Hour=="20",4, ifelse(zoop_epi_tows$Hour=="21",5,
                                                                                                                  ifelse(zoop_epi_tows$Hour=="0:" | zoop_epi_tows$Hour=="23",6,
                                                                                                                         ifelse(zoop_epi_tows$Hour=="4:" | zoop_epi_tows$Hour=="3:",7,
                                                                                                                                ifelse(zoop_epi_tows$Hour=="5:",8,
                                                                                                                                       ifelse(zoop_epi_tows$Hour=="6:",9,10)))))))))

#add order 11 for noon2
zoop_epi_tows$order[zoop_epi_tows$order==1 & (zoop_epi_tows$collect_date=="2019-07-10" | zoop_epi_tows$collect_date=="2019-07-24" |
                                                zoop_epi_tows$collect_date=="2020-08-12" | zoop_epi_tows$collect_date=="2021-06-15" |
                                                zoop_epi_tows$collect_date=="2021-07-07")] <- 11

#now specify whether it is noon1 or noon2
zoop_epi_tows$time[zoop_epi_tows$order==1] <- "noon1"
zoop_epi_tows$time[zoop_epi_tows$order==11] <- "noon2"

#average by MSN, site, then hour
zoop_avg <- zoop_epi_tows %>% group_by(groups,site,order) %>%
  summarise_at(vars(Calanoida_density_NopL_1:Conochilidae_density_NopL_1), list(mean = mean))

#pelagic vs littoral dfs
zoop_pel <- zoop_avg[zoop_avg$site=="pel",]
zoop_lit <- zoop_avg[zoop_avg$site=="lit",]

#only select data cols
zoop_temporal_avg_dens <- zoop_avg[,c(grepl("mean",colnames(zoop_avg)))] 

zoop_pel_dens <- zoop_pel[,c(grepl("mean",colnames(zoop_pel)))]
zoop_lit_dens <- zoop_lit[,c(grepl("mean",colnames(zoop_lit)))]

#transforming data - hellinger transformation because gives low weight to low/zero values
#converts species abundances from absolute to relative - use w/ bray curtis
zoop_temporal_dens_avg_trans <- hellinger(zoop_temporal_avg_dens)

zoop_pel_dens_trans <- hellinger(zoop_pel_dens)
zoop_lit_dens_trans <- hellinger(zoop_lit_dens)

#turn transformed community data into Euclidean distance matrix
zoop_euc <- as.matrix(vegdist(zoop_temporal_dens_avg_trans, method='euclidean'))
#zoop_bray <- as.matrix(vegdist(zoop_temporal_dens_avg_trans, method='bray'))

#scree plot to choose dimension #
dimcheckMDS(zoop_euc, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)

set.seed(3)

#now do NMDS using averages w/ 4 dimensions for consistency
NMDS_temporal_avg_bray <- metaMDS(zoop_euc, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_avg_bray$stress

#-------------------------------------------------------------------------------#
#                                 NMDS ms figs                                  #
#-------------------------------------------------------------------------------#

ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),choices = c(1,2),type = "n")
sites <- gg_ordiplot(ord, zoop_avg$site, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_site <- sites$plot + geom_point() + theme_bw() + 
  geom_polygon(data = sites$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=sites$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=c("#882255","#3399CC")) +
  # xlim(c(-0.63, 0.6)) + ylim(c(-0.7,0.64)) +
  #scale_x_continuous(labels=c('-0.4', '-0.2', '0.0', '0.2', '')) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
  annotate("text", x=-0.27, y=0.4, label= "a: sites", 
           fontface = "italic", size = 3) +
  scale_fill_manual("",values=c("#882255","#3399CC"))+
  scale_color_manual("",values=c("#882255","#3399CC"),
                     label=c('littoral','pelagic')) 

days <- gg_ordiplot(ord, zoop_avg$groups, kind = "ehull", 
                    ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day <- days$plot + geom_point() + theme_bw() + geom_path() + ylab(NULL) +
  geom_polygon(data = days$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=days$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=c("#EFECBF","#DB9B5A","#C7522B")) +
  # xlim(c(-0.63, 0.6)) + ylim(c(-0.7,0.64)) +
#  scale_x_continuous(labels=c('-0.4', '-0.2', '0.0', '0.2', '')) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,-0.1), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  annotate("text", x=-0.05, y=0.4, label= "b: sampling campaigns", 
           fontface = "italic", size = 3) +
  guides(fill="none", color = guide_legend(ncol=2)) +
  scale_fill_manual("",values=c("#EFECBF","#DB9B5A","#C7522B"))+
  scale_color_manual("",values=c("#EFECBF","#DB9B5A","#C7522B"),
                     label=c('12-13 Aug 2020','15-16 Jun 2021', '7-8 Jul 2021'))


hours <- gg_ordiplot(ord, zoop_avg$order, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 
#order hours properly
hours$df_hull$Group <- factor(hours$df_hull$Group, levels = c(unique(hours$df_hull$Group)))

NMDS_hour <- hours$plot + geom_point() + theme_bw() + geom_path() + ylab(NULL) +
  geom_polygon(data = hours$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=hours$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=hcl.colors(11,"sunset")) +
  #scale_x_continuous(labels=c('-0.4', '-0.2', '0.0', '0.2', '')) +
  theme(text = element_text(size=7), axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0,0,-0.1), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
  annotate("text", x=-0.1, y=0.4, label= "c: hours of the day",
           fontface = "italic", size=3) +
  scale_fill_manual("",values=hcl.colors(11,"sunset"))+
  scale_color_manual("",values=hcl.colors(11,"sunset"),
                     label=c('12pm','6pm','7pm','8pm','9pm','12am',
                             '4am','5am','6am','7am','12pm'))

fig5 <- egg::ggarrange(NMDS_site, NMDS_day, NMDS_hour, nrow=1)
#ggsave("figures/NMDS_multipanel_2v1.jpg", fig5, width=5, height=3) 


#------------------------------------------------------------------------------#
#now read in daily thermistor data

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/725/3/a9a7ff6fe8dc20f7a8f89447d4dc2038" 
infile1 <- tempfile()
options(timeout=3000)
try(download.file(inUrl1,infile1,method="auto"))

sensors <-read.csv(infile1) %>%
  mutate(date = as.Date(DateTime)) %>%
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-$d %H:%M:%S")) %>%
  mutate(hour = lubridate::hour(DateTime)) %>%
  filter(Reservoir =="BVR" & Site == 50 &
    date %in% c(as.Date("2020-08-12"), as.Date("2020-08-13"),
                as.Date("2021-06-15"), as.Date("2021-06-16"), as.Date("2021-07-07"), 
                as.Date("2021-07-08")) &  #can only do TN/TP for all 5 MSNs
          hour %in% c(0, 4, 5, 6, 7, 12, 18, 19, 20, 21))

#average by hour and date
sensors_avg <- sensors %>% group_by(date, hour) %>% summarise_at(vars(ThermistorTemp_C_1:EXODepth_m), mean, na.rm=T)

#only select columns needed for env fit
sensors_avg <- sensors_avg %>% select(date, hour, ThermistorTemp_C_3, ThermistorTemp_C_10, 
                                      EXOSpCond_uScm_1.5, EXOChla_ugL_1.5, EXOTDS_mgL_1.5,
                                      EXOfDOM_RFU_1.5)

#only select the 11 hours for each MSN
sensors_avg <- sensors_avg[c(6:16,26:36,46:56),]

#add new column for order to join below
sensors_avg$order <- rep(1:11,3)

#add groups
sensors_avg$groups <- as.character(ifelse(sensors_avg$date=="2020-08-12" | sensors_avg$date=="2020-08-13", 3,
                             ifelse(sensors_avg$date=="2021-06-15" | sensors_avg$date=="2021-06-16",4,5)))


#convert from wide to long
sensors_avg_long <- sensors_avg %>% pivot_longer(cols = ThermistorTemp_C_3:EXOfDOM_RFU_1.5, names_to = "variable")

#order df so I can match NMDS data to rows easier
sensors_avg_long <- sensors_avg_long[order(sensors_avg_long$variable, sensors_avg_long$order),]

#add NMDS2 col - just pelagic site
sensors_avg_long$NMDS2 <- ifelse(sensors_avg_long$order==1, hours$df_ord$y[hours$df_ord$Group==1][c(2,4,6)],
                                 ifelse(sensors_avg_long$order==2, hours$df_ord$y[hours$df_ord$Group==2][c(2,4,6)],
                                        ifelse(sensors_avg_long$order==3, hours$df_ord$y[hours$df_ord$Group==3][c(2,4,6)],
                                               ifelse(sensors_avg_long$order==4, hours$df_ord$y[hours$df_ord$Group==4][c(2,4,6)],
                                                      ifelse(sensors_avg_long$order==5, hours$df_ord$y[hours$df_ord$Group==5][c(2,4,6)],
                                                             ifelse(sensors_avg_long$order==6, hours$df_ord$y[hours$df_ord$Group==6][c(2,4,6)],
                                                                    ifelse(sensors_avg_long$order==7, hours$df_ord$y[hours$df_ord$Group==7][c(2,4,6)],
                                                                           ifelse(sensors_avg_long$order==8, hours$df_ord$y[hours$df_ord$Group==8][c(2,4,6)],
                                                                                  ifelse(sensors_avg_long$order==9, hours$df_ord$y[hours$df_ord$Group==9][c(2,4,6)],
                                                                                         ifelse(sensors_avg_long$order==10, hours$df_ord$y[hours$df_ord$Group==10][c(2,4,6)],
                                                                                                hours$df_ord$y[hours$df_ord$Group==11][c(2,4,6)]))))))))))

#multipanel plot
hourly_driver_NMDS <- ggplot(data=sensors_avg_long, aes(NMDS2, value, color=as.factor(hour))) + geom_point() +
  facet_wrap(~variable, scales = "free_y") +
  scale_color_manual("Hour",values=hcl.colors(11,"sunset"), 
                     labels=c("12pm", "6pm", "7pm", "8pm", "9pm", "12am", 
                              "4am", "5am", "6am", "7am", "12pm"), guide=guide_legend(order=1)) +
  theme(text = element_text(size=5), axis.text = element_text(size=4, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), legend.position = "top", 
        legend.spacing = unit(-0.5, 'cm'), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line"))
#ggsave("figures/hourly_driver_vs_NMDS2.jpg", hourly_driver_NMDS, width=3, height=3) 

hourly_driver_NMDS_byday <- ggplot(data=sensors_avg_long, aes(NMDS2, value, color=as.factor(groups))) + geom_point() +
  facet_wrap(~variable, scales = "free_y") +
  scale_color_manual("Sampling Campaign",values=c("#F2E2B0","#DEA868","#C7522B"), 
                     labels=c("12-13 Aug 2020", "15-16 Jun 2021", "7-8 Jul 2021"), 
                     guide=guide_legend(order=1)) +
  theme(text = element_text(size=5), axis.text = element_text(size=4, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), legend.position = "top", 
        legend.spacing = unit(-0.5, 'cm'), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line"))
#ggsave("figures/hourly_driver_vs_NMDS2_byday.jpg", hourly_driver_NMDS_byday, width=3, height=3) 

