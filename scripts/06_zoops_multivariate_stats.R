#multivariate stats for 2019 - 2021 zoop data
#created 25Nov2021

source("scripts/01_install.R")

#select density cols to keep
zoop_summary <- zoop_summary |>  
  select("Reservoir","Site","DateTime","StartDepth_m","EndDepth_m", "Rep",
         "CollectionMethod","Taxon", "Density_IndPerL")

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
                    
#convert from long to wide for NMDS
zoop_summary_wide <- zoop_summary |> 
  group_by(Reservoir, Site, DateTime, CollectionMethod, Taxon) |> 
  pivot_wider(names_from = Taxon, values_from = Density_IndPerL,
              names_glue = "{Taxon}_density_NopL")

#create df for temporal epi tows
zoop_epi_tows <- zoop_summary_wide |> 
  filter(Site != 49, StartDepth_m %in% c(2,4), CollectionMethod == "Tow") |> 
  mutate(Hour=hour(DateTime)) |> 
  select(!c(StartDepth_m, EndDepth_m, CollectionMethod))

#manually remove collection method and reservoir
zoop_epi_tows <- zoop_epi_tows[,!colnames(zoop_epi_tows) %in%
                                     c("CollectionMethod","Reservoir")]

#change a couple hours so they group properly
zoop_epi_tows$Hour[as.character(zoop_epi_tows$DateTime)=="2019-07-10 11:30:00"] <- 12
zoop_epi_tows$Hour[as.character(zoop_epi_tows$DateTime)=="2019-07-10 19:04:00"] <- 18
zoop_epi_tows$Hour[as.character(zoop_epi_tows$DateTime)=="2021-07-08 23:54:00"] <- 0

#now average the reps and drop rep and date cols
zoop_epi_tows <- zoop_epi_tows |> 
  select(!c(Rep)) |> 
  group_by(Site, date, Hour)  |> 
  summarise_at(vars(Copepoda_density_NopL:Lecane_density_NopL),list(mean=mean))

#new col for time
zoop_epi_tows$time <-ifelse(zoop_epi_tows$Hour=="12" | 
                            zoop_epi_tows$Hour=="11", "noon", ifelse(
                            zoop_epi_tows$Hour =="0" | zoop_epi_tows$Hour =="23", 
                            "midnight",ifelse(zoop_epi_tows$Hour=="18"|
                            zoop_epi_tows$Hour=="19" | zoop_epi_tows$Hour=="20" | 
                            zoop_epi_tows$Hour=="21", "sunset", "sunrise")))

zoop_epi_tows$site <- ifelse(zoop_epi_tows$Site==51,"lit","pel")

#also add columns to group times by sampling event and time and sites together
zoop_epi_tows$groups <- 
  ifelse(zoop_epi_tows$date=="2019-07-10" | 
         zoop_epi_tows$date=="2019-07-11","1",
  ifelse(zoop_epi_tows$date=="2019-07-24" | 
         zoop_epi_tows$date=="2019-07-25","2",
  ifelse(zoop_epi_tows$date=="2020-08-12" | 
         zoop_epi_tows$date=="2020-08-13","3",
  ifelse(zoop_epi_tows$date=="2021-06-15" | 
         zoop_epi_tows$date=="2021-06-16","4","5"))))


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
                              ifelse(zoop_epi_tows$Hour=="18",2, 
                                     ifelse(zoop_epi_tows$Hour=="19",3,
                              ifelse(zoop_epi_tows$Hour=="20",4, 
                                     ifelse(zoop_epi_tows$Hour=="21",5,
                              ifelse(zoop_epi_tows$Hour=="0" | 
                                       zoop_epi_tows$Hour=="23",6,
                              ifelse(zoop_epi_tows$Hour=="4" | 
                                       zoop_epi_tows$Hour=="3",7,
                              ifelse(zoop_epi_tows$Hour=="5",8,
                              ifelse(zoop_epi_tows$Hour=="6",9,10)))))))))

#add order 11 for noon2
zoop_epi_tows$order[zoop_epi_tows$order==1 & 
                   (zoop_epi_tows$date=="2019-07-10" | 
                   zoop_epi_tows$date=="2019-07-24" |
                   zoop_epi_tows$date=="2020-08-12" | 
                   zoop_epi_tows$date=="2021-06-15" |
                   zoop_epi_tows$date=="2021-07-07")] <- 11

#now specify whether it is noon1 or noon2
zoop_epi_tows$time[zoop_epi_tows$order==1] <- "noon1"
zoop_epi_tows$time[zoop_epi_tows$order==11] <- "noon2"

#average by MSN, site, then hour
zoop_avg <- zoop_epi_tows %>% group_by(groups,site,order) %>%
  summarise_at(vars(Copepoda_density_NopL_mean:Lecane_density_NopL_mean), 
               list(mean = mean))

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

#scree plot to choose dimension # (Figure S2)
#jpeg("figures/scree.jpg") 
dimcheckMDS(zoop_euc, distance = "bray", k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(15)

#now do NMDS using averages w/ 4 dimensions for consistency
NMDS_temporal_avg_bray <- metaMDS(zoop_euc, distance='bray', k=4, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_avg_bray$stress

#-------------------------------------------------------------------------------#
#                                 NMDS ms figs                                  #
#-------------------------------------------------------------------------------#

ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),
                choices = c(1,2),type = "n")
sites <- gg_ordiplot(ord, zoop_avg$site, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_site <- sites$plot + geom_point() + theme_bw() + 
  geom_polygon(data = sites$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=sites$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=c("#882255","#3399CC")) +
  theme(text = element_text(size=7), 
        axis.text = element_text(size=7, color="black"), 
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
        annotate("text", x=-0.15, y=0.5, label= "a: sites", 
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
             fill=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B")) +
  theme(text = element_text(size=7), 
        axis.text = element_text(size=7, color="black"), 
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
        annotate("text", x=-0.02, y=0.5, label= "b: sampling days", 
                 fontface = "italic", size = 3) +
        guides(fill="none", color = guide_legend(ncol=2)) +
        scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
        scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"),
                     label=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020',
                             '15-16 Jun 2021', '7-8 Jul 2021'))


hours <- gg_ordiplot(ord, zoop_avg$order, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 
#order hours properly
hours$df_hull$Group <- factor(hours$df_hull$Group, levels = 
                                c(unique(hours$df_hull$Group)))

NMDS_hour <- hours$plot + geom_point() + theme_bw() + geom_path() + ylab(NULL) +
  geom_polygon(data = hours$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=hours$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=hcl.colors(11,"sunset")) +
  theme(text = element_text(size=7), 
        axis.text = element_text(size=7, color="black"), 
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
        annotate("text", x=0, y=0.5, label= "c: hours of the day",
                 fontface = "italic", size=3) +
        scale_fill_manual("",values=hcl.colors(11,"sunset"))+
        scale_color_manual("",values=hcl.colors(11,"sunset"),
                     label=c('12pm','6pm','7pm','8pm','9pm','12am',
                             '4am','5am','6am','7am','12pm'))

fig5 <- egg::ggarrange(NMDS_site, NMDS_day, NMDS_hour, nrow=1)
#ggsave("figures/NMDS_multipanel_2v1.jpg",fig5, width=5, height=3.5) 

#-------------------------------------------------------------------------------#
#                     Calculating euclidean distance                            #
#-------------------------------------------------------------------------------#
#step 1: use Euclidean distance matrix from transformed community data (zoop_euc)

#order zoop epi tows by hour, MSN, and site
zoop_epi_tows <- zoop_epi_tows |> dplyr::arrange(site, groups, order)

#convert ED matrix back into distance structure for next steps
zoop_euc <- vegdist(zoop_temporal_dens_avg_trans, method='euclidean', 
                    upper = TRUE)

#Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
centroids_sites <- betadisper(zoop_euc, group = as.factor(zoop_epi_tows$site), 
                              type="centroid")
centroids_hours <- betadisper(zoop_euc, group = as.factor(zoop_epi_tows$order), 
                              type="centroid")
centroids_days <-  betadisper(zoop_euc, group = as.factor(zoop_epi_tows$groups), 
                              type="centroid")

#-------------------------------------------------------------------------------#
#METHOD 1:average distance of each point to polygon centroid (dispersion approach)
#METHOD 2: average distance between all combinations of centroids (pairwise approach)

#"bootstrapping" or randomly select 10 points per group and 2 groups and then calculating values using methods above

#create a df to store all the different variability values
var_results <- data.frame("site_disp"=rep(NA,500))

set.seed(1)

for (i in 1:500){ 
  #randomly select 10 points in each group
  ord_sub <- sample(unique(zoop_epi_tows$order), 2)
  groups_sub <-  sample(unique(zoop_epi_tows$groups), 2)
  
  zoop_sub <-  zoop_epi_tows |> group_by(order) |> 
               filter(order %in% c(ord_sub), groups %in% c(groups_sub)) |>
               slice_sample(n=10)
  
  #only select data cols
  zoop_dens_sub <- zoop_sub[,c(grepl("density",colnames(zoop_sub)))] 
  
  #hellinger transformation
  zoop_dens_sub_trans <- hellinger(zoop_dens_sub)
  
  #convert ED matrix back into distance structure for next steps
  zoop_euc_sub <- vegdist(zoop_dens_sub_trans, method='euclidean', 
                          upper = TRUE)
  
  #Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
  centroids_sites_sub <- betadisper(zoop_euc_sub, group = as.factor(zoop_sub$site), 
                                    type="centroid")
  centroids_days_sub <-  betadisper(zoop_euc_sub, group = as.factor(zoop_sub$groups), 
                                    type="centroid")
  centroids_hours_sub <- betadisper(zoop_euc_sub, group = as.factor(zoop_sub$order), 
                                    type="centroid")
  
  #distance calcs
  var_results$site_disp[i] <- mean(centroids_sites_sub$group.distances)
  var_results$day_disp[i] <- mean(centroids_days_sub$group.distances)
  var_results$hour_disp[i] <- mean(centroids_hours_sub$group.distances)
  
  var_results$site_pair[i] <- mean(dist(centroids_sites_sub$centroids))
  var_results$day_pair[i] <- mean(dist(centroids_days_sub$centroids))
  var_results$hour_pair[i] <- mean(dist(centroids_hours_sub$centroids))
  
  }
  
#-------------------------------------------------------------------------------#
#Kruskal-wallis test to determine if group means are significant

#first convert wide to long
disp_df <- var_results[,grepl("disp",colnames(var_results))] |>  
  pivot_longer(everything(), names_to="group")
pair_df <- var_results[,grepl("pair",colnames(var_results))] |> 
  pivot_longer(everything(), names_to="group")

#now kw test
kw_disp <- kruskal.test(value ~ group, data = disp_df) #significant
kw_pair <- kruskal.test(value ~ group, data = pair_df) #significant

#now dunn test to determine which groups are different from each other
dunnTest(value ~ as.factor(group),
         data=disp_df,
         method="bonferroni")

dunnTest(value ~ as.factor(group),
         data=pair_df,
         method="bonferroni")

#sites, days, and hours are all different from each other!

disp_box <- ggboxplot(disp_df, x = "group", y = "value", 
          fill = "group", palette = c("#A4C6B8", "#81858B", "#5E435D"),
          order = c("site_disp", "day_disp", "hour_disp"),
          ylab = "Dispersion", xlab = "") +
  theme(text = element_text(size=7),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.2,0,-0.5,0), 'lines')) +
  annotate("text",label=c("b","c","a"), x=c(1.1,2.1,3.1),
           y=c(mean(disp_df$value[disp_df$group=="site_disp"]) + 
                 sd(disp_df$value[disp_df$group=="site_disp"]),
               mean(disp_df$value[disp_df$group=="day_disp"]) + 
                 sd(disp_df$value[disp_df$group=="day_disp"]),
               mean(disp_df$value[disp_df$group=="hour_disp"]) + 
                 sd(disp_df$value[disp_df$group=="hour_disp"]))) +
  guides(fill = "none") 

pair_box <- ggboxplot(pair_df, x = "group", y = "value", 
          fill = "group", palette = c("#A4C6B8", "#81858B", "#5E435D"),
          order = c("site_pair", "day_pair", "hour_pair"),
          ylab = "Location", xlab = "") +
  scale_x_discrete(labels = c("sites", "sampling \n\ days", "hours of \n\ the day")) +
  theme(text = element_text(size=7),
        plot.margin = unit(c(-0.5,0,0,0), 'lines')) +
  annotate("text",label=c("b","a","c"), x=c(1.1,2.1,3.1),
           y=c(mean(pair_df$value[pair_df$group=="site_pair"]) + 
                 sd(pair_df$value[pair_df$group=="site_pair"]),
               mean(pair_df$value[pair_df$group=="day_pair"]) + 
                 sd(pair_df$value[pair_df$group=="day_pair"]),
               mean(pair_df$value[pair_df$group=="hour_pair"]) + 
                 sd(pair_df$value[pair_df$group=="hour_pair"]))) +
  guides(fill = "none") 

among_scales <- egg::ggarrange(disp_box, pair_box, nrow=2)
#ggsave("figures/among_variability_boxplots.jpg",among_scales, width=3, height=4) 

#create table for kw test results
kw_results <- data.frame("Variability" = c(" ", "Dispersion", " ", " ", "Location", " "),
                         "Scale" = rep(c("Sites", "Sampling campaigns", "Hours of the day"),2),
                         "n" = c(rep(500,6)),
                         "mean" = c(round(mean(var_results$site_disp),2), 
                                    round(mean(var_results$day_disp),2),
                                    round(mean(var_results$hour_disp),2), 
                                    round(mean(var_results$site_pair),2),
                                    round(mean(var_results$day_pair),2), 
                                    round(mean(var_results$hour_pair),2)),
                         "sd" = c(round(sd(var_results$site_disp),2), 
                                  round(sd(var_results$day_disp),2),
                                  round(sd(var_results$hour_disp),2), 
                                  round(sd(var_results$site_pair),2),
                                  round(sd(var_results$day_pair),2), 
                                  round(sd(var_results$hour_pair),2)),
                         "df" = c(" ",kw_disp$parameter, " ", " ", 
                                  kw_pair$parameter, " "),
                         "χ2" = c(" ",round(kw_disp$statistic,3), " ", " ", 
                                  round(kw_pair$statistic,3), " "),
                         "p-value" = c(" ",kw_disp$p.value, " ", " ", 
                                       kw_pair$p.value, " "))
#write.csv(kw_results, "output/Euclidean_distances_bootstrapped_kw_results.csv", row.names = FALSE)

#-------------------------------------------------------------------------------#
#dfs to calculate significance within sites, days, and hours
within_site_dist <- data.frame("group" = c(rep("lit",55),rep("pel",55)),
                                "dist" = c(centroids_sites$distances[zoop_epi_tows$site=="lit"],
                                          centroids_sites$distances[zoop_epi_tows$site=="pel"]))

within_day_dist <- data.frame("group" = c(rep("day1",22),rep("day2",22),rep("day3",22),
                                          rep("day4",22),rep("day5",22)),
                              "dist" = c(centroids_days$distances[zoop_epi_tows$groups==1],
                                         centroids_days$distances[zoop_epi_tows$groups==2],
                                         centroids_days$distances[zoop_epi_tows$groups==3],
                                         centroids_days$distances[zoop_epi_tows$groups==4],
                                         centroids_days$distances[zoop_epi_tows$groups==5]))

within_hour_dist <- data.frame("group" = c(rep("hour1",10),rep("hour2",10),rep("hour3",10),
                                           rep("hour4",10), rep("hour5",10), rep("hour6",10),
                                           rep("hour7",10), rep("hour8",10), rep("hour9",10),
                                           rep("hour10",10), rep("hour11",10)),
                               "dist" = c(centroids_hours$distances[zoop_epi_tows$order==1],
                                          centroids_hours$distances[zoop_epi_tows$order==2],
                                          centroids_hours$distances[zoop_epi_tows$order==3],
                                          centroids_hours$distances[zoop_epi_tows$order==4],
                                          centroids_hours$distances[zoop_epi_tows$order==5],
                                          centroids_hours$distances[zoop_epi_tows$order==6],
                                          centroids_hours$distances[zoop_epi_tows$order==7],
                                          centroids_hours$distances[zoop_epi_tows$order==8],
                                          centroids_hours$distances[zoop_epi_tows$order==9],
                                          centroids_hours$distances[zoop_epi_tows$order==10],
                                          centroids_hours$distances[zoop_epi_tows$order==11]))


#now kw tests for significance 
kw_sites <- kruskal.test(dist ~ group, data = within_site_dist) #sig! - so pelagic and littoral are different
kw_days <- kruskal.test(dist ~ group, data = within_day_dist) #not sig
kw_hours <- kruskal.test(dist ~ group, data = within_hour_dist) #nope

#dunn post-hoc test
dunn_within_days <- dunnTest(dist ~ as.factor(group),
                        data=within_day_dist,
                        method="bonferroni")

cldList(P.adj ~ Comparison, data=dunn_within_days$res, threshold = 0.05)

#-------------------------------------------------------------------------------
#within boxplot for sites, days, and hours

site_box <- ggboxplot(within_site_dist, x = "group", y = "dist", 
                     fill = "group", palette = c("#882255","#3399CC"),
                     order = c("lit", "pel"),
                     ylab = "Dispersion", xlab = "") +
           ylim(c(0,1)) +
           scale_x_discrete(labels = c("littoral", "pelagic")) +
           theme(text = element_text(size=8),
                 plot.margin = unit(c(0,-0.4,0,0), 'lines'),
                 axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8)) +
           annotate("text",label=c("a","b"), x=c(1.2,2.2), size=4,
                    y=c(0.6, 0.437)) +
           annotate("text", x=1.3, y=1, label= "a: sites",
                    fontface = "italic", size=3) +
           guides (fill = "none")

day_box <- ggboxplot(within_day_dist, x = "group", y = "dist", 
              fill = "group", palette = c("#008585", "#9BBAA0", "#F2E2B0","#DEA868","#C7522B"),
              order = c("day1", "day2", "day3","day4","day5"),
              ylab = "", xlab = "") + ylim(c(0,1)) +
          scale_x_discrete(labels = c("10-11 Jul 2019", "24-25 Jul 2019", 
                              "12-13 Aug 2020", "15-16 Jun 2021", "7-8 Jul 2021")) +
          theme(text = element_text(size=8),
                plot.margin = unit(c(0,-0.4,0,-0.4), 'lines'),
                axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8),
                axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          annotate("text", x=2.2, y=1, label= "b: sampling days",
                    fontface = "italic", size=3) +
          guides (fill = "none")

hour_box <- ggboxplot(within_hour_dist, x = "group", y = "dist", 
                      fill = "group", palette = hcl.colors(11,"sunset"),
                      order = c("hour1", "hour2","hour3","hour4","hour5",
                                "hour6","hour7","hour8","hour9","hour10","hour11"),
                      ylab = "", xlab = "") + ylim(c(0,1)) +
           scale_x_discrete(labels = c("12pm", "6pm", "7pm", "8pm", "9pm", "12am", 
                                       "4am", "5am", "6am", "7am", "12pm")) +
           theme(text = element_text(size=8),
                 plot.margin = unit(c(0,0,0,-0.4), 'lines'),
                 axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8),
                 axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           annotate("text", x=2.7, y=1, label= "c: hours of the day",
                    fontface = "italic", size=3) +
           guides (fill = "none")

within_scales <- egg::ggarrange(site_box, day_box, hour_box, nrow=1, widths = c(0.9,2.2,4))
#ggsave("figures/within_variability_boxplots.jpg", within_scales, width=5, height=4) 


#create table for within scale kw test results
kw_results_disp <- data.frame("Group" = c("Littoral", "Pelagic", "10-11 Jul 2019", 
                                          "24-25 Jul 2019", "12-13 Aug 2020", 
                                          "15-16 Jun 2021", "7-8 Jul 2021", "12pm",
              "6pm", "7pm", "8pm", "9pm", "12am", "4am", "5am", "6am", "7am", "12pm"),
  "n" = c(55, 55, 22, 22, 22, 22, 22, rep(10,11)),
  "mean" = c(round(mean(within_site_dist$dist[within_site_dist$group=="lit"]),2), 
             round(mean(within_site_dist$dist[within_site_dist$group=="pel"]),2),
             round(mean(within_day_dist$dist[within_day_dist$group=="day1"]),2),
             round(mean(within_day_dist$dist[within_day_dist$group=="day2"]),2),
             round(mean(within_day_dist$dist[within_day_dist$group=="day3"]),2),
             round(mean(within_day_dist$dist[within_day_dist$group=="day4"]),2),
             round(mean(within_day_dist$dist[within_day_dist$group=="day5"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour1"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour2"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour3"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour4"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour5"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour6"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour7"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour8"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour9"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour10"]),2),
             round(mean(within_hour_dist$dist[within_hour_dist$group=="hour11"]),2)),
  "sd" = c(round(sd(within_site_dist$dist[within_site_dist$group=="lit"]),2), 
           round(sd(within_site_dist$dist[within_site_dist$group=="pel"]),2),
           round(sd(within_day_dist$dist[within_day_dist$group=="day1"]),2),
           round(sd(within_day_dist$dist[within_day_dist$group=="day2"]),2),
           round(sd(within_day_dist$dist[within_day_dist$group=="day3"]),2),
           round(sd(within_day_dist$dist[within_day_dist$group=="day4"]),2),
           round(sd(within_day_dist$dist[within_day_dist$group=="day5"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour1"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour2"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour3"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour4"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour5"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour6"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour7"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour8"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour9"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour10"]),2),
           round(sd(within_hour_dist$dist[within_hour_dist$group=="hour11"]),2)),
  "df" = c(kw_sites$parameter, " ", rep(" ",2), kw_days$parameter, rep(" ",2),
           rep(" ", 5), kw_hours$parameter, rep(" ", 5)),
  "χ2" = c(round(kw_sites$statistic,3), " ", rep(" ",2), round(kw_days$statistic,3), rep(" ",2),
           rep(" ", 5), round(kw_hours$statistic,3), rep(" ", 5)),
  "p-value" = c(kw_sites$p.value, " ", rep(" ",2), kw_days$p.value, rep(" ",2),
           rep(" ", 5), kw_hours$p.value, rep(" ", 5)))
#write.csv(kw_results_disp, "output/within_group_dispersion_kw_results.csv",row.names = FALSE)

#------------------------------------------------------------------------------#
#pull in driver data --> can't use ysi or fp because only have 2 or 3 casts out of 5

msn_dates <- c("2019-07-10", "2019-07-11", "2019-07-24",
               "2019-07-25", "2020-08-12", "2020-08-13",
               "2021-06-15", "2021-06-16", "2021-07-07", 
               "2021-07-08")

#chem
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

chem <-read.csv(infile1,header=T) |>
  mutate(DateTime = as.Date(DateTime)) |>
  filter(Reservoir =="BVR" & Site==50 & 
           DateTime %in% c(msn_dates)) #can only do TN/TP for all 5 MSNs

#add msn #
chem$msn <- ifelse(chem$DateTime=="2019-07-10" | chem$DateTime=="2019-07-11","1",
                        ifelse(chem$DateTime=="2019-07-24" | 
                                 chem$DateTime=="2019-07-25","2",
                               ifelse(chem$DateTime=="2020-08-12" | 
                                        chem$DateTime=="2020-08-13","3",
                                      ifelse(chem$DateTime=="2021-06-15" | 
                                               chem$DateTime=="2021-06-16","4","5"))))

#average across msn and depth
chem_avg <- chem |> group_by(Reservoir, Depth_m, msn) |>  
  summarise(across(everything(), list(mean)))

#--------------------------------------#
#Secchi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,timeout = max(300, getOption("timeout"))))

secchi <-read.csv(infile1,header=T) |>
  mutate(DateTime = as.Date(DateTime)) |>
  filter(Reservoir =="BVR" & Site==50 & 
           DateTime %in% c(msn_dates)) #can only do TN/TP for all 5 MSNs

#add secchi for 5th MSN - 1.8m
secchi[5,] <- c("BVR",50,"2021-07-08",1.8,0,0)

#--------------------------------------#
#ctd
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |>
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(Reservoir =="BVR" & Depth_m >0 &
           DateTime %in% c(msn_dates, "2021-07-12")) #ctd cast from 4 days after MSN 5

#select every 0.5m from casts
ctd_final_temp <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) |>
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) |>
  dplyr::summarise(value = mean(Temp_C)) |>
  dplyr::rename(depth = rdepth) 

#do the same for DO values now
ctd_final_DO <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) |>
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) |>
  dplyr::summarise(value = mean(DO_mgL)) |>
  dplyr::rename(depth = rdepth) |>
  dplyr::group_by(Reservoir,Site, DateTime) |>
  dplyr::mutate(oxy = min(depth[value<=2], na.rm=TRUE))

depths <- ctd_final_DO$oxy[(ctd_final_DO$DateTime=="2019-07-10" | 
                              ctd_final_DO$DateTime=="2019-07-24" |
                              ctd_final_DO$DateTime=="2020-08-12" | 
                              ctd_final_DO$DateTime=="2021-06-16" | 
                              ctd_final_DO$DateTime=="2021-07-12") & 
                              ctd_final_DO$Reservoir=="BVR" & 
                              ctd_final_DO$Site==50 & ctd_final_DO$depth>0]

#oxygen plots for all 5 MSNs (Figure S1)
ggplot(subset(ctd_final_DO, depth > 0 & Reservoir=="BVR" & Site==50 & 
                DateTime %in% c(as.Date("2019-07-10"), as.Date("2019-07-24"),
                as.Date("2020-08-12"), as.Date("2021-06-16"), as.Date("2021-07-12"))), 
       aes(value,depth,color=as.factor(DateTime))) + 
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin= depths, ymax=Inf), fill="red",alpha=0.03) +
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), guide="none") +  
  xlab("DO (mg/L)") + ylab("Depth (m)") +
  ylim(10,0) + geom_point() + geom_path() + 
  theme_bw() + facet_wrap(~DateTime, ncol=5) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        legend.position = c(0.76,0.02), legend.background = element_blank(),
        legend.direction = "horizontal", legend.title = element_text(size = 4.5),
        panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,0,0), "cm"),legend.key.size = unit(0.5, "lines"),
        panel.grid.major = element_blank(), axis.text.y = element_text(size=6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), 
        legend.text  = element_text(size = 4.5), panel.spacing=unit(0, "cm")) 
#ggsave("figures/DO_profiles.jpg"), width=4, height=3) 


#calculate thermocline depth and oxycline depth by date
ctd_thermo_depth <- ctd_final_temp |> group_by(DateTime) |> 
  mutate(therm_depth = thermo.depth(value,depth))

ctd_oxy_depth <- ctd_final_DO %>% group_by(DateTime) |>
  mutate(oxy_depth = thermo.depth(value,depth))

#add msn # 
ctd_thermo_depth$msn <- ifelse(ctd_thermo_depth$DateTime=="2019-07-10" | 
                            ctd_thermo_depth$DateTime=="2019-07-11","1",
                        ifelse(ctd_thermo_depth$DateTime=="2019-07-24" | 
                            ctd_thermo_depth$DateTime=="2019-07-25","2",
                        ifelse(ctd_thermo_depth$DateTime=="2020-08-12" | 
                            ctd_thermo_depth$DateTime=="2020-08-13","3",
                        ifelse(ctd_thermo_depth$DateTime=="2021-06-15" | 
                            ctd_thermo_depth$DateTime=="2021-06-16","4","5"))))

ctd_oxy_depth$msn <- ifelse(ctd_oxy_depth$DateTime=="2019-07-10" | 
                            ctd_oxy_depth$DateTime=="2019-07-11","1",
                     ifelse(ctd_oxy_depth$DateTime=="2019-07-24" | 
                            ctd_oxy_depth$DateTime=="2019-07-25","2",
                     ifelse(ctd_oxy_depth$DateTime=="2020-08-12" | 
                            ctd_oxy_depth$DateTime=="2020-08-13","3",
                     ifelse(ctd_oxy_depth$DateTime=="2021-06-15" | 
                            ctd_oxy_depth$DateTime=="2021-06-16","4","5"))))

depths = c(0.1, seq(1, 10, by = 1))
ctd.final<-data.frame() 

for (i in 1:length(depths)){
  ctd_layer <- ctd |> group_by(DateTime) |> slice(which.min(abs(
    as.numeric(Depth_m) - depths[i])))
  # Bind each of the data layers together.
  ctd.final = bind_rows(ctd.final, ctd_layer)
}

#round to one decimal place
ctd.final$Depth_m <- round(ctd.final$Depth_m,1)

#add msn # 
ctd.final$msn <- ifelse(ctd.final$DateTime=="2019-07-10" | 
                            ctd.final$DateTime=="2019-07-11","1",
                 ifelse(ctd.final$DateTime=="2019-07-24" | 
                            ctd.final$DateTime=="2019-07-25","2",
                 ifelse(ctd.final$DateTime=="2020-08-12" | 
                            ctd.final$DateTime=="2020-08-13","3",
                 ifelse(ctd.final$DateTime=="2021-06-15" | 
                            ctd.final$DateTime=="2021-06-16","4","5"))))

#average across msn and depth
ctd.final_avg <- ctd.final |> group_by(Reservoir, Depth_m, msn) |> 
  summarise(across(everything(), list(mean)))

#read in migration metrics
migration_metrics <- read.csv("output/migration_metrics.csv",header=T)

#read in met data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/7/02d36541de9088f2dd99d79dc3a7a853" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

#list hours in day vs night
day <- c(6:20)
night <- c(21:23,0:5)

#monthly met for jul 2019, Aug 2020, Jun 2021, Jul 2021
monthly_met <- read.csv(infile1,header=T) |> 
  dplyr::mutate(my = paste0(format(as.Date(DateTime), "%m"),
                      format(as.Date(DateTime), "%y"))) |> 
  dplyr::filter(my %in% c("0719","0820","0621","0721")) |> 
  dplyr::select(my, AirTemp_C_Average,
                WindSpeed_Average_m_s, WindDir_degrees) |> 
  dplyr::group_by(my) |> 
  dplyr::summarise(AirTemp_month = mean(AirTemp_C_Average),
                   WindSpeed_month = mean(WindSpeed_Average_m_s),
                   WindDir_month = mean(WindDir_degrees)) 
  
#avg obs for each 24 hr period and day vs. night
met <-read.csv(infile1,header=T) |> 
  dplyr::filter(as.Date(DateTime) %in% c(msn_dates)) |> 
  dplyr::select(DateTime, AirTemp_C_Average,
                  WindSpeed_Average_m_s, WindDir_degrees) |> 
  dplyr::mutate(hour = hour(DateTime)) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |> 
  dplyr::group_by(DateTime) |> 
  dplyr::summarise(AirTemp_C_24hrs = mean(AirTemp_C_Average),
                   WindSpeed_24hrs = mean(WindSpeed_Average_m_s),
                   WindDir_24hrs = mean(WindDir_degrees),
                   AirTemp_day = mean(AirTemp_C_Average[hour %in% day]),
                   AirTemp_night = mean(AirTemp_C_Average[hour %in% night]),
                   WindSpeed_day = mean(WindSpeed_Average_m_s[hour %in% day]),
                   WindSpeed_night = mean(WindSpeed_Average_m_s[hour %in% night]),
                   WindDir_day = mean(WindDir_degrees[hour %in% day]),
                   WindDir_night = mean(WindDir_degrees[hour %in% night]))
  
#make an environmental driver df
msn_drivers <- data.frame("groups" = as.character(1:5), #epi is 0.1m, hypo is 10m - avg for both noons of msn when available
     "epilimnetic_temperature" = c(ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                    ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                    ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                    ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                    ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]), 
     "hypolimnetic_temperature" = c(ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                     ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                     ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                     ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                     ctd.final_avg$Temp_C_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
     "epilimnetic_sp_cond" = c(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                      ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                      ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                      ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                      ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]),
     "hypolimnetic_sp_cond" = c(ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                       ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                       ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                       ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                       ctd.final_avg$SpCond_uScm_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
     "epilimnetic_chlorophyll" = c(ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                   ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                   ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                   ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                   ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]),
     "hypolimnetic_chlorophyll" = c(ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==1],
                    ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==2],
                    ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==3],
                    ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==4],
                    ctd.final_avg$Chla_ugL_1[ctd.final_avg$Depth_m==10 & ctd.final_avg$msn==5]),
     "epilimnetic_PAR" = c(ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==1],
                   ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==2],
                   ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==3],
                   ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==4],
                   ctd.final_avg$PAR_umolm2s_1[ctd.final_avg$Depth_m==0.1 & ctd.final_avg$msn==5]),
     "epilimnetic_TN" = c(chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==1],
                  chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==2],
                  chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==3],
                  chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==4],
                  chem_avg$TN_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==5]), 
     "epilimnetic_TP" = c(chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==1],
                  chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==2],
                  chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==3],
                  chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==4],
                  chem_avg$TP_ugL_1[chem_avg$Depth_m==0.1 & chem_avg$msn==5]), 
     "Secchi" = as.numeric(secchi$Secchi_m),
     "thermocline_depth" = c(mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==1]),
                       mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==2]),
                       mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==3]),
                       mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==4]),
                       mean(ctd_thermo_depth$therm_depth[ctd_thermo_depth$msn==5])),
     "Monthly_wind_speed" = c(monthly_met$WindSpeed_month[monthly_met$my=="0719"],
                              monthly_met$WindSpeed_month[monthly_met$my=="0719"],
                              monthly_met$WindSpeed_month[monthly_met$my=="0820"],
                              monthly_met$WindSpeed_month[monthly_met$my=="0621"],
                              monthly_met$WindSpeed_month[monthly_met$my=="0721"]),
     "Monthly_wind_dir" = c(monthly_met$WindDir_month[monthly_met$my=="0719"],
                            monthly_met$WindDir_month[monthly_met$my=="0719"],
                            monthly_met$WindDir_month[monthly_met$my=="0820"],
                            monthly_met$WindDir_month[monthly_met$my=="0621"],
                            monthly_met$WindDir_month[monthly_met$my=="0721"]),
     "Monthly_air_temp" = c(monthly_met$AirTemp_month[monthly_met$my=="0719"],
                            monthly_met$AirTemp_month[monthly_met$my=="0719"],
                            monthly_met$AirTemp_month[monthly_met$my=="0820"],
                            monthly_met$AirTemp_month[monthly_met$my=="0621"],
                            monthly_met$AirTemp_month[monthly_met$my=="0721"]),
     "24hr_wind_speed" = c(mean(met$WindSpeed_24hrs[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                           mean(met$WindSpeed_24hrs[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                           mean(met$WindSpeed_24hrs[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                           mean(met$WindSpeed_24hrs[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                           mean(met$WindSpeed_24hrs[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "24hr_wind_dir" = c(mean(met$WindDir_24hrs[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                         mean(met$WindDir_24hrs[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                         mean(met$WindDir_24hrs[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                         mean(met$WindDir_24hrs[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                         mean(met$WindDir_24hrs[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "24hr_air_temp" = c(mean(met$AirTemp_C_24hrs[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                         mean(met$AirTemp_C_24hrs[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                         mean(met$AirTemp_C_24hrs[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                         mean(met$AirTemp_C_24hrs[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                         mean(met$AirTemp_C_24hrs[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "Day_wind_speed" = c(mean(met$WindSpeed_day[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                          mean(met$WindSpeed_day[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                          mean(met$WindSpeed_day[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                          mean(met$WindSpeed_day[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                          mean(met$WindSpeed_day[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "Day_wind_dir" = c(mean(met$WindDir_day[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                        mean(met$WindDir_day[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                        mean(met$WindDir_day[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                        mean(met$WindDir_day[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                        mean(met$WindDir_day[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "Day_air_temp" = c(mean(met$AirTemp_day[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                        mean(met$AirTemp_day[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                        mean(met$AirTemp_day[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                        mean(met$AirTemp_day[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                        mean(met$AirTemp_day[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "Night_wind_speed" = c(mean(met$WindSpeed_night[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                          mean(met$WindSpeed_night[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                          mean(met$WindSpeed_night[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                          mean(met$WindSpeed_night[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                          mean(met$WindSpeed_night[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "Night_wind_dir" = c(mean(met$WindDir_night[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                        mean(met$WindDir_night[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                        mean(met$WindDir_night[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                        mean(met$WindDir_night[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                        mean(met$WindDir_night[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     "Night_air_temp" = c(mean(met$AirTemp_night[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                        mean(met$AirTemp_night[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                        mean(met$AirTemp_night[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                        mean(met$AirTemp_night[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                        mean(met$AirTemp_night[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
     
    "Cladocera_DVM" = c(migration_metrics$value[
      migration_metrics$metric=="Cladocera_density_NopL" & 
      migration_metrics$migration=="DVM_avg"]),
    "Copepoda_DVM" = c(migration_metrics$value[
      migration_metrics$metric=="Copepoda_density_NopL" & 
      migration_metrics$migration=="DVM_avg"]),
    "Rotifera_DVM" = c(migration_metrics$value[
      migration_metrics$metric=="Rotifera_density_NopL" & 
      migration_metrics$migration=="DVM_avg"]),
    "Cladocera_DHM" = c(migration_metrics$value[
      migration_metrics$metric=="Cladocera_density_NopL" & 
      migration_metrics$migration=="DHM_avg"]),
    "Copepoda_DHM" = c(migration_metrics$value[
      migration_metrics$metric=="Copepoda_density_NopL" & 
      migration_metrics$migration=="DHM_avg"]),
    "Rotifera_DHM" = c(migration_metrics$value[
      migration_metrics$metric=="Rotifera_density_NopL" & 
      migration_metrics$migration=="DHM_avg"]))

#join driver and env dfs
zoops_plus_drivers <- left_join(zoop_avg, msn_drivers)

#fit environmental drivers onto ordination
fit_env <- envfit(ord$sites, zoops_plus_drivers[,c(25:47)])

#pull out vectors - need to multiply by the sqrt of r2 to get magnitude!
scores <- data.frame((fit_env$vectors)$arrows * sqrt(fit_env$vectors$r), 
                     pvals=(fit_env$vectors)$pvals)
scores <- cbind(scores, env = rownames(scores))


#plot drivers w/ NMDS
ord <- ordiplot(NMDS_temporal_avg_bray,display = c('sites','species'),
                choices = c(1,2),type = "n")
days <- gg_ordiplot(ord, zoop_avg$groups, kind = "ehull", 
                    ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day_env <- days$plot + geom_point() + theme_bw() + geom_path() +
  geom_polygon(data = days$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=days$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B")) +
  theme(text = element_text(size=6), axis.text = element_text(size=6, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"), 
        legend.key = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(ncol=3)) +
  scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
  scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"),
                     label=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020',
                             '15-16 Jun 2021', '7-8 Jul 2021')) +
  geom_segment(data = scores,
              aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
              arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  #geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = env), size = 3)
  
  #2 vs. 1
  annotate(geom="text", x= 0.22, y= 0.63, label="hypolimnetic sp. cond.", size=1.5) +
  annotate(geom="text", x= 0.38, y= 0.56, label="epilimnetic chl.", size=1.5) +
  annotate(geom="text", x= 0.35, y= 0.51, label="daily wind dir.", size=1.5) +
  annotate(geom="text", x= 0.09, y= 0.46, label="24-hr wind dir.", size=1.5) +
  annotate(geom="text", x= 0.37, y= 0.41, label="epilimnetic sp. cond.", size=1.5) +
  annotate(geom="text", x= 0.32, y= 0.16, label="monthly air temp.", size=1.5) +
  annotate(geom="text", x= -0.07, y= 0.39, label="nightly wind dir.", size=1.5) +
  annotate(geom="text", x= -0.36, y= 0.35, label="monthly wind dir.", size=1.5) +
  annotate(geom="text", x= -0.42, y= 0.29, label="nightly wind speed", size=1.5) +
  annotate(geom="text", x= -0.26, y= 0.12, label="24-hr wind speed", size=1.5) +
  annotate(geom="text", x= -0.7, y= 0.25, label="monthly wind speed", size=1.5) +
  annotate(geom="text", x= -0.57, y= 0.08, label="Secchi", size=1.5) +
  annotate(geom="text", x= -0.69, y= 0.01, label="thermocline depth", size=1.5) +
  annotate(geom="text", x= -0.44, y= -0.05, label="daily wind speed", size=1.5) +
  annotate(geom="text", x= 0.81, y= 0.25, label="epilimnetic PAR", size=1.5) +
  annotate(geom="text", x= 0.55, y= 0.04, label="epilimnetic temp.", size=1.5) +
  annotate(geom="text", x= 0.7, y= 0, label="hypolimnetic chl.", size=1.5) +
  annotate(geom="text", x= 0.8, y= -0.08, label="hypolimnetic temp.", size=1.5) +
  annotate(geom="text", x= 0.8, y= -0.14, label="epilimnetic TN", size=1.5) +
  annotate(geom="text", x= 0.8, y= -0.21, label="epilimnetic TP", size=1.5) +
  annotate(geom="text", x= 0.48, y= -0.13, label="daily air temp.", size=1.5) +
  annotate(geom="text", x= 0.6, y= -0.24, label="24-hr air temp.", size=1.5) +
  annotate(geom="text", x= 0.52, y= -0.33, label="nightly air temp.", size=1.5)
  
  #3 vs. 1
  #annotate(geom="text", x= 0.51, y= -0.33, label="hypolimnetic sp. cond.", size=1.5) +
  #geom_segment(aes(x= 0.18, y= 0, xend= 0.4, yend= -0.31), linewidth=0.2)+
  #annotate(geom="text", x= 0.23, y= -0.4, label="epilimnetic chl.", size=1.5) +
  #geom_segment(aes(x= 0.2, y= 0.05, xend= 0.24, yend= -0.39), linewidth=0.2)+
  #annotate(geom="text", x= 0.25, y= -0.08, label="daily wind dir.", size=1.5) +
  #annotate(geom="text", x= 0.14, y= -0.5, label="24-hr wind dir.", size=1.5) +
  #geom_segment(aes(x= 0.09, y= 0.03, xend= 0.04, yend= -0.48), linewidth=0.2)+
  #annotate(geom="text", x= 0.4, y= 0.24, label="epilimnetic sp. cond.", size=1.5) +
  #geom_segment(aes(x= 0.19, y= 0.15, xend= 0.2, yend= 0.22), linewidth=0.2)+
  #annotate(geom="text", x= 0.3, y= 0.3, label="monthly air temp.", size=1.5) +
  #annotate(geom="text", x= -0.25, y= 0.28, label="nightly wind dir.", size=1.5) +
  #geom_segment(aes(x=-0.04, y=0.12, xend=-0.24, yend=0.26), linewidth=0.2)+
  #annotate(geom="text", x= -0.35, y= 0.17, label="monthly wind dir.", size=1.5) +
  #annotate(geom="text", x= -0.41, y= -0.03, label="nightly wind speed", size=1.5) +
  #annotate(geom="text", x= -0.45, y= -0.07, label="24-hr wind speed", size=1.5) +
  #annotate(geom="text", x= -0.67, y= 0, label="monthly wind speed", size=1.5) +
  #annotate(geom="text", x= -0.6, y= -0.25, label="Secchi", size=1.5) +
  #annotate(geom="text", x= -0.55, y= 0.07, label="thermocline depth", size=1.5) +
  #annotate(geom="text", x= -0.47, y= -0.1, label="daily wind speed", size=1.5) +
  #annotate(geom="text", x= 0.8, y= -0.15, label="epilimnetic PAR", size=1.5) +
  #annotate(geom="text", x= 0.4, y= 0.18, label="epilimnetic temp.", size=1.5) +
  #annotate(geom="text", x= 0.77, y= 0.04, label="hypolimnetic chl.", size=1.5) +
  #annotate(geom="text", x= 0.78, y= -0.04, label="hypolimnetic temp.", size=1.5) +
  #annotate(geom="text", x= 0.8, y= 0.15, label="epilimnetic TN", size=1.5) +
  #annotate(geom="text", x= 0.8, y= 0.08, label="epilimnetic TP", size=1.5) +
  #annotate(geom="text", x= 0.44, y= -0.15, label="daily air temp.", size=1.5) +
  #geom_segment(aes(x= 0.43, y= 0.09, xend= 0.42, yend= -0.13), linewidth=0.2) +
  #annotate(geom="text", x= 0.6, y= -0.21, label="24-hr air temp.", size=1.5) +
  #geom_segment(aes(x= 0.46, y= 0.10, xend= 0.5, yend= -0.19), linewidth=0.2) +
  #annotate(geom="text", x= 0.78, y= 0.31, label="nightly air temp.", size=1.5) +
  #geom_segment(aes(x= 0.49, y= 0.11, xend= 0.78, yend= 0.29), linewidth=0.2)

  #4 vs. 1
  #annotate(geom="text", x= 0.5, y= -0.22, label="hypolimnetic sp. cond.", size=1.5) +
  #annotate(geom="text", x= 0.43, y= -0.26, label="epilimnetic chl.", size=1.5) +
  #annotate(geom="text", x= 0.49, y= -0.33, label="daily wind dir.", size=1.5) +
  #annotate(geom="text", x= 0.27, y= -0.37, label="24-hr wind dir.", size=1.5) +
  #annotate(geom="text", x= 0.39, y= -0.14, label="epilimnetic sp. cond.", size=1.5) +
  #annotate(geom="text", x= 0.28, y= -0.05, label="monthly air temp.", size=1.5) +
  #annotate(geom="text", x= 0.08, y= -0.42, label="nightly wind dir.", size=1.5) +
  #geom_segment(aes(x= 0.08, y= -0.32, xend= 0.08, yend= -0.4), linewidth=0.2) +
  #annotate(geom="text", x= -0.4, y= -0.16, label="monthly wind dir.", size=1.5) +
  #geom_segment(aes(x= -0.11, y= -0.12, xend= -0.23, yend= -0.16), linewidth=0.2) +
  #annotate(geom="text", x= -0.29, y= -0.35, label="nightly wind speed", size=1.5) +
  #annotate(geom="text", x= -0.36, y= -0.24, label="24-hr wind speed", size=1.5) +
  #annotate(geom="text", x= -0.67, y= -0.03, label="monthly wind speed", size=1.5) +
  #annotate(geom="text", x= -0.51, y= 0.12, label="Secchi", size=1.5) +
  #annotate(geom="text", x= -0.68, y= 0.04, label="thermocline depth", size=1.5) +
  #annotate(geom="text", x= -0.42, y= -0.11, label="daily wind speed", size=1.5) +
  #annotate(geom="text", x= 0.8, y= -0.13, label="epilimnetic PAR", size=1.5) +
  #annotate(geom="text", x= 0.44, y= 0.07, label="epilimnetic temp.", size=1.5) +
  #annotate(geom="text", x= 0.68, y= 0.14, label="hypolimnetic chl.", size=1.5) +
  #geom_segment(aes(x= 0.74, y= 0.07, xend= 0.68, yend= 0.13), linewidth=0.2) +
  #annotate(geom="text", x= 0.79, y= 0.03, label="hypolimnetic temp.", size=1.5) +
  #annotate(geom="text", x= 0.66, y= -0.07, label="epilimnetic TN", size=1.5) +
  #annotate(geom="text", x= 0.84, y= 0.1, label="epilimnetic TP", size=1.5) +
  #annotate(geom="text", x= 0.59, y= 0.19, label="daily air temp.", size=1.5) +
  #annotate(geom="text", x= 0.61, y= 0.24, label="24-hr air temp.", size=1.5) +
  #annotate(geom="text", x= 0.63, y= 0.3, label="nightly air temp.", size=1.5)
#ggsave("figures/NMDS_2v1_envfit.jpg", NMDS_day_env, width=3, height=2) 

#------------------------------------------------------------------------------#
#multiplot of env variable vs NMDS1 centroids for SI

#convert msn_drivers from wide to long
msn_drivers_long <- msn_drivers %>% pivot_longer(cols = epilimnetic_temperature:Night_air_temp, 
                                                 names_to = "variable")

#add NMDS1 col
msn_drivers_long$NMDS1 <- ifelse(msn_drivers_long$groups==1, days$df_mean.ord$x[1],
                          ifelse(msn_drivers_long$groups==2, days$df_mean.ord$x[2],
                          ifelse(msn_drivers_long$groups==3, days$df_mean.ord$x[3],
                          ifelse(msn_drivers_long$groups==4, days$df_mean.ord$x[4],
                                                 days$df_mean.ord$y[5]))))

#multipanel plot (Figure S2)
driver_NMDS <- ggplot(data=msn_drivers_long, aes(NMDS1, value, color=groups)) + 
  geom_point() + facet_wrap(~variable, scales = "free_y") + scale_color_manual(
    "",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), 
       labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020",
                 "15-16 Jun 2021","7-8 Jul 2021"), guide=guide_legend(order=1)) +
                    theme(text = element_text(size=5), 
                          axis.text = element_text(size=5, color="black"), 
                          legend.background = element_blank(), 
                          legend.key = element_blank(), 
                          legend.key.height=unit(0.3,"line"),
                          axis.text.x = element_text(angle = 90, 
                                                     vjust = 0.5, hjust=1), 
                          strip.background = element_rect(fill = "transparent"), 
                          legend.position = c(0.88,0.09), 
                          legend.spacing = unit(-0.5, 'cm'), 
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), 
                          legend.key.width =unit(0.7,"line"))
#ggsave("figures/driver_vs_NMDS1.jpg", driver_NMDS, width=5, height=5) 

#-----------------------------------------------------------------------------------------------
#migration metrics across taxa vs dispersion to better link migration and variability (Fig. S10)

#add dispersion for sampling days to df
msn_drivers$disp <- c(mean(within_day_dist$dist[within_day_dist$group=="day1"]),
                      mean(within_day_dist$dist[within_day_dist$group=="day2"]),
                      mean(within_day_dist$dist[within_day_dist$group=="day3"]),
                      mean(within_day_dist$dist[within_day_dist$group=="day4"]),
                      mean(within_day_dist$dist[within_day_dist$group=="day5"]))


migration_drivers <- msn_drivers |>  
  pivot_longer(cols = Cladocera_DVM:Rotifera_DHM, names_to = "variable")

#add NMDS2 col
migration_drivers$NMDS2 <- ifelse(migration_drivers$groups==1, days$df_mean.ord$y[1],
                           ifelse(migration_drivers$groups==2, days$df_mean.ord$y[2],
                           ifelse(migration_drivers$groups==3, days$df_mean.ord$y[3],
                           ifelse(migration_drivers$groups==4, days$df_mean.ord$y[4],
                                                      days$df_mean.ord$y[5]))))

migration_vs_disp <- ggplot(data=migration_drivers, aes(disp, value, color=groups)) + 
  geom_point() + facet_wrap(~variable, scales = "free_y") + ylab("Migration metric") +
  scale_color_manual("",values=c("#008585","#9BBAA0","#F2E2B0","#DEA868","#C7522B"), 
                     labels=c("10-11 Jul 2019","24-25 Jul 2019","12-13 Aug 2020",
                              "15-16 Jun 2021","7-8 Jul 2021"), 
                     guide=guide_legend(order=1)) +
  guides(color = guide_legend(nrow=2,byrow=TRUE)) +
  theme(text = element_text(size=6), axis.text = element_text(size=5, color="black"), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-4,-4,-4,-4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", 
        legend.spacing = unit(-0.5, 'cm'), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.key.width =unit(0.7,"line")) 
#ggsave("figures/migration_metrics_vs_dispersion.jpg", migration_vs_disp, width=3, height=3) 
