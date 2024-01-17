#first need to run lines 1-112 in 06_zoops_multivariate_stats.R

#day vs night ordinations for reviewer 1
zoop_day <- zoop_avg |> filter(order %in% c(1:4,9:11)) #6am to 8pm
zoop_night <- zoop_avg |> filter(!order %in% c(1:4,9:11)) #9pm to 5am

#only select data cols
zoop_temporal_day <- zoop_day[,c(grepl("mean",colnames(zoop_day)))] 
zoop_temporal_night <- zoop_night[,c(grepl("mean",colnames(zoop_night)))] 

#hellinger transformation
zoop_temporal_day_trans <- hellinger(zoop_temporal_day)
zoop_temporal_night_trans <- hellinger(zoop_temporal_night)

#create Euclidean distance matrix
zoop_euc_day<- as.matrix(vegdist(zoop_temporal_day_trans, method='euclidean'))
zoop_euc_night<- as.matrix(vegdist(zoop_temporal_night_trans, method='euclidean'))

#now do NMDS using averages w/ 4 dimensions for consistency
set.seed(15)
NMDS_temporal_day_bray <- metaMDS(zoop_euc_day, distance='bray', k=4, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_day_bray$stress

set.seed(15)
NMDS_temporal_night_bray <- metaMDS(zoop_euc_night, distance='bray', k=4, trymax=20, 
                                    autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_temporal_night_bray$stress

#------------------------------------------------------------------------------#
# day and night ords
ord_day <- ordiplot(NMDS_temporal_day_bray,display = c('sites','species'),
                    choices = c(1,2),type = "n")
sites_day <- gg_ordiplot(ord_day, zoop_day$site, kind = "ehull", 
                         ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_site_day <- sites_day$plot + geom_point() + theme_bw() + xlab(NULL) +
  geom_polygon(data = sites_day$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=sites_day$df_mean.ord, aes(x, y), 
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
        plot.margin = unit(c(t=1,r=-0.2,b=0,l=0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
  annotate("text", x=-0.15, y=0.5, label= "sites", 
           fontface = "italic", size = 3) +
  scale_fill_manual("",values=c("#882255","#3399CC"))+
  scale_color_manual("",values=c("#882255","#3399CC"),
                     label=c('littoral','pelagic')) 

days_day <- gg_ordiplot(ord_day, zoop_day$groups, kind = "ehull", 
                        ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day_day <- days_day$plot + geom_point() + theme_bw() + 
  geom_path() + ylab(NULL) + xlab(NULL) +
  geom_polygon(data = days_day$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=days_day$df_mean.ord, aes(x, y), 
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
        plot.margin = unit(c(t=1,r=-0.2,b=0,l=-0.2), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  annotate("text", x=-0.02, y=0.5, label= "sampling dates", 
           fontface = "italic", size = 3) +
  guides(fill="none", color = guide_legend(ncol=2)) +
  scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
  scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"),
                     label=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020',
                             '15-16 Jun 2021', '7-8 Jul 2021'))


hours_day <- gg_ordiplot(ord_day, zoop_day$order, kind = "ehull", 
                         ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 
#order hours properly
hours_day$df_hull$Group <- factor(hours_day$df_hull$Group, levels = 
                                    c(unique(hours_day$df_hull$Group)))

NMDS_hour_day <- hours_day$plot + geom_point() + theme_bw() + 
  geom_path() + ylab(NULL) + xlab(NULL) +
  geom_polygon(data = hours_day$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=hours_day$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=hcl.colors(7,"sunset")) +
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
        plot.margin = unit(c(t=1,r=0,b=0,l=-0.2), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
  annotate("text", x=0, y=0.5, label= "hours of the day",
           fontface = "italic", size=3) +
  scale_fill_manual("",values=hcl.colors(7,"sunset"))+
  scale_color_manual("",values=hcl.colors(7,"sunset"),
                     label=c('12pm','6pm','7pm','8pm','6am','7am','12pm'))



ord_night <- ordiplot(NMDS_temporal_night_bray,display = c('sites','species'),
                      choices = c(1,2),type = "n")
sites_night <- gg_ordiplot(ord_night, zoop_night$site, kind = "ehull", 
                           ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_site_night <- sites_night$plot + geom_point() + theme_bw() + xlim(c(-0.5,0.5)) +
  geom_polygon(data = sites_night$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=sites_night$df_mean.ord, aes(x, y), 
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
        legend.position = "none", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(t=0,r=-0.2,b=0,l=0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
  scale_fill_manual("",values=c("#882255","#3399CC"))+
  scale_color_manual("",values=c("#882255","#3399CC"),
                     label=c('littoral','pelagic')) 

days_night <- gg_ordiplot(ord_night, zoop_night$groups, kind = "ehull", 
                          ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

NMDS_day_night <- days_night$plot + geom_point() + theme_bw() + 
  geom_path() + ylab(NULL) + xlim(c(-0.5,0.5)) +
  geom_polygon(data = days_night$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=days_night$df_mean.ord, aes(x, y), 
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
        legend.position = "none", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(t=0,r=-0.2,b=0,l=-0.2), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + 
  guides(fill="none", color = guide_legend(ncol=2)) +
  scale_fill_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"))+
  scale_color_manual("",values=c("#008585","#89B199","#EFECBF","#DB9B5A","#C7522B"),
                     label=c('10-11 Jul 2019', '24-25 Jul 2019','12-13 Aug 2020',
                             '15-16 Jun 2021', '7-8 Jul 2021'))


hours_night <- gg_ordiplot(ord_night, zoop_night$order, kind = "ehull", 
                           ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 
#order hours properly
hours_night$df_hull$Group <- factor(hours_night$df_hull$Group, levels = 
                                      c(unique(hours_night$df_hull$Group)))

NMDS_hour_night <- hours_night$plot + geom_point() + theme_bw() + 
  geom_path() + ylab(NULL) + xlim(c(-0.5,0.5)) +
  geom_polygon(data = hours_night$df_hull, aes(x = x, y = y, fill = Group), alpha=0.2) +
  geom_point(data=hours_night$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, fill=hcl.colors(4,"sunset")) +
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
        legend.position = "none", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(t=0,r=0,b=0,l=-0.2), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill="none") +
  scale_fill_manual("",values=hcl.colors(4,"sunset"))+
  scale_color_manual("",values=hcl.colors(4,"sunset"),
                     label=c('9pm','12am','4am','5am'))


fig5 <- egg::ggarrange(NMDS_site_day, NMDS_day_day, NMDS_hour_day, 
                       NMDS_site_night, NMDS_day_night, NMDS_hour_night,
                       nrow=2,labels = c('a', '','','b','',''),
                       heights = c(1,1))
  
#ggsave("figures/DAYNIGHT_NMDS_multipanel_2v1.jpg",fig5, width=5, height=3.5) 

#-------------------------------------------------------------------------------#
#                     Calculating euclidean distance  (DAY)                          #
#-------------------------------------------------------------------------------#
#step 1: use Euclidean distance matrix from transformed community data (zoop_euc)

#order zoop epi tows by hour, MSN, and site
zoop_day <- zoop_day |> dplyr::arrange(site, groups, order)

#convert ED matrix back into distance structure for next steps
zoop_euc_day <- vegdist(zoop_temporal_day_trans, method='euclidean', 
                    upper = TRUE)

#Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
centroids_sites_day <- betadisper(zoop_euc_day, group = as.factor(zoop_day$site), 
                              type="centroid")
centroids_hours_day <- betadisper(zoop_euc_day, group = as.factor(zoop_day$order), 
                              type="centroid")
centroids_days_day <-  betadisper(zoop_euc_day, group = as.factor(zoop_day$groups), 
                              type="centroid")

#-------------------------------------------------------------------------------#
#METHOD 1:average distance of each point to polygon centroid (dispersion approach)
#METHOD 2: average distance between all combinations of centroids (pairwise approach)

#"bootstrapping" or randomly select 10 points per group and 2 groups and then calculating values using methods above

#create a df to store all the different variability values
var_results_day <- data.frame("site_disp"=rep(NA,500))

set.seed(1)

for (i in 1:500){ 
  #randomly select 10 points in each group
  ord_sub <- sample(unique(zoop_day$order), 2)
  groups_sub <-  sample(unique(zoop_day$groups), 2)
  
  zoop_sub <-  zoop_day |> group_by(order) |> 
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
  var_results_day$site_disp[i] <- mean(centroids_sites_sub$group.distances)
  var_results_day$day_disp[i] <- mean(centroids_days_sub$group.distances)
  var_results_day$hour_disp[i] <- mean(centroids_hours_sub$group.distances)
  
  var_results_day$site_pair[i] <- mean(dist(centroids_sites_sub$centroids))
  var_results_day$day_pair[i] <- mean(dist(centroids_days_sub$centroids))
  var_results_day$hour_pair[i] <- mean(dist(centroids_hours_sub$centroids))
  
}

#-------------------------------------------------------------------------------#
#Kruskal-wallis test to determine if group means are significant

#first convert wide to long
disp_df_day <- var_results_day[,grepl("disp",colnames(var_results_day))] |>  
  pivot_longer(everything(), names_to="group")
pair_df_day <- var_results_day[,grepl("pair",colnames(var_results_day))] |> 
  pivot_longer(everything(), names_to="group")

#now kw test
kw_disp_day <- kruskal.test(value ~ group, data = disp_df_day) #significant
kw_pair_day <- kruskal.test(value ~ group, data = pair_df_day) #significant

#now dunn test to determine which groups are different from each other
dunnTest(value ~ as.factor(group),
         data=disp_df_day,
         method="bonferroni")

dunnTest(value ~ as.factor(group),
         data=pair_df_day,
         method="bonferroni")

#sites, days, and hours are all different from each other!

disp_box_day <- ggboxplot(disp_df_day, x = "group", y = "value", 
                      fill = "group", palette = c("#A4C6B8", "#81858B", "#5E435D"),
                      order = c("site_disp", "day_disp", "hour_disp"),
                      ylab = "Dispersion", xlab = "") +
  theme(text = element_text(size=7),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.2,0,-0.5,0), 'lines')) +
  annotate("text",label=c("b","c","a"), x=c(1.1,2.1,3.1),
           y=c(mean(disp_df_day$value[disp_df_day$group=="site_disp"]) + 
                 sd(disp_df_day$value[disp_df_day$group=="site_disp"]),
               mean(disp_df_day$value[disp_df_day$group=="day_disp"]) + 
                 sd(disp_df_day$value[disp_df_day$group=="day_disp"]),
               mean(disp_df_day$value[disp_df_day$group=="hour_disp"]) + 
                 sd(disp_df_day$value[disp_df_day$group=="hour_disp"]))) +
  annotate("text", x=2, y=0.5, label= "Day",
           fontface = "italic", size=5) +
  guides(fill = "none") 

pair_box_day <- ggboxplot(pair_df_day, x = "group", y = "value", 
                      fill = "group", palette = c("#A4C6B8", "#81858B", "#5E435D"),
                      order = c("site_pair", "day_pair", "hour_pair"),
                      ylab = "Location", xlab = "") +
  scale_x_discrete(labels = c("sites", "sampling \n\ dates", "hours of \n\ the day")) +
  theme(text = element_text(size=7),
        plot.margin = unit(c(-0.5,0,0,0), 'lines')) +
  annotate("text",label=c("b","a","c"), x=c(1.1,2.1,3.1),
           y=c(mean(pair_df_day$value[pair_df_day$group=="site_pair"]) + 
                 sd(pair_df_day$value[pair_df_day$group=="site_pair"]),
               mean(pair_df_day$value[pair_df_day$group=="day_pair"]) + 
                 sd(pair_df_day$value[pair_df_day$group=="day_pair"]),
               mean(pair_df_day$value[pair_df_day$group=="hour_pair"]) + 
                 sd(pair_df_day$value[pair_df_day$group=="hour_pair"]))) +
  guides(fill = "none") 

among_scales_day <- egg::ggarrange(disp_box_day, pair_box_day, nrow=2)
#ggsave("figures/DAYONLY_among_variability_boxplots.jpg",among_scales_day, width=3, height=4) 


#------------------------------------------------------------------------------#
#                                  NIGHT
#------------------------------------------------------------------------------#
#step 1: use Euclidean distance matrix from transformed community data (zoop_euc)

#order zoop epi tows by hour, MSN, and site
zoop_night <- zoop_night |> dplyr::arrange(site, groups, order)

#convert ED matrix back into distance structure for next steps
zoop_euc_night <- vegdist(zoop_temporal_night_trans, method='euclidean', 
                        upper = TRUE)

#Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
centroids_sites_night <- betadisper(zoop_euc_night, group = as.factor(zoop_night$site), 
                                  type="centroid")
centroids_hours_night <- betadisper(zoop_euc_night, group = as.factor(zoop_night$order), 
                                  type="centroid")
centroids_days_night <-  betadisper(zoop_euc_night, group = as.factor(zoop_night$groups), 
                                  type="centroid")

#-------------------------------------------------------------------------------#
#METHOD 1:average distance of each point to polygon centroid (dispersion approach)
#METHOD 2: average distance between all combinations of centroids (pairwise approach)

#"bootstrapping" or randomly select 10 points per group and 2 groups and then calculating values using methods above

#create a df to store all the different variability values
var_results_night <- data.frame("site_disp"=rep(NA,500))

set.seed(1)

for (i in 1:500){ 
  #randomly select 10 points in each group
  ord_sub <- sample(unique(zoop_night$order), 2)
  groups_sub <-  sample(unique(zoop_night$groups), 2)
  
  zoop_sub <-  zoop_night |> group_by(order) |> 
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
  var_results_night$site_disp[i] <- mean(centroids_sites_sub$group.distances)
  var_results_night$day_disp[i] <- mean(centroids_days_sub$group.distances)
  var_results_night$hour_disp[i] <- mean(centroids_hours_sub$group.distances)
  
  var_results_night$site_pair[i] <- mean(dist(centroids_sites_sub$centroids))
  var_results_night$day_pair[i] <- mean(dist(centroids_days_sub$centroids))
  var_results_night$hour_pair[i] <- mean(dist(centroids_hours_sub$centroids))
  
}

#-------------------------------------------------------------------------------#
#Kruskal-wallis test to determine if group means are significant

#first convert wide to long
disp_df_night <- var_results_night[,grepl("disp",colnames(var_results_night))] |>  
  pivot_longer(everything(), names_to="group")
pair_df_night <- var_results_night[,grepl("pair",colnames(var_results_night))] |> 
  pivot_longer(everything(), names_to="group")

#now kw test
kw_disp_night <- kruskal.test(value ~ group, data = disp_df_night) #significant
kw_pair_night <- kruskal.test(value ~ group, data = pair_df_night) #significant

#now dunn test to determine which groups are different from each other
dunnTest(value ~ as.factor(group),
         data=disp_df_night,
         method="bonferroni")

dunnTest(value ~ as.factor(group),
         data=pair_df_night,
         method="bonferroni")

#sites, days, and hours are all different from each other!

disp_box_night <- ggboxplot(disp_df_night, x = "group", y = "value", 
                          fill = "group", palette = c("#A4C6B8", "#81858B", "#5E435D"),
                          order = c("site_disp", "day_disp", "hour_disp"),
                          ylab = "Dispersion", xlab = "") +
  theme(text = element_text(size=7),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.2,0,-0.5,0), 'lines')) +
  annotate("text",label=c("b","c","a"), x=c(1.1,2.1,3.1),
           y=c(mean(disp_df_night$value[disp_df_night$group=="site_disp"]) + 
                 sd(disp_df_night$value[disp_df_night$group=="site_disp"]),
               mean(disp_df_night$value[disp_df_night$group=="day_disp"]) + 
                 sd(disp_df_night$value[disp_df_night$group=="day_disp"]),
               mean(disp_df_night$value[disp_df_night$group=="hour_disp"]) + 
                 sd(disp_df_night$value[disp_df_night$group=="hour_disp"]))) +
  annotate("text", x=2, y=0.58, label= "Night",
           fontface = "italic", size=5) +
  guides(fill = "none") 

pair_box_night <- ggboxplot(pair_df_night, x = "group", y = "value", 
                          fill = "group", palette = c("#A4C6B8", "#81858B", "#5E435D"),
                          order = c("site_pair", "day_pair", "hour_pair"),
                          ylab = "Location", xlab = "") +
  scale_x_discrete(labels = c("sites", "sampling \n\ dates", "hours of \n\ the day")) +
  theme(text = element_text(size=7),
        plot.margin = unit(c(-0.5,0,0,0), 'lines')) +
  annotate("text",label=c("b","a","c"), x=c(1.1,2.1,3.1),
           y=c(mean(pair_df_night$value[pair_df_night$group=="site_pair"]) + 
                 sd(pair_df_night$value[pair_df_night$group=="site_pair"]),
               mean(pair_df_night$value[pair_df_night$group=="day_pair"]) + 
                 sd(pair_df_night$value[pair_df_night$group=="day_pair"]),
               mean(pair_df_night$value[pair_df_night$group=="hour_pair"]) + 
                 sd(pair_df_night$value[pair_df_night$group=="hour_pair"]))) +
  guides(fill = "none") 

among_scales_night <- egg::ggarrange(disp_box_night, pair_box_night, nrow=2)
#ggsave("figures/NIGHTONLY_among_variability_boxplots.jpg",among_scales_night, width=3, height=4) 





