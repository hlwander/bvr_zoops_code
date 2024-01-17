#revisions - temporal wind and air temp aggregations for nmds driver plot

#monthly met for jul 2019, Aug 2020, Jun 2021, Jul 2021
monthly_met <- read.csv(infile1,header=T) |> 
  dplyr::mutate(my = paste0(format(as.Date(DateTime), "%m"),
                            format(as.Date(DateTime), "%y"))) |> 
  dplyr::filter(my %in% c("0719","0820","0621","0721")) |> 
  dplyr::select(my, AirTemp_C_Average,
                WindSpeed_Average_m_s, WindDir_degrees) |> 
  dplyr::group_by(my) |> 
  dplyr::summarise(AirTemp_month = mean(AirTemp_C_Average),
                   WindSpeed_month = mean(WindSpeed_Average_m_s)) 

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
                   AirTemp_day = mean(AirTemp_C_Average[hour %in% day]),
                   AirTemp_night = mean(AirTemp_C_Average[hour %in% night]),
                   WindSpeed_day = mean(WindSpeed_Average_m_s[hour %in% day]),
                   WindSpeed_night = mean(WindSpeed_Average_m_s[hour %in% night]))

#avg obs for each 24 hr period
met_24 <-read.csv(infile1,header=T) |> 
  dplyr::mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S")) |> 
  dplyr::filter((DateTime >= as.POSIXct("2019-07-10 12:00:00") & 
                   DateTime <= as.POSIXct("2019-07-11 12:00:00")) |
                  (DateTime >= as.POSIXct("2019-07-24 12:00:00") & 
                     DateTime <= as.POSIXct("2019-07-25 12:00:00")) |
                  (DateTime >= as.POSIXct("2020-08-12 12:00:00") & 
                     DateTime <= as.POSIXct("2020-08-13 12:00:00")) |
                  (DateTime >= as.POSIXct("2021-06-15 12:00:00") & 
                     DateTime <= as.POSIXct("2021-06-16 12:00:00")) |
                  (DateTime >= as.POSIXct("2021-07-07 12:00:00") & 
                     DateTime <= as.POSIXct("2021-07-08 12:00:00"))) |> 
  dplyr::select(DateTime, AirTemp_C_Average,
                WindSpeed_Average_m_s, WindDir_degrees) |> 
  dplyr::mutate(hour = hour(DateTime)) |> 
  dplyr::mutate(msn = ifelse(as.Date(DateTime) == as.Date(msn_dates[1]) | 
                               as.Date(DateTime) == as.Date(msn_dates[2]), 1,
                             ifelse(as.Date(DateTime) == as.Date(msn_dates[3]) |
                                      as.Date(DateTime) == as.Date(msn_dates[4]), 2,
                                    ifelse(as.Date(DateTime) == as.Date(msn_dates[5]) |
                                             as.Date(DateTime) == as.Date(msn_dates[6]), 3,
                                           ifelse(as.Date(DateTime) == as.Date(msn_dates[7]) |
                                                    as.Date(DateTime) == as.Date(msn_dates[8]), 4, 5))))) |> 
  dplyr::group_by(msn) |> 
  dplyr::summarise(AirTemp_C = mean(AirTemp_C_Average),
                   WindSpeed = mean(WindSpeed_Average_m_s))

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
                          "Monthly_air_temp" = c(monthly_met$AirTemp_month[monthly_met$my=="0719"],
                                                 monthly_met$AirTemp_month[monthly_met$my=="0719"],
                                                 monthly_met$AirTemp_month[monthly_met$my=="0820"],
                                                 monthly_met$AirTemp_month[monthly_met$my=="0621"],
                                                 monthly_met$AirTemp_month[monthly_met$my=="0721"]),
                          "24hr_wind_speed" = c(mean(met_24$WindSpeed[met_24$msn==1]),
                                                mean(met_24$WindSpeed[met_24$msn==2]),
                                                mean(met_24$WindSpeed[met_24$msn==3]),
                                                mean(met_24$WindSpeed[met_24$msn==4]),
                                                mean(met_24$WindSpeed[met_24$msn==5])),
                          "24hr_air_temp" = c(mean(met_24$AirTemp_C[met_24$msn==1]),
                                              mean(met_24$AirTemp_C[met_24$msn==2]),
                                              mean(met_24$AirTemp_C[met_24$msn==3]),
                                              mean(met_24$AirTemp_C[met_24$msn==4]),
                                              mean(met_24$AirTemp_C[met_24$msn==5])),
                          "Day_wind_speed" = c(mean(met$WindSpeed_day[met$DateTime %in% c("2019-07-10", "2019-07-11")]),
                                               mean(met$WindSpeed_day[met$DateTime %in% c("2019-07-24", "2019-07-25")]),
                                               mean(met$WindSpeed_day[met$DateTime %in% c("2020-08-12", "2020-08-13")]),
                                               mean(met$WindSpeed_day[met$DateTime %in% c("2021-06-15", "2021-06-16")]),
                                               mean(met$WindSpeed_day[met$DateTime %in% c("2021-07-07", "2021-07-08")])),
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
fit_env <- envfit(ord$sites, zoops_plus_drivers[,c(23:41)]) 

#pull out vectors - need to multiply by the sqrt of r2 to get magnitude!
scores <- data.frame((fit_env$vectors)$arrows * sqrt(fit_env$vectors$r), 
                     pvals=(fit_env$vectors)$pvals)
scores <- cbind(scores, env = rownames(scores))

#supplmental table w/ r2 and p values for ms
driver_correlation <- data.frame("variable" = scores$env,
                                 "R2" = fit_env$vectors$r,
                                 "p-value" = fit_env$vectors$pvals)
#write.csv(driver_correlation, "output/driver_correlation_temporal_agg.csv", row.names=F)

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
  geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = env), size = 1.5)
  
  #2 vs. 1
  #annotate(geom="text", x= 0.35, y= 0.21, label="hypolimnetic sp. cond.", size=1.5) +
  #annotate(geom="text", x= 0.37, y= 0.12, label="epilimnetic chl.", size=1.5) +
  #annotate(geom="text", x= 0.25, y= 0.03, label="epilimnetic sp. cond.", size=1.5) +
  #annotate(geom="text", x= 0, y= 0.3, label="wind speed", size=1.5) +
  #annotate(geom="text", x= -0.25, y= 0.5, label="Secchi", size=1.5) +
  #annotate(geom="text", x= -0.53, y= 0.5, label="thermocline depth", size=1.5) +
  #annotate(geom="text", x= 0.73, y= -0.32, label="epilimnetic PAR", size=1.5) +
  #annotate(geom="text", x= 0.25, y= -0.29, label="epilimnetic temp.", size=1.5) +
  #annotate(geom="text", x= 0.71, y= -0.45, label="hypolimnetic chl.", size=1.5) +
  #annotate(geom="text", x= 0.7, y= -0.6, label="hypolimnetic temp.", size=1.5) +
  #annotate(geom="text", x= 0.6, y= -0.65, label="epilimnetic TN", size=1.5) +
  #annotate(geom="text", x= 0.55, y= -0.7, label="epilimnetic TP", size=1.5) +
  #annotate(geom="text", x= 0, y= -0.39, label="air temp.", size=1.5) 
#ggsave("figures/NMDS_2v1_envfit_temporal_aggs.jpg", NMDS_day_env, width=3, height=2.5) 
