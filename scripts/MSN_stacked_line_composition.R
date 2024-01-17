# stacked line plot for each MSN
# created 13 Jan 2023

#read in packages
pacman::p_load(zoo, dplR, dplyr, tidyverse, scales, 
               ggplot2, ggpubr, sf, lubridate, NatParksPalettes)

#read in zoop data from EDI
inUrl1  <-  "https://pasta-s.lternet.edu/package/data/eml/edi/1090/19/c7a04035b0a99adc489f5b6daec1cd52" 
infile1 <-  tempfile()
download.file(inUrl1,infile1,method="curl")

#list msn dates
msn_dates <- c("2019-07-10", "2019-07-11", "2019-07-24",
               "2019-07-25", "2020-08-12", "2020-08-13",
               "2021-06-15", "2021-06-16", "2021-07-07", 
               "2021-07-08")

#list of 17 "dominant" taxa (> 0.1%)
taxa <- c("Calanoida","Cyclopoida","Nauplius", "Bosmina",
          "Ceriodaphnia", "Daphnia", "Conochilus", "Keratella","Trichocerca") 
#need to remove n=8 bc they were not observed at all sampling time
# "Pompholyx", "Conochiloides", "Lepadella", "Monostyla",  "Polyarthra", "Kellicottia", "Gastropus", "Collotheca"

zoops <- read.csv(infile1, header=T) |>
  mutate(DateTime = as_datetime(DateTime),
           date = as.Date(DateTime)) |> 
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m == 4 & 
           date %in% c(msn_dates) &
           Taxon %in% c(taxa)) |> 
  select(-c(Site,EndDepth_m,CollectionMethod)) |> 
  mutate(hour = ifelse(hour(DateTime) == 11, "12:00",
                       ifelse(hour(DateTime) == 23 |
                                hour(DateTime) == 0, "00:00",
                              ifelse(hour(DateTime) == 03, "4:00", 
                                     ifelse(DateTime %in% c("2019-07-10 19:04:00"), "18:00",
                                     format(paste0(hour(DateTime),":00"), 
                                            format='%H:%M')))))) |> 
  group_by(date, hour, Taxon) |> 
  summarise(Density_IndPerL = mean(Density_IndPerL)) |>    #average reps when appropriate
  mutate(msn = ifelse(date=="2019-07-10" | 
                        date=="2019-07-11", "10-11 Jul 2019",
                      ifelse(date=="2019-07-24" | 
                               date=="2019-07-25", "24-25 Jul 2019",
                             ifelse(date=="2020-08-12" | 
                                      date=="2020-08-13","12-13 Aug 2020",
                                    ifelse(date=="2021-06-15" | 
                                             date=="2021-06-16",
                                           "15-16 Jun 2021","7-8 Jul 2021")))))
  
zoops$date_temp <- ifelse(zoops$date=="2019-07-10" | 
                            zoops$date=="2019-07-24"| 
                            zoops$date=="2020-08-12" | 
                            zoops$date=="2021-06-15" |
                            zoops$date=="2021-07-07",
                           "2024-01-13","2024-01-14")

#combine hour and date 
zoops$hour <- strptime(paste0(as.character(zoops$date_temp), 
                              zoops$hour),format="%Y-%m-%d %H:%M")
zoops$hour <- as.POSIXct(zoops$hour)

#order by hour for plotting
zoops <- zoops[order(zoops$hour),]

#look at doy on x and year by color
ggplot(zoops, aes(hour, Density_IndPerL, color=as.factor(msn))) + 
  geom_point() + theme_bw() + geom_line() +
  scale_x_datetime(expand = c(0,0),labels = date_format("%H:%M",tz="EST5EDT")) +
  facet_wrap(~forcats::fct_relevel(Taxon, "Bosmina", "Calanoida", "Conochilus",
                                   "Ceriodaphnia", "Cyclopoida", "Keratella",
                                   "Daphnia", "Nauplius", "Trichocerca"), 
             scales="free_y", nrow=3) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6), 
                     labels=c("10-11 Jul 2019","24-25 Jul 2019",
                              "12-13 Aug 2020","15-16 Jun 2021","7-8 Jul 2021"))+ 
  coord_cartesian(clip = 'off') 

#calculate relative abundance by hour for each group
zoops_relabund <- zoops |> group_by(date, hour) |> 
  mutate(total = sum(Density_IndPerL)) |> 
    ungroup() |> 
  group_by(date, hour, Taxon) |> 
  mutate(rel_abund = Density_IndPerL / total)

zoops_relabund$Taxon <- factor(zoops_relabund$Taxon, 
                         levels = c("Bosmina","Ceriodaphnia","Daphnia",
                                    "Calanoida","Cyclopoida","Nauplius",
                                    "Conochilus","Keratella","Trichocerca"))

#stacked line plot time
ggplot(zoops_relabund, aes(x = hour, y = rel_abund)) + 
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "stack", stat="identity",
            linewidth=3) +
  facet_wrap(~msn, scales = "free_x")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 9))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 9))+
  ylab("Relative abundance")+
  #labs(fill = "Taxon", color = "Taxon")+
  scale_x_datetime(expand = c(0,0),labels = date_format("%H:%M",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  xlab("Hour") +
  guides(color=guide_legend(ncol=2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.88,0.23),
        text = element_text(size=8), 
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = NatParksPalettes::natparks.pals("KingsCanyon", 1)),
        panel.spacing = unit(0.5, "lines"))
ggsave("figures/BVR_stacked_composition.jpg", width=5, height=4) 



