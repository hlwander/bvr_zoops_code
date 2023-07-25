#### Zoop Summary code modified from DCR - Summer 2020 Density and Biomass in VA reservoirs

##This runs the code for zooplankton files to generate summary stats
#Summarize data and then output a file that can be merged with other lake data

##make sure to drop BVR_schind and BVR_trap samples because these are NOT tows!!!!!
#########################################

#######################START OF READING IN DATA#######################################
###Function to get the median if odd number in a column, or one entry below median if even number####
median.zoop<-function(vector){
  #Remove nas and then sort
  vector2<-vector[!is.na(vector)]
  vector2<-sort(vector2)
  #If the vector contains an odd number of entries, return the median
  if(length(vector2)%%2==1){return(median(vector2))}else{
    return(vector2[length(vector2)/2]) #If the vector is even #, return the entry immediately below the median
  }
}

###############################End of function

#Read in the csv file that contains all the taxa from the UNH zooplankton key ####
#http://cfb.unh.edu/cfbkey/html/index.html
ZoopTaxonomy<-read.csv("conversions/Zooplankton-TaxonomicIdentification.csv", fill = TRUE,as.is=TRUE)
#record which rows are NA
SpeciesNA.rows<-which(is.na(ZoopTaxonomy$Species)) 
#Replace the species with Genus and species, e.g., D. magna
ZoopTaxonomy$Species<-paste0(substring(ZoopTaxonomy$Genus,1,1),". ",ZoopTaxonomy$Species)
#If there are any NA, then fill those back in
ZoopTaxonomy$Species[SpeciesNA.rows]<-NA   

#Read in the csv file with the occular lens conversions####
ZoopOccularLensConversions<-read.csv("conversions/Zooplankton-OccularLensConversions.csv", fill = TRUE,as.is=TRUE)

#Read in csv file for size to biomass conversions for rotifers with all taxa ID'ed####
RotiferBiomassConversions<-read.csv("conversions/ZooplanktonLengthBiomassConversion-Rotifers.csv", fill = TRUE,as.is=TRUE) 
#Any of the blank species get NA 
RotiferBiomassConversions$Species[RotiferBiomassConversions$Species==""]<-NA
#Replace the species with Genus and species, e.g., K. earlinae
RotiferBiomassConversions$Species[!is.na(RotiferBiomassConversions$Species)]<-paste0(substring(RotiferBiomassConversions$Genus[!is.na(RotiferBiomassConversions$Species)],1,1),". ",RotiferBiomassConversions$Species[!is.na(RotiferBiomassConversions$Species)])

#Read in csv file for size to biomass conversions for crustaceans with all taxa ID'ed####
CrustaceanBiomassConversions<-read.csv("conversions/ZooplanktonLengthBiomassConversion-Crustaceans.csv", fill = TRUE,as.is=TRUE) 
#Any of the blank species get NA 
CrustaceanBiomassConversions$Species[CrustaceanBiomassConversions$Species==""]<-NA
#Replace the species with Genus and species, e.g., D. magna
CrustaceanBiomassConversions$Species[!is.na(CrustaceanBiomassConversions$Species)]<-paste0(substring(CrustaceanBiomassConversions$Genus[!is.na(CrustaceanBiomassConversions$Species)],1,1),". ",CrustaceanBiomassConversions$Species[!is.na(CrustaceanBiomassConversions$Species)])
#Move all Cladacora ids to suborder
CrustaceanBiomassConversions$Suborder[CrustaceanBiomassConversions$Order=="Cladocera"]<-"Cladocera"
CrustaceanBiomassConversions$Order[CrustaceanBiomassConversions$Order=="Cladocera"]<-""
#Make sure all Daphniidae are subOrder cladocera
CrustaceanBiomassConversions$Suborder[CrustaceanBiomassConversions$Family=="Daphniidae"]<-"Cladocera"

#For loop through all the years(or specify year here)
year<-2020

#Names of the two files for a particular year
dataFiles<-c(paste("FCR",year,"_ZooplanktonCounting_Density_DataEntry.csv",sep=""),paste("FCR",year,"_ZooplanktonCounting_SizeID_DataEntry.csv",sep=""))

#Create Density and SizeID data frames for that year
Density<-read.csv(paste0(getwd(),"/zoop_data/",dataFiles[1]), fill = TRUE,as.is=TRUE)
SizeID<-read.csv(paste0(getwd(),"/zoop_data/",dataFiles[2]), fill = TRUE,as.is=TRUE)

#If any SizeID are missing, add NA
SizeID$LowestTaxonomicLevelOfID[SizeID$LowestTaxonomicLevelOfID==""]<-NA
SizeID$TaxaID[SizeID$TaxaID==""]<-NA
SizeID$MarksInOcularMicrometer_No.[SizeID$MarksInOcularMicrometer_No.==""]<-NA
SizeID$ObjectiveMagnification[SizeID$ObjectiveMagnification==""]<-NA

#Radius of the sampling net in meters --> FOR 80um NET!
RadiusOfZoopNet_m<-0.1524
AreaOfZoopNet_m2<-(RadiusOfZoopNet_m^2) *pi

#Radius of the sampling net in meters --> FOR 20um NET!
RadiusOfZoopNet_m_20<-0.099
AreaOfZoopNet_m2_20<-(RadiusOfZoopNet_m_20^2) *pi

#Calculate the volume counted for each sample using each subsample volume
#conditional for each of the mesh types/methods; schindler vol= proportional vol * 30L schindler trap; trap volume = proportional vol * .25L
Density$Volume_L <- ifelse(Density$site_no=="BVR_trap", Density$Volume_L<- (Density$SubsampleVolume_mL/Density$InitialSampleVolume_mL) * .25,
                           ifelse(Density$mesh_size_μm==61, Density$Volume_L<- (Density$SubsampleVolume_mL/Density$InitialSampleVolume_mL) * 30,
                                  ifelse(Density$mesh_size_μm==80, 
###################   This is the 1-5mL in the counting chamber relative to the ######## This is the volume of the zoop tow    #####  This is conversion 
##################   concentrated/diluted total volume of the sample (200-1000mL) ######## depth of tow times area of circle (in m3) #####   from m3 to L      
Density$Volume_L <- (Density$SubsampleVolume_mL/Density$InitialSampleVolume_mL)    *      (Density$DepthOfTow_m*pi*RadiusOfZoopNet_m^2)   *    1000,                
ifelse(Density$mesh_size_μm==20, Density$Volume_L<- (Density$SubsampleVolume_mL/Density$InitialSampleVolume_mL) * (Density$DepthOfTow_m*pi*RadiusOfZoopNet_m_20^2)   *    1000,  Density$Volume_L<-NA))))

#ignoring volume counted for hypo calcs, as this will be accounted for in numerator (proportional_vol)
Density$Volume_unadj <- ifelse(Density$site_no=="BVR_trap", Density$Volume_unadj <-.25,
                               ifelse(Density$mesh_size_μm==80, Density$Volume_unadj <- (Density$DepthOfTow_m*pi*RadiusOfZoopNet_m^2)*1000, 
                               ifelse(Density$mesh_size_μm==20, Density$Volume_unadj<- (Density$DepthOfTow_m*pi*RadiusOfZoopNet_m_20^2)*1000,
                                      ifelse(Density$mesh_size_μm==61, Density$Volume_unadj<- 30, Density$Volume_unadj<-NA))))

#proportion of sample counted (used to scale raw density/biomass for hypo calcs in numerator)
Density$proportional_vol<- Density$SubsampleVolume_mL/Density$InitialSampleVolume_mL

#REPLACE plomida with ploima (they are synonymous orders!)
SizeID$TaxaID[SizeID$TaxaID=="Plomida"]<- "Ploima"

#Aggregate by sample no
#Create a vector of the sample ID
uniqueSampleID<-unique(Density$sample_ID)

#Initialize a data frame with the correct columns
ZoopDensityCalcs<-data.frame(matrix(ncol = length(names(Density)), nrow = length(uniqueSampleID)))
names(ZoopDensityCalcs)<-names(Density)
#Store the sample_no
ZoopDensityCalcs$sample_ID<-uniqueSampleID[order(uniqueSampleID)]
#Aggregate volume and # of zoop by sum (mean for vol_unadj because not dependent on sample volume)
ZoopDensityCalcs$Volume_L<-aggregate(Volume_L~sample_ID,data=Density,FUN=sum,na.action=na.omit)[,2]
ZoopDensityCalcs$Volume_unadj<-aggregate(Volume_unadj~sample_ID,data=Density,FUN=mean,na.action=na.omit)[,2]
ZoopDensityCalcs$proportional_vol<-aggregate(proportional_vol~sample_ID,data=Density,FUN=sum,na.action=na.omit)[,2]
ZoopDensityCalcs$Zooplankton_No.<-aggregate(Zooplankton_No.~sample_ID,data=Density,FUN=sum,na.action=na.omit)[,2]

#Pull the unique entries for some of the string columns
ZoopDensityCalcs$Project<-aggregate(Project~sample_ID,data=Density,FUN=unique,na.action=na.omit)[,2]
ZoopDensityCalcs$site_no<-aggregate(site_no~sample_ID,data=Density,FUN=unique,na.action=na.omit)[,2]
ZoopDensityCalcs$collect_date<-as.Date(aggregate(collect_date~sample_ID,data=Density,FUN=unique,na.action=na.omit)[,2],"%d-%b-%y")
ZoopDensityCalcs$Hour<-aggregate(Hour~sample_ID,data=Density,FUN=unique,na.action=na.omit)[,2]

ZoopDensityCalcs$INT<-ifelse(sum(!is.na(Density$INT))==0,NA,aggregate(INT~sample_ID,data=Density,FUN=unique,na.action=na.omit)[,2])

ZoopDensityCalcs$InitialSampleVolume_mL<-aggregate(InitialSampleVolume_mL~sample_ID,data=Density,FUN=unique,na.action=na.omit)[,2]
#convert to integer so can export at the end
ZoopDensityCalcs$InitialSampleVolume_mL<- as.numeric(as.character(ZoopDensityCalcs$InitialSampleVolume_mL))

#Average depth of tows and mesh size
ZoopDensityCalcs$DepthOfTow_m<-aggregate(DepthOfTow_m~sample_ID,data=Density,FUN=mean,na.action=na.pass)[,2]
ZoopDensityCalcs$mesh_size_μm<-aggregate(mesh_size_μm~sample_ID,data=Density,FUN=mean,na.action=na.pass)[,2]


#Keep important columns
keep<-c("Project","site_no","collect_date","Hour","sample_ID","DepthOfTow_m","InitialSampleVolume_mL","Zooplankton_No.","INT","Volume_L","Volume_unadj","proportional_vol","mesh_size_μm")
ZoopDensityCalcs<-ZoopDensityCalcs[,keep]

#Net efficiencies from NetEfficiencyCalcs script (note that bvr value is averaged from 2020 and 2021 data)
NetEfficiency2020 <- c(0.03664618, 0.05279169) #fcr is only for 2020 for now

#Calculate the zooplankton density 2 different ways (tows vs. schindler/horiz traps)
#multiplying # by net efficiency ratio (calculated from tow (apparent dens) / schindler (actual dens) counts in 2020 (n=2)) 
ZoopDensityCalcs$ZoopDensity_No.pL<- ifelse(ZoopDensityCalcs$site_no=="BVR_trap" | ZoopDensityCalcs$site_no=="BVR_schind" | ZoopDensityCalcs$site_no=="FCR_schind", ZoopDensityCalcs$ZoopDensity_No.pL <- ZoopDensityCalcs$Zooplankton_No./ZoopDensityCalcs$Volume_unadj,
                                            ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="BVR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopDensityCalcs$ZoopDensity_No.pL <- ZoopDensityCalcs$Zooplankton_No.* (1/NetEfficiency2020[1]) / ZoopDensityCalcs$Volume_L,
                                                   ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="FCR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopDensityCalcs$ZoopDensity_No.pL <- ZoopDensityCalcs$Zooplankton_No.* (1/NetEfficiency2020[2]) / ZoopDensityCalcs$Volume_L, 
                                                          ZoopDensityCalcs$ZoopDensity_No.pL <- ZoopDensityCalcs$Zooplankton_No. / ZoopDensityCalcs$Volume_L)))
                                                         

#convert from character to numeric
SizeID$MarksInOcularMicrometer_No.<- as.numeric(SizeID$MarksInOcularMicrometer_No.)

#Calculate the size based on the marks in the ocular micrometer and conversion
#Checks for non-NA and non-numeric strings
#Prints out the row to examine
for(i in 1:length(SizeID$MarksInOcularMicrometer_No.)){
  if(!is.na(SizeID$MarksInOcularMicrometer_No.[i])&is.na(as.numeric(SizeID$MarksInOcularMicrometer_No.[i]))){print(paste("Row",i))}
}

#Calculates the size of each zooplankton
SizeID$Size_mm<-as.numeric(SizeID$MarksInOcularMicrometer_No.)*(ZoopOccularLensConversions[match(SizeID$ObjectiveMagnification,ZoopOccularLensConversions$Objective),2]/ZoopOccularLensConversions[match(SizeID$ObjectiveMagnification,ZoopOccularLensConversions$Objective),3])

#Create a new column for each level of taxonomic richness
SizeID[,names(ZoopTaxonomy)]<-NA

#Check if order name has a first letter capital
substr(SizeID$LowestTaxonomicLevelOfID, 1, 1) <- toupper(substr(SizeID$LowestTaxonomicLevelOfID, 1, 1))

#Figure out which column the taxa ID is
for(j in 1:length(SizeID$TaxaID)){
    #Pick out only values with Taxa Level and ID
    if(!is.na(SizeID$LowestTaxonomicLevelOfID[j])&!is.na(SizeID$TaxaID[j])){
    #Get the column that the lowest level of ID
    TaxaColumn<-match(SizeID$LowestTaxonomicLevelOfID[j],names(ZoopTaxonomy))
    #Get the row for the correct taxa (it will output the first match)
    TaxaRow<-match(SizeID$TaxaID[j],ZoopTaxonomy[,TaxaColumn])
    #Run through each of the levels up to the lowest level and store it in new columns in the SizeID data frame
  #if nauplius, then record lowest taxanomic level as Copepoda (insert NA in order-species columns)
    #Note: actually might want to keep the cal vs cyc naup difference - I feel pretty confident in IDs this year (plus I actually do see a difference between them)
  for(k in 1:TaxaColumn){
      SizeID[j,names(ZoopTaxonomy)[k]]<-ZoopTaxonomy[TaxaRow,k]
  }
    }
      for(rot.i in 5:10){
      if((SizeID$TaxaID[j]=="Cyclopoida" | SizeID$TaxaID[j]=="Calanoida") & SizeID$Nauplius[j]!=""){
        SizeID[j,names(ZoopTaxonomy)[rot.i]]<- NA}
     }
  }
    
###################################################
#     Biomass conversion in SizeID data frame     #
###################################################

#For crustaceans, prepare the conversion equation
CrustaceanBiomassConversions$a<-exp(CrustaceanBiomassConversions$ln.a..intercept.)
#Replace all blanks with NA
CrustaceanBiomassConversions$Subclass[CrustaceanBiomassConversions$Subclass==""]<-NA
CrustaceanBiomassConversions$Order[CrustaceanBiomassConversions$Order==""]<-NA
CrustaceanBiomassConversions$Suborder[CrustaceanBiomassConversions$Suborder==""]<-NA
CrustaceanBiomassConversions$Family[CrustaceanBiomassConversions$Family==""]<-NA
CrustaceanBiomassConversions$Subfamily[CrustaceanBiomassConversions$Subfamily==""]<-NA
CrustaceanBiomassConversions$Genus[CrustaceanBiomassConversions$Genus==""]<-NA
CrustaceanBiomassConversions$Species[CrustaceanBiomassConversions$Species==""]<-NA

#Loop through each individual
for(k in 1:length(SizeID$TaxaID)){
  #If to find if the individual is a rotifer
  if(!is.na(SizeID[k,"Phylum"])&SizeID[k,"Phylum"]=="Rotifera"){
    #loop through each taxa starting with species and ending with Phylum
    for(taxa in length(names(ZoopTaxonomy)):1){
      #If there is an id there and that matches an id in the conversion dataframe, 
      #then do the conversion using the median value 
      if(!is.na(SizeID[k,names(ZoopTaxonomy)[taxa]])&SizeID[k,names(ZoopTaxonomy)[taxa]]  %in% as.matrix(RotiferBiomassConversions[names(ZoopTaxonomy)[taxa]])){
        #Find the median of all values of that taxa
        conversion.parameter.a<-median(RotiferBiomassConversions$Parameter_vol.a3[RotiferBiomassConversions[,names(ZoopTaxonomy)[taxa]]==SizeID[k,names(ZoopTaxonomy)[taxa]]],na.rm=T)
        #Convert to biovolume using the Volume (mm3) = parameter * length^3
        Volume_mm3<-conversion.parameter.a * (SizeID[k,"Size_mm"])^3
        #Convert from volume to fresh weight using a specific gravity of 1 kg/m3 (1000 ug/mm3)
        Freshweight_ug<-Volume_mm3*1000
        #Convert to dry weight assuming ratio of 0.1 (dry:wet) - dry weight = wet * 0.1
        SizeID[k,"Biomass_ug"]<- Freshweight_ug*0.1
        #break the for loop
        break
        #End of if loop to find which id matches
      }
      #End of looping through taxonomy
    }
    
    #End of Rotifer if else
  }else if(!is.na(SizeID[k,"Subphylum"])&SizeID[k,"Subphylum"]=="Crustacea"){
    #If it is a copepod nauplius
    if(SizeID$Nauplius[k]=="nauplius"&SizeID$Subclass[k]=="Copepoda"){  
      #Find the conversion for Copepod nauplius
      conversion.parameter.b<-CrustaceanBiomassConversions$b[is.na(CrustaceanBiomassConversions$Genus)]
      #choose the a that matches that b
      conversion.parameter.a<-CrustaceanBiomassConversions$a[is.na(CrustaceanBiomassConversions$Genus)]
      #Convert to dry mass using the equation w=e^(a+b(ln(L)) where L is length in mm
      SizeID[k,"Biomass_ug"]<-conversion.parameter.a*(SizeID[k,"Size_mm"]^conversion.parameter.b)
    }else{
      #Else for non-nauplii  
      #loop through each taxa starting with species and ending with Phylum
      for(taxa in length(names(ZoopTaxonomy)):1){
        #If there is an id there and that matches an id in the conversion dataframe, 
        #then do the conversion using the median value 
        if(!is.na(SizeID[k,names(ZoopTaxonomy)[taxa]]) & SizeID[k,names(ZoopTaxonomy)[taxa]]%in% as.matrix(CrustaceanBiomassConversions[names(ZoopTaxonomy)[taxa]])){
          #Remove NAs from the data frame
          temp.CrustaceanBiomassConversions<-CrustaceanBiomassConversions[!is.na(CrustaceanBiomassConversions[,names(ZoopTaxonomy)[taxa]]),]
          #Find the median of all values of that taxa for b
          conversion.parameter.b<-median.zoop(temp.CrustaceanBiomassConversions$b[temp.CrustaceanBiomassConversions[,names(ZoopTaxonomy)[taxa]]==SizeID[k,names(ZoopTaxonomy)[taxa]]])
          #choose the a that matches that b, use median.zoop in case there are several of that taxa with the same slope
          conversion.parameter.a<-median.zoop(temp.CrustaceanBiomassConversions$a[temp.CrustaceanBiomassConversions$b==conversion.parameter.b&temp.CrustaceanBiomassConversions[,names(ZoopTaxonomy)[taxa]]==SizeID[k,names(ZoopTaxonomy)[taxa]]])
          #Convert to dry mass using the equation w=e^(a+b(ln(L)) where L is length in mm
          SizeID[k,"Biomass_ug"]<-conversion.parameter.a*(SizeID[k,"Size_mm"]^conversion.parameter.b)
          
          #break the for loop
          break
          #End of if loop to find which id matches
        }
        #End of looping through taxonomy
      }
      
      #End of else for non-nauplii
    }
    
    
    #End of crustacea if else  
  }else{
    #All others get NA biomass
    SizeID$Biomass_ug[k]<-NA} 
  #End of loop through the individuals
}

#asplanchna WW:DW is 0.039 according to 2016 EPA SOP so multiplying this value by .39 (0.039 = 0.1 * 0.39) for all asplanchna biomass values
SizeID$Biomass_ug[SizeID$Genus=="Asplanchna"& !is.na(SizeID$Genus)] <-  SizeID$Biomass_ug[SizeID$Genus=="Asplanchna"& !is.na(SizeID$Genus)]*0.39

####################################################
####THIS SECTION IS USED TO CONFIRM THE BIOMASS CONVERSION IS WORKING####

hist(SizeID$Biomass_ug[SizeID$Subphylum=="Crustacea"])
max(SizeID$Biomass_ug[SizeID$Phylum=="Rotifera"],na.rm=T)
head(SizeID[order(SizeID$Biomass_ug,decreasing=TRUE),c("LowestTaxonomicLevelOfID","TaxaID","Size_mm","Biomass_ug")],15)
temp<-SizeID[order(SizeID$Biomass_ug,decreasing=TRUE),]
head(temp[temp$Phylum=="Rotifera"&!is.na(temp$Phylum),c("LowestTaxonomicLevelOfID","TaxaID","Size_mm","Biomass_ug")],15)
head(temp[temp$Subphylum=="Crustacea"&!is.na(temp$Phylum),c("LowestTaxonomicLevelOfID","TaxaID","Size_mm","Biomass_ug")],15)
median(temp[temp$Subphylum=="Crustacea"&!is.na(temp$Phylum),c("Biomass_ug")],na.rm=T)
length(temp[temp$Subphylum=="Crustacea"&!is.na(temp$Phylum),c("Biomass_ug")])
length(temp[temp$Phylum=="Rotifera"&!is.na(temp$Phylum),c("Biomass_ug")])
hist(SizeID$Biomass_ug[SizeID$Phylum=="Rotifera"])
SizeID[!is.na(SizeID$Biomass_ug>100),c("LowestTaxonomicLevelOfID","TaxaID","Biomass_ug")]
names(ZoopTaxonomy)

#Check total lengths to see how many in each taxa of rotifer vs. crustacean vs. neither
length(SizeID$Biomass_ug)
sum(SizeID$Phylum=="Rotifera"&!is.na(SizeID$Phylum))
sum(SizeID$Subphylum=="Crustacea"&!is.na(SizeID$Subphylum))
sum(!is.na(SizeID$Biomass_ug))
#SizeID[which(is.na(SizeID$MarksInOcularMicrometer_No.)& is.na(SizeID$Phylum)),]
#SizeID[which(!is.na(SizeID$MarksInOcularMicrometer_No.)& is.na(SizeID$Biomass_ug)),]
##########################################################

#Create new dataframe that is incorporating all the data from density calc and size
ZoopFinal<-ZoopDensityCalcs

#Overall size stats for each sample (moved biomass to after taxa calcs)
#Average size by site_ID
ZoopFinal<- merge(ZoopFinal,setNames(aggregate(SizeID$Size_mm~SizeID$sample_ID,FUN=mean,na.rm=TRUE),c("sample_ID","OverallMeanSize_mm")),by=c("sample_ID"),all.x=T)
#SD of size by site_ID
ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Size_mm~sample_ID,data=SizeID,FUN=sd,na.action=na.omit),c("sample_ID","OverallSESize_mm")),by=c("sample_ID"),all.x=T)
#count of size by site_ID
ZoopFinal<-merge(ZoopFinal,setNames(aggregate(TaxaID~sample_ID,data=SizeID,FUN=length,na.action=na.pass),c("sample_ID","OverallCount_n")),by=c("sample_ID"),all.x=T)

#total biomass is calculated after taxa biomass (see below)
ZoopFinal$TotalBiomass_ug <- NA
ZoopFinal$BiomassConcentration_ugpL <- NA
###original calc was biom_unadj * count/zoops_measured, but this means that smaller zoop biomass is largely overestimated because I generally measure a lot more of the large zoops (i.e., total biomass/zoops measured is skewed towards the laeger zoops). 

#Create columns broken down by the following taxa
taxaOfInterest<-c("Daphniidae","Copepoda","Calanoida","Cladocera","Cyclopoida","Rotifera","Keratella","Kellicottia","Crustacea","Bosminidae","nauplius","Ceriodaphnia","Daphnia","Bosmina","Ploima","Gastropidae","Collothecidae","Conochilidae","Synchaetidae","Trichocercidae","Lepadella","Monostyla","Lecane") #no holopedium this year
CorrespondingTaxaLevel<-c("Family","Subclass","Order","Suborder","Order","Phylum","Genus","Genus","Subphylum","Family","Nauplius","Genus","Genus","Genus","Order","Family","Family","Family","Family","Family","Genus","Genus","Genus")
#Here crustacean is the sum of copepoda and cladocera; 

#For loop that runs through all the taxa of interest and generates output (including nauplius!)
for(taxa.index in 1:length(taxaOfInterest)){
  taxa<-taxaOfInterest[taxa.index]
  temp<-SizeID[SizeID[,CorrespondingTaxaLevel[taxa.index]]==taxa &!is.na(SizeID[,CorrespondingTaxaLevel[taxa.index]]),]
 
  #If there are no individuals in that taxa, then export as NA
  if(nrow(temp)==0){
    ZoopFinal[,c(paste(taxa,"MeanSize_mm",sep=""),paste(taxa,"SESize_mm",sep=""))]<-NA
    ZoopFinal[,c(paste(taxa,"Count_n",sep=""),paste(taxa,"_PercentOfTotal",sep=""),paste(taxa,"_totalbiomass_ug",sep=""),paste(taxa,"_PercentOftotalbiomassConcentration",sep=""))]<-0
  }else{
    #Average size by sample_ID
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Size_mm~sample_ID,data=temp,FUN=mean,na.action=na.omit),c("sample_ID",paste(taxa,"MeanSize_mm",sep=""))),by=c("sample_ID"),all.x=T)
    #SD of size by sample_ID
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Size_mm~sample_ID,data=temp,FUN=sd,na.action=na.omit),c("sample_ID",paste(taxa,"SESize_mm",sep=""))),by=c("sample_ID"),all.x=T)
    #Count by sample_ID
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(TaxaID~sample_ID,data=temp,FUN=length),c("sample_ID",paste(taxa,"Count_n",sep=""))),by=c("sample_ID"),all.x=T)
    #Percent that are that taxa
    ZoopFinal[,paste(taxa,"_PercentOfTotal",sep="")]<-ZoopFinal[,paste(taxa,"Count_n",sep="")]*100/ZoopFinal$OverallCount_n
    #biomass - unadjusted biomass because does not account for individuals w/o size measurements ( total biomass * count/num.measured)
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Biomass_ug~sample_ID,data=temp,FUN=sum,na.action=na.omit),c("sample_ID",paste(taxa,"_biomass_unadj",sep=""))),by=c("sample_ID"),all.x=T)
    #this is number of individuals (within each taxa) measured per sample
    ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Size_mm~sample_ID,data=temp,FUN=length),c("sample_ID",paste(taxa,"_zoops_measured",sep=""))),by=c("sample_ID"),all.x=T)
    #now calculate total biomass by multiplying biom_unadj by count/zoops_measured
    ZoopFinal[,paste(taxa,"_totalbiomass_ug",sep="")]<- ZoopFinal[,paste(taxa,"_biomass_unadj",sep="")] * (ZoopFinal[,paste(taxa,"Count_n",sep="")] / ZoopFinal[,paste(taxa,"_zoops_measured",sep="")])
    #Taxa density
    ZoopFinal[,paste(taxa,"_density_NopL",sep="")]<-ifelse(ZoopFinal$site_no=="BVR_trap" | ZoopFinal$site_no=="BVR_schind" | ZoopFinal$site_no=="FCR_schind", ZoopFinal[,paste(taxa,"_density_NopL",sep="")] <- ZoopFinal[,paste(taxa,"Count_n",sep="")]/ZoopFinal$Volume_unadj,
                                                       ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="BVR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopFinal[,paste(taxa,"_density_NopL",sep="")] <- ZoopFinal[,paste(taxa,"Count_n",sep="")]* (1/NetEfficiency2020[1])/ ZoopFinal$Volume_L,
                                                              ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="FCR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopFinal[,paste(taxa,"_density_NopL",sep="")] <- ZoopFinal[,paste(taxa,"Count_n",sep="")]* (1/NetEfficiency2020[2])/ ZoopFinal$Volume_L, 
                                                              ZoopFinal[,paste(taxa,"_density_NopL",sep="")] <- ZoopFinal[,paste(taxa,"Count_n",sep="")] / ZoopFinal$Volume_L)))
    #Taxa biomass concentration
    ZoopFinal[,paste(taxa,"_BiomassConcentration_ugpL",sep="")]<-ifelse(ZoopFinal$site_no=="BVR_trap" | ZoopFinal$site_no=="BVR_schind" | ZoopFinal$site_no=="FCR_schind", ZoopFinal[,paste(taxa,"_BiomassConcentration_ugpL",sep="")]<- ZoopFinal[,paste(taxa,"_totalbiomass_ug",sep="")] / ZoopFinal$Volume_unadj,
                                                                ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="BVR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopFinal[,paste(taxa,"_BiomassConcentration_ugpL",sep="")]<- ZoopFinal[,paste(taxa,"_totalbiomass_ug",sep="")] * (1/NetEfficiency2020[1])/ ZoopFinal$Volume_L,
                                                                       ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="FCR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopFinal[,paste(taxa,"_BiomassConcentration_ugpL",sep="")]<- ZoopFinal[,paste(taxa,"_totalbiomass_ug",sep="")] * (1/NetEfficiency2020[2])/ ZoopFinal$Volume_L, 
                                                                              ZoopFinal[,paste(taxa,"_BiomassConcentration_ugpL",sep="")]<- ZoopFinal[,paste(taxa,"_totalbiomass_ug",sep="")] / ZoopFinal$Volume_L)))
    #Percent biomass concentration of total
    ZoopFinal[,paste(taxa,"_PercentOftotalbiomassConcentration",sep="")]<-ZoopFinal[,paste(taxa,"_BiomassConcentration_ugpL",sep="")]*100/ZoopFinal$BiomassConcentration_ugpL
  }   
  #End of loop through taxa of interest  
}

#now drop biomass unadj and zoop_measured cols
ZoopFinal<- ZoopFinal[, -grep("biomass_unadj", colnames(ZoopFinal))]
ZoopFinal<- ZoopFinal[, -grep("zoops_measured", colnames(ZoopFinal))]

#taxa used to calculate total biomass and biomass concentration
TaxaforTotalCalcs<- c("Calanoida","Cyclopoida","nauplius","Cladocera","Rotifera") #not sure if this should be copepod instead of cyc and cals bc not all were IDed to order (some are subclass cope)
#Calculate raw biomass 
ZoopFinal$TotalBiomass_ug<- rowSums(ZoopFinal[,paste0(TaxaforTotalCalcs,"_totalbiomass_ug")], na.rm=TRUE)
#Calculate biomass concentration as biomass relative to total volume counted
ZoopFinal$BiomassConcentration_ugpL<-ifelse(ZoopFinal$site_no=="BVR_trap" | ZoopFinal$site_no=="BVR_schind" | ZoopFinal$site_no=="FCR_schind", ZoopFinal$BiomassConcentration_ugpL <- ZoopFinal$TotalBiomass_ug / ZoopFinal$Volume_unadj,
  ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="BVR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopFinal$BiomassConcentration_ugpL <- ZoopFinal$TotalBiomass_ug *(1/NetEfficiency2020[1])/ ZoopFinal$Volume_L,
         ifelse(substr(ZoopDensityCalcs$site_no,1,3)=="FCR" & ZoopDensityCalcs$DepthOfTow_m > 8, ZoopFinal$BiomassConcentration_ugpL <- ZoopFinal$TotalBiomass_ug *(1/NetEfficiency2020[2])/ ZoopFinal$Volume_L, 
                ZoopFinal$BiomassConcentration_ugpL <- ZoopFinal$TotalBiomass_ug / ZoopFinal$Volume_L)))
  
#Replace NA with 0's for the relevant cells (count, density, biomass, biomass concentration) NOTE - if this doesn't work, then at least one taxa is not found across all samples and should be removed!
ZoopFinal[c(paste(taxaOfInterest,"Count_n",sep=""),paste(taxaOfInterest,"_PercentOfTotal",sep=""),paste(taxaOfInterest,"_totalbiomass_ug",sep=""),paste(taxaOfInterest,"_density_NopL",sep=""),paste(taxaOfInterest,"_BiomassConcentration_ugpL",sep=""),paste(taxaOfInterest,"_PercentOftotalbiomassConcentration",sep=""))][is.na(ZoopFinal[c(paste(taxaOfInterest,"Count_n",sep=""),paste(taxaOfInterest,"_PercentOfTotal",sep=""),paste(taxaOfInterest,"_totalbiomass_ug",sep=""),paste(taxaOfInterest,"_density_NopL",sep=""),paste(taxaOfInterest,"_BiomassConcentration_ugpL",sep=""),paste(taxaOfInterest,"_PercentOftotalbiomassConcentration",sep=""))])] <- 0

#Convert the initials to characters for export
ZoopFinal$INT<-as.character(ZoopFinal$INT)

#Export the final Zoop table
write.csv(ZoopFinal,paste("output/FCR_ZooplanktonSummary",year,".csv",sep=""),row.names = FALSE)

ZoopFinal$BiomassConcentration_ugpL
ZoopFinal$Cladocera_BiomassConcentration_ugpL
ZoopFinal$Calanoida_BiomassConcentration_ugpL
ZoopFinal$Cladocera_PercentOftotalbiomassConcentration
ZoopFinal$Calanoida_PercentOftotalbiomassConcentration

#####Diversity#####
##IGNORE THIS FOR EXPORT!!

#Overall diversity by different levels here
#Create a function that can aggregate by levels, genus, order/suborder
#e.g., aggregate(order_no ~ name, myvec, function(x) length(unique(x)))
#For Genus, calculate richness, ShannonDiversity, and Evenness
taxa<-"Genus"
if(sum(!is.na(SizeID[,c(taxa)]))>0){
  #Richness
  ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Genus~sample_ID,data=SizeID,FUN=function(x) length(unique(x))),c("sample_ID",paste(taxa,"Richness",sep=""))),by=c("sample_ID"),all.x=T)
  #ShannonDiversity
  ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Genus~sample_ID,data=SizeID,FUN=function(x) -sum(prop.table(table(x[!is.na(x)]))*log(prop.table(table(x[!is.na(x)]))))),c("sample_ID",paste(taxa,"ShannonDiversity",sep=""))),by=c("sample_ID"),all.x=T)
  #Evenness
  ZoopFinal<-merge(ZoopFinal,setNames(aggregate(Genus~sample_ID,data=SizeID,FUN=function(x) -sum(prop.table(table(x[!is.na(x)]))*log(prop.table(table(x[!is.na(x)]))))/log(length(unique(x)))),c("sample_ID",paste(taxa,"Evenness",sep=""))),by=c("sample_ID"),all.x=T)
  
}else{
  ZoopFinal[,c(paste(taxa,"Richness",sep=""),paste(taxa,"ShannonDiversity",sep=""),paste(taxa,"Evenness",sep=""))]<-NA
}

#For Order/Suborder, calculate richness, ShannonDiversity, and Evenness
taxa<-"OrderSubOrder"
#Create dataframe that is the order and suborder with sample_no
#Merge the columns and remove the NAs, only the non-NA column is selected OR NA if both are NA
SizeID.Order<-SizeID[,c("sample_ID","Order","Suborder")]
SizeID.Order$OrderSubOrder<-apply(SizeID.Order[,c("Order","Suborder")], 1, function(x) paste(na.omit(x),collapse=""))
SizeID.Order$OrderSubOrder[SizeID.Order$OrderSubOrder==""]<-NA
if(sum(!is.na(SizeID.Order[,c(taxa)]))>0){
  #Richness
  ZoopFinal<-merge(ZoopFinal,setNames(aggregate(OrderSubOrder~sample_ID,data=SizeID.Order,FUN=function(x) length(unique(x))),c("sample_ID",paste(taxa,"Richness",sep=""))),by=c("sample_ID"),all.x=T)
  #ShannonDiversity
  ZoopFinal<-merge(ZoopFinal,setNames(aggregate(OrderSubOrder~sample_ID,data=SizeID.Order,FUN=function(x) -sum(prop.table(table(x[!is.na(x)]))*log(prop.table(table(x[!is.na(x)]))))),c("sample_ID",paste(taxa,"ShannonDiversity",sep=""))),by=c("sample_ID"),all.x=T)
  #Evenness
  ZoopFinal<-merge(ZoopFinal,setNames(aggregate(OrderSubOrder~sample_ID,data=SizeID.Order,FUN=function(x) -sum(prop.table(table(x[!is.na(x)]))*log(prop.table(table(x[!is.na(x)]))))/log(length(unique(x)))),c("sample_ID",paste(taxa,"Evenness",sep=""))),by=c("sample_ID"),all.x=T)
  
}else{
  ZoopFinal[,c(paste(taxa,"Richness",sep=""),paste(taxa,"ShannonDiversity",sep=""),paste(taxa,"Evenness",sep=""))]<-NA
}


