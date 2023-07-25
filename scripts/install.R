#Install packages

install.packages("pacman")
library(pacman)

#read in libraries
pacman::p_load(dplyr, vegan, labdsv, goeveg, rLakeAnalyzer, 
               ggplot2,tidyr,lubridate, scales, colorblindcheck, 
               viridis, RColorBrewer, plotrix, BiocManager, ggpubr,
               remotes, egg, splancs, FSA, rcompanion, ggrepel)

remotes::install_github("jfq3/ggordiplots")

#------------------------------------------------------------------#
#functions

#Count characters starting at the end of the string
substrEnd <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Calculates the standard error####
stderr <- function(x) {
  sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
}

#calculating SE using difference between epi mean and hypo mean
SE.diffMean<- function(x,y){
  sqrt((sd(x,na.rm=TRUE)^2/length(na.omit(x))) + 
         (sd(y,na.rm=TRUE)^2/length(na.omit(y))))
}

#Function to get the median if odd number in a column, or one entry below median if even number
median.zoop<-function(vector){
  #Remove nas and then sort
  vector2<-vector[!is.na(vector)]
  vector2<-sort(vector2)
  #If the vector contains an odd number of entries, return the median
  if(length(vector2)%%2==1){return(median(vector2))}else{
    return(vector2[length(vector2)/2]) #If the vector is even #, return the entry immediately below the median
  }
}

