##Cluster data
## EM SOGIN 
## June 6 2017
##--------------------------------------------------------------------------------------------------
# Clustering
#
## Description: Goal of script is to cluster dataset with cardinal functions using an unsupervised, but sptially aware technique
#
rm(list=ls())
#
dir<- "/opt/extern/bremen/symbiosis/sogin/Data/MaldiFish"
setwd(file.path(dir,'Scripts'))
#
## Load libraries
library(Cardinal)
#
## Load in dataset
setwd(file.path(dir,'Results'))
load('Cardinal_Processed_Data_June2017.RData')
##--------------------------------------------------------------------------------------------------
## Cluster Analysis
#
# Spatial Aware K-means using the shrunken centriods method
#
ssc.a<-spatialShrunkenCentroids(maldifishmz,r=c(1,2),k=c(10,15,20),s=c(0,6,12,18,24),method='adaptive')
#
setwd(file.path(dir, 'Results'))
save(list=c('ssc.a'), file='cluster-analysis-final.RData')
## END

