# Clustering 

## Description: Goal of script is to cluster dataset with cardinal functions using both supervised and unsuprvised tools. To be done on high mem cluster to avoid crashing local machine. Output will be as .RData files in the results directory

# rm(list=ls())
# 
# dir<- "/opt/extern/bremen/symbiosis/sogin/Data/MaldiFish"
# setwd(file.path(dir,'Scripts'))
# 
# ## Load libraries
# library(Cardinal)
# 
# ## Load in dataset 
# setwd(file.path(dir,'Results'))
# load('Cardinal_Processed_Data_June2017.RData')
# 
# ## Take out Mixed Category 
# maldifishmz_3grps<-maldifishmz[,maldifishmz$Class %in% c('Tissue','Red','Green')]
# 
# ## Add in crossvalidation group
# maldifishmz_3grps$cvgroup<-cut(maldifishmz_3grps$y, breaks=6, labels = F) #based on y-location value in dataset
# 
# ## Generate Clusters with K-Means spatial clustering algrithm 
# 
# ## Unspervised
# skmg<-spatialKMeans(maldifishmz_3grps, r = c(1, 2), k = c(3,5,7), method = "gaussian")
# skma<-spatialKMeans(maldifishmz_3grps, r = c(1, 2), k = c(3,5,7), method = "adaptive")
# 
# setwd(file.path(dir,'Results'))
# save(list=c('skmg','skma'), file='unsupervised-clustering-analysis.RData')
# 
# 
# ## Do Spatial shrunken centriods clustering method
# ssc.a<-spatialShrunkenCentroids(maldifishmz_3grps,r=c(1,2),k=c(15,20),s=c(0,3,6,9,15,20),method='adaptive')
# save(list=c('ssc.a'), file='ssc-adaptive-unsupervised-clustering-analysis.RData')
# 
# ##--------------------------------------------------------------------------------------------------
# ## Supervised Clustering Method 
# 
# 
# ## Add in Ciliated Edge Group for supervised analysis 
# grps<-skma$cluster$`r = 2, k = 7`
# ciliatedgrps<-grps[which(grps==1)]
# grp1<-names(ciliatedgrps)
# msidata<-maldifishmz_3grps
# pData(msidata)[rownames(pData(msidata)) %in% grp1,'Class']<-'CiliatedEdge'
# 
# table(msidata$Class, msidata$cvgroup)
# 
# ssca.cv<-cvApply(msidata, .y=msidata$Class, .fun="spatialShrunkenCentroids",method='adaptive', .fold=cvgroup, r=c(1,2,3),s=c(0, 4, 8, 12, 16, 20, 24, 28))
# 
# setwd(file.path(dir,'Results'))
# save(list=c('ssca.cv', 'msidata'), file='supervised-clustering-anlaysis.RData')

## transfer files to local machine

## END



# Clustering
#
## Description: Goal of script is to cluster dataset with cardinal functions using both supervised and unsuprvised$
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
#load('Cardinal_Processed_Data.RData')
load('Cardinal_Processed_Data_June2017.RData')
##--------------------------------------------------------------------------------------------------
## Cluster Analysis
#
# Spatial Aware K-means using the shrunken centriods method
#
ssc.a<-spatialShrunkenCentroids(maldifishmz,r=c(1,2),k=c(20,40),s=c(0,6,18,24),method='adaptive')
#
## Spatial aware K Means without centriods
skma<-spatialKMeans(maldifishmz, r = c(1, 2), k = c(6,12,18), method = "adaptive")
#
setwd(file.path(dir, 'Results'))
save(list=c('ssc.a', 'skma'), file='cluster-analysis-final.RData')
## END



