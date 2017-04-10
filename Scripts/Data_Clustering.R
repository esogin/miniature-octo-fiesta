# Clustering 

## Description: Goal of script is to cluster dataset with cardinal functions using both supervised and unsuprvised tools. To be done on high mem cluster to avoid crashing local machine. Output will be as .RData files in the results directory

rm(list=ls())

dir<- "/opt/extern/bremen/symbiosis/sogin/Data/MaldiFish"
setwd(file.path(dir,'Scripts'))

## Load libraries
library(Cardinal)

## Load in dataset 
setwd(file.path(dir,'Results'))
load('Cardinal_Processed_Data.RData')

## Take out Mixed Category 
maldifishmz_3grps<-maldifishmz[,maldifishmz$Class %in% c('Tissue','Red','Green')]

## Add in crossvalidation group
maldifishmz_3grps$cvgroup<-cut(maldifishmz_3grps$y, breaks=10, labels = F) #based on y-location value in dataset

## Generate Clusters with K-Means spatial clustering algrithm 

## Unspervised
skmg<-spatialKMeans(maldifishmz_3grps, r = c(1, 2), k = c(3,5,7), method = "gaussian")
skma<-spatialKMeans(maldifishmz_3grps, r = c(1, 2), k = c(3,5,7), method = "adaptive")

setwd(file.path(dir,'Results'))
save(list=c('skmg','skma'), file='unsupervised-clustering-analysis.RData')


## Do Spatial shrunken centriods clustering method

ssc.a<-spatialShrunkenCentroids(maldifishmz_3grps,r=c(1,2),k=c(10,15,20),s=c(0,5,10,15),method='adaptive')
save(list=c('ssc.a'), file='ssc-adaptive-unsupervised-clustering-analysis.RData')

## Supervised
ssca.cv<-cvApply(maldifishmz_3grps, .y=maldifishmz_3grps$Class, .fun="spatialShrunkenCentroids",method='adaptive', .fold=cvgroup, r=c(1,2),s=c(0,3,5,10,15))

setwd(file.path(dir,'Results'))
save(list=c('ssca.cv'), file='supervised-clustering-anlaysis.RData')

## transfer files to local machine

## END