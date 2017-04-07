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
maldifishmz<-maldifishmz[,maldifishmz$Class %in% c('Tissue','Red','Green')]

## Add in crossvalidation group
pixelNo<-nrow(pData(maldifishmz))
randNum<-sample(1:10,size = pixelNo,replace = T) ## For now, just split the data 10 ways, increase to see how changes in overnight run
pData(maldifishmz)$cvgroup<-as.vector(randNum)

## Generate Clusters with K-Means spatial clustering algrithm 

## Unspervised
skmg<-spatialKMeans(maldifishmz, r = c(1, 2), k = c(3,5,7), method = "gaussian")
skma<-spatialKMeans(maldifishmz, r = c(1, 2), k = c(3,5,7), method = "gaussian")

setwd(file.path(dir,'Results'))
save(list=c('skmg','skma'), file='unsupervised-clustering-analysis.RData')

ssc.g<-spatialShrunkenCentroids(maldifishmz,r=c(1,2),kl=c(15,20),method='gaussian')
ssc.a<-spatialShrunkenCentroids(maldifishmz,r=c(1,2),kl=c(15,20),method='adaptive')
save(list=c('ssc.a','ssc.g'), file='ssc-unsupervised-clustering-analysis.RData')

## Supervised
sscg.cv<-cvApply(maldifishmz, .y=maldifishmz$Class, .fun="spatialShrunkenCentroids",method='gaussian', .fold=cvgroup, r=c(1,2,3),s=c(0,2,3,5,10))
ssca.cv<-cvApply(maldifishmz, .y=maldifishmz$Class, .fun="spatialShrunkenCentroids",method='adaptive', .fold=cvgroup, r=c(1,2,3),s=c(0,2,3,5,10))

setwd(file.path(dir,'Results'))
save(list=c('sscg.cv', 'ssca.cv'), file='supervised-clustering-anlaysis.RData')

## transfer files to local machine

## END