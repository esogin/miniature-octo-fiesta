##Data Processing
## Description: This script will import data from peak-picking pipeline and reduce the dimensions by subsetting by pixel category and then reducing the number of m/z values in dataset
## All data dimensionality is done in cardinal 

dir<- "/opt/extern/bremen/symbiosis/sogin/Data/MaldiFish"
#dir<- "/Users/esogin/Documents/Projects/maldifish/RAnalysis"

setwd(file.path(dir,'Scripts'))

## Load libraries
library(Cardinal)

## Load Maldi Data
setwd(file.path(dir, 'Data'))
load('FullData.PeakPicked.RData') ## Peak Picked Dataset

## Load fish data
FISH<-read.csv('RGBchannelIntensities.csv') ## 
FISH$XY<-paste('x = ',FISH$x, ', y = ',FISH$y,sep='')
FISH$Annotation<-NA

## Set up annotation rules for FISH signals 
FISH[FISH$mox==0 & FISH$sox ==  0 & FISH$host_dapi==0,'Annotation']<-'Background' ## Identify Background Pixels
FISH[FISH$mox > 20,'Annotation']<-'Mox'
FISH[FISH$sox > 20,'Annotation']<-'Sox'
FISH[FISH$mox > 20 & FISH$sox > 20,'Annotation']<- 'Mixed'
FISH[FISH$host_dapi > 20 & FISH$mox < 20 & FISH$sox < 20,'Annotation']<-'Host'

table(FISH$Annotation)
FISH<-FISH[order(FISH$y, FISH$x),]

## Comind FISH and MALDI 
pData(msdataset)$Class<-as.vector(FISH$Annotation)
identical(pData(msdataset)$Class, as.vector(FISH$Annotation)) ## another sanity check 

## Remove Background signals
maldifish<-msdataset[,msdataset$Class %in% c('Sox', 'Mox', 'Host','Mixed')]  ## Select Only Annotated Pixels
#maldifish<-msdataset
# maldifish dimensions 275279 x 17273

## Subset data by frequency of mz value
threshold<-min(table(maldifish$Class))*0.1
peaks<-pData(mzData(imageData(maldifish)))
mz<-unlist(peaks)
tab<-as.data.frame(table(mz))
peaks_in_data<-tab[tab$Freq >=threshold,'mz']
maldifishmz<-maldifish[mz(maldifish) %in% peaks_in_data,]
## dimensions 2322 x 17273 Hashmat 

mat<-as.data.frame(spectra(maldifishmz)[,]) ## This will take awhile
rownames(mat)<-maldifishmz@imageData@dimnames[[1]] ## gets rownames
colnames(mat)<-maldifishmz@imageData@dimnames[[2]]
df<-as.data.frame(t(mat))
int_matrix<-data.frame(Class=maldifishmz$Class,df)

setwd(file.path(dir,'Results'))
save(list=c('maldifishmz','peaks_in_data','int_matrix','df'),file='Cardinal_Processed_Data_June2017.RData')



## END CODE


