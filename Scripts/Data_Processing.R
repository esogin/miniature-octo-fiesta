#Data Processing
## Description: This script will import data from peak-picking pipeline and reduce the dimensions by subsetting by pixel category and then reducing the number of m/z values in dataset
## All data dimensionality is done in cardinal 

dir<- "/opt/extern/bremen/symbiosis/sogin/Data/MaldiFish"
setwd(file.path(dir,'Scripts'))

## Load libraries
library(Cardinal)

## Load Maldi Data
setwd(file.path(dir, 'Data'))
load('FullData.PeakPicked.RData') ## Peak Picked Dataset

## Load fish data
FISH.All<-read.csv('pixels.csv') ## 
FISH.All<-FISH.All[order(FISH.All$y, FISH.All$x),]

## Comind FISH and MALDI 
pData(msdataset)$Class<-as.vector(FISH.All$Category)
identical(pData(msdataset)$Class, as.vector(FISH.All$Category)) ## another sanity check 

## Remove Background signals 
maldifish<-msdataset[,msdataset$Class %in% c('Green', 'Red', 'Tissue','Mixed')]  ## Select Only Annotated Pixels
#maldifish<-msdataset

## Subset data by frequency of mz value
threshold<-min(table(maldifish$Class))*0.1 
peaks<-pData(mzData(imageData(maldifish)))
mz<-unlist(peaks)
tab<-as.data.frame(table(mz))
peaks_in_data<-tab[tab$Freq >=threshold,'mz']
maldifishmz<-maldifish[mz(maldifish) %in% peaks_in_data,] 

mat<-as.data.frame(spectra(maldifishmz)[,]) ## This will take awhile
rownames(mat)<-maldifishmz@imageData@dimnames[[1]] ## gets rownames
colnames(mat)<-maldifishmz@imageData@dimnames[[2]]
df<-as.data.frame(t(mat))
int_matrix<-data.frame(Class=maldifishmz$Class,df)


setwd(file.path(dir,'Results'))
save(list=c('maldifishmz','peaks_in_data','int_matrix','df'),file='Cardinal_Processed_Data.RData')
#save(list=c('maldifishmz','peaks_in_data','int_matrix','df'),file='Cardinal_Processed_Data_withbkg.RData')


## END CODE











## PCA Analysis 
pdf('PCAAnalysis.pdf')
pca <- PCA(maldifish, ncomp = 5)
summary(pca)
image(pca, column = c("PC1", "PC2", "PC3"), superpose = FALSE,layout = c(3, 1))
plot(pca, column = c("PC1", "PC2", "PC3"), superpose = FALSE, layout = c(3,1))
dev.off()
##------------------------------------------------------------------------------------------------------------------
## Spatial K Means 
pdf('KMeansClustering.pdf')
clust<-spatialKMeans(maldifishmz,r=c(1,2), k=c(3,4,8), method='adaptive')
summary(clust)
image(clust, key = T, layout = c(2,3))
plot(clust, key = T, layout = c(2,3), type = c("p", "h"), pch = 20)
dev.off()
##------------------------------------------------------------------------------------------------------------------

## Spatial shrunken centriods

ssc<-spatialShrunkenCentroids(maldifishmz, r=c(1,2),k=c(5,10),s = c(0, 1,3, 5, 7), method = "adaptive")



pdf('spatialShrunkenCentriods.pdf')
image(ssc, model = list(r =2, s = c(0,3)), key = FALSE, layout = c(1,2))
dev.off()

##------------------------------------------------------------------------------------------------------------------

## Discriminate Analysis 

## Even out groups 




pdf('OPLSAnalysis.pdf')

pData(maldifishmz)$pix<-as.vector(paste('x = ',maldifishmz$x,', y = ',maldifishmz$y))
opls.mod<-OPLS(maldifishmz, y=maldifishmz$Class, ncomp=1:2)

opls <- cvApply(maldifishmz, .y = maldifishmz$Class, .fun = "OPLS", ncomp = 1:2, keep.Xnew = FALSE, .fold=pix)
plot(summary(opls))
(summary(opls))
image(opls, model = list(ncomp = 12), layout = c(4, 2))
(topLabels(opls))

dev.off()

save.image(file="Classification.RData")