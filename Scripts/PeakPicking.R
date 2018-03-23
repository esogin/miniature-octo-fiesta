# peak picking 
# Description: Code to select peaks from maldifish data set. Must be performed on high mem node to avoid crashing the machine. 

library(Cardinal)

dir<- "/opt/extern/bremen/symbiosis/sogin/Data/MaldiFish"
setwd(file.path(dir,'Data'))

file<-'20161206_MPIBremen_Bputeoserpentis_MALDI-FISH8_Sl16_s1_DHB_233x233_3um.imzML'
data <- readMSIData(file, mass.accuracy=1)
norm_data<-normalize(data,method='tic')
msdataset<-peakPick(norm_data, method="adaptive", SNR=10)
save(msdataset,file='FullData.PeakPicked.RData')

## END 