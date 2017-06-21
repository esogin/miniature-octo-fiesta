## Explore Data
## EM SOGIN 
## June 17, 2017

## Description: This script will import data processed by the data processing script, it will bin metabolomes based on FISH annotations, and it will run exploratory data analyses including PCA analysis and volcano plots

## Need: (1) Decision on threshold for metabolome binning, how to decide what is noise and what is signal (mean column intensities per group, range is typically from very low to very high, median values aroun 8, mean values around 100)

## (2) Need to add text and further thresholding to volcano plots 



##-------------------------------------------------------
## Set up working space
rm(list=ls())
library(Cardinal)
library(VennDiagram)
library(ggplot2)
library(gplots)
library(gridExtra)
library(dplyr)

## Load data
load('Data/Cardinal_Processed_Data_June2017.RData')
ls()
maldifishmz ## Dimensions: 2322 peaks by 17,273 pixels

##-------------------------------------------------------
## Metabolome Binning 
##-------------------------------------------------------
table(maldifishmz$Class) # Get class information 
trsh<-10 # set a value for the signal to noise intensity for the symbiont binning

# Tissue
channel<-'Host'
df<-int_matrix[which(int_matrix$Class==channel),]
mz.means<-colMeans(df[,grep('m.z.', colnames(df))])
Host<-names(mz.means[which(mz.means > trsh)])

# Mox
channel<-'Mox'
df<-int_matrix[which(int_matrix$Class==channel),]
mz.means<-colMeans(df[,grep('m.z.', colnames(df))])
Mox<-names(mz.means[which(mz.means > trsh)])

# Sox
channel<-'Sox'
df<-int_matrix[which(int_matrix$Class==channel),]
mz.means<-colMeans(df[,sapply(df,is.numeric)])
Sox<-names(mz.means[which(mz.means > trsh)])

# Combine channels into a list 
chns<-list(Host, Mox, Sox)

# Plot venn diagram 
venn<-venn.diagram(x = chns, category.names = c('Host','Mox','Sox'),filename = NULL, fill=c('skyblue', 'lightcoral','mediumseagreen'))
grid.arrange(gTree(children=venn), ncol=1,newpage = T)

## Get unique values per channel 
metlist<-chns
names(metlist)<-c('Host', 'Mox', 'Sox')
a<-venn(metlist, show.plot=F)
inters <- attr(a,"intersections")
m<-data.frame(mz=inters$Mox, label='mox')
s<-data.frame(mz=inters$Sox, label='sox')
h<-data.frame(mz=inters$Host, label='host')
c<-data.frame(mz=inters$`Host:Mox:Sox`, label='common')
vlists<-rbind(m,s,h,c)
write.csv(vlists, file='Results/VennDiagram-Seperated-Ion-List.csv')



##-------------------------------------------------------
## PCA Analysis
##-------------------------------------------------------

## PCA in Cardinal
set.seed(seed = 1)
pca.mod<-PCA(maldifishmz,ncomp=3)
summary(pca.mod)
plot(summary(pca.mod))

## Plot PCA scores across tissue section
mycols<-gradient.colors(10,start='royalblue', end='goldenrod')
pdf('Results/ExploreData_Results/PCA.pdf', height=12, width=18)
par(mfrow=c(1,2))
image(pca.mod, column=c('PC1','PC2'), superpose=F, col.regions=mycols)
dev.off()

## Plot in normal 2D Plot
pca.scores.red<-as.data.frame(pca.mod[[1]]$scores[maldifishmz$Class=='Mox',])
pca.scores.red$Class<-'mox'
pca.scores.green<-as.data.frame(pca.mod[[1]]$scores[maldifishmz$Class=='Sox',])
pca.scores.green$Class<-'sox'
pca.scores.tissue<-as.data.frame(pca.mod[[1]]$scores[maldifishmz$Class=='Host',])
pca.scores.tissue$Class<-'host'
pca.scores<-rbind(pca.scores.green,pca.scores.red,pca.scores.tissue)
pca.scores.symbionts<-pca.scores[pca.scores$Class %in% c('mox','sox'),]

ggplot(pca.scores, aes(x=PC1, y=PC2, color=Class)) + geom_point() + geom_point(data=pca.scores.symbionts,aes( x=PC1, y=PC2, color=Class)) + theme_bw() + scale_color_manual(values=c('skyblue', 'lightcoral','mediumseagreen'))
ggsave('Results/ExploreData_Results/PCA_scatter_PC1_PC2.eps')

ggplot(pca.scores, aes(x=PC2, y=PC3, color=Class)) + geom_point() + geom_point(data=pca.scores.symbionts,aes( x=PC1, y=PC2, color=Class)) + theme_bw() + scale_color_manual(values=c('skyblue', 'lightcoral','mediumseagreen'))
ggsave('Results/ExploreData_Results/PCA_scatter_PC3_PC2.eps')

ggplot(pca.scores, aes(x=PC1, y=PC3, color=Class)) + geom_point() + geom_point(data=pca.scores.symbionts,aes( x=PC1, y=PC2, color=Class)) + theme_bw() + scale_color_manual(values=c('skyblue', 'lightcoral','mediumseagreen'))
ggsave('Results/ExploreData_Results/PCA_scatter_PC1_PC3.eps')


## PCA on Mox pixels

## PCA on Host pixels 


##-------------------------------------------------------
## Volcano Plots
##-------------------------------------------------------
## Get group means
host<-colMeans(int_matrix[int_matrix$Class %in% 'Host',grep('m.z', colnames(int_matrix))], na.rm=T)
mox<-colMeans(int_matrix[int_matrix$Class %in% 'Mox',grep('m.z', colnames(int_matrix))], na.rm=T)
sox<-colMeans(int_matrix[int_matrix$Class %in% 'Sox',grep('m.z', colnames(int_matrix))], na.rm=T)

## Host vs. Mox
## Fold Change B/A
df<-data.frame(host, mox)
df$Metabolite<-rownames(df)
df$FC<-df$mox/df$host

## T-test
Groups<-c('Host', 'Mox')
X<-int_matrix[int_matrix$Class %in% Groups,]
ttest<-lapply(X[-1], function(x) t.test(x ~ X$Class))
pvals<-data.frame(unlist(lapply(ttest, function(x) x$p.value)))
colnames(pvals)<-'pvalue'
pvals$padjust<-p.adjust(p = pvals$pvalue, method='BH')
pvals$Metabolite<-rownames(pvals)

## Combine
df.m<-merge(pvals[,colnames(pvals) %in% c('Metabolite','padjust')], df, by='Metabolite')
df.m$threshold <- as.factor(df.m$padjust < 0.05)

## Plot Volcano 
g1 <- ggplot(data=df.m, aes(x=log2(FC), y =-log10(padjust), 
                colour=threshold)) +
  geom_point(size=1.75) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() +
  theme(legend.position="none") + ggtitle(label = 'Mox vs. Host Volcano Plot')
g1
ggsave(filename = 'Results/ExploreData_Results/Host_Mox_Volcano.eps')

## Mox vs. Sox
## Fold Change B/A
df<-data.frame(sox, mox)
df$Metabolite<-rownames(df)
df$FC<-df$mox/df$sox

## T-test
Groups<-c('Sox', 'Mox')
X<-int_matrix[int_matrix$Class %in% Groups,]
ttest<-lapply(X[-1], function(x) t.test(x ~ X$Class))
pvals<-data.frame(unlist(lapply(ttest, function(x) x$p.value)))
colnames(pvals)<-'pvalue'
pvals$padjust<-p.adjust(p = pvals$pvalue, method='BH')
pvals$Metabolite<-rownames(pvals)

## Combine
df.m<-merge(pvals[,colnames(pvals) %in% c('Metabolite','padjust')], df, by='Metabolite')
df.m$threshold <- as.factor(df.m$padjust < 0.05)

## Plot Volcano 
g2 <- ggplot(data=df.m, aes(x=log2(FC), y =-log10(padjust), 
                           colour=threshold)) +
  geom_point(size=1.75) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() +
  theme(legend.position="none") + ggtitle(label = 'Mox vs. Sox Volcano Plot')
g2
ggsave(filename = 'Results/ExploreData_Results/Mox_Sox_Volcano.eps')

## Sox vs. Host
## Fold Change B/A
df<-data.frame(sox, host)
df$Metabolite<-rownames(df)
df$FC<-df$sox/df$host

## T-test
Groups<-c('Sox', 'Host')
X<-int_matrix[int_matrix$Class %in% Groups,]
ttest<-lapply(X[-1], function(x) t.test(x ~ X$Class))
pvals<-data.frame(unlist(lapply(ttest, function(x) x$p.value)))
colnames(pvals)<-'pvalue'
pvals$padjust<-p.adjust(p = pvals$pvalue, method='BH')
pvals$Metabolite<-rownames(pvals)

## Combine
df.m<-merge(pvals[,colnames(pvals) %in% c('Metabolite','padjust')], df, by='Metabolite')
df.m$threshold <- as.factor(df.m$padjust < 0.05)

## Plot Volcano 
g3 <- ggplot(data=df.m, aes(x=log2(FC), y =-log10(padjust), 
                            colour=threshold)) +
  geom_point(size=1.75) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() +
  theme(legend.position="none") + ggtitle(label = 'Host vs. Sox Volcano Plot') 
g3
ggsave(filename = 'Results/ExploreData_Results/Host_Sox_Volcano.eps')


##-------------------------------------------------------
## Do Clustering Analysis?
##-------------------------------------------------------









## Select Mox Pixels
mox_int<-int_matrix[int_matrix$Class=='Mox',]
head(mox_int[,1:10])

x<-as.matrix(mox_int[sample(x=1:nrow(mox_int), size = 250, replace = F),grep('m.z.', colnames(mox_int))])
x<-x[,colSums(x) > range(colSums(x))[2]*0.0001]
dim(x)

## Visualize heatmap
heatmap.2(log(x+1),trace = 'none')

## Visualize in a pca

## Select Sox Pixels
sox_int<-int_matrix[int_matrix$Class=='Sox',]
head(sox_int[,1:10])

x<-as.matrix(sox_int[sample(x=1:nrow(sox_int), size = 250, replace = F),grep('m.z.', colnames(sox_int))])
x<-x[,colSums(x) > range(colSums(x))[2]*0.001]
dim(x)

## Visualize heatmap
heatmap.2(log(x+1),trace = 'none')

## PCA
sox_int<-int_matrix[int_matrix$Class=='Sox',]
x<-sox_int[,grep('m.z',colnames(sox_int))]
x<-x[,!colSums(x) == 0]

mpca<-prcomp(x, scale=T, center=T)
scores<-data.frame(mpca$x[,1:3])
plot(scores)

scores$rn<-as.vector(rownames(scores))
sox_int$rn<-as.vector(rownames(sox_int))

## merge hopanoid value with PC scores
sox_hop<-data.frame(rn=sox_int$rn, hop1=sox_int$m.z...546.526, ce=sox_int$m.z...756.554)
sc.m<-merge(scores, sox_hop, by= 'rn')
ggplot(sc.m, aes(x=PC1, y=PC2, color=ce)) + geom_point() 

## Host
## Select host Pixels
host_int<-int_matrix[int_matrix$Class=='Host',]
head(host_int[,1:10])

x<-as.matrix(host_int[sample(x=1:nrow(host_int), size = 250, replace = F),grep('m.z.', colnames(host_int))])
x<-x[,colSums(x) > range(colSums(x))[2]*0.001]
dim(x)

## Visualize heatmap
heatmap.2(log(x+1),trace = 'none')

## PCA
x<-host_int[,grep('m.z',colnames(host_int))]
x<-x[,!colSums(x) == 0]

mpca<-prcomp(x, scale=T, center=T)
scores<-data.frame(mpca$x[,1:3])
plot(scores)
scores$rn<-as.vector(rownames(scores))
host_int$rn<-as.vector(rownames(host_int))

## merge hopanoid value with PC scores
host_hop<-data.frame(rn=host_int$rn, hop1=host_int$m.z...546.526, ce=host_int$m.z...756.554)
sc.m<-merge(scores, host_hop, by= 'rn')
ggplot(sc.m, aes(x=PC1, y=PC2, color=ce)) + geom_point() 



## All
## Select host Pixels
x<-int_matrix[sample(x=1:nrow(int_matrix), size = 250, replace = F),]
x<-x[,colSums(x) > range(colSums(x))[2]*0.001]
dim(x)

## Visualize heatmap
heatmap.2(log(as.matrix(x[,2:ncol(x)])+1),trace = 'none')

## PCA
mpca<-prcomp(x, scale=T, center=T)
scores<-data.frame(mpca$x[,1:3])
plot(scores)
scores$rn<-as.vector(rownames(scores))
host_int$rn<-as.vector(rownames(host_int))

## merge hopanoid value with PC scores
host_hop<-data.frame(rn=host_int$rn, hop1=host_int$m.z...546.526, ce=host_int$m.z...756.554)
sc.m<-merge(scores, host_hop, by= 'rn')
ggplot(sc.m, aes(x=PC1, y=PC2, color=ce)) + geom_point() 







