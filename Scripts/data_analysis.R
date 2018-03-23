## Description: This script will import data processed by the data processing script, it will bin metabolomes based on FISH annotations, and it will run exploratory data analyses including PCA analysis and volcano plots

##-------------------------------------------------------
## Set up working space
rm(list=ls())
library(Cardinal)
library(VennDiagram)
library(ggplot2)
library(gplots)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(ggdendro)

dir<-'~/Documents/Projects/maldifish/RAnalysis/'

## Import Data
setwd(file.path(dir))
load('Data/Cardinal_Processed_Data_June2017.RData')
load('Results/cluster-analysis-final.RData')
ls()
maldifishmz ## Dimensions: 2322 peaks by 17,273 pixels

setwd(file.path(dir,'Results'))
## Extra functions
##-------------------------------------------------------
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
summaryPlots <- function(dataset, results, model, segment, name, col) {
  image(results, model = model, key = FALSE, column = segment, main = name,
        layout = c(3, 2), col = col)
  plot(results, model = model, key = FALSE, column = segment, mode = "centers",
       main = "Shrunken mean spectrum", col = col)
  plot(results, model = model, key = FALSE, column = segment, mode = "tstatistics",
       main = "Shrunken t-statistics", col = col)
  top <- topLabels(results, n = 3, model = model, filter = list(classes = segment))
  image(dataset, mz = top$mz, contrast.enhance = "histogram")
}

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
summary(mz.means)
Mox<-names(mz.means[which(mz.means > trsh)])

# Sox
channel<-'Sox'
df<-int_matrix[which(int_matrix$Class==channel),]
mz.means<-colMeans(df[,sapply(df,is.numeric)])
summary(mz.means)
Sox<-names(mz.means[which(mz.means > trsh)])

# Combine channels into a list 
chns<-list(Host, Mox, Sox)

# Plot venn diagram 
venn<-venn.diagram(x = chns, category.names = c('Host','Mox','Sox'),filename = NULL, fill=c('skyblue', 'lightcoral','mediumseagreen'))

pdf('metabolome_binning.pdf')
grid.arrange(gTree(children=venn), ncol=1,newpage = T)
dev.off()

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
write.csv(vlists, file='VennDiagram-Seperated-Ion-List.csv')


##-------------------------------------------------------
## Cluster Results
##-------------------------------------------------------
## Cluster analysis was ran on the high mem node in colonge following the cluster_data.R script

## Visualize ssc.a cluster results and choose "best" fitting model
plot.new()
plot(summary(ssc.a))
abline(h=7, col='blue',lty=3)

summary(ssc.a)

# All models converge on approximately 9-10 clusters describing the data-set. In order to retain as much information as possible, the model with r=1, k=20 and s=18 will be kept for down stream analysis. 
mod<-list(r=2, k=10, s=18) ## Choose a model -- 10 groups, 124 features 

pdf('cluster_analysis.pdf', height=10, width = 10)
mycol<-c('olivedrab3', 'mediumslateblue','mediumblue','coral4','darkslategray','hotpink','mediumseagreen','darkgoldenrod','gray25','mediumblue',rep('gray',11))
mycol<-c('olivedrab3','mediumblue','lightcyan3','indianred3','deepskyblue','hotpink','mediumseagreen',rep('gray',11))
image(ssc.a, model=mod,col=mycol,strip=F)
dev.off()

## Plot clusters independently 
pdf('Cluster-Distribution-Plots.pdf', height=6, width=11)
for (i in 1:7){
  summaryPlots(maldifishmz,ssc.a,model=mod,segment=i, name='', col=mycol)
}
dev.off()

##_______________ Dendrogram of Clusters _________________________
modResults<-ssc.a@resultData$`r = 2, k = 10, s = 18`

## Add in cluster group to spectra information 
maldifishmz$SSCA_Clusters_pixeln<-names(modResults$classes)
maldifishmz$SSCA_Clusters<-paste('ssca_cluster',modResults$classes,sep='-')

## Get mean spectra for each group
mean.sp<-data.frame(featureApply(maldifishmz, mean, .pixel.groups=maldifishmz$SSCA_Clusters))
mean.sp[,1:10]

## Run hclust using bray curtis distances 
d<-vegdist(mean.sp,method = 'bray')
clusters<-hclust(d, method='average')
plot(clusters,main='bray')

dhc <- as.dendrogram(clusters)
ddata <- dendro_data(dhc, type = "rectangle")
ddata$labels$label<-as.numeric(gsub('ssca_cluster-',' ',ddata$labels$label))
ddata$labels<-ddata$labels[order(ddata$labels$label),]

ggplot() + geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = ddata$labels, aes(x = x, y = y, label = label, size = 3, vjust = 0.5, hjust=-0.5)) + 
  geom_point(data=label(ddata), aes(x=x, y=y, color=as.factor(label),size=3),shape=15) + 
  theme_dendro() + scale_y_reverse() + coord_flip() + scale_color_manual(values=mycol)  + guides(color=F, size=F) 
ggsave('Cluster_dendrogram_segmentationmap.eps')


##_______________   Relate to FISH Data  _________________________
## Get pixel data 
pDf<-data.frame(pData(maldifishmz))
pDf$rowname<-rownames(pDf)

## Get cluster data plus coordinates
cDf<-data.frame(Cluster=ssc.a$classes$`r = 2, k = 10, s = 18`)
cDf$rowname<-rownames(cDf)

## merge cluster annotations and pixel data
df.m<-merge(pDf, cDf, by = 'rowname')

## Call all Red and Green annotations symbionts
df.m[df.m$Class %in% c('Mox','Sox','Mixed'),'Class']<-'Symbiont'


## Class membership across clusters
tab<-table(df.m$Class , df.m$Cluster)

## Get Symbiont ratio per cluster group 
sym.ratio<-t(tab)

sr<-data.frame(host=sym.ratio[,1], symbiont=sym.ratio[,2], Cluster=rownames(sym.ratio))
sr$ratio_s_to_h<-sr$symbiont/sr$host
sr
setwd(file.path(dir,'Results'))
write.csv(sr, 'Symbiont_to_host_ratio.csv')

total<-rowSums(tab)
tab<-data.frame(tab)
tab$Total<-total
tab$percent_of_group<-tab$Freq/total*100
tab<-tab[!tab$Freq==0,]
tab$Var2<-factor(tab$Var2, levels=c(3,6,4,1,5,2,7))

head(tab)


ggplot(tab, aes(x=as.factor(Var2), y=percent_of_group, fill=Var1)) + geom_bar(stat='identity', position='dodge')+ scale_fill_manual(values=c('gray80','black')) + ggtitle('% of all symbiont or host tissue identification signals across clusters') + ylab('% of  "FISH" Group') + xlab('Cluster Group') + guides(fill=guide_legend(title='FISH signal')) + theme_bw() + coord_flip()
ggsave('Cluster_Membership.eps')

## _______________  Get Important Ions   _________________________

tL<-topLabels(ssc.a, model=mod, n=10000) # just set a really high number
tL.sig<-tL[which(tL$adj.p.values < 0.2),]
dim(tL.sig)
table(tL.sig$classes)
write.csv(tL.sig, file='ssc-ions-clusterassignment.csv')
