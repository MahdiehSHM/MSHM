############################################
# desert endophytes
############################################
library(vegan)
library(mvabund)
library(MASS)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(tidytree)
library(gtable)
library(grid)
library(magrittr)
library(scales)
library(reshape)
library (lattice)


data <- read.csv("MSHM-1.csv",header = T, row.names = 1)
otu <- read.csv(file="OTU.csv",header = T, row.names = 1)
str(data)
str(otu)
summary(otu)
summary(data)

######## Temperature as fector:
data$TEMP<-as.factor(data$TEMP)
######## MERG the OTU abundance data by 4
S1<-seq(1,121920,4)
S2<-seq(4,121920,4)
O1<-matrix(0,length(S1),133)
for (i in 1:length(S1)) {
  O1[i,]<-colSums(otu[S1[i]:S2[i],])}

########now convert to data frame and rename the columns
OTUabund<-data.frame(O1)
class(OTUabund)
colnames(OTUabund)=colnames(otu)
name<- row.names(data)
row.names(OTUabund)<- name[1:30480]
summary(OTUabund)
########the abundance data is now saved under object name: OTUabund

######## Merging METADATA
S1<-seq(1,121920,4)
S2<-seq(4,121920,4)
D<-matrix(0,length(S1),7)
for (i in 1:length(S1)) {
  D[i,1]<-noquote(paste(data[S2[i],1]))
  D[i,2]<-noquote(paste(data[S2[i],2]))
  D[i,3]<-noquote(paste(data[S2[i],3]))
  D[i,4]<-noquote(paste(data[S2[i],4]))
  D[i,5]<-noquote(paste(data[S2[i],5]))
  D[i,6]<-noquote(paste(data[S2[i],6]))
  D[i,7]<-noquote(paste(data[S2[i],7]))
}
MetaData<-data.frame(D)
class(MetaData)
colnames(MetaData)=colnames(data)
row.names(MetaData)<- name[1:30480]
summary(MetaData)
######## the merged metadata is now saved under object named: MetaData

######## calculating Isolation Rate for each sample, We need this for diversity estimations
isolationRate= apply(OTUabund,1, sum)
MetaData = cbind(MetaData, IR = isolationRate)
######## now check the MetaData object to see that a new variable is added "IR"
###Exporting as csv files

write.csv(MetaData,file="MataDataMerg.csv")
write.csv(OTUabund, file="OTUabundMerg.csv")
##############################################
##### IR histograms for each variable/factor
##############################################

hist(MetaData $ IR) 
hist(log(MetaData $ IR))
hist(log(MetaData $ IR), prob=TRUE)
hist(log(MetaData $ IR), prob=TRUE, breaks=20)

boxplot(IR ~ SITE, data = MetaData)
boxplot(IR ~ SOIL, data = MetaData)
boxplot(IR ~ TIME, data = MetaData)
boxplot(IR ~ HOST, data = MetaData)
boxplot(IR ~ TISSUE, data = MetaData)
boxplot(IR ~ TEMP, data = MetaData)
boxplot(IR ~ MEDIA, data = MetaData)



aggregate(OTUabund $SRE.ss.4.Neocamarosporium.sp. ~HOST, data = MetaData, sum)
aggregate(OTUabund $SRE.ss.4.Neocamarosporium.sp. ~TISSUE, data = MetaData, sum)
aggregate(OTUabund $SRE.ss.4.Neocamarosporium.sp. ~SOIL, data = MetaData, sum)
aggregate(OTUabund $SRE.ss.4.Neocamarosporium.sp. ~SITE, data = MetaData, sum)
aggregate(OTUabund $SRE.ss.4.Neocamarosporium.sp. ~TIME, data = MetaData, sum)
aggregate(OTUabund $SRE.ss.4.Neocamarosporium.sp. ~TEMP, data = MetaData, sum)
aggregate(OTUabund $SRE.ss.4.Neocamarosporium.sp. ~MEDIA, data = MetaData, sum)

summary(MetaData)







