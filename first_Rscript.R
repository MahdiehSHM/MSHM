
###################################
###################################
## Pakages we need for the analysis
library(mvabund)
library(vegan)
library(coda)
library(rjags)
library(boral)
library(effects)
library(MASS)
library(reshape)
library(jpeg)
##################################
##################################
## DATA Input

data <- read.csv("MSHM.csv",header = T, row.names = 1)
otu <- read.csv(file="OTU.csv",header = T, row.names = 1)
str(data)
str(otu)
summary(otu)
summary(data)
# Temperature as fector:
data$TEMP<-as.factor(data$TEMP)
######## MERG the OTU abundance data by 4
S1<-seq(1,121920,4)
S2<-seq(4,121920,4)
O1<-matrix(0,length(S1),133)
for (i in 1:length(S1)) {
  
O1[i,]<-colSums(otu[S1[i]:S2[i],])}

#now convert to data frame and rename the columns
OTUabund<-data.frame(O1)
class(OTUabund)
colnames(OTUabund)=colnames(otu)
name<- row.names(data)
row.names(OTUabund)<- name[1:30480]
summary(OTUabund)
# the abundance data is now saved under object name: OTUabund

# one OTU has 0 in all obs "APE.se5.Staphylotrichum.coccosporum". check this
# if you want to delet this OTU in OTUabund run the code below:
#OTUabund$APE.se5.Staphylotrichum.coccosporum<-NULL

## Merging METADATA
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
# the merged metadata is now saved under object named: MetaData

# calculating Isolation Rate for each sample, We need this for diversity estimations
isolationRate= apply(OTUabund,1, sum)
MetaData = cbind(MetaData, IR = isolationRate)
# now check the MetaData object to see that a new variable is added "IR"


##############################################
##### IR histograms for each variable/factor
##############################################






##############################################################
##### Diversity Indices
##############################################################

metaDATA <- read.csv("MSHM.csv")
OTUabund <- read.csv("OTU.csv")

