#Dear Mahdie
#Please put comments before each code and specify what you whant to do otherwise I wont understand 
#what is the purpos of the codes!


<<<<<<< HEAD
str(MSHM)
str(OTU)

MSHM1 <- read.csv("MSHM.csv", stringsAsFactors = FALSE)
str(MSHM1)

OTU1 <- read.csv("OTU.csv", stringsAsFactors = FALSE)
str(OTU1)

summary(MSHM) 
summary(OTU)
=======
#I changed the names of the datasets so they would be meaningfull.
# this is a good start keep going on...
# I commented out the codes you wrote becuse they are not relly helpfull... I suggest you delete them all
# instead: try writing the codes to install the neccesary packages that we need... like vegan and mvabund...
# go on step by step on my Rscript and ask if you have problems...
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



<<<<<<< HEAD
metaDATA <- read.csv("MSHM.csv")
OTUabund <- read.csv("OTU.csv")
>>>>>>> f150812809de455a3f981bc9aa47915baa111858


str(MSHM)
str(OTU)

<<<<<<< HEAD
head(OTU1)
head(OTU1, n = 10)
tail(OTU1)
tail(OTU1, n = 10)
names(OTU1)
nrow(OTU1)
ncol(OTU1)
dim(OTU1)

MSHM$Soil.Types
OTU$Sample.ID.ORIGINAL
MSHM$Soil.Types
OTU$Sample.ID.ORIGINAL


subset(OTU,TPEsh28.Humicola.fuscoatra > 0)
subset(OTU,TPEsh28.Humicola.fuscoatra > 1)
subset(OTU,TPEsh28.Humicola.fuscoatra > 2)

subset(OTU,TPEsh28.Humicola.fuscoatra > 3)
View(subset(OTU, TPEsh28.Humicola.fuscoatra > 2))
TPEsh28.Humicola.fuscoatra <-subset(OTU,TPEsh28.Humicola.fuscoatra > 0)
nrow(TPEsh28.Humicola.fuscoatra)

subset(MSHM, Soil.Types == "Gypsum Soil")
View(subset(MSHM, Soil.Types == "Gypsum Soil"))
GypsumSoil <- subset(MSHM, Soil.Types == "Gypsum Soil")

subset(MSHM, Soil.Types == "Gypsum Soil" & Plant.segments == "Root")
subset(MSHM, Soil.Types == "Gypsum Soil" & Plant.segments == "Root")
subset(MSHM, Soil.Types %in% c("Gypsum Soil","Arid soil") & Plant.segments == "Root")

table(MSHM$Soil.Types)
summary(MSHM$Soil.Types)
table (MSHM$Soil.Types , MSHM$Plant.segments)


median(OTU$TPEsh28.Humicola.fuscoatra)
sd(OTU$TPEsh28.Humicola.fuscoatra)
range(OTU$TPEsh28.Humicola.fuscoatra)
sum(OTU$TPEsh28.Humicola.fuscoatra)
mean(OTU$TPEsh28.Humicola.fuscoatra)


table(MSHM$Soil.Types)
Soiltable <-table(MSHM$Soil.Types)
prop.table(Soiltable, margin = 1)
prop.table(Soiltable, margin = 1)

aggregate(ZEE.se11.Periconia.macrospinosa ~ TPEsh28.Humicola.fuscoatra, data=OTU, median)

table (OTU$ZEE.se11.Periconia.macrospinosa)
table (OTU$TPEsh28.Humicola.fuscoatra)
table (OTU$TPEsh28.Humicola.fuscoatra + OTU$ZEE.se11.Periconia.macrospinosa)


plot(x = MSHM$Soil.Types, y = MSHM$ Host.plant.species)
plot(Host.plant.species ~ Soil.Types, data=MSHM)

plot(Host.plant.species ~ Soil.Types, data=MSHM, 
     main = "Frist plot",
     ylab = "Host plant species",
     xlab = "Soil Types")


hist(OTU$ZEE.se11.Periconia.macrospinosa) 
hist(log(OTU$ZEE.se11.Periconia.macrospinosa)) 
hist(log(OTU$ZEE.se11.Periconia.macrospinosa), prob=TRUE) 
hist(log(OTU$ZEE.se11.Periconia.macrospinosa), prob=TRUE, breaks=20) 


boxplot(Temperature ~ Soil.Types , data=MSHM, 
        main="boxplot1",
        ylab = "Host plant species",
        xlab = "Soil Types")

library(ggplot2)
=======
#you want the string to be saved as factors. so you dont need the rest of the codes... DELETE them
#MSHM1 <- read.csv("MSHM.csv", stringsAsFactors = FALSE)
# str(MSHM1)
# 
# OTU1 <- read.csv("OTU.csv", stringsAsFactors = FALSE)
# str(OTU1)
# summary(MSHM) 
# summary(OTU)
# summary(MSHM1) 
# summary(OTU1)
# 
# head(MSHM)
# head(MSHM, n = 10)
# tail(MSHM)
# tail(MSHM, n = 10)
# names(MSHM)
# nrow(MSHM)
# ncol(MSHM)
# dim(MSHM)
# 
# head(MSHM1)
# head(MSHM1, n = 10)
# tail(MSHM1)
# tail(MSHM1, n = 10)
# names(MSHM1)
# nrow(MSHM1)
# ncol(MSHM1)
# dim(MSHM1)
# 
# head(OTU)
# head(OTU, n = 10)
# tail(OTU)
# tail(OTU, n = 10)
# names(OTU)
# nrow(OTU)
# ncol(OTU)
# dim(OTU)
# 
# head(OTU1)
# head(OTU1, n = 10)
# tail(OTU1)
# tail(OTU1, n = 10)
# names(OTU1)
# nrow(OTU1)
# ncol(OTU1)
# dim(OTU1)
>>>>>>> f150812809de455a3f981bc9aa47915baa111858
=======
>>>>>>> 824a45db65c6ca42b85ab44cceacfb8ae75e6c7b
