
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


##############################################################
##### Diversity Indices
##############################################################
##Fixing Time levels:
levels(MetaData$TIME)<- list("W2015"=c("2015-winter","2015-Winter"),"S2016"="2016-Summer",
                             "W2016"="2016-Winter","S2017"="2017-Summer")
#Check and see if it worked?
summary(MetaData)


############## Richness (Species number)#################
# Richness is the nmber of culture observations in the samples. 
# Some samples had more, than one observed species.
hist(isolationRate)

# For the evaluation of richness we need to remove the samples with zero observations.
NotZero = isolationRate > 0 #filter for zero-observation samples

# Keep only the samples with at least one observed species
AbundNotZero=OTUabund[NotZero,]

# Richness in the samples
Richness = specnumber(AbundNotZero)
hist(Richness)
hist(log(Richness))

# Remove the samples with zero observation from the metadata
MetaRich = MetaData[NotZero,]

########### Shannon and Simpson indices####################
# Keep only samples with at least two OTUs
RichNotOne = Richness > 1
AbundNotOne=AbundNotZero[RichNotOne,]

# This keeps observations with at least two OTUs
MetaNotOne = MetaRich[RichNotOne,] 

# Calculate diversity indices
shannon = diversity(AbundNotOne,index = "shannon")
simpson = diversity(AbundNotOne,index = "simpson")

hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

#################################################
####MODELS
#################################################
#Step 1: Write your reserach questions here:

#Step 2: Find the right kind of model (lm or glm) for your data:

# you may need to subset your data based on your questions and goals!!

#Step 3: Choose the best model based on model diagnstics plots and AIC(read your workshop scripts) 

aggregate(OTUabund $ ZEE.se11.Periconia.macrospinosa ~HOST, data = MetaData, sum)
aggregate(OTUabund $ ZEE.se11.Periconia.macrospinosa ~TISSUE, data = MetaData, sum)
aggregate(OTUabund $ ZEE.se11.Periconia.macrospinosa ~SOIL, data = MetaData, sum)
aggregate(OTUabund $ ZEE.se11.Periconia.macrospinosa ~SITE, data = MetaData, sum)
aggregate(OTUabund $ ZEE.se11.Periconia.macrospinosa ~TIME, data = MetaData, sum)
aggregate(OTUabund $ ZEE.se11.Periconia.macrospinosa ~TEMP, data = MetaData, sum)
aggregate(OTUabund $ ZEE.se11.Periconia.macrospinosa ~MEDIA, data = MetaData, sum)


aggregate(OTUabund $ LREwh64..Neocamarosporium.chichastianum ~HOST, data = MetaData, sum)
aggregate(OTUabund $ LREwh64..Neocamarosporium.chichastianum ~TISSUE, data = MetaData, sum)
aggregate(OTUabund $ LREwh64..Neocamarosporium.chichastianum ~SOIL, data = MetaData, sum)
aggregate(OTUabund $ LREwh64..Neocamarosporium.chichastianum ~SITE, data = MetaData, sum)
aggregate(OTUabund $ LREwh64..Neocamarosporium.chichastianum ~TIME, data = MetaData, sum)
aggregate(OTUabund $ LREwh64..Neocamarosporium.chichastianum ~TEMP, data = MetaData, sum)
aggregate(OTUabund $ LREwh64..Neocamarosporium.chichastianum ~MEDIA, data = MetaData, sum)


aggregate(OTUabund $ SREwh19.Neocamarosporium.goegapense ~HOST, data = MetaData, sum)
aggregate(OTUabund $ SREwh19.Neocamarosporium.goegapense ~TISSUE, data = MetaData, sum)
aggregate(OTUabund $ SREwh19.Neocamarosporium.goegapense ~SOIL, data = MetaData, sum)
aggregate(OTUabund $ SREwh19.Neocamarosporium.goegapense ~SITE, data = MetaData, sum)
aggregate(OTUabund $ SREwh19.Neocamarosporium.goegapense ~TIME, data = MetaData, sum)
aggregate(OTUabund $ SREwh19.Neocamarosporium.goegapense ~TEMP, data = MetaData, sum)
aggregate(OTUabund $ SREwh19.Neocamarosporium.goegapense ~MEDIA, data = MetaData, sum)



Article1Meta= subset (MetaData, HOST %in% c("Alhagi persarum", "Artemisia sieberi", "Haloxylon ammodendron", 
                                        "Launaea acunthodes",
                                        "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                        "Tamrix hispida"))


Article1OTU = subset (OTUabund, MetaData$HOST %in% c("Alhagi persarum", "Artemisia sieberi", "Haloxylon ammodendron", 
                                               "Launaea acunthodes",
                                               "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                               "Tamrix hispida"))
# check to see if it worked 
rownames(Article1OTU)==rownames(Article1Meta)
class(Article1Meta)
class(Article1OTU)
View(Article1OTU)
View(Article1Meta)
