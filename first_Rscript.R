##################################
## Pakages we need for the analysis
library(mvabund)
library(vegan)
library(coda)
library(jags)
library(rjags)
library(boral)
library(effects)
library(MASS)
library(ggplot2)

##################################
## DATA Input

data <- read.csv("MSHM.csv",header = T, row.names = 1)
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


##############################################
##### Diversity Indices
##############################################
########Fixing Time levels:
levels(MetaData$TIME)<- list("W2015"=c("2015-winter","2015-Winter"),"S2016"="2016-Summer",
                             "W2016"="2016-Winter","S2017"="2017-Summer")
#Check and see if it worked?
summary(MetaData)

######## Richness (Species number)#################
######## Richness is the nmber of culture observations in the samples. 
######## Some samples had more than one observed species.
hist(isolationRate)

########  For the evaluation of richness we need to remove the samples with zero observations.
NotZero = isolationRate > 0 
########filter for zero-observation samples
######## Keep only the samples with at least one observed species
AbundNotZero=OTUabund[NotZero,]

######## Richness in the samples
Richness = specnumber(AbundNotZero)
hist(Richness)
hist(log(Richness))

######## Remove the samples with zero observation from the metadata
MetaRich = MetaData[NotZero,]

######## Shannon and Simpson indices####################
######## Keep only samples with at least two OTUs
RichNotOne = Richness > 1
AbundNotOne=AbundNotZero[RichNotOne,]

######## This keeps observations with at least two OTUs
MetaNotOne = MetaRich[RichNotOne,] 

######## Calculate diversity indices
shannon = diversity(AbundNotOne,index = "shannon")
simpson = diversity(AbundNotOne,index = "simpson")
hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

######## Aggregate of 3 Fungi(Periconia macrospinosa,Neocamarosporium chichastianum and Neocamarosporium.goegapense)
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

#################################################
#### Article 1
#################################################

########Step 1: Write your reserach questions for the first article here:

# What is the effect of saline and arid soil on root endophyte communities in contrast 
#to the type of host plant?
# The effects of climate parameters on diversity of root endophytes in different soils?
# Do dominant endophytes in saline and arid soils can induce salinity and 
#drought resistance in model plants?

########Step 2: We can do these analyses

# Perecentage of root pieces yielding 
# Yield root endophytic fungal growth
# Overal colonization perecentage 
# An averaged colonization per population 
# OTU richness, frequency, diversity and Community structure 
#(VEGAN v2.2-1 (Oksanenet al., 2015)) across every Hosts, soils and climate parameters
# Investigating general patterns of variation in the endophytic fungi 
#impact on model plants growth across abiotic conditionsfactor levels using 
#the Kruskal–Wallis rank sum test
# Taxonomic classification of isolates

  
####### Subsetting the data for ARTICLE1
Article1M = subset (MetaData, HOST%in%c("Alhagi persarum","Artemisia sieberi", "Haloxylon ammodendron", 
                                            "Launaea acunthodes",
                                            "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                            "Tamrix hispida"))
# some how this still shows the 40 hists!! I am trying another way to fix it!
levels(Article1M$SITE)
levels(Article1M$HOST)
write.csv(Article1M, file = "testmetadata.csv")

Article1Meta<-read.csv("testmetadata.csv",header = T, row.names = 1)

# remane the long variable levels
levels(Article1Meta$SITE)
levels(Article1Meta$HOST)
levels(Article1Meta$HOST)<- list("A.pers"="Alhagi persarum","A.sieb"="Artemisia sieberi",
                                 "H.ammo"="Haloxylon ammodendron","L.acun"="Launaea acunthodes",
                                 "P.step"="Prosopis stephaniana","S.inca"="Salsola incanescens",
                                 "S.rosm"="Seidlitzia rosmarinus","T.hisp"="Tamrix hispida")

levels(Article1Meta$SOIL)<-list("Arid"="Arid soil","Saline"="Saline Soil")
levels(Article1Meta$TISSUE)
levels(Article1Meta$MEDIA)<-list("PDA"= "PDA", "PDA+Plant"="PDA+Plant extract")
## Subset OTU frequency dataframe
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
colnames(Article1OTU)
#FROM NOW ON: WE ONLY USE THESE TWO DATA FRAMES FOR ARTICLE 1: Article1OTU & Article1Meta

# What OTUs are in this project?

OTUsINarticl1<-colSums(Article1OTU)
Article1OTU<-Article1OTU[, colSums(Article1OTU != 0) > 0]
# OTU frequencies Article 1
colSums(Article1OTU)
#OTU list Article 1
colnames(Article1OTU)

### chosen OTU for other experiments:LREwh64..Neocamarosporium.chichastianum
# isolated from?
aggregate(LREwh64..Neocamarosporium.chichastianum ~ Article1Meta$SOIL, Article1OTU, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum ~ Article1Meta$HOST, Article1OTU, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum ~ Article1Meta$TISSUE, Article1OTU, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum ~ Article1Meta$TIME, Article1OTU, sum)


#############################################################################
######## Step 3: Find the right kind of analysis for each research questions:

#########################################################
#########################################################
##### Diversity indices
#For diversity we are using model based aproaches:
## I am using the codes that we wrote earlier:
## first delete the IR in metadata.not sure if it is ok
Article1Meta$IR<-NULL

IR.art1= apply(Article1OTU,1, sum)
Article1Meta = cbind(Article1Meta, IR = IR.art1)

######## Richness (Species number)#################

hist(IR.art1)

########  For the evaluation of richness we need to remove the samples with zero observations.
NotZero.art1= IR.art1 > 0 
########filter for zero-observation samples
######## Keep only the samples with at least one observed species
AbundNotZero.art1=Article1OTU[NotZero.art1,]

######## Richness in the samples
Richness.art1 = specnumber(AbundNotZero.art1)
hist(Richness.art1)
hist(log(Richness.art1))

######## Remove the samples with zero observation from the metadata
MetaRich.ART1 = Article1Meta[NotZero.art1,]
row.names(AbundNotZero.art1)==row.names(MetaRich.ART1)
######## Shannon and Simpson indices####################
######## Keep only samples with at least two OTUs
RichNotOne.art1 = Richness.art1 > 1
AbundNotOne.art1=AbundNotZero.art1[RichNotOne.art1,]

######## This keeps observations with at least two OTUs
MetaNotOne.art1 = MetaRich.ART1[RichNotOne.art1,] 

######## Calculate diversity indices
shannon.art1 = diversity(AbundNotOne.art1,index = "shannon")
simpson.art1 = diversity(AbundNotOne.art1,index = "simpson")
hist(shannon.art1)
hist(simpson.art1)
hist(log(shannon.art1))
hist(log(simpson.art1))

#Diversity index objects: shannon.art1 & simpson.art1 & Richness.art1
###################
# Richness model
###################
hist(Richness.art1)
hist(sqrt(Richness.art1))
hist(log10(Richness.art1))
mean(Richness.art1)
var(Richness.art1)
# Since variance is lower than mean: We use poisson model instead of negbin
# poisson model
# r.m<-glm(formula =Richness.art1~SOIL*TISSUE+ SOIL*HOST+ TIME+ SITE,data = MetaRich.ART1,
#          family=poisson(link = "log"))
# summary(r.m) ##NA??
# AIC(r.m)
# par(mfrow=c(2,2))
# plot(r.m)
# ######
# r.m1<-glm(formula =Richness.art1~SOIL*TISSUE*HOST+ TIME+ SITE,data = MetaRich.ART1,
#          family=poisson(link = "log"))
# summary(r.m1)
# plot(r.m1)
# AIC(r.m1)
# #####
# r.m2<-glm(formula =Richness.art1~SOIL*TISSUE+HOST+ TIME+ SITE,data = MetaRich.ART1,
#           family=poisson(link = "log"))
# summary(r.m2)
# plot(r.m2)
# AIC(r.m2)
# ####
# r.m3<-glm(formula =Richness.art1~SOIL+TISSUE+HOST+ TIME+ SITE,data = MetaRich.ART1,
#           family=poisson(link = "log"))
# summary(r.m3)
# plot(r.m3)
# AIC(r.m3)
# #########
# r.m4<-glm(formula =Richness.art1~SOIL*HOST+TISSUE+ TIME+ SITE,data = MetaRich.ART1,
#           family=poisson(link = "log"))
# summary(r.m4)
# plot(r.m4)
# AIC(r.m4)
# # # negbin models
# r.m5<-glm(formula =Richness.art1~SOIL*TISSUE+HOST+ TIME+ SITE,data = MetaRich.ART1,
#           family=negative.binomial(theta = 0.5))
# summary(r.m5)
# plot(r.m5)
# AIC(r.m5)
# dev.off()
# # r.m6<-glm(formula =Richness.art1~SOIL*TISSUE+HOST+ TIME+ SITE,data = MetaRich.ART1,
#           family=negative.binomial(theta = 0.2))
# summary(r.m6)
# plot(r.m6)
# AIC(r.m6)
# r.m7<-glm(formula =Richness.art1~SOIL*TISSUE+HOST+ TIME+ SITE,data = MetaRich.ART1,
#           family=negative.binomial(theta =10000 ))
# summary(r.m7)
# anova(r.m7, test = "Chisq")
# plot(r.m7)
# AIC(r.m7)
# r.m.nb<-glm.nb(formula =Richness.art1~SOIL*TISSUE+ HOST+ TIME+
#                  SITE,data = MetaRich.ART1,link = "log")
# summary(r.m.nb)
# anova(r.m.nb,test = "Chisq")
# AIC(r.m.nb)
# 
# plot(r.m.nb)
# ## model comparison with anova test
# aa.m<-anova(r.m, r.m1, r.m2,r.m3,r.m4,r.m5,r.m6,r.m7,r.m.nb, test = "Chisq")
# summary(aa.m)

### stepwise model selection with all variables
R.m<-glm(formula =Richness.art1~SOIL+HOST+TISSUE+TIME+SITE+MEDIA+TEMP ,data = MetaRich.ART1,
        family=poisson(link = "log"))

stepAIC(R.m,direction="backward")
# try couple of interactions
# R.m1<-glm(formula =Richness.art1~SOIL*HOST+TISSUE+TIME+SITE,data = MetaRich.ART1,
#          family=poisson(link = "log"))
# summary(R.m1)
# anova(R.m1, test = "Chisq")
# AIC(R.m1)
R.m2<-glm(formula =Richness.art1~SOIL*HOST+TISSUE+TIME+SITE,data = MetaRich.ART1,
          family=poisson(link = "log"))
summary(R.m2)
anova(R.m2, test = "Chisq")
AIC(R.m2)
################################
# SELECTED MODEL FOR RICHNESS
################################
Richness.m<-glm(formula =Richness.art1~SOIL+HOST+TISSUE+TIME+SITE,data = MetaRich.ART1,
                family=poisson(link = "log"))
#MODEL SUMMARY FOR REACHNESS
Rich.summ<-summary(Richness.m)
#ANOVA RESULTS FOR RICHNESS
Rich.anova<-anova(Richness.m, test = "Chisq")
AIC(Richness.m)
par(mfrow=(c(2,2)))
plot(Richness.m)
dev.off()

################
# SHANNON MODEL
################

sh.m<-lm(formula =log(shannon.art1)~SOIL*TISSUE+HOST+TIME+SITE,data = MetaNotOne.art1)
# sh.m2<-lm(formula =log(shannon.art1)~SOIL*HOST+TISSUE+TIME+SITE,data = MetaNotOne.art1)
# sh.m3<-lm(formula =log(shannon.art1)~HOST*SOIL+TISSUE+TIME+SITE,data = MetaNotOne.art1)

#MODEL SUMMARY FOR SHANNON
shannon.sum<-summary(sh.m)
# shannon.sum2<-summary(sh.m2)
# shannon.sum3<-summary(sh.m3)

#ANOVA RESULS FOR SHANON
shannon.anov<-anova(sh.m, test="F")
AIC(sh.m)
# shannon.anov2<-anova(sh.m2, test="F")
# AIC(sh.m2)
# 
# shannon.anov2<-anova(sh.m3, test="F")
# AIC(sh.m3)
###############
#SIMPSON MODEL
###############
simp.m<-lm(formula =log(simpson.art1)~SOIL*TISSUE+HOST+TIME+SITE,data = MetaNotOne.art1)
# simp.m2<-lm(formula =log(simpson.art1)~SOIL*HOST+TISSUE+TIME+SITE,data = MetaNotOne.art1)
# simp.m3<-lm(formula =log(simpson.art1)~HOST*SOIL+TISSUE+TIME+SITE,data = MetaNotOne.art1)

#MODEL SUMMARY FOR SIMPSON
simp.sum<-summary(simp.m)
# simp.sum<-summary(simp.m2)
# simp.sum<-summary(simp.m3)
#ANOVA RESULTS FOR SIMPSON
simp.anov<-anova(simp.m, test="F")
AIC(simp.m)
# simp.anov<-anova(simp.m2, test="F")
# AIC(simp.m2)
# simp.anov<-anova(simp.m3, test="F")
# AIC(simp.m3)
######## Step 4: Choose the best way to visualize the results



####################################################
####################################################
# COMMUNITY COMPOSITION ANALYSIS
####################################################
####################################################

# for this we are using Permutational multivariate analysis of variance (PERMANOVA)
#this function from vegan does it
memory.limit()
memory.limit(size=30000)
# data transformation:
tran.abund.notzero<-decostand(AbundNotZero.art1,method="hellinger")
#PERMANOVA
comm.anova<-adonis(formula=tran.abund.notzero~SOIL*TISSUE+HOST+TIME+SITE, data= MetaRich.ART1,
       permutations = 999, method = "bray",by=NULL)
# if you get an error about the memory allocation run the above lines:memory.limit


## Coefficients
comm.anova.coef = as.data.frame(comm.anova$coefficients)

## which OTUs significantly affected by any variables:



#####################
####### NMDS PLOTS 
#####################
# first I am computing an NMDS matrix for all of the OTUs
# nmds.art1<-metaMDS(AbundNotZero.art1, distance = "bray", k= 2, trymax = 20)
#NMDS with square root transformed data:
nmds.art2<-metaMDS(tran.abund.notzero, distance = "bray", k= 2, trymax = 20)
plot(nmds.art1)
# red + are species and black dots are sites (samples)

###########################################
#this is the base of our plots: 
# NM.pl<-ordiplot(nmds.art2,type = "none")
# you can show species ot samples or both
# points(NM.pl,"sites",pch = 19, cex= 0.5, col="grey30")
# p.spe<-points(NM.pl,"species",pch = 2, col= "black", cex= 0.8)

###########################################
# NMDS soil:
###########################################
dev.off()
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5))
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.8)

ordiellipse(nmds.art2, MetaRich.ART1$SOIL,cex=1,lwd = 3,alpha = 200, 
            draw="polygon", col= "Blue",border="Blue",
             kind="se", conf=0.95,show.groups=(c("Saline")))

ordiellipse(nmds.art2, MetaRich.ART1$SOIL,cex=1,alpha = 200, 
            draw="polygon", col= "red",lwd = 3,border="red",
             kind="se", conf=0.95,show.groups=(c("Arid")))
legend("bottomright", c("Saline soil","Arid soil"), 
       fill=c("blue","red"),
       border="white", bty="n")

###########################################
# NMDS Hosts
###########################################
levels(MetaRich.ART1$HOST)
dev.off()
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5))
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.8)
ordiellipse(nmds.art2, MetaRich.ART1$HOST,cex=1,alpha = 200, 
            draw="polygon", col= 1:8,border=1:8,lwd=3, kind="se", conf=0.95)

legend("bottomright", c("A.pers" ,"A.sieb", "H.ammo" ,"L.acun" ,"P.step" ,"S.inca", "S.rosm" ,"T.hisp"), 
       fill= 1:8,border="white", bty="n")
###########################################
#NMDS sampling sites:
###########################################
levels(MetaRich.ART1$SITE)
dev.off()
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5))
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.8)
ordiellipse(nmds.art2, MetaRich.ART1$SITE,cex=1,alpha = 200, 
            draw="polygon", col= 1:5,border= 1:5,lwd=3, kind="se", conf=0.95)

legend("bottomright", c("Garmsar","Haj Ali Gholi Lake","Hoze Soltan Lake","Maranjab Desert","Rig-Boland Desert"), 
       fill= 1:5, border="white", bty="n")
###########################################
#NMDS ORGAN:
###########################################
levels(MetaRich.ART1$TISSUE)
dev.off()
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5))
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.8)
ordiellipse(nmds.art2, MetaRich.ART1$TISSUE,cex=1,
            draw="polygon", col= 1:3,border=1:3,lwd=3,alpha = 200, kind="se", conf=0.95)

legend("bottomright", c("Branch","Leaf","Root" ), 
       fill= 1:3,border="white", bty="n")


#################
#plot organ and soil toghether:
dev.off()
NM.pl<-ordiplot(nmds.art1,type = "none")
# points(NM.pl,"sites",pch = 19, cex= 0.5, col="grey30")
points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.6)
# show soil communities:
ordihull(nmds.art1, MetaRich.ART1$SOIL,cex=1.5,
         draw="line", col= "grey20",
         lwd = 2, lty = 2,
         show.groups=(c("Arid")))#Arid soil endophyte community
ordihull(nmds.art1, MetaRich.ART1$SOIL,cex=1.5,
         draw="line", col= "grey20",
         lwd = 2,lty = 9,
         show.groups=(c("Saline")))# Saline soil endophyte community
# show different organs: leaf, twig and root:

ordiellipse(nmds.art1, MetaRich.ART1$TISSUE,cex=1,
            draw="polygon", col= "green",
            alpha=200, kind="se", conf=0.95,
            show.groups=(c("Leaf")))
ordiellipse(nmds.art1, MetaRich.ART1$TISSUE,cex=1.5,
            draw="polygon", col= "blue",
            alpha=200, kind="se", conf=0.95,
            show.groups=(c("Branch")))
ordiellipse(nmds.art1, MetaRich.ART1$TISSUE,cex=1,
            draw="polygon", col= "red",
            alpha=200, kind="se", conf=0.95,
            show.groups=(c("Root")))
legend("topleft", c("Arid soil","Saline soil","Leaf", "Branch",
                    "Root"), 
       col=c("grey20","grey20","green","blue","red"),
       lty = c(2,9,1,1,1), border="white", bty="n")

########################################################
#### Individul reactions of OTUs to our variables
########################################################
####use glm.nb for each OTU to figur out each of their frequency recation to every variables

OTU1.model<-glm.nb(formula =AbundNotZero.art1$APE.se5.Staphylotrichum.coccosporum~SOIL*TISSUE+ HOST+ TIME+
                           SITE,data = MetaRich.ART1,link = "log")
OTU1.anov<-anova(OTU1.model,test = "Chisq")

OTU2.model<-glm.nb(formula =AbundNotZero.art1$TPEsh28.Humicola.fuscoatra~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU2.anov<-anova(OTU2.model,test = "Chisq")

OTU3.model<-glm.nb(formula =AbundNotZero.art1$PFE.sh7..Rosellinia.limonispora~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU3.anov<-anova(OTU3.model,test = "Chisq")

OTU4.model<-glm.nb(formula =AbundNotZero.art1$LREwh64..Neocamarosporium.chichastianum~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU4.anov<-anova(OTU4.model,test = "Chisq")

OTU5.model<-glm.nb(formula =AbundNotZero.art1$SREwh22.Preussia.minimoides~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU5.anov<-anova(OTU5.model,test = "Chisq")

OTU6.model<-glm.nb(formula =AbundNotZero.art1$THE.we10..Aporospora.terricola~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU6.anov<-anova(OTU6.model,test = "Chisq")

OTU7.model<-glm.nb(formula =AbundNotZero.art1$THE.ss3.Fusarium.sp.~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU7.anov<-anova(OTU7.model,test = "Chisq")

OTU8.model<-glm.nb(formula =AbundNotZero.art1$SIE.sh1.Briansuttonomyces.eucalypti~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU8.anov<-anova(OTU8.model,test = "Chisq")

OTU9.model<-glm.nb(formula =AbundNotZero.art1$RAE.sh12.Acrocalymma.vagum~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU9.anov<-anova(OTU9.model,test = "Chisq")

OTU10.model<-glm.nb(formula =AbundNotZero.art1$PSE.wh14.Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = MetaRich.ART1,link = "log")
OTU10.anov<-anova(OTU10.model,test = "Chisq")

OTU11.model<-glm.nb(formula =AbundNotZero.art1$PSE.wh40.Dimorphosporicola.tragani~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU11.anov<-anova(OTU11.model,test = "Chisq")

OTU12.model<-glm.nb(formula =AbundNotZero.art1$PSE.wh66.Comoclathris.italica~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU12.anov<-anova(OTU12.model,test = "Chisq")

OTU13.model<-glm.nb(formula =AbundNotZero.art1$PSE.ss7.Penicillium.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU13.anov<-anova(OTU13.model,test = "Chisq")

OTU14.model<-glm.nb(formula =AbundNotZero.art1$PSE.we4..Chaetomium.globosum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU14.anov<-anova(OTU14.model,test = "Chisq")

OTU15.model<-glm.nb(formula =AbundNotZero.art1$PSE.se8.Podospora.minicauda~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU15.anov<-anova(OTU15.model,test = "Chisq")

OTU16.model<-glm.nb(formula =AbundNotZero.art1$PSE.we8.Alternaria.chlamydospora~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU16.anov<-anova(OTU16.model,test = "Chisq")

OTU17.model<-glm.nb(formula =AbundNotZero.art1$LDE.se7.Coniothyrium.aleuritis~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU17.anov<-anova(OTU17.model,test = "Chisq")

OTU18.model<-glm.nb(formula =AbundNotZero.art1$LAE.se5.Sordaria.humana~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU18.anov<-anova(OTU18.model,test = "Chisq")

OTU19.model<-glm.nb(formula =AbundNotZero.art1$HAE.se5.Camarosporomyces.flavigenus~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU19.anov<-anova(OTU19.model,test = "Chisq")

OTU20.model<-glm.nb(formula =AbundNotZero.art1$HAE.ss4.Cytospora.chrysosperma~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU20.anov<-anova(OTU20.model,test = "Chisq")

OTU21.model<-glm.nb(formula =AbundNotZero.art1$HAE.we5.Coniolariella.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU21.anov<-anova(OTU21.model,test = "Chisq")

OTU22.model<-glm.nb(formula =AbundNotZero.art1$HAE.wh26.Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU22.anov<-anova(OTU22.model,test = "Chisq")

OTU23.model<-glm.nb(formula =AbundNotZero.art1$HAE.wh10.Raffaelea.montetyi~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU23.anov<-anova(OTU23.model,test = "Chisq")

OTU24.model<-glm.nb(formula =AbundNotZero.art1$HAE.wh65.Coniophora.marmorata~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU24.anov<-anova(OTU24.model,test = "Chisq")

OTU25.model<-glm.nb(formula =AbundNotZero.art1$HAE.se9.Chaetomium.nigricolor~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU25.anov<-anova(OTU25.model,test = "Chisq")

OTU26.model<-glm.nb(formula =AbundNotZero.art1$HAE.se1.Acrocalymma.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU26.anov<-anova(OTU26.model,test = "Chisq")

OTU27.model<-glm.nb(formula =AbundNotZero.art1$SREwh19.Neocamarosporium.goegapense~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU27.anov<-anova(OTU27.model,test = "Chisq")

OTU28.model<-glm.nb(formula =AbundNotZero.art1$APEsh6.Dictyosporium.digitatum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU28.anov<-anova(OTU28.model,test = "Chisq")

OTU29.model<-glm.nb(formula =AbundNotZero.art1$APE.sh8.Pestalotiopsis.vismiae~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU29.anov<-anova(OTU29.model,test = "Chisq")

OTU30.model<-glm.nb(formula =AbundNotZero.art1$APE.se3.Dactylonectria.macrodidyma~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU30.anov<-anova(OTU30.model,test = "Chisq")

OTU31.model<-glm.nb(formula =AbundNotZero.art1$APE.sh5.Nigrospora.sphaerica~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU31.anov<-anova(OTU31.model,test = "Chisq")

OTU32.model<-glm.nb(formula =AbundNotZero.art1$SREwh18...Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU32.anov<-anova(OTU32.model,test = "Chisq")

OTU33.model<-glm.nb(formula =AbundNotZero.art1$LAEsh5.Coniolariella.ershadii~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1 = "log")
OTU33.anov<-anova(OTU33.model,test = "Chisq")

OTU34.model<-glm(formula =AbundNotZero.art1$LAE.se3.Neosetophoma.lunariae~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU34.anov<-anova(OTU34.model,test = "Chisq")

OTU35.model<-glm.nb(formula =AbundNotZero.art1$LAE.SH.7.Muriphaeosphaeria.viburni~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU35.anov<-anova(OTU35.model,test = "Chisq")

OTU36.model<-glm.nb(formula =AbundNotZero.art1$LAE.sh1.Acrocalymma.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU36.anov<-anova(OTU36.model,test = "Chisq")

OTU37.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh30..Preussia.grandispora~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU37.anov<-anova(OTU37.model,test = "Chisq")

OTU38.model<-glm.nb(formula =AbundNotZero.art1$SRE.ss.4.Neocamarosporium.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU38.anov<-anova(OTU38.model,test = "Chisq")

OTU39.model<-glm.nb(formula =AbundNotZero.art1$SRE.ws8.Botryotrichum.murorum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU39.anov<-anova(OTU39.model,test = "Chisq")

OTU40.model<-glm.nb(formula =AbundNotZero.art1$SRE.ws10.Sarocladium.kiliense~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU40.anov<-anova(OTU40.model,test = "Chisq")

OTU41.model<-glm.nb(formula =AbundNotZero.art1$SRE.wh16.Paracamarosporium.hawaiiense~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data =  MetaRich.ART1,link = "log")
OTU41.anov<-anova(OTU41.model,test = "Chisq")

OTU42.model<-glm.nb(formula =AbundNotZero.art1$SRE.wh13.Ovatospora.senegalensis~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU42.anov<-anova(OTU42.model,test = "Chisq")

OTU43.model<-glm.nb(formula =AbundNotZero.art1$SRE.we6.Fusariella.sinensis~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU43.anov<-anova(OTU43.model,test = "Chisq")

OTU44.model<-glm.nb(formula =AbundNotZero.art1$SRE.we.10.pichia.kudriavzevii~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU44.anov<-anova(OTU44.model,test = "Chisq")

OTU45.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh7.Chaetomium.cucumericola~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU45.anov<-anova(OTU45.model,test = "Chisq")

OTU46.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh5.Fusarium.redolens~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU46.anov<-anova(OTU46.model,test = "Chisq")

OTU47.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh9.Preussia.intermedia~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU47.anov<-anova(OTU47.model,test = "Chisq")

OTU48.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh4.Penicillium.vinaceum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU48.anov<-anova(OTU48.model,test = "Chisq")

OTU49.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh3.Trichoderma.rifaii~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1,link = "log")
OTU49.anov<-anova(OTU49.model,test = "Chisq")

### MAHDIEH!!!!
# what are these numbers 
########## SOIL 
# 1.3.4.5.6.7.10.11.12.13.14.15.16.17.18.19.20.21.22.23.24.25.26.27.28.29.30.31.
# 32.33.34.35.36.37.38.39.4043.44.45.46.47.48.49

########## Tissue
# 1.2.3.4.5.6.7.8.9.10.11.12.13.14.16.17.18.19.20.21.22.24.25.26.27.28.29.30.31.32.33.34.35.36.37.39.40.41.42.43.
# 45.46.47.48.49
########## HOST
# 1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21.22.23.24.25.26.27.28.29.30.31.32.33.34.35.36.37.38.39.40.
# 41.42.43.44.45.46.47.48.49
# 
########## TIME
# 2.5.9.11.12.13.14.15.1618.20.21.22.23.24.25.27.32.33.34.37.38.39.40.41.44.45.46.47.48.49

########## SITE
#2.3.4.13.14.1520.2125.33.38.39.40.43.44.45.46.47.48.49

########## SITE:Tissue
#2.8.

##############################################
##############################################
# Variation partitioning
##############################################
##############################################

?varpart
var.art1<-varpart(Article1OTU,Article1Meta$SOIL,
                 Article1Meta$HOST, Article1Meta$SITE,Article1Meta$TISSUE)
plot(var.art1)

# >0 samles only
var.art2<-varpart(AbundNotZero.art1,MetaRich.ART1$SOIL,
                  MetaRich.ART1$HOST, MetaRich.ART1$SITE,MetaRich.ART1$TISSUE, 
                  transfo="hellinger")
plot(var.art2)
#without tissue
var.art3<-varpart(AbundNotZero.art1,MetaRich.ART1$TIME,
                  MetaRich.ART1$HOST, MetaRich.ART1$SITE,MetaRich.ART1$SOIL, 
                  transfo="hellinger")
plot(var.art3)
##varpart with models
var.part4<- varpart(AbundNotZero.art1,~MetaRich.ART1$TIME+
                    MetaRich.ART1$HOST+ MetaRich.ART1$SITE+MetaRich.ART1$TISSUE,
                   ~env.Rich.ART1$Ece+env.Rich.ART1$pHe+env.Rich.ART1$Cle+env.Rich.ART1$EC,
                   transfo = "hellinger")
plot(var.part4)

#with soil*tissue
var.part5<-varpart(AbundNotZero.art1, MetaRich.ART1$HOST,MetaRich.ART1$TIME,
                  MetaRich.ART1$SITE,~MetaRich.ART1$SOIL*MetaRich.ART1$TISSUE, 
                  transfo="hellinger")
plot(var.part5,Xnames = c("Host","Time","Site","Soil*Organ"),
     bg=c("grey20","green","blue","red"),alpha=100)

##############################################
##############################################
#soil data analysis
##############################################
##############################################
# data input
env.data<-read.csv("MataDataMergSoil.csv", header = T, row.names = 1)

#subset for this article:
env.data.art1<- subset (env.data, MetaData$HOST%in%c("Alhagi persarum","Artemisia sieberi", "Haloxylon ammodendron", 
                                           "Launaea acunthodes",
                                           "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                           "Tamrix hispida"))
row.names(env.data.art1)==row.names(Article1Meta)

# remove thesamples with zero obs for PCA
env.Rich.ART1 = env.data.art1[NotZero.art1,]
row.names(env.Rich.ART1)==row.names(MetaRich.ART1)
#remove soil and site variables
env.Rich.ART1$SITE<-NULL
env.Rich.ART1$SOIL<-NULL
env.Rich.ART1$Cle<-as.numeric(env.Rich.ART1$Cle)
##### PCA (Principal Components Analysis) analysis 
#use rda fun in vegan for PCA:
?rda()
rda.art1<-rda(AbundNotZero.art1~.,data =env.Rich.ART1)
head(summary(rda.art1))

#Plot PCA results
plot(rda.art1, xlim=c(-2,2), ylim=c(-1,2.5))# this is a simple plot with both sites and species

# #coloring:
# soilfactor<-factor(MetaRich.ART1$SOIL)
# colvec <-  c("red","blue")
# cols<-colvec[soilfactor]
# #scors:
rda.scores <- scores(rda.art1, display = 'bp')
mul <- ordiArrowMul(rda.scores, fill = 0.75)

dev.off()
####### FINAL RDA PLOT
plot(rda.art1, type = "n")
points(rda.art1, display = "sites", col = colvec[soilfactor], 
       pch = (16:17)[soilfactor],cex=0.85)#salin=blue
#points(rda.art1, display = "species", pch = "+")
arrows(0, 0, mul * rda.scores[,1], mul * rda.scores[,2],
       length = 0.05, col ="black", lwd=2)
labs <- rownames(rda.scores)
#labs<-c("Ece","EC","Cle","pHe")
text(ordiArrowTextXY(mul * rda.scores, labs), labs)
legend("topleft", c("Arid soil","Saline soil"), 
             col=c("red","blue"),
       pch = c(16,17), border="white", bty="n")

#### SOIL DATA VARATION PARTITIONING

var.soil1<-varpart(AbundNotZero.art1, env.Rich.ART1$Ece+env.Rich.ART1$EC+env.Rich.ART1$Cle+env.Rich.ART1$pHe ,
                   env.Rich.ART1$pH.1.5+env.Rich.ART1$OM+env.Rich.ART1$OC+env.Rich.ART1$Nt+env.Rich.ART1$P+
                     env.Rich.ART1$P2O5+env.Rich.ART1$M+env.Rich.ART1$SP+env.Rich.ART1$M.SP,
                   transfo="hellinger")
plot(var.soil1)
#This is not good for us rather ignore it

##############################################
##############################################
# Greenhouse data analysis
##############################################
##############################################
# data input
GH.data<-read.csv("Green house data.csv", header = T, row.names = 1)
View(GH.data)

GH.data$Biomass<-GH.data$DWshoot+GH.data$DWroot

# factors: Fungi, Drought, Salinity
# ?aov
# GH.1<-aov(Lshoot~Fungi*Drought*Salinity,data = GH.data)
# summary(GH.1)
hist(GH.data$Lshoot)
hist(GH.data$Wshoot)
hist(GH.data$DWshoot)
hist(GH.data$Lroot)
hist(GH.data$Wroot)
hist(GH.data$DWroot)
hist(GH.data$Photosyntesis)
hist(GH.data$Biomass)


boxplot(GH.data$Lshoot)
boxplot(GH.data$Wshoot)
boxplot(GH.data$DWshoot)
boxplot(GH.data$Lroot)
boxplot(GH.data$Wroot)
boxplot(GH.data$DWroot)
boxplot(GH.data$Photosyntesis)
boxplot(GH.data$Biomass)


#outlier fixed
dev.off()
# try GLM
glm.Lshoot<-glm(Lshoot~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.Lshoot)
anova(glm.Lshoot, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.Lshoot)

glm.Wshoot<-glm(Wshoot~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.Wshoot)
anova(glm.Wshoot, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.Wshoot)

glm.DWshoot<-glm(DWshoot~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.DWshoot)
anova(glm.DWshoot, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.DWshoot)

glm.Lroot<-glm(Lroot~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.Lroot)
anova(glm.Lroot, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.Lroot)

glm.Wroot<-glm(Wroot~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.Wroot)
anova(glm.Wroot, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.Wroot)

glm.DWroot<-glm(DWroot~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.DWroot)
anova(glm.DWroot, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.DWroot)

glm.Photosyntesis<-glm(Photosyntesis~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.Photosyntesis)
anova(glm.Photosyntesis, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.Photosyntesis)

glm.Biomass<-glm(Biomass~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.Biomass)
anova(glm.Biomass, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.Biomass)

# how to show this?

ggplot(GH.data, aes(x=Drought,y=Lshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()

ggplot(GH.data, aes(x=Drought,y=Wshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()

ggplot(GH.data, aes(x=Drought,y=DWshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()

ggplot(GH.data, aes(x=Drought,y=Lroot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()

ggplot(GH.data, aes(x=Drought,y=Wroot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()

ggplot(GH.data, aes(x=Drought,y=DWroot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()

ggplot(GH.data, aes(x=Drought,y=Photosyntesis, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()


ggplot(GH.data, aes(x=Drought,y=Biomass, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()



##############################################
##############################################
# Antifungal test
##############################################
##############################################
# data input
#########
Antifungal.data<-read.csv("Antifungal.csv", header = T, row.names = 1)
# View(Antifungal.data)
# 
# boxplot(Growth ~ Pathogen *Fungi , data=Antifungal.data)
# 
# interaction.plot(x.factor = Antifungal.data$Fungi,
#                  trace.factor = Antifungal.data$Pathogen,
#                  response = Antifungal.data$Growth)

# P. oryzea
Antifungal.PO.data<-read.csv("Antifungal.PO.csv", header = T, row.names = 1)
View(Antifungal.PO.data)

t.test(Growth ~ Fungi,data = Antifungal.PO.data)
ggplot(Antifungal.PO.data, aes(x = Fungi, y = Growth)) + 
  geom_boxplot() 

# par(mfrow = c(2, 3))
# plot(Antifungal.PO.data)
# hist(Antifungal.PO.data$Growth)
# boxplot(Antifungal.PO.data$Growth)
# plot(Growth ~ Fungi, data=Antifungal.PO.data)
# Antifungal.PO.lm <- lm(Growth ~ Fungi, data=Antifungal.PO.data)

# A.conoides
Antifungal.AC.data<-read.csv("Antifungal.AC.csv", header = T, row.names = 1)
View(Antifungal.AC.data)

t.test(Growth ~ Fungi,data = Antifungal.AC.data)
ggplot(Antifungal.AC.data, aes(x = Fungi, y = Growth)) + 
  geom_boxplot() 

# par(mfrow = c(2, 3))
# plot(Antifungal.AC.data)
# hist(Antifungal.AC.data$Growth)
# boxplot(Antifungal.AC.data$Growth)

# P. graminea
Antifungal.PG.data<-read.csv("Antifungal.PG.csv", header = T, row.names = 1)
View(Antifungal.PG.data)

t.test(Growth ~ Fungi,data = Antifungal.PG.data)
ggplot(Antifungal.PG.data, aes(x = Fungi, y = Growth)) + 
  geom_boxplot() 

# par(mfrow = c(2, 3))
# plot(Antifungal.PG.data)
# hist(Antifungal.PG.data$Growth)
# boxplot(Antifungal.PG.data$Growth)

## Order DATA Input

# Order <- read.csv(file="OrderArticle.csv",header = T, row.names = 1)
Order <- read.csv(file="OrderArticle.csv",header = T, row.names = 1)
str(Order)
summary(Order)


######## MERG the ORDER abundance data by 4
S1<-seq(1,121920,4)
S2<-seq(4,121920,4)
Order1<-matrix(0,length(S1),11)
for (i in 1:length(S1)) {
  Order1[i,]<-colSums(Order[S1[i]:S2[i],])}

########now convert to data frame and rename the columns and rows
orderabund<-data.frame(Order1)
class(orderabund)
colnames(orderabund)=colnames(Order)
name<- row.names(data)
row.names(orderabund)<- name[1:30480]
summary(orderabund)

####### Subsetting the Order abundance data for ARTICLE1(keep only 8 hosts)
####### Subsetting the data for ARTICLE1
Article1M = subset (MetaData, HOST%in%c("Alhagi persarum","Artemisia sieberi", "Haloxylon ammodendron", 
                                        "Launaea acunthodes",
                                        "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                        "Tamrix hispida"))
# some how this still shows the 40 hists!! I am trying another way to fix it!
levels(Article1M$SITE)
levels(Article1M$HOST)
write.csv(Article1M, file = "testmetadata.csv")

Article1Meta<-read.csv("testmetadata.csv",header = T, row.names = 1)

# remane the long variable levels
levels(Article1Meta$SITE)
levels(Article1Meta$HOST)
levels(Article1Meta$HOST)<- list("A.pers"="Alhagi persarum","A.sieb"="Artemisia sieberi",
                                 "H.ammo"="Haloxylon ammodendron","L.acun"="Launaea acunthodes",
                                 "P.step"="Prosopis stephaniana","S.inca"="Salsola incanescens",
                                 "S.rosm"="Seidlitzia rosmarinus","T.hisp"="Tamrix hispida")

levels(Article1Meta$SOIL)<-list("Arid"="Arid soil","Saline"="Saline Soil")
levels(Article1Meta$TISSUE)
levels(Article1Meta$MEDIA)<-list("PDA"= "PDA", "PDA+Plant"="PDA+Plant extract")
## Subset OTU frequency dataframe
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
colnames(Article1OTU)


Article1Order = subset (orderabund, MetaData$HOST %in% c("Alhagi persarum", "Artemisia sieberi", "Haloxylon ammodendron", 
                                                     "Launaea acunthodes",
                                                     "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                                   "Tamrix hispida"))
# check to see if it worked 
rownames(Article1Order)==rownames(Article1Meta)
class(Article1Order)
colnames(Article1Order)
View(Article1Order)
# creat a Pie chart 
order.colsum<-colSums(Article1Order)
order.slic<- c(265,170,589,301,1300,443,4944,88,1736,89,741)  #get the valus from order.colsum
order.lbls<- c("Eurotiales","Diaporthales","Ophiostomatales","Boletales","Xylariales","Hypocreales"," Pleosporales",
           "Saccharomycetales","Sordariales ", "Amphisphaeriales","Unknown ")

order.Percent<-round(order.slic/sum(order.slic)*100, digits=2)
order.lbls <- paste(order.lbls, order.Percent)
order.lbls<-paste(order.lbls,"%",sep="")
pie(order.slic,labels =order.lbls, col = c("firebrick","indianred1","skyblue1","magenta",
                                   "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                   "yellow", "springgreen2") , main = "Order", cex=1,border = NA,cex.main= 1.5, radius = 0.7)

# aggregate the Order abundance for all variables
aggregateHOST = aggregate (.~Article1Meta$HOST,Article1Order, sum)
aggregateSITE = aggregate (.~Article1Meta$SITE,Article1Order, sum)
aggregateSOIL = aggregate (.~Article1Meta$SOIL,Article1Order, sum)
aggregateTIME = aggregate (.~Article1Meta$TIME,Article1Order, sum)
aggregateTISSUE = aggregate (.~Article1Meta$TISSUE,Article1Order, sum)
View(aggregateSITE)

####### plot porpotion of each order  (HOST)
dat.Order.Host <- read.table (text = "    A.persarum A.sieberi  H.ammodendron L.acunthodes P.stephaniana S.incanescens S.rosmarinus T.hispida
                  Eurotiales  0 0  170 0 95 0 0 0
                              Diaporthales   0 0 170 0 0 0 0 0 
                              Ophiostomatales   0 0 589 0 0 0 0 0
                              Boletales   0 0 301 0 0 0 0 0
                              Amphisphaeriales   89 0 0 0 0 0 0 0
                              Saccharomycetales 0 0 0 0 0 0 0 0 
                              Hypocreales 82 0 0 0 0 0 195 166
                              Pleosporales 103 572 1513 564 886 264 1042 0 
                              Xylariales 72 0 580 0 0 201 447 0
                              Sordariales 0 0 418 272 387 0 546 113
                              Unknown 221 0 0 0 0 0 348 172",sep = "",header = TRUE)
View(dat.Order.Host)
library(reshape)
datH <- melt(cbind(dat.Order.Host, ind = rownames(dat.Order.Host)), id.vars = c('ind'))
View(datH)
library(scales)
ggplot(datH,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())+
  xlab("Host plant species")+ ylab("Proportional frequency")+
  labs(fill = "Order")

# plot porpotion of each order  (SITE)
dat.Order.Site <- read.table (text = " Garmsar  HajAli  HozeSoltan  Maranjab  RigBoland
                  Eurotiales  265 0 0 0 0
                  Diaporthales 170 0 0 0 0   
                  Ophiostomatales  0 0 589 0 0  
                  Boletales 0 0 301 0 0  
                  Amphisphaeriales  0 0 89 0 0 
                  Saccharomycetales 0 0 0 0 88
                  Hypocreales 0 196 165 0 82
                  Pleosporales 180 639 2959 1025 141
                  Xylariales 0 0 1031 269 0 
                  Sordariales 0 135 625 796 180 
                  Unknown 0 1313 98 0 515 ",sep = "",header = TRUE)

View(dat.Order.Site)
datS <- melt(cbind(dat.Order.Site, ind = rownames(dat.Order.Site)), id.vars = c('ind'))
View(datS)

ggplot(datS,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())


# plot porpotion of each order  (TIME)

dat.Order.Time <- read.table (text = " Summer2016 Summer2017 Winter2015 Winter2016
                              Eurotiales  0 265 0 0
                              Diaporthales 0 170 0 0  
                              Ophiostomatales  235 0 354 0
                              Boletales 186 0 115 0  
                              Amphisphaeriales  89 0 0 0 
                              Saccharomycetales 0 0 0 88
                              Hypocreales 0 349 0 94
                              Pleosporales 1511 1807 1171 455
                              Xylariales 436 541 323 0
                              Sordariales 186 1066 257 227
                              Unknown 0 475 58 208 ",sep = "",header = TRUE)


datT <- melt(cbind(dat.Order.Time, ind = rownames(dat.Order.Time)), id.vars = c('ind'))


ggplot(datT,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())


# plot porpotion of each order  (SOIL)

dat.Order.Soil <- read.table (text = " Arid  Saline
                              Eurotiales  265 0
                              Diaporthales 170 0  
                              Ophiostomatales  0 589 
                              Boletales 0 301 
                              Amphisphaeriales  0 89
                              Saccharomycetales 88 0
                              Hypocreales 112 331
                              Pleosporales 1385 3559 
                              Xylariales 269 1031
                              Sordariales 1056 680
                              Unknown 577 164",sep = "",header = TRUE)


datSO <- melt(cbind(dat.Order.Soil , ind = rownames(dat.Order.Soil )), id.vars = c('ind'))


ggplot(datSO,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())


# plot porpotion of each order  (TISSUE)

dat.Order.Tissue <- read.table (text = "  Branch  Leaf  Root
                              Eurotiales  120 145 0
                              Diaporthales 120 50 0  
                              Ophiostomatales  153 136 300 
                              Boletales 301 0 0 
                              Amphisphaeriales  1 88 0
                              Saccharomycetales 18 31 39
                              Hypocreales 38 24 381
                              Pleosporales 606 1867 2471 
                              Xylariales 110 323 867
                              Sordariales 802 671 263
                              Unknown 453 285 3",sep = "",header = TRUE)


datTI <- melt(cbind(dat.Order.Tissue , ind = rownames(dat.Order.Tissue)), id.vars = c('ind'))


ggplot(datTI,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())

#aggregate 3 Varibles 

aggregate(Article1Order $  Eurotiales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Diaporthales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Ophiostomatales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Boletales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Amphisphaeriales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Saccharomycetales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Hypocreales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Pleosporales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Xylariales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Sordariales  ~ HOST + SITE , data = Article1Meta, sum)
aggregate(Article1Order $  Unknown  ~ HOST + SITE , data = Article1Meta, sum)


myd1 <- data.frame(  host  = c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8),   
                     site = c("A","B","c","d","E","A","B","c","d","E","A","B","c","d","E","A","B","c","d","E","A","B","c","d","E","A","B","c","d","E","A","B","c","d","E","A","B","c","d","E"),   
                     Eurotiales = c(170,0,0,0,0,95,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                     Diaporthales = c(170,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     Ophiostomatales = c(0,0,589,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     Boletales = c(0,0,301,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     Amphisphaeriales = c(0,0,0,0,0,0,0,0,0,0, 0,0,89,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     Saccharomycetales = c(0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,88,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     Hypocreales = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,0,30,165,0,0,0,166,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     Pleosporales = c(180,0,914,419,0,0,0,771,115,0,0,0,103,0,0,0,202,699,0,141,0,0,0,0,0,0,437,0,135,0,0,0,331,233,0,0,0,141,123,0),
                     Xylariales = c(0,0,311,269,0,0,0,0,0,0,0,0,72,0,0,0,0,447,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,201,0,0),
                     Sordariales = c(0,0,95,323,0,0,0,0,387,0,0,0,0,0,0,0,135,344,0,67,0,0,0,0,113,0,0,0,0,0,0,0,186,86,0,0,0,0,0,0),
                     Unknown = c(0,0,0,0,0,0,0,0,0,0,0,0,0,221,0,0,133,98,0,117,0,0,0,0,172,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) 
View(myd1)
# rshaping data to long form for ggplot2 

meltd<- melt(myd1 , id.vars=1:2) 
View (meltd )
#plot 

ggplot(meltd, aes(x=host, y=value, fill=variable)) +
  geom_bar(stat="identity") + facet_wrap(~site) + theme_bw()


##### plot soil& host&orders

#Eurotiales data frame
host.s.e <-aggregate(Article1Order $  Eurotiales ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.e)
View(host.s.e)
host.s.e$value=host.s.e$`Article1Order$Eurotiales`
host.s.e$`Article1Order$Eurotiales`<-NULL
host.s.e$Order<-c("Eurotiales")
View(host.s.e)

#Diaporthales data frame
host.s.d <-aggregate(Article1Order $  Diaporthales  ~ HOST + SOIL , data = Article1Meta, sum)
View(host.s.d)
host.s.d$value=host.s.d$`Article1Order$Diaporthales`
host.s.d$`Article1Order$Diaporthales`<-NULL
host.s.d$Order<-c("Diaporthales")
View(host.s.d)


#Ophiostomatales data frame
host.s.o <-aggregate(Article1Order $  Ophiostomatales  ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.o)
View (host.s.o)
host.s.o$value=host.s.o$`Article1Order$Ophiostomatales`
host.s.o$`Article1Order$Ophiostomatales`<-NULL
host.s.o$Order<-c("Ophiostomatales")
View (host.s.o)


#Boletales data frame
host.s.b <-aggregate(Article1Order $  Boletales  ~ HOST + SOIL , data = Article1Meta, sum)
View(host.s.b)
host.s.b$value=host.s.b$`Article1Order$Boletales`
host.s.b$`Article1Order$Boletales`<-NULL
host.s.b$Order<-c("Boletales")
View(host.s.b)


#Amphisphaeriales data frame
host.s.a<-aggregate(Article1Order $  Amphisphaeriales  ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.a)
View(host.s.a)
host.s.a$value=host.s.a$`Article1Order$Amphisphaeriales`
host.s.a$`Article1Order$Amphisphaeriales`<-NULL
host.s.a$Order<-c("Amphisphaeriales")
View(host.s.a)


#Saccharomycetales data frame
host.s.s<-aggregate(Article1Order $  Saccharomycetales  ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.s)
View(host.s.s)
host.s.s$value=host.s.s$`Article1Order$Saccharomycetales`
host.s.s$`Article1Order$Saccharomycetales`<-NULL
host.s.s$Order<-c("Saccharomycetales")
View(host.s.s)

#Hypocreales data frame
host.s.h<-aggregate(Article1Order $  Hypocreales  ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.h)
View(host.s.h)
host.s.h$value=host.s.h$`Article1Order$Hypocreales`
host.s.h$`Article1Order$Hypocreales`<-NULL
host.s.h$Order<-c("Hypocreales")
View(host.s.h)

#Pleosporales data frame
host.s.p<-aggregate(Article1Order $  Pleosporales  ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.p)
View(host.s.p)
host.s.p$value=host.s.p$`Article1Order$Pleosporales`
host.s.p$`Article1Order$Pleosporales`<-NULL
host.s.p$Order<-c("Pleosporales")
View(host.s.p)

#Xylariales data frame
host.s.x<-aggregate(Article1Order $  Xylariales  ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.x)
View(host.s.x)
host.s.x$value=host.s.x$`Article1Order$Xylariales`
host.s.x$`Article1Order$Xylariales`<-NULL
host.s.x$Order<-c("Xylariales")
View(host.s.x)


# Sordariales data frame
host.s.ss<-aggregate(Article1Order $  Sordariales  ~ HOST + SOIL , data = Article1Meta, sum)
class(host.s.ss)
View(host.s.ss)
host.s.ss$value=host.s.ss$`Article1Order$Sordariales`
host.s.ss$`Article1Order$Sordariales`<-NULL
host.s.ss$Order<-c("Sordariales")
View(host.s.ss)


# Unknown data frame
host.s.u <-aggregate(Article1Order $  Unknown  ~ HOST + SOIL , data = Article1Meta, sum)
View(host.s.u)
host.s.u$value=host.s.u$`Article1Order$Unknown`
host.s.u$`Article1Order$Unknown`<-NULL
host.s.u$Order<-c("Unknown")
View(host.s.u)

# a new dataframe for ploting
eample1<-rbind (host.s.d, host.s.e, host.s.o,host.s.b,host.s.a, host.s.s, host.s.h, host.s.p, host.s.x, host.s.ss, host.s.u)
View(eample1 )

ggplot(eample1, aes(x=HOST, y=value, fill=Order)) +
  geom_bar(stat="identity") + facet_wrap(~SOIL) + theme_bw()

#### PLOT!
library(reshape)
eample1<-rbind (host.s.d, host.s.e, host.s.o,host.s.b,host.s.a, host.s.s, host.s.h, host.s.p, host.s.x, host.s.ss, host.s.u)
View(eample1 )
library(scales)
ggplot(eample1,aes(x = HOST, y = value,fill = Order)) + 
  geom_bar(position = "fill",stat = "identity") +facet_wrap(~SOIL) + theme_bw()+
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())+
  xlab("Host plant species")+ ylab("Proportional frequency")+
  labs(fill = "Order")





















