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
Article1Meta = subset (MetaData, HOST%in%c("Alhagi persarum","Artemisia sieberi", "Haloxylon ammodendron", 
                                            "Launaea acunthodes",
                                            "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                            "Tamrix hispida"))
# some how this still shows the 40 hists!! I am trying another way to fix it!
write.csv(Article1Meta, file = "testmetadata.csv")

Article1Meta<-read.csv("testmetadata.csv",header = T, row.names = 1)
levels(Article1Meta$HOST)

# remane the long variable levels
levels(Article1Meta$SITE)
levels(Article1Meta$SITE)<- list("Garmsar"="Garmsar","HajAli"="Haj Ali Gholi Lake",
                             "Hoze"="Hoze Soltan Lake","Rig"="Rig-Boland Desert", "Sorkhe"="Sorkhe")

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

aggregate(LREwh64..Neocamarosporium.chichastianum ~ Article1Meta$HOST, Article1OTU, sum)
aggregate(SREwh19.Neocamarosporium.goegapense ~ Article1Meta$HOST, Article1OTU, sum)

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
# R.m2<-glm(formula =Richness.art1~SOIL*TISSUE+HOST+TIME+SITE,data = MetaRich.ART1,
#           family=poisson(link = "log"))
# summary(R.m2)
# anova(R.m2, test = "Chisq")
# AIC(R.m2)
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
#MODEL SUMMARY FOR SHANNON
shannon.sum<-summary(sh.m)
#ANOVA RESULS FOR SHANON
shannon.anov<-anova(sh.m, test="F")
AIC(sh.m)

###############
#SIMPSON MODEL
###############
simp.m<-lm(formula =log(simpson.art1)~SOIL*TISSUE+HOST+TIME+SITE,data = MetaNotOne.art1)
#MODEL SUMMARY FOR SIMPSON
simp.sum<-summary(simp.m)
#ANOVA RESULTS FOR SIMPSON
simp.anov<-anova(simp.m, test="F")
AIC(simp.m)

######## Step 4: Choose the best way to visualize the results








####################################################
####################################################
# COMMUNITY COMPOSITION ANALYSIS
####################################################
####################################################

# for this we are using Permutational multivariate analysis of variance (PERMANOVA)
#this function from vegan does it

comm.anova<-adonis(formula=AbundNotZero.art1~SOIL*TISSUE+HOST+TIME+SITE, data= MetaRich.ART1,
       permutations = 999, method = "bray",by=NULL)
# if you get an error about the memory allocation run the following lines:
# memory.limit()
# memory.limit(size=30000)

## Coefficients
comm.anova.coef = as.data.frame(comm.anova$coefficients)

## which OTUs significantly affected by any variables:
#####################
####### NMDS PLOTS 
#####################
# first I am computing an NMDS matrix for all of the OTUs
nmds.art1<-metaMDS(AbundNotZero.art1, distance = "bray", k= 2, trymax = 20)

plot(nmds.art1)
# red + are species and black dots are sites (samples)
dev.off()
###########################################
#this is the base of our plots: 
NM.pl<-ordiplot(nmds.art1,type = "none")
# you can show species ot samples or both
# points(NM.pl,"sites",pch = 19, cex= 0.5, col="grey30")
points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.6)
###########################################
# show soil communities:
ordihull(nmds.art1, MetaRich.ART1$SOIL,cex=1.5,
         draw="line", col= "grey20",
         lwd = 2, lty = 2,
         show.groups=(c("Arid")))#Arid soil endophyte community
ordihull(nmds.art1, MetaRich.ART1$SOIL,cex=1.5,
         draw="line", col= "grey20",
         lwd = 2,lty = 9,
         show.groups=(c("Saline")))# Saline soil endophyte community
# show different sampling sites:
ordibar(nmds.art1, MetaRich.ART1$SITE,
            col= "green",lwd = 2,
            kind="se", conf=0.99,
            show.groups=(c("Garmsar")))
ordibar(nmds.art1, MetaRich.ART1$SITE,
             col= "blue",lwd = 2,
            kind="se", conf=0.99,
            show.groups=(c("Haj Ali Gholi Lake")))
ordibar(nmds.art1, MetaRich.ART1$SITE,
            col= "red",lwd = 2,
           kind="se", conf=0.99,
            show.groups=(c("Hoze Soltan Lake")))
ordibar(nmds.art1, MetaRich.ART1$SITE,
             col= "yellow",lwd = 2,
             kind="se", conf=0.99,
            show.groups=(c("Maranjab Desert")))
ordibar(nmds.art1, MetaRich.ART1$SITE,
        col= "violetred4",lwd = 2,
         kind="se", conf=0.99, 
        show.groups=(c("Rig-Boland Desert")))
pl.legend = legend("topleft", c("Arid soil","Saline soil","Garmsar", "Haj Ali Gholi Lake",
                                "Hoze Soltan Lake","Maranjab Desert","Rig-Boland Desert"), 
                  col=c("grey20","grey20","green","blue","red","yellow","violetred4"),
                  lty = c(2,9,1,1,1,1,1), border="white", bty="n")
#or:
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1,
            draw="polygon", col= "green",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Garmsar")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1.5,
            draw="polygon", col= "blue",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Haj Ali Gholi Lake")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1,
            draw="polygon", col= "red",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Hoze Soltan Lake")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1.5,
            draw="polygon", col= "yellow",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Maranjab Desert")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1.5,
            draw="polygon", col= "violetred4",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Rig-Boland Desert")))

?legend
?points
?ordiellipse
?ordiplot
?ordibar
# in case we  decided to show the data in other shapes:
# try a new plot for sites with different x & y limites
dev.off()
site.pl<-ordiplot(nmds.art1,type = "none", xlim = c(-8,8), ylim = c(-4,4))
points(site.pl,"species",pch = 2, col= "grey20", cex= 0.6)
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1,
            draw="polygon", col= "green",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Garmsar")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1.5,
            draw="polygon", col= "blue",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Haj Ali Gholi Lake")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1,
            draw="polygon", col= "red",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Hoze Soltan Lake")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1.5,
            draw="polygon", col= "yellow",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Maranjab Desert")))
ordiellipse(nmds.art1, MetaRich.ART1$SITE,cex=1.5,
            draw="polygon", col= "violetred4",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Rig-Boland Desert")))

dev.off()
?points
levels(MetaRich.ART1$SITE)
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
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Leaf")))
ordiellipse(nmds.art1, MetaRich.ART1$TISSUE,cex=1.5,
            draw="polygon", col= "blue",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Branch")))
ordiellipse(nmds.art1, MetaRich.ART1$TISSUE,cex=1,
            draw="polygon", col= "red",
            alpha=200, kind="se", conf=0.99,
            show.groups=(c("Root")))
legend("topleft", c("Arid soil","Saline soil","Leaf", "Branch",
                    "Root"), 
       col=c("grey20","grey20","green","blue","red"),
       lty = c(2,9,1,1,1), border="white", bty="n")


levels(MetaRich.ART1$TISSUE)



