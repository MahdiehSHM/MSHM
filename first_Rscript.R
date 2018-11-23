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
memory.limit()
memory.limit(size=30000)
comm.anova<-adonis(formula=AbundNotZero.art1~SOIL*TISSUE+HOST+TIME+SITE, data= MetaRich.ART1,
       permutations = 999, method = "bray",by=NULL)
# if you get an error about the memory allocation run the following lines:


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

####use glm.nb for each OTU to figur out each of their frequency recation to every variables
# example: first OTU

OTU1.model<-glm.nb(formula =AbundNotZero.art1$APE.se5.Staphylotrichum.coccosporum~SOIL*TISSUE+ HOST+ TIME+
                           SITE,data = MetaRich.ART1,link = "log")
OTU1.anov<-anova(OTU1.model,test = "Chisq")

OTU2.model<-glm.nb(formula =Article1OTU$TPEsh28.Humicola.fuscoatra~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU2.anov<-anova(OTU2.model,test = "Chisq")

OTU3.model<-glm.nb(formula =Article1OTU$PFE.sh7..Rosellinia.limonispora~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU3.anov<-anova(OTU3.model,test = "Chisq")

OTU4.model<-glm.nb(formula =Article1OTU$LREwh64..Neocamarosporium.chichastianum~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU4.anov<-anova(OTU4.model,test = "Chisq")

OTU5.model<-glm.nb(formula =Article1OTU$SREwh22.Preussia.minimoides~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU5.anov<-anova(OTU5.model,test = "Chisq")

OTU6.model<-glm.nb(formula =Article1OTU$THE.we10..Aporospora.terricola~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU6.anov<-anova(OTU6.model,test = "Chisq")

OTU7.model<-glm.nb(formula =Article1OTU$THE.ss3.Fusarium.sp.~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU7.anov<-anova(OTU7.model,test = "Chisq")

OTU8.model<-glm.nb(formula =Article1OTU$SIE.sh1.Briansuttonomyces.eucalypti~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU8.anov<-anova(OTU8.model,test = "Chisq")

OTU9.model<-glm.nb(formula =Article1OTU$RAE.sh12.Acrocalymma.vagum~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU9.anov<-anova(OTU9.model,test = "Chisq")

OTU10.model<-glm.nb(formula =Article1OTU$PSE.wh14.Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
                     SITE,data = Article1Meta,link = "log")
OTU10.anov<-anova(OTU10.model,test = "Chisq")

OTU11.model<-glm.nb(formula =Article1OTU$PSE.wh40.Dimorphosporicola.tragani~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU11.anov<-anova(OTU11.model,test = "Chisq")

OTU12.model<-glm.nb(formula =Article1OTU$PSE.wh66.Comoclathris.italica~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU12.anov<-anova(OTU12.model,test = "Chisq")

OTU13.model<-glm.nb(formula =Article1OTU$PSE.ss7.Penicillium.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU13.anov<-anova(OTU13.model,test = "Chisq")

OTU14.model<-glm.nb(formula =AbundNotZero.art1$PSE.we4..Chaetomium.globosum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = MetaRich.ART1)
OTU14.anov<-anova(OTU14.model,test = "Chisq")

OTU15.model<-glm(formula =AbundNotZero.art1$PSE.se8.Podospora.minicauda~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU15.anov<-anova(OTU15.model,test = "Chisq")

OTU16.model<-glm(formula =AbundNotZero.art1$PSE.we8.Alternaria.chlamydospora~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU16.anov<-anova(OTU16.model,test = "Chisq")

OTU17.model<-glm(formula =AbundNotZero.art1$LDE.se7.Coniothyrium.aleuritis~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU17.anov<-anova(OTU17.model,test = "Chisq")

OTU18.model<-glm.nb(formula =Article1OTU$LAE.se5.Sordaria.humana~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU18.anov<-anova(OTU18.model,test = "Chisq")

OTU19.model<-glm.nb(formula =Article1OTU$HAE.se5.Camarosporomyces.flavigenus~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU19.anov<-anova(OTU19.model,test = "Chisq")

OTU20.model<-glm.nb(formula =Article1OTU$HAE.ss4.Cytospora.chrysosperma~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU20.anov<-anova(OTU20.model,test = "Chisq")

OTU21.model<-glm(formula =AbundNotZero.art1$HAE.we5.Coniolariella.sp.~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU21.anov<-anova(OTU21.model,test = "Chisq")

OTU22.model<-glm.nb(formula =Article1OTU$HAE.wh26.Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU22.anov<-anova(OTU22.model,test = "Chisq")

OTU23.model<-glm.nb(formula =Article1OTU$HAE.wh10.Raffaelea.montetyi~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU23.anov<-anova(OTU23.model,test = "Chisq")

OTU24.model<-glm.nb(formula =Article1OTU$HAE.wh65.Coniophora.marmorata~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU24.anov<-anova(OTU24.model,test = "Chisq")

OTU25.model<-glm(formula =AbundNotZero.art1$HAE.se9.Chaetomium.nigricolor~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU25.anov<-anova(OTU25.model,test = "Chisq")

OTU26.model<-glm(formula =AbundNotZero.art1$HAE.se1.Acrocalymma.sp.~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU26.anov<-anova(OTU26.model,test = "Chisq")

OTU27.model<-glm.nb(formula =Article1OTU$SREwh19.Neocamarosporium.goegapense~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU27.anov<-anova(OTU27.model,test = "Chisq")

OTU28.model<-glm.nb(formula =Article1OTU$APEsh6.Dictyosporium.digitatum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU28.anov<-anova(OTU28.model,test = "Chisq")

OTU29.model<-glm.nb(formula =Article1OTU$APE.sh8.Pestalotiopsis.vismiae~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU29.anov<-anova(OTU29.model,test = "Chisq")

OTU30.model<-glm.nb(formula =Article1OTU$APE.se3.Dactylonectria.macrodidyma~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU30.anov<-anova(OTU30.model,test = "Chisq")

OTU31.model<-glm.nb(formula =Article1OTU$APE.sh5.Nigrospora.sphaerica~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU31.anov<-anova(OTU31.model,test = "Chisq")

OTU32.model<-glm.nb(formula =Article1OTU$SREwh18...Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU32.anov<-anova(OTU32.model,test = "Chisq")

OTU33.model<-glm.nb(formula =Article1OTU$LAEsh5.Coniolariella.ershadii~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU33.anov<-anova(OTU33.model,test = "Chisq")

OTU34.model<-glm(formula =AbundNotZero.art1$LAE.se3.Neosetophoma.lunariae~SOIL*TISSUE+ HOST+ TIME+
                   SITE,data = MetaRich.ART1, family=poisson(link = "log"))
OTU34.anov<-anova(OTU34.model,test = "Chisq")

OTU35.model<-glm.nb(formula =Article1OTU$LAE.SH.7.Muriphaeosphaeria.viburni~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU35.anov<-anova(OTU35.model,test = "Chisq")

OTU36.model<-glm.nb(formula =Article1OTU$LAE.sh1.Acrocalymma.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU36.anov<-anova(OTU36.model,test = "Chisq")

OTU37.model<-glm.nb(formula =Article1OTU$SRE.sh30..Preussia.grandispora~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU37.anov<-anova(OTU37.model,test = "Chisq")

OTU38.model<-glm.nb(formula =Article1OTU$SRE.ss.4.Neocamarosporium.sp.~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU38.anov<-anova(OTU38.model,test = "Chisq")

OTU39.model<-glm.nb(formula =Article1OTU$SRE.ws8.Botryotrichum.murorum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU39.anov<-anova(OTU39.model,test = "Chisq")

OTU40.model<-glm.nb(formula =Article1OTU$SRE.ws10.Sarocladium.kiliense~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU40.anov<-anova(OTU40.model,test = "Chisq")

OTU41.model<-glm.nb(formula =Article1OTU$SRE.wh16.Paracamarosporium.hawaiiense~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU41.anov<-anova(OTU41.model,test = "Chisq")

OTU42.model<-glm.nb(formula =Article1OTU$SRE.wh13.Ovatospora.senegalensis~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU42.anov<-anova(OTU42.model,test = "Chisq")

OTU43.model<-glm.nb(formula =Article1OTU$SRE.we6.Fusariella.sinensis~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU43.anov<-anova(OTU43.model,test = "Chisq")

OTU44.model<-glm.nb(formula =Article1OTU$SRE.we.10.pichia.kudriavzevii~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU44.anov<-anova(OTU44.model,test = "Chisq")

OTU45.model<-glm.nb(formula =Article1OTU$SRE.sh7.Chaetomium.cucumericola~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU45.anov<-anova(OTU45.model,test = "Chisq")

OTU46.model<-glm.nb(formula =Article1OTU$SRE.sh5.Fusarium.redolens~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU46.anov<-anova(OTU46.model,test = "Chisq")

OTU47.model<-glm.nb(formula =Article1OTU$SRE.sh9.Preussia.intermedia~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU47.anov<-anova(OTU47.model,test = "Chisq")

OTU48.model<-glm.nb(formula =Article1OTU$SRE.sh4.Penicillium.vinaceum~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU48.anov<-anova(OTU48.model,test = "Chisq")

OTU49.model<-glm.nb(formula =Article1OTU$SRE.sh3.Trichoderma.rifaii~SOIL*TISSUE+ HOST+ TIME+
                      SITE,data = Article1Meta,link = "log")
OTU49.anov<-anova(OTU49.model,test = "Chisq")


########## SOIL 
# OTU1-OTU2-OTU3-OTU4-OTU5-OTU6-OTU7-OTU8-OTU9-OTU10-OTU11-OTU12-OTU13-OTU18-OTU19-OTU20-OTU22-OTU23-OTU24-OTU27-
# OTU28-OTU29-OTU30-OTU31-OTU32-OTU33-OTU35-OTU36-OTU37-OTU38-OTU39-OTU40-OTU41-OTU42-OTU43-OTU44-OTU45-OTU46-
# OTU47-OTU48-OTU49

########## Tissue
# OTU1-OTU2-OTU3-OTU4-OTU5-OTU6-OTU7-OTU9-OTU10-OTU11-OTU12-OTU13-OTU18-OTU19-OTU20-OTU22-OTU24-OTU27-
# OTU28-OTU29-OTU30-OTU31-OTU32-OTU33-OTU35-OTU36-OTU37-OTU39-OTU40-OTU41-OTU42-OTU43-OTU45-OTU46-OTU47-OTU48-OTU49

########## HOST
# OTU1-OTU2-OTU3-OTU4-OTU5-OTU6-OTU7-OTU8-OTU9-OTU10-OTU11-OTU12-OTU13-OTU18-OTU19-OTU20-OTU22-OTU23-OTU24-OTU27-OTU28-
# OTU29-OTU30-OTU31-OTU32-OTU33-OTU35-OTU36-OTU37-OTU38-OTU39-OTU40-OTU41-OTU43-OTU44-OTU45-OTU46-OTU47-OTU48-OTU49

########## TIME
# OTU2-OTU3-OTU4-OTU5-OTU7-OTU8-OTU9-OTU11-OTU12-OTU13-OTU20-OTU22-OTU23-OTU24-OTU27-OTU32-OTU33-OTU37-
# OTU38-OTU39-OTU40-OTU44-OTU45-OTU46-OTU47-OTU48-OTU49

########## SITE
#OTU2-OTU3-OTU4-OTU33-OTU38-OTU39-OTU40-OTU44-OTU45-OTU46-OTU47-OTU48-OTU49

########## SITE:Tissue
#OTU2

##############################################
##############################################
# Variation partitioning
##############################################
##############################################

?varpart
var.art1<-varpart(Article1OTU,Article1Meta$SOIL,
                 Article1Meta$HOST, Article1Meta$SITE,Article1Meta$TISSUE)
plot(var.art1)

row.names(AbundNotZero.art1)==row.names(MetaRich.ART1)

View(AbundNotZero.art1)
View(MetaRich.ART1)
# >0 samles only
var.art2<-varpart(AbundNotZero.art1,MetaRich.ART1$SOIL,
                  MetaRich.ART1$HOST, MetaRich.ART1$SITE,MetaRich.ART1$TISSUE, 
                  transfo="hellinger")
plot(var.art2)
#without soil
var.art3<-varpart(AbundNotZero.art1,MetaRich.ART1$TIME,
                  MetaRich.ART1$HOST, MetaRich.ART1$SITE,MetaRich.ART1$SOIL, 
                  transfo="hellinger")
plot(var.art3)


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

#coloring:
soilfactor<-factor(MetaRich.ART1$SOIL)
colvec <-  c("red","blue")
cols<-colvec[soilfactor]
#scors:
rda.scores <- scores(rda.art1, display = 'bp')
mul <- ordiArrowMul(rda.scores, fill = 0.75)
 
?scores()
dev.off()
plot(rda.art1, type = "n")
points(rda.art1, display = "sites", col = colvec[soilfactor], 
       pch = (16:17)[soilfactor],cex=0.85)#salin=blue
#points(rda.art1, display = "species", pch = "+")
arrows(0, 0, mul * rda.scores[,1], mul * rda.scores[,2],
       length = 0.05, col ="black", lwd=2)
labs <- rownames(rda.scores)
# labs<-c("Ece","EC","Cle","pHe")
text(ordiArrowTextXY(mul * rda.scores, labs), labs)
legend("topleft", c("Arid soil","Saline soil"), 
             col=c("red","blue"),
       pch = c(16,17), border="white", bty="n")
?arrows
?legend
?points()


##############################################
##############################################
# Greenhouse data analysis
##############################################
##############################################
# data input
GH.data<-read.csv("Green house data.csv", header = T, row.names = 1)
# factors: Fungi, Drought, Salinity

?aov
GH.1<-aov(Lshoot~Fungi*Drought*Salinity,data = GH.data)
summary(GH.1)
par(mfrow = c(2, 3))
plot(GH.1)
hist(GH.data$Wshoot)
hist(GH.data$Lshoot)
hist(GH.data$DWshoot)
hist(GH.data$Lroot)
hist(GH.data$Wroot)
hist(GH.data$DWroot)
boxplot(GH.data$DWshoot)
boxplot(GH.data$Lshoot)

#outlier fixed
dev.off()
# try GLM
glm.Lshoot<-glm(Lshoot~Fungi*Drought*Salinity,data = GH.data, family =poisson(link = "log"))
summary(glm.Lshoot)
anova(glm.Lshoot, test = "Chisq")
par(mfrow = c(2, 2))
plot(glm.Lshoot)

# how to show this?

ggplot(GH.data, aes(x=Drought,y=Lshoot, fill=Salinity)) + 
facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()
# please try to make the lables presentable for the paper. you can find the right codes in workshop ggplot2 script
# do the same for the rest of variables 





