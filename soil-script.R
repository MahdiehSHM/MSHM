
############################################
# desert endophytes
############################################
library(vegan)
library(rjags)
library(MASS)
library(ggplot2)
library(ggtree)
library(tidyverse)
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


####### new subset : remove GARMSAR data from Article1Meta&Article1OTU


# subset OTU abundance data FIRST
levels(Article1Meta$SITE)
Article1OTU.n = subset (Article1OTU, Article1Meta$SITE!="Garmsar")
Article1OTU<-Article1OTU.n[, colSums(Article1OTU.n != 0) > 0]
#export the abundance data of article 1:
write.csv(Article1OTU,file = "art1.abundance.cav")
#this is the abundance data for article 1 :art1.abundance.cav
#SUBSET METADATA
Article1data <- subset(Article1Meta, Article1Meta$SITE!="Garmsar")
write.csv(Article1data, file = "art1data.csv")
#this is the metadata for article 1 :art1data.csv

# I am using the same name for metada object this way I only need to run everything again 
#to have the new results. this overwrites all the objects we had in the old script
# if you need the old results you should run the codes in first-Rscript again

Article1Meta<-read.csv(file = "art1data.csv",header = TRUE, row.names = 1)
#now we have 4 sites: 2 have salin soil and the other 2 have dry soil
levels(Article1Meta$SOIL)
levels(Article1Meta$SOIL)<-list("Dry soil" = "Arid", "Saline soil"= "Saline")
#############################################################################
######## Step 3: Find the right kind of analysis for each research questions:


#######################################
#IR model
#######################################
hist(Article1Meta$IR)
IR.m<-glm(formula =IR~SOIL*TISSUE+HOST+season+SITE,data = Article1Meta,
                family=poisson(link = "log"))
#MODEL SUMMARY FOR REACHNESS
IR.summ<-summary(IR.m)
#ANOVA RESULTS FOR RICHNESS
IR.anova<-anova(IR.m, test = "Chisq")
boxplot(IR~SOIL,data = Article1Meta)
boxplot(IR~HOST,data = Article1Meta)
aggregate(IR~HOST,data= Article1Meta,mean )
aggregate(IR~TISSUE,data= Article1Meta,mean )


#plot the interaction between soil and organ
ggplot (Article1Meta, aes(x=TISSUE,y=IR))+
  facet_wrap(~SOIL)+ geom_boxplot(position = "dodge", width=0.5) + theme_bw()+
xlab("Organ type")+ylab("Isolation rate")+ scale_x_discrete(labels = c("Twig", "Leaf","Root"))+ theme(legend.position = "non")



#########################################################
#########################################################
##### Diversity indices
#########################################################
#########################################################
#For diversity we are using model based aproaches:
## I am using the codes that we wrote earlier:
## first delete the IR in metadata.not sure if it is ok
Article1Meta$IR<-NULL

IR.art1= apply(Article1OTU,1, sum)
Article1Meta = cbind(Article1Meta, IR = IR.art1)
### add season as varible

Article1Meta$season<-Article1Meta$TIME
levels(Article1Meta$season) <- list(summer=c( "S2016","S2017"), winter=c("W2015","W2016"))
levels(Article1Meta$season) 
       
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

################################
# SELECTED DIVERSITY MODELs 
################################

################
# RICHNESS MODEL
################
Richness.m<-glm(formula =Richness.art1~SOIL*TISSUE+HOST+season+SITE,data = MetaRich.ART1,
                family=poisson(link = "log"))
#MODEL SUMMARY FOR REACHNESS
Rich.summ<-summary(Richness.m)
#ANOVA RESULTS FOR RICHNESS
Rich.anova<-anova(Richness.m, test = "Chisq")
AIC(Richness.m)
par(mfrow=(c(2,2)))
plot(Richness.m)
dev.off()
boxplot(Richness.art1~MetaRich.ART1$HOST)
MetaRich.ART1$richness.art1<-Richness.art1
aggregate(richness.art1~SITE, data=MetaRich.ART1,mean)



################
# SHANNON MODEL
################
sh.m<-lm(formula =log(shannon.art1)~SOIL*TISSUE+HOST+season+SITE,data = MetaNotOne.art1)
#MODEL SUMMARY FOR SHANNON
shannon.sum<-summary(sh.m)
#ANOVA RESULS FOR SHANON
shannon.anov<-anova(sh.m, test="F")
###############
#SIMPSON MODEL
###############
simp.m<-lm(formula =log(simpson.art1)~SOIL*TISSUE+HOST+season+SITE,data = MetaNotOne.art1)
#MODEL SUMMARY FOR SIMPSON
simp.sum<-summary(simp.m)
#ANOVA RESULTS FOR SIMPSON
simp.anov<-anova(simp.m, test="F")
AIC(simp.m)


#plot the interaction
dev.off()
?geom_step()
MetaNotOne.art1$shannon<-shannon.art1
MetaNotOne.art1$simpson.art1<-simpson.art1

ggplot (MetaNotOne.art1)+
  geom_boxplot(aes(x=TISSUE,y=simpson.art1),  width=0.5, position = "dodge") +
  geom_boxplot(aes(x=TISSUE,y=shannon), width=0.5, position = "dodge")+
  facet_wrap(~SOIL)+ 
  xlab("Organ type")+ylab("Diversity")+ scale_x_discrete(labels = c("Twig", "Leaf","Root"))


# فعلا جدا جدا میکشیم تا وقتیکه مشکل حل بشه
#simpson
ggplot (MetaNotOne.art1)+
  geom_boxplot(aes(x=TISSUE,y=simpson.art1),  width=0.5, position = "dodge") +
  facet_wrap(~SOIL)+ 
  xlab("Organ type")+ylab("Simpson Diversity")+ scale_x_discrete(labels = c("Twig", "Leaf","Root"))
aggregate(simpson.art1~SOIL+TISSUE, data = MetaNotOne.art1, mean)

#shanon
ggplot (MetaNotOne.art1)+
  geom_boxplot(aes(x=TISSUE,y=shannon), width=0.5, position = "dodge")+
  facet_wrap(~SOIL)+ 
  xlab("Organ type")+ylab("Diversity")+ scale_x_discrete(labels = c("Twig", "Leaf","Root"))
aggregate(shannon~SOIL+TISSUE, data = MetaNotOne.art1, mean)

?position_dodge()
?geom_boxplot

##############################
######## visualize the results
##############################


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
comm.anova<-adonis(formula=tran.abund.notzero~SOIL*TISSUE+HOST+season+SITE, data= MetaRich.ART1,
                   permutations = 999, method = "bray",by=NULL)
# if you get an error about the memory allocation run the above lines:memory.limit

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
levels(MetaRich.ART1$SOIL)
###########################################
# NMDS soil:
###########################################
dev.off()
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5))
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$SOIL,cex=1,alpha = 200, 
            draw="polygon", col= "red",lwd = 3,border="red",
            kind="se", conf=0.95,show.groups=(c("Arid soil")))
ordiellipse(nmds.art2, MetaRich.ART1$SOIL,cex=1,lwd = 3,alpha = 200, 
            draw="polygon", col= "Blue",border="Blue",
            kind="se", conf=0.95,show.groups=(c("Saline Soil")))
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
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$SITE,cex=1,alpha = 200, 
            draw="polygon", col= 1:5,border= 1:5,lwd=3, kind="se", conf=0.95)

legend("bottomright", c("Haj Ali Gholi Lake","Hoze Soltan Lake","Maranjab Desert","Rig-Boland Desert"), 
       fill= 1:5, border="white", bty="n")
###########################################
#NMDS ORGAN:
###########################################
levels(MetaRich.ART1$TISSUE)
dev.off()
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5))
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$TISSUE,cex=1,
            draw="polygon", col= c("blue","green","red"),border=c("blue","green","red"),lwd=3,alpha = 200, kind="se", conf=0.95)

legend("bottomright", c("Branch","Leaf","Root" ), 
       fill= c("blue","green","red"),border="white", bty="n")
###########################################
#NMDS time:
###########################################
levels(MetaRich.ART1$season)
dev.off()
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5))
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$season,cex=1,
            draw="polygon", col= c("red","blue"),border=c("red","blue"),lwd=3,alpha = 200, kind="se", conf=0.95)

legend("bottomright", c("Summer","Winter"), 
       fill= c("red","blue"),border="white", bty="n")
#################

########################################################
#### Individul reactions of OTUs to our variables
########################################################

#useing manyglm function from mvabund package
library(mvabund)
fun.Mvabund = mvabund(Article1OTU)

fun.Mvabund.m = manyglm (fun.Mvabund ~ SOIL*TISSUE+HOST+season+SITE, data= Article1Meta,
                         family="negative.binomial", show.residuals=T)

fun.Mvabund.m.sum = summary.manyglm (fun.Mvabund.m, nBoot=100, test="LR",p.uni="adjusted", 
                                     resamp="montecarlo")

## Analysis of variance explained by the predictors
# fun.Mvabund.m.anova = anova.manyglm (fun.Mvabund.m, nBoot=300, test="LR", p.uni="adjusted", 
#                                 resamp="montecarlo")
#options(max.print=25000)
fun.Mvabund.m.anova100 = anova.manyglm (fun.Mvabund.m, nBoot=100, test="LR", p.uni="adjusted", 
                                        resamp="montecarlo")

## Coefficients
fun.Mvabund.m.coef = as.data.frame(fun.Mvabund.m$coefficients)
## mean-centering the contrasts
fun.Mvabund.m.coef.contrast = fun.Mvabund.m.coef - apply (fun.Mvabund.m.coef,2,mean)

# ## OTUs significantly affected by the source??
mvabund.m.anova <- as.data.frame(fun.Mvabund.m.anova100$uni.p)
OTU.INTER<-colnames(mvabund.m.anova)[mvabund.m.anova["SOIL:TISSUE",]<= 0.05]#2otus affected
OTU.HOST<-colnames(mvabund.m.anova)[mvabund.m.anova["HOST",]<= 0.05]#47otus affected
OTU.SITE<-colnames(mvabund.m.anova)[mvabund.m.anova["SITE",]<= 0.05]#14otus affected
OTU.SOIL<-colnames(mvabund.m.anova)[mvabund.m.anova["SOIL",]<= 0.05]#45otus affected
OTU.SEASON<-colnames(mvabund.m.anova)[mvabund.m.anova["season",]<= 0.05]#26otus affected

# ####use glm.nb for each OTU to figur out each of their frequency recation to every variables
# 
# OTU1.model<-glm.nb(formula =AbundNotZero.art1$APE.se5.Staphylotrichum.coccosporum ~ SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU1.anov<-anova(OTU1.model,test = "Chisq")
# 
# OTU2.model<-glm.nb(formula =AbundNotZero.art1$TPEsh28.Humicola.fuscoatra~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU2.anov<-anova(OTU2.model,test = "Chisq")
# 
# OTU3.model<-glm.nb(formula =AbundNotZero.art1$PFE.sh7..Rosellinia.limonispora~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU3.anov<-anova(OTU3.model,test = "Chisq")
# 
# OTU4.model<-glm.nb(formula =AbundNotZero.art1$LREwh64..Neocamarosporium.chichastianum~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU4.anov<-anova(OTU4.model,test = "Chisq")
# 
# OTU5.model<-glm.nb(formula =AbundNotZero.art1$SREwh22.Preussia.minimoides~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU5.anov<-anova(OTU5.model,test = "Chisq")
# 
# OTU6.model<-glm.nb(formula =AbundNotZero.art1$THE.we10..Aporospora.terricola~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU6.anov<-anova(OTU6.model,test = "Chisq")
# 
# OTU7.model<-glm.nb(formula =AbundNotZero.art1$THE.ss3.Fusarium.sp.~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU7.anov<-anova(OTU7.model,test = "Chisq")
# 
# OTU8.model<-glm.nb(formula =AbundNotZero.art1$SIE.sh1.Briansuttonomyces.eucalypti~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU8.anov<-anova(OTU8.model,test = "Chisq")
# 
# OTU9.model<-glm.nb(formula =AbundNotZero.art1$RAE.sh12.Acrocalymma.vagum~SOIL*TISSUE+ HOST+ TIME+
#                      SITE,data = MetaRich.ART1,link = "log")
# OTU9.anov<-anova(OTU9.model,test = "Chisq")
# 
# OTU10.model<-glm.nb(formula =AbundNotZero.art1$PSE.wh14.Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU10.anov<-anova(OTU10.model,test = "Chisq")
# 
# OTU11.model<-glm.nb(formula =AbundNotZero.art1$PSE.wh40.Dimorphosporicola.tragani~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU11.anov<-anova(OTU11.model,test = "Chisq")
# 
# OTU12.model<-glm.nb(formula =AbundNotZero.art1$PSE.wh66.Comoclathris.italica~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU12.anov<-anova(OTU12.model,test = "Chisq")
# 
# OTU14.model<-glm.nb(formula =AbundNotZero.art1$PSE.we4..Chaetomium.globosum~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU14.anov<-anova(OTU14.model,test = "Chisq")
# 
# OTU15.model<-glm.nb(formula =AbundNotZero.art1$PSE.se8.Podospora.minicauda~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU15.anov<-anova(OTU15.model,test = "Chisq")
# 
# OTU16.model<-glm.nb(formula =AbundNotZero.art1$PSE.we8.Alternaria.chlamydospora~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU16.anov<-anova(OTU16.model,test = "Chisq")
# 
# OTU17.model<-glm.nb(formula =AbundNotZero.art1$LDE.se7.Coniothyrium.aleuritis~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU17.anov<-anova(OTU17.model,test = "Chisq")
# 
# OTU18.model<-glm.nb(formula =AbundNotZero.art1$LAE.se5.Sordaria.humana~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU18.anov<-anova(OTU18.model,test = "Chisq")
# 
# OTU19.model<-glm.nb(formula =AbundNotZero.art1$HAE.se5.Camarosporomyces.flavigenus~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU19.anov<-anova(OTU19.model,test = "Chisq")
# 
# OTU21.model<-glm.nb(formula =AbundNotZero.art1$HAE.we5.Coniolariella.sp.~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU21.anov<-anova(OTU21.model,test = "Chisq")
# 
# OTU22.model<-glm.nb(formula =AbundNotZero.art1$HAE.wh26.Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU22.anov<-anova(OTU22.model,test = "Chisq")
# 
# OTU23.model<-glm.nb(formula =AbundNotZero.art1$HAE.wh10.Raffaelea.montetyi~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU23.anov<-anova(OTU23.model,test = "Chisq")
# 
# OTU24.model<-glm.nb(formula =AbundNotZero.art1$HAE.wh65.Coniophora.marmorata~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU24.anov<-anova(OTU24.model,test = "Chisq")
# 
# OTU25.model<-glm.nb(formula =AbundNotZero.art1$HAE.se9.Chaetomium.nigricolor~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU25.anov<-anova(OTU25.model,test = "Chisq")
# 
# OTU26.model<-glm.nb(formula =AbundNotZero.art1$HAE.se1.Acrocalymma.sp.~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU26.anov<-anova(OTU26.model,test = "Chisq")
# 
# OTU27.model<-glm.nb(formula =AbundNotZero.art1$SREwh19.Neocamarosporium.goegapense~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU27.anov<-anova(OTU27.model,test = "Chisq")
# 
# OTU28.model<-glm.nb(formula =AbundNotZero.art1$APEsh6.Dictyosporium.digitatum~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU28.anov<-anova(OTU28.model,test = "Chisq")
# 
# OTU29.model<-glm.nb(formula =AbundNotZero.art1$APE.sh8.Pestalotiopsis.vismiae~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU29.anov<-anova(OTU29.model,test = "Chisq")
# 
# OTU30.model<-glm.nb(formula =AbundNotZero.art1$APE.se3.Dactylonectria.macrodidyma~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU30.anov<-anova(OTU30.model,test = "Chisq")
# 
# OTU31.model<-glm.nb(formula =AbundNotZero.art1$APE.sh5.Nigrospora.sphaerica~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU31.anov<-anova(OTU31.model,test = "Chisq")
# 
# OTU32.model<-glm.nb(formula =AbundNotZero.art1$SREwh18...Preussia.sp.~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU32.anov<-anova(OTU32.model,test = "Chisq")
# 
# OTU33.model<-glm.nb(formula =AbundNotZero.art1$LAEsh5.Coniolariella.ershadii~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1 = "log")
# OTU33.anov<-anova(OTU33.model,test = "Chisq")
# 
# OTU34.model<-glm(formula =AbundNotZero.art1$LAE.se3.Neosetophoma.lunariae~SOIL*TISSUE+ HOST+ TIME+
#                    SITE,data = MetaRich.ART1, family=poisson(link = "log"))
# OTU34.anov<-anova(OTU34.model,test = "Chisq")
# 
# OTU35.model<-glm.nb(formula =AbundNotZero.art1$LAE.SH.7.Muriphaeosphaeria.viburni~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU35.anov<-anova(OTU35.model,test = "Chisq")
# 
# OTU36.model<-glm.nb(formula =AbundNotZero.art1$LAE.sh1.Acrocalymma.sp.~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU36.anov<-anova(OTU36.model,test = "Chisq")
# 
# OTU37.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh30..Preussia.grandispora~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU37.anov<-anova(OTU37.model,test = "Chisq")
# 
# OTU38.model<-glm.nb(formula =AbundNotZero.art1$SRE.ss.4.Neocamarosporium.sp.~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU38.anov<-anova(OTU38.model,test = "Chisq")
# 
# OTU39.model<-glm.nb(formula =AbundNotZero.art1$SRE.ws8.Botryotrichum.murorum~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU39.anov<-anova(OTU39.model,test = "Chisq")
# 
# OTU40.model<-glm.nb(formula =AbundNotZero.art1$SRE.ws10.Sarocladium.kiliense~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU40.anov<-anova(OTU40.model,test = "Chisq")
# 
# OTU41.model<-glm.nb(formula =AbundNotZero.art1$SRE.wh16.Paracamarosporium.hawaiiense~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data =  MetaRich.ART1,link = "log")
# OTU41.anov<-anova(OTU41.model,test = "Chisq")
# 
# OTU42.model<-glm.nb(formula =AbundNotZero.art1$SRE.wh13.Ovatospora.senegalensis~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU42.anov<-anova(OTU42.model,test = "Chisq")
# 
# OTU43.model<-glm.nb(formula =AbundNotZero.art1$SRE.we6.Fusariella.sinensis~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU43.anov<-anova(OTU43.model,test = "Chisq")
# 
# OTU44.model<-glm.nb(formula =AbundNotZero.art1$SRE.we.10.pichia.kudriavzevii~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU44.anov<-anova(OTU44.model,test = "Chisq")
# 
# OTU45.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh7.Chaetomium.cucumericola~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU45.anov<-anova(OTU45.model,test = "Chisq")
# 
# OTU46.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh5.Fusarium.redolens~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU46.anov<-anova(OTU46.model,test = "Chisq")
# 
# OTU47.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh9.Preussia.intermedia~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU47.anov<-anova(OTU47.model,test = "Chisq")
# 
# OTU48.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh4.Penicillium.vinaceum~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU48.anov<-anova(OTU48.model,test = "Chisq")
# 
# OTU49.model<-glm.nb(formula =AbundNotZero.art1$SRE.sh3.Trichoderma.rifaii~SOIL*TISSUE+ HOST+ TIME+
#                       SITE,data = MetaRich.ART1,link = "log")
# OTU49.anov<-anova(OTU49.model,test = "Chisq")


########################################
########################################
# creat a Pie chart (Order)
OTU.colsum<-colSums(AbundNotZero.art1)

OTU.slic<- c(866,352,303,589,434,348,535,877,500,394,301, 4603)  #get the valus from OTU.colsum
OTU.lbls<- c("Rosellinia limonispora","Acrocalymma vagum","Dimorphosporicola tragani","Raffaelea montetyi",
               "Paracamarosporium hawaiiense","Fusariella sinensis","Humicola fuscoatra",
               "Neocamarosporium chichastianum", "Camarosporomyces flavigenus",
               "Preussia sp.","Coniophora marmorata", "Other")

OTU.Percent<-round(OTU.slic/sum (OTU.slic)*100, digits=2)
OTU.lbls <- paste(OTU.lbls, OTU.Percent)
order.lbls<-paste(OTU.lbls,"%",sep="")
pie(OTU.slic,labels =OTU.lbls, col = c("red","skyblue1","magenta",
                                           "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                           "yellow", "springgreen2", "pink","green" ) , main = "OTU", cex=1,border = NA,cex.main= 1.5, radius = 0.7)
########################################
######################################## 2 the most frequently isolated species
# 
# aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~HOST, data = Article1Meta, sum)
# aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~TISSUE, data = Article1Meta, sum)
# aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~SOIL, data = Article1Meta, sum)
# aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~SITE, data = Article1Meta, sum)
# aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~TIME, data = Article1Meta, sum)
# aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~TEMP, data = Article1Meta, sum)
# aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~MEDIA, data = Article1Meta, sum)
# 
# aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~HOST, data = Article1Meta, sum)
# aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~TISSUE, data = Article1Meta, sum)
# aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~SOIL, data = Article1Meta, sum)
# aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~SITE, data = Article1Meta, sum)
# aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~TIME, data = Article1Meta, sum)
# aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~TEMP, data = Article1Meta, sum)
# aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~MEDIA, data = Article1Meta, sum)

### MAHDIEH!!!!
# what are these numbers???? 
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

# >0 samples only
var.art2<-varpart(AbundNotZero.art1,MetaRich.ART1$SOIL,
                  MetaRich.ART1$HOST, MetaRich.ART1$SITE,MetaRich.ART1$TISSUE, 
                  transfo="hellinger")
plot(var.art2)
#without tissue
var.art3<-varpart(AbundNotZero.art1,MetaRich.ART1$season,
                  MetaRich.ART1$HOST, MetaRich.ART1$SITE,MetaRich.ART1$SOIL, 
                  transfo="hellinger")
plot(var.art3)
##varpart with models
var.part4<- varpart(AbundNotZero.art1,~MetaRich.ART1$TIME+
                      MetaRich.ART1$HOST+ MetaRich.ART1$SITE+ MetaRich.ART1$TISSUE,
                    ~env.Rich.ART1$Ece+env.Rich.ART1$pHe+env.Rich.ART1$Cle+env.Rich.ART1$EC,
                    transfo = "hellinger")
plot(var.part4)

#with soil*tissue
var.part5<-varpart(AbundNotZero.art1, MetaRich.ART1$HOST,MetaRich.ART1$season,
                   MetaRich.ART1$SITE,~MetaRich.ART1$SOIL*MetaRich.ART1$TISSUE, 
                   transfo="hellinger")
plot(var.part5,Xnames = c("Host","Time","Site","Soil*Organ"),
     bg=c("grey20","green","blue","red"),alpha=100)

###remove season
var.part6<-varpart(AbundNotZero.art1, MetaRich.ART1$HOST,
                   MetaRich.ART1$SITE,MetaRich.ART1$SOIL,MetaRich.ART1$TISSUE, 
                   transfo="hellinger")
plot(var.part6,Xnames = c("Host","Site","Soil","Organ"),
     bg=c("grey20","green","blue","red"),alpha=100)

# time and soil*tissue
var.part7<-varpart(AbundNotZero.art1, MetaRich.ART1$HOST,MetaRich.ART1$TIME,
                   MetaRich.ART1$SITE,~MetaRich.ART1$SOIL*MetaRich.ART1$TISSUE, 
                   transfo="hellinger")
plot(var.part7,Xnames = c("Host","Season","Site","Soil*Organ"),
     bg=c("green","blue","red","black"),alpha=100)

##############################################
##############################################
#soil data analysis
##############################################
##############################################
# data input
env.data<-read.csv("MataDataMergSoil.csv", header = T, row.names = 1)
view (env.data)

#subset for this article:
env.data.ar1<- subset (env.data, MetaData$HOST%in%c("Alhagi persarum","Artemisia sieberi", "Haloxylon ammodendron", 
                                                     "Launaea acunthodes",
                                                     "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                                     "Tamrix hispida"))
write.csv(env.data.ar1,file = "soildatagarmsar.csv")

# removed garmsar from file nd then read the new file
env.data.art1<-read.csv("soil.data.art1.csv",header = T)

# remove thesamples with zero obs for PCA
env.Rich.ART1 = env.data.art1[NotZero.art1,]
row.names(env.Rich.ART1)=row.names(MetaRich.ART1)
#remove soil and site variables
env.Rich.ART1$SITE<-NULL
env.Rich.ART1$SOIL<-NULL
env.Rich.ART1$Cle<-as.numeric(env.Rich.ART1$Cle)
##### PCA (Principal Components Analysis) analysis 
#use rda fun in vegan for PCA:
?rda()
rda.art1<-rda(AbundNotZero.art1~Ece + EC + pHe +Cle +  pH.1.5 + OM + OC + Nt + P + P2O5 + M + SP +
                M.SP,data =env.Rich.ART1)
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
#labs<-c("Ece","EC","Cle")
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

aggregate(glm.Lshoot)

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

Article1Order = subset (orderabund, MetaData$HOST %in% c("Alhagi persarum", "Artemisia sieberi", "Haloxylon ammodendron", 
                                                         "Launaea acunthodes",
                                                         "Prosopis stephaniana","Salsola incanescens","Seidlitzia rosmarinus",
                                                         "Tamrix hispida")&MetaData$SITE%in% c("Haj Ali Gholi Lake","Hoze Soltan Lake",
                                                                                               "Maranjab Desert","Rig-Boland Desert" ))
levels(Article1Meta$SITE)
# check to see if it worked 
rownames(Article1Order)==rownames(Article1Meta)
class(Article1Order)
colnames(Article1Order)
colSums(Article1Order)
# remove zero observed orders
Article1Order$Diaporthales<-NULL
View(Article1Order)
##################################################
####################################################
#Mahdie pls Check this! and make sure of the coloring

# creat a Pie chart (Order)
order.colsum<-colSums(Article1Order)
order.slic<- c(51, 589,301,1300,443,4764,88,1736,89,741)  #get the valus from order.colsum
order.lbls<- c("Eurotiales","Ophiostomatales","Boletales","Xylariales","Hypocreales"," Pleosporales",
               "Saccharomycetales","Sordariales ", "Amphisphaeriales","Unknown ")

order.Percent<-round(order.slic/sum(order.slic)*100, digits=2)
order.lbls <- paste(order.lbls, order.Percent)
order.lbls<-paste(order.lbls,"%",sep="")
pie(order.slic,labels =order.lbls, col = c("red","skyblue1","magenta",
                                           "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                           "yellow", "springgreen2") , main = "Order", cex=1,border = NA,cex.main= 1.5, radius = 0.7)


##############################################
################################################
### this plots needs to be re run
##### plot soil& host&orders
## first creat the data frames we need for the plot

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
eample1<-rbind (host.s.o,host.s.b,host.s.a, host.s.s, host.s.h, host.s.p, host.s.x, host.s.ss, host.s.u)
View(eample1 )

ggplot(eample1, aes(x=HOST, y=value, fill=Order)) +
  geom_bar(stat="identity") + facet_wrap(~SOIL) + theme_bw()

#### PLOT!
library(reshape)
eample1<-rbind (host.s.o,host.s.b,host.s.a, host.s.s, host.s.h, host.s.p, host.s.x, host.s.ss, host.s.u)
View(eample1 )
library(scales)
ggplot(eample1,aes(x = HOST, y = value,fill = Order)) + 
  geom_bar(position = "fill",stat = "identity") +facet_wrap(~SOIL) + theme_bw()+
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())+
  xlab("Host plant species")+ ylab("Proportional frequency")+
  labs(fill = "Order")
### lables should be fixed and readable


#########################################################
#########################################################
#########################################################
# phylogenic tree for the article


#### ggtree!
#install ggtree
#source("https://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade") # you may need this
#biocLite("ggtree")
library(ggtree)
library(tidyverse)

#loaded your tree
treeDP <- read.tree("A.tree")

ggplot(treeDP) + geom_tree() + theme_tree()
#Add a tree scale
ggtree(treeDP) + geom_treescale()
ggtree(treeDP) + theme_tree2()

#turn your tree into a cladogram
ggtree(treeDP, branch.length="none")

ggtree(treeDP, branch.length="none", color="red", size=1, linetype=1)

#turn your tree into circular layout
ggtree (treeDP, layout="circular") + ggtitle("(Phylogram) circular layout")

p <- ggtree(treeDP)
p + geom_nodepoint()
p + geom_tippoint()
p + geom_tiplab()
p + geom_tiplab()+ geom_nodepoint() 

ggtree(treeDP)+  geom_text(aes(label=node), hjust=-3)


# label and color for every branch
ggtree(treeDP) + geom_cladelabel(node=17, label="APEsh6", color="blue")

ggtree(treeDP) + geom_tiplab() +  geom_cladelabel(node=17, label="APEsh6", color="red", offset=1.5)

ggtree(treeDP) + geom_tiplab() + geom_cladelabel(node=17, label="Some random clade", color="red2", offset=1.5) + 
  geom_cladelabel(node=80, label="A different clade",  color="blue", offset=1.5)


ggtree(treeDP) + geom_tiplab() + geom_cladelabel(node=17, label="ABC", color="red2", offset=5, align=TRUE) + 
  geom_cladelabel(node=80, label="DEF", color="blue", offset=5, align=TRUE) + theme_tree2() + 
  xlim(0, 15) + theme_tree()

ggtree(treeDP) + geom_tiplab() + geom_hilight(node=17, fill="gold") + geom_hilight(node=80, fill="purple")

#Plot tree with other data
#read the help first!!!!!!!!!!!!!!!!
?facet_plot

#this shows the tip lables
d = fortify(treeDP)
d = subset(d, isTip)
d2<-with(d, label[order(y, decreasing=T)])

#extract and use in data for ploting
write.csv(d2, file = "tiplab.csv")

# modified data import
Treedata1 <- read.csv("Tree.csv")
View(Treedata1)

tree.p<-ggtree(treeDP)
# now plot toghether
tree.p2<-facet_plot(tree.p, panel='branch', data=Treedata1, geom=geom_segment, 
                    aes(x=0, xend=val, y=y, yend=y), size=3, color='blue')

tree.p3<-facet_plot(tree.p2, panel='leaf', data=Treedata1, geom=geom_segment, 
                    aes(x=0, xend=val, y=y, yend=y), size=3, color='green')
facet_plot(tree.p3, panel='root', data=Treedata1, geom=geom_segment, 
           aes(x=0, xend=val, y=y, yend=y), size=3, color='red')

### fixed the points size
facet_plot(tree.p, panel='data', data=Treedata1, geom=geom_point, 
           
           aes(x=0,size=val), color='blue')+ theme(legend.position="right")

#load the new tree
TREEARTICLE1 <- read.tree("Tree.nxs.tree")

ggplot(TREEARTICLE1) + geom_tree() + theme_tree()
#Add a tree scale
ggtree(TREEARTICLE1) + geom_treescale()

#turn your tree into a cladogram
ggtree(TREEARTICLE1, branch.length="none")

ggtree(TREEARTICLE1, branch.length="none", color="red", size=1, linetype=1)

#turn your tree into circular layout
ggtree (TREEARTICLE1, layout="circular") + ggtitle("(Phylogram) circular layout")

p <- ggtree(TREEARTICLE1)
p + geom_nodepoint()
p + geom_tippoint()
p + geom_tiplab()
p + geom_tiplab()+ geom_nodepoint() 



