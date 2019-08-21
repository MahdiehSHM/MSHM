
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

>>>>>>> c271edb1cb5f4adfc063c1c1864518a99de3425a
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
aggregate(shannon~season, data=MetaNotOne.art1,mean)
aggregate(simpson~TISSUE, data=MetaNotOne.art1,mean)


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

####################################################################
######## visualize the significant interaction results for the paper
####################################################################
### barplot
IR.interact<-aggregate(IR~SOIL+TISSUE,Article1Meta,mean)
MetaNotOne.art1$shannon<-shannon.art1
MetaNotOne.art1$simpson<-simpson.art1
shanon.interact<-aggregate(shannon~SOIL+TISSUE,MetaNotOne.art1,mean)
simpson.interact<-aggregate(simpson~SOIL+TISSUE,MetaNotOne.art1,mean)
inter.data<-cbind(IR.interact,shanon.interact,simpson.interact)
inter.data$SOIL<-NULL#run 2 times
inter.data$TISSUE<-NULL
#inter.data
#write.csv(inter.data, "interaction.csv")
# I manualy fixed the data fram
# now read it agin
interact.data<-read.csv("interaction.csv", header = TRUE)
levels(interact.data$Organ)
interact.data$Organ<-factor(interact.data$Organ, levels = c("Leaf" ,"Twig","Root" ))


# barplot
ggplot(data = interact.data,aes(x=Organ,y=value, fill= Index))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Soil)+ theme_bw()+ ylab("Mean value")+
  theme(legend.position="top",axis.text=element_text(size=12),
  axis.title=element_text(size=14),legend.text=element_text(size=12),
strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


####################################################
####################################################
# COMMUNITY COMPOSITION ANALYSIS
####################################################
####################################################
# for this we are using Permutational multivariate analysis of variance (PERMANOVA)
#this function from vegan does it
# memory.limit()
# memory.limit(size=30000)
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
OTU.organ<-colnames(mvabund.m.anova)[mvabund.m.anova["TISSUE",]<= 0.05]#40otus affected

###############################################
## plot affected OTUs by SOIL*TISSUE for the paper
################################################
# LAE.se5.Sordaria.humana data frame
LAEse5 <-aggregate(Article1OTU $  LAE.se5.Sordaria.humana  ~ TISSUE + SOIL , data = Article1Meta, sum)
LAEse5$value=LAEse5$`Article1OTU$LAE.se5.Sordaria.humana`
LAEse5$`Article1OTU$LAE.se5.Sordaria.humana`<-NULL
LAEse5$OTU<-c("LAE.se5.Sordaria.humana")

# TPEsh28.Humicola.fuscoatra data frame
TPEsh28 <-aggregate(Article1OTU $ TPEsh28.Humicola.fuscoatra  ~ TISSUE + SOIL , data = Article1Meta, sum)
TPEsh28$value=TPEsh28$`Article1OTU$TPEsh28.Humicola.fuscoatra`
TPEsh28$`Article1OTU$TPEsh28.Humicola.fuscoatra`<-NULL
TPEsh28$OTU<-c("TPEsh28.Humicola.fuscoatra")

# a new dataframe for ploting
SOILTISSUE<-rbind (TPEsh28,LAEse5)
<<<<<<< HEAD
View(SOILTISSUE )
write.csv(SOILTISSUE, file = "H.csv")
dev.off()
read.csv ("H.csv")

# the plot
ggplot(SOILTISSUE, aes(x=TISSUE, y=value, fill=OTU)) + xlab("Tissue")+ ylab("Frequency")+
  labs(fill = "OTU") + geom_bar(stat="identity") + facet_wrap(~SOIL) + theme_bw()

=======
# View(SOILTISSUE )
# write.csv(SOILTISSUE, file = "SOILTISSUE.csv")
# read.csv("SOILTISSUE.csv")
ST<- read.csv ("ST.csv")
ST$Organ<-factor(ST$Organ,levels = c("Leaf", "Twig","Root" ))
levels(ST$Organ)
##################### the soil:organ interaction plot for 2 otus
ggplot(ST, aes(x=Organ, y=Frequency, fill=OTU)) +
  geom_bar(stat="identity", width = 0.5) + facet_wrap(~Soil) + theme_bw()+
  xlab("Organ")+ ylab("Frequency")+ theme(legend.position="top",axis.text=element_text(size=12),
                                          axis.title=element_text(size=14),legend.text=element_text(size=12),
                                          strip.text = element_text(size = 12))+
  guides(fill=guide_legend(title=NULL))
>>>>>>> 7b58a5ad70b6ff45745f15e9d3f15190be98c65c

########################################
########################################
# Pie chart (OTUs!)

OTU.colsum<-colSums(AbundNotZero.art1)
# اعداد درصد ها رو از آبجکت بگیر و مثل این مثال اول وارد کن برای همه
OTU.slic<- c(866,352,303,589,434,348,535,877,500,394,301, 4603)  #get the valus from OTU.colsum
<<<<<<< HEAD
OTU.lbls<- c(expression(italic("Rosellinia sp.")),expression(italic("Acrocalymma vagum")),
                         expression(italic("Dimorphosporicola tragani")),
                         expression(italic("Raffaelea montetyi")),expression(italic("Paracamarosporium hawaiiense")),
                         expression(italic("Fusariella sinensis")),expression(italic("Humicola fuscoatra")),
                         expression(italic("Neocamarosporium chichastianum")),
                         expression(italic("Camarosporomyces flavigenus")),
                         expression(italic("Preussia sp.")),expression(italic("Coniophora marmorata")),"Other")
                         

OTU.Percent<-round(OTU.slic/sum (OTU.slic)*100, digits=2)
OTU.lbls <- paste(OTU.lbls, OTU.Percent)
order.lbls<-paste(OTU.lbls,"%",sep="")
pie(OTU.slic,labels =OTU.lbls, col = c("red","skyblue1","magenta",
                                           "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                           "yellow", "springgreen2", "pink","green" ) , main = "OTU", 
    cex=1,border = NA,cex.main= 0.5, radius = 1)

OTU.Percent<-round(OTU.slic/sum (OTU.slic)*100, digits=2)
OTU.lbls<- c (expression (italic("Rosellinia")*' sp. (8.57%)' ), expression (italic("Acrocalymma vagum")*' (3.48%)'), expression(italic("Dimorphosporicola")*' sp. (3%)'),
             expression(italic("Raffaelea montetyi")*' (5.83%)'), expression(italic("Paracamarosporium")*' sp. (4.3%)'), expression(italic("Fusariella")*' sp. (3.44%)'),
             expression(italic("Humicola fuscoatra")*' (5.3%)'),expression(italic("Neocamarosporium chichastianum")*' (8.68%)'),expression(italic("Chaetosphaeronema")*' sp. (4.95%)'),
             expression(italic("Preussia")*' sp. (3.9%)'),expression (italic("Coniophora")*' sp. (2.98%)'), "Other (45.57%)")

pie (OTU.slic,labels =OTU.lbls, col = c("red","skyblue1","magenta",
                                       "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                       "yellow", "springgreen2", "pink","green" ) , main = "OTUs frequency", cex=0.8,border = NA,cex.main= 1.1, radius =0.85)


########################################
######################################## Aggregate for ggtree
?par
aggregate(Article1OTU $ APE.se5.Staphylotrichum.coccosporum ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.se5.Staphylotrichum.coccosporum ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.se5.Staphylotrichum.coccosporum ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ TPEsh28.Humicola.fuscoatra ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ TPEsh28.Humicola.fuscoatra ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ TPEsh28.Humicola.fuscoatra ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ PFE.sh7..Rosellinia.limonispora ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ LREwh64..Neocamarosporium.chichastianum ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ SREwh22.Preussia.minimoides ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SREwh22.Preussia.minimoides ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SREwh22.Preussia.minimoides ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ THE.we10..Aporospora.terricola ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ THE.we10..Aporospora.terricola ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ THE.we10..Aporospora.terricola ~SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ THE.ss3.Fusarium.sp. ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ THE.ss3.Fusarium.sp. ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ THE.ss3.Fusarium.sp. ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ SIE.sh1.Briansuttonomyces.eucalypti ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SIE.sh1.Briansuttonomyces.eucalypti ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SIE.sh1.Briansuttonomyces.eucalypti ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ RAE.sh12.Acrocalymma.vagum ~HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ RAE.sh12.Acrocalymma.vagum ~TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ RAE.sh12.Acrocalymma.vagum ~SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ PSE.wh40.Dimorphosporicola.tragani ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.wh40.Dimorphosporicola.tragani ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.wh40.Dimorphosporicola.tragani ~ SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ PSE.wh66.Comoclathris.italica ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.wh66.Comoclathris.italica ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.wh66.Comoclathris.italica ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ PSE.we4..Chaetomium.globosum ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.we4..Chaetomium.globosum ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.we4..Chaetomium.globosum ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ PSE.we8.Alternaria.chlamydospora ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.we8.Alternaria.chlamydospora ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.we8.Alternaria.chlamydospora ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ LDE.se7.Coniothyrium.aleuritis ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ LDE.se7.Coniothyrium.aleuritis ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ LDE.se7.Coniothyrium.aleuritis ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ LAE.se5.Sordaria.humana ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.se5.Sordaria.humana ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.se5.Sordaria.humana ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ HAE.se5.Camarosporomyces.flavigenus ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.se5.Camarosporomyces.flavigenus ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.se5.Camarosporomyces.flavigenus ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ HAE.we5.Coniolariella.sp. ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.we5.Coniolariella.sp. ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.we5.Coniolariella.sp. ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ HAE.wh26.Preussia.sp. ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.wh26.Preussia.sp. ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.wh26.Preussia.sp. ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ HAE.wh10.Raffaelea.montetyi ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.wh10.Raffaelea.montetyi ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.wh10.Raffaelea.montetyi ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ HAE.wh65.Coniophora.marmorata ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.wh65.Coniophora.marmorata ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.wh65.Coniophora.marmorata ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ HAE.se9.Chaetomium.nigricolor ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.se9.Chaetomium.nigricolor ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.se9.Chaetomium.nigricolor ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ HAE.se1.Acrocalymma.sp. ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.se1.Acrocalymma.sp. ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ HAE.se1.Acrocalymma.sp. ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SREwh19.Neocamarosporium.goegapense ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SREwh19.Neocamarosporium.goegapense ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SREwh19.Neocamarosporium.goegapense ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ APEsh6.Dictyosporium.digitatum ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ APEsh6.Dictyosporium.digitatum ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ APEsh6.Dictyosporium.digitatum ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ APE.sh8.Pestalotiopsis.vismiae ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.sh8.Pestalotiopsis.vismiae ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.sh8.Pestalotiopsis.vismiae ~ SOIL, data = Article1Meta, sum) 

aggregate(Article1OTU $ APE.se3.Dactylonectria.macrodidyma ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.se3.Dactylonectria.macrodidyma ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.se3.Dactylonectria.macrodidyma ~ SOIL, data = Article1Meta, sum) 

aggregate(Article1OTU $ APE.sh5.Nigrospora.sphaerica ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.sh5.Nigrospora.sphaerica ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ APE.sh5.Nigrospora.sphaerica ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SREwh18...Preussia.sp. ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SREwh18...Preussia.sp. ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SREwh18...Preussia.sp. ~ SOIL, data = Article1Meta, sum) 

aggregate(Article1OTU $ LAEsh5.Coniolariella.ershadii ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ LAEsh5.Coniolariella.ershadii ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ LAEsh5.Coniolariella.ershadii ~ SOIL, data = Article1Meta, sum) 

aggregate(Article1OTU $ LAE.se3.Neosetophoma.lunariae ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.se3.Neosetophoma.lunariae ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.se3.Neosetophoma.lunariae ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ LAE.SH.7.Muriphaeosphaeria.viburni ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.SH.7.Muriphaeosphaeria.viburni ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.SH.7.Muriphaeosphaeria.viburni ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ LAE.sh1.Acrocalymma.sp. ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.sh1.Acrocalymma.sp. ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ LAE.sh1.Acrocalymma.sp. ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.sh30..Preussia.grandispora ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh30..Preussia.grandispora ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh30..Preussia.grandispora ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.ss.4.Neocamarosporium.sp. ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.ss.4.Neocamarosporium.sp. ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.ss.4.Neocamarosporium.sp. ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.ws8.Botryotrichum.murorum ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.ws8.Botryotrichum.murorum ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.ws8.Botryotrichum.murorum ~ SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ SRE.ws10.Sarocladium.kiliense ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.ws10.Sarocladium.kiliense ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.ws10.Sarocladium.kiliense ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.wh16.Paracamarosporium.hawaiiense ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.wh16.Paracamarosporium.hawaiiense ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.wh16.Paracamarosporium.hawaiiense ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.wh13.Ovatospora.senegalensis ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.wh13.Ovatospora.senegalensis ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.wh13.Ovatospora.senegalensis ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.we6.Fusariella.sinensis ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.we6.Fusariella.sinensis ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.we6.Fusariella.sinensis ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.we.10.pichia.kudriavzevii ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.we.10.pichia.kudriavzevii ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.we.10.pichia.kudriavzevii ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.sh9.Preussia.intermedia ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh9.Preussia.intermedia ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh9.Preussia.intermedia ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.sh4.Penicillium.vinaceum ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh4.Penicillium.vinaceum ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh4.Penicillium.vinaceum ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.sh3.Trichoderma.rifaii ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh3.Trichoderma.rifaii ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh3.Trichoderma.rifaii ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.sh7.Chaetomium.cucumericola ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh7.Chaetomium.cucumericola ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh7.Chaetomium.cucumericola ~ SOIL, data = Article1Meta, sum)

aggregate(Article1OTU $ SRE.sh5.Fusarium.redolens ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh5.Fusarium.redolens ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ SRE.sh5.Fusarium.redolens ~ SOIL, data = Article1Meta, sum)


aggregate(Article1OTU $ PSE.wh14.Preussia.sp. ~ HOST, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.wh14.Preussia.sp. ~ TISSUE, data = Article1Meta, sum)
aggregate(Article1OTU $ PSE.wh14.Preussia.sp. ~ SOIL, data = Article1Meta, sum)

PSE.wh14.Preussia.sp.

aggregate (IR~ TISSUE, data = Article1Meta, sum)

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
var.part7<-varpart(AbundNotZero.art1, MetaRich.ART1$HOST,
                   MetaRich.ART1$SITE,~MetaRich.ART1$SOIL*MetaRich.ART1$TISSUE, 
                   transfo="hellinger")
plot(var.part7,Xnames = c("Host","Site","Soil*Organ"),
     bg=c("green","blue","red"),alpha=150)

##############################################
##############################################
#soil data analysis
##############################################
##############################################
# data input
env.data<-read.csv("MataDataMergSoil.csv", header = T, row.names = 1)


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
# plot(rda.art1, xlim=c(-2,2), ylim=c(-1,2.5))# this is a simple plot with both sites and species
# 
# # #coloring:
# # soilfactor<-factor(MetaRich.ART1$SOIL)
# # colvec <-  c("red","blue")
# # cols<-colvec[soilfactor]
# # #scors:
# rda.scores <- scores(rda.art1, display = 'bp')
# mul <- ordiArrowMul(rda.scores, fill = 0.75)
# 
# dev.off()
# ####### FINAL RDA PLOT
# plot(rda.art1, type = "n")
# points(rda.art1, display = "sites", col = colvec[soilfactor], 
#        pch = (16:17)[soilfactor],cex=0.85)#salin=blue
# #points(rda.art1, display = "species", pch = "+")
# arrows(0, 0, mul * rda.scores[,1], mul * rda.scores[,2],
#        length = 0.09, col ="black", lwd=2)
# labs <- rownames(rda.scores)
# #labs<-c("Ece","EC","Cle")
# text(ordiArrowTextXY(mul * rda.scores, labs), labs)
# legend("topleft", c("Arid soil","Saline soil"), 
#        col=c("red","blue"),
#        pch = c(16,17), border="white", bty="n")


# new plot for the paper

#coloring:
soilfactor<-factor(MetaRich.ART1$SOIL)
colvec <-  c("red","blue")
cols<-colvec[soilfactor]

rda.art1 %>% plot(type="n", corr=TRUE) %>% 
  points("sites", pch=(16:17)[soilfactor], col=colvec[soilfactor]) %>% 
  text( "biplot", arrows=TRUE, col= "black", lwd=2, cex=1.2, xpd=TRUE, len=0.05)
legend("topleft", c("Dry soil","Saline soil"), 
       col=c("red","blue"),
       pch = c(16,17), border="white", bty="n",cex = 1.2)
dev.off()
##############################################
##############################################
# Greenhouse data analysis
##############################################
##############################################
# data input
GH.data<-read.csv("Green house data.csv", header = T, row.names = 1)
# fix order of levels of variable
levels(GH.data$Drought)
GH.data$Drought<-factor(GH.data$Drought, levels = c("25%","50%","75%","100%"))

levels(GH.data$Salinity)
GH.data$Salinity<-factor(GH.data$Salinity, levels = c( "8 dS/m","12 dS/m","16 dS/m"))

GH.data$Biomass<-GH.data$DWshoot+GH.data$DWroot

# 3-way anova
# factors: Fungi, Drought, Salinity
# ?aov
# shoot lenght
GH.lshoot<-aov(Lshoot~Fungi*Drought*Salinity,data = GH.data)
ghlshoot.sum<-summary(GH.lshoot)
# shoot wet weight
GH.Wshoot<-aov(Wshoot~Fungi*Drought*Salinity,data = GH.data)
ghWshoot.sum<-summary(GH.Wshoot)
# shoot dry weight
GH.DWshoot<-aov(DWshoot~Fungi*Drought*Salinity,data = GH.data)
ghDWshoot.sum<-summary(GH.DWshoot)
# root lenght
GH.Lroot<-aov(Lroot~Fungi*Drought*Salinity,data = GH.data)
ghLroot.sum<-summary(GH.Lroot)
# root wet weight
GH.Wroot<-aov(Wroot~Fungi*Drought*Salinity,data = GH.data)
ghWroot.sum<-summary(GH.Wroot)
# root dry weight
GH.DWroot<-aov(DWroot~Fungi*Drought*Salinity,data = GH.data)
ghDWroot.sum<-summary(GH.DWroot)
#Chlorophyll concentration
GH.photo<-aov(Photosyntesis~Fungi*Drought*Salinity,data = GH.data)
ghphoto.sum<-summary(GH.photo)
# total biomass
GH.Biomass<-aov(Biomass~Fungi*Drought*Salinity,data = GH.data)
ghBiomass.sum<-summary(GH.Biomass)
################################
# plots for the paper
#################################
# how to show this?
# only the significant interactions
# shoot lenght
ggplot(GH.data, aes(x=Drought,y=Lshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Shoot lenght (cm)")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# shoot wet weight
ggplot(GH.data, aes(x=Drought,y=Wshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+
  ylab("Shoot wet weight (g)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# shoot dry weight
ggplot(GH.data, aes(x=Drought,y=DWshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+
  ylab("Shoot dry weight (g)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# root lenght
ggplot(GH.data, aes(x=Drought,y=Lroot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Root lenght (cm)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# root wet weight
ggplot(GH.data, aes(x=Drought,y=Wroot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+
  ylab("Root wet weight (g)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))


# Chlorophyll concentration
ggplot(GH.data, aes(x=Drought,y=Photosyntesis, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()
# total biomass
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
png("orderplot1.png",width = 1700, height = 1700,res=600)
order.colsum<-colSums(Article1Order)
order.slic<- c(51, 589,301,1300,443,4764,88,1736,89,741)  #get the valus from order.colsum
Order.Percent<- round(order.slic/sum(OTU.slic)*100, digits = 2 )
order.lbls<- c("Eurotiales (0.5%)","Ophiostomatales (5.83%)","Boletales (2.98%)","Xylariales (12.87%)",
               "Hypocreales (4.39%)"," Pleosporales (47.16%)",
               "Saccharomycetales (0.87%)","Sordariales (17.18%) ", "Amphisphaeriales (0.88%)","Unknown (7.34%)")

<<<<<<< HEAD
order.Percent<-round(order.slic/sum(order.slic)*100, digits=2)
order.lbls <- paste(order.lbls, order.Percent)
order.lbls<-paste(order.lbls,"%",sep="")
pie(order.slic,labels =order.lbls, col = c("red","skyblue1","magenta",
                                           "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                           "yellow", "springgreen2") , main = "Order", cex=0.2,border = NA,cex.main= 0.3, radius = 0.7)
dev.off()
=======
pie(order.slic,labels =order.lbls, col = c("red","yellowgreen","dodgerblue",
                                           "plum2","blue","orchid4","yellow","tomato2",
                                           "goldenrod", "aquamarine4") , main = "Orders frequency", 
    cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)

>>>>>>> 7b58a5ad70b6ff45745f15e9d3f15190be98c65c

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
# eample1<-rbind (host.s.o,host.s.b,host.s.a, host.s.s, host.s.h, host.s.p, host.s.x, host.s.ss, host.s.u)
# View(eample1 )
# 
# ggplot(eample1, aes(x=HOST, y=value, fill=Order)) +
#   geom_bar(stat="identity") + facet_wrap(~SOIL) + theme_bw()

#### PLOT!

eample1<-rbind (host.s.o,host.s.b,host.s.a, host.s.s, host.s.h, host.s.p, host.s.x, host.s.ss, host.s.u)
View(eample1 )

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

# biocLite("ggtree")
library(ggtree)
# BiocManager::install("ggtree")
library(tidyverse)
library(tidytree)

#biocLite("ggtree")

#loaded your tree

# treeDP1 <- read.tree("Final.nxs.tree")

# ggplot(treeDP1) + geom_tree() + theme_tree()
#Add a tree scale
# ggtree(treeDP1) + geom_treescale()
# ggtree(treeDP1) + theme_tree2()

#turn your tree into a cladogram
# ggtree(treeDP, branch.length="none")

# ggtree(treeDP, branch.length="none", color="red", size=1, linetype=1)

#turn your tree into circular layout
# ggtree (treeDP, layout="circular") + ggtitle("(Phylogram) circular layout")

# p <- ggtree(treeDP)
# p + geom_nodepoint()
# p + geom_tippoint()
# p + geom_tiplab()
# p + geom_tiplab()+ geom_nodepoint() 

# ggtree(treeDP)+  geom_text(aes(label=node), hjust=-3)


# label and color for every branch
# ggtree(treeDP1) + geom_cladelabel(node=17, label="APEsh6", color="blue")
# 
# ggtree(treeDP) + geom_tiplab() +  geom_cladelabel(node=17, label="APEsh6", color="red", offset=1.5)
# 
# ggtree(treeDP) + geom_tiplab() + geom_cladelabel(node=17, label="Some random clade", color="red2", offset=1.5) + 
#   geom_cladelabel(node=80, label="A different clade",  color="blue", offset=1.5)
# 
# 
# ggtree(treeDP) + geom_tiplab() + geom_cladelabel(node=17, label="ABC", color="red2", offset=5, align=TRUE) + 
#   geom_cladelabel(node=80, label="DEF", color="blue", offset=5, align=TRUE) + theme_tree2() + 
#   xlim(0, 15) + theme_tree()
# 
# ggtree(treeDP) + geom_tiplab() + geom_hilight(node=17, fill="gold") + geom_hilight(node=80, fill="purple")

#Plot tree with other data
#read the help first!!!!!!!!!!!!!!!!
?facet_plot

# s <- read.tree("AfterEdit.tree")
# ggtree(tree, ladderize = FALSE)
# ggplot(s) + geom_tree() + theme_tree()
# #Add a tree scale
# ggtree(s) + geom_treescale()
# ggtree(s) + theme_tree2()
# 
# #turn your tree into a cladogram
# ggtree(s, branch.length="none",size=0.5, linetype=1, hjust=-3) + geom_tiplab()+ geom_nodepoint() 
# 
# ggtree(s, branch.length="none", color="red", size=1, linetype=1)
# 
# #turn your tree into circular layout
# ggtree (s, layout="circular") + ggtitle("(Phylogram) circular layout")
# 
# p <- ggtree(s)
# p + geom_nodepoint()
# p + geom_tippoint()
# p + geom_tiplab()
# p + geom_tiplab()+ geom_nodepoint() 
# 
# ggtree(s)+  geom_text(aes(label=node), hjust=-3)
# 
# 
# # label and color for every branch
# ggtree(s) + geom_cladelabel(node=17, label="APEsh6", color="blue")
# 
# ggtree(s) + geom_tiplab() +  geom_cladelabel(node=17, label="APEsh6", color="red", offset=1.5)
# 
# ggtree(s) + geom_tiplab() + geom_cladelabel(node=17, label="Some random clade", color="red2", offset=1.5) + 
#   geom_cladelabel(node=80, label="A different clade",  color="blue", offset=1.5)
# 
# 
# ggtree(s) + geom_tiplab() + geom_cladelabel(node=17, label="ABC", color="red2", offset=5, align=TRUE) + 
#   geom_cladelabel(node=80, label="DEF", color="blue", offset=5, align=TRUE) + theme_tree2() + 
#   xlim(0, 15) + theme_tree()
# 
# ggtree(s) + geom_tiplab() + geom_hilight(node=17, fill="gold") + geom_hilight(node=80, fill="purple")
# 
# #Plot tree with other data
# #read the help first!!!!!!!!!!!!!!!!
# ?facet_plot
# 
# #this shows the tip lables
# d = fortify(treeDP)
# d = subset(d, isTip)
# d2<-with(d, label[order(y, decreasing=T)])
# 
# #extract and use in data for ploting
# write.csv(d2, file = "tiplab.csv")


############# THE TREE!

# modified data import

#Treedata1 <- read.csv("Tree.csv")
#View(Treedata1)

#tree.p<-ggtree(treeDP1)
# # now plot toghether
# tree.p2<-facet_plot(tree.p, panel='branch', data=Treedata1, geom=geom_segment, 
#                     aes(x=0, xend=val, y=y, yend=y), size=3, color='blue')
# 
# tree.p3<-facet_plot(tree.p2, panel='leaf', data=Treedata1, geom=geom_segment, 
#                     aes(x=0, xend=val, y=y, yend=y), size=3, color='green')
# facet_plot(tree.p3, panel='root', data=Treedata1, geom=geom_segment, 
#            aes(x=0, xend=val, y=y, yend=y), size=3, color='red')
# 
# ### fixed the points size
# facet_plot(tree.p, panel='data', data=Treedata1, geom=geom_point, 
#            
#            aes(x=0,size=val), color='blue')+ theme(legend.position="right")

# #load the new tree
# TREEARTICLE1 <- read.tree("Tree.nxs.tree")
# 
# ggplot(TREEARTICLE1) + geom_tree() + theme_tree()
# #Add a tree scale
# ggtree(TREEARTICLE1) + geom_treescale()
# 
# #turn your tree into a cladogram
# ggtree(TREEARTICLE1, branch.length="none")
# 
# ggtree(TREEARTICLE1, branch.length="none", color="red", size=1, linetype=1)
# 
# #turn your tree into circular layout
# ggtree (TREEARTICLE1, layout="circular") + ggtitle("(Phylogram) circular layout")
# 
# p <- ggtree(TREEARTICLE1)
# p + geom_nodepoint()
# p + geom_tippoint()
# p + geom_tiplab()
# p + geom_tiplab()+ geom_nodepoint() 


#############################################################
#Fig 4- article 
#############################################################

############# THE TREE!
library(ggplot2)
library(ggtree)
library(gtable)
library(grid)
library(tidytree)

# import the tree 
treetmu<-read.tree("BeforeEdit.tree")
############
#get the tip lables from the tree
# d = fortify(treetmu)
# d = subset(d, isTip)
# d2<-with(d, label[order(y, decreasing=T)])
# #extract and use in data for ploting
# write.csv(d2, file = "tiplabtmu.csv")
# dev.off()

# data import for ploting
Treedata.tmu<- read.csv("tiplabtmu.csv")
View(Treedata.tmu)

####plotig the tree

tree.tmua<-ggtree(treetmu,branch.length= "none")
tree.tmu1<-tree.tmua + geom_tiplab()

# now plot toghether with soil data
tree.tmu2<-facet_plot(tree.tmu1+xlim_tree(50), panel='Dry soil', data=Treedata.tmu, geom=geom_point, 
                    aes(x=0, size=Dry), color='blue')

tree.tmu3<-facet_plot(tree.tmu2, panel='Saline soil', data=Treedata.tmu, geom=geom_point,
                    aes(x=0, size=Saline), color='blue')+ theme(legend.position="right",legend.title = element_blank())
# fixig the pael size for better visualization
#panel 1
gtree1$layout$l[grep("panel-1",gtree1$layout$name)]
gtree1$widths[5]=2.5*gtree1$widths[5]
grid.draw(gtree1)
#panel 2
gtree1<-ggplot_gtable(ggplot_build(tree.tmu3))
gtree1$layout$l[grep("panel-2",gtree1$layout$name)]
gtree1$widths[7]=0.3*gtree1$widths[7]
grid.draw(gtree1)
#panel 3
gtree1$layout$l[grep("panel-3",gtree1$layout$name)]
gtree1$widths[9]=0.3*gtree1$widths[9]
##########################################
####plot final tree with lables
grid.draw(gtree1)
# export file
png("tree-Fig 4.png", height = 3000, width = 2500, res= 300)
grid.draw(gtree1)
dev.off()

# try to change the colors of different branch based on soils: saline:yellowe, Dry: red, both: green
#new data set for coloring of tiplabs
color.tiplab<-read.csv("tiplabsoils.csv")

tree.tmua %<+% color.tiplab + 
  geom_tiplab(aes(fill = factor(soilfactor)) , color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.08, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  theme(legend.position = "right", 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
?geom_tiplab()
?`ggtree-package`
?geom_cladelabel
?geom_hilight
# label and color for every branch
# ggtree(s) + geom_cladelabel(node=17, label="APEsh6", color="blue")
# 
# ggtree(s) + geom_tiplab() +  geom_cladelabel(node=17, label="APEsh6", color="red", offset=1.5)
# 
# ggtree(s) + geom_tiplab() + geom_cladelabel(node=17, label="Some random clade", color="red2", offset=1.5) + 
#   geom_cladelabel(node=80, label="A different clade",  color="blue", offset=1.5)
# 
# 
# ggtree(s) + geom_tiplab() + geom_cladelabel(node=17, label="ABC", color="red2", offset=5, align=TRUE) + 
#   geom_cladelabel(node=80, label="DEF", color="blue", offset=5, align=TRUE) + theme_tree2() + 
#   xlim(0, 15) + theme_tree()
# 
# ggtree(s) + geom_tiplab() + geom_hilight(node=17, fill="gold") + geom_hilight(node=80, fill="purple")

####################################

##### fig 1-article



### data prep

otu.sum.organ<-aggregate(.~MetaRich.ART1$TISSUE,AbundNotZero.art1,sum)



otu.count<-colSums(AbundNotZero.art1)



OTU.250<-subset(otu.count, otu.count>250)

class(OTU.250)

OTU.250<-as.data.frame(OTU.250)

OTU.250<-t(OTU.250)

# 15 OTUs selected



selectedOTU<-colnames(OTU.250)



#subset the sum aggregate

selected.otu.count.organ<-otu.sum.organ[,names(otu.sum.organ) %in% colnames(OTU.250)]

#add the organ column to it

selected.otu.count.organ$organ<-otu.sum.organ$`MetaRich.ART1$TISSUE`

View(selected.otu.count.organ)

# now you have a data fram of 15 selected OTUs with more than 250 observations in selected.otu.count.organ

# total observation per organ?

total.obs.organ<-aggregate(IR~TISSUE,MetaRich.ART1,sum)



#######################

#leaf OTU pie chart 

leaf.pie.otu<- subset(selected.otu.count.organ, selected.otu.count.organ$organ=="Leaf")

leaf.pie.otu$organ<-NULL

leaf.pie.otu$Others<-1432

# more than zero?

leaf.pie.otu.0<- leaf.pie.otu[, colSums(leaf.pie.otu != 0) > 0]

lables.pie.leaf<-colnames(leaf.pie.otu.0)

numbers.pie.leaf<-as.numeric(leaf.pie.otu.0[1,])

3245-1813



#plot the leaf pie chart

p.leaf<-pie(numbers.pie.leaf,labels =lables.pie.leaf, col = c("chartreuse4","yellowgreen","dodgerblue2",
                                                              
                                                              "yellow1","orchid4","maroon2","darkcyan","plum2",
                                                              
                                                              "cyan1", "goldenrod","tomato2") , main = "Leaf", 
            
            cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)

#######################

#Twig OTU pie chart 

Twig.pie.otu<- subset(selected.otu.count.organ, selected.otu.count.organ$organ=="Branch")

Twig.pie.otu$organ<-NULL

Twig.pie.otu$Others<-985

# more than zero?

Twig.pie.otu.0<- Twig.pie.otu[, colSums(Twig.pie.otu != 0) > 0]

lables.pie.Twig<-colnames(Twig.pie.otu.0)

numbers.pie.Twig<-as.numeric(Twig.pie.otu.0[1,])

2533-1548



#plot the Twig pie chart

p.twig<-pie(numbers.pie.Twig,labels =lables.pie.Twig, col = c("chartreuse4","yellowgreen","dodgerblue2","yellow1",
                                                              
                                                              "orchid4","darkcyan","plum2","cyan1","slategray",
                                                              
                                                              "goldenrod","tomato2") , main = "Twig", 
            
            cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)



#######################

#Root OTU pie chart

Root.pie.otu<- subset(selected.otu.count.organ, selected.otu.count.organ$organ=="Root")

Root.pie.otu$organ<-NULL

Root.pie.otu$Others<-1129

# more than zero?

Root.pie.otu.0<- Root.pie.otu[, colSums(Root.pie.otu != 0) > 0]

lables.pie.Root<-colnames(Root.pie.otu.0)

numbers.pie.Root<-as.numeric(Root.pie.otu.0[1,])

4324-3195



#plot the Root pie chart

p.root<-pie(numbers.pie.Root,labels =lables.pie.Root, col = c("chartreuse4","Violetred4","darkturquoise",
                                                              
                                                              "yellowgreen","plum4","yellow1","cyan1","darkorange",
                                                              
                                                              "goldenrod","tomato2") , main = "Root", 
            
            cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)

##############

# Orders pie chart

p.order<-pie(order.slic,labels =order.lbls, col = c("red","yellowgreen","dodgerblue",
                                                    
                                                    "plum2","blue","orchid4","yellow","tomato2",
                                                    
                                                    "goldenrod", "aquamarine4") , main = "Orders frequency", 
             
             cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)



### host otu barplot

otu.sum.Host<-aggregate(.~MetaRich.ART1$HOST,AbundNotZero.art1,sum)

#subset the sum aggregate

selected.otu.count.Host<-otu.sum.Host[,names(otu.sum.Host) %in% colnames(OTU.250)]



#add the Host column to it

selected.otu.count.Host$Host<-otu.sum.Host$`MetaRich.ART1$HOST`

View(selected.otu.count.Host)

# total observation per Host?

total.obs.Host<-aggregate(IR~HOST,MetaRich.ART1,sum)

#manually prep the data frame for the plot

#write.csv(selected.otu.count.Host,"selected.otu.host.csv")

# import new data

host.otu.data<-read.csv("selected.otu.host.csv", header = TRUE)

otu.cols<-c("plum4","orchid4","yellowgreen","yellow1","maroon2","darkcyan","slategray","dodgerblue2",
            
            "goldenrod","chartreuse4","darkturquoise","tomato2",'darkorange',"plum2","cyan1","Violetred4")



#barplot

host.barplot<-ggplot(host.otu.data,aes(x = Host, y = Frequency,fill = OTU)) + 
  
  geom_bar(position = "fill",stat = "identity", width = 0.5)+ theme_bw()+scale_y_continuous(labels = percent_format())+
  
  xlab("Host plant species")+ ylab("Proportional frequency")+
  
  labs(fill = "OTU")+ scale_fill_manual(values=otu.cols)+theme(legend.position="top")

### export all of the plots together  

margin1<-matrix(c(0,2,0,0,2,0,1,3,5,1,3,5,0,4,0,0,4,0))  

layout(margin1)

dev.off()

par(mfrow=c(3,2))

p.leaf

p.twig

p.root

p.order

host.barplot




Treedata.final <- read.csv("finalTREEdata.csv")
View(Treedata.final)
# import the tree
tree1 <- read.tree("FinalTree.tree")
tree.p<-ggtree(tree1,branch.length= "none")

# now plot toghether with soil data
tree.p2<-facet_plot(tree.p, panel='Dry soil', data=Treedata.final, geom=geom_point, 
                    aes(x=0, size=Dry), color='blue')

tree.p3<-facet_plot(tree.p2, panel='Saline soil', data=Treedata.final, geom=geom_point, 
                    aes(x=0, size=Saline), color='blue')+
  theme(legend.position="right",legend.title = element_blank())
  
####################################
##### fig 1-article
###################################

### data prep
otu.sum.organ<-aggregate(.~MetaRich.ART1$TISSUE,AbundNotZero.art1,sum)

otu.count<-colSums(AbundNotZero.art1)

OTU.250<-subset(otu.count, otu.count>250)
class(OTU.250)
OTU.250<-as.data.frame(OTU.250)
OTU.250<-t(OTU.250)
# 15 OTUs selected

selectedOTU<-colnames(OTU.250)

#subset the sum aggregate
selected.otu.count.organ<-otu.sum.organ[,names(otu.sum.organ) %in% colnames(OTU.250)]
#add the organ column to it
selected.otu.count.organ$organ<-otu.sum.organ$`MetaRich.ART1$TISSUE`
View(selected.otu.count.organ)
# now you have a data fram of 15 selected OTUs with more than 250 observations in selected.otu.count.organ
# total observation per organ?
total.obs.organ<-aggregate(IR~TISSUE,MetaRich.ART1,sum)

#######################
#leaf OTU pie chart 
leaf.pie.otu<- subset(selected.otu.count.organ, selected.otu.count.organ$organ=="Leaf")
leaf.pie.otu$organ<-NULL
leaf.pie.otu$Others<-1432
# more than zero?
leaf.pie.otu.0<- leaf.pie.otu[, colSums(leaf.pie.otu != 0) > 0]
lables.pie.leaf<-colnames(leaf.pie.otu.0)
numbers.pie.leaf<-as.numeric(leaf.pie.otu.0[1,])
3245-1813

#plot the leaf pie chart
p.leaf<-pie(numbers.pie.leaf,labels =lables.pie.leaf, col = c("chartreuse4","yellowgreen","dodgerblue2",
                                           "yellow1","orchid4","maroon2","darkcyan","plum2",
                                           "cyan1", "goldenrod","tomato2") , main = "Leaf", 
    cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)
#without lables
p.leaf<-pie(numbers.pie.leaf,labels ="", col = c("chartreuse4","yellowgreen","dodgerblue2",
                                                              "yellow1","orchid4","maroon2","darkcyan","plum2",
                                                              "cyan1", "goldenrod","tomato2") , main = "Leaf (b)", 
            cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)

#######################
#Twig OTU pie chart 
Twig.pie.otu<- subset(selected.otu.count.organ, selected.otu.count.organ$organ=="Branch")
Twig.pie.otu$organ<-NULL
Twig.pie.otu$Others<-985
# more than zero?
Twig.pie.otu.0<- Twig.pie.otu[, colSums(Twig.pie.otu != 0) > 0]
lables.pie.Twig<-colnames(Twig.pie.otu.0)
numbers.pie.Twig<-as.numeric(Twig.pie.otu.0[1,])
2533-1548

#plot the Twig pie chart
p.twig<-pie(numbers.pie.Twig,labels =lables.pie.Twig, col = c("chartreuse4","yellowgreen","dodgerblue2","yellow1",
                                                      "orchid4","darkcyan","plum2","cyan1","slategray",
                                                      "goldenrod","tomato2") , main = "Twig", 
    cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)
#without the labs
p.twig<-pie(numbers.pie.Twig,labels ="", col = c("chartreuse4","yellowgreen","dodgerblue2","yellow1",
                                                              "orchid4","darkcyan","plum2","cyan1","slategray",
                                                              "goldenrod","tomato2") , main = "Twig (c)", 
            cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)

#######################
#Root OTU pie chart
Root.pie.otu<- subset(selected.otu.count.organ, selected.otu.count.organ$organ=="Root")
Root.pie.otu$organ<-NULL
Root.pie.otu$Others<-1129
# more than zero?
Root.pie.otu.0<- Root.pie.otu[, colSums(Root.pie.otu != 0) > 0]
lables.pie.Root<-colnames(Root.pie.otu.0)
numbers.pie.Root<-as.numeric(Root.pie.otu.0[1,])
4324-3195

#plot the Root pie chart
p.root<-pie(numbers.pie.Root,labels =lables.pie.Root, col = c("chartreuse4","Violetred4","darkturquoise",
                                                      "yellowgreen","plum4","yellow1","cyan1","darkorange",
                                                      "goldenrod","tomato2") , main = "Root", 
    cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)
#without the labs
p.root<-pie(numbers.pie.Root,labels ="", col = c("chartreuse4","Violetred4","darkturquoise",
                                                              "yellowgreen","plum4","yellow1","cyan1","darkorange",
                                                              "goldenrod","tomato2") , main = "Root (d)", 
            cex=0.8,border = NA,cex.main= 1.1, radius = 0.9)

##############
# Orders pie chart
p.order<-pie(order.slic,labels =order.lbls, col = c("red","yellowgreen","dodgerblue",
                                                    "plum2","blue","orchid4","yellow","tomato2",
                                                    "goldenrod", "aquamarine4") , main = "Taxonomic Orders (a)", 
             cex=1.2,border = NA,cex.main= 1.1, radius = 0.9)

### host otu barplot
otu.sum.Host<-aggregate(.~MetaRich.ART1$HOST,AbundNotZero.art1,sum)
#subset the sum aggregate
selected.otu.count.Host<-otu.sum.Host[,names(otu.sum.Host) %in% colnames(OTU.250)]

#add the Host column to it
selected.otu.count.Host$Host<-otu.sum.Host$`MetaRich.ART1$HOST`
View(selected.otu.count.Host)
# total observation per Host?
total.obs.Host<-aggregate(IR~HOST,MetaRich.ART1,sum)
#manually prep the data frame for the plot
#write.csv(selected.otu.count.Host,"selected.otu.host.csv")
# import new data
host.otu.data<-read.csv("selected.otu.host.csv", header = TRUE)
otu.cols<-c("plum4","orchid4","yellowgreen","yellow1","maroon2","darkcyan","slategray","dodgerblue2",
"goldenrod","chartreuse4","darkturquoise","tomato2",'darkorange',"plum2","cyan1","Violetred4")
<<<<<<< HEAD

#barplot
host.barplot<-ggplot(host.otu.data,aes(x = Host, y = Frequency,fill = OTU)) + 
  geom_bar(position = "fill",stat = "identity", width = 0.5)+ theme_bw()+scale_y_continuous(labels = percent_format())+
  xlab("Host plant species (e)")+ ylab("Proportional frequency")+
  labs(fill = "OTU")+ scale_fill_manual(values=otu.cols)+
  theme(legend.position="top",legend.text=element_text(size=16),
        legend.title=element_text(size=16),axis.text=element_text(size=14,face = "italic"),
        axis.title=element_text(size=14,face="bold"))


#########################
#### Fig 2- article
#########################
# barplot
ggplot(data = interact.data,aes(x=Organ,y=value, fill= Index))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Soil)+ theme_bw()+ ylab("Mean value")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))



#############################
##Fig 3- article
#############################
dev.off()
# NMDS Hosts
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5), cex.axis = 1.5, cex.lab = 1.5)
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.8)
ordiellipse(nmds.art2, MetaRich.ART1$HOST,cex=1,alpha = 200, 
            draw="polygon", col= 1:8,border=1:8,lwd=3, kind="se", conf=0.95)

legend("bottomright", c("A.pers" ,"A.sieb", "H.ammo" ,"L.acun" ,"P.step" ,"S.inca", "S.rosm" ,"T.hisp"), 
       fill= 1:8,border="white", bty="n",title = "a", cex = 1.5)

#NMDS ORGAN:
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5), cex.axis = 1.5, cex.lab = 1.5)
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$TISSUE,cex=1,
            draw="polygon", col= c("blue","green","red"),border=c("blue","green","red"),lwd=3,alpha = 200, kind="se", conf=0.95)

legend("bottomright", c("Branch","Leaf","Root" ), 
       fill= c("blue","green","red"),border="white", bty="n",title = "b", cex = 1.5)

#NMDS sampling sites:
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5), cex.axis = 1.5, cex.lab = 1.5)
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$SITE,cex=1,alpha = 200, 
            draw="polygon", col= 1:5,border= 1:5,lwd=3, kind="se", conf=0.95)

legend("bottomright", c("Haj Ali Gholi Lake","Hoze Soltan Lake","Maranjab Desert","Rig-Boland Desert"), 
       fill= 1:5, border="white", bty="n", title = "c", cex = 1.5)

# NMDS soil:
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5), cex.axis = 1.5, cex.lab = 1.5)
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$SOIL,cex=1,alpha = 200, 
            draw="polygon", col= "red",lwd = 3,border="red",
            kind="se", conf=0.95,show.groups=(c("Arid soil")))
ordiellipse(nmds.art2, MetaRich.ART1$SOIL,cex=1,lwd = 3,alpha = 200, 
            draw="polygon", col= "Blue",border="Blue",
            kind="se", conf=0.95,show.groups=(c("Saline Soil")))
legend("bottomright", c("Saline soil","Arid soil"), 
       fill=c("blue","red"),
       border="white", bty="n",title = "d", cex = 1.5)
#NMDS time
NM.pl<-ordiplot(nmds.art2,type = "none", xlim = c(-5,5),ylim = c(-5,5), cex.axis = 1.5, cex.lab = 1.5)
p.spe<-points(NM.pl,"species",pch = 2, col= "grey20", cex= 0.7)
ordiellipse(nmds.art2, MetaRich.ART1$season,cex=1,
            draw="polygon", col= c("red","blue"),border=c("red","blue"),lwd=3,alpha = 200, kind="se", conf=0.95)

legend("bottomright", c("Summer","Winter"), 
       fill= c("red","blue"),border="white", bty="n",title = "e", cex = 1.5)
dev.off()
#################

#######################################
# Fig 5

ggplot(ST, aes(x=Organ, y=Frequency, fill=OTU)) +
  geom_bar(stat="identity", width = 0.5) + facet_wrap(~Soil) + theme_bw()+
  xlab("Organ")+ ylab("Frequency")+ theme(legend.position="top",axis.text=element_text(size=12),
                                          axis.title=element_text(size=14),legend.text=element_text(size=12),
                                          strip.text = element_text(size = 12))+
  guides(fill=guide_legend(title=NULL))


###########################################
#Fig 6
# pca plot

#coloring:
soilfactor<-factor(MetaRich.ART1$SOIL)
colvec <-  c("red","blue")
cols<-colvec[soilfactor]

rda.art1 %>% plot(type="n", corr=TRUE) %>% 
  points("sites", pch=(16:17)[soilfactor], col=colvec[soilfactor]) %>% 
  text( "biplot", arrows=TRUE, col= "black", lwd=2, cex=1.2, xpd=TRUE, len=0.05)
legend("topleft", c("Dry soil","Saline soil"), 
       col=c("red","blue"),
       pch = c(16,17), border="white", bty="n",cex = 1.2)
# varpart
plot(var.part7,Xnames = c("Host","Site","Soil*Organ"),
     bg=c("green","blue","red"),alpha=150)

##############################
#Fig 7
par(mfrow=c(1,4))

# total biomass
ggplot(GH.data, aes(x=Drought,y=Biomass, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Biomass (g)")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# Chlorophyll concentration
ggplot(GH.data, aes(x=Drought,y=Photosyntesis, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Chlorophyll concentration (%)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.text = element_text(size = 12))

###################
# Fig S1
# shoot lenght
ggplot(GH.data, aes(x=Drought,y=Lshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Shoot lenght (cm)")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# shoot wet weight
ggplot(GH.data, aes(x=Drought,y=Wshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+
  ylab("Shoot wet weight (g)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# shoot dry weight
ggplot(GH.data, aes(x=Drought,y=DWshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+
  ylab("Shoot dry weight (g)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# root lenght
ggplot(GH.data, aes(x=Drought,y=Lroot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Root lenght (cm)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
# root wet weight
ggplot(GH.data, aes(x=Drought,y=Wroot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+
  ylab("Root wet weight (g)")+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))


#barplot
host.barplot<-ggplot(host.otu.data,aes(x = Host, y = Frequency,fill = OTU)) + 
  geom_bar(position = "fill",stat = "identity", width = 0.5)+ theme_bw()+scale_y_continuous(labels = percent_format())+
  xlab("Host plant species")+ ylab("Proportional frequency")+
  labs(fill = "OTU")+ scale_fill_manual(values=otu.cols)+theme(legend.position="top")
### export all of the plots together  
margin1<-matrix(c(0,2,0,0,2,0,1,3,5,1,3,5,0,4,0,0,4,0))  
layout(margin1)
dev.off()
par(mfrow=c(3,2))
p.leaf
p.twig
p.root
p.order
host.barplot
