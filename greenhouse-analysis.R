
##############################################
##############################################
# Greenhouse data analysis
##############################################
##############################################

library(vegan)
library(mvabund)
library(rjags)
library(MASS)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(mvabund)
library(reshape)
library(scales)

# data input
Hordeum2<-read.csv("Hordeum2.csv", header = T, row.names = 1)
colnames(Hordeum2)
#delete extra columns
Hordeum2$X.4<-NULL
#delet rows
Hordeum.new <- Hordeum2[-c(145),] 
#output
write.csv(Hordeum.new, "new.h.csv")
#input
Hordeum2<-read.csv("new.h.csv",header = T, row.names = 1 )

Hordeum2$Biomass<-Hordeum2$Wdshoot+Hordeum2$Wdroot

# 3-way anova

# ?aov
# shoot lenght
Hordeum2Lshoot<-aov(Lshoot~Fungi*Drought*salinity,data = Hordeum2)
Hordeum2Lshoot.sum<-summary(Hordeum2Lshoot)
# shoot wet weight
GH.Wshoot<-aov(Wshoot~Fungi*Drought*salinity,data = Hordeum2)
ghWshoot.sum<-summary(GH.Wshoot)
# shoot dry weight
GH.DWshoot<-aov(DWshoot~Fungi*Drought*Salinity,data = Hordeum2)
ghDWshoot.sum<-summary(GH.DWshoot)
# root lenght
GH.Lroot<-aov(Lroot~Fungi*Drought*Salinity,data = Hordeum2)
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
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()
# shoot wet weight
ggplot(GH.data, aes(x=Drought,y=Wshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()
# shoot dry weight
ggplot(GH.data, aes(x=Drought,y=DWshoot, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()
# Chlorophyll concentration
ggplot(GH.data, aes(x=Drought,y=Photosyntesis, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()
# total biomass
ggplot(GH.data, aes(x=Drought,y=Biomass, fill=Salinity)) + 
  facet_wrap(~Fungi)+geom_boxplot(position = "dodge")+theme_bw()

