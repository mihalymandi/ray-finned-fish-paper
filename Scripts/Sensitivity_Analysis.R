#Sensitivity Analysis for Tables 1-3:
#Table 1. Care from associations
#Table 2. Care bias in each stage
#Table 3. Care bias and fertilization modes

#Contact for this script: oscar.garcia.miranda@outlook.com
#November 2023
#R version 4.3.2.

#Directory====
setwd("C:/Documentos/University of Debrecen/Balazs Vagi/Fish materials")

#Load Libraries====
library(ape)
library(caper)
library(beepr)
library(corHMM)
library(diversitree)
library(dplyr)
library(geiger)
library(ggplot2)
library(MuMIn)
library(phylolm)
library(phytools)
library(tidyr)
#*====
#1. Care form associations====
#*====

#Gergely Tree====
##Load Data=====
data<-read.csv2(file = "Fish_care_noNA.csv")
rownames(data)<-data$Species

##Load Tree====

teleost.tree<-read.tree(file="teleost.tre")

##Prune Tree====
fish.tree<-treedata(teleost.tree,data,sort=T,warnings=T)$phy
fish.data<-treedata(teleost.tree,data,sort=T,warnings=T)$data
fish.data<-as.data.frame(fish.data)

name.check(fish.tree,fish.data) #OK

length(fish.tree$tip.label) #7581 tip labels
length(fish.data$Species) #7581 species

fish.data$Fertilization.mode_ext<-as.factor(fish.data$Fertilization.mode_ext)
summary(fish.data$Fertilization.mode_ext)
###phylolm====

####1a. Male Egg Care ~ Male Nest building====
fish.data$Male_egg_care<-as.factor(fish.data$Male_egg_care)
fish.data$Male_nest_building<-as.factor(fish.data$Male_nest_building)

fish.data$Species[is.na(fish.data$Male_egg_care)] 
#Apistogramma_bitaeniata has no Male egg care data

summary(fish.data$Male_egg_care)
summary(fish.data$Male_nest_building)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m1a <- phylolm(as.numeric(Male_egg_care) ~ Male_nest_building,
                     data=fish.data,phy=fish.tree2,
                     boot=100,method = c("logistic_MPLE"))
m.egg.care.nest<-summary(m1a)

####1b. Male Egg Care ~ Male Fry Care====
fish.data$Male_egg_care<-as.factor(fish.data$Male_egg_care)
fish.data$Male_fry_care<-as.factor(fish.data$Male_fry_care)

fish.data$Species[is.na(fish.data$Male_egg_care)] 
#Apistogramma_bitaeniata has no Male egg care data

fish.data$Species[is.na(fish.data$Male_fry_care)] 
#Apistogramma_bitaeniata has no Male fry care data

summary(fish.data$Male_egg_care)
summary(fish.data$Male_fry_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m1b <- phylolm(as.numeric(Male_egg_care) ~ Male_fry_care,
              data=fish.data,phy=fish.tree2,
              boot=100,method = c("logistic_MPLE"))
m.egg.care.fry<-summary(m1b)

####1c. Female Egg Care ~ Female Nest Building====
fish.data$Female_egg_care<-as.factor(fish.data$Female_egg_care)
fish.data$Female_nest_building<-as.factor(fish.data$Female_nest_building)

summary(fish.data$Female_egg_care)
summary(fish.data$Female_nest_building)

m1c <- phylolm(as.numeric(Female_egg_care) ~ Female_nest_building,
               data=fish.data,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
f.egg.care.nest<-summary(m1c)

####1d. Female Egg Care ~ Female Fry Care====
fish.data$Female_egg_care<-as.factor(fish.data$Female_egg_care)
fish.data$Female_fry_care<-as.factor(fish.data$Female_fry_care)

summary(fish.data$Female_egg_care)
summary(fish.data$Female_fry_care)

m1d <- phylolm(as.numeric(Female_egg_care) ~ Female_fry_care,
               data=fish.data,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
f.egg.care.fry<-summary(m1d)

####Write phylolm table====
r1<-c("Male egg care", "Male nest building", length(m.egg.care.nest$residuals),
      round(m.egg.care.nest$coefficients[2],3), round(m.egg.care.nest$coefficients[4],3),
      round(m.egg.care.nest$coefficients[6],3), "<0.0001")
r2<-c("", "Male fry care", "", 
      round(m.egg.care.fry$coefficients[2],3), round(m.egg.care.fry$coefficients[4],3), 
      round(m.egg.care.fry$coefficients[6],3), "<0.0001")
r3<-c("Female egg care", "Female nest building", length(f.egg.care.nest$residuals), 
      round(f.egg.care.nest$coefficients[2],3), round(f.egg.care.nest$coefficients[4],3), 
      round(f.egg.care.nest$coefficients[6],3), "<0.0001")
r4<-c("", "Female fry care", "", 
      round(f.egg.care.fry$coefficients[2],3), round(f.egg.care.fry$coefficients[4],3), 
      round(f.egg.care.fry$coefficients[6],3), "<0.0001")

t1<-rbind(r1,r2,r3,r4)
colnames(t1)<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
rownames(t1)<-NULL
write.csv(t1,file="Table1_phylolm.csv",row.names = F)

##pgls

#I tried to run this but it took to long with Gergely's Tree. 

#Fishtreeoflife.org====

##Load Data=====
data<-read.csv2(file = "Fish_care_noNA.csv")
rownames(data)<-data$Species

##Load Tree====
teleost.tree<-read.tree(file="actinopt_12k_treePL.tre")

##Prune Tree====
fish.tree<-treedata(teleost.tree,data,sort=T,warnings=T)$phy
fish.data<-treedata(teleost.tree,data,sort=T,warnings=T)$data
fish.data<-as.data.frame(fish.data)

name.check(fish.tree,fish.data) #OK

length(fish.tree$tip.label) #4205 tip labels
length(fish.data$Species) #4205 species

fish.data$Fertilization.mode_ext<-as.factor(fish.data$Fertilization.mode_ext)
summary(fish.data$Fertilization.mode_ext)

###phylolm====

####1a. Male Egg Care ~ Male Nest building====
fish.data$Male_egg_care<-as.factor(fish.data$Male_egg_care)
fish.data$Male_nest_building<-as.factor(fish.data$Male_nest_building)

fish.data$Species[is.na(fish.data$Male_egg_care)] 
#Apistogramma_bitaeniata has no Male egg care data

summary(fish.data$Male_egg_care)
summary(fish.data$Male_nest_building)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m1a <- phylolm(as.numeric(Male_egg_care) ~ Male_nest_building,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
m.egg.care.nest<-summary(m1a)

####1b. Male Egg Care ~ Male Fry Care====
fish.data$Male_egg_care<-as.factor(fish.data$Male_egg_care)
fish.data$Male_fry_care<-as.factor(fish.data$Male_fry_care)

fish.data$Species[is.na(fish.data$Male_egg_care)] 
#Apistogramma_bitaeniata has no Male egg care data

fish.data$Species[is.na(fish.data$Male_fry_care)] 
#Apistogramma_bitaeniata has no Male fry care data

summary(fish.data$Male_egg_care)
summary(fish.data$Male_fry_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m1b <- phylolm(as.numeric(Male_egg_care) ~ Male_fry_care,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
m.egg.care.fry<-summary(m1b)

####1c. Female Egg Care ~ Female Nest Building====
fish.data$Female_egg_care<-as.factor(fish.data$Female_egg_care)
fish.data$Female_nest_building<-as.factor(fish.data$Female_nest_building)

summary(fish.data$Female_egg_care)
summary(fish.data$Female_nest_building)

m1c <- phylolm(as.numeric(Female_egg_care) ~ Female_nest_building,
               data=fish.data,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
f.egg.care.nest<-summary(m1c)

####1d. Female Egg Care ~ Female Fry Care====
fish.data$Female_egg_care<-as.factor(fish.data$Female_egg_care)
fish.data$Female_fry_care<-as.factor(fish.data$Female_fry_care)

summary(fish.data$Female_egg_care)
summary(fish.data$Female_fry_care)

m1d <- phylolm(as.numeric(Female_egg_care) ~ Female_fry_care,
               data=fish.data,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
f.egg.care.fry<-summary(m1d)

####Write phylolm table====

r1<-c("Male egg care", "Male nest building", length(m.egg.care.nest$residuals),
      round(m.egg.care.nest$coefficients[2],3), round(m.egg.care.nest$coefficients[4],3),
      round(m.egg.care.nest$coefficients[6],3), "<0.0001")
r2<-c("", "Male fry care", "", 
      round(m.egg.care.fry$coefficients[2],3), round(m.egg.care.fry$coefficients[4],3), 
      round(m.egg.care.fry$coefficients[6],3), "<0.0001")
r3<-c("Female egg care", "Female nest building", length(f.egg.care.nest$residuals), 
      round(f.egg.care.nest$coefficients[2],3), round(f.egg.care.nest$coefficients[4],3), 
      round(f.egg.care.nest$coefficients[6],3), "<0.0001")
r4<-c("", "Female fry care", "", 
      round(f.egg.care.fry$coefficients[2],3), round(f.egg.care.fry$coefficients[4],3), 
      round(f.egg.care.fry$coefficients[6],3), "<0.0001")

t1<-rbind(r1,r2,r3,r4)
colnames(t1)<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
rownames(t1)<-NULL
write.csv(t1,file="Table1_phylolm-fishtreeoflife.csv",row.names = F)

###pgls=====

####1a. Male Egg Care ~ Male Nest building====
fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
comp1<- comparative.data(fish.tree2, fish.data[,c("Species","Male_egg_care", "Male_nest_building")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls1a<-pgls(as.numeric(Male_egg_care)~Male_nest_building, data=comp1,lambda = "ML")
beep("mario")
m.egg.care.nest<-summary(pgls1a)
r1<-c("Male egg care", "Male nest building", length(m.egg.care.nest$residuals),
      round(m.egg.care.nest$coefficients[2],3), round(m.egg.care.nest$coefficients[4],3),
      round(m.egg.care.nest$coefficients[6],3), "<0.0001")
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t1<-rbind(r2,r1)
write.csv(t1,file="1a_PGLS-fishtreeoflife.csv",row.names = F)

####1b. Male Egg Care ~ Male Fry Care====
fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
comp2<- comparative.data(fish.tree2, fish.data[,c("Species","Male_egg_care", "Male_fry_care")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls1b<-pgls(as.numeric(Male_egg_care)~Male_fry_care, data=comp2,lambda = "ML")
beep("mario")
m.egg.care.fry<-summary(pgls1b)
r1<-c("Male egg care", "Male fry care", length(m.egg.care.fry$residuals),
      round(m.egg.care.fry$coefficients[2],3), round(m.egg.care.fry$coefficients[4],3),
      round(m.egg.care.fry$coefficients[6],3), "<0.0001")
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t2<-rbind(r2,r1)
write.csv(t2,file="1b_PGLS-fishtreeoflife.csv",row.names = F)

####1c. Female Egg Care ~ Female Nest Building====
comp3<- comparative.data(fish.tree, fish.data[,c("Species","Female_egg_care", "Female_nest_building")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls1c<-pgls(as.numeric(Female_egg_care)~Female_nest_building, data=comp3,lambda = "ML")
f.egg.care.nest<-summary(pgls1c)
r1<-c("Female egg care", "Female nest building", length(f.egg.care.nest$residuals),
      round(f.egg.care.nest$coefficients[2],3), round(f.egg.care.nest$coefficients[4],3),
      round(f.egg.care.nest$coefficients[6],3), "<0.0001")
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t3<-rbind(r2,r1)
write.csv(t3,file="1c_PGLS-fishtreeoflife.csv",row.names = F)

####1d. Female Egg Care ~ Female Fry Care====
comp4<- comparative.data(fish.tree, fish.data[,c("Species","Female_egg_care", "Female_fry_care")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls1d<-pgls(as.numeric(Female_egg_care)~Female_fry_care, data=comp4,lambda = "ML")
f.egg.care.fry<-summary(pgls1d)
r1<-c("Female egg care", "Female fry care", length(f.egg.care.fry$residuals),
      round(f.egg.care.fry$coefficients[2],3), round(f.egg.care.fry$coefficients[4],3),
      round(f.egg.care.fry$coefficients[6],3), "<0.0001")
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t4<-rbind(r2,r1)
write.csv(t4,file="1d_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
