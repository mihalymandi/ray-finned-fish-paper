#*====
#2. Care bias in each stage====
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

####All species====
#####2a. Care_bias_score ~ Nest_building ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Nest_building<-as.factor(fish.data$Nest_building)

fish.data$Species[is.na(fish.data$Care_bias_score)] 

summary(fish.data$Care_bias_score)
summary(fish.data$Nest_building)

m2a <- phylolm(as.numeric(Care_bias_score) ~ Nest_building,
               data=fish.data,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
a.care.score.nest<-summary(m2a)


#####2b. Care_bias_score ~ Egg_care====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Egg_care<-as.factor(fish.data$Egg_care)

fish.data$Species[is.na(fish.data$Care_bias_score)] 

fish.data$Species[is.na(fish.data$Egg_care)] 
#Apistogramma_bitaeniata has no Egg care data

summary(fish.data$Care_bias_score)
summary(fish.data$Egg_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2b <- phylolm(as.numeric(Care_bias_score) ~ Egg_care,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
a.care.score.egg<-summary(m2b)

#####2c. Care_bias_score ~ Fry_care====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fry_care<-as.factor(fish.data$Fry_care)

fish.data$Species[is.na(fish.data$Fry_care)] 
#Apistogramma_bitaeniata has no Fry care data

summary(fish.data$Care_bias_score)
summary(fish.data$Fry_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2c <- phylolm(as.numeric(Care_bias_score) ~ Fry_care,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
a.care.score.fry<-summary(m2c)

####Caring only====
data<-read.csv2(file = "Fish_care_noNA.csv")
rownames(data)<-data$Species
data$Care_presence<-as.factor(data$Care_presence)
dataP<-subset(data,data$Care_presence == "Present")
dataP<-droplevels(dataP)
teleost.tree<-read.tree(file="teleost.tre")

fish.tree<-treedata(teleost.tree,dataP,sort=T,warnings=T)$phy
fish.dataP<-treedata(teleost.tree,dataP,sort=T,warnings=T)$data
fish.dataP<-as.data.frame(fish.dataP)

name.check(fish.tree,fish.dataP) #OK

length(fish.tree$tip.label) #2363 tip labels
length(fish.dataP$Species) #2363 species

#####2d. Care_bias_score ~ Nest_building ====
fish.dataP$Care_bias_score<-as.factor(fish.dataP$Care_bias_score)
fish.dataP$Nest_building<-as.factor(fish.dataP$Nest_building)

fish.dataP$Species[is.na(fish.dataP$Care_bias_score)] 

summary(fish.dataP$Care_bias_score)
summary(fish.dataP$Nest_building)

m2d <- phylolm(as.numeric(Care_bias_score) ~ Nest_building,
               data=fish.dataP,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
p.care.score.nest<-summary(m2d)


#####2e. Care_bias_score ~ Egg_care====
fish.dataP$Care_bias_score<-as.factor(fish.dataP$Care_bias_score)
fish.dataP$Egg_care<-as.factor(fish.dataP$Egg_care)

fish.dataP$Species[is.na(fish.dataP$Care_bias_score)] 

fish.dataP$Species[is.na(fish.dataP$Egg_care)] 
#Apistogramma_bitaeniata has no Egg care data

summary(fish.dataP$Care_bias_score)
summary(fish.dataP$Egg_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2e <- phylolm(as.numeric(Care_bias_score) ~ Egg_care,
               data=fish.dataP,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
p.care.score.egg<-summary(m2e)

#####2f. Care_bias_score ~ Fry_care====
fish.dataP$Care_bias_score<-as.factor(fish.dataP$Care_bias_score)
fish.dataP$Fry_care<-as.factor(fish.dataP$Fry_care)

fish.dataP$Species[is.na(fish.dataP$Fry_care)] 
#Apistogramma_bitaeniata has no Fry care data

summary(fish.dataP$Care_bias_score)
summary(fish.dataP$Fry_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2e <- phylolm(as.numeric(Care_bias_score) ~ Fry_care,
               data=fish.dataP,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
p.care.score.fry<-summary(m2e)


####Write phylolm table====
r1<-c("Care bias (all species)", "Nest building", length(a.care.score.nest$residuals),
      round(a.care.score.nest$coefficients[2],3), round(a.care.score.nest$coefficients[4],3),
      round(a.care.score.nest$coefficients[6],3), "<0.0001")
r2<-c("Care bias (caring only)", "Nest building", length(p.care.score.nest$residuals), 
      round(p.care.score.nest$coefficients[2],3), round(p.care.score.nest$coefficients[4],3),
      round(p.care.score.nest$coefficients[6],3), "<0.0001")
r3<-c("Care bias (all species)", "Egg care", length(a.care.score.egg$residuals),
      round(a.care.score.egg$coefficients[2],3), round(a.care.score.egg$coefficients[4],3),
      round(a.care.score.egg$coefficients[6],3), "<0.0001")
r4<-c("Care bias (caring only)", "Egg care", length(p.care.score.egg$residuals), 
      round(p.care.score.egg$coefficients[2],3), round(p.care.score.egg$coefficients[4],3),
      round(p.care.score.egg$coefficients[6],3), round(p.care.score.egg$coefficients[12],3))
r5<-c("Care bias (all species)", "Fry care", length(a.care.score.fry$residuals),
      round(a.care.score.fry$coefficients[2],3), round(a.care.score.fry$coefficients[4],3),
      round(a.care.score.fry$coefficients[6],3), "<0.0001")
r6<-c("Care bias (caring only)", "Fry care", length(p.care.score.fry$residuals), 
      round(p.care.score.fry$coefficients[2],3), round(p.care.score.fry$coefficients[4],3),
      round(p.care.score.fry$coefficients[6],3), "<0.0001")
t2<-rbind(r1,r2,r3,r4,r5,r6)
colnames(t2)<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
rownames(t2)<-NULL
write.csv(t2,file="Table2_phylolm.csv",row.names = F)

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

###phylolm====

####All species====
#####2a. Care_bias_score ~ Nest_building ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Nest_building<-as.factor(fish.data$Nest_building)

fish.data$Species[is.na(fish.data$Care_bias_score)] 

summary(fish.data$Care_bias_score)
summary(fish.data$Nest_building)

m2a <- phylolm(as.numeric(Care_bias_score) ~ Nest_building,
               data=fish.data,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
a.care.score.nest<-summary(m2a)


#####2b. Care_bias_score ~ Egg_care====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Egg_care<-as.factor(fish.data$Egg_care)

fish.data$Species[is.na(fish.data$Care_bias_score)] 

fish.data$Species[is.na(fish.data$Egg_care)] 
#Apistogramma_bitaeniata has no Egg care data

summary(fish.data$Care_bias_score)
summary(fish.data$Egg_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2b <- phylolm(as.numeric(Care_bias_score) ~ Egg_care,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
a.care.score.egg<-summary(m2b)

#####2c. Care_bias_score ~ Fry_care====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fry_care<-as.factor(fish.data$Fry_care)

fish.data$Species[is.na(fish.data$Fry_care)] 
#Apistogramma_bitaeniata has no Fry care data

summary(fish.data$Care_bias_score)
summary(fish.data$Fry_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2c <- phylolm(as.numeric(Care_bias_score) ~ Fry_care,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
a.care.score.fry<-summary(m2c)

####Caring only====
data<-read.csv2(file = "Fish_care_noNA.csv")
rownames(data)<-data$Species
data$Care_presence<-as.factor(data$Care_presence)
dataP<-subset(data,data$Care_presence == "Present")
dataP<-droplevels(dataP)
teleost.tree<-read.tree(file="actinopt_12k_treePL.tre")

fish.tree<-treedata(teleost.tree,dataP,sort=T,warnings=T)$phy
fish.dataP<-treedata(teleost.tree,dataP,sort=T,warnings=T)$data
fish.dataP<-as.data.frame(fish.dataP)

name.check(fish.tree,fish.dataP) #OK

length(fish.tree$tip.label) #1285 tip labels
length(fish.dataP$Species) #1285 species

#####2d. Care_bias_score ~ Nest_building ====
fish.dataP$Care_bias_score<-as.factor(fish.dataP$Care_bias_score)
fish.dataP$Nest_building<-as.factor(fish.dataP$Nest_building)

fish.dataP$Species[is.na(fish.dataP$Care_bias_score)] 

summary(fish.dataP$Care_bias_score)
summary(fish.dataP$Nest_building)

m2d <- phylolm(as.numeric(Care_bias_score) ~ Nest_building,
               data=fish.dataP,phy=fish.tree,
               boot=100,method = c("logistic_MPLE"))
p.care.score.nest<-summary(m2d)


#####2e. Care_bias_score ~ Egg_care====
fish.dataP$Care_bias_score<-as.factor(fish.dataP$Care_bias_score)
fish.dataP$Egg_care<-as.factor(fish.dataP$Egg_care)

fish.dataP$Species[is.na(fish.dataP$Care_bias_score)] 

fish.dataP$Species[is.na(fish.dataP$Egg_care)] 
#Apistogramma_bitaeniata has no Egg care data

summary(fish.dataP$Care_bias_score)
summary(fish.dataP$Egg_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2e <- phylolm(as.numeric(Care_bias_score) ~ Egg_care,
               data=fish.dataP,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
p.care.score.egg<-summary(m2e)

#####2f. Care_bias_score ~ Fry_care====
fish.dataP$Care_bias_score<-as.factor(fish.dataP$Care_bias_score)
fish.dataP$Fry_care<-as.factor(fish.dataP$Fry_care)

fish.dataP$Species[is.na(fish.dataP$Fry_care)] 
#Apistogramma_bitaeniata has no Fry care data

summary(fish.dataP$Care_bias_score)
summary(fish.dataP$Fry_care)

fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
m2e <- phylolm(as.numeric(Care_bias_score) ~ Fry_care,
               data=fish.dataP,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
p.care.score.fry<-summary(m2e)


####Write phylolm table====
r1<-c("Care bias (all species)", "Nest building", length(a.care.score.nest$residuals),
      round(a.care.score.nest$coefficients[2],3), round(a.care.score.nest$coefficients[4],3),
      round(a.care.score.nest$coefficients[6],3), "<0.0001")
r2<-c("Care bias (caring only)", "Nest building", length(p.care.score.nest$residuals), 
      round(p.care.score.nest$coefficients[2],3), round(p.care.score.nest$coefficients[4],3),
      round(p.care.score.nest$coefficients[6],3), "<0.0001")
r3<-c("Care bias (all species)", "Egg care", length(a.care.score.egg$residuals),
      round(a.care.score.egg$coefficients[2],3), round(a.care.score.egg$coefficients[4],3),
      round(a.care.score.egg$coefficients[6],3), "<0.0001")
r4<-c("Care bias (caring only)", "Egg care", length(p.care.score.egg$residuals), 
      round(p.care.score.egg$coefficients[2],3), round(p.care.score.egg$coefficients[4],3),
      round(p.care.score.egg$coefficients[6],3), round(p.care.score.egg$coefficients[12],3))
r5<-c("Care bias (all species)", "Fry care", length(a.care.score.fry$residuals),
      round(a.care.score.fry$coefficients[2],3), round(a.care.score.fry$coefficients[4],3),
      round(a.care.score.fry$coefficients[6],3), "<0.0001")
r6<-c("Care bias (caring only)", "Fry care", length(p.care.score.fry$residuals), 
      round(p.care.score.fry$coefficients[2],3), round(p.care.score.fry$coefficients[4],3),
      round(p.care.score.fry$coefficients[6],3), "<0.0001")
t2<-rbind(r1,r2,r3,r4,r5,r6)
colnames(t2)<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
rownames(t2)<-NULL
write.csv(t2,file="Table2_phylolm-fishtreeoflife.csv",row.names = F)

###pgls=====
####All species====
data<-read.csv2(file = "Fish_care_noNA.csv")
rownames(data)<-data$Species
teleost.tree<-read.tree(file="actinopt_12k_treePL.tre")
fish.tree<-treedata(teleost.tree,data,sort=T,warnings=T)$phy
fish.data<-treedata(teleost.tree,data,sort=T,warnings=T)$data
fish.data<-as.data.frame(fish.data)
name.check(fish.tree,fish.data) #OK
length(fish.tree$tip.label) #4205 tip labels
length(fish.data$Species) #4205 species

#####2a. Care_bias_score ~ Nest_building ====
comp1<- comparative.data(fish.tree, fish.data[,c("Species","Care_bias_score", "Nest_building")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls2a<-pgls(as.numeric(Care_bias_score)~Nest_building, data=comp1,lambda = "ML")
a.care.score.nest<-summary(pgls2a)
r1<-c("Care bias (all species)", "Nest building", length(a.care.score.nest$residuals),
      round(a.care.score.nest$coefficients[2],3), round(a.care.score.nest$coefficients[4],3),
      round(a.care.score.nest$coefficients[6],3), round(a.care.score.nest$coefficients[8],5))
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t1<-rbind(r2,r1)
write.csv(t1,file="2a_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
#####2b. Care_bias_score ~ Egg_care ====
fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
comp1<- comparative.data(fish.tree2, fish.data[,c("Species","Care_bias_score", "Egg_care")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls2b<-pgls(as.numeric(Care_bias_score)~Egg_care, data=comp1,lambda = "ML")
a.care.score.egg<-summary(pgls2b)
r1<-c("Care bias (all species)", "Egg care", length(a.care.score.egg$residuals),
      round(a.care.score.egg$coefficients[2],3), round(a.care.score.egg$coefficients[4],3),
      round(a.care.score.egg$coefficients[6],3), round(a.care.score.egg$coefficients[8],5))
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t1<-rbind(r2,r1)
write.csv(t1,file="2b_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
#####2c. Care_bias_score ~ Fry_care ====
fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
comp1<- comparative.data(fish.tree2, fish.data[,c("Species","Care_bias_score", "Fry_care")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls2c<-pgls(as.numeric(Care_bias_score)~Fry_care, data=comp1,lambda = "ML")
a.care.score.fry<-summary(pgls2c)
r1<-c("Care bias (all species)", "Fry care", length(a.care.score.fry$residuals),
      round(a.care.score.fry$coefficients[2],3), round(a.care.score.fry$coefficients[4],3),
      round(a.care.score.fry$coefficients[6],3), round(a.care.score.fry$coefficients[8],5))
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t1<-rbind(r2,r1)
write.csv(t1,file="2c_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
####Caring only====
data<-read.csv2(file = "Fish_care_noNA.csv")
rownames(data)<-data$Species
data$Care_presence<-as.factor(data$Care_presence)
dataP<-subset(data,data$Care_presence == "Present")
dataP<-droplevels(dataP)
teleost.tree<-read.tree(file="actinopt_12k_treePL.tre")

fish.tree<-treedata(teleost.tree,dataP,sort=T,warnings=T)$phy
fish.dataP<-treedata(teleost.tree,dataP,sort=T,warnings=T)$data
fish.dataP<-as.data.frame(fish.dataP)

name.check(fish.tree,fish.dataP) #OK

length(fish.tree$tip.label) #1285 tip labels
length(fish.dataP$Species) #1285 species

#####2d. Care_bias_score ~ Nest_building ====
comp1<- comparative.data(fish.tree, fish.dataP[,c("Species","Care_bias_score", "Nest_building")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls2d<-pgls(as.numeric(Care_bias_score)~Nest_building, data=comp1,lambda = "ML")
p.care.score.nest<-summary(pgls2d)
r1<-c("Care bias (caring only)", "Nest building", length(p.care.score.nest$residuals),
      round(p.care.score.nest$coefficients[2],3), round(p.care.score.nest$coefficients[4],3),
      round(p.care.score.nest$coefficients[6],3), round(p.care.score.nest$coefficients[8],5))
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t1<-rbind(r2,r1)
write.csv(t1,file="2d_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
#####2e. Care_bias_score ~ Egg_care ====
fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
comp1<- comparative.data(fish.tree2, fish.dataP[,c("Species","Care_bias_score", "Egg_care")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls2e<-pgls(as.numeric(Care_bias_score)~Egg_care, data=comp1,lambda = "ML")
p.care.score.egg<-summary(pgls2e)
r1<-c("Care bias (caring only)", "Egg care", length(p.care.score.egg$residuals),
      round(p.care.score.egg$coefficients[2],3), round(p.care.score.egg$coefficients[4],3),
      round(p.care.score.egg$coefficients[6],3), round(p.care.score.egg$coefficients[8],5))
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t1<-rbind(r2,r1)
write.csv(t1,file="2e_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
#####2f. Care_bias_score ~ Fry_care ====
fish.tree2<-drop.tip(fish.tree,"Apistogramma_bitaeniata")
comp1<- comparative.data(fish.tree2, fish.dataP[,c("Species","Care_bias_score", "Fry_care")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls2f<-pgls(as.numeric(Care_bias_score)~Fry_care, data=comp1,lambda = "ML")
p.care.score.fry<-summary(pgls2f)
r1<-c("Care bias (caring only)", "Fry care", length(p.care.score.fry$residuals),
      round(p.care.score.fry$coefficients[2],3), round(p.care.score.fry$coefficients[4],3),
      round(p.care.score.fry$coefficients[6],3), round(p.care.score.fry$coefficients[8],5))
r2<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t1<-rbind(r2,r1)
write.csv(t1,file="2f_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
