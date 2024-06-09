#*====
#3. Care bias and fertilization modes====
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

###phylolm====

#####3a. Care_bias_score ~ Fertilization_code_ext ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fertilization_code_ext<-as.factor(fish.data$Fertilization_code_ext)

fish.data$Species[is.na(fish.data$Care_bias_score)] 
spp<-fish.data$Species[is.na(fish.data$Fertilization_code_ext)]


summary(fish.data$Care_bias_score)
summary(fish.data$Fertilization_code_ext)

fish.tree2<-drop.tip(fish.tree, spp)
m3a <- phylolm(as.numeric(Care_bias_score) ~ Fertilization_code_ext,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
care.score.ext<-summary(m3a)

#####3b. Care_bias_score ~ Fertilization_code_restr ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fertilization_code_restr<-as.factor(fish.data$Fertilization_code_restr)

fish.data$Species[is.na(fish.data$Care_bias_score)] 
spp<-fish.data$Species[is.na(fish.data$Fertilization_code_restr)]


summary(fish.data$Care_bias_score)
summary(fish.data$Fertilization_code_restr)

fish.tree2<-drop.tip(fish.tree, spp)
m3b <- phylolm(as.numeric(Care_bias_score) ~ Fertilization_code_restr,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
care.score.restr<-summary(m3b)

####Write phylolm table====
r1<-c("Care bias (pouch definition extended)", "Fertilization"," "," "," "," "," ")
r2<-c(" ", "external-mouth", length(care.score.ext$residuals),
      round(care.score.ext$coefficients[2],3), round(care.score.ext$coefficients[6],3),
      round(care.score.ext$coefficients[10],3), round(care.score.ext$coefficients[22],5))
r3<-c(" ", "external-oviduct", " ",
      round(care.score.ext$coefficients[3],3), round(care.score.ext$coefficients[7],3),
      round(care.score.ext$coefficients[11],3), round(care.score.ext$coefficients[23],5))
r4<-c(" ", "external-pouch", " ",
      round(care.score.ext$coefficients[4],3), round(care.score.ext$coefficients[8],3),
      round(care.score.ext$coefficients[12],3), round(care.score.ext$coefficients[24],5))

r5<-c("Care bias (pouch definition restricted)", "Fertilization"," "," "," "," "," ")
r6<-c(" ", "external-mouth", length(care.score.restr$residuals),
      round(care.score.restr$coefficients[2],3), round(care.score.restr$coefficients[6],3),
      round(care.score.restr$coefficients[10],3), round(care.score.restr$coefficients[22],5))
r7<-c(" ", "external-oviduct", " ",
      round(care.score.restr$coefficients[3],3), round(care.score.restr$coefficients[7],3),
      round(care.score.restr$coefficients[11],3), round(care.score.restr$coefficients[23],5))
r8<-c(" ", "external-pouch", " ",
      round(care.score.restr$coefficients[4],3), round(care.score.restr$coefficients[8],3),
      round(care.score.restr$coefficients[12],3), round(care.score.restr$coefficients[24],5))
t3<-rbind(r1,r2,r3,r4,r5,r6,r7,r8)
colnames(t3)<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
rownames(t3)<-NULL
write.csv(t3,file="Table3_phylolm.csv",row.names = F)

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

#####3a. Care_bias_score ~ Fertilization_code_ext ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fertilization_code_ext<-as.factor(fish.data$Fertilization_code_ext)

fish.data$Species[is.na(fish.data$Care_bias_score)] 
spp<-fish.data$Species[is.na(fish.data$Fertilization_code_ext)]


summary(fish.data$Care_bias_score)
summary(fish.data$Fertilization_code_ext)

fish.tree2<-drop.tip(fish.tree, spp)
m3a <- phylolm(as.numeric(Care_bias_score) ~ Fertilization_code_ext,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
care.score.ext<-summary(m3a)

#####3b. Care_bias_score ~ Fertilization_code_restr ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fertilization_code_restr<-as.factor(fish.data$Fertilization_code_restr)

fish.data$Species[is.na(fish.data$Care_bias_score)] 
spp<-fish.data$Species[is.na(fish.data$Fertilization_code_restr)]


summary(fish.data$Care_bias_score)
summary(fish.data$Fertilization_code_restr)

fish.tree2<-drop.tip(fish.tree, spp)
m3b <- phylolm(as.numeric(Care_bias_score) ~ Fertilization_code_restr,
               data=fish.data,phy=fish.tree2,
               boot=100,method = c("logistic_MPLE"))
care.score.restr<-summary(m3b)

####Write phylolm table====
r1<-c("Care bias (pouch definition extended)", "Fertilization"," "," "," "," "," ")
r2<-c(" ", "external-mouth", length(care.score.ext$residuals),
      round(care.score.ext$coefficients[2],3), round(care.score.ext$coefficients[6],3),
      round(care.score.ext$coefficients[10],3), round(care.score.ext$coefficients[22],5))
r3<-c(" ", "external-oviduct", " ",
      round(care.score.ext$coefficients[3],3), round(care.score.ext$coefficients[7],3),
      round(care.score.ext$coefficients[11],3), round(care.score.ext$coefficients[23],5))
r4<-c(" ", "external-pouch", " ",
      round(care.score.ext$coefficients[4],3), round(care.score.ext$coefficients[8],3),
      round(care.score.ext$coefficients[12],3), round(care.score.ext$coefficients[24],5))

r5<-c("Care bias (pouch definition restricted)", "Fertilization"," "," "," "," "," ")
r6<-c(" ", "external-mouth", length(care.score.restr$residuals),
      round(care.score.restr$coefficients[2],3), round(care.score.restr$coefficients[6],3),
      round(care.score.restr$coefficients[10],3), round(care.score.restr$coefficients[22],5))
r7<-c(" ", "external-oviduct", " ",
      round(care.score.restr$coefficients[3],3), round(care.score.restr$coefficients[7],3),
      round(care.score.restr$coefficients[11],3), round(care.score.restr$coefficients[23],5))
r8<-c(" ", "external-pouch", " ",
      round(care.score.restr$coefficients[4],3), round(care.score.restr$coefficients[8],3),
      round(care.score.restr$coefficients[12],3), round(care.score.restr$coefficients[24],5))
t3<-rbind(r1,r2,r3,r4,r5,r6,r7,r8)
colnames(t3)<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
rownames(t3)<-NULL
write.csv(t3,file="Table3_phylolm-fishtreeoflife.csv",row.names = F)

###pgls=====
data<-read.csv2(file = "Fish_care_noNA.csv")
rownames(data)<-data$Species
teleost.tree<-read.tree(file="actinopt_12k_treePL.tre")
fish.tree<-treedata(teleost.tree,data,sort=T,warnings=T)$phy
fish.data<-treedata(teleost.tree,data,sort=T,warnings=T)$data
fish.data<-as.data.frame(fish.data)
name.check(fish.tree,fish.data) #OK
length(fish.tree$tip.label) #4205 tip labels
length(fish.data$Species) #4205 species

#####3a. Care_bias_score ~ Fertilization_code_ext ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fertilization_code_ext<-as.factor(fish.data$Fertilization_code_ext)

fish.data$Species[is.na(fish.data$Care_bias_score)] 
spp<-fish.data$Species[is.na(fish.data$Fertilization_code_ext)]


summary(fish.data$Care_bias_score)
summary(fish.data$Fertilization_code_ext)


fish.tree2<-drop.tip(fish.tree, spp)

comp1<- comparative.data(fish.tree2, fish.data[,c("Species","Care_bias_score", "Fertilization_code_ext")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls3a<-pgls(as.numeric(Care_bias_score)~Fertilization_code_ext, data=comp1,lambda = "ML")
care.score.ext<-summary(pgls3a)
r1<-c("Care bias (pouch definition extended)", "Fertilization",length(care.score.ext$residuals),
      " "," "," "," ")
r2<-c(" ", "external-mouth", " ",
      round(care.score.ext$coefficients[2],3), round(care.score.ext$coefficients[6],3),
      round(care.score.ext$coefficients[10],3), round(care.score.ext$coefficients[14],5))
r3<-c(" ", "external-oviduct", " ",
      round(care.score.ext$coefficients[3],3), round(care.score.ext$coefficients[7],3),
      round(care.score.ext$coefficients[11],3), round(care.score.ext$coefficients[15],5))
r4<-c(" ", "external-pouch", " ",
      round(care.score.ext$coefficients[4],3), round(care.score.ext$coefficients[8],3),
      round(care.score.ext$coefficients[12],3), round(care.score.ext$coefficients[16],5))
r5<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t3<-rbind(r5,r1,r2,r3,r4)
write.csv(t3,file="3a_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")

#####3b. Care_bias_score ~ Fertilization_code_restr ====
fish.data$Care_bias_score<-as.factor(fish.data$Care_bias_score)
fish.data$Fertilization_code_restr<-as.factor(fish.data$Fertilization_code_restr)

fish.data$Species[is.na(fish.data$Care_bias_score)] 
spp<-fish.data$Species[is.na(fish.data$Fertilization_code_restr)]


summary(fish.data$Care_bias_score)
summary(fish.data$Fertilization_code_restr)

fish.tree2<-drop.tip(fish.tree, spp)

comp1<- comparative.data(fish.tree2, fish.data[,c("Species","Care_bias_score", "Fertilization_code_restr")], 
                         names.col=Species, vcv = T, vcv.dim = 2)
pgls3b<-pgls(as.numeric(Care_bias_score)~Fertilization_code_restr, data=comp1,lambda = "ML")
care.score.restr<-summary(pgls3b)
r1<-c("Care bias (pouch definition restricted)", "Fertilization",length(care.score.restr$residuals),
      " "," "," "," ")
r2<-c(" ", "external-mouth", " ",
      round(care.score.restr$coefficients[2],3), round(care.score.restr$coefficients[6],3),
      round(care.score.restr$coefficients[10],3), round(care.score.restr$coefficients[14],5))
r3<-c(" ", "external-oviduct", " ",
      round(care.score.restr$coefficients[3],3), round(care.score.restr$coefficients[7],3),
      round(care.score.restr$coefficients[11],3), round(care.score.restr$coefficients[16],5))
r4<-c(" ", "external-pouch", " ",
      round(care.score.restr$coefficients[4],3), round(care.score.restr$coefficients[8],3),
      round(care.score.restr$coefficients[12],3), round(care.score.restr$coefficients[16],5))
r5<-c("Response","Predictor","N species", "Estimate", "StErr", "t", "p")
t3<-rbind(r5,r1,r2,r3,r4)
write.csv(t3,file="3b_PGLS-fishtreeoflife.csv",row.names = F)
beep("mario")
