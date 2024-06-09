
#Sensitivity Analysis for Figure 2,3b, 4, 5:
#Figure 2. Transition rates in type of care, care bias, and care duration in male
#Figure 3b. Transition rates in fertilization modes
#Figure 4. Co-evolution of the presence of parental care and internal/external fertilization
#Figure 5. Co-evolution of parental care types and fertilization modes

#Contact for this script: oscar.garcia.miranda@outlook.com
#November 2023
#R version 4.3.2.

#Directory====
setwd("C:/Documentos/University of Debrecen/Balazs Vagi/Fish materials")

#Load Libraries====
library(ape)
library(beepr)
library(caper)
library(geiger)
library(ggpubr)
library(phytools)
library(lmtest)
library(corHMM)

##Load Data=====
data <- read.csv2(file = "Fish_care_noNA.csv")
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

class(fish.data)


##*====
##a. Fertilization.mode_ext====
fish.data$Fertilization.mode_ext<-as.factor(fish.data$Fertilization.mode_ext)
summary(fish.data$Fertilization.mode_ext)
levels(fish.data$Fertilization.mode_ext)<-c("E","M","O","P")

fish.data2 <- fish.data[complete.cases(fish.data$Fertilization.mode_ext), ]
summary(fish.data2$Fertilization.mode_ext)
length(fish.data2$Species) #4086

##Prune Tree
fish.tree2<-treedata(fish.tree,fish.data2,sort=T,warnings=T)$phy
fish.data3<-treedata(fish.tree,fish.data2,sort=T,warnings=T)$data
fish.data3<-as.data.frame(fish.data3)

name.check(fish.tree2,fish.data3) #OK

length(fish.tree2$tip.label) #4086 tip labels
length(fish.data3$Species) #4086 species

fish.data3$Fertilization.mode_ext<-as.factor(fish.data3$Fertilization.mode_ext)
levels(fish.data3$Fertilization.mode_ext)

mode<-setNames(fish.data3$Fertilization.mode_ext, rownames(fish.data3))
mode

#equal rates model
fitER<-fitDiscrete(fish.tree2,mode,model="ER")
beep("mario")
print(fitER,digits=3)

#symmetric model
fitSYM<-fitDiscrete(fish.tree2,mode,model="SYM")
beep("mario")
print(fitSYM,digits=3)

#all rates different model
fitARD<-fitDiscrete(fish.tree2,mode,model="ARD")
beep("mario")
print(fitARD,digits=3)

###Comparing models====
a1<-lrtest(fitER,fitSYM) #p=0.006 - means that SYM is better
a2<-lrtest(fitER,fitARD) #p=<2.2e-16 - means that ARD is better
a3<-lrtest(fitSYM,fitARD) #p=<2.2e-16 - means that ARD is better
A<-data.frame(matrix(c(round(a1$LogLik,3),round(a2$LogLik,3),round(a3$LogLik,3),
                       a1$Df,a2$Df,a3$Df,
                       round(a1$Chisq,3),round(a2$Chisq,3),round(a3$Chisq,3),
                       NA,"<0.0001",NA,"<0.0001",NA,"<0.0001"),ncol=4))
colnames(A)<-c("LogLik","Df","Chisq","Pr(>Chisq)")
A$Model<-c("ER","SYM","ER","ARD","SYM","ARD")
A$Variable<-c("Fertilization mode",NA,NA,NA,NA,NA)
A<-A[,c(6,5,1:4)]
A

aic<-setNames(c(AIC(fitER),AIC(fitSYM),AIC(fitARD)),
              c("ER","SYM","ARD"))
aic
aic.w(aic)
###Akaike====
t<-round(data.frame(
  k=c(fitER$opt$k,fitSYM$opt$k,fitARD$opt$k),
  logL=c(round(logLik(fitER),3),round(logLik(fitSYM),3),round(logLik(fitARD),3)),
  AIC=round(aic,3),Akaike.w=as.vector(aic.w(aic))),3)
t$Model<-rownames(t)
t$Trait<-c("Fertilization mode",NA,NA)
t<-t[,c(6,5,1:3)]
At<-t
write.csv(t,file="a. Transition_Fertilization_Mode.csv",row.names = T)

###Plots====
dev.off()
plot(fitER,color=FALSE,width=TRUE,min.lwd=1,max.lwd=1,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(a) Equal rate model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitER),2)),cex.main=0.8,line=-1)

plot(fitSYM,color=F,width=T,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(b) Symmetric model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitSYM),2)),cex.main=0.8,line=-1)

plot(fitARD,color=TRUE,width=TRUE,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=F,xlim=c(-3,1),ylim=c(-1,0))
title(main="(c) All rates different model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitARD),2)),cex.main=0.8,line=-1)

##*====
##b. Parental.care_1====

fish.data$Parental.care_1<-as.factor(fish.data$Parental.care_1)
summary(fish.data$Parental.care_1)
levels(fish.data$Parental.care_1)<-c("B","F","M","N")

##Prune Tree
fish.treeb<-treedata(fish.tree,fish.data,sort=T,warnings=T)$phy
fish.datab<-treedata(fish.tree,fish.data,sort=T,warnings=T)$data
fish.datab<-as.data.frame(fish.datab)

name.check(fish.treeb,fish.datab) #OK

length(fish.treeb$tip.label) #4205 tip labels
length(fish.datab$Species) #4205 species

fish.datab$Parental.care_1<-as.factor(fish.datab$Parental.care_1)
levels(fish.datab$Parental.care_1)

care<-setNames(fish.datab$Parental.care_1, rownames(fish.datab))
care

#equal rates model
fitERb<-fitDiscrete(fish.treeb,care,model="ER")
beep("mario")
print(fitERb,digits=3)

#symmetric model
fitSYMb<-fitDiscrete(fish.treeb,care,model="SYM")
beep("mario")
print(fitSYMb,digits=3)

#all rates different model
fitARDb<-fitDiscrete(fish.treeb,care,model="ARD")
beep("mario")
print(fitARDb,digits=3)

###Comparing models====
b1<-lrtest(fitERb,fitSYMb) #p = < 2.2e-16 SYM is better
b2<-lrtest(fitERb,fitARDb) #p = < 2.2e-16 ARD is better
b3<-lrtest(fitSYMb,fitARDb) #p = 3.125e-09 ARD is better
B<-data.frame(matrix(c(round(b1$LogLik,3),round(b2$LogLik,3),round(b3$LogLik,3),
                       b1$Df,b2$Df,b3$Df,
                       round(b1$Chisq,3),round(b2$Chisq,3),round(b3$Chisq,3),
                       NA,"<0.0001",NA,"<0.0001",NA,"<0.0001"),ncol=4))
colnames(B)<-c("LogLik","Df","Chisq","Pr(>Chisq)")
B$Model<-c("ER","SYM","ER","ARD","SYM","ARD")
B$Variable<-c("Type of care",NA,NA,NA,NA,NA)
B<-B[,c(6,5,1:4)]
B


aicb<-setNames(c(AIC(fitERb),AIC(fitSYMb),AIC(fitARDb)),
              c("ER","SYM","ARD"))
aicb
aic.w(aicb)
###Alaike====
t<-round(data.frame(
  k=c(fitERb$opt$k,fitSYMb$opt$k,fitARDb$opt$k),
  logL=c(round(logLik(fitERb),3),round(logLik(fitSYMb),3),round(logLik(fitARDb),3)),
  AIC=round(aicb,3),Akaike.w=as.vector(aic.w(aicb))),3)
t$Model<-rownames(t)
t$Trait<-c("Type of care",NA,NA)
t<-t[,c(6,5,1:3)]
Bt<-t
write.csv(t,file="b. Transition_Parental_care.csv",row.names = T)

###Plots====
dev.off()
plot(fitERb,color=FALSE,width=TRUE,min.lwd=1,max.lwd=1,
     signif=5,show.zeros=FALSE,text=TRUE,cex.size=)
title(main="(a) Equal rate model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitERb),2)),cex.main=0.8,line=-1)

plot(fitSYMb,color=F,width=T,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(b) Symmetric model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitSYMb),2)),cex.main=0.8,line=-1)

plot(fitARDb,color=TRUE,width=TRUE,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=F,xlim=c(-3,1),ylim=c(-1,0))
title(main="(c) All rates different model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitARDb),2)),cex.main=0.8,line=-1)

##*====
##c. Care_bias_transitions====

data <- read.csv2(file = "Fish_care_noNA_updated.csv")
rownames(data)<-data$Species

teleost.tree<-read.tree(file="actinopt_12k_treePL.tre")

##Prune Tree====
fish.tree<-treedata(teleost.tree,data,sort=T,warnings=T)$phy
fish.data<-treedata(teleost.tree,data,sort=T,warnings=T)$data
fish.data<-as.data.frame(fish.data)

name.check(fish.tree,fish.data) #OK

length(fish.tree$tip.label) #4205 tip labels
length(fish.data$Species) #4205 species

fish.data$Care_bias_transitions<-as.factor(fish.data$Care_bias_transitions)
summary(fish.data$Care_bias_transitions)
levels(fish.data$Care_bias_transitions)<-c("E","FB","MB","NC")

##Prune Tree
fish.treec<-treedata(fish.tree,fish.data,sort=T,warnings=T)$phy
fish.datac<-treedata(fish.tree,fish.data,sort=T,warnings=T)$data
fish.datac<-as.data.frame(fish.datac)

name.check(fish.treec,fish.datac) #OK

length(fish.treec$tip.label) #4205 tip labels
length(fish.datac$Species) #4205 species

fish.datac$Care_bias_transitions<-as.factor(fish.datac$Care_bias_transitions)
levels(fish.datac$Care_bias_transitions)

careb<-setNames(fish.datac$Care_bias_transitions, rownames(fish.datac))
careb

#equal rates model
fitERc<-fitDiscrete(fish.treec,careb,model="ER")
beep("mario")
print(fitERc,digits=3)

#symmetric model
fitSYMc<-fitDiscrete(fish.treec,careb,model="SYM")
beep("mario")
print(fitSYMc,digits=3)

#all rates different model
fitARDc<-fitDiscrete(fish.treec,careb,model="ARD")
beep("mario")
print(fitARDc,digits=3)

###Comparing models====
c1<-lrtest(fitERc,fitSYMc) #p = < 2.2e-16 SYM is better
c2<-lrtest(fitERc,fitARDc) #p = < 2.2e-16 ARD is better
c3<-lrtest(fitSYMc,fitARDc) #p = 3.216e-11 ARD is better
C<-data.frame(matrix(c(round(c1$LogLik,3),round(c2$LogLik,3),round(c3$LogLik,3),
                       c1$Df,c2$Df,c3$Df,
                       round(c1$Chisq,3),round(c2$Chisq,3),round(c3$Chisq,3),
                       NA,"<0.0001",NA,"<0.0001",NA,"<0.0001"),ncol=4))
colnames(C)<-c("LogLik","Df","Chisq","Pr(>Chisq)")
C$Model<-c("ER","SYM","ER","ARD","SYM","ARD")
C$Variable<-c("Care bias",NA,NA,NA,NA,NA)
C<-C[,c(6,5,1:4)]
C

aicc<-setNames(c(AIC(fitERc),AIC(fitSYMc),AIC(fitARDc)),
               c("ER","SYM","ARD"))
aicc
aic.w(aicc)
###Akaike====
t<-round(data.frame(
  k=c(fitERc$opt$k,fitSYMc$opt$k,fitARDc$opt$k),
  logL=c(round(logLik(fitERc),3),round(logLik(fitSYMc),3),round(logLik(fitARDc),3)),
  AIC=round(aicc,3),Akaike.w=as.vector(aic.w(aicc))),3)
t$Model<-rownames(t)
t$Trait<-c("Care bias",NA,NA)
t<-t[,c(6,5,1:3)]
Ct<-t
write.csv(t,file="c. Transition_Care_bias.csv",row.names = T)

###Plots====
dev.off()
plot(fitERc,color=FALSE,width=TRUE,min.lwd=1,max.lwd=1,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(a) Equal rate model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitERc),2)),cex.main=0.8,line=-1)

plot(fitSYMc,color=F,width=T,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(b) Symmetric model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitSYMc),2)),cex.main=0.8,line=-1)

plot(fitARDc,color=TRUE,width=TRUE,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=F,xlim=c(-3,1),ylim=c(-1,0))
title(main="(c) All rates different model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitARDc),2)),cex.main=0.8,line=-1)

##*====
##d. Male_care_duration====
fish.data$Male_care_duration<-as.factor(fish.data$Male_care_duration)
summary(fish.data$Male_care_duration)
levels(fish.data$Male_care_duration)<-c("EC","FC","NB","NC")

##Prune Tree
fish.treed<-treedata(fish.tree,fish.data,sort=T,warnings=T)$phy
fish.datad<-treedata(fish.tree,fish.data,sort=T,warnings=T)$data
fish.datad<-as.data.frame(fish.datad)

name.check(fish.treed,fish.datad) #OK

length(fish.treed$tip.label) #4205 tip labels
length(fish.datad$Species) #4205 species

fish.datad$Male_care_duration<-as.factor(fish.datad$Male_care_duration)
levels(fish.datad$Male_care_duration)

male<-setNames(fish.datad$Male_care_duration, rownames(fish.datad))
male

#equal rates model
fitERd<-fitDiscrete(fish.treed,male,model="ER")
beep("mario")
print(fitERd,digits=3)

#symmetric model
fitSYMd<-fitDiscrete(fish.treed,male,model="SYM")
beep("mario")
print(fitSYMd,digits=3)

#all rates different model
fitARDd<-fitDiscrete(fish.treed,male,model="ARD")
beep("mario")
print(fitARDd,digits=3)


###Comparing models====
d1<-lrtest(fitERd,fitSYMd) #p = < 2.2e-16 SYM is better
d2<-lrtest(fitERd,fitARDd) #p = < 2.2e-16 ARD is better
d3<-lrtest(fitSYMd,fitARDd) #p = 3.216e-11 ARD is better
D<-data.frame(matrix(c(round(d1$LogLik,3),round(d2$LogLik,3),round(d3$LogLik,3),
                       d1$Df,d2$Df,d3$Df,
                       round(d1$Chisq,3),round(d2$Chisq,3),round(d3$Chisq,3),
                       NA,"<0.0001",NA,"<0.0001",NA,"<0.0001"),ncol=4))
colnames(D)<-c("LogLik","Df","Chisq","Pr(>Chisq)")
D$Model<-c("ER","SYM","ER","ARD","SYM","ARD")
D$Variable<-c("Care duration in males",NA,NA,NA,NA,NA)
D<-D[,c(6,5,1:4)]
D

aicd<-setNames(c(AIC(fitERd),AIC(fitSYMd),AIC(fitARDd)),
               c("ER","SYM","ARD"))
aicd
aic.w(aicd)
###Akaike====
t<-round(data.frame(
  k=c(fitERd$opt$k,fitSYMd$opt$k,fitARDd$opt$k),
  logL=c(round(logLik(fitERd),3),round(logLik(fitSYMd),3),round(logLik(fitARDd),3)),
  AIC=round(aicd,3),Akaike.w=as.vector(aic.w(aicd))),3)
t$Model<-rownames(t)
t$Trait<-c("Care duration in males",NA,NA)
t<-t[,c(6,5,1:3)]
Dt<-t
write.csv(t,file="d. Transition_Male_care_duration.csv",row.names = T)

###Plots====
dev.off()
plot(fitERd,color=FALSE,width=TRUE,min.lwd=1,max.lwd=1,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(a) Equal rate model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitERd),2)),cex.main=0.8,line=-1)

plot(fitSYMd,color=F,width=T,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(b) Symmetric model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitSYMd),2)),cex.main=0.8,line=-1)

plot(fitARDd,color=TRUE,width=TRUE,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=F,xlim=c(-3,1),ylim=c(-1,0))
title(main="(c) All rates different model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitARDd),2)),cex.main=0.8,line=-1)
##*====
##e. Female_care_duration====
fish.data$Female_care_duration<-as.factor(fish.data$Female_care_duration)
summary(fish.data$Female_care_duration)
levels(fish.data$Female_care_duration)<-c("EC","FC","NB","NC")

##Prune Tree
fish.treee<-treedata(fish.tree,fish.data,sort=T,warnings=T)$phy
fish.datae<-treedata(fish.tree,fish.data,sort=T,warnings=T)$data
fish.datae<-as.data.frame(fish.datae)

name.check(fish.treee,fish.datae) #OK

length(fish.treee$tip.label) #4205 tip labels
length(fish.datae$Species) #4205 species

fish.datae$Female_care_duration<-as.factor(fish.datae$Female_care_duration)
levels(fish.datae$Female_care_duration)

female<-setNames(fish.datae$Female_care_duration, rownames(fish.datae))
female

#equal rates model
fitERe<-fitDiscrete(fish.treee,female,model="ER")
beep("mario")
print(fitERe,digits=3)

#symmetric model
fitSYMe<-fitDiscrete(fish.treee,female,model="SYM")
beep("mario")
print(fitSYMe,digits=3)

#all rates different model
fitARDe<-fitDiscrete(fish.treee,female,model="ARD")
beep("mario")
print(fitARDe,digits=3)

###Comparing models====
e1<-lrtest(fitERe,fitSYMe) #p = < 2.2e-16 SYM is better
e2<-lrtest(fitERe,fitARDe) #p = < 2.2e-16 ARD is better
e3<-lrtest(fitSYMe,fitARDe) #p = 3.216e-11 ARD is better
E<-data.frame(matrix(c(round(e1$LogLik,3),round(e2$LogLik,3),round(e3$LogLik,3),
                       e1$Df,e2$Df,e3$Df,
                       round(e1$Chisq,3),round(e2$Chisq,3),round(e3$Chisq,3),
                       NA,"<0.0001",NA,"<0.0001",NA,"<0.0001"),ncol=4))
colnames(E)<-c("LogLik","Df","Chisq","Pr(>Chisq)")
E$Model<-c("ER","SYM","ER","ARD","SYM","ARD")
E$Variable<-c("Care duration in females",NA,NA,NA,NA,NA)
E<-E[,c(6,5,1:4)]
E

aice<-setNames(c(AIC(fitERe),AIC(fitSYMe),AIC(fitARDe)),
               c("ER","SYM","ARD"))
aice
aic.w(aice)
###Akaike====
t<-round(data.frame(
  k=c(fitERe$opt$k,fitSYMe$opt$k,fitARDe$opt$k),
  logL=c(round(logLik(fitERe),3),round(logLik(fitSYMe),3),round(logLik(fitARDe),3)),
  AIC=round(aice,3),Akaike.w=as.vector(aic.w(aice))),3)
t$Model<-rownames(t)
t$Trait<-c("Care duration in females",NA,NA)
t<-t[,c(6,5,1:3)]
Et<-t
write.csv(t,file="e. Transition_Female_care_duration.csv",row.names = T)

###Plots====
dev.off()
plot(fitERe,color=FALSE,width=TRUE,min.lwd=1,max.lwd=1,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(a) Equal rate model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitERe),2)),cex.main=0.8,line=-1)

plot(fitSYMe,color=F,width=T,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=TRUE)
title(main="(b) Symmetric model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitSYMe),2)),cex.main=0.8,line=-1)

plot(fitARDe,color=TRUE,width=TRUE,min.lwd=1,max.lwd=3,
     signif=5,show.zeros=FALSE,text=F,xlim=c(-3,1),ylim=c(-1,0))
title(main="(c) All rates different model",cex.main=1.2,line=0)
title(main=paste("AIC = ",round(AIC(fitARDe),2)),cex.main=0.8,line=-1)


##abcde.TABLES ====
ABCDE<-rbind(A,B,C,D,E)
write.csv(ABCDE,file="Table S4 Comparing Transition Models.csv",row.names = F)

ABCDEt<-rbind(At,Bt,Ct,Dt,Et)
write.csv(ABCDEt,file="Table S5 Comparing Transition Models-AIC.csv",row.names = F)

##f. Fert.mode_ext-Care.bias.binary====
fish.data$Fertilization.mode_ext  <- as.factor(fish.data$Fertilization.mode_ext )
summary(fish.data$Fertilization.mode_ext)

fish.data$Care_bias_binary <- as.factor(fish.data$Care_bias_binary)
summary(fish.data$Care_bias_binary)

fish.data2 <- fish.data[complete.cases(fish.data$Fertilization.mode_ext), ]
summary(fish.data2$Fertilization.mode_ext)
length(fish.data2$Species) #4086

##Prune Tree
fish.tree2<-treedata(fish.tree,fish.data2,sort=T,warnings=T)$phy
fish.data3<-treedata(fish.tree,fish.data2,sort=T,warnings=T)$data
fish.data3<-as.data.frame(fish.data3)

name.check(fish.tree2,fish.data3) #OK

length(fish.tree2$tip.label) #4086 tip labels
length(fish.data3$Species) #4086 species

fish.data3$Fertilization.mode_ext<-as.factor(fish.data3$Fertilization.mode_ext)
summary(fish.data3$Fertilization.mode_ext)
fish.data3$Care_bias_binary <- as.factor(fish.data3$Care_bias_binary)
summary(fish.data3$Care_bias_binary)

fish.data5<-subset(fish.data3, select= c(Species, Fertilization.mode_ext, Care_bias_binary))
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "oviduct"]) #304
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "pouch"]) #167
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "mouth"]) #46
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "external"]) #3569

length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "oviduct" & fish.data5$Care_bias_binary == "-1"]) #291
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "oviduct" & fish.data5$Care_bias_binary == " 0"]) #9
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "oviduct" & fish.data5$Care_bias_binary == " 1"]) #4

length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "pouch" & fish.data5$Care_bias_binary == "-1"]) #2
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "pouch" & fish.data5$Care_bias_binary == " 0"]) #44
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "pouch" & fish.data5$Care_bias_binary == " 1"]) #121

length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "mouth" & fish.data5$Care_bias_binary == "-1"]) #43
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "mouth" & fish.data5$Care_bias_binary == " 0"]) #2
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "mouth" & fish.data5$Care_bias_binary == " 1"]) #1

length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "external" & fish.data5$Care_bias_binary == "-1"]) #140
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "external" & fish.data5$Care_bias_binary == " 0"]) #2921
length(fish.data5$Species[fish.data5$Fertilization.mode_ext == "external" & fish.data5$Care_bias_binary == " 1"]) #508


phy <- fish.tree2
phy <- multi2di(phy)
data <- fish.data5
plot(phy, show.tip.label = FALSE)
colnames(data)
data.sort <- data.frame(data[, 2], data[, 3], row.names = data[, 1])
data.sort <- data.sort[phy$tip.label, ]
tiplabels(pch = 16, col = data.sort[, 1] , cex = 0.5)
tiplabels(pch = 16, col = data.sort[, 2] , cex = 0.5, offset = 1.2)

MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
beep("mario")
load("corHMMResults.Rsave")
MK_3state
plotMKmodel(MK_3state)
plot(as.Qmatrix(MK_3state),show.zeros=FALSE,lwd=1,
     cex.traits=0.1)

fig5_table <- data.frame(
  Legend = c("external_F", "external_N", "external_M", "mouth_F", "mouth_N", "mouth_M", "oviduct_F", "oviduct_N", "oviduct_M", "pouch_F", "pouch_N", "pouch_M"),
  external_F = c(NA, 0.0001, 0.0002, 0.0521, NA, NA, 0.0000, NA, NA, 0.0000, NA, NA),
  external_N = c(0.0297, NA, 0.0013, NA, 0.5878, NA, NA, 0.0019, NA, NA, 0.0000, NA),
  external_M = c(0.0000, 0.0004, NA, NA, NA, 1.0312, NA, NA, 0.0000, NA, NA, 0.0003),
  mouth_F = c(0.0070, NA, NA, NA, 4.6891, 14.0, 0.0000, NA, NA, 0.0000, 0.0000, NA),
  mouth_N = c(NA, 0.0000, NA, 0.0416, NA, 10.8574, NA, 0.0000, 0.0000, NA, NA, NA),
  mouth_M = c(NA, NA, 0.0000, 0.5014, 0.0000, NA, NA, 0.0033, 0.0574, NA, 0.0067, 0.0000),
  oviduct_F = c(0.0009, NA, NA, 0.0000, NA, NA, NA, NA, 0.0254, NA, NA, NA),
  oviduct_N = c(NA, 0.0001, NA, NA, 0.0100, NA, 0.0005, NA, NA, 0.0000, NA, NA),
  oviduct_M = c(NA, NA, 0.0000, NA, NA, NA, 0.0000, NA, 0.0000, NA, NA, 0.0000),
  pouch_F = c(0.0000, NA, 0.0000, 0.0000, NA, 0.0473, 0.0000, 0.0000, NA, 0.0229, 0.0067, 0.0000),
  pouch_N = c(NA, 0.0000, NA, NA, NA, 0.0000, 0.0000, NA, NA, NA, NA, 0.0005),
  pouch_M = c(NA, NA, 0.0000, NA, 0.0478, NA, NA, NA, 0.0000, NA, NA, 0.0011)
)
write.csv(fig5_table,file="Fig5_transitionrates.csv",row.names = T)

##g. Fert.binary-Care.presence====
fish.data$Fertilization.mode_ext  <- as.factor(fish.data$Fertilization.mode_ext )
summary(fish.data$Fertilization.mode_ext )
fish.data$Fertilization.mode_binary<-NA
fish.data$Fertilization.mode_binary[fish.data$Fertilization.mode_ext == "external"] <- "External"
fish.data$Fertilization.mode_binary[fish.data$Fertilization.mode_ext == "mouth"] <- "Internal"
fish.data$Fertilization.mode_binary[fish.data$Fertilization.mode_ext == "oviduct"] <- "Internal"
fish.data$Fertilization.mode_binary[fish.data$Fertilization.mode_ext == "pouch"] <- "Internal"
fish.data$Fertilization.mode_binary<-as.factor(fish.data$Fertilization.mode_binary)
summary(fish.data$Fertilization.mode_binary)

fish.data$Care_presence <- as.factor(fish.data$Care_presence)
summary(fish.data$Care_presence)

fish.data6 <- fish.data[complete.cases(fish.data$Fertilization.mode_binary), ]
summary(fish.data6$Fertilization.mode_binary)
length(fish.data6$Fertilization.mode_binary) #4086

##Prune Tree
fish.tree4<-treedata(fish.tree,fish.data6,sort=T,warnings=T)$phy
fish.data7<-treedata(fish.tree,fish.data6,sort=T,warnings=T)$data
fish.data7<-as.data.frame(fish.data7)

name.check(fish.tree4,fish.data7) #OK

length(fish.tree4$tip.label) #4086 tip labels
length(fish.data7$Species) #4086 species

fish.data7$Fertilization.mode_binary<-as.factor(fish.data7$Fertilization.mode_binary)
summary(fish.data7$Fertilization.mode_binary)
fish.data7$Care_presence <- as.factor(fish.data7$Care_presence)
summary(fish.data7$Care_presence)

fish.data8<-subset(fish.data7, select= c(Species, Fertilization.mode_binary, Care_presence))

phy <- fish.tree4
phy <- multi2di(phy)
data <- fish.data8
plot(phy, show.tip.label = FALSE)
colnames(data)
data.sort <- data.frame(data[, 2], data[, 3], row.names = data[, 1])
data.sort <- data.sort[phy$tip.label, ]
tiplabels(pch = 16, col = data.sort[, 1] , cex = 0.5)
tiplabels(pch = 16, col = data.sort[, 2] , cex = 0.5, offset = 1.2)

MK_4state <- corHMM(phy = phy, data = data, rate.cat = 1)
beep("mario")
load("corHMMResults.Rsave")
MK_4state
plotMKmodel(MK_4state)


##h. Fert.binary.restr-Care.presence====
fish.data$Fertilization.mode_ext  <- as.factor(fish.data$Fertilization.mode_ext )
summary(fish.data$Fertilization.mode_ext)
fish.data$Fertilization.mode_restr  <- as.factor(fish.data$Fertilization.mode_restr )
summary(fish.data$Fertilization.mode_restr)
fish.data$Fertilization.mode_binary_restr<-NA
fish.data$Fertilization.mode_binary_restr[fish.data$Fertilization.mode_restr == "external"] <- "External"
fish.data$Fertilization.mode_binary_restr[fish.data$Fertilization.mode_restr == "mouth"] <- "Internal"
fish.data$Fertilization.mode_binary_restr[fish.data$Fertilization.mode_restr == "oviduct"] <- "Internal"
fish.data$Fertilization.mode_binary_restr[fish.data$Fertilization.mode_restr == "pouch"] <- "Internal"
fish.data$Fertilization.mode_binary_restr<-as.factor(fish.data$Fertilization.mode_binary_restr)
summary(fish.data$Fertilization.mode_binary_restr)

fish.data$Care_presence <- as.factor(fish.data$Care_presence)
summary(fish.data$Care_presence)

fish.data9 <- fish.data[complete.cases(fish.data$Fertilization.mode_binary_restr), ]
summary(fish.data9$Fertilization.mode_binary_restr)
length(fish.data9$Species) #4086

##Prune Tree
fish.tree5<-treedata(fish.tree,fish.data9,sort=T,warnings=T)$phy
fish.data10<-treedata(fish.tree,fish.data9,sort=T,warnings=T)$data
fish.data10<-as.data.frame(fish.data10)

name.check(fish.tree5,fish.data10) #OK

length(fish.tree5$tip.label) #4086 tip labels
length(fish.data10$Species) #4086 species

fish.data10$Fertilization.mode_binary_restr<-as.factor(fish.data10$Fertilization.mode_binary_restr)
summary(fish.data10$Fertilization.mode_binary_restr)
fish.data10$Care_presence <- as.factor(fish.data10$Care_presence)
summary(fish.data10$Care_presence)

fish.data11<-subset(fish.data10, select= c(Species, Fertilization.mode_binary_restr, Care_presence))

phy <- fish.tree5
phy <- multi2di(phy)
data <- fish.data11
plot(phy, show.tip.label = FALSE)
colnames(data)
data.sort <- data.frame(data[, 2], data[, 3], row.names = data[, 1])
data.sort <- data.sort[phy$tip.label, ]
tiplabels(pch = 16, col = data.sort[, 1] , cex = 0.5)
tiplabels(pch = 16, col = data.sort[, 2] , cex = 0.5, offset = 1.2)

MK_5state <- corHMM(phy = phy, data = data, rate.cat = 1)
beep("mario")
load("corHMMResults.Rsave")
MK_5state
plotMKmodel(MK_5state)
