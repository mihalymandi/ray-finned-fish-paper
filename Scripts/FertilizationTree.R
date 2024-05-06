library(caper)
library(ape)
library(diversitree)

m.df <- read.csv2(file = "Fish_care_noNA.csv")
names(m.df)
nrow(m.df) #7604
summary(m.df)

a.tre <- read.tree(file="teleost.tre")

#Fanplots
m.df$Male_nest_building  <- as.factor(m.df$Male_nest_building )
m.df$Male_egg_care     <- as.factor(m.df$Male_egg_care    )
m.df$Male_fry_care      <- as.factor(m.df$Male_fry_care     )
m.df$Female_nest_building  <- as.factor(m.df$Female_nest_building )
m.df$Female_egg_care   <- as.factor(m.df$Female_egg_care  )
m.df$Female_fry_care    <- as.factor(m.df$Female_fry_care   )
m.df$Fertilization_code <- as.factor(m.df$Fertilization_code)
m.df$Fertilization.mode <- as.factor(m.df$Fertilization.mode)
m.df$Care_bias_positive <- as.factor(m.df$Care_bias_positive)
m.df$Care_presence <- as.factor(m.df$Care_presence)

#Care_bias
m2.df <- m.df[!is.na(m.df$Care_bias_positive) & 
                !is.na(m.df$Fertilization_code_ext),]
nrow(m2.df) 
summary(m2.df)


#onlycaring
m2.df <- m.df[(m2.df$Care_presence %in% "Present"),]
nrow(m2.df) 

m3.df=subset(m2.df, select= c(Species, Fertilization_code_ext))
summary(m3.df)


(a.lab <- a.tre$tip.label)
(t <- table(a.lab))
t[t>1]
a.lab[order(a.lab)]
attributes(a.tre)


#Drop missing species
todrop <-a.lab[!a.lab %in% m3.df$Species]
(a.tre2 <- drop.tip(a.tre, todrop))

#Remove missing species from data frame
(a2.lab <- a.tre2$tip.label)

(sp <- unique(m3.df$Species))

t <- table(a2.lab)
t[t>1]
t.sp <- table(sp)
t.sp[t.sp>1]
(intree <- sp[sp %in% a2.lab])  
(notintree <- sp[!(sp %in% a2.lab)])  
(notintree <- notintree[order(notintree)])
m4.df <- m3.df[m3.df$Species %in% intree,]
summary(m4.df)


#Remove duplicates
(t <- table(m4.df$Species))
t[t>1]
m.df <- m.df[!m.df$Species_in_tree %in% names(t[t>1]),]
names(m.df)

#Branch length 
a.tre5 <- compute.brlen(a.tre2, method = "Grafen", power = 0.5)

row.names(m4.df) <- m4.df$Species
summary(m4.df)

m4.df$Care_bias_positive  <- as.numeric(m4.df$Care_bias_positive )
info$Fertilization_code_ext <- as.factor(info$Fertilization_code_ext)

#Calibrate data
info <- m4.dfc
tree <- a.tre5

pfert<-setNames(info$Fertilization_code_ext,
                rownames(info))
rn= rownames(info)

info <- as.data.frame(sapply(info, as.character))
rownames(info) <- rn


#Set color palette
cols<-setNames(c("black","palegreen", "lightpink","paleturquoise"),levels(pfert))

#Write csv
write.csv(m4.df,"test2.csv",row.names=FALSE)
m4.dfc<-read.csv("test2.csv",header=TRUE,row.names=1)

#Make simmap
mtrees2_ARD<-make.simmap(tree,pfert,model="ARD",nsim=100)
ARD<-summary(mtrees2_ARD)

plot(mtrees2_ARD,colors=cols,ftype="i",type="fan",fsize=0.01)
