# Fertilization mode drives the evolution of parental care in ray-finned fishes
This repository contains all materials relating to the publication: Fertilization mode drives the evolution of parental care in ray-finned fishes.
The publication is based on a comprehensive dataset, that contains fertilization and parental care data for more than 7000 ray-finned fishes.

## About
This study includes several analysis on the connections between ferilization and parental care. We used several statistical and visual methods in this study.
We carried out stochastic character reconstruction for parental care types, episodes, and modes. We also calculated evolutionary transition rates between character combinations of parental care and the mode of fertilization. We also used Phylogenetic Generalised Least Squares (PGLS) model in the comparisions.

## Example Codes
To visualize the transition rates and character reconstructions we used .....



To visualize our data on phylogenetic scale, we used the `caper`, `ggtree`, `ape` and `diversitree` R packages.

```ruby
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

```


To visualize fertelization modes and parental care on the phylogeny we used `diversitree` R package.

```ruby
library(caper)
library(ape)
library(diversitree)

m.df <- read.csv2(file = "Fish_care_noNA.csv")
names(m.df)
nrow(m.df) #7604
summary(m.df)

a.tre <- read.tree(file="teleost.tre")

#fanplots
m.df$Male_nest_building  <- as.factor(m.df$Male_nest_building )
m.df$Male_egg_care     <- as.factor(m.df$Male_egg_care    )
m.df$Male_fry_care      <- as.factor(m.df$Male_fry_care     )
m.df$Female_nest_building  <- as.factor(m.df$Female_nest_building )
m.df$Female_egg_care   <- as.factor(m.df$Female_egg_care  )
m.df$Female_fry_care    <- as.factor(m.df$Female_fry_care   )
m.df$Fertilization_code_ext <- as.factor(m.df$Fertilization_code_ext)
m.df$Fertilization_code_restr <- as.factor(m.df$Fertilization_code_restr)
m.df$Fertilization.mode_ext <- as.factor(m.df$Fertilization.mode_ext)
m.df$Fertilization.mode_restr <- as.factor(m.df$Fertilization.mode_restr)
m.df$Care_bias_positive <- as.factor(m.df$Care_bias_positive)
m.df$Care_bias_positive_careonly <- as.factor(m.df$Care_bias_positive_careonly)
m.df$Care_presence <- as.factor(m.df$Care_presence)

#Care_bias
m2.df <- m.df[!is.na(m.df$Care_bias_positive) & 
                !is.na(m.df$Fertilization_code_ext),]
nrow(m2.df) 
summary(m2.df)


#Onlycaring
m2.df <- m.df[(m2.df$Care_presence %in% "Present"),]
nrow(m2.df) 

m3.df=subset(m2.df, select= c(Species, Male_nest_building,Male_egg_care,Male_fry_care,Female_nest_building,Female_egg_care,Female_fry_care))
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

m4.df$Male_nest_building  <- as.numeric(m4.df$Male_nest_building )
m4.df$Male_egg_care <- as.numeric(m4.df$Male_egg_care)
m4.df$Male_fry_care <- as.numeric(m4.df$Male_fry_care)
m4.df$Female_nest_building  <- as.numeric(m4.df$Female_nest_building )
m4.df$Female_egg_care <- as.numeric(m4.df$Female_egg_care)
m4.df$Female_fry_care <- as.numeric(m4.df$Female_fry_care)


col2.tree <- list(
  Male_nest_building=c("lightblue","blue"),
  Male_egg_care=c("lightblue","blue"),
  Male_fry_care=c("lightblue","blue"),
  Female_nest_building=c("lightpink","red1"),
  Female_egg_care=c("lightpink","red1"),
  Female_fry_care=c("lightpink","red1") )


col2.tree <- list(
  Fertilization_code=c("gray40", "palegreen", "lightpink", "paleturquoise"),
  Care_bias_positive_careonly=c("red", "indianred1", "pink", "grey90", "deepskyblue", "dodgerblue", "blue", "black"))

#(clade <- sapply(a.tre5$tip.label, function(x) m4.df[m4.df$Species==x,]$Class))

trait.plot(a.tre5, dat=m4.df, 
           #class = as.character(clade), 
           cols=col2.tree, quiet=TRUE, w=1/15, margin=1/1.4,
           cex.lab = 0.001, cex.legend = 0.001)

```

To analyse our data we used PGLS and corHMM models



## [Figures folder](Figures)
- this folder has PNG files which contain the figures that were shown in the supplementary files
- these figures were made with the R scrips shown above, and later expanded with original drawings by Mihály Mándi using [GIMP](https://www.gimp.org/)
- the figure were made with `ggplot2` and `ggthemes` R packages
- phylogenetic trees were made with `ape`,`phytools` and `diversitree`

## [Scripts folder](Scripts)
- this folder contains some of our R scripts, that we used in this study
- all scripts are in a separate file
- the scripts in the folder are examples, these codes were used in every class
- the sripts use several R packages like: `ggplot2`,`ggthemes`,`caper`,`phylolm`,`ape` and `diversitree` .....

## [Tables folder](Tables)
- this folder has several tables from the publication

## [Stats folder](Stats)
- this folder has several statistical results used in this study

## Dataset
The original dataset will be uploaded to an OpenAccess Repository, and linked in the original publication.

