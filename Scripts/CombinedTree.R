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
