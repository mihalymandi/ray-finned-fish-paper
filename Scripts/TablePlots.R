#*====
#Plot for Table S1====
#*====
##Load Data=====
df<-read.csv(file = "Table S1.csv")
df[sapply(df,is.integer)]<-lapply(df[sapply(df,is.integer)],as.numeric)
df$Orders<-as.factor(df$Orders)
str(df)
df<-df[-c(64),]
df2 <- df %>%
  pivot_longer(cols = c(No.care, Female.care, Male.care, Biparental.care), names_to = "Care Type", values_to = "Count")
df2<-as.data.frame(df2)
str(df2)
df2$`Care Type`<-as.factor(df2$`Care Type`)
levels(df2$`Care Type`)<-c("Biparental care", "Female care", "Male care", "No care")
##Plot=====
ggplot(df2, aes(x = Orders, y = Count, fill = `Care Type`)) +
  geom_bar(stat = "identity") +
  labs(y = "Count")+
  theme(legend.position=c(0.06,0.9),
        panel.spacing = unit(0.1,"lines"),
        panel.background = element_rect(fill=NA),
        panel.grid=element_blank(),
        strip.background = element_rect(),
        axis.line = element_line(size = 0.8, colour = "black"),
        strip.text.x = element_text(hjust=0, size=11),
        plot.title = element_text(face="plain", size=12),
        axis.title.x = element_blank(),
        legend.key=element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("No care" = "grey", 
                               "Female care" = "red", 
                               "Male care" = "blue", 
                               "Biparental care" = "purple"))
#*====
#Plot for Table S2====
#*====
##Load Data=====
df<-read.csv(file = "Table S2.csv")
df[sapply(df,is.integer)]<-lapply(df[sapply(df,is.integer)],as.numeric)
df$Orders<-as.factor(df$Orders)
str(df)
df<-df[-c(64),]
df2 <- df %>%
  pivot_longer(cols = c(external, mouth, oviduct, pouch), names_to = "Fertilization modes", values_to = "Count")
df2<-as.data.frame(df2)
str(df2)
df2$`Fertilization modes`<-as.factor(df2$`Fertilization modes`)
levels(df2$`Fertilization modes`)<-c("External", "Mouth", "Oviduct", "Pouch")
##Plot=====
ggplot(df2, aes(x = Orders, y = Count, fill = `Fertilization modes`)) +
  geom_bar(stat = "identity") +
  labs(y = "Count")+
  theme(legend.position=c(0.06,0.9),
        panel.spacing = unit(0.1,"lines"),
        panel.background = element_rect(fill=NA),
        panel.grid=element_blank(),
        strip.background = element_rect(),
        axis.line = element_line(size = 0.8, colour = "black"),
        strip.text.x = element_text(hjust=0, size=11),
        plot.title = element_text(face="plain", size=12),
        axis.title.x = element_blank(),
        legend.key=element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("External" = "#ED7953FF", 
                               "Mouth" = "#9C179EFF", 
                               "Oviduct" = "#F0F921FF", 
                               "Pouch" = "#0D0887FF"))
