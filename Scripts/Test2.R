#remove missing species from data frame
(a2.lab <- a.tre2$tip.label)

(sp <- unique(m3.df$Species))
t <- table(a2.lab)
t[t>1]
t.sp <- table(sp)
t.sp[t.sp>1]
(faban <- sp[sp %in% a2.lab])  #6540  egg care_7439  7608_juv care 8694_nest
(nincsfaban <- sp[!(sp %in% a2.lab)])  #4375
(nincsfaban <- nincsfaban[order(nincsfaban)])
m4.df <- m3.df[m3.df$Species %in% faban,]
summary(m4.df)
