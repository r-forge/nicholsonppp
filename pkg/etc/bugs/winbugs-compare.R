load("c.Rdata")
library(ggplot2)
ppp <- sapply(res4,function(L)L$ppp)
molt <- melt(ppp)
names(molt) <- c("locus","generation","ppp")
molt$generation <- molt$generation*25
type=factor(c(rep("none",800),rep("positive",100),rep("balancing",100)))
df <- data.frame(molt,type)
library(lattice)
densityplot(~ppp|generation,df,groups=type,layout=c(1,4))

load("bugs/all.models.Rdata")
ppp <- melt(cbind(old=R$ppp$PPPVALtot[1:1000],r1=fit1$ppp,r0.7=fit0.7$ppp))
names(ppp) <- c("locus","program","ppp")
df <- data.frame(ppp,type)
densityplot(~ppp|program,df,groups=type,layout=c(1,3))
densityplot(~ppp|program,df,layout=c(1,3))

