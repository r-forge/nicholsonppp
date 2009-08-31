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

source("../R/sim.R")
source("../R/plot.R")
library(ggplot2)
library(lattice)
## initial study of which s values approach 
s <- 10^c(-3,-2,-1.5,-1,0)
sim <- sim.drift.selection(loci=100,generations=100,s=s,p.neutral=0.1)
df <- sim2df(sim)
sim.summary.plot(df)

sims <- mlply(data.frame(s),
              function(s)
              sim.drift.selection(loci=100,generations=100,s=s))
models <- lapply(sims,function(L)nicholsonppp(L$sim[,,L$p$gen]))
save.image(file="s.models.Rdata")

compare <- function(v){
  ppp5 <- melt(sapply(models,function(L)L[[v]]))
  names(ppp5) <- c("locus","sim","value")
  densityplot(~value,ppp5,groups=sim,main=v)
}
