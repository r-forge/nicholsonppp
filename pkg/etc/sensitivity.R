source("../R/sim.R")
source("../R/plot.R")
source("../R/model.R")
library(ggplot2)
library(lattice)
dyn.load("../src/0mers_twist.so")
## initial study of which s values approach 
s <- 10^c(-3,-2,-1.5,-1,0)
sim <- sim.drift.selection(loci=100,generations=100,s=s,p.neutral=0.1)
df <- sim2df(sim)
sim.summary.plot(df)

sims <- mlply(data.frame(s),
              function(s)
              sim.drift.selection(loci=100,generations=100,s=s))
models <- lapply(sims,function(L)nicholsonppp(L$sim[,,L$p$gen]))
save(sims,models,file="sims.models.Rdata")

load(file="s.models.Rdata")
res.df <- ldply(1:length(models),
                function(i)data.frame(ppp=models[[i]]$ppp,
                                      type=sims[[i]]$s$type,
                                      s=s[i]))
model5 <- nicholsonppp(sims[[5]]$sim[,,100])
library(latticedl)
dl(densityplot,res.df,~ppp|s,type,
   layout=c(1,length(s)),
   strip=strip.custom(strip.levels=c(TRUE,TRUE)),
   n=500)
densityplot(~ppp|s,res.df,layout=c(1,length(s)),n=500)

sim <- sims[[5]]
write.table(round(sim$sim[,,sim$p$gen]*100),
            file="/home/thocking/Desktop/testnich/Y_SIM_HMN",
            quote=F,row.names=F,col.names=F)
ppp <- read.table("~/Desktop/testnich/PPPval.out",header=TRUE)
ppp <- ppp[ppp$POP==1,4]
   
compare <- function(v){
  ppp5 <- melt(sapply(models,function(L)L[[v]]))
  names(ppp5) <- c("locus","sim","value")
  densityplot(~value,ppp5,groups=sim,main=v)
}




