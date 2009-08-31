source("../R/sim.R")
source("../R/plot.R")
source("../R/model.R")
dyn.load("../src/0mers_twist.so")
dyn.unload("../src/0mers_twist.so")
library(ggplot2)
sim <- sim.drift.selection(generations=200,s=0.05,loci=100,popsize=c(500,1000,2000))
df <- sim2df(sim)
res8 <- sapply((1:8)*25,function(g)nicholsonppp(sim$sim[,,g]),simplify=FALSE)
save.image("c.models.8.Rdata")


## Compare old and new fortran program estimates:
compare <- function(v,m){
  ppp5 <- melt(sapply(m,function(L)L[[v]]))
  names(ppp5) <- c("locus","sim","value")
  densityplot(~value,ppp5,groups=sim,main=v)
}
res1 <- nicholsonppp(sim$sim[,,25])
write.table(round(sim$sim[,,25]*100),
            file="/home/thocking/Desktop/testnich/Y_SIM_HMN",
            quote=F,row.names=F,col.names=F)
vars <- c("c","alpha","pi")
old.est <- lapply(vars,function(v)read.table(paste("~/Desktop/testnich/summary_",v,".out",sep=""),header=TRUE))
names(old.est) <- vars
plot(old.est$pi$MOY,res8[[1]]$p)
plot(old.est$alpha$MOY,res8[[1]]$a)
plot(old.est$c$MOY,res8[[1]]$c)

##load("c.Rdata")
load("c.models.8.Rdata")
library(lattice)
cc <- melt(sapply(res4,function(L)L$c))
names(cc) <- c("population","generation","c.est")
cc$generation <- cc$generation*25
g <- (1:4)*25
cc <- rbind(cc,data.frame(population=0,
                          generation=g,
                          c.est=g/(2*1000)))
cc$c.est.inv <- 1/cc$c.est
cc$population <- factor(cc$population)
levels(cc$population)[1] <- "g/2N"
library(latticedl,lib="~/lib")

pdf("c.est.pdf",paper="a4",h=0,w=0)
p <- dl(xyplot,cc,c.est~generation,population,type='l',
           main="c estimates do not depend on number of generations")
plog <- dl(xyplot,cc,log(c.est)~generation,population,type='l',
        sub="lines represent different populations")
plot(p,split=c(1,1,1,2))
plot(plog,split=c(1,2,1,2),newpage=FALSE)
dev.off()

qplot(generation,log(c.est),data=cc,group=population,geom="line")


load("bugs/all.models.Rdata")
cc <- data.frame(cbind(old=apply(R$c,2,mean),r1=fit1$c,r0.7=fit0.7$c))
cpts <- melt(R$c)
names(cpts) <- c("population","c")
dl(densityplot,cpts,~c,population)
splom(cc)
