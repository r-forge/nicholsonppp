source("../R/sim.R")
source("../R/plot.R")
source("../R/model.R")
source("../R/c.R")
dyn.load("../src/0mers_twist.so")
dyn.unload("../src/0mers_twist.so")
library(ggplot2)
library(lattice)

sim.diff <- sim.drift.selection(generations=200,s=0.05,loci=100,popsize=c(500,1000,2000))
nppp.diff <- nppp.over.time(sim.diff)
save(sim.diff,nppp.diff,file="c.models.4.diff.Rdata")
c.diff <- c.df(nppp.diff)
pdf("c.over.time.diff.pdf",h=0,w=0,paper="a4")
c.over.time(c.diff,main="Differentiation parameter c depends on population size")
dev.off()
## Compare several model fits:
c.over.time(rbind(do.fit(c.est~generation:popsize+generation,c.diff),
                  do.fit(c.est~generation:popsize,c.diff,FALSE,"only.int"),
                  do.fit(c.est~generation*popsize,c.diff,FALSE,"fitfull"),
                  do.fit(c.est~generation:popsize+generation-1,c.diff,FALSE,"no.intercept")
                  ))
## Try log-fit, similar to theoretical model:
log.data <- do.fit(log(c.est)~log(generation)+log(popsize),c.diff,FALSE,"log")
log.data$c.est <- exp(log.fit$c.est)
c.over.time(rbind(do.fit(c.est~generation:popsize+generation,
                         c.diff,label="linear"),log.data))
## The original seems to be the best:
c.over.time(do.fit(c.est~generation:popsize+generation,c.diff))



sim.same <- sim.drift.selection(generations=200,s=0.05,loci=100,popsize=1000)
nppp.same <- nppp.over.time(sim.same)
save(sim.same,nppp.same,file="c.models.4.same.Rdata")
c.same <- c.df(nppp.same)
pdf("c.over.time.same.pdf",h=0,w=0,paper="a4")
c.over.time(c.same,main="Differentiation parameter c increases linearly over time")
dev.off()
do.fit <- function(f,data,keep=TRUE,label="fit"){
  data.same <- subset(data,type=="simulated")
  fit.same <- lm(f,data.same)
  print(summary(fit.same))
  data.same$c.est <- predict(fit.same)
  data.same$population <- factor(paste(label,data.same$popsize))
  data.same$type <- factor(label)
  if(keep)data.same <- rbind(data,data.same)
  unique(data.same)
}
c.over.time(do.fit(c.est~generation,c.same))

c.all <- rbind(data.frame(c.diff,popsizes="different.population.sizes"),
               data.frame(c.same,popsizes="uniform.population.size"))
pdf("c.over.time.all.pdf",h=0,w=0,paper="a4")
c.over.time(c.all,main="Differentiation parameter c depends on population size and time",facets=~popsizes)+theme_bw()
dev.off()



load("c.models.8.Rdata")
cc <- c.df(res8)
c.over.time(cc,main="Differentiation parameter c changes with time but disagrees with theory")


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
