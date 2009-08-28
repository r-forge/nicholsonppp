source("../R/sim.R")
source("../R/plot.R")
source("../R/model.R")
dyn.load("../src/0mers_twist.so")
sim <- sim.drift.selection(generations=200,s=0.05,loci=100)
df <- sim2df(sim)
res20 <- sapply((1:20)*10,function(g)nicholsonppp(sim$sim[,,g]),simplify=FALSE)
load("c.Rdata")
save.image("allc.Rdata")
library(ggplot2)
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
plog <- dl(xyplot,cc,log(c.est)~generation,population,type='l',
           main="c estimates do not depend on number of generations")
p <- dl(xyplot,cc,c.est~generation,population,type='l',
        sub="lines represent different populations")
plot(p,split=c(1,1,1,2))
plot(plog,split=c(1,2,1,2),newpage=FALSE)
dev.off()

qplot(generation,log(c.est),data=cc,group=population,geom="line")
