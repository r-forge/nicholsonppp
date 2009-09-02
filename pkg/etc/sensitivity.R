source("../R/sim.R")
source("../R/plot.R")
source("../R/model.R")
library(ggplot2)
library(lattice)
dyn.load("../src/0mers_twist.so")
## initial study of which s values show a good range of selection
## behavior
s <- 10^c(-3,-2,-1.5,-1,0)
sim <- sim.drift.selection(loci=100,generations=100,s=s,p.neutral=0.1)
df <- sim2df(sim)
sim.summary.plot(df)

sim.several.s <- function
### Do several simulations each with only 1 s value.
(s,
### Vector of s values.
 ...
### Other arguments to sim.drift.selection.
 ){
  mlply(data.frame(s),function(s)sim.drift.selection(s=s,...))
### List of selection simulations.
}
nicholsonppp.list <- function
### Fit a model for each simulation.
(sims
### List of simulations from nicholsonppp.
 ){
  lapply(sims,function(L)nicholsonppp(L$sim[,,L$p$gen]))
### List of nicholson model fits.
}

sims <- sim.several.s(s,loci=100,generations=100)
models <- nicholsonppp.list(sims)
##save(sims,models,file="sims.models.Rdata")
##load(file="sims.models.Rdata")

sims.neu <- sim.several.s(s,loci=100,generations=100,p.neutral=0.99)
models.neu <- nicholsonppp.list(sims.neu)
##save(sims,models,file="sims.models.Rdata")
##load(file="sims.models.Rdata")

sims.few <- sim.several.s(s,loci=100,generations=100,populations=4)
models.few <- nicholsonppp.list(sims.few)
##save(sims,models,file="sims.models.Rdata")
##load(file="sims.models.Rdata")


ppp.df <- function
### Convert simulation and model fit lists to a data frame that can be
### used for plotting diagnostics for the ppp-value classifier.
(sims,
### List of simulations.
 models
### List of nicholson model fits.
 ){
  res.df <- ldply(1:length(models),
                  function(i)data.frame(ppp=models[[i]]$ppp,
                                        type=sims[[i]]$s$type,
                                        s=factor(s[i])))
  levels(res.df$s) <- format(s,digits=2)
  res.df
### Data frame suitable for plotting densities or sensitivity.
}



## Plotting code, to be encapsulated in R functions:
library(latticedl)
subt <- deduce.param.label(do.call("cbind",sims[[1]]$p[display.params]))
pdf("ppp-density-several-s.pdf",paper="a4",h=0,w=0)
p1 <- dl(densityplot,res.df,~ppp|s,type,
         layout=c(1,length(s)),n=500,sub=subt,
         strip=strip.custom(strip.levels=c(TRUE,TRUE),strip.names=TRUE),
         main="Small PPP-values indicate strong positive selection")
plot(p1)
dev.off()
densityplot(~ppp|s,res.df,layout=c(1,length(s)),n=500)

xmax <- 0.6
ymax <- 21
classify.group <- function(sdf){
  i <- which(sdf$ppp<xmax)
  cutoff <- c(0,sort(sdf$ppp[c(i,max(i)+1)]))
  classify1 <- function(cutoff){
    summarise(transform(sdf,guess=ppp<cutoff,positive=type=="positive"),
              true.positive=sum(guess==TRUE&positive==TRUE),
              false.positive=sum(guess==TRUE&positive==FALSE),
              true.negative=sum(guess==FALSE&positive==FALSE),
              false.negative=sum(guess==FALSE&positive==TRUE))
  }
  mdply(data.frame(cutoff),classify1)
}
cl <- ddply(res.df,.(s),classify.group)
cl2 <- transform(cl,
                 correct=true.positive+true.negative,
                 incorrect=false.positive+false.negative)
molt <- transform(melt(cl2,id=1:2),percent=value/sims[[1]]$p$n.loc*100)
minrisk.panel <- function(x,y,group.number,...){
  panel.xyplot(x=x,y=y,group.number=group.number,...)
  best.y <- min(y)
  best.x <- x[which(y==min(y))]
  if(group.number==2 & y[1]!=best.y){
    panel.abline(v=best.x,col="grey")
    panel.abline(h=best.y,col="grey")
    ltext(x=best.x,y=ymax,round(best.x,2),srt=90,adj=c(1,0),col="grey")
    ltext(x=xmax,y=best.y,round(best.y,2),adj=c(1,0),col="grey")
  }
}
## cool but useless:
##dl(xyplot,molt,value~cutoff,variable,type='l')
p2 <- dl(xyplot,subset(molt,variable%in%c("true.positive","incorrect")),
         percent~cutoff|s,variable,
         panel=minrisk.panel,
         type='l',layout=c(1,nlevels(molt$s)),
         xlim=c(-0.15,xmax),ylim=c(-2,ymax),
         xlab="Percent of loci in simulation",
         ylab="Cutoff for PPP-value decision rule",
         main="Prediction rates change with PPP-value cutoffs",
         sub=subt)
pdf("ppp-classify-several-s.pdf",paper="a4",h=0,w=0)
plot(p2)
dev.off()

plot(p1,split=c(1,1,2,1))
plot(p2,split=c(2,1,2,1),newpage=FALSE)











