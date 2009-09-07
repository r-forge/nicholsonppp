source("../R/sim.R")
source("../R/plot.R")
source("../R/model.R")
library(ggplot2)
library(lattice)
library(ff)
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
                                        s=max(sims[[i]]$s$s)))
  res.df$s <- factor(res.df$s,labels=format(unique(res.df$s),digits=2))
  attr(res.df,"N") <- sims[[1]]$p$n.loc
  res.df
### Data frame suitable for plotting densities or sensitivity.
}
density.several.s <- function
### Plot density curves for ppp values for each selection type, to
### visualize the extent to which PPP-values can be used to indicate
### selection state in a simulation.
(df,
### Data frame with columns ppp type s to be plotted.
 ...
### To be passed to dl.
 ){
  p <- densityplot(~ppp|s,df,groups=type,
                   layout=c(1,nlevels(df$s)),n=500,
                   strip=strip.custom(
                     strip.levels=c(TRUE,TRUE),strip.names=TRUE),
                   main="Small PPP-values indicate strong positive selection",
                   ...)
  direct.label(p)
### The lattice plot.
}


sims <- sim.several.s(s,loci=100,generations=100)
models <- nicholsonppp.list(sims)
##save(sims,models,file="sims.models.Rdata")
##load(file="sims.models.Rdata")
sims.neu <- sim.several.s(s,loci=100,generations=100,p.neutral=0.99,array.fun=ff)
models.neu <- nicholsonppp.list(sims.neu)
save(sims.neu,models.neu,file="sims.neu.models.Rdata")
##load(file="sims.models.Rdata")
library(latticedl)
neu.df <- ppp.df(sims.neu,models.neu)
subt <- deduce.param.label(do.call("cbind",sims.neu[[1]]$p[display.params]))
density.several.s(neu.df,sub=subt)

sims.few <- sim.several.s(s,loci=100,generations=100,populations=4)
models.few <- nicholsonppp.list(sims.few)
##save(sims.few,models.few,file="sims.few.models.Rdata")
load(file="sims.few.models.Rdata")
dfew <- ppp.df(sims.few,models.few)
library(latticedl)
subt <- deduce.param.label(do.call("cbind",sims.few[[1]]$p[display.params]))
density.several.s(dfew,sub=subt)
classify.loci(dfew,sub=subt)


pdf("ppp-density-several-s.pdf",paper="a4",h=0,w=0)
dev.off()

pdf("ppp-classify-several-s.pdf",paper="a4",h=0,w=0)
dev.off()

classify.loci <- function
### Classify loci into selection or not.
(res.df,
### Result of ppp.df.
 xmax=0.6,
 xlim=c(-0.15,xmax),
 ymax=21,
 ylim=c(-2,ymax),
 ...
### Other arguments for xyplot.
 ){
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
  molt <- transform(melt(cl2,id=1:2),percent=value/attr(res.df,"N")*100)
  minrisk.panel <- function(x,y,group.number,...){
    panel.xyplot(x=x,y=y,group.number=group.number,...)
    best.y <- min(y)
    best.x <- x[which(y==min(y))]
    if(group.number==2 & y[1]!=best.y){
      panel.abline(v=best.x,col="grey")
      panel.abline(h=best.y,col="grey")
      ## to use grid eventually:
      ltext(x=best.x,y=ymax,round(best.x,2),srt=90,adj=c(1,0),col="grey")
      ltext(x=xmax,y=best.y,round(best.y,2),adj=c(1,0),col="grey")
    }
  }
  ## cool but useless:
  ##dl(xyplot,molt,value~cutoff,variable,type='l')
  p2 <- xyplot(percent~cutoff|s,
               subset(molt,variable%in%c("true.positive","incorrect")),
               groups=variable,
               panel=minrisk.panel,
               type='l',layout=c(1,nlevels(molt$s)),
               xlim=xlim,ylim=ylim,
               xlab="Percent of loci in simulation",
               ylab="Cutoff for PPP-value decision rule",
               main="Prediction rates change with PPP-value cutoffs",
               ...)
  direct.label(p2)
### The lattice plot.
}


