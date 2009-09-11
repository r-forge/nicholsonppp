source("../R/sim.R")
source("../R/plot.R")
source("../R/sensitivity.R")
source("../R/model.R")
library(ggplot2)
library(lattice)
library(ff)
library(latticedl)
dyn.load("../src/0mers_twist.so")
## initial study of which s values show a good range of selection
## behavior
s <- 10^c(-3,-2,-1.5,-1,0)
sim <- sim.drift.selection(loci=100,generations=100,s=s,p.neutral=0.1)
df <- sim2df(sim)
sim.summary.plot(df)



sims <- sim.several.s(s,loci=100,generations=100)
models <- nicholsonppp.list(sims)
##save(sims,models,file="sims.models.Rdata")
load(file="sims.models.Rdata")



sims.neu <- sim.several.s(s,loci=100,generations=100,
                          p.neutral=0.99,array.fun=ff)
models.neu <- nicholsonppp.list(sims.neu)
##save(sims.neu,models.neu,file="sims.neu.models.Rdata")

loadpdf <- function
### Load a data set and plot it.
(desc="",
### Data set description.
 viewer="xpdf",
### program to view output.
 ...
### To be passed to classify.loci.
 ){
  if(desc!="")desc <- paste('.',desc,sep="")
  lfile <- paste("sims",desc,".models.Rdata",sep="")
  simsname <- paste("sims",desc,sep="")
  modelsname <- paste("models",desc,sep="")
  load(file=lfile)
  sims <- get(simsname)
  models <- get(modelsname)
  df <- ppp.df(sims,models)
  subt <- deduce.param.label(do.call("cbind",sims[[1]]$p[display.params]))
  print(subt)

  makepdf <- function(fun.name){
    plot.fun <- get(fun.name)
    outf <- paste(fun.name,desc,".pdf",sep="")
    pdf(outf,paper="a4",h=0,w=0)
    print(plot.fun(df,sub=subt,...))
    dev.off()
    cmd <- paste(viewer,outf,"&")
    system(cmd)
  }

  ##makepdf("dens.several.s")
  ##makepdf("classify.loci")

  df
}


## loadpdf()
## loadpdf("few")
## loadpdf("neu",ymax=2)

a <- mdply(data.frame(desc=c("","few","neu"),viewer="none"),loadpdf)
levels(a$desc) <- c("12 populations, 1000 loci",
                    "4 populations, 1000 loci",
                    "12 populations, 19999 loci")
acl <- classify.loci(a)

pdf("roc-desc.pdf",h=6,w=8.5)
roc.loci(acl,layout=c(3,1),aspect=1)
dev.off()

pdf("roc-s.pdf",h=0,w=0,paper="a4")
xyplot(sensitivity~1-specificity|s,acl,
       type='l',
       groups=desc,
       panel=function(...){
         panel.abline(0,1,col="grey")
         panel.xyplot(...)
       },
       main="ROCs vary with selection strength and number of populations",
       auto.key=list(title="Simulation",space="right",lines=TRUE,points=FALSE),
       layout=c(1,5),
       aspect=1)
dev.off()
