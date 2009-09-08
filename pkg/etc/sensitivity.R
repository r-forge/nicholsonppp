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
  print(desc)
  desc <- if(desc=="orig")"" else paste('.',desc,sep="")
  print(desc)
  lfile <- paste("sims",desc,".models.Rdata",sep="")
  simsname <- paste("sims",desc,sep="")
  modelsname <- paste("models",desc,sep="")
  load(file=lfile)
  sims <- get(simsname)
  models <- get(modelsname)
  df <- ppp.df(sims,models)
  subt <- deduce.param.label(do.call("cbind",sims[[1]]$p[display.params]))

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

a <- mdply(data.frame(desc=c("orig","few","neu")),loadpdf)
##TODO: ROC curves between the groups.
