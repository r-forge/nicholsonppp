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
dens.several.s <- function
### Plot density curves for ppp values for each selection type, to
### visualize the extent to which PPP-values can be used to indicate
### selection state in a simulation.
(df,
### Data frame with columns ppp type s to be plotted.
 ...
### To be passed to densityplot.
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
classify.loci <- function
### Classify loci into selection or not.
(res.df,
### Result of ppp.df.
 xmax=0.6,
 xlim=c(-0.15,xmax),
 ymax=21,
 ylim=c(-2,ymax),
 cutoff=seq(0,xmax,l=200)
### Vector of cutoff values to use for the classifier.
 ){
  classify.group <- function(sdf,cutoff=NULL){
    i <- which(sdf$ppp<xmax)
    if(is.null(cutoff))cutoff <- c(0,sort(sdf$ppp[c(i,max(i)+1)]))
    classify1 <- function(cutoff){
      summarise(transform(sdf,guess=ppp<cutoff,positive=type=="positive"),
                true.positive=sum(guess==TRUE&positive==TRUE),
                false.positive=sum(guess==TRUE&positive==FALSE),
                true.negative=sum(guess==FALSE&positive==FALSE),
                false.negative=sum(guess==FALSE&positive==TRUE),
                sensitivity=sum(guess==TRUE&positive==TRUE)/sum(positive),
                specificity=1-sum(guess==TRUE&positive==FALSE)/sum(!positive))
    }
    mdply(data.frame(cutoff),classify1)
  }
  ddply(res.df,.(s,desc),classify.group,cutoff)
}

cutoff.plot <- function
### Plot prediction counts against cutoff values.
(cl
### Data frame from classify.loci.
 ){
  N <- sum(cl[1,3:6])
  cl2 <- transform(cl,
                   correct=true.positive+true.negative,
                   incorrect=false.positive+false.negative)
  molt <- transform(melt(cl2,id=1:2),percent=value/N*100)
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
roc.loci <- function
### Plot ROC curves of loci classification for several s values.
(cl,
### Result of classify.loci.
 ...
### Arguments for xyplot.
 ){
  roc <- xyplot(sensitivity~1-specificity|desc,cl,
                type='l',
                groups=s,
                panel=function(...){
                  panel.abline(0,1,col="grey")
                  panel.xyplot(...)
                },
                main="ROCs for several selection strengths",
                ...)
  direct.label(roc,method=function(...){
    data.frame(perpendicular.lines(...),hjust=0)})
### The lattice plot.
}

