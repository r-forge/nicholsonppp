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
  splitby <- c("s")
  if("desc"%in%names(res.df))splitby <- c(splitby,"desc")
  ddply(res.df,splitby,classify.group,cutoff)
}

hilite.best <- function
### panel.groups function for highlighting the best points by drawing
### grey lines and labeling the actual values.
(x,
 y
 ){
  panel.abline(v=x,col="grey")
  panel.abline(h=y,col="grey")

  first.x <- x[1]
  grid.text(round(first.x,2),unit(first.x,"native"),1,
            rot=90,just=c("right","bottom"),gp=gpar(col="grey"))
  
  if(length(x)>1){
    last.x <- x[length(x)]
    grid.text(round(last.x,2),unit(last.x,"native"),1,
              rot=90,just=c("right","top"),gp=gpar(col="grey"))
  }
  
  grid.text(round(y,2),1,unit(y,"native"),
            just=c("right","bottom"),gp=gpar(col="grey"))
}

panel.densityplot.offset <- function
(x, darg = list(n = 30), plot.points = "jitter", ref = FALSE, 
 groups = NULL, weights = NULL,
 jitter.amount = 0.01 * diff(current.panel.limits()$ylim), 
 type = "p",
 ...){
  if (ref) {
    reference.line <- trellis.par.get("reference.line")
    panel.abline(h = 0, col = reference.line$col, lty = reference.line$lty, 
                 lwd = reference.line$lwd)
  }
  plot.line <- trellis.par.get("plot.line")
  superpose.line <- trellis.par.get("superpose.line")
  if (!is.null(groups)) {
    panel.superpose(x, darg = darg, plot.points = plot.points, 
                    ref = FALSE, groups = groups, weights = weights, 
                    panel.groups = panel.densityplot,
                    jitter.amount = jitter.amount, 
                    type = type, ...)
  }
  else {
    switch(as.character(plot.points),
           `TRUE`=panel.xyplot(x = x,y = rep(0, length(x)), type = type, ...),
           rug = panel.rug(x = x,start = 0, end = 0, x.units = c("npc", "native"), 
                                                                                    type = type, ...), jitter = panel.xyplot(x = x, y = jitter(rep(0, 
                                                                                                                                      length(x)), amount = jitter.amount), type = type, 
                                                                                                         ...))
    density.fun <- function(x, weights, subscripts = TRUE, 
                            darg, ...) {
      do.call("density", c(list(x = x, weights = weights[subscripts]), 
                           darg))
    }
    if (sum(!is.na(x)) > 1) {
      h <- density.fun(x = x, weights = weights, ..., darg = darg)
      lim <- current.panel.limits()$xlim
      id <- h$x > min(lim) & h$x < max(lim)
      panel.lines(x = h$x[id], y = h$y[id], ...)
    }
  }
}


cutoff.plot <- function
### Plot prediction counts against cutoff values.
(cl,
### Data frame from classify.loci.
 ylim=c(-3,21),
### Limits for y axis.
 xlim=c(-0.15,0.6),
### Limits for x axis.
 dens=NULL,
### Optional density data to plot.
 ...
### args for xyplot.
 ){
  N <- sum(cl[1,3:6])
  cl2 <- transform(cl,
                   correct=true.positive+true.negative,
                   incorrect=false.positive+false.negative)
  molt <- transform(melt(cl2,id=1:2),percent=value/N*100)
  minrisk.panel <- function(x,y,subscripts,group.number,...){
    panel.xyplot(x=x,y=y,group.number=group.number,...)
    if(group.number==2){
      best.y <- min(y)
      if(y[1]!=best.y){
        best.x <- x[which(y==min(y))]
        hilite.best(best.x,best.y)
      }
    }
  }
  ## cool but useless:
  ##dl(xyplot,molt,value~cutoff,variable,type='l')
  ss <- subset(molt,variable%in%c("true.positive","incorrect"))
  ss$variable <- factor(ss$variable)
  p2 <- xyplot(percent~cutoff|s,ss,
               groups=variable,
               panel=function(subscripts,groups,...){
                 panel.superpose(subscripts=subscripts,groups=groups,...)
                 if(!is.null(dens)){
                   d <- subset(dens,s==ss[subscripts,"s"][1])
                   panel.superpose(d$ppp,y=NULL,1:nrow(d),d$type,
                                   panel.groups=panel.rug,
                                   end=0.05,
                                   col=selection.colors.default)
                 }
               },
               panel.groups=minrisk.panel,
               type='l',layout=c(1,nlevels(molt$s)),
               xlim=xlim,
               ylim=ylim,
               ylab="Percent of loci in simulation",
               xlab="Cutoff for PPP-value decision rule",
               main="Prediction rates change with PPP-value cutoffs",
               par.settings=list(superpose.line=list(col=c("brown","orange"))),
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
                main="ROCs vary with selection strength",
                ...)
  direct.label(roc,method=function(...){
    x <- data.frame(perpendicular.lines(...),hjust=0)
    x[2,c("x","y")] <- c(0.6,0.8)
    x
  })
### The lattice plot.
}

