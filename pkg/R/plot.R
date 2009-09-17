### Text size for lattice legend titles.
cex.title <- 1

### Default colors used for representing population types.
pop.colors.default <- c("blue","turquoise","red")

### Default colors used for selection types (balancing, neutral, positive).
selection.colors.default <- c("#0080ff","#ff00ff","green")

### Default symbols used for selection types (balancing, neutral, positive).
selection.symbols.default <- c("B","N","P")

### Default simulation parameters summarized in plot subtitles.
display.params <- c("populations",
                    'popsize',
                    'generations',
                    'n.locus',
                    'p.neutral',
                    "s",
                    'loci.per.s.value')

deduce.param.label <- function # Deduce parameters for plot subtitles
### Summarize parameter values used in simulations.
(lt,
### Data frame of results of simulation.
 imp=display.params
### Vector of column names to be searched and reported.
 ){
  lt <- as.data.frame(lt)
  imp <- imp[imp%in%names(lt)]
  L <- sapply(imp,function(cn)if(cn%in%names(lt)){
    tmp <- unique(lt[[cn]])
    if(is.numeric(tmp))tmp <- round(tmp,digits=2)
    tmp
  } else NULL)
  lab <- paste(sapply(names(L),function(N)paste(N,
                                         ": ",
                                         paste(L[[N]],collapse=", "),
                                         "  ",
                                         sep='')),
        collapse='')
  if(lab=="")NULL else lab
### Text string, or NULL, to use as sub= argument for a high-level
### lattice function.
}

interesting.loci <- function
### Find a small subset of loci for each selection type, from a
### variety of ancestral frequencies. This will be used to show how
### allele frequency evolution depends on selection type and
### population color.
(fr
### Data frame of all the simulated allele frequencies.
 ){
  to.df <- function(subs){
    subs <- subset(subs,s==max(s))
    rows <- with(subs,c(which(d<0.1)[1],which(d>0.3&d<0.4)[1]))
    rows <- rows[!is.na(rows)]
    if(length(rows)<2)rows <- 1:2
    loci <- subs[rows,"locus"]
    subset(all.fr,locus%in%loci)
  }
  all.fr <- fr
  n.locus <- attr(fr,"parameters")$n.locus
  ## Optimization when fr is a data frame with repeated loci:
  if(!is.null(n.locus))fr <- head(fr,n.locus)
  anc.diff <- transform(fr,d=abs(ancestral-0.5))
  res <- ddply(anc.diff,.(type),to.df)
  res
### Data frame, subset of input data.
}

loci.over.time <- function # Allele frequency time series plots
### Plot allele frequency evolution over time for each locus, grouping
### by population color. This shows the difference between selection
### types and population colors.
(fr,
### Data frame of simulated allele frequencies.
 pop.colors=pop.colors.default,
### Colors to use to distinguish populations (blue, neutral, red)
 generation=NULL,
### Generation to emphasize, or NULL for no emphasis.
 m="Allele frequency evolution varies with selection type and population color",
### Main title.
 popsize.enc=NULL
### Aesthetic to encode population size. lty or lwd work well. NULL
### means do not show population size.
 ){
  p <- ggplot(fr,aes(generation,simulated,group=population,colour=color))+
    geom_line()+
    ylim(c(0,1))+
    facet_wrap(S~locus,nrow=1)+
    scale_colour_manual(values=pop.colors)+
    labs(y="Simulated blue allele frequency",
         colour="Population color")+
    opts(title=m,
         panel.background=theme_rect(colour=NA),
         strip.background=theme_rect(colour=NA),
         legend.key=theme_rect(colour=NA))
  if(!is.null(generation)){
    p <- p +
      geom_vline(xintercept=generation)+
      geom_hline(aes(yintercept=simulated[1]),colour="grey")
  }
  if(!is.null(popsize.enc)){
    mf <- function(fun.name,argval){
      L <- list(as.name(fun.name),argval)
      names(L)[2] <- popsize.enc
      cl <- as.call(L)
      eval(cl)
    }
    p <- p+
      mf("aes",call("factor",as.name("popsize")))+
      mf("labs","Population size")
  }
  p
### The ggplot2 plot.
}

fixation.endpoints <- function # Plot allele frequency by selection type
### Plot gene frequencies for all loci and populations for a given
### generation, stratified by selection type and coefficient.
(lf,
### Subset of simulated gene frequency data frame, with just 1
### generation.
 main="Loci fixation depends on selection type and population color",
### Plot title.
 pop.colors=c(pop.colors.default,"black"),
### Population color scheme (blue, neutral, red, ancestral).
 par.settings=list(superpose.symbol=list(col=pop.colors,pch=20)),
### Plot settings as described in trellis.par.get.
 hilite.locus=NULL,
### Locus to highlight with a vertical line, or NULL to highlight
### nothing.
 selection.colors=selection.colors.default,
### List with element "col" which contains a vector of colors to label
### the selection types (balancing, none, positive)
 sub=deduce.param.label(attr(lf,"parameters")),
### Subtitle for the plot.
 ...
### Other parameters to pass to xyplot.
 ){
  lf$locus.ord <- reorder(as.factor(lf$locus),lf$ancestral)
  mapfrom <- levels(lf$locus.ord)
  levels(lf$locus.ord) <- 1:nlevels(lf$locus.ord)
  map <- structure(levels(lf$locus.ord),names=mapfrom)
  anc <- ddply(lf,.(locus),summarise,
               simulated=ancestral[1],
               color="ancestral",
               S=S[1],
               type=type[1],
               locus.ord=locus.ord[1])
  anc.sim <- rbind(anc,lf[,names(anc)])
  anc.sim$color <- reorder(anc.sim$color,anc.sim$color,
                           function(x)if(x[1]=="ancestral")1 else 0)
  hilite.locus.ord <- as.integer(map[paste(hilite.locus)])
  mydot.panel <- function(x,...){
    panel.xyplot(x=x,...)
    if((!is.null(hilite.locus))&&hilite.locus.ord%in%x)
      panel.abline(v=hilite.locus.ord)
  }
  xyplot(simulated~locus.ord|S,anc.sim,
         alpha=1,
         panel=mydot.panel,
         auto.key=list(space='right',title="Population type",
           cex.title=cex.title),
         xlab="Locus (a dot is drawn for each population and locus)",
         ylab="Simulated blue allele frequency",
         ylim=c(0,1),
         groups=color,
         main=main,
         sub=sub,
         par.settings=par.settings,
         scales=list(x=list(draw=FALSE)),
         layout=c(nlevels(anc.sim$S),1),
         strip=function(which.panel,factor.levels,bg,...){
           level <- gsub(" .*$","",factor.levels[which.panel])
           supcol <- level==levels(anc.sim$type)
           strip.default(which.panel=which.panel,
                         factor.levels=factor.levels,
                         bg=selection.colors[supcol],
                         strip.names=c(TRUE,TRUE),
                         strip.levels=c(TRUE,TRUE),
                         ...)
         },
         ...)
### The lattice plot.
}

anc.est.naive <- function # Naive ancestral allele frequency estimates
### Estimate ancestral allele frequencies by simply taking the mean of
### all the current observations. This is a simple statistic we can
### use to reality check the model estimates, and which are fast
### enough to calculate for each generation in the simulation,
### facilitating animations.
(fr
### Subset of simulation data frame containing the generation of
### interest from which we estimate the ancestral frequency.
 ){
  est.df <- ddply(fr,.(locus),summarise,
                  ancestral.est=mean(simulated),
                  ancestral=ancestral[1],
                  type=type[1])
  attr(est.df,"parameters") <- attr(fr,"parameters")
  est.df
### Data frame with a row for each locus, summarizing naive ancestral
### allele frequencies in the ancestral.est column.
}

anc.est.nicholson <- function
### Ancestral allele frequency estimate given by Nicholson model.
(sim,
 model
 ){
  data.frame(sim$s,ancestral.est=model$p)
}

anc.est.plot <- function # Ancestral estimate plot
### Plot naive estimates of ancestral allele frequency versus actual
### values from the simulation, to see if they agree.
(est.df,
### Data frame describing loci and ancestral estimates. Need columns
### type ancestral ancestral.est.
 hilite.locus=NULL,
### Locus number to highlight on the plot with a circle. NULL means do
### not highlight.
 sub=deduce.param.label(attr(est.df,"parameters")),
### Subtitle for the plot.
 selection.colors=selection.colors.default,
### Colors for the different selection types (balancing, neutral, positive)
 selection.symbols=selection.symbols.default,
### Symbols for the different selection types (balancing, neutral, positive)
 f=ancestral.est~ancestral,
### Plot formula for xyplot.
 xlab="Simulated blue allele frequency",
 ylab="Estimated blue allele frequency",
 main="Allele frequency estimates vary with selection type"
 ){
  xyplot(f,est.df,
         alpha=1,
         cols=selection.colors,
         groups=type,
         panel=function(...){
           panel.xyplot(...)
           if(!is.null(hilite.locus)){
             pp <- est.df[hilite.locus,]
             lpoints(pp$ancestral,pp$ancestral.est,pch=1,cex=2,col="black")
           }
         },
         par.settings=list(superpose.symbol=list(
                             pch=selection.symbols,
                             cex=1.2,
                             col=selection.colors)),
         strip=strip.custom(strip.names=TRUE),
         sub=sub,
         ylim=c(0,1),
         xlim=c(0,1),
         aspect=1,
         xlab=xlab,
         ylab=ylab,
         main=main,
         auto.key=list(space="right",
           title="Selection type",cex.title=cex.title))
### The lattice plot.
}

sim.summary.plot <- function # Simulation summary plot
### Draw 3 simulation summary plots on the same screen
### (loci.over.time, anc.est.plot, fixation.endpoints).
(fr,
### Data frame from a simulation.
 g=NULL,
### Generation to plot in anc.est.plot and fixation.endpoints, and
### generation to emphasize in loci.over.time. NULL means the last
### generation in the simulation.
 hilite.locus=NULL
### Locus to highlight in the plots. NULL means choose an interesting
### locus under positive selection.
 ){
  if(is.null(g))g <- attr(fr,"parameters")$gen
  if(is.null(hilite.locus)){
    fr.i <- interesting.loci(fr)
    hilite.locus <- tail(fr.i,1)$locus
  }
  parameters <- attr(fr,"parameters")
  ss <- fr[fr$generation==g,]
  vpl <- function(x,y)viewport(layout.pos.row=x,layout.pos.col=y)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,2)))
  p <- loci.over.time(fr[fr$locus==hilite.locus,],
                      generation=g,
                      m="Allele frequency evolution for 1 locus")
  print(p,vp=vpl(1,1))
  aep <- anc.est.plot(anc.est.naive(ss),hilite.locus=hilite.locus,sub=NULL)
  pushViewport(vpl(1,2));print(aep,newpage=FALSE);popViewport()
  fe <- fixation.endpoints(ss,hilite.locus=hilite.locus)
  pushViewport(vpl(2,1:2));print(fe,newpage=FALSE);popViewport()
}

evolution.animation <- function
### Create an animation that summarizes a simulation, using
### sim.summary.plot for every generation of the simulation.
(df,
### Result of sim2df.
 outdir=tempfile(),
### Subdirectory for plot files, to be created.
 tit="Allele frequency and ancestral estimate evolution"
### Title for the animation.
 ){
  dir.create(outdir)
  ## convert relative path to full path (animation package bug)
  outdir <- tools::file_path_as_absolute(outdir)
  gens <- unique(df$generation)
  ani.start(nmax=length(gens),
            title=tit,
            outdir=outdir,
            ani.width=1000,
            ani.height=800)
  on.exit(ani.stop())
  loci <- interesting.loci(df)
  for(g in gens){
    cat(outdir,g,"\n")
    sim.summary.plot(df,g)
  }
}

