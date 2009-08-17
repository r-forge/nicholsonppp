### Text size for lattice legend titles.
cex.title <- 1

### Default colors used for representing population types.
pop.colors.default <- c("blue","turquoise","red")

### Default colors used for selection types.
selection.colors.default <- trellis.par.get("superpose.symbol")$col
selection.colors.default[3] <- "green"

### Default symbols used for selection types.
selection.symbols.default <- c("B","N","P")

### Default simulation parameters summarized in plot subtitles.
display.params <- c("populations",
                    'loci.per.s.value',
                    'popsize',
                    'generations',
                    'p.neutral',
                    'fixed')

deduce.param.label <- function
### Summarize parameter values used in simulations.
(lt,
### Data frame of results of simulation.
 imp=display.params
### Vector of column names to be searched and reported.
 ){
  L <- sapply(imp,function(cn)if(cn%in%names(lt)){
    tmp <- unique(lt[,cn])
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
}

interesting.loci <- function
### Find a small subset of loci for each selection type, from a
### variety of ancestral frequencies. This will be used to show how
### allele frequency evolution depends on selection type on population
### color.
(fr
### Data frame of all the simulated allele frequencies.
 ){
  to.df <- function(subs){
    subs <- subset(subs,s==max(s))
    rows <- c(which(subs$d<0.1)[1],which(subs$d>0.3 & subs$d<0.4)[1])
    rows <- rows[!is.na(rows)]
    loci <- subs[rows,"locus"]
    subset(subs,locus%in%loci)
  }
  res <- ddply(transform(fr,d=abs(ancestral-0.5)),.(type),to.df)
  res[order(res$generation),]
### Data frame, subset of input data.
}

loci.over.time <- function
### Plot allele frequency evolution over time for each locus, grouping
### by population color. This shows the difference between selection
### types and population colors.
(fr,
### Data frame of simulated allele frequencies.
 pop.colors=pop.colors.default,
### Colors to use to distinguish populations (blue, neutral, red)
 generation=NULL,
### Generation to emphasize, or NULL for no emphasis.
 m="Allele frequency evolution varies with selection type and population color"
### Main title.
 ){
  p <- ggplot(fr,aes(generation,simulated,group=population,colour=color))+
    geom_line()+
    ylim(c(0,1))+
    facet_wrap(S~locus,nrow=1)+
    scale_colour_manual(values=pop.colors)+
    labs(y="Simulated blue allele frequency",colour="Population color")+
    opts(title=m)
  if(!is.null(generation))p <- p + geom_vline(xintercept=generation)
  p
}

fixation.endpoints <- function
### Plot gene frequencies for all loci and populations for a given
### generation.
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
 sub=deduce.param.label(lf),
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
         scales=list(x=list(draw=F)),
         layout=c(5,1),
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
}

evolution.animation <- function(f,tit,df,...){
  naivedir <- paste(getwd(),f,sep="/")
  dir.create(naivedir)
  ani.start(nmax=max(df$generations),
            title=tit,
            outdir=naivedir,
            ani.width=1000,
            ani.height=800)
  for(g in 1:ani.options("nmax")){
    cat(naivedir,g,"\n")
    ## this actually plots for the less complicated plots
    ## but does nothing for the more complicated bigplot()
    bigplot(df,g,...)
  }
  ani.stop()
}


anc.est.plot <- function
### Plot naive estimates of ancestral allele frequency versus actual
### values from the simulation, to see if they agree.
(fr,
### Subset of simulation data frame containing the generation of
### interest to estimate and plot.
 hilite.locus=NULL,
### Locus number to highlight on the plot with a circle. NULL means do
### not highlight.
 sub=deduce.param.label(fr),
### Subtitle for the plot.
 selection.colors=selection.colors.default,
### Colors for the different selection types (balancing, neutral, positive)
 selection.symbols=selection.symbols.default
### Symbols for the different selection types (balancing, neutral, positive)
 ){
  est.df <- ddply(fr,.(locus),summarise,
                  ancestral.est=mean(simulated),
                  ancestral=ancestral[1],
                  type=type[1])
  xyplot(ancestral.est~ancestral,est.df,
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
         sub=sub,
         ylim=c(0,1),
         xlim=c(0,1),
         aspect=1,
         xlab="Simulated blue allele frequency",
         ylab="Estimated blue allele frequency",
         main="Allele frequency estimates vary with selection type",
         auto.key=list(space="right",
           title="Selection type",cex.title=cex.title))
}

bigplot <- function
### Draw 3 simulation summary plots on the same screen
### (loci.over.time, anc.est.plot, fixation.endpoints).
(fr,
### Data frame from a simulation.
 g=1,
### Generation to plot in anc.est.plot and fixation.endpoints, and
### generation to emphasize in loci.over.time.
 hilite.locus=1
### Locus to highlight in the plots.
 ){
  ss <- subset(fr,generation==g)
  vpl <- function(x,y)viewport(layout.pos.row=x,layout.pos.col=y)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,2)))
  p <- loci.over.time(subset(fr,locus==hilite.locus),
                      generation=g,
                      m="Allele frequency evolution for 1 locus")
  print(p,vp=vpl(1,1))
  aep <- anc.est.plot(ss,hilite.locus=hilite.locus,sub=NULL)
  pushViewport(vpl(1,2));print(aep,newpage=FALSE);popViewport()
  fe <- fixation.endpoints(ss,hilite.locus=hilite.locus)
  pushViewport(vpl(2,1:2));print(fe,newpage=FALSE);popViewport()
}

