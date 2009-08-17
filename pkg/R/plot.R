cex.title <- 1
pop.colors.default <- c("blue","turquoise","red")
selection.colors.default <- trellis.par.get("superpose.symbol")$col
selection.colors.default[3] <- "green"
selection.symbols.default <- c("B","N","P")
display.params <- c("populations",
                    'loci.per.s.value',
                    'popsize',
                    'generations',
                    'p.neutral',
                    'fixed')
deduce.param.label <- function(lt,imp=display.params){
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
 pop.colors=pop.colors.default
### Colors to use to distinguish populations (blue, neutral, red)
 ){
  ggplot(fr,aes(generation,simulated,group=population,colour=color))+
    geom_line()+
    facet_wrap(S~locus,nrow=1)+
    scale_colour_manual(values=pop.colors)+
    labs(y="Simulated blue allele frequency",colour="Population color")+
    opts(title="Allele frequency evolution varies with selection type and population color")
}


## conversion to format compatible with endpoint allele frequency
## fixation plots
convert.df <- function(d,g){
  ## with lattice the easiest way to do legends is just with auto.key
  ## so in order to have ancestral allele frequencies lets add some
  ## rows to the data frame
  d <- subset(d,generation==g)
  ancest <- d
  ancest$color <- 'ancestral'
  ancest$simulated <- ancest$ancestral
  ancest$population <- 0
  ancest <- unique(ancest)
  d <- rbind(ancest,d)
  d$locus <- factor(d$locus)
  print(levels(d$locus))
  d$locus <- reorder(d$locus,d$ancestral,median)
  print(levels(d$locus))
  d
}
fixation.endpoints <- function
### Plot gene frequencies for all loci and populations for a given
### generation.
(lf,
### Subset of simulated gene frequency data frame, with just 1
### generation.
 main="Loci are not always fixed for the same allele in each population",
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
  mydot.panel <- function(...,x,subscripts){
    panel.xyplot(x=x,subscripts=subscripts,...)
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


myanim <- function(f,tit,pfun){
  naivedir <- paste(getwd(),subdir,f,sep="/")
  dir.create(naivedir)
  ani.start(nmax=generations,
            title=tit,
            outdir=naivedir,
            ani.width=1200,
            ani.height=800)
  for(g in 1:ani.options("nmax")){
    cat(naivedir,g,"\n")
    ## this actually plots for the less complicated plots
    ## but does nothing for the more complicated bigplot()
    print(pfun(g)) 
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
         auto.key=list(space="right",title="Selection type",cex.title=cex.title))
}
