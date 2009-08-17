simulate.drift.selection <- function
### Simulate allele frequency evolution in several populations
### according to a simple drift and selection model.
(populations=8,
### Number of populations.
 popsize=1000,
### Size of each simulated population.
 generations=25,
### Generations of evolution to simulate.
 loci.per.s.value=50,
### Number of loci simulated for each s value.
 beta1=0.7,
### Parameter for beta distribution of initial allele frequencies.
 beta2=0.7,
### Parameter for beta distribution of initial allele frequencies.
 p.neutral=1/3,
### Proportion of neutral alleles to simulate.
 s,
### Vector of selection strength values to simulate.
 adapt.pop=c(blue=0.4,red=0.4,neutral=0.2)
### Probabilities of color adaptation in populations.
 ){
  timestr <- gsub('[ :]','-',format(Sys.time()))
  subdir <- paste('results',timestr,sep='/')
  rfile <- function(name)paste(subdir,'/',name,'.txt',sep='')
  rpdf <- function(name)
    pdf(paste(subdir,'/',timestr,'-',name,'.pdf',sep=''),h=10,w=7)
  results <- function(name)
    write.table(get(name),rfile(name),quote=FALSE,row.names=TRUE,col.names=TRUE)
  dir.create(subdir)

  srep <- rep(s,each=loci.per.s.value)
  Nloc <- length(srep)
  n.neutral <- 2*Nloc*p.neutral/(1-p.neutral)
  s.LOC=c(rep(0,n.neutral),srep,srep)
  Type.LOC=factor(c(
    rep("none",n.neutral),
    rep("positive",Nloc),
    rep("balancing",Nloc)))
  Nlocus=length(Type.LOC)
  ## choice: diminishes the visible effects of the simulation
  ## three primary choices: bell curve, U shape, uniform.
  ## show the probability density used to create initial F values
  rpdf('beta-density')

  ##Nlocus=1000;beta1=0.7;beta2=0.7
  F.init=F.init.noncoup=rbeta(Nlocus,beta1,beta2)
  x <- seq(0,1,l=1e3)
  x <- x[2:(length(x)-1)]
  y <- dbeta(x,beta1,beta2)
  par(xpd=NA,las=1,bty='n',xaxs='i',mfrow=c(2,1),mar=c(2,4,4,2)+0.1)
  plot(x,y,type='n',ylim=c(0,max(y)),
       ylab="Probability density",xlab='',xlim=c(0,1)
       )
  title(substitute(paste(beta,
                         group("(", list(beta1, beta2), ")"),
                         " probabilities"),
                   list(beta1 = beta1, beta2 = beta2)))
  cutoff <- 0.05
  F.init[F.init<cutoff] <- cutoff
  F.init[F.init>(1-cutoff)] <- 1-cutoff
  hist(F.init,add=T,prob=T,col='grey',breaks=20)
  hist(F.init.noncoup,add=T,prob=T,col='white',breaks=20)
  lines(x,y)
  mtext(paste(length(F.init),'simulated values of allele frequency'))

  par(mar=c(3,4,3,2)+0.1)
  plot(x <- seq(cutoff,1-cutoff,l=1e2),
       y <- pbeta(x,beta1,beta2),
       type='l',
       xlim=c(0,1),
       ylim=c(0,1),
       xlab="",
       xaxt='n',
       yaxs='i',
       ylab="Cumulative distribution function")
  segments(0,0,x[1],0)
  segments(x[length(x)],1,1,1)
  par(pch=1)
  points(x[1],0)
  points(x[length(x)],y[length(y)])
  par(pch=20)
  points(x[1],y[1])
  points(x[length(x)],1)
  par(las=3)
  axis(1,c(0,x[1]),lwd=0)
  axis(3,c(1,x[length(x)]),lwd=0)
  sorted <- sort(F.init)
  points(cbind(x=unique(F.init),
               y=sapply(unique(F.init),
                 function(v)sum(sorted<=v)/length(F.init))),
         pch=".",
         col='grey',
         cex=2)
  legend('bottomright',
         legend=c('truncated','usual'),
         fill=c('grey','white'),
         bty='n',bg='white')
  legend('top',
         legend='theoretical',
         lty=1,
         bty='n',
         bg='white')

  
  dev.off()

  s <- data.frame(s=s.LOC,type=Type.LOC,locus=1:length(s.LOC))
  sfac <- factor(paste(s$type,format(s$s)))
  levs <- levels(sfac)
  ids <- rev(grep("balancing",levs))
  levs <- levs[c(ids,(ids[1]+1):length(levs))]
  s$S <- factor(sfac,levs)

  ## each population/locus combination has a different "color"
  Type.POP <- matrix(sample(names(adapt.pop),populations*Nlocus,T,adapt.pop),
                     Nlocus,populations)

  ## this will store all of the simulated allele frequencies
  fr <- array(F.init,c(Nlocus,populations,generations))
  ## these will be intermediate data structures to support fast selection
  selected.loci <- matrix(as.character(Type.LOC)!="none",Nlocus,populations)
  selective.pops <- Type.POP!='neutral'
  selected <- selected.loci & selective.pops
  smat <- matrix(s.LOC,Nlocus,populations)[selected]
  type <- matrix(Type.LOC,Nlocus,populations)[selected]
  popcol <- Type.POP[selected]
  ## Genotype weights for each locus and population under selection:
  ## (not including neutral loci)
  BB <- ifelse(type=="positive"&popcol=="blue",1+smat  ,1)
  BR <- ifelse(type=="positive",               1+smat/2,1+smat)
  RR <- ifelse(type=="positive"&popcol=="red", 1+smat  ,1)
  ## diagnostic data frame
  ## data.frame(BB,BR,RR,popcol,s,type)
  for(t in 2:generations){
    cat("Generation: ",t,"\n")
    ## update all loci for drift
    fr[,,t] <- rbinom(populations*Nlocus,popsize,fr[,,t-1])/popsize
    ft <- fr[,,t][selected]
    ## update for selection (just loci under selection)
    fr[,,t][selected] <- (BB*ft^2 + BR*ft*(1-ft))/
      (BB*ft^2 + BR*2*ft*(1-ft) + RR*(1-ft)^2)
  }

  RES.F.fin <- fr[,,generations]
  parameters <- data.frame(populations,popsize,generations,
                           loci.per.s.value,beta1,beta2,p.neutral)
  results('parameters')
  results("F.init")
  results("RES.F.fin")
  results("s")

  labels <- data.frame(generations=generations,
                       loci=loci.per.s.value,
                       fixed='all',
                       populations=populations,
                       popsize=popsize,
                       p.neutral=p.neutral)

  pops <- data.frame(id=1:nrow(Type.POP),color=Type.POP)
  poplong <- reshape(pops,dir="long",varying=2:ncol(pops),timevar="population")
  names(poplong)[1] <- "locus"

  
  ## Ancestral allele frequency estimation plots
  ##myanim("naive-pi","Naive estimation of ancestral allele frequencies degrades for selected loci",anc.est.plot)
  ## Allele frequency variation animation
  ##myanim("fixanim","Allele frequency evolution varies with selection type and strength",freq.var.plot)

  ##   mypng('fixation-endpoints',
  ##         fixation.endpoints(convert.df(generations)),
  ##         paste(subdir,'/',sep=''),show='png')


  ##   grab <- function(type){
  ##     thistype <- s$type==type
  ##     maxs <- max(s$s[thistype])
  ##     cond <- s$s==maxs & thistype
  ##     c(which(cond & abs(F.init - 0.5) < 0.1)[1],
  ##       which(cond & abs(F.init - 0.5) > 0.3 & abs(F.init - 0.5) < 0.4)[1])
  ##   }
  ##   browser()
  ##   big <- combine(lapply(c(grab("none"),posloci <- grab("positive"),grab("balancing")),to.df))
  ## reshaping the data for single allele frequency over time plots
  ## to.df <- function(i){
  ##   x <- data.frame(fr[i,,])
  ##   names(x) <- paste('frequency',1:ncol(x),sep='.')
  ##   y <- reshape(x,dir='long',varying=names(x),idvar="population")
  ##   data.frame(y,locus=i)
  ## }
##   big <- something
##   z <- merge(merge(big,s),poplong)
##   ##print(unique(z[,c('s','type')]))
##   vals <- unique(z$locus)
##   locus.to.hilite <- posloci[2]
##   map <- rep(1:(length(vals)/3),3)
##   names(map) <- vals
##   #IDEA source("/home/thocking/directlabels/pkg/latticedirectlabels/R/direct.labels.R");dl(xyplot,z,frequency~time|locus+type,population,method="last.points",type="l",plot.symbol='.')
##   z$fakelocus <- map[as.character(z$locus)]
##   zna <- ddply(z,.(population,locus,type),
##                function(d){
##                  rbind(d[sort(d$time,i=TRUE)$i,],
##                        {tmp <- d[1,];tmp$frequency <- NA;tmp})})
  ##   xyplot(frequency~time|locus+type,zna,
  ##          groups=color,
  ##          par.settings=list(
  ##            superpose.line=list(col=c("blue","turquoise","red"))),
  ##          auto.key=list(points=F,lines=T,space="top",columns=3),
  ##          type="l")


  ##   mypng('allele-freq-evolution',
  ##   latticeplot(frequency~time|locus+type,zna,xyplot,
  ##               groups=color,
  ##               type='l',
  ##               main="Allele frequency variance over time can be used to distinguish different selection types",
  ##          par.settings=list(
  ##            superpose.line=list(col=c("blue","turquoise","red"))),
  ##          auto.key=list(points=F,lines=T,space="top",columns=3),
  ##               xlab="Generation",
  ##               alpha=1,
  ##               ylab="Blue allele frequency"),
  ##         paste(subdir,'/',sep=''),show='')

  bigplot <- function(g){
    line.xyplot <- function(x,y,...){
      z <- y[x==min(x)]
      panel.xyplot(x=x,y=y,...)
      panel.abline(v=g)
      ancest <- z[!is.na(z)][1]
      panel.abline(h=ancest,col="grey")
      ##lpoints(min(x),ancest,pch=2,col="black")
    }
    p <- xyplot(frequency~time,subset(zna,locus==locus.to.hilite),
                panel=line.xyplot,
                groups=color,
                type='l',
                ylim=c(0,1),
                main="Evolution of 1 locus over time in several populations",
                par.settings=list(
                  superpose.line=list(col=c("blue","turquoise","red"))),
                auto.key=list(points=F,lines=T,space="right",
                  title="Population type",cex.title=CEX.LEG),
                xlab="Generation",
                alpha=1,
                xlim=c(0,generations+1),
                ylab="Simulated blue allele frequency")
    plot(p,split=c(1,1,2,2))
    plot(anc.est.plot(g),split=c(2,1,2,2),newpage=F)
    p2 <- freq.var.plot(g)
    plot(p2,split=c(1,2,1,2),newpage=F)
    g
  }

  ##bigplot(100)
  ## animate both at the same time ??
  ##if(animate)
  ##  myanim("both","Allele frequency and ancestral estimate evolution",bigplot)
  molt <- melt(fr)
  names(molt) <- c("locus","population","generation","simulated")
  molt <- merge(molt,s)
  molt <- merge(molt,data.frame(ancestral=F.init,locus=1:Nlocus))
  molt <- merge(molt,poplong)
  list(id=timestr,
       s=s,
       freq=fr,
       fin=data.frame(molt,parameters),
       parameters=parameters)
}
