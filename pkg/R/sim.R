## Question: doesn't this result in doing the migration calculation
## several times for each population? That is, when 1 population gains
## members, another population has to lose members? Is that taken into
## account with the calcul?
update.freq.mig<-function(Freq.all1,F.pop,Taille.POP){ 
  ##Freq.all1=vecteur des frequences alleliques dans chaque population
  ##F.pop=vecteur des F pour chaque pop
  Freq.new=Freq.all1
  Npop=length(Freq.all1)
  m.pop=(1-F.pop)/(2*F.pop*Taille.POP)
  m.pop[m.pop>1]=1
  for(i in 1:Npop){
    ##migration
    n.imm=rbinom(1,Taille.POP,m.pop[i])
    if(n.imm>0){
      n1.lost=rbinom(1,n.imm,Freq.new[i])
      imm.orig=rmultinom(1,n.imm,rep(1/(Npop-1),Npop-1))
      n1.gain=0
      for(j in 1:(Npop-1)){
        n1.gain=n1.gain+rbinom(1,imm.orig[j],Freq.new[-i][j])
      }
      Freq.new[i]=Freq.all1[i]+(n1.gain-n1.lost)/Taille.POP
    }
    Freq.new[i]=min(1,max(0,Freq.new[i])) 
  }
  list(Freq.new=Freq.new)
}
##debug(update.freq.mig)

## apply-ified
update.freq.mig2 <- function(Freq.all1,F.pop,Taille.POP){ 
  ##Freq.all1=vecteur des frequences alleliques dans chaque population
  ##F.pop=vecteur des F pour chaque pop
  Freq.new <- Freq.all1
  Npop <- length(Freq.all1)
  m.pop <- (1-F.pop)/(2*F.pop*Taille.POP)
  m.pop[m.pop>1] <- 1
  ## number of immigrants in each population ??
  n.immi <- rbinom(Npop,Taille.POP,m.pop)
  lost <- rbinom(Npop,n.immi,Freq.new)
  imm.origi <- sapply(n.immi,function(ni)rmultinom(1,ni,rep(1/(Npop-1),Npop-1)))
  subfreqs <- sapply(1:Npop,function(i)Freq.new[-i])
  n1.gaini <- matrix(rbinom(length(subfreqs),imm.origi,subfreqs),
                     nrow=nrow(subfreqs),ncol=ncol(subfreqs))
  gain <- colSums(n1.gaini)
  fnew <- Freq.all1 + (gain-lost)/Taille.POP
  sapply(fnew,function(x)min(1,max(0,x)))
}
##debug(update.freq.mig2)

update.freq.sel<-function(Freq.all1,w.coef,type.pop){ 
  ##Freq.all1=vecteur des frequences alleliques dans chaque population
  ##w.coef=coefficient de seleciton pour chaque genotype
  ##type.pop=1 (bleu),2 (rouge) ou 3 (blanc)
  Freq.new=Freq.all1
  Npop=length(Freq.all1)
  for(i in 1:Npop){   
    if(type.pop[i]!=3){
      tmp.w=w.coef 
      tmp.freq=Freq.new[i]
      tmp.freq=c(tmp.freq**2,2*tmp.freq*(1-tmp.freq),(1-tmp.freq)**2)
      if(type.pop[i]==2){tmp.w=rev(w.coef)}
      w.bar=t(tmp.freq)%*%tmp.w
      w.num=tmp.freq[1]*w.coef[1]+(tmp.freq[2]*w.coef[2])/2
      Freq.new[i]=min(1,max(0,w.num/w.bar))     
    }
  }
  Freq.new
}
##debug(update.freq.sel)

FST.WC<-function(DATA.SIZE,DATA.FREQ){
  Nrace=dim(DATA.SIZE)[2]
  SumNi=rowSums(DATA.SIZE)
  Nic=DATA.SIZE-(DATA.SIZE**2)/SumNi
  Nc=rowSums(Nic)/(Nrace-1)

  MSG=(rowSums(DATA.FREQ*(1-DATA.FREQ)*DATA.SIZE)) /(SumNi-1)
  PA=rowSums(DATA.FREQ*DATA.SIZE)/SumNi
  MSP=(rowSums(DATA.SIZE*((DATA.FREQ-PA)**2)))/(Nrace-1)

  FST=(MSP-MSG)/(MSP+(Nc-1)*MSG)
  FST[FST<0]=0 ##mise � zero des FST negatifs

  ##calcul pour des Bi estime FSTi pour cahque pop relativement a ThetaA: Equation 8 Weir et Hill 2002
  ##Attention cas bi-allelique
  SumNic=rowSums(Nic)
  Pmoy.u=rowSums(DATA.SIZE*DATA.FREQ)/SumNi
  Num=SumNic*DATA.SIZE*DATA.FREQ*(1-DATA.FREQ)/(DATA.SIZE-1)
  Den=rowSums(DATA.SIZE*((DATA.FREQ-Pmoy.u)**2) + Nic*DATA.FREQ*(1-DATA.FREQ))
  BETA.I=1-(Num/Den)

  list(FST=FST,FSTmoy=mean(FST,na.rm=T),FST.I=BETA.I)
}


sim <- function(Npop=8,Taille.POP=1000,Ngeneration=25,N.s.rep=50,
                migrate=F,beta1=0.7,beta2=0.7,p.neutral=1/3,s,scale=1,shape=1){
  timestr <- gsub('[ :]','-',format(Sys.time()))
  subdir <- paste('results',timestr,sep='/')
  rfile <- function(name)paste(subdir,'/',name,'.txt',sep='')
  rpdf <- function(name)
    pdf(paste(subdir,'/',timestr,'-',name,'.pdf',sep=''),h=10,w=7)
  results <- function(data,name='')
    write.table(data,rfile(name),quote=F,row.names=T,col.names=T)
  dir.create(subdir)

  srep <- rep(s,each=N.s.rep)
  Nloc <- length(srep)
  n.neutral <- 2*Nloc*p.neutral/(1-p.neutral)
  s.LOC=c(rep(0,n.neutral),srep,srep)
  Type.LOC=factor(c(rep("NEU",n.neutral),rep("POS",Nloc),rep("BAL",Nloc)))
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
  
  NOM.LOC=paste(seq(1:length(Type.LOC)),Type.LOC,sep="")
  names(s.LOC)=NOM.LOC
  names(Type.LOC)=NOM.LOC
  d <- data.frame(f=F.init,s=s.LOC,type=Type.LOC,id=1:length(s.LOC))
  usual.colors <- trellis.par.get("superpose.symbol")$col
  usual.colors[3] <- "green"

  ## probabilities of color adaptation in populations
  adapt.pop <- c(blue=0.4,red=0.4,neutral=0.2)
  ## each population/locus combination has a different "color"
  Type.POP <- matrix(sample(names(adapt.pop),Npop*Nlocus,T,adapt.pop),
                     Nlocus,Npop)

  ## this will store all of the simulated allele frequencies
  fr <- array(F.init,c(Nlocus,Npop,Ngeneration))
  ## these will be intermediate data structures to support fast selection
  selected.loci <- matrix(as.character(Type.LOC)!="NEU",Nlocus,Npop)
  selective.pops <- Type.POP!='neutral'
  selected <- selected.loci & selective.pops
  s <- matrix(s.LOC,Nlocus,Npop)[selected]
  type <- matrix(Type.LOC,Nlocus,Npop)[selected]
  popcol <- Type.POP[selected]
  ## Genotype weights for each locus and population under selection:
  ## (not including neutral loci)
  BB <- ifelse(type=="POS"&popcol=="blue",1+s  ,1)
  BR <- ifelse(type=="POS",               1+s/2,1+s)
  RR <- ifelse(type=="POS"&popcol=="red", 1+s  ,1)
  ## diagnostic data frame
  ## data.frame(BB,BR,RR,popcol,s,type)
  for(t in 2:Ngeneration){
    cat("Generation: ",t,"\n")
    ## update all loci for drift
    fr[,,t] <- rbinom(Npop*Nlocus,Taille.POP,fr[,,t-1])/Taille.POP
    ft <- fr[,,t][selected]
    ## update for selection (just loci under selection)
    fr[,,t][selected] <- (BB*ft^2 + BR*ft*(1-ft))/(BB*ft^2 + BR*2*ft*(1-ft) + RR*(1-ft)^2)
  }

  RES.F.fin <- fr[,,Ngeneration]

  rownames(RES.F.fin) <- NOM.LOC

  results(data.frame(Npop,Taille.POP,Ngeneration,migrate,N.s.rep,beta1,beta2,p.neutral),
          'parameters')
  results(F.init,"F.init")
  results(RES.F.fin,"RES.F.fin")
  results(s.LOC,"S.LOC")

  labels <- data.frame(generations=Ngeneration,
                       loci=N.s.rep,
                       fixed='all',
                       populations=Npop,
                       popsize=Taille.POP,
                       p.neutral=p.neutral)

  pops <- data.frame(id=1:nrow(Type.POP),color=Type.POP)
  poplong <- reshape(pops,dir="long",varying=2:ncol(pops),timevar="population")
  names(poplong)[1] <- "locus"
  ## conversion to format compatible with endpoint allele frequency
  ## fixation plots
  convert.df <- function(g){
    x <- data.frame(fr[,,g],fr[,1,1])
    names(x) <- paste('alpha',1:ncol(x),sep='.')
    y <- reshape(x,dir='long',varying=names(x),timevar='population')
    names(y)[3] <- "locus"
    init <- y$pop==max(y$pop)
    z <- merge(data.frame(y[!init,],pi=y[init,'alpha']),d)
    xx <- merge(poplong,z)
    final <- data.frame(xx,labels)
    final$locus <- factor(final$locus)
    ## with lattice the easiest way to do legends is just with auto.key
    ## so in order to have ancestral allele frequencies lets add some
    ## rows to the data frame
    ancest <- final
    ancest$color <- 'ancestral'
    ancest$alpha <- ancest$pi
    ancest$population <- 0
    ancest <- unique(ancest)
    rbind(ancest,final)
  }
  ##debug(convert.df)

  myanim <- function(f,tit,pfun){
    naivedir <- paste(getwd(),subdir,f,sep="/")
    dir.create(naivedir)
    ani.start(nmax=Ngeneration,
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

  ## reshaping the data for single allele frequency over time plots
  to.df <- function(i){
    x <- data.frame(fr[i,,])
    names(x) <- paste('frequency',1:ncol(x),sep='.')
    y <- reshape(x,dir='long',varying=names(x),idvar="population")
    data.frame(y,locus=i)
  }
  ##debug(to.df)
  
  anc.est.plot <- function(g){
    est.df <- data.frame(estimate=rowMeans(fr[,,g]),d)
    pp <- est.df[locus.to.hilite,]
    pi.xyplot(estimate~f,est.df,alpha=1,
              cols=usual.colors,
              panel=function(...){
                panel.xyplot(...)
                lpoints(pp$f,pp$est,pch=1,cex=2,col="black")
              },
              xlab="Simulated blue allele frequency",
              ylab="Estimated blue allele frequency",
              main="Allele frequency estimates vary with selection type",
              auto.key=list(space="right",title="Selection type",cex.title=CEX.LEG))
  }
  
  freq.var.plot <- function(g){
    fixation.endpoints(convert.df(g),
                       hilite.locus=locus.to.hilite,
                       main="Allele frequency variance depends on selection type and strength",
                       par.settings=list(fontsize=list(text=12,points=8)))
  }

  
  ## Ancestral allele frequency estimation plots
  ##myanim("naive-pi","Naive estimation of ancestral allele frequencies degrades for selected loci",anc.est.plot)
  ## Allele frequency variation animation
  ##myanim("fixanim","Allele frequency evolution varies with selection type and strength",freq.var.plot)

  ##   mypng('fixation-endpoints',
  ##         fixation.endpoints(convert.df(Ngeneration)),
  ##         paste(subdir,'/',sep=''),show='png')


  grab <- function(type){
    thistype <- d$type==type
    maxs <- max(d$s[thistype])
    c(which(d$s==maxs & thistype & abs(d$f - 0.5) < 0.1)[1],
      which(d$s==maxs & thistype & abs(d$f - 0.5) > 0.3 & abs(d$f - 0.5) < 0.4)[1])
  }
  big <- combine(lapply(c(grab("NEU"),posloci <- grab("POS"),grab("BAL")),to.df))
  names(d)[names(d)=="id"] <- "locus"
  z <- merge(merge(big,d),poplong)
  ##print(unique(z[,c('s','type')]))
  vals <- unique(z$locus)
  locus.to.hilite <- posloci[2]
  map <- rep(1:(length(vals)/3),3)
  names(map) <- vals
  #IDEA source("/home/thocking/directlabels/pkg/latticedirectlabels/R/direct.labels.R");dl(xyplot,z,frequency~time|locus+type,population,method="last.points",type="l",plot.symbol='.')
  z$fakelocus <- map[as.character(z$locus)]
  zna <- ddply(z,.(population,locus,type),function(d){rbind(d[sort(d$time,i=T)$i,],{tmp <- d[1,];tmp$frequency <- NA;tmp})})
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
                xlim=c(0,Ngeneration+1),
                ylab="Simulated blue allele frequency")
    plot(p,split=c(1,1,2,2))
    plot(anc.est.plot(g),split=c(2,1,2,2),newpage=F)
    p2 <- freq.var.plot(g)
    plot(p2,split=c(1,2,1,2),newpage=F)
    g
  }

  ##bigplot(100)
  ## animate both at the same time ??
  myanim("both","Allele frequency and ancestral estimate evolution",bigplot)
  
  timestr # return the id of the simulation just created
}
##debug(sim)


