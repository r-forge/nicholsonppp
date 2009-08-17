simulate.drift.selection <- function
### Simulate allele frequency evolution in several populations
### according to a simple drift and selection model.
(populations=12,
### Number of populations.
 popsize=1000,
### Size of each simulated population.
 generations=200,
### Generations of evolution to simulate.
 loci.per.s.value=50,
### Number of loci simulated for each s value.
 beta1=0.7,
### Parameter for beta distribution of initial allele frequencies.
 beta2=0.7,
### Parameter for beta distribution of initial allele frequencies.
 p.neutral=0.8,
### Proportion of neutral alleles to simulate.
 s=c(1,5)/100,
### Vector of selection strength values to simulate.
 adapt.pop=c(blue=0.4,red=0.4,neutral=0.2)
### Probabilities of color adaptation in populations.
 ){
  timestr <- gsub('[ :]','-',format(Sys.time()))
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
  F.init <- rbeta(Nlocus,beta1,beta2)
  cutoff <- 0.05
  F.init[F.init<cutoff] <- cutoff
  F.init[F.init>(1-cutoff)] <- 1-cutoff
  
  s <- data.frame(s=s.LOC,type=Type.LOC,locus=1:length(s.LOC),ancestral=F.init)
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
  cat("Generation: ")
  for(t in 2:generations){
    cat(" ",t,sep="")
    ## update all loci for drift
    fr[,,t] <- rbinom(populations*Nlocus,popsize,fr[,,t-1])/popsize
    ft <- fr[,,t][selected]
    ## update for selection (just loci under selection)
    fr[,,t][selected] <- (BB*ft^2 + BR*ft*(1-ft))/
      (BB*ft^2 + BR*2*ft*(1-ft) + RR*(1-ft)^2)
  }
  cat("\n")

  cat("Converting simulated allele frequencies to molten data\n")
  molt <- melt(fr)
  names(molt) <- c("locus","population","generation","simulated")
  
  cat("Merging selection coefficients\n")
  molt <- merge(molt,s)

  cat("Merging population colors\n")
  pops <- data.frame(id=1:nrow(Type.POP),color=Type.POP)
  poplong <- reshape(pops,dir="long",varying=2:ncol(pops),timevar="population")
  names(poplong)[1] <- "locus"
  molt <- merge(molt,poplong)

  cat("Adding simulation parameters\n")
  parameters <- data.frame(populations,popsize,generations,
                           loci.per.s.value,beta1,beta2,p.neutral,id=timestr)
  data.frame(molt,parameters)
### Data frame containing all simulated data and parameter values, 1
### line per (population,locus,generation) tuple.
}
