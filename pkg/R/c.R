nppp.over.time <- function
### Fit a series of models over time.
(sim,
### Simulation to fit models on.
 generation=seq(0,sim$p$gen,l=5)[-1]
### Vector of generation numbers to fit on. By default, fit 4 equally
### spaced models.
 ){
  r <- mlply(data.frame(generation),function(g)nicholsonppp(sim$sim[,,g]))
  attr(r,"popsize") <- sim$p$popsize
  r
### List of model fit result lists, with attributes from mlply.
}

estc.df <- function
### Convert list of model fit lists from nppp.over.time to a data
### frame suitable for plotting c values.
(fit.list
### List of model fit lists from nppp.over.time.
 ){
  cc <- melt(ldply(fit.list,function(L)L$c),id="generation")
  names(cc) <- c("generation","population","c.est")
  g <- attr(fit.list,"split_labels")$g
  cc$population <- factor(cc$population)
  cc <- cc[order(cc$generation),]
  cc$popsize <- attr(fit.list,"popsize")
  cc$type <- factor("simulated")
  add.level <- function(popsize){
    denom <- popsize
    rbind(cc,data.frame(population=factor(paste("g/",denom,sep="")),
                        generation=g,c.est=g/denom,popsize=popsize,
                        type=factor("theoretical")))
  }
  for(ps in unique(cc$popsize))cc <- add.level(ps)
  cc
### Data frame for plotting the simulated and theoretical values of c.
}

estc.over.time <- function
### Plot evolution of C estimates over time. This is the
### differentiation parameter and is expected to increase linearly
### over time.
(data,
### Data frame to plot, result of estc.df.
 ...
### Other arguments to qplot.
 ){
  qplot(generation,c.est,data=data,group=population,
        size=factor(popsize),geom="line",colour=type,...)
### The ggplot2 plot.
}
