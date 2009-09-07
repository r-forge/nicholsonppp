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
density.several.s <- function
### Plot density curves for ppp values for each selection type, to
### visualize the extent to which PPP-values can be used to indicate
### selection state in a simulation.
(df,
### Data frame with columns ppp type s to be plotted.
 ...
### To be passed to dl.
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
