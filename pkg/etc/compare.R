## We did a bunch of simulations, now we need to compare them to see
## what happened when we change 1 parameter and leave the rest fixed.

## this file contains functions primarily for reading the database

## to be executed in projects/sim

## Database subset definition -- just read a subset of sims
gstrs <- c("2009-05-14","2009-05-26")

source("db.R")
pfiles <- dir('results',pattern='parameters.txt',recursive=T,full.names=T)
pfiles <- pfiles[unique(unlist(lapply(gstrs,function(gs)grep(gs,pfiles))))]
names(pfiles) <- sapply(strsplit(pfiles,split='/'),function(x)x[2])
result.dirs <- sub('parameters.txt','',pfiles)
params <- combine(lapply(pfiles,read.params))
rownames(params) <- names(pfiles)

##uncomment for subset based on parameter values
ss <- rownames(params[params$N.s.rep==50,])
result.dirs <- result.dirs[ss]
params <- params[ss,]


results <- read.results(result.dirs)


class(results) <- c('simresult','list')
source("nicholson.run.R")
