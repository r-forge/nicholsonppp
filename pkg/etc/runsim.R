library(plyr)
source("db.R")
source("sim.R")
s1 <- seq(0.0025,0.01,0.0025)
s2 <- seq(0.02,0.05,0.01)
s <- c(s1,s2,0.075,0.1,0.5,0.75,1,1.5,2)
length(s) # 15
s <- c(0.01,0.05)
ids <- c()
ids <- sim(migrate=F,Npop=12,Ngeneration=200,p.neutral=0.8,s=s)


## ids is now a vector of simulation ids
result.dirs <- paste("results/",ids,"/",sep='')
names(result.dirs) <- ids
source("db.R") # for read.res.tables
source("nicholson.run.R")
result.list <- lapply(result.dirs,read.res.tables,c(alpha.sim="RES.F.fin.txt"))
fit.nicholson(ids)

## read results of model fit and make fixation plots
results <- read.results(result.dirs)
pfiles <- paste(result.dirs,"parameters.txt",sep='/')
params <- combine(lapply(pfiles,read.params))
rownames(params) <- ids
sapply(ids,fixation.plot)

source("nicholson.analyze.R")
source("naive.R")
