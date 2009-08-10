library(R2WinBUGS)
ids <- "2009-06-30-15-56-18"
result.dirs <- paste("results/",ids,"/",sep='')
names(result.dirs) <- ids
source("db.R") # for read.res.tables
r <- read.res.tables(result.dirs,c(Y="Y_SIM_HMN",N="N_SIM_HMN"))
r$I <- nrow(r$Y)
r$J <- ncol(r$Y)
r$N <- as.matrix(r$N)
r$Y <- as.matrix(r$Y)
pfiles <- paste(result.dirs,"parameters.txt",sep='/')
params <- combine(lapply(pfiles,read.params))
rownames(params) <- ids
## Variables in model:
## data: I J Y N
## init: alpha p tau c
init <- function(){
  list(p=
       c=)
}
inits <- list(list(p=rep(0.5,r$I),c=rep(0.1,r$J))) #,alpha=as.matrix(Y/N))
vars <- c("alpha","p","c")
res <- bugs(r,inits,vars,"/home/thocking/projects/bugs/nicholson.bugs",n.chains=1,n.iter=500,n.burnin=100)
