library(R2WinBUGS)
ids <- "2009-06-30-15-56-18"
result.dirs <- paste("~/projects/sim/results/",ids,"/",sep='')
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
inits <- list(list(p=rep(0.5,r$I),c=rep(0.1,r$J))) #,alpha=as.matrix(Y/N))
vars <- c("alpha","p","c")
res <- bugs(r,inits,vars,"/home/thocking/projects/bugs/nicholson.bugs",n.chains=1,n.iter=500,n.burnin=100)
R <- read.res.tables(result.dirs)
est <- data.frame(simulated=R$pi.sim[,1],
                  fortran=R$pi.est$MOY,
                  winbugs=res$mean$p,
                  locus=factor(1:r$I),
                  ppp=R$ppp$PPPVALtot[1:r$I],
                  R$s)
molten <- melt(est,c(1,4),2:3)
names(molten)[3:4] <- c("estimator","estimate")
xyplot(estimate~simulated,molten,aspect=1,groups=estimator,auto.key=TRUE)
xyplot(fortran~winbugs,est,aspect=1,groups=type,auto.key=TRUE,panel=function(...){panel.abline(0,1);panel.xyplot(...)})

pdf("winbugs-fortran-disagree.pdf",h=8,w=8)
dl(xyplot,est,fortran~winbugs,type,aspect=1,
   panel=function(...){panel.abline(0,1);panel.xyplot(...)},
   method="empty.grid",
   main="2 programs yield different ancestral allele frequency estimates for the same model")
dev.off()

d <- est$fortran-est$winbugs
den <- density(d)
sorted <- sort(d)
cutoff <- sorted[which.max(diff(sorted))]
cutoff <- -min(d)
densityplot(~fortran-winbugs,est,
            aspect=1,groups=type,auto.key=TRUE,
            panel=function(...){
              panel.densityplot(...)
              panel.abline(v=cutoff)
  })
histogram(~fortran-winbugs,est,n=100)
est
