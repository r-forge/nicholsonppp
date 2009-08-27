library(R2WinBUGS)
ids <- "2009-06-30-15-56-18"
result.dirs <- paste("~/projects/sim/results/",ids,"/",sep='')
names(result.dirs) <- ids
source("~/projects/sim/db.R") # for read.res.tables
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

betafit <- function(beta_pi=0.7){
  bugsfile <- paste("~/nicholsonppp/pkg/etc/bugs/nicholson.bugs.",
                    beta_pi,sep="")
  res <- bugs(r,inits,vars,bugsfile,n.chains=1,
              n.iter=6000,n.burnin=5000,n.thin=50)
  res$mean$p
}

fit0.7 <- nicholsonppp(r$Y,beta_pi=0.7)
fit1 <- nicholsonppp(r$Y,beta_pi=1)
winbugs0.7 <- betafit(0.7)
winbugs1 <- betafit(1)
R <- read.res.tables(result.dirs)
est <- data.frame(simulated=R$pi.sim[,1],
                  fortran.old=R$pi.est$MOY,
                  fortranr0.7=fit0.7$p,
                  fortranr1=fit1$p,
                  winbugs0.7,
                  winbugs1,
                  locus=factor(1:r$I),
                  ppp=R$ppp$PPPVALtot[1:r$I],
                  R$s)
pdf("2009-08-27-notbeta.pdf",paper="a4",h=0,w=0)
splom(~est[1:6],est,groups=type,auto.key=T,main="Beta parameters do not explain discrepancy between estimates",sub="Similar programs yield similar estimates")
dev.off()


library(ggplot2)
molten <- melt(est,c(1,4),2:3)
names(molten)[3:4] <- c("estimator","estimate")
##xyplot(estimate~simulated,molten,aspect=1,groups=estimator,auto.key=TRUE)
##xyplot(fortran~winbugs,est,aspect=1,groups=type,auto.key=TRUE,panel=function(...){panel.abline(0,1);panel.xyplot(...)})

library(latticedl)
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
