## prereqs to running this script:
## params, results


source("db.R")
subdir <- paste("results/",ids,"/",sep="") ## assume 1 sim
x <- combine(lapply(names(results),make.lattice.table))
x$p.neutral <- factor(round(x$p.neutral,digits=2))
x$error <- x$actual-x$estimated
mypng2 <- function(pre,P)mypng(pre,P,subdir,show='png')
##Comparing loci exclusion rules
mypng2('pi-neutral-fixed',
      pi.xyplot(estimated~actual|fixed,
                x,
                aspect='iso',
                main="Ancestral allele frequency is harder to estimate for fixed loci")
      )
z <- x[as.character(x$fixed)=="all",]
mypng2('pi-error-fixed',
pi.xyplot(error^2~log(s)|fixed+p.neutral,
          z[as.character(z$type)!="NEU",],
          aspect='fill',
          panel=panel.xyplot,
          main='Error in estimating ancestral allele frequency does not depend on proportion of neutral alleles')
      )



##debug(my.xyplot)
## loci independent of correlation value
## this is expected since each loci is simulated independent of each other
##pdf(paste(format(Sys.time(),"%F"),"loci-indep.pdf",sep='-'),h=10,w=10)
mypng2('pi',
      pi.xyplot(estimated~actual,
                z,
                aspect='iso',
                main="Ancestral allele frequency estimate depends on selection type")
      )
mypng2('ppp',
      latticeplot(~ppp,z,densityplot,
            group=type,
                  alpha=1,
            main="PPP values depend on selection type",
            auto.key=list(space='top',columns=3))
      )
## ppp-value dependent on s-value?
y <- z[as.character(z$type)!="NEU",]
mypng2('ppp-s',
      pi.xyplot(ppp~log(s),y,aspect='fill',
                panel=panel.xyplot,
                main='PPP values depend on selection strength')
      )

