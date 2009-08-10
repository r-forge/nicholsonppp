source("db.R")
com <- combine(lapply(names(pfiles),compare.alleles))
params2 <- params
names(params2) <- rename(names(params2),
                         c(N.s.rep='loci',
                           Ngeneration='generations',
                           Taille.POP='popsize',
                           Npop='populations'))
params2$id <- rownames(params2)
com <- merge(com,params2,'id')
nonzero <- !is.nan(com$prop)
nonall <- !as.character(com$fixed)=='all'
prange <- range(com$prop[!(is.nan(com$prop)|com$prop%in%c(0,1))])
##pdf(pdfname <- paste(format(Sys.time(),"%F"),"selection-fixation.pdf",sep='-'),h=12,w=10)
mypng('fixation-selection',
latticeplot(prop~S|generations+populations,
            com[nonzero&nonall&com$loci==50&com$popsize==1000,],
            xyplot,
            alpha=1,
            groups=factor(fixed),
            ylab='proportion of loci "not fixed"',
            xlab='selection type and coefficient',
            main='Loci fixation depends on selection type and coefficient',
            scales=list(
              x=list(rot=90,cex=0.65,fontfamily='mono'),
              y=list(axs='i',rot=0,at=c(0,1,0.25,0.5,0.75))),
            auto.key=list(space='right'),
            panel=function(...){grid.lines(c(0.5,0.5),c(0,1),gp=gpar(col='grey'));panel.xyplot(...)},
            ylim=c(0,1))
      )
##dev.off()


results <- lapply(result.dirs,read.res.tables,c(s="S.LOC.txt"))
repvals <- unique(params$N.s.rep)
names(repvals) <- repvals
rep.description <- sapply(repvals,describe.alleles,simplify=F)
