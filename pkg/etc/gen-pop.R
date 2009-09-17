source("../R/sim.R")
source("../R/model.R")
source("../R/plot.R")
source("../R/sensitivity.R")
library(ggplot2)
library(lattice)
library(latticedl)
dyn.load("../src/nicholsonppp.so")


pars <- expand.grid(generations=c(25,50,100),populations=c(4,8,12))
sims.gp <- mlply(pars,sim.drift.selection,loci.per.s.value=100)
models.gp <- llply(sims.gp,function(L)nicholsonppp(L$sim[,,L$p$gen]))
gp <- cbind(sims.gp,models.gp)

est.gp <- adply(gp,1,
                function(r)
                data.frame(anc.est.nicholson(r$sim,r$model),
                           populations=factor(r$sim$p$populations),
                           generations=factor(r$sim$p$generations)))
plist <- sims.gp[[1]]$p
plist <- plist[-which(names(plist)%in%names(pars))]
attr(est.gp,"parameters") <- plist
pdf("gen-pop.pdf")
anc.est.plot(est.gp,f=ancestral.est~ancestral|populations+generations,
             main="More generations and fewer populations increases estimate variability")
dev.off()

ppp.gp <- adply(gp,1,function(r)
                data.frame(r$sim$s,
                           ppp=r$model$ppp,
                           populations=factor(r$sim$p$populations),
                           generations=factor(r$sim$p$generations)))
pdf("gen-pop-dens.pdf",h=10,w=10)
dens.several.s(ppp.gp,f=~ppp|populations+generations,
               sub=deduce.param.label(plist),
               main="Number of populations affects PPP-values")
dev.off()

cl.gp <- classify.loci(ppp.gp,splitby=c("populations","generations"))
pdf("gen-pop-roc.pdf")
qplot(1-specificity,sensitivity,data=cl.gp,geom="line",
      group=.(populations,generations),
      colour=generations,
      size=populations)+theme_bw()+scale_colour_brewer(pal="Blues")+
  opts(title="PPP-value classifier depends on populations and generations")
dev.off()


pdf("anc-est-plot.pdf",h=10,w=10)
anc.est.plot(anc.est.naive(df[df$gen==sim$p$gen,]))
dev.off()

pdf("fixation-endpoints.pdf",h=10,w=10)
fixation.endpoints(df[df$gen==sim$p$gen,])
dev.off()

pdf("sim-summary-plot.pdf",h=10,w=12)
sim.summary.plot(df,g=80)
dev.off()
system("xpdf sim-summary-plot.pdf")

sim.summary.plot(sim2df(sims.gp[[1]]))
