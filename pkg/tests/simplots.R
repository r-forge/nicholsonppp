library(nicholsonppp)
sim <- simulate.drift.selection(populations=12,
                                generations=50,
                                loci=5,
                                p.neutral=0.8,
                                s=c(1,5)/10)
loci.over.time(subset(sim$fin,locus<10),pop.colors=c('black','grey','white'))
loci.over.time(subset(sim$fin,locus==1))
anc.est.plot(subset(sim$fin,generation==generations),1)
fixation.endpoints(subset(sim$fin,generation==generations),hilite.locus=1)
bigplot(sim$fin)
evolution.animation("newtest","Allele frequency and ancestral estimate evolution",sim$fin)
