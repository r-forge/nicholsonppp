library(nicholsonppp)

## A really simple test data set:
observed <- matrix(rbeta(100,0.7,0.7)*100,nrow=20,ncol=5)
fit <- nicholsonppp(observed)
df <- cbind(melt(fit$a),melt(observed)$value)
names(df) <- c("locus","population","estimated","observed")
xyplot(estimated~observed,df)

## A data set that comes from an evolution simulation:
sim <- sim.drift.selection(populations=4,loci=1)
ev.obs <- sim$sim[,,sim$p$gen]
ev.fit <- nicholsonppp(ev.obs)
## Compare alpha estimates to simulated values:
head(round(ev.obs*100,2))
head(round(ev.fit$a*100,2))
ev.df <- sim2df(sim)
ev.df <- ev.df[ev.df$generation==sim$p$gen,]
ev.df$estimated <- as.vector(ev.fit$a)
plot(xyplot(simulated~estimated,ev.df))
## Plot ancestral estimates given by model fit:
s <- anc.est.nicholson(sim,ev.fit)
plot(anc.est.plot(s))
