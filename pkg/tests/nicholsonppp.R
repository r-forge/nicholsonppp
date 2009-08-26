library(nicholsonppp)

## A really simple test data set:
observed <- matrix(rbeta(100,0.7,0.7)*100,nrow=20,ncol=5)
fit <- nicholsonppp(observed)
df <- cbind(melt(fit$a),melt(observed)$value)
names(df) <- c("locus","population","estimated","observed")
xyplot(estimated~observed,df)

## A data set that comes from an evolution simulation:
## FIXME: Too long of model fit.
## sim <- sim.drift.selection()
## ev.obs <- sim$sim[,,sim$p$gen]
## ev.fit <- nicholsonppp(ev.obs)
## head(round(ev.obs*100))
## head(round(ev.fit$a*10000))
## ev.df <- sim2df(sim)
## ev.df <- ev.df[ev.df$generation==sim$p$gen,]
## ev.df$estimated <- as.vector(ev.fit$a)
## xyplot(simulated~estimated,ev.df)
## sim$param
