library(nicholsonppp)
sim <- sim.drift.selection()
df <- sim2df(sim)
anc.est.plot(df[df$generation==sim$p$gen,])
