library(nicholsonppp)
sim <- sim.drift.selection()
df <- sim2df(sim)
sim.summary.plot(df,ceiling(sim$p$gen/3*2))
