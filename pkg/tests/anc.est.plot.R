sim <- sim.drift.selection()
df <- sim2df(sim)
anc.est.plot(subset(df,generation==generations))


