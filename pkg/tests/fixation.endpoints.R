sim <- sim.drift.selection()
df <- sim2df(sim)
fixation.endpoints(subset(df,generation==generations))
