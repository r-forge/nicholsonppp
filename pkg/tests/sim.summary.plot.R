library(nicholsonppp)
sim <- simulate.drift.selection()
df <- sim2df(sim)
sim.summary.plot(df)
