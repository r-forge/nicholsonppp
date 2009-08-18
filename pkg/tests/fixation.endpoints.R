library(nicholsonppp)
sim <- sim.drift.selection()
df <- sim2df(sim)
fixation.endpoints(df[df$generation==sim$p$gen,])
