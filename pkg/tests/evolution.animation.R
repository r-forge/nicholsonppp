library(nicholsonppp)
sim <- sim.drift.selection()
df <- sim2df(sim)
evolution.animation(df[df$generation%%10==1,])
