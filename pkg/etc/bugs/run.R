library(R2WinBUGS)
ss<- bugs(list("J", "y", "sigma.y"), inits, c("theta", "mu.theta", "sigma.theta"), model.file)
