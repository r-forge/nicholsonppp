system("R CMD SHLIB *.f90")
try(dyn.unload("0mers_twist.so"))
dyn.load("0mers_twist.so")
observed <- matrix(rbeta(100,0.7,0.7)*100,nrow=20,ncol=5)
source("../R/model.R")
fit <- nicholsonppp(observed)

