system("R CMD SHLIB fitnicholsonppp.f90")
dyn.load("fitnicholsonppp.so")
source("../R/model.R")
nicholsonppp(matrix(rbeta(100,0.7,0.7)*100,nrow=20,ncol=5))

