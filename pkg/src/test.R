system("R CMD SHLIB test.f90");try(dyn.unload("test.so"));dyn.load("test.so");.Fortran("update_alpha",Y=as.integer(c(9,100)),n=as.integer(2),foo="bar")
