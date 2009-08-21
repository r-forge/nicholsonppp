system("R CMD SHLIB test.f90")
try(dyn.unload("test.so"))
dyn.load("test.so")
.Fortran("update_alpha",as.integer(c(9,100)),2,foo="bar")
