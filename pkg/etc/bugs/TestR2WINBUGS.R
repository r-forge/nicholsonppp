library(R2WinBUGS)
# Copied from Prof Andrew Gelman's example
model.file <- system.file(package = "R2WinBUGS", "model", "schools.txt")
file.show(model.file)
data(schools)
J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list ("J", "y", "sigma.y")
inits <- function(){
    list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100),
        sigma.theta = runif(1, 0, 100))
}
parameters <- c("theta", "mu.theta", "sigma.theta")

schools.sim<-
  bugs(data, inits, parameters, model.file,
       n.chains = 3,
       n.iter = 5,debug=T,
       bugs.directory = "/home/thocking/.wine/drive_c/Program\ Files/WinBUGS14/",
       working.directory = "/tmp", clearWD = FALSE, useWINE=T, newWINE=T,
       WINE="/usr/bin/wine",WINEPATH="/usr/bin/winepath")

################################################################

library(fUtilities)
setwd("/home/mathieu/Desktop/ANALYSES/TEST_NICHOLSON/TEST_R")
DATA_Y=read.table("R_Yfile.txt",sep="\t",header=T,row.names=1)
DATA_N=read.table("R_Nfile.txt",sep="\t",header=T,row.names=1)
LISTE_RACE=read.table("LISTE_RACE_FDIST")[,1]
Nmrk=dim(DATA_Y)[1];Npop=dim(DATA_Y)[2]

modele<-"/home/mathieu/Desktop/ANALYSES/TEST_NICHOLSON/TEST_R/Test_model.txt"
file.show(modele)
Y<-as.matrix(read.table("R_Yfile.txt",sep="\t",skip=1)[,-1])
N<-as.matrix(read.table("R_Nfile.txt",sep="\t",skip=1)[,-1])
I<-Nmrk;J<-Npop
data <- list (I=I, J=J, Y=Y,N=N)
inits <- list(list(p=rep(0.5,Nmrk),c=rep(0.1,Npop))) #,alpha=as.matrix(Y/N))
parameters <- c("p","c","alpha")
WINEPATH="/usr/bin/winepath"

tst.sim <- bugs(data, inits, parameters, modele, n.chains = 1,
n.iter = 500,n.burnin=100,bugs.directory = "/home/mathieu/.wine/drive_c/Program\ Files/WinBUGS14/",
working.directory = "/tmp", clearWD = FALSE, useWINE=T, newWINE=T,
WINE="/usr/bin/wine",WINEPATH=WINEPATH)

#####Test extractions aleatoirs de N markers
Nechantillons=300
Ntirages=500
RES_TIRAGES=list()
MRK_TIRES=list()
inits <- list(list(p=rep(0.5,Nechantillons),c=rep(0.1,Npop))) 
parameters <- c("p","c","alpha")
WINEPATH="/usr/bin/winepath"
for(i in 1:Ntirages){
  tmp_mrk_sel=round(runif(Nechantillons,1,Nmrk))
  MRK_TIRES[[i]]=tmp_mrk_sel
  tmp_Y=Y[tmp_mrk_sel,];tmp_N=N[tmp_mrk_sel,];tmp_I=Nechantillons;tmp_J=Npop
  data <- list(I=tmp_I, J=tmp_J, Y=tmp_Y,N=tmp_N)
RES_TIRAGES[[i]]<-bugs(data, inits, parameters, modele, n.chains = 1,n.iter = 2000,n.burnin=1000, #n.thin=1,
                       bugs.directory = "/home/mathieu/.wine/drive_c/Program\ Files/WinBUGS14/",
                       working.directory = "/tmp", clearWD = FALSE, useWINE=T, newWINE=T,WINE="/usr/bin/wine",WINEPATH=WINEPATH)  
cat(paste("TIRAGE ",i," REALISE\n",sep=""))
}
