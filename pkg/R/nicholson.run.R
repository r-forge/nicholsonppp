## takes result.list as input, fits nicholson model

## STRUCTURE OF NICHOLSON INPUT FILE (to create)
##    -nicholson.inp:
##       -ligne1 (4elements):
##             -Nbre_de_Marquers
##             -Nbre_Populations
##             -FILENAME_N_allelesREFERENCES_Observes FILENAME_N_allelesTOT_Observes
##             -FILENAME_N_allelesTOT_Observes
##       -ligne2 (5elements):
##             -nvaleurs à generer (ex=5000)
##             -thin rate (ex=25)
##             -burn_in= (ex 5000)
##             -Netudes pilots maximum (ex=100)
##             -pilot_length (ex=1000)
##       -ligne3 (7 elements):
##             -delta_a (ex=0.6), delta_pi (ex=0.5), delta_c (ex=0.1):  valeur pour la marche aleaoitre   
##             -ajustement du taux de changement lors des etudes pilotes (ex=1.25)
##             -limite inferieure pour le taux d'acceptation (ex=0.25)
##             -limite superieure pour le taux d'acceptation (ex=0.40)
##             -graine pour le mersenne twister (ex: 4357)
##       -ligne4 (1 element):
##             -si 0: on fait les summary stat pour les alphas: ATTENTION necessite d'allouer une enorme matrice Nmrk*Npop*Nvaleurs 
##             -si 1: on calcule les summary stats que la moyenne et la variance
## 	    -si 2: pour les enormes jeux de donnees: on sort un fichier RES.out (tres volumineux) qui contient ITER POP MRK Alpha_ij FST_ij Pi_i Cj (il pourra ensuite etre traite avec un autre porgramme en dynamique)  
## EXAMPLE
## 55961 3 Y_SIM_HMN N_SIM_HMN
## 1000 50 5000 50 1000
## 0.4 1.25 0.06 1.25 0.20 0.45 4501
## 1

write.nicholson <- function(d,f)
  write.table(round(d),f,col.names=F,row.names=F,quote=F)
run.nicholson <- function(d,N,f,trt.fun.name=NULL){
  yfile <- "Y_SIM_HMN"
  nfile <- "N_SIM_HMN"
  sfile <- "S.LOC.txt"
  ## if treatment specified we create a subdirectory with the relevant
  ## files and run nicholson there
  if(!is.null(trt.fun.name)){
    trt <- loci.subsets[[trt.fun.name]]
    d <- d[trt(d),]
    print(dim(d))
    olds <- read.table(paste(f,sfile,sep='/'))
    news <- data.frame(x=olds[rownames(d),],row.names=rownames(d))
    f <- paste(f,trt.fun.name,sep='/')
    ##f <- paste(f,deparse(substitute(trt)),sep='/')
    system(paste("mkdir",f))
    write.table(news,paste(f,sfile,sep='/'))
  }
  write.nicholson(d*N,paste(f,yfile,sep='/'))
  write.nicholson(d*0+N,paste(f,nfile,sep='/'))
  cat(paste(nrow(d),ncol(d),yfile,nfile,collapse=' '),
      paste(1000,50,5000,50,1000,collapse=' '),
      paste(0.4,1.25,0.06,1.25,0.20,0.45,4501),
      1,
      sep='\n',
      file=paste(f,'nicholson.inp',sep='/'))
  cmd <- paste("cd",f,"&& nicholson")
  if(trt.fun.name=="all")cmd <- paste(cmd,"&&cd ..&& ln -s all/* .")
  print(cmd)
  system(cmd)
}
##debug(run.nicholson)

fit.nicholson <- function(ids){
  for(N in ids)
    for(FUN in names(loci.subsets))
      run.nicholson(result.list[[N]]$alpha.sim,
                    100,
                    paste("results/",N,sep=''),
                    FUN)
}
