nicholsonppp <- function
### Fit the Nicholson model and calculate PPP-values, using MCMC
### implented using fast fortran code.
(Y_OBS,
### Matrix of allele frequencies, one for each (locus,population)
### pair. Used for nmrk, npop, Y_OBS.
 seed=4501,
### Random seed for mersenne twister.
 nvaleurs=1000,
### Number of observations to sample in MCMC.
 thin=50,
### Thinning rate.
 burn_in=5000,
### Burn in period length.
 npilot=50,
### Max number of pilot runs.
 pilot_length=1000,
### Pilot run length.
 delta_a_init=0.4,
### Random walk delta for alpha.
 delta_p_init=1.25,
### Random walk delta for pi.
 delta_c_init=0.06,
### Random walk delta for c.
 rate_adjust=1.25,
### Pilot run adjustment rate.
 acc_inf=0.2,
### Target acceptance rate for rejection sampling, min value.
 acc_sup=0.45,
### Target acceptance rate for rejection sampling, max value.
 out_option=1,
### How much detail to print?
 N_OBS=100
### Population size.
 ){
  nmrk <- nrow(Y_OBS)
  npop <- ncol(Y_OBS)
  return_ppp <- rep(-2.5,nmrk)
  return_a <- rep(-2.5,nmrk*npop)
  return_c <- rep(-2.5,npop)
  return_p <- rep(-2.5,nmrk)
  N_OBS <- matrix(N_OBS,nrow=nrow(Y_OBS),ncol=ncol(Y_OBS))
  YY <- Y_OBS
  NN <- N_OBS
  ## order of argument names here should match fortran code
  fargs <- list(integer=c("npop","nmrk","seed","nvaleurs","thin","burn_in",
                  "npilot","pilot_length","out_option","YY","NN"),
                single=c("delta_a_init","delta_p_init","delta_c_init",
                  "rate_adjust","acc_inf","acc_sup",
                  "return_ppp","return_a","return_c","return_p"))
  fc <- list(as.name(".Fortran"),"fitnicholsonppp")
  for(N in names(fargs))for(a in fargs[[N]]){
    fc[[L <- length(fc)+1]] <- call(paste("as.",N,sep=""),as.name(a))
    names(fc)[L] <- a
  }
  names(fc)[2] <- ""
  cc <- as.call(fc)
  ## debugging:
  ##print(cc)
  ##print(sapply(cc,function(x)if(!is.character(x))eval(x) else x))
  res <- eval(cc)
  ## post-process...
  pre <- "return_"
  rnames <- grep(pre,names(res),val=TRUE)
  ret <- res[rnames]
  names(ret) <- sub(pre,"",rnames)
  ret$a <- matrix(ret$a,nrow=nmrk,ncol=npop)
  ret
}
