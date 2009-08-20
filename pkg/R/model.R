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
 out_option=1
 ){
  nmrk <- nrow(Y_OBS)
  npop <- ncol(Y_OBS)
  .Fortran("fitnicholsonppp",
           as.integer(npop),
           as.integer(nmrk),
           as.integer(seed),
           as.integer(nvaleurs),
           as.integer(thin),
           as.integer(burn_in),
           as.integer(npilot),
           as.integer(pilot_length),
           as.integer(out_option),
           as.double(Y_OBS),
           as.double(N_OBS),
           as.double(delta_a_init),
           as.double(delta_p_init),
           as.double(delta_c_init),
           as.double(rate_adjust),
           as.double(acc_inf),
           as.double(acc_sup)
           )
}
