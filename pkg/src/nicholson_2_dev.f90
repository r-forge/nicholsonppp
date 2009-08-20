!version 2.0 30/01/09
!on fait les acceptance rate (et les delta) par locus, par pop... plutot qu'une moyenne generale
!Chnages: on vire le coefficient binomial dans l'update des alphas (ca se simplifie: a priori on doit gagner du temps et limiter l'overflow
!Changement: calcul des PPPvalues
include 'fitnicholsonppp.f90'
program nicholson

  use fitnicholsonpppmod

  implicit none

  integer:: mrk, err, npop, nmrk, seed, nvaleurs, thin, burn_in, npilot, pilot_length, out_option
  integer, allocatable :: Y_OBS(:,:), N_OBS(:,:)
  real :: delta_a_init,delta_p_init,delta_c_init, rate_adjust, acc_inf, acc_sup
  character*30 :: Y_file, N_file

!!!!!!!!!!!!!!!!!!!!!
!!!INPUT
!!!!!!!!!!!!!!!!!!!

  open(1,file='nicholson.inp',status='old')
  read(1,*) nmrk,npop,Y_file,N_file
  print *,' No de Marqueurs declares           = ',nmrk
  print *,' No de Populations                  = ',npop
  print *,' Fichiers N_allele Ref              = ',Y_file
  print *,' Fichiers N_allele TOT (2*Nindgen)  = ',N_file
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'           PARAMETRES MCMC'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  read(1,*) nvaleurs,thin,burn_in,npilot,pilot_length
  print *,' Nbre Valeurs Desirees              = ',nvaleurs
  print *,' Thinning Rate                      = ',thin
  print *,' Burn in Period Length              = ',burn_in
  print *,' Max Number of Pilot runs           = ',npilot
  print *,' Pilot run Length                   = ',pilot_length
  print *,''
  read(1,*) delta_a_init,delta_p_init,delta_c_init,rate_adjust, acc_inf , acc_sup, seed
  print *,' Random Walk deltas (alpha,pi and c)= ',delta_a_init,delta_p_init,delta_c_init
  print *,' Pilot Run Adjustment Rate Pilot    = ',rate_adjust
  print *,' Targeted Rejection/Acceptation     = ',acc_inf,'-',acc_sup
  print *,' Seed (Mersenne-Twister)            = ',seed 
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  read(1,*) out_option
  if(out_option==1) print *,'Only Mean and Var will be computed in summary_* files (no output res for iterations except for c)'
  if(out_option==0) print *,'Summary stats will be computed in summary_* files'
  if(out_option==2) print *,'A (huge) res MCMC file will be printed out: no summary stats'
  print *,''
  close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !allocation des valeurs initiales
  allocate(Y_OBS(nmrk,npop), N_OBS(nmrk,npop) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !lecture des observations

  open(1,file=Y_file,status='old')
  open(2,file=N_file,status='old')
  do mrk=1,nmrk
     read(1,*,end=7,iostat=err) Y_OBS(mrk,1:Npop)
     read(2,*,end=7,iostat=err) N_OBS(mrk,1:Npop)
  end do
7   continue
  close(1)
  close(2)

  print *,'Premiere ligne (fichier N) ',N_OBS(1,:)
  print *,'Derniere ligne (fichier N) ',N_OBS(nmrk,:)
  print *,'Premiere ligne (fichier Y) ',Y_OBS(1,:)
  print *,'Derniere ligne (fichier Y) ',Y_OBS(nmrk,:)

  call fitnicholsonppp(npop, nmrk, seed, nvaleurs, thin, burn_in, npilot, pilot_length, out_option, Y_OBS, N_OBS, delta_a_init,delta_p_init,delta_c_init, rate_adjust, acc_inf, acc_sup)
end program nicholson
