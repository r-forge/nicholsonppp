!version 2.0 30/01/09
!on fait les acceptance rate (et les delta) par locus, par pop... plutot qu'une moyenne generale
!Chnages: on vire le coefficient binomial dans l'update des alphas (ca se simplifie: a priori on doit gagner du temps et limiter l'overflow
!Changement: calcul des PPPvalues
include 'prob_module.f90'
include 'utils_stat.f90'
program nicholson

  use prob_mod !fonctions diverse pour nombre aleatoires et pdf
  use utils_stat !fonctions diverses pour les summary statistics

  implicit none

  integer:: err, i_thin, npop, nmrk, pop, mrk, tmp, tst, tst_a, tst_p,tst_c, iter , pilot, &
       seed, nvaleurs, thin,burn_in,npilot,pilot_length , out_option
  integer, allocatable :: Y_OBS(:,:), N_OBS(:,:) , RANGS(:)          
  real, allocatable :: INITS_A(:,:),INITS_P_I(:),INITS_C(:) , & 
       RES_A(:,:,:),RES_PI(:,:),RES_C(:,:) , PPPval(:,:) , &! BPval(:,:,:) , &
       delta_p(:),delta_c(:),delta_a(:,:), acc_a(:,:),acc_p(:),acc_c(:) , &
       mean_cij(:,:,:), mean_a(:,:,:),mean_c(:,:),mean_p(:,:)
  real :: tmp_mod,tmp_mean,delta_a_init,delta_p_init,delta_c_init,a_up, p_up, c_up , &
       rate_adjust , acc_tol=0.005, acc_inf , acc_sup , tmp_min, tmp_max ,beta_pi=0.7 , tmp_chi2 ,tmp_pval ,&
       tmp_t_obs, tmp_t_rep, tmp_y_rep,tmp_t,tmp_vij
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call sgrnd(seed) !seed du mersenne twister
  allocate(INITS_A(nmrk,npop) , INITS_P_I(nmrk) , INITS_C(npop) )

  allocate(delta_p(nmrk),delta_c(npop),delta_a(nmrk,npop), acc_a(nmrk,npop),acc_p(nmrk),acc_c(npop)) 
  !initialisation des delta
  delta_p(:)=delta_p_init
  delta_c(:)=delta_c_init
  delta_a(:,:)=delta_a_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!etudes pilotes pour ajuster les delta (acceptance rate entre 0.25 et 0.45)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !open(100,file='adj_acc_a.out',status='unknown')
  !open(101,file='adj_acc_p.out',status='unknown')
  !open(102,file='adj_acc_c.out',status='unknown')

  !write (100,*) 'ITER MRK POP DELTA RATE'
  !write (101,*) 'ITER MRK DELTA RATE'
  !write (102,*) 'ITER POP DELTA RATE'

  pilot=1 ; tst=1 !on elimine le burn-in a vide car on reinitialise chaque fois
  do while (pilot<=npilot .and. tst /= 0)

!!!!Reinitialisation des inits
     INITS_P_I(:)=0.0 !on initialise à la moyenne des A_IJ
     do mrk=1,nmrk
        do pop=1,npop
           INITS_A(mrk,pop)=(Y_OBS(mrk,pop)+0.0)/(N_OBS(mrk,pop)+0.0)
           ! INITS_P_I(mrk)=INITS_P_I(mrk) + INITS_A(mrk,pop)*N_OBS(mrk,pop)
        end do
        INITS_P_I(mrk)=sum(INITS_A(mrk,:))/npop       !INITS_P_I(mrk) / sum(N_OBS(mrk,:))
        INITS_P_I(mrk)=max(0.001,INITS_P_I(mrk))
        INITS_P_I(mrk)=min(0.999,INITS_P_I(mrk))
     end do

     !initialisation des c
     do pop=1,npop
        INITS_C(pop)=0.
        do mrk=1,nmrk
           tmp_min=max(0.001,INITS_P_I(mrk))

           INITS_C(pop)= INITS_C(pop)+ ((INITS_A(mrk,pop)-INITS_P_I(mrk))**2)/(INITS_P_I(mrk)*(1-INITS_P_I(mrk)))
        end do
        INITS_C(pop)= INITS_C(pop)/nmrk
     end do

     print *,'###########INITS VALUES#############'
     print *,'C inits values',INITS_C
     tmp_min=REAL_MIN(INITS_P_I(:)) ; tmp_max=REAL_MAX(INITS_P_I(:)) 
     print *,'min des PI=',tmp_min,'max des PI=',tmp_max,'mean des PI=',sum(INITS_P_I(:))/nmrk
     print *,'####################################'
     !INITS_C(:)=0.05

     acc_a(:,:)=0. ; acc_p(:)=0. ; acc_c(:)=0. 

     print *,'ETUDE PILOTE: ',pilot
     do iter=1,pilot_length
        if(mod(iter,pilot_length/10)==0) print *,'  iteration=',iter

        do mrk=1,nmrk
           do pop=1,npop   
              call update_alpha(Y_OBS(mrk,pop),N_OBS(mrk,pop),INITS_A(mrk,pop),INITS_C(pop),INITS_P_I(mrk),delta_a(mrk,pop),tmp,a_up)
              acc_a(mrk,pop)=acc_a(mrk,pop)+tmp
              INITS_A(mrk,pop)=a_up
           end do
        end do

        do mrk=1,nmrk
           call update_p(INITS_P_I(mrk),INITS_A(mrk,:),INITS_C,delta_p(mrk),tmp,p_up) 
           acc_p(mrk)=acc_p(mrk)+tmp
           INITS_P_I(mrk)=p_up
        end do

        do pop=1,npop
           call update_c(INITS_C(pop),INITS_A(:,pop),INITS_P_I,delta_c(pop),tmp,c_up) 
           acc_c(pop)=acc_c(pop)+tmp 
           INITS_C(pop)=c_up

        end do

     end do

     if(pilot>0) then ! on fait une etude pilote à vide mini-burn-in
        tst=0
        tst_a=0 ; tst_p=0 ; tst_c=0

        do mrk=1,nmrk

           !ajustement des pi
           ! write(101,'(2(i5,1x),2(f8.6,1x))') pilot,mrk,delta_p(mrk),acc_p(mrk)/pilot_length
           if(acc_p(mrk)/pilot_length>acc_sup) then
              delta_p(mrk)=delta_p(mrk)*rate_adjust
              tst_p=tst_p+1
           end if
           if(acc_p(mrk)/pilot_length<acc_inf) then
              delta_p(mrk)=delta_p(mrk)/rate_adjust
              tst_p=tst_p+1
           end if

           do pop=1,npop
              !ajustement des c
              if(mrk == 1)  then
                 !   write(102,'(2(i5,1x),2(f8.6,1x))') pilot,pop,delta_c(pop),acc_c(pop)/pilot_length 
                 if(acc_c(pop)/pilot_length>acc_sup) then
                    delta_c(pop)=delta_c(pop)*rate_adjust
                    tst_c=tst_c+1
                 end if
                 if(acc_c(pop)/pilot_length<acc_inf) then
                    delta_c(pop)=delta_c(pop)/rate_adjust
                    tst_c=tst_c+1
                 end if
              end if


              !  write(100,'(3(i5,1x),2(f8.6,1x))') pilot,mrk,pop,delta_a(mrk,pop),acc_a(mrk,pop)/pilot_length 
              if(acc_a(mrk,pop)/pilot_length>acc_sup) then 
                 delta_a(mrk,pop)=delta_a(mrk,pop)*rate_adjust
                 tst_a=tst_a+1
              end if
              if(acc_a(mrk,pop)/pilot_length<acc_inf) then
                 delta_a(mrk,pop)=delta_a(mrk,pop)/rate_adjust
                 tst_a=tst_a+1
              end if

           end do
        end do

        print *,'  Acceptance Rate not achieved for ',tst_a,' out of ',nmrk*npop ,' alphas' 
        print *,'  Acceptance Rate not achieved for ',tst_p,' out of ',nmrk ,' p' 
        print *,'  Acceptance Rate not achieved for ',tst_c,' out of ',npop ,' c'

        if( ((tst_a+0.0)/(nmrk*npop))>acc_tol .or.  ((tst_p+0.0)/nmrk)>acc_tol .or. tst_c>0 ) tst=1 !!condition pour continuer à ajuster

        print *,'      Mean Acceptance Rate Alpha= ',sum(acc_a(:,:))/(pilot_length*nmrk*npop),' mean delta_a= ',sum(delta_a(:,:))/(nmrk*npop)
        print *,'      Mean Acceptance Rate Pi   = ',sum(acc_p(:))/(pilot_length*nmrk),' mean delta_p= ',sum(delta_p(:))/nmrk
        print *,'      Mean Acceptance Rate C    = ',sum(acc_c(:))/(pilot_length*npop),' mean delta_c= ',sum(delta_c(:))/npop

     end if

     pilot=pilot+1


     print *,'C node values',INITS_C
     tmp_min=REAL_MIN(INITS_P_I(:)) ; tmp_max=REAL_MAX(INITS_P_I(:)) 
     print *,'min des PI=',tmp_min,'max des PI=',tmp_max,'mean des PI=',sum(INITS_P_I(:))/nmrk

  end do

  print *,'FIN AJUSTEMENT DES Delta '
  print *,'  Mean Delta_a= ',sum(delta_a(:,:))/(nmrk*npop)
  print *,'  Mean Delta_p= ',sum(delta_p(:))/nmrk
  print *,'  Mean Delta_c= ',sum(delta_c(:))/npop

  ! close(100) ; close(101) ; close(102)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! BURN IN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print *,'BEGIN BURN_IN PERIOD'
  do iter=1,burn_in

     do mrk=1,nmrk
        do pop=1,npop   
           call update_alpha(Y_OBS(mrk,pop),N_OBS(mrk,pop),INITS_A(mrk,pop),INITS_C(pop),INITS_P_I(mrk),delta_a(mrk,pop),tmp,a_up)
           INITS_A(mrk,pop)=a_up
        end do
     end do

     do mrk=1,nmrk
        call update_p(INITS_P_I(mrk),INITS_A(mrk,:),INITS_C,delta_p(mrk),tmp,p_up) 
        INITS_P_I(mrk)=p_up
     end do

     do pop=1,npop
        call update_c(INITS_C(pop),INITS_A(:,pop),INITS_P_I,delta_c(pop),tmp,c_up) 
        INITS_C(pop)=c_up
     end do

  end do

  print *,'END BURN_IN PERIOD'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! CHAINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  acc_a(:,:)=0. ; acc_p(:)=0. ; acc_c(:)=0.
  if(out_option==0) then
     open(3,file='res_alpha.out',status='unknown')
     open(4,file='res_pi.out',status='unknown')
  end if

  open(5,file='res_c.out',status='unknown') !quelle que soit l'option on sort les c car pas possible avec le sum_stat sur res_mcmc trie de faire les stats sur les c
  !open(6,file='res_mcmc.out',status='unknown')

  print *,'BEGIN CHAIN'

  if(out_option/=2) then
     allocate(mean_cij(nmrk,npop,1:2)) ; mean_cij(:,:,:)=0.0
     allocate(mean_a(nmrk,nvaleurs,1:2),mean_c(npop,1:2),mean_p(nmrk,1:2))
     mean_a(:,:,:)=0.0 ;  mean_p(:,:)=0.0 ; mean_c(:,:)=0.0 
  end if

  if(out_option==0) then
     allocate(RES_A(nmrk,npop,nvaleurs))
     allocate(RES_PI(nmrk,nvaleurs) , RES_C(npop,nvaleurs) )
  end if

  allocate(PPPval(nmrk,npop+1)) !le dernier correspond à al totale
  PPPval(:,:)=0.0
  !allocate(BPval(nmrk,npop,1:6)) ! BPval à 3.841 et 6.635 , mean_PVALUE(CHI2) , sd_Pvalue , min,max
  !BPval(:,:,:)=0.0
  !BPval(:,:,5)=10.0

  !open(1000,file='tt',status='unknown')

  do iter=1,nvaleurs

     if(mod(iter,nvaleurs/100)==0) print *,'  iteration=',iter
!!!!!!!!!!!!!!!!!!!!!
     !impression resultats et calucl des fstij
     if(out_option==0) write(4,'(100000(f12.8,1x))') INITS_P_I(1:nmrk)
     write(5,'(100(f12.8,1x))') INITS_C(1:npop)
     do mrk=1,nmrk
        if(out_option==0) write(3,'(100(f12.8,1x))') INITS_A(mrk,:)
        !calcul des moyennes et moy2 des fstij et des PPPVal
        tmp_t=0.0
        do pop=1,npop
           tmp_vij=INITS_P_i(mrk)*(1-INITS_P_i(mrk))*((N_OBS(mrk,pop)-1)*INITS_C(pop)+1)/N_OBS(mrk,pop)
           !*(1-INITS_P_i(mrk))*(INITS_C(pop)+1/N_OBS(mrk,pop))
           tmp_t_obs=(((Y_OBS(mrk,pop)+0.0)/(N_OBS(mrk,pop)+0.0) -INITS_P_i(mrk))**2)/tmp_vij

           tmp_y_rep=random_binomial2(N_OBS(mrk,pop),INITS_A(mrk,pop),.TRUE.)
           tmp_t_rep=(((tmp_y_rep/N_OBS(mrk,pop)) -INITS_P_i(mrk))**2)/tmp_vij
           ! write(1000,*) tmp_y_rep, Y_OBS(mrk,pop)/N_OBS(mrk,pop), tmp_y_rep/N_OBS(mrk,pop)

           if(tmp_t_rep > tmp_t_obs) PPPval(mrk,pop)=PPPval(mrk,pop)+1.0
           tmp_t=tmp_t_rep-tmp_t_obs + tmp_t

           ! tmp_chi2= ((INITS_A(mrk,pop)-INITS_P_i(mrk))**2)/(INITS_P_i(mrk)*(1-INITS_P_i(mrk))*INITS_C(pop))
           ! call chi_square_cdf(tmp_chi2,1.0,tmp_pval)
           ! if(tmp_pval >= 0.95) BPval(mrk,pop,1)=BPval(mrk,pop,1)+1.0
           ! if(tmp_pval >= 0.99) BPval(mrk,pop,2)=BPval(mrk,pop,2)+1.0                
           ! BPval(mrk,pop,3)=(BPval(mrk,pop,3)*(iter-1)+1-tmp_pval)/iter
           ! BPval(mrk,pop,4)=(BPval(mrk,pop,4)*(iter-1)+(1-tmp_pval)**2)/iter
           ! if(tmp_pval >BPval(mrk,pop,6)) BPval(mrk,pop,6)=tmp_pval
           ! if(tmp_pval <BPval(mrk,pop,5)) BPval(mrk,pop,5)=tmp_pval                 
           !write(6,'(2(i5,1x),3(f8.6,1x))') mrk,pop,INITS_A(mrk,pop),INITS_P_I(mrk), INITS_C(pop) !, 1-tmp_pval &
           ! ((INITS_A(mrk,pop)-INITS_P_i(mrk))**2)/(INITS_P_i(mrk)*(1-INITS_P_i(mrk))), tmp_chi2 , 
           if(out_option/=2) then 
              mean_cij(mrk,pop,1)=(mean_cij(mrk,pop,1)*(iter-1)+((INITS_A(mrk,pop)-INITS_P_i(mrk))**2)/(INITS_P_i(mrk)*(1-INITS_P_i(mrk))))/iter
              mean_cij(mrk,pop,2)=(mean_cij(mrk,pop,2)*(iter-1)+(((INITS_A(mrk,pop)-INITS_P_i(mrk))**2)/(INITS_P_i(mrk)*(1-INITS_P_i(mrk))))**2)/iter
           end if
        end do
        if(tmp_t > 0.0) PPPval(mrk,npop+1)=PPPval(mrk,npop+1)+1.0
     end do
!!!!!!!!!!!!!!!!!!!!!
     if(out_option==0)  then
        RES_A(:,:,iter)=INITS_A(:,:)
        RES_PI(:,iter)=INITS_P_I(:)
        RES_C(:,iter)=INITS_C(:)
     end if

     if(out_option==1) then
        do mrk=1,nmrk
           mean_p(mrk,1)=(mean_p(mrk,1)*(iter-1)+INITS_P_I(mrk))/iter
           mean_p(mrk,2)=(mean_p(mrk,2)*(iter-1)+(INITS_P_I(mrk))**2)/iter
           do pop=1,npop
              if(mrk==1) then
                 mean_c(pop,1)=(mean_c(pop,1)*(iter-1)+INITS_C(pop))/iter
                 mean_c(pop,2)=(mean_c(pop,2)*(iter-1)+(INITS_C(pop))**2)/iter
              end if
              mean_a(mrk,pop,1)=(mean_a(mrk,pop,1)*(iter-1)+INITS_A(mrk,pop))/iter
              mean_a(mrk,pop,2)=(mean_a(mrk,pop,2)*(iter-1)+(INITS_A(mrk,pop))**2)/iter                  
           end do
        end do
     end if
!!!!!!!!!!!!!!!!!!!!!!

     do i_thin=1,thin
        do mrk=1,nmrk
           do pop=1,npop   
              call update_alpha(Y_OBS(mrk,pop),N_OBS(mrk,pop),INITS_A(mrk,pop),INITS_C(pop),INITS_P_I(mrk),delta_a(mrk,pop),tmp,a_up)
              acc_a(mrk,pop)=acc_a(mrk,pop)+tmp
              INITS_A(mrk,pop)=a_up
           end do
        end do

        do mrk=1,nmrk
           call update_p(INITS_P_I(mrk),INITS_A(mrk,:),INITS_C,delta_p(mrk),tmp,p_up) 
           acc_p(mrk)=acc_p(mrk)+tmp
           INITS_P_I(mrk)=p_up
        end do

        do pop=1,npop
           call update_c(INITS_C(pop),INITS_A(:,pop),INITS_P_I,delta_c(pop),tmp,c_up) 
           acc_c(pop)=acc_c(pop)+tmp 
           INITS_C(pop)=c_up
        end do

     end do

  end do

  close(3) ; close(4) ; close(5)  !close(6)

  deallocate(INITS_A,INITS_C,INITS_P_I,Y_OBS,N_OBS)

  print *,'Mean Final Acceptance Rate Alpha= ',sum(acc_a(:,:))/((nvaleurs+0.0)*(thin+0.0)*(nmrk+0.0)*(npop+0.0)),' mean delta_a= ',sum(delta_a(:,:))/(nmrk*npop)
  print *,'Mean Final Acceptance Rate Pi   = ',sum(acc_p(:))/((nvaleurs+0.0)*(thin+0.0)*(nmrk+0.0)),' mean delta_p= ',sum(delta_p(:))/nmrk
  print *,'Mean Final Acceptance Rate C    = ',sum(acc_c(:))/((nvaleurs+0.0)*(thin+0.0)*(npop+0.0)),' mean delta_c= ',sum(delta_c(:))/npop

!!!!!
  !!impression des acceptance rate par locus...
!!!!!

  open(100,file='fin_acc_a.out',status='unknown')
  open(101,file='fin_acc_p.out',status='unknown')
  open(102,file='fin_acc_c.out',status='unknown')

  write (100,*) 'MRK POP DELTA RATE'
  write (101,*) 'MRK DELTA RATE'
  write (102,*) 'POP DELTA RATE'

  do mrk=1,nmrk
     write(101,'(1(i5,1x),2(f8.6,1x))') mrk,delta_p(mrk),acc_p(mrk)/(nvaleurs*thin)

     do pop=1,npop

        if(mrk == 1)  write(102,'(1(i5,1x),2(f8.6,1x))') pop,delta_c(pop),acc_c(pop)/(nvaleurs*thin)

        write(100,'(2(i5,1x),2(f8.6,1x))') mrk,pop,delta_a(mrk,pop),acc_a(mrk,pop)/(nvaleurs*thin)

     end do
  end do

  close(100) ; close(101) ; close(102)


  !!SUMMARY STATISTICS: NOM_code, moyenne, variance, mediane, mode, quantile1%, 5%,25%,75%,95%,99%
  if(out_option==0) then
     allocate(RANGS(nvaleurs))
     open(1,file='summary_alpha.out',status='unknown')
     write (1,*) 'POP MRK MOY VAR MED MODE Q_0.01 Q_0.05 Q_0.25 Q_0.75 Q_0.99'
     do pop=1,npop
        do mrk=1,nmrk
           tmp_mean=sum(RES_A(mrk,pop,:))/nvaleurs
           call TRIIND(RES_A(mrk,pop,:),RANGS)
           call CALMOD (RES_A(mrk,pop,:),RANGS, tmp_mod, tmp) 
           write(1,'(2i5,9(f12.8,1x))') pop,mrk,tmp_mean,sum((RES_A(mrk,pop,:)-tmp_mean)**2)/(nvaleurs-1),FRCTIL(RES_A(mrk,pop,:), RANGS, 0.5) ,&
                tmp_mod, FRCTIL(RES_A(mrk,pop,:), RANGS, 0.01),FRCTIL(RES_A(mrk,pop,:), RANGS, 0.05), &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.25),FRCTIL(RES_A(mrk,pop,:), RANGS, 0.75),&
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.95),FRCTIL(RES_A(mrk,pop,:), RANGS, 0.99)
        end do
     end do
     close(1) 
     deallocate(RES_A)

     open(1,file='summary_pi.out',status='unknown')
     write (1,*) 'MRK MOY VAR MED MODE Q_0.01 Q_0.05 Q_0.25 Q_0.75 Q_0.99'
     do mrk=1,nmrk
        tmp_mean=sum(RES_PI(mrk,:))/nvaleurs !on recycle acc_a pour la moyenne
        call TRIIND(RES_PI(mrk,:),RANGS)
        call CALMOD (RES_PI(mrk,:),RANGS, tmp_mod, tmp) !on recycle acc_p pour le mode
        write(1,'(i5,9(f12.8,1x))') mrk,tmp_mean,sum((RES_PI(mrk,:)-tmp_mean)**2)/(nvaleurs-1),FRCTIL(RES_PI(mrk,:), RANGS, 0.5) ,&
             tmp_mod, FRCTIL(RES_PI(mrk,:), RANGS, 0.01),FRCTIL(RES_PI(mrk,:), RANGS, 0.05), &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.25),FRCTIL(RES_PI(mrk,:), RANGS, 0.75),&
             FRCTIL(RES_PI(mrk,:), RANGS, 0.95),FRCTIL(RES_PI(mrk,:), RANGS, 0.99)
     end do
     close(1)     
     deallocate(RES_PI)  

     open(1,file='summary_c.out',status='unknown')
     write (1,*) 'POP MOY VAR MED MODE Q_0.01 Q_0.05 Q_0.25 Q_0.75 Q_0.99'
     do pop=1,npop
        tmp_mean=sum(RES_C(pop,:))/nvaleurs !on recycle acc_a pour la moyenne
        call TRIIND(RES_C(pop,:),RANGS)
        call CALMOD (RES_C(pop,:),RANGS, tmp_mod, tmp) !on recycle acc_p pour le mode
        write(1,'(i5,9(f12.8,1x))') pop,tmp_mean,sum((RES_C(pop,:)-tmp_mean)**2)/(nvaleurs-1),FRCTIL(RES_C(pop,:), RANGS, 0.5) ,&
             tmp_mod, FRCTIL(RES_C(pop,:), RANGS, 0.01),FRCTIL(RES_C(pop,:), RANGS, 0.05), &
             FRCTIL(RES_C(pop,:), RANGS, 0.25),FRCTIL(RES_C(pop,:), RANGS, 0.75),&
             FRCTIL(RES_C(pop,:), RANGS, 0.95),FRCTIL(RES_C(pop,:), RANGS, 0.99)
     end do
     close(1)     
     deallocate(RES_C) 
  end if


  if(out_option==1) then ! on ecrit juste les moyennes et les sd
     open(1,file='summary_alpha.out',status='unknown')
     write (1,*) 'POP MRK MOY VAR'

     open(2,file='summary_pi.out',status='unknown')
     write (2,*) 'MRK MOY VAR'

     open(3,file='summary_c.out',status='unknown')
     write (3,*) 'POP MOY VAR'

     do pop=1,npop
        write(3,'(1i5,2(f12.8,1x))') pop,mean_c(pop,1),mean_c(pop,2)-(mean_c(pop,1))**2

        do mrk=1,nmrk
           write(1,'(2i5,2(f12.8,1x))') pop,mrk,mean_a(mrk,pop,1),mean_a(mrk,pop,2)-(mean_a(mrk,pop,1))**2

           if(pop==1) write(2,'(1i5,2(f12.8,1x))') mrk,mean_p(mrk,1),mean_p(mrk,2)-(mean_p(mrk,1))**2    

        end do
     end do
     close(1) ; close(2) ; close(3)

  end if


  if(out_option/=2) then !
     !!impressions des FST_ij
     open(1,file='summary_cij.out',status='unknown')
     write (1,*) 'POP MRK MOY VAR'
     do pop=1,npop
        do mrk=1,nmrk
           write(1,'(2i5,2(f12.8,1x))') pop,mrk,mean_cij(mrk,pop,1),mean_cij(mrk,pop,2)-(mean_cij(mrk,pop,1))**2
        end do
     end do
     close(1)
  end if
!!!!!!!!!!!!!!!!
  !impression des PPvalues (par pop/marker et par marker
  PPPval(:,:)=PPPval(:,:)/nvaleurs
  open(1,file='PPPval.out',status='unknown')
  write (1,*) 'POP MRK PPPVAL PPPVALtot'
  do pop=1,npop
     do mrk=1,nmrk
        write(1,'(2(i5,1x),2(f12.8,1x))') pop,mrk,PPPval(mrk,pop),PPPval(mrk,npop+1)
     end do
  end do
  close(1)



!!!!!!!!!!!!!!!!



  !open(1,file='BPval.out',status='unknown')
  !write (1,*) 'POP MRK BP5 BP1 Mean_BP VAR_P MIN_P MAX_P '
  !  do pop=1,npop
  !   do mrk=1,nmrk
  !   write(1,'(2(i5,1x),6(f12.8,1x))') pop,mrk,BPval(mrk,pop,1)/nvaleurs,BPval(mrk,pop,2)/nvaleurs,BPval(mrk,pop,3), &
  !     BPval(mrk,pop,4)-BPval(mrk,pop,3)**2,BPval(mrk,pop,5),BPval(mrk,pop,6)
  !   end do
  !  end do
  ! close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!update des alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_alpha (Y_obs,N_obs,a_cur,C_j,P_i,delta_a,accept,a_out) !delta_a=taille de l'intervalle autour de la proposal de a 
    implicit none

    integer, intent(in) :: Y_obs, N_obs !locus pop consideres
    real, intent(in) :: C_j, P_i , delta_a ,a_cur

    integer, intent(out) :: accept !accept=0 si accepte, 1 sinon
    real, intent(out) :: a_out

    real :: b_inf, b_sup, a_prop, inv_fwd, inv_bwd, rnd_prop , rnd_accep , diff_log, log_cur , log_new !rnd_prop pour la proposal e rnd accep pour l'acceptance

    rnd_prop=grnd()

    b_inf=max(0.,a_cur-delta_a/2) ; b_sup=min(1.,a_cur+delta_a/2)
    a_prop=b_inf+rnd_prop*(b_sup - b_inf)

    !1/p(q,q') et 1/p(q',q): forward et backward prop probability
    inv_fwd=b_sup - b_inf
    inv_bwd=min(1.,a_prop+delta_a/2) - max(0.,a_prop-delta_a/2)

    !pour evaluer on regarde juste le changement au locus: ca sert à rien de calculer toute la vraisemblance
    !log_cur=log_lik(Y_obs,N_obs,A_i,C_j,P_i)
    !A_i(i_cur,j_cur)=a_prop
    !log_new=log_lik(Y_obs,N_obs,A_i,C_j,P_i)

    log_cur=log_norm_pdf(a_cur,P_i,P_i*(1-P_i)*C_j) + Y_obs*log(a_cur) + (N_obs-Y_obs)*log(1-a_cur) 
    !log(norm_pdf(a_cur,P_i,P_i*(1-P_i)*C_j))+log(binomial_pdf(Y_obs,N_obs,a_cur))   
    log_new= log_norm_pdf(a_prop,P_i,P_i*(1-P_i)*C_j) + Y_obs*log(a_prop) + (N_obs-Y_obs)*log(1-a_prop) 
    !log(norm_pdf(a_prop,P_i,P_i*(1-P_i)*C_j))+log(binomial_pdf(Y_obs,N_obs,a_prop))   


    !diff_log=log(P_new*P_bwd/P_old*P_fwd) 
    diff_log=log_new - log_cur + log(inv_fwd) - log(inv_bwd) 

    !on refuse le saut si (rnd_accept)>min(1,P) soit log(rnd_accep)>min(0,diff_log) 
    ! or log(rnd_accep)<=0 car rnd_accep<=1 donc meme si diff_log>0 on peut garder la condition sur 0

    rnd_accep=grnd()
    if(log(rnd_accep)>diff_log) then 
       accept=0
       a_out=a_cur
    else
       accept=1 ! on avait deja initialise à a_prop
       a_out=a_prop
    end if

  end subroutine update_alpha

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!update des C_j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_c (c_cur,A_ij,P_i,delta_c,accept,c_out) !delta_a=taille de l'intervalle autour de la proposal de a 
    implicit none

    !integer, intent(in) :: j_cur !pop consideree
    real, intent(in), dimension(:) :: A_ij , P_i 
    real, intent(in) :: delta_c , c_cur

    real, intent(out) :: c_out
    integer, intent(out) :: accept !0 si accepte, 1 sinon

    integer :: tmp_nmrk , tmp_i
    real :: b_inf, b_sup, c_prop, inv_fwd, inv_bwd, rnd_prop , rnd_accep ,&
         diff_log, log_cur , log_new !rnd_prop pour la proposal e rnd accep pour l'acceptance

    log_cur=0. ; log_new=0. !uniforme sur les c

    tmp_nmrk=size(P_i)

    rnd_prop=grnd()
    b_inf=max(0.,c_cur-delta_c/2) ; b_sup=min(1.,c_cur+delta_c/2)
    c_prop=b_inf+rnd_prop*(b_sup - b_inf)

    !1/p(q,q') et 1/p(q',q): forward et backward prop probability
    inv_fwd=b_sup - b_inf
    inv_bwd=min(1.,c_prop+delta_c/2) - max(0.,c_prop-delta_c/2)

    do tmp_i=1, tmp_nmrk
       log_cur=log_cur + log_norm_pdf(A_ij(tmp_i),P_i(tmp_i),P_i(tmp_i)*(1-P_i(tmp_i))*c_cur) 
       !log(norm_pdf(A_ij(tmp_i),P_i(tmp_i),P_i(tmp_i)*(1-P_i(tmp_i))*c_cur)) 
       log_new=log_new + log_norm_pdf(A_ij(tmp_i),P_i(tmp_i),P_i(tmp_i)*(1-P_i(tmp_i))*c_prop)
       !log(norm_pdf(A_ij(tmp_i),P_i(tmp_i),P_i(tmp_i)*(1-P_i(tmp_i))*c_prop)) 

    end do

    !diff_log=log(P_new*P_bwd/P_old*P_fwd) 
    diff_log=log_new - log_cur + log(inv_fwd) - log(inv_bwd) 


    ! write (*,'(4f14.8)') log_new, log_cur , log(inv_fwd) , log(inv_bwd) 



    !on refuse le saut si (rnd_accept)>min(1,P) soit log(rnd_accep)>min(0,diff_log) 
    ! or log(rnd_accep)<=0 car rnd_accep<=1 donc meme si diff_log>0 on peut garder la condition sur 0

    rnd_accep=grnd()
    if(log(rnd_accep)>diff_log) then 
       accept=0
       c_out=c_cur 
    else
       accept=1 ! on avait deja initialise à a_prop
       c_out=c_prop
    end if

  end subroutine update_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!update des P_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_p (p_cur,A_ij,C_j,delta_p,accept,p_out) !delta_a=taille de l'intervalle autour de la proposal de a 
    implicit none

    !integer, intent(in) :: i_cur !locus considere
    real, intent(in), dimension(:) :: A_ij , C_j
    real, intent(in) :: delta_p , p_cur

    integer, intent(out) :: accept !0 si accepte, 1 sinon
    real, intent(out) :: p_out

    real :: b_inf, b_sup, p_prop, inv_fwd, inv_bwd, rnd_prop , rnd_accep , diff_log, log_cur , log_new !rnd_prop pour la proposal e rnd accep pour l'acceptance
    integer :: tmp_npop , tmp_i

    tmp_npop=size(C_j)

    rnd_prop=grnd()
    b_inf=max(0.001,p_cur-delta_p/2) ; b_sup=min(0.999,p_cur+delta_p/2)
    p_prop=b_inf+rnd_prop*(b_sup - b_inf)


    !1/p(q,q') et 1/p(q',q): forward et backward prop probability
    inv_fwd=b_sup - b_inf
    inv_bwd=min(0.999,p_prop+delta_p/2) - max(0.001,p_prop-delta_p/2)

    log_cur=(beta_pi - 1.0)*(log(p_cur) + log(1-p_cur))
    !log(beta_pdf(p_cur,beta_pi,beta_pi))
    log_new=(beta_pi - 1.0)*(log(p_prop) + log(1-p_prop))
    !log(beta_pdf(p_prop,beta_pi,beta_pi))

    do tmp_i=1,tmp_npop
       log_cur=log_cur + log_norm_pdf(A_ij(tmp_i),p_cur,p_cur*(1-p_cur)*C_j(tmp_i))
       log_new=log_new + log_norm_pdf(A_ij(tmp_i),p_prop,p_prop*(1-p_cur)*C_j(tmp_i))
    end do



    !diff_log=log(P_new*P_bwd/P_old*P_fwd) 
    diff_log=log_new - log_cur + log(inv_fwd) - log(inv_bwd) 

    !on refuse le saut si (rnd_accept)>min(1,P) soit log(rnd_accep)>min(0,diff_log) 
    ! or log(rnd_accep)<=0 car rnd_accep<=1 donc meme si diff_log>0 on peut garder la condition sur 0

    rnd_accep=grnd()
    if(log(rnd_accep)>diff_log) then 
       accept=0
       p_out=p_cur 
    else
       accept=1 ! on avait deja initialise à a_prop
       p_out=p_prop
    end if

  end subroutine update_p


end program nicholson




