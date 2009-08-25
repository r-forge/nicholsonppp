  close(3) ; close(4) ; close(5)  !close(6)

  deallocate(INITS_A,INITS_C,INITS_P_I,Y_OBS,N_OBS)

  !print *,'Mean Final Acceptance Rate Alpha= ', &
  !     sum(acc_a(:,:))/((nvaleurs+0.0)*(thin+0.0)*(nmrk+0.0)*(npop+0.0)), &
  !     ' mean delta_a= ',sum(delta_a(:,:))/(nmrk*npop)
  !print *,'Mean Final Acceptance Rate Pi   = ', &
  !     sum(acc_p(:))/((nvaleurs+0.0)*(thin+0.0)*(nmrk+0.0)), &
  !     ' mean delta_p= ',sum(delta_p(:))/nmrk
  !print *,'Mean Final Acceptance Rate C    = ', &
  !     sum(acc_c(:))/((nvaleurs+0.0)*(thin+0.0)*(npop+0.0)), &
  !     ' mean delta_c= ',sum(delta_c(:))/npop

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
     write(101,'(1(i5,1x),2(f8.6,1x))') mrk,delta_p(mrk), &
          acc_p(mrk)/(nvaleurs*thin)

     do pop=1,npop

        if(mrk == 1)  write(102,'(1(i5,1x),2(f8.6,1x))') pop, &
             delta_c(pop),acc_c(pop)/(nvaleurs*thin)

        write(100,'(2(i5,1x),2(f8.6,1x))') mrk,pop, &
             delta_a(mrk,pop),acc_a(mrk,pop)/(nvaleurs*thin)

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
           write(1,'(2i5,9(f12.8,1x))') pop,mrk,tmp_mean, &
                sum((RES_A(mrk,pop,:)-tmp_mean)**2)/(nvaleurs-1), &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.5), &
                tmp_mod, &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.01), &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.05), &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.25), &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.75), &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.95), &
                FRCTIL(RES_A(mrk,pop,:), RANGS, 0.99)
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
        write(1,'(i5,9(f12.8,1x))') mrk,tmp_mean, &
             sum((RES_PI(mrk,:)-tmp_mean)**2)/(nvaleurs-1), &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.5) ,&
             tmp_mod, &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.01), &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.05), &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.25), &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.75), &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.95), &
             FRCTIL(RES_PI(mrk,:), RANGS, 0.99)
     end do
     close(1)     
     deallocate(RES_PI)  

     open(1,file='summary_c.out',status='unknown')
     write (1,*) 'POP MOY VAR MED MODE Q_0.01 Q_0.05 Q_0.25 Q_0.75 Q_0.99'
     do pop=1,npop
        tmp_mean=sum(RES_C(pop,:))/nvaleurs !on recycle acc_a pour la moyenne
        call TRIIND(RES_C(pop,:),RANGS)
        call CALMOD (RES_C(pop,:),RANGS, tmp_mod, tmp) !on recycle acc_p pour le mode
        write(1,'(i5,9(f12.8,1x))') pop,tmp_mean, &
             sum((RES_C(pop,:)-tmp_mean)**2)/(nvaleurs-1), &
             FRCTIL(RES_C(pop,:), RANGS, 0.5) ,&
             tmp_mod, &
             FRCTIL(RES_C(pop,:), RANGS, 0.01), &
             FRCTIL(RES_C(pop,:), RANGS, 0.05), &
             FRCTIL(RES_C(pop,:), RANGS, 0.25), &
             FRCTIL(RES_C(pop,:), RANGS, 0.75), &
             FRCTIL(RES_C(pop,:), RANGS, 0.95), &
             FRCTIL(RES_C(pop,:), RANGS, 0.99)
     end do
     close(1)     
     deallocate(RES_C) 
  end if


  end if


  if(out_option/=2) then !
     !!impressions des FST_ij
     open(1,file='summary_cij.out',status='unknown')
     write (1,*) 'POP MRK MOY VAR'
     do pop=1,npop
        do mrk=1,nmrk
           write(1,'(2i5,2(f12.8,1x))') pop,mrk,mean_cij(mrk,pop,1), &
                mean_cij(mrk,pop,2)-(mean_cij(mrk,pop,1))**2
        end do
     end do
     close(1)
  end if


  PPPval(:,:)=PPPval(:,:)/nvaleurs
  do mrk=1,nmrk
     !print *,mrk,PPPvec(mrk)
     PPPvec(mrk)=PPPval(mrk,npop+1)
     !print *,mrk,PPPvec(mrk)
  end do

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

  !print *,"right before the contains"

