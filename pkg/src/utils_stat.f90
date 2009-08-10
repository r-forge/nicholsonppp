module utils_stat

!!!ce module contient un ensemble dutilite pour calucler
!!les moyennes, quantiles,modes,variances...
!!!Il contient egalemetn une fonction de tri

!!!!LISTE:

!     -sous rountine: TRIIND(vecteur des valeurs real , vectur de rang entier): donne le rang de chaque valeur du vecteur
!     -fonction: FRCTIL (XVALT, IRNGT, XORD) : renvoie le fractile d'ordre XORD (0.1,0.05...à choisir) pour le vecteur de reals XVALT dont les rangs (calucles precedemment avec TRIIND) sont dans le vectur IRNGT
!     -sous rountine: CALMOD (XVALT, IRNGT, XMOD, NMOD): renvoie le mode dans XMOD si NMOD=1 sinon renvoie le nombre de modes
!     -fonction: LISOUB (XVALT, FREN) : lissage par moyenne courante avec oubli (determine par le reel FREN) des valuer de XVALT
!     -fonction: LISFEN (XVALT) : lissage par moyenne courante des valuers de XVALT


contains

SUBROUTINE TRIIND (XVALT, IRNGT)
!   triind = tri par interclassement suivant xvalt croissant
REAL, DIMENSION (:)            :: XVALT 
INTEGER, DIMENSION (:)         :: IRNGT
! __________________________________________________________
    INTEGER, DIMENSION (:), ALLOCATABLE :: JWRKT
!
    NVAL = MIN (SIZE (XVALT), SIZE (IRNGT))
    IF (NVAL <= 0) THEN
        RETURN
    ENDIF
!
! .... initialisation du tableau d'indices
!      (on cree des monotonies de longueurs 2)
!
    DO  IIND = 2, NVAL, 2
        IF (XVALT (IIND - 1) < XVALT (IIND)) THEN
            IRNGT (IIND - 1) = IIND - 1
            IRNGT (IIND) = IIND
        ELSE
            IRNGT (IIND - 1) = IIND
            IRNGT (IIND) = IIND - 1
        ENDIF
    ENDDO
    IF (MOD (NVAL, 2) /= 0) THEN
        IRNGT (NVAL) = NVAL
    ENDIF
!
! .... initialisation des monotonies 'a' et 'c'
!
    ALLOCATE (JWRKT (1:NVAL))
    LMTNC = 2
    LMTNA = 2
!
! .... initialisation nouvelles monotonies 
!
    DO
        IF (LMTNA >= NVAL) EXIT
        IWRKF = 0
        LMTNC = 2 * LMTNC
        IWRK = 0
!
! .... passage aux monotonies suivantes 'a', 'b' et 'c'
! on regarde si on deborde, auquel cas on regarde
! si on a encore 2 monotonies a interclasser
!
        DO
          IINDA = IWRKF
          IWRKD = IWRKF + 1
          IWRKF = IINDA + LMTNC
          JINDA = IINDA + LMTNA
          IF (IWRKF >= NVAL) THEN
            IF (JINDA >= NVAL) EXIT
            IWRKF = NVAL
          ENDIF
          IINDB = JINDA
!
! .... on compare la derniere valeur de 'a' et la 1ere de 'b',
!      pour sauter directement si a < b
!  
          IVALA = IRNGT (JINDA)
          IVALB = IRNGT (JINDA + 1)
          IF (XVALT (IVALA) <= XVALT (IVALB)) THEN
              IWRK = IWRKF
              CYCLE    
          ENDIF
!
! .... boucle sur la nouvelle monotonie 'c'
!      qu'on cree dans le tableau wrk
!
          DO
              IF (IWRK >= IWRKF) THEN
!
! .... on recopie le tableau de travail dans le tableau d'indices
!
                  IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
                  EXIT
              ENDIF
!
              IWRK = IWRK + 1
!
! .... monotonie 'a'  et 'b' ne sont pas epuisees
!
              IF (IINDA < JINDA) THEN
                 IF (IINDB < IWRKF) THEN
                    IVALA = IRNGT (IINDA + 1)
                    IVALB = IRNGT (IINDB + 1)
                    IF (XVALT (IVALA) > XVALT (IVALB)) THEN
                       IINDB = IINDB + 1
                       JWRKT (IWRK) = IRNGT (IINDB)
                    ELSE
                       IINDA = IINDA + 1
                       JWRKT (IWRK) = IRNGT (IINDA)
                    ENDIF
                 ELSE
!
! .... monotonie 'b' epuisee
!
                    IINDA = IINDA + 1
                    JWRKT (IWRK) = IRNGT (IINDA)
                 ENDIF
              ELSE
!
! .... monotonie 'a' epuisee
!
                 IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
                 IWRK = IWRKF
                 EXIT
              ENDIF
!
           ENDDO
        ENDDO
!
! ... on double les monotonies de depart et on recommence
!
        LMTNA = 2 * LMTNA
    ENDDO
!
! ... fin
!
    DEALLOCATE (JWRKT)
    RETURN
END SUBROUTINE TRIIND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION FRCTIL (XVALT, IRNGT, XORD)
!   frctil = fractile d'ordre XORD de l'ensemble XVALT dont
!            les rangs sont donnes par IRNGT 
REAL, DIMENSION (:), INTENT (IN)     :: XVALT 
INTEGER, DIMENSION (:), INTENT (IN)  :: IRNGT
REAL, INTENT (IN)                    :: XORD 
REAL                                 :: FRCTIL 
! __________________________________________________________
!
! ... nombre de valeurs rangees au total
!
    NVAL = SIZE (IRNGT)
    IF (NVAL <= 0 .OR. NVAL > SIZE (XVALT)) THEN
        WRITE (*, *) 'Dimensions incorrectes'
        FRCTIL = 0.0
        RETURN
    ENDIF
    IF (XORD <= 0.0 .OR. XORD >= 1.0) THEN
        WRITE (*, *) 'Ordre incorrect ( ]0.0, 1.0[ )'
        FRCTIL = 0.0
        RETURN
    ENDIF
!
! ... encadrement de la valeur
!
    IINF = 1 + FLOOR (XORD * REAL (NVAL - 1))
    ISUP = 1 + CEILING (XORD * REAL (NVAL - 1))
    ITST = MIN ((IINF - 1), (NVAL - ISUP))
    IF (ITST == 0) THEN
        WRITE (*, *) 'Significativite douteuse'
    ENDIF
!
    XVAL1 = XVALT (IRNGT (IINF))
    XVAL2 = XVALT (IRNGT (ISUP))
    IF (XVAL1 == XVAL2) THEN
        FRCTIL = XVAL1
    ELSE
!
! ... interpolation
!
        XORD1 = REAL (IINF - 1) / REAL (NVAL - 1)
        XORD2 = REAL (ISUP - 1) / REAL (NVAL - 1)
        FRCTIL = ((XORD - XORD1) * XVAL2 +  &
                  (XORD2 - XORD) * XVAL1) / &
                       (XORD2 - XORD1)
    ENDIF
    RETURN
END FUNCTION FRCTIL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALMOD (XVALT, IRNGT, XMOD, NMOD)
!   calcul de la valeur modale de l'ensemble XVALT dont
!            les rangs sont donnes par IRNGT 
REAL, DIMENSION (:), INTENT (IN)     :: XVALT 
INTEGER, DIMENSION (:), INTENT (IN)  :: IRNGT
REAL, INTENT (OUT)                   :: XMOD ! le mode si NMOD=1
INTEGER, INTENT (OUT)                :: NMOD ! nbre de modes 
! __________________________________________________________
    INTEGER, PARAMETER   :: NITRM = 16 ! maximum iterations
    REAL, DIMENSION (1:SIZE(IRNGT))     :: XWRKT
    REAL, DIMENSION (1:(SIZE(IRNGT)+2)) :: XDIST
    REAL, DIMENSION (1:SIZE(IRNGT))     :: XCCVT
    LOGICAL, DIMENSION (1:SIZE(IRNGT)-3):: IFMODT
!    INTERFACE
!        SUBROUTINE LISOUB (XVALT, FREN)
!  Lisse XVALT par moyenne courante a oubli progressif 
!        REAL, DIMENSION (:), INTENT (INOUT) :: XVALT 
!        REAL, INTENT (IN)                   :: FREN  
!                              Facteur de renouvellement 
!        END SUBROUTINE LISOUB
!        SUBROUTINE LISFEN (XVALT)
!  Lisse XVALT par moyenne courante a poids 1 - 2 - 1 
!        REAL, DIMENSION (:), INTENT (INOUT) :: XVALT 
!        END SUBROUTINE LISFEN
!    END INTERFACE
!
! ... nombre de valeurs rangees au total
!
    NVAL = SIZE (IRNGT)
    IF (NVAL <= 0 .OR. NVAL > SIZE (XVALT)) THEN
        WRITE (*, *) 'Dimensions incorrectes'
        NMOD = 0
        RETURN
    ENDIF
!
! ... estimation de la fonction de repartition
!
    XDIST (2:NVAL+1) = XVALT (IRNGT (1:NVAL))
    XDIST (1)  = XDIST (2)      &
               - (XDIST (NVAL+1) - XDIST (2)) / REAL (NVAL)
    XDIST (NVAL+2) = XDIST (NVAL+1) &
               + (XDIST (NVAL+1) - XDIST (2)) / REAL (NVAL)
    CALL LISFEN (XDIST (2:NVAL+1))
!
! ... concavite
!
    XCCVT (1:NVAL) = 2.0 * XDIST (2:NVAL+1) &
                         - XDIST (1:NVAL)   &
                         - XDIST (3:NVAL+2)
    DO JPUI = 1, NITRM
        FOUB = 2.0 ** (-JPUI)
        CALL LISOUB (XCCVT, FOUB)
        IFMODT(1:NVAL-3) =  XCCVT (2:NVAL-2)  > 0.0 .AND. &
                            XCCVT (3:NVAL-1) <= 0.0
        NMOD = COUNT (IFMODT (1:NVAL-3))
        XWRKT (1:NMOD) = PACK (XDIST (3:NVAL-1), IFMODT (:))
        SELECT CASE (NMOD)
        CASE (1)
            XMOD = XWRKT (1)
            EXIT
        CASE (0)
            XMOD = XDIST (2)
            EXIT
        CASE DEFAULT
            XMOD = XWRKT (1)
            CYCLE
        END SELECT
    ENDDO
    RETURN
END SUBROUTINE CALMOD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE LISOUB (XVALT, FREN)
!  Lisse XVALT par moyenne courante a oubli progressif 
REAL, DIMENSION (:), INTENT (INOUT) :: XVALT 
REAL, INTENT (IN)                   :: FREN  
!                              Facteur de renouvellement 
! __________________________________________________________
REAL, DIMENSION (SIZE (XVALT)) :: XWRKT
!
! ... nombre de valeurs au total
!
    NVAL = SIZE (XVALT)
    XWRKT (:) = 0.0
!
! ... initialisations
!
    FVIE = 1.0 / (1.0 + FREN)
    FNOV = FREN / (1.0 + FREN)
    ICRS = 1
    IDCR = NVAL
    XCRS = XVALT (ICRS)
    XDCR = XVALT (IDCR)
!
! ... Faire progresser les deux moyennes en sens contraire
!
    DO
       XWRKT (ICRS) = XWRKT (ICRS) + XCRS
       XWRKT (IDCR) = XWRKT (IDCR) + XDCR
       IF (ICRS >= NVAL) EXIT
       ICRS = ICRS + 1
       IDCR = IDCR - 1
       XCRS = FVIE * XCRS + FNOV * XVALT (ICRS)
       XDCR = FVIE * XDCR + FNOV * XVALT (IDCR)
    ENDDO
!
! ... Recopier dans le tableau d'origine
!
    XVALT (:) = 0.5 * XWRKT (:)
    RETURN
END SUBROUTINE LISOUB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE LISFEN (XVALT)
!  Lisse XVALT par moyenne courante a poids 1 - 2 - 1 
REAL, DIMENSION (:), INTENT (INOUT) :: XVALT 
! __________________________________________________________
!
! ... nombre de valeurs au total
!
    NVAL = SIZE (XVALT)
!
! ... initialisations
!
    XVAL1 = 0.5 * (XVALT (1) + XVALT (2))
    XVALN = 0.5 * (XVALT (NVAL-1) + XVALT (NVAL))
!
! ... La grosse formule: (Xi-1 + 2 Xi + Xi+1) / 4
!
    XVALT (2:NVAL-1) = 0.5 * XVALT (2:NVAL-1) + &
                     0.25 * (XVALT (1:NVAL-2) + XVALT (3:NVAL))
    XVALT (1) = XVAL1
    XVALT (NVAL) = XVALN
    RETURN
END SUBROUTINE LISFEN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION REAL_MAX(XVALT)

REAL, DIMENSION (:), INTENT (IN)     :: XVALT 
REAL                   :: tmp_max
INTEGER :: tmp_i,tmp_n

tmp_n=size(XVALT)
tmp_max=XVALT(1)
DO tmp_i=2,tmp_n
 if(XVALT(tmp_i)>tmp_max) tmp_max=XVALT(tmp_i)
END DO

 REAL_MAX=tmp_max

end FUNCTION REAL_MAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION REAL_MIN(XVALT)

REAL, DIMENSION (:), INTENT (IN)     :: XVALT 
REAL                    :: tmp_min
INTEGER :: tmp_i,tmp_n

tmp_n=size(XVALT)
tmp_min=XVALT(1)
DO tmp_i=2,tmp_n
 if(XVALT(tmp_i)<tmp_min) tmp_min=XVALT(tmp_i)
END DO

 REAL_MIN=tmp_min

RETURN
end FUNCTION REAL_MIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module utils_stat
