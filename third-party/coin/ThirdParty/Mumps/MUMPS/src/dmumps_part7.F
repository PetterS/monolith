C
C  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
C
C
C  This version of MUMPS is provided to you free of charge. It is public
C  domain, based on public domain software developed during the Esprit IV
C  European project PARASOL (1996-1999). Since this first public domain
C  version in 1999, research and developments have been supported by the
C  following institutions: CERFACS, CNRS, ENS Lyon, INPT(ENSEEIHT)-IRIT,
C  INRIA, and University of Bordeaux.
C
C  The MUMPS team at the moment of releasing this version includes
C  Patrick Amestoy, Maurice Bremond, Alfredo Buttari, Abdou Guermouche,
C  Guillaume Joslin, Jean-Yves L'Excellent, Francois-Henry Rouet, Bora
C  Ucar and Clement Weisbecker.
C
C  We are also grateful to Emmanuel Agullo, Caroline Bousquet, Indranil
C  Chowdhury, Philippe Combes, Christophe Daniel, Iain Duff, Vincent Espirat,
C  Aurelia Fevre, Jacko Koster, Stephane Pralet, Chiara Puglisi, Gregoire
C  Richard, Tzvetomila Slavova, Miroslav Tuma and Christophe Voemel who
C  have been contributing to this project.
C
C  Up-to-date copies of the MUMPS package can be obtained
C  from the Web pages:
C  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
C
C
C   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
C
C
C  User documentation of any code that uses this software can
C  include this complete notice. You can acknowledge (using
C  references [1] and [2]) the contribution of this package
C  in any scientific publication dependent upon the use of the
C  package. You shall use reasonable endeavours to notify
C  the authors of the package of this publication.
C
C   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
C   A fully asynchronous multifrontal solver using distributed dynamic
C   scheduling, SIAM Journal of Matrix Analysis and Applications,
C   Vol 23, No 1, pp 15-41 (2001).
C
C   [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
C   S. Pralet, Hybrid scheduling for the parallel solution of linear
C   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
C
      SUBROUTINE DMUMPS_635(N,KEEP,ICNTL,MPG)
      IMPLICIT NONE
          INTEGER N, KEEP(500), ICNTL(40), MPG
          KEEP(19)=0
          RETURN
      END SUBROUTINE DMUMPS_635
      SUBROUTINE DMUMPS_634(ICNTL,KEEP,MPG,INFO)
      IMPLICIT NONE
      INTEGER KEEP(500), MPG, INFO(40), ICNTL(40)
      IF (KEEP(19).EQ.0.AND.KEEP(110).EQ.0) THEN
        IF (KEEP(111).NE.0) THEN
         INFO(1) = -37
         INFO(2) = 16  
         IF (KEEP(110).EQ.0) INFO(2) = 24 
          IF(MPG.GT.0) THEN
           WRITE( MPG,'(A)')
     &'** ERROR  : Null space computation requirement'
          WRITE( MPG,'(A)')
     &'** not consistent with factorization options'
         ENDIF
         GOTO 333
        ENDIF
      ENDIF
      IF (ICNTL(9).NE.1) THEN
        IF (KEEP(111).NE.0) THEN
         INFO(1) = -37
         INFO(2) = 9  
         IF (MPG.GT.0) THEN
          WRITE(MPG,'(A)')
     &'** ERROR  ICNTL(25) incompatible with '
          WRITE( MPG,'(A)')
     &'** option transposed system (ICNLT(9)=1) '
         ENDIF
        ENDIF
        GOTO 333
      ENDIF
 333  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_634
      SUBROUTINE DMUMPS_637(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) id
      NULLIFY(id%root%QR_TAU)
      RETURN
      END SUBROUTINE DMUMPS_637
      SUBROUTINE DMUMPS_636(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) id
      IF (associated(id%root%QR_TAU))  THEN
        DEALLOCATE(id%root%QR_TAU)
        NULLIFY(id%root%QR_TAU)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_636
      SUBROUTINE DMUMPS_146( MYID, root, N, IROOT,
     &           COMM, IW, LIW, IFREE,
     &           A, LA, PTRAST, PTLUST_S, PTRFAC,
     &           STEP, INFO, LDLT, QR,
     &           WK, LWK, KEEP,KEEP8,DKEEP)
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      TYPE ( DMUMPS_ROOT_STRUC ) :: root
      INTEGER N, IROOT, COMM, LIW, MYID, IFREE
      INTEGER(8) :: LA
      INTEGER(8) :: LWK
      DOUBLE PRECISION WK( LWK )
      INTEGER KEEP(500)
      DOUBLE PRECISION    DKEEP(30)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PTLUST_S(KEEP(28)), STEP(N), IW( LIW )
      INTEGER INFO( 2 ), LDLT, QR
      DOUBLE PRECISION A( LA )
      INTEGER IOLDPS
      INTEGER(8) :: IAPOS
      INTEGER LOCAL_M, LOCAL_N, LPIV, IERR, allocok
      INTEGER FWD_LOCAL_N_RHS, FWD_MTYPE 
      INCLUDE 'mumps_headers.h'
      EXTERNAL numroc
      INTEGER numroc
        IF ( .NOT. root%yes ) RETURN
        IF ( KEEP(60) .NE. 0 ) THEN
          IF ((LDLT == 1 .OR. LDLT == 2) .AND. KEEP(60) == 3 ) THEN
            CALL DMUMPS_320( WK, root%MBLOCK,
     &      root%MYROW, root%MYCOL, root%NPROW, root%NPCOL,
     &      root%SCHUR_POINTER(1),
     &      root%SCHUR_LLD, root%SCHUR_NLOC,
     &      root%TOT_ROOT_SIZE, MYID, COMM )
          ENDIF
        RETURN
        ENDIF
        IOLDPS  = PTLUST_S(STEP(IROOT))+KEEP(IXSZ)
        IAPOS   = PTRAST(STEP(IROOT))
        LOCAL_M = IW( IOLDPS + 2 )
        LOCAL_N = IW( IOLDPS + 1 )
        IAPOS = PTRFAC(IW ( IOLDPS + 4 ))
        IF ( LDLT.EQ.0 .OR. LDLT.EQ.2 .OR. QR.ne.0 ) THEN
         LPIV = LOCAL_M + root%MBLOCK
        ELSE
         LPIV = 1
        END IF
        IF (associated( root%IPIV )) DEALLOCATE(root%IPIV)
        root%LPIV = LPIV
        ALLOCATE( root%IPIV( LPIV ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = LPIV
          WRITE(*,*) MYID,': problem allocating IPIV(',LPIV,') in root'
          CALL MUMPS_ABORT()
        END IF
        CALL DESCINIT( root%DESCRIPTOR(1), root%TOT_ROOT_SIZE,
     &      root%TOT_ROOT_SIZE, root%MBLOCK, root%NBLOCK,
     &      0, 0, root%CNTXT_BLACS, LOCAL_M, IERR )
        IF ( LDLT.EQ.2 ) THEN
            IF(root%MBLOCK.NE.root%NBLOCK) THEN
              WRITE(*,*) ' Error: symmetrization only works for'
              WRITE(*,*) ' square block sizes, MBLOCK/NBLOCK=',
     &        root%MBLOCK, root%NBLOCK
              CALL MUMPS_ABORT()
            END IF
            IF ( LWK .LT. min(
     &           int(root%MBLOCK,8) * int(root%NBLOCK,8),
     &           int(root%TOT_ROOT_SIZE,8)* int(root%TOT_ROOT_SIZE,8 )
     &         )) THEN
               WRITE(*,*) 'Not enough workspace for symmetrization.'
               CALL MUMPS_ABORT()
            END IF
            CALL DMUMPS_320( WK, root%MBLOCK,
     &      root%MYROW, root%MYCOL, root%NPROW, root%NPCOL,
     &      A( IAPOS ), LOCAL_M, LOCAL_N,
     &      root%TOT_ROOT_SIZE, MYID, COMM )
        END IF
        IF (LDLT.EQ.0.OR.LDLT.EQ.2) THEN
          CALL pdgetrf( root%TOT_ROOT_SIZE, root%TOT_ROOT_SIZE,
     &      A( IAPOS ),
     &      1, 1, root%DESCRIPTOR(1), root%IPIV(1), IERR )
          IF ( IERR .GT. 0 ) THEN
              INFO(1)=-10
              INFO(2)=IERR-1
          END IF
        ELSE
          CALL pdpotrf('L',root%TOT_ROOT_SIZE,A(IAPOS),
     &      1,1,root%DESCRIPTOR(1),IERR)
            IF ( IERR .GT. 0 ) THEN
              INFO(1)=-40
              INFO(2)=IERR-1
            END IF
        END IF
        IF (KEEP(258).NE.0) THEN
          IF (root%MBLOCK.NE.root%NBLOCK) THEN
            write(*,*) "Internal error in DMUMPS_146:",
     &      "Block size different for rows and columns",
     &      root%MBLOCK, root%NBLOCK
            CALL MUMPS_ABORT() 
          ENDIF
          CALL DMUMPS_763(root%MBLOCK, root%IPIV(1),root%MYROW,
     &         root%MYCOL, root%NPROW, root%NPCOL, A(IAPOS), LOCAL_M,
     &         LOCAL_N, root%TOT_ROOT_SIZE, MYID, DKEEP(6), KEEP(259),
     &         LDLT)
        ENDIF
        IF (KEEP(252) .NE. 0) THEN
          FWD_LOCAL_N_RHS = numroc(KEEP(253), root%NBLOCK,
     &                      root%MYCOL, 0, root%NPCOL)
          FWD_LOCAL_N_RHS = max(1,FWD_LOCAL_N_RHS)
          FWD_MTYPE       = 1 
          CALL DMUMPS_768(
     &         root%TOT_ROOT_SIZE,
     &         KEEP(253),          
     &         FWD_MTYPE,
     &         A(IAPOS),
     &         root%DESCRIPTOR(1),
     &         LOCAL_M, LOCAL_N, FWD_LOCAL_N_RHS,
     &         root%IPIV(1), LPIV,
     &         root%RHS_ROOT(1,1), LDLT,
     &         root%MBLOCK, root%NBLOCK,
     &         root%CNTXT_BLACS, IERR)
        ENDIF
        RETURN
      END SUBROUTINE DMUMPS_146
      SUBROUTINE DMUMPS_556(
     &     N,PIV,FRERE,FILS,NFSIZ,IKEEP,
     &     NCST,KEEP,KEEP8,id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) :: id
      INTEGER N,NCST
      INTEGER PIV(N),FRERE(N),FILS(N),NFSIZ(N),IKEEP(N,3)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER I,P11,P1,P2,K1,K2,NLOCKED
      LOGICAL V1,V2
      NCST = 0
      NLOCKED = 0
      P11 = KEEP(93)
      DO I=KEEP(93)-1,1,-2
         P1 = PIV(I)
         P2 = PIV(I+1)
         K1 = IKEEP(P1,1)
         IF(K1 .GT. 0) THEN
            V1 = (abs(id%A(K1))*(id%ROWSCA(P1)**2).GE.1.0D-1)
         ELSE
            V1 = .FALSE.
         ENDIF
         K2 = IKEEP(P2,1)
         IF(K2 .GT. 0) THEN
            V2 = (abs(id%A(K2))*(id%ROWSCA(P2)**2).GE.1.0D-1)
         ELSE
            V2 = .FALSE.
         ENDIF
         IF(V1 .AND. V2) THEN
            PIV(P11) = P1
            P11 = P11 - 1
            PIV(P11) = P2
            P11 = P11 - 1
         ELSE IF(V1) THEN
            NCST = NCST+1
            FRERE(NCST) = P1
            NCST = NCST+1
            FRERE(NCST) = P2
         ELSE IF(V2) THEN
            NCST = NCST+1
            FRERE(NCST) = P2                
            NCST = NCST+1
            FRERE(NCST) = P1
         ELSE
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P1
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P2                   
         ENDIF
      ENDDO
      DO I=1,NLOCKED
         PIV(I) = FILS(I)
      ENDDO
      KEEP(94) = KEEP(94) + KEEP(93) - NLOCKED
      KEEP(93) = NLOCKED
      DO I=1,NCST
         NLOCKED = NLOCKED + 1
         PIV(NLOCKED) = FRERE(I)
      ENDDO
      DO I=1,KEEP(93)/2
         NFSIZ(I) = 0
      ENDDO
      DO I=(KEEP(93)/2)+1,(KEEP(93)/2)+NCST,2
         NFSIZ(I) = I+1
         NFSIZ(I+1) = -1
      ENDDO
      DO I=(KEEP(93)/2)+NCST+1,(KEEP(93)/2)+KEEP(94)
         NFSIZ(I) = 0
      ENDDO
      END SUBROUTINE DMUMPS_556
      SUBROUTINE DMUMPS_550(N,NCMP,N11,N22,PIV,
     &     INVPERM,PERM)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER, intent(in) :: PIV(N),PERM(N)
      INTEGER, intent (out):: INVPERM(N)
      INTEGER CMP_POS,EXP_POS,I,J,N2,K
      N2 = N22/2
      EXP_POS = 1
      DO CMP_POS=1,NCMP
         J = PERM(CMP_POS)
         IF(J .LE. N2) THEN
            K = 2*J-1
            I = PIV(K)
            INVPERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
            K = K+1
            I = PIV(K)
            INVPERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ELSE
            K = N2 + J
            I = PIV(K)
            INVPERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ENDIF
      ENDDO
      DO K=N22+N11+1,N
         I = PIV(K)
         INVPERM(I) = EXP_POS
         EXP_POS = EXP_POS+1
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_550
      SUBROUTINE DMUMPS_547(
     &     N,NZ, IRN, ICN, PIV,
     &     NCMP, IW, LW, IPE, LEN, IQ, 
     &     FLAG, ICMP, IWFR,
     &     IERROR, KEEP,KEEP8, ICNTL)
      IMPLICIT NONE
      INTEGER N,NZ,NCMP,LW,IWFR,IERROR
      INTEGER ICNTL(40),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IRN(NZ),ICN(NZ),IW(LW),PIV(N),IPE(N+1)
      INTEGER LEN(N),IQ(N),FLAG(N),ICMP(N) 
      INTEGER MP,N11,N22,NDUP
      INTEGER I,K,J,N1,LAST,K1,K2,L
      MP = ICNTL(2)
      IERROR = 0
      N22 = KEEP(93)
      N11 = KEEP(94)
      NCMP = N22/2 + N11
      DO I=1,NCMP
         IPE(I) = 0
      ENDDO
      K = 1
      DO I=1,N22/2
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
      ENDDO
      K = N22/2 + 1
      DO I=N22+1,N22+N11
         J = PIV(I)
         ICMP(J) = K
         K = K + 1
      ENDDO
      DO I=N11+N22+1,N
         J = PIV(I)
         ICMP(J) = 0
      ENDDO
      DO K=1,NZ
         I = IRN(K)
         J = ICN(K)
         I = ICMP(I)
         J = ICMP(J)
         IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     &        .OR.(J.LT.1)) THEN
            IERROR = IERROR + 1
         ELSE
            IF (I.NE.J) THEN
               IPE(I) = IPE(I) + 1
               IPE(J) = IPE(J) + 1
            ENDIF
         ENDIF
      ENDDO
      IQ(1) = 1
      N1 = NCMP - 1
      IF (N1.GT.0) THEN
         DO I=1,N1
            IQ(I+1) = IPE(I) + IQ(I)
         ENDDO 
      ENDIF
      LAST = max(IPE(NCMP)+IQ(NCMP)-1,IQ(NCMP))
      DO I = 1,NCMP
         FLAG(I) = 0
         IPE(I)  = IQ(I)
      ENDDO
      DO K=1,LAST
        IW(K) = 0
      ENDDO
      IWFR = LAST + 1
      DO K=1,NZ
         I = IRN(K)
         J = ICN(K)
         I = ICMP(I)
         J = ICMP(J)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
         ENDIF
      ENDDO
      NDUP = 0
      DO I=1,NCMP
         K1 = IPE(I) 
         K2 = IQ(I) -1
         IF (K1.GT.K2) THEN
            LEN(I) = 0
            IQ(I)  = 0
         ELSE
            DO K=K1,K2
               J     = -IW(K)
               IF (J.LE.0) GO TO 250
               L     = IQ(J) 
               IQ(J) = L + 1
               IF (FLAG(J).EQ.I) THEN
                  NDUP = NDUP + 1
                  IW(L) = 0
                  IW(K) = 0
               ELSE
                  IW(L)   = I
                  IW(K)   = J
                  FLAG(J) = I
               ENDIF
            ENDDO
 250        IQ(I) = IQ(I) - IPE(I)
            IF (NDUP.EQ.0) LEN(I) = IQ(I)
         ENDIF
      ENDDO
      IF (NDUP.NE.0) THEN
         IWFR = 1
         DO I=1,NCMP
            K1 = IPE(I) 
            IF (IQ(I).EQ.0) THEN
               LEN(I) = 0
               IPE(I) = IWFR
               CYCLE
            ENDIF
            K2 = K1 + IQ(I) - 1
            L = IWFR
            IPE(I) = IWFR
            DO K=K1,K2
               IF (IW(K).NE.0) THEN
                  IW(IWFR) = IW(K)
                  IWFR     = IWFR + 1
               ENDIF
            ENDDO
            LEN(I) = IWFR - L 
         ENDDO
      ENDIF
      IPE(NCMP+1) = IPE(NCMP) + LEN(NCMP)
      IWFR = IPE(NCMP+1)
      RETURN
      END SUBROUTINE DMUMPS_547
      SUBROUTINE DMUMPS_551(
     &     N, NE, IP, IRN, SCALING,LSC,CPERM, DIAG,
     &     ICNTL, WEIGHT,MARKED,FLAG,
     &     PIV_OUT, INFO)
      IMPLICIT NONE
      INTEGER N, NE, ICNTL(10), INFO(10),LSC
      INTEGER CPERM(N),PIV_OUT(N),IP(N+1),IRN(NE), DIAG(N)
      DOUBLE PRECISION SCALING(LSC),WEIGHT(N+2)
      INTEGER MARKED(N),FLAG(N)
      INTEGER NUM1,NUM2,NUMTOT,PATH_LENGTH,NLAST
      INTEGER I,CUR_EL,CUR_EL_PATH,CUR_EL_PATH_NEXT,BEST_BEG
      INTEGER L1,L2,PTR_SET1,PTR_SET2,TUP,T22
      DOUBLE PRECISION BEST_SCORE,CUR_VAL,TMP,VAL
      DOUBLE PRECISION INITSCORE, DMUMPS_739, 
     &     DMUMPS_740, DMUMPS_741
      LOGICAL VRAI,FAUX,MAX_CARD_DIAG,USE_SCALING
      INTEGER SUM
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (SUM = 1, VRAI = .TRUE., FAUX = .FALSE.)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      MAX_CARD_DIAG = .TRUE.
      NUM1 = 0
      NUM2 = 0
      NUMTOT = 0
      NLAST = N
      INFO = 0      
      MARKED = 1
      FLAG = 0
      VAL = ONE
      IF(LSC .GT. 1) THEN
         USE_SCALING = .TRUE.
      ELSE
         USE_SCALING = .FALSE.
      ENDIF
      TUP = ICNTL(2)
      IF(TUP .EQ. SUM) THEN
        INITSCORE = ZERO
      ELSE
        INITSCORE = ONE
      ENDIF
      IF(ICNTL(2) .GT. 2 .OR. ICNTL(2) .LE. 0) THEN
         WRITE(*,*)
     &        'ERROR: WRONG VALUE FOR ICNTL(2) = ',ICNTL(2)
         INFO(1) = -1
         RETURN
      ENDIF
      T22 = ICNTL(1)
      IF(ICNTL(1) .LT. 0 .OR. ICNTL(1) .GT. 2) THEN
         WRITE(*,*)
     &        'ERROR: WRONG VALUE FOR ICNTL(1) = ',ICNTL(1)
         INFO(1) = -1
         RETURN
      ENDIF
      DO CUR_EL=1,N
         IF(MARKED(CUR_EL) .LE. 0) THEN
            CYCLE
         ENDIF
         IF(CPERM(CUR_EL) .LT. 0) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         PATH_LENGTH = 2
         CUR_EL_PATH = CPERM(CUR_EL)
         IF(CUR_EL_PATH .EQ. CUR_EL) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         MARKED(CUR_EL) = 0
         WEIGHT(1) = INITSCORE
         WEIGHT(2) = INITSCORE
         L1 = IP(CUR_EL+1)-IP(CUR_EL)
         L2 = IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH)
         PTR_SET1 = IP(CUR_EL)
         PTR_SET2 = IP(CUR_EL_PATH)
         IF(USE_SCALING) THEN
            VAL = -SCALING(CUR_EL_PATH) - SCALING(CUR_EL+N)
         ENDIF
         CUR_VAL = DMUMPS_741(
     &        CUR_EL,CUR_EL_PATH,
     &        IRN(PTR_SET1),IRN(PTR_SET2),
     &        L1,L2,
     &        VAL,DIAG,N,FLAG,FAUX,T22)
         WEIGHT(PATH_LENGTH+1) = 
     &        DMUMPS_739(WEIGHT(1),CUR_VAL,TUP)
         DO
            IF(CUR_EL_PATH .EQ. CUR_EL) EXIT
            PATH_LENGTH = PATH_LENGTH+1
            MARKED(CUR_EL_PATH) = 0
            CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
            L1 = IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH)
            L2 = IP(CUR_EL_PATH_NEXT+1)-IP(CUR_EL_PATH_NEXT)
            PTR_SET1 = IP(CUR_EL_PATH)
            PTR_SET2 = IP(CUR_EL_PATH_NEXT)
            IF(USE_SCALING) THEN
               VAL = -SCALING(CUR_EL_PATH_NEXT) 
     &              - SCALING(CUR_EL_PATH+N)
            ENDIF
            CUR_VAL = DMUMPS_741(
     &           CUR_EL_PATH,CUR_EL_PATH_NEXT,
     &           IRN(PTR_SET1),IRN(PTR_SET2),
     &           L1,L2,
     &           VAL,DIAG,N,FLAG,VRAI,T22)
            WEIGHT(PATH_LENGTH+1) = 
     &           DMUMPS_739(WEIGHT(PATH_LENGTH-1),CUR_VAL,TUP)
            CUR_EL_PATH = CUR_EL_PATH_NEXT
         ENDDO
         IF(mod(PATH_LENGTH,2) .EQ. 1) THEN
            IF(WEIGHT(PATH_LENGTH+1) .GE. WEIGHT(PATH_LENGTH)) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
            ELSE
               CUR_EL_PATH = CUR_EL
            ENDIF
            DO I=1,(PATH_LENGTH-1)/2
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 1
         ELSE
            IF(MAX_CARD_DIAG) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
               IF(DIAG(CUR_EL) .NE. 0) THEN 
                  BEST_BEG = CUR_EL_PATH
                  GOTO 1000
               ENDIF
               DO I=1,(PATH_LENGTH/2)
                  CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
                  IF(DIAG(CUR_EL_PATH) .NE. 0) THEN 
                     BEST_BEG = CUR_EL_PATH_NEXT
                     GOTO 1000
                  ENDIF
               ENDDO
            ENDIF
            BEST_BEG = CUR_EL
            BEST_SCORE = WEIGHT(PATH_LENGTH-1)
            CUR_EL_PATH = CPERM(CUR_EL)
            DO I=1,(PATH_LENGTH/2)-1
               TMP = DMUMPS_739(WEIGHT(PATH_LENGTH),
     &              WEIGHT(2*I-1),TUP)
               TMP = DMUMPS_740(TMP,WEIGHT(2*I),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               TMP = DMUMPS_739(WEIGHT(PATH_LENGTH+1),
     &              WEIGHT(2*I),TUP)
               TMP = DMUMPS_740(TMP,WEIGHT(2*I+1),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
 1000       CUR_EL_PATH = BEST_BEG
            DO I=1,(PATH_LENGTH/2)-1
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 2
            MARKED(CUR_EL_PATH) = -1
         ENDIF
      ENDDO
      DO I=1,N
         IF(MARKED(I) .LT. 0) THEN
            IF(DIAG(I) .EQ. 0) THEN
               PIV_OUT(NLAST) = I
               NLAST = NLAST - 1
            ELSE
               NUM1 = NUM1 + 1
               PIV_OUT(NUM2+NUM1) = I
               NUMTOT = NUMTOT + 1
            ENDIF
         ENDIF
      ENDDO
      INFO(2) = NUMTOT
      INFO(3) = NUM1
      INFO(4) = NUM2
      RETURN
      END SUBROUTINE DMUMPS_551
      FUNCTION DMUMPS_739(A,B,T)
      IMPLICIT NONE
      DOUBLE PRECISION DMUMPS_739
      DOUBLE PRECISION A,B
      INTEGER T
      INTEGER SUM
      PARAMETER(SUM = 1)
      IF(T .EQ. SUM) THEN
         DMUMPS_739 = A+B
      ELSE
         DMUMPS_739 = A*B
      ENDIF
      END FUNCTION DMUMPS_739
      FUNCTION DMUMPS_740(A,B,T)
      IMPLICIT NONE
      DOUBLE PRECISION DMUMPS_740
      DOUBLE PRECISION A,B
      INTEGER T
      INTEGER SUM
      PARAMETER(SUM = 1)
      IF(T .EQ. SUM) THEN
         DMUMPS_740 = A-B
      ELSE
         DMUMPS_740 = A/B
      ENDIF
      END FUNCTION DMUMPS_740
      FUNCTION DMUMPS_741(CUR_EL,CUR_EL_PATH,
     &     SET1,SET2,L1,L2,VAL,DIAG,N,FLAG,FLAGON,T)
      IMPLICIT NONE
      DOUBLE PRECISION DMUMPS_741
      INTEGER CUR_EL,CUR_EL_PATH,L1,L2,N
      INTEGER SET1(L1),SET2(L2),DIAG(N),FLAG(N)
      DOUBLE PRECISION VAL
      LOGICAL FLAGON
      INTEGER T
      INTEGER I,INTER,MERGE
      INTEGER STRUCT,MA47
      PARAMETER(STRUCT=0,MA47=1)
      IF(T .EQ. STRUCT) THEN
         IF(.NOT. FLAGON) THEN
            DO I=1,L1
               FLAG(SET1(I)) = CUR_EL
            ENDDO            
         ENDIF
         INTER = 0
         DO I=1,L2
            IF(FLAG(SET2(I)) .EQ. CUR_EL) THEN
               INTER = INTER + 1
               FLAG(SET2(I)) = CUR_EL_PATH
            ENDIF
         ENDDO
         MERGE = L1 + L2 - INTER
         DMUMPS_741 = dble(INTER) / dble(MERGE)
      ELSE IF (T .EQ. MA47) THEN
         MERGE = 3
         IF(DIAG(CUR_EL) .NE. 0) MERGE = 2
         IF(DIAG(CUR_EL_PATH) .NE. 0) MERGE = MERGE - 2
         IF(MERGE .EQ. 0) THEN
            DMUMPS_741 = dble(L1+L2-2)
            DMUMPS_741 = -(DMUMPS_741**2)/2.0D0
         ELSE IF(MERGE .EQ. 1) THEN
            DMUMPS_741 = - dble(L1+L2-4) * dble(L1-2)
         ELSE IF(MERGE .EQ. 2) THEN
            DMUMPS_741 = - dble(L1+L2-4) * dble(L2-2)
         ELSE
            DMUMPS_741 = - dble(L1-2) * dble(L2-2)
         ENDIF
      ELSE
         DMUMPS_741 = VAL
      ENDIF
      RETURN
      END FUNCTION 
      SUBROUTINE DMUMPS_622(NA, NCMP,
     &      INVPERM,PERM, 
     &      LISTVAR_SCHUR, SIZE_SCHUR, AOTOA)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN):: NA, NCMP
      INTEGER, INTENT(IN):: AOTOA(NCMP), PERM(NCMP)
      INTEGER, INTENT(OUT):: INVPERM(NA) 
      INTEGER CMP_POS, IO, I, K, IPOS
      DO CMP_POS=1, NCMP
        IO              = PERM(CMP_POS)
        INVPERM(AOTOA(IO)) = CMP_POS
      ENDDO
      IPOS = NCMP
      DO K =1,  SIZE_SCHUR
        I       = LISTVAR_SCHUR(K)
        IPOS    = IPOS+1
        INVPERM(I) = IPOS
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_622
      SUBROUTINE DMUMPS_623
     & (NA,N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     & IQ, FLAG, IWFR,
     & NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     & symmetry, SYM, MedDens, NBQD, AvgDens, 
     & LISTVAR_SCHUR, SIZE_SCHUR, ATOAO, AOTOA)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: NA,N,NZ,LW
      INTEGER, INTENT(IN)  :: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN)  :: IRN(NZ), ICN(NZ) 
      INTEGER, INTENT(IN)  :: ICNTL(40), SYM
      INTEGER, INTENT(INOUT) :: IFLAG
      INTEGER, INTENT(OUT) :: IERROR,NRORM,NIORM,IWFR
      INTEGER, INTENT(OUT) :: AOTOA(N)
      INTEGER, INTENT(OUT) :: ATOAO(NA)  
      INTEGER, INTENT(OUT) :: LEN(N), IPE(N+1)
      INTEGER, INTENT(OUT) :: symmetry, 
     &                        MedDens, NBQD, AvgDens
      INTEGER, INTENT(OUT)  :: FLAG(N), IW(LW), IQ(N)
      INTEGER MP, MPG
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR, THRESH, IAO
      INTEGER NZOFFA, NDIAGA
      DOUBLE PRECISION RSYM
      INTRINSIC nint
      ATOAO(1:NA) = 0
      DO I = 1, SIZE_SCHUR
        ATOAO(LISTVAR_SCHUR(I)) = -1
      ENDDO
      IAO = 0  
      DO I= 1, NA
        IF (ATOAO(I).LT.0) CYCLE
        IAO = IAO +1   
        ATOAO(I)   = IAO
        AOTOA(IAO) = I
      ENDDO
      MP = ICNTL(2)
      MPG= ICNTL(3)
      NIORM  = 3*N
      NDIAGA = 0
      IERROR = 0
      IPE(1:N+1) = 0
      DO K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ELSE
           NDIAGA = NDIAGA + 1
          ENDIF
        ENDIF
      ENDDO
      NZOFFA  = NIORM - 3*N
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (mod(IFLAG,2).EQ.0) IFLAG  = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (mod(K,10).GT.3 .OR. mod(K,10).EQ.0 .OR.
     &             (10.LE.K .AND. K.LE.20)) THEN
                 WRITE (MP,'(I8,A,I8,A,I8,A)')
     &             K,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (mod(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     &             K,'st entry (in row',I,' and column',J,') ignored'
                 IF (mod(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     &             K,'nd entry (in row',I,' and column',J,') ignored'
                 IF (mod(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     &             K,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2*N
      IQ(1) = 1
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO 110 I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
  110   CONTINUE
      ENDIF
      LAST = max(IPE(N)+IQ(N)-1,IQ(N))
      FLAG(1:N) = 0
      IPE(1:N)  = IQ(1:N)
      IW(1:LAST) = 0
      IWFR = LAST + 1
      DO 200 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                          .OR.(J.LT.1)) CYCLE
         I = ATOAO(I)
         J = ATOAO(J)
         IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
          ELSE
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
          ENDIF
         ENDIF
  200 CONTINUE
      NDUP = 0
      DO 260 I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) -1
        IF (K1.GT.K2) THEN
         LEN(I) = 0
         IQ(I)  = 0
        ELSE
         DO 240 K=K1,K2
           J     = -IW(K)
           IF (J.LE.0) GO TO 250
           L     = IQ(J) 
           IQ(J) = L + 1
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1
            IW(L) = 0
            IW(K) = 0
           ELSE
            IW(L)   = I
            IW(K)   = J
            FLAG(J) = I
           ENDIF
  240    CONTINUE
  250    IQ(I) = IQ(I) - IPE(I)
         IF (NDUP.EQ.0) LEN(I) = IQ(I)
        ENDIF
  260 CONTINUE
      IF (NDUP.NE.0) THEN
       IWFR = 1
       DO 280 I=1,N
         IF (IQ(I).EQ.0) THEN
             LEN(I) = 0
            IPE(I) = IWFR
            GOTO 280
         ENDIF
         K1 = IPE(I) 
         K2 = K1 + IQ(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K=K1,K2
           IF (IW(K).NE.0) THEN
            IW(IWFR) = IW(K)
            IWFR     = IWFR + 1
           ENDIF
  270    CONTINUE
         LEN(I) = IWFR - L 
  280  CONTINUE
      ENDIF
      IPE(N+1) = IPE(N) + LEN(N)
      IWFR = IPE(N+1)
      IF (SYM.EQ.0) THEN
      RSYM =  dble(NDIAGA+2*NZOFFA - (IWFR-1))/
     &            dble(NZOFFA+NDIAGA) 
      symmetry = nint (100.0D0*RSYM)
         IF (MPG .GT. 0)
     &  write(MPG,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
      ELSE
       symmetry = 100
      ENDIF
      AvgDens = nint(dble(IWFR-1)/dble(N))
      THRESH  = AvgDens*50 - AvgDens/10 + 1
      NBQD    = 0
      IF (N.GT.2) THEN
        IQ(1:N) = 0
        DO I= 1, N
          K = max(LEN(I),1)
          IQ(K) = IQ(K) + 1
          IF (K.GT.THRESH) NBQD = NBQD+1
        ENDDO
        K = 0
        MedDens = 0
        DO WHILE (K .LT. (N/2))
         MedDens = MedDens + 1
         K       = K+IQ(MedDens)
        ENDDO
      ELSE
        MedDens = AvgDens
      ENDIF
         IF (MPG .GT. 0)
     &  write(MPG,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE DMUMPS_623
      SUBROUTINE DMUMPS_549(N,PE,INVPERM,NFILS,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),INVPERM(N),NFILS(N),WORK(N)
      INTEGER I,FATHER,STKLEN,STKPOS,PERMPOS,CURVAR
      NFILS = 0
      DO I=1,N
         FATHER = -PE(I)
         IF(FATHER .NE. 0) NFILS(FATHER) = NFILS(FATHER) + 1
      ENDDO
      STKLEN = 0
      PERMPOS = 1
      DO I=1,N
         IF(NFILS(I) .EQ. 0) THEN
            STKLEN = STKLEN + 1
            WORK(STKLEN) = I
            INVPERM(I) = PERMPOS
            PERMPOS = PERMPOS + 1
         ENDIF
      ENDDO
      DO STKPOS = 1,STKLEN
         CURVAR = WORK(STKPOS)
         FATHER = -PE(CURVAR)
         DO
            IF(FATHER .EQ. 0) EXIT
            IF(NFILS(FATHER) .EQ. 1) THEN
               INVPERM(FATHER) = PERMPOS
               FATHER = -PE(FATHER)
               PERMPOS = PERMPOS + 1
            ELSE
               NFILS(FATHER) = NFILS(FATHER) - 1
               EXIT
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_549
      SUBROUTINE DMUMPS_548(N,PE,NV,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),NV(N),WORK(N)
      INTEGER I,FATHER,LEN,NEWSON,NEWFATHER
      DO I=1,N
         IF(NV(I) .GT. 0) CYCLE
         LEN = 1
         WORK(LEN) = I
         FATHER = -PE(I)
         DO
            IF(NV(FATHER) .GT. 0) THEN
               NEWSON = FATHER
               EXIT
            ENDIF
            LEN = LEN + 1
            WORK(LEN) = FATHER
            NV(FATHER) = 1
            FATHER = -PE(FATHER)
         ENDDO
         NEWFATHER = -PE(FATHER)
         PE(WORK(LEN)) = -NEWFATHER
         PE(NEWSON) = -WORK(1)
      ENDDO      
      END SUBROUTINE DMUMPS_548
