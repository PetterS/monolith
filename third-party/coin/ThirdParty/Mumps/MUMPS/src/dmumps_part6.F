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
      SUBROUTINE DMUMPS_324(A, LDA, NPIV, NBROW, K50 )
      IMPLICIT NONE
      INTEGER LDA, NPIV, NBROW, K50
      DOUBLE PRECISION A(int(LDA,8)*int(NBROW+NPIV,8))
      INTEGER(8) :: IOLD, INEW, J8
      INTEGER I , ILAST
      INTEGER NBROW_L_RECTANGLE_TO_MOVE
      IF ((NPIV.EQ.0).OR.(LDA.EQ.NPIV)) GOTO 500
      IF ( K50.NE.0 ) THEN
        IOLD = int(LDA  + 1,8)
        INEW = int(NPIV + 1,8)
        IF (IOLD .EQ. INEW ) THEN
          INEW = INEW + int(NPIV,8) * int(NPIV - 1,8)
          IOLD = IOLD + int(LDA,8) * int(NPIV - 1,8)
        ELSE
          DO I = 1, NPIV - 1
            IF ( I .LE. NPIV-2 ) THEN
              ILAST = I+1
            ELSE
              ILAST = I
            ENDIF
            DO J8 = 0_8, int(ILAST,8)
              A( INEW + J8 ) = A( IOLD + J8 )
            END DO
            INEW = INEW + int(NPIV,8)
            IOLD = IOLD + int(LDA,8)
          END DO
        ENDIF
        NBROW_L_RECTANGLE_TO_MOVE = NBROW
      ELSE 
        INEW = 1_8 + int(NPIV,8) * int(LDA + 1,8)
        IOLD = 1_8 + int(LDA,8) * int(NPIV +1,8)
        NBROW_L_RECTANGLE_TO_MOVE = NBROW - 1
      ENDIF
      DO I = 1, NBROW_L_RECTANGLE_TO_MOVE
         DO J8 = 0_8, int(NPIV - 1,8)
           A( INEW + J8 ) = A( IOLD + J8 )
         END DO
         INEW = INEW + int(NPIV,8)
         IOLD = IOLD + int(LDA,8)
      ENDDO
 500  RETURN
      END SUBROUTINE DMUMPS_324
      SUBROUTINE DMUMPS_651(A, LDA, NPIV, NCONTIG )
      IMPLICIT NONE
      INTEGER NCONTIG, NPIV, LDA
      DOUBLE PRECISION A(NCONTIG*LDA)
      INTEGER I, J
      INTEGER(8) :: INEW, IOLD
      INEW = int(NPIV+1,8)
      IOLD = int(LDA+1,8)
      DO I = 2, NCONTIG
        DO J = 1, NPIV
          A(INEW)=A(IOLD)
          INEW = INEW + 1_8
          IOLD = IOLD + 1_8
        ENDDO
        IOLD = IOLD + int(LDA - NPIV,8)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_651
      SUBROUTINE DMUMPS_652( A, LA, LDA, POSELT,
     &           IPTRLU, NPIV,
     &           NBCOL_STACK, NBROW_STACK,
     &           NBROW_SEND, SIZECB, KEEP, COMPRESSCB,
     &           LAST_ALLOWED, NBROW_ALREADY_STACKED )
      IMPLICIT NONE
      INTEGER(8), intent (in) :: POSELT, IPTRLU, LA, SIZECB
      LOGICAL, intent (in) :: COMPRESSCB
      DOUBLE PRECISION A(LA)
      INTEGER, intent(in):: LDA, NPIV, NBCOL_STACK, NBROW_STACK,
     &                      NBROW_SEND
      INTEGER, intent(inout) :: NBROW_ALREADY_STACKED
      INTEGER(8), intent(in)    :: LAST_ALLOWED
      INTEGER(8) :: APOS, NPOS
      INTEGER NBROW
      INTEGER(8) :: J
      INTEGER I, KEEP(500)
#if ! defined(ALLOW_NON_INIT)
      DOUBLE PRECISION ZERO
        PARAMETER( ZERO = 0.0D0 )
#endif
      NBROW = NBROW_STACK + NBROW_SEND
      IF (NBROW_STACK .NE. 0 ) THEN
        NPOS = IPTRLU + SIZECB         
        APOS = POSELT + int(NPIV+NBROW,8) * int(LDA,8) - 1_8 
        IF ( KEEP(50) .EQ. 0 .OR. .NOT. COMPRESSCB ) THEN
          APOS = APOS - int(LDA,8) * int(NBROW_ALREADY_STACKED,8)
          NPOS = NPOS
     &         - int(NBCOL_STACK,8) * int(NBROW_ALREADY_STACKED,8)
        ELSE
          APOS = APOS - int(LDA - 1,8) * int(NBROW_ALREADY_STACKED,8)
          NPOS = NPOS - ( int(NBROW_ALREADY_STACKED,8) *
     &                    int(NBROW_ALREADY_STACKED+1,8) ) / 2_8
        ENDIF
        DO I = NBROW - NBROW_ALREADY_STACKED, NBROW_SEND+1, -1
          IF (KEEP(50).EQ.0) THEN
            IF ( NPOS - int(NBCOL_STACK,8) + 1_8 .LT.
     &                                  LAST_ALLOWED ) THEN
              EXIT
            ENDIF
            DO J= 1_8,int(NBCOL_STACK,8)
              A(NPOS-J+1_8) = A(APOS-J+1_8)
            ENDDO
            NPOS = NPOS - int(NBCOL_STACK,8)
          ELSE
            IF (.NOT. COMPRESSCB) THEN
              IF ( NPOS - int(NBCOL_STACK,8) + 1_8 .LT.
     &                                  LAST_ALLOWED ) THEN
                EXIT
              ENDIF
#if ! defined(ALLOW_NON_INIT)
              DO J = 1_8, int(NBCOL_STACK - I,8)
                A(NPOS - J + 1_8) = ZERO
              END DO
#endif
              NPOS = NPOS + int(- NBCOL_STACK + I,8)
            ENDIF
            IF ( NPOS - int(I,8) + 1_8 .LT. LAST_ALLOWED ) THEN
              EXIT
            ENDIF
            DO J =1_8, int(I,8)
              A(NPOS-J+1_8) = A(APOS-J+1_8)
            ENDDO
            NPOS = NPOS - int(I,8)
          ENDIF
          IF (KEEP(50).EQ.0) THEN
            APOS = APOS - int(LDA,8)
          ELSE
            APOS = APOS - int(LDA + 1,8)
          ENDIF
          NBROW_ALREADY_STACKED = NBROW_ALREADY_STACKED + 1
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_652
      SUBROUTINE DMUMPS_705( A, LA, LDA, POSELT,
     &           IPTRLU, NPIV,
     &           NBCOL_STACK, NBROW_STACK,
     &           NBROW_SEND, SIZECB, KEEP, COMPRESSCB)
      IMPLICIT NONE
      INTEGER(8), intent (in) :: POSELT, IPTRLU, LA, SIZECB
      LOGICAL, intent (in) :: COMPRESSCB
      DOUBLE PRECISION A(LA)
      INTEGER, intent(in):: LDA, NPIV, NBCOL_STACK, NBROW_STACK,
     &                      NBROW_SEND
      INTEGER(8) :: APOS, NPOS, APOS_ini, NPOS_ini
      INTEGER I, KEEP(500)
      INTEGER(8) :: J, LDA8
#if ! defined(ALLOW_NON_INIT)
      DOUBLE PRECISION ZERO
        PARAMETER( ZERO = 0.0D0 )
#endif
      LDA8 = int(LDA,8)
      NPOS_ini = IPTRLU + 1_8
      APOS_ini = POSELT + int(NPIV+NBROW_SEND,8)* LDA8 + int(NPIV,8)
      DO I = 1, NBROW_STACK
         IF (COMPRESSCB) THEN
            NPOS = NPOS_ini + int(I-1,8) * int(I,8)/2_8 +
     &             int(I-1,8) * int(NBROW_SEND,8)
         ELSE
            NPOS = NPOS_ini + int(I-1,8) * int(NBCOL_STACK,8)
        ENDIF
        APOS  =  APOS_ini + int(I-1,8) * LDA8
        IF (KEEP(50).EQ.0) THEN
          DO J = 1_8, int(NBCOL_STACK,8)
            A(NPOS+J-1_8) = A(APOS+J-1_8)
          ENDDO
        ELSE
          DO J  = 1_8, int(I + NBROW_SEND,8)
            A(NPOS+J-1_8)=A(APOS+J-1_8)
          ENDDO
#if ! defined(ALLOW_NON_INIT)
          IF (.NOT. COMPRESSCB) THEN
            A(NPOS+int(I+NBROW_SEND,8):
     &        NPOS+int(NBCOL_STACK-1,8))=ZERO
          ENDIF
#endif
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_705
      SUBROUTINE DMUMPS_140( N, INODE, IW, LIW, A, LA,
     &                           IOLDPS, POSELT, IFLAG,
     &                           UU, NNEG, NPVW,
     &                           KEEP,KEEP8,
     &                           MYID, SEUIL, AVOID_DELAYED, ETATASS,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, IWPOS )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, NNEG, NPVW
      INTEGER MYID, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A( LA )
      INTEGER, TARGET :: IW( LIW )
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER INOPV, IFINB, NFRONT, NPIV, NBOLKJ,
     &        NBTLKJ,IBEG_BLOCK
      INTEGER NASS, NEL1, IFLAG_OOC
      INTEGER :: LDA
      DOUBLE PRECISION UUTEMP
      INCLUDE 'mumps_headers.h'
      EXTERNAL DMUMPS_222, DMUMPS_234, 
     &         DMUMPS_230, DMUMPS_226, 
     &         DMUMPS_237
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      INTEGER PIVSIZ,IWPOSP2
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, NextPiv2beWritten, 
     &        IDUMMY
      LOGICAL POSTPONE_COL_UPDATE, IS_MAXFROMM_AVAIL
      DOUBLE PRECISION MAXFROMM
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER PP_FIRST2SWAP_L
      INTEGER PP_LastPIVRPTRFilled  
      IS_MAXFROMM_AVAIL = .FALSE.
      INOPV = 0
      SEUIL_LOC = SEUIL
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
        UUTEMP=UU
      ENDIF
      POSTPONE_COL_UPDATE = (UUTEMP == 0.0D0 .AND. KEEP(201).NE.1)
      IBEG_BLOCK = 1
      NFRONT = IW(IOLDPS+KEEP(IXSZ))
      LDA    = NFRONT
      NASS   = iabs(IW(IOLDPS+2+KEEP(IXSZ)))
      IF (NASS .GT. KEEP(3)) THEN
        NBOLKJ = min( KEEP(6), NASS )
      ELSE
        NBOLKJ = min( KEEP(5), NASS )
      ENDIF
      NBTLKJ = NBOLKJ
      IF (KEEP(201).EQ.1) THEN 
          IDUMMY    = -8765
          CALL MUMPS_729(LAFAC,IW(IOLDPS+XXR))
          LIWFAC    = IW(IOLDPS+XXI)
          TYPEFile  = TYPEF_L
          NextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = NextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          PP_LastPIVRPTRFilled       = 0 
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 1
          MonBloc%NROW     = NFRONT
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -77777    
          MonBloc%INDICES  =>
     &              IW(IOLDPS+6+NFRONT+KEEP(IXSZ):
     &                 IOLDPS+5+NFRONT+KEEP(IXSZ)+NFRONT)
      ENDIF
      IW(IOLDPS+3+KEEP(IXSZ)) = min0(NASS,NBTLKJ)
      UUTEMP = UU
 50   CONTINUE
      CALL DMUMPS_222(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &                INOPV, NNEG, IFLAG,IOLDPS,POSELT,UUTEMP,
     &                SEUIL_LOC,KEEP,KEEP8,PIVSIZ,
     &     DKEEP(1),PIVNUL_LIST(1),LPN_LIST, KEEP(IXSZ),
     &     PP_FIRST2SWAP_L, MonBloc%LastPanelWritten_L,
     &     PP_LastPIVRPTRFilled, MAXFROMM, IS_MAXFROMM_AVAIL)
      IF (IFLAG.LT.0) GOTO 500
      IF(KEEP(109).GT. 0) THEN
         IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
            IWPOSP2 = IOLDPS+IW(IOLDPS+1+KEEP(IXSZ))+6+KEEP(IXSZ)
     &              +IW(IOLDPS+5+KEEP(IXSZ))
            PIVNUL_LIST(KEEP(109)) = IW(IWPOSP2)
         ENDIF
      ENDIF
      IF (INOPV.EQ.1) THEN
         IF(STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
         CALL DMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &        LDA, IOLDPS,POSELT, KEEP,KEEP8, POSTPONE_COL_UPDATE,
     &        ETATASS,
     &         TYPEFile, LAFAC, MonBloc, NextPiv2beWritten,
     &         LIWFAC, MYID, IFLAG)
         GOTO 500
      END IF
      IF (INOPV.EQ.2) THEN
         CALL DMUMPS_234(IBEG_BLOCK,
     &            NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &            LDA, IOLDPS,POSELT, NBOLKJ, NBTLKJ,KEEP(4),
     &            POSTPONE_COL_UPDATE,
     &            KEEP,KEEP8)
         GOTO 50
      ENDIF
      NPVW = NPVW + PIVSIZ
      IF (NASS.LE.1) THEN
       CALL DMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     &                 IOLDPS,POSELT)
       IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+1+KEEP(IXSZ)) + 1
       GO TO 500
      ENDIF
       CALL DMUMPS_226(IBEG_BLOCK,
     &             NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     &             LDA, POSTPONE_COL_UPDATE, IOLDPS,
     &             POSELT,IFINB,
     &             NBTLKJ,PIVSIZ, KEEP(IXSZ),MAXFROMM,
     &             IS_MAXFROMM_AVAIL, (UUTEMP.NE.0.0D0),
     &             KEEP(253) )
       IF(PIVSIZ .EQ. 2) THEN
          IWPOSP2 = IOLDPS+IW(IOLDPS+1+KEEP(IXSZ))+6
          IW(IWPOSP2+NFRONT+KEEP(IXSZ)) = -IW(IWPOSP2+NFRONT+KEEP(IXSZ))
       ENDIF
       IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+1+KEEP(IXSZ)) + PIVSIZ
       IF (IFINB.EQ.0) GOTO 50
       NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
       NEL1   = NASS - NPIV
        IF (KEEP(201).EQ.1) THEN
          IF (IFINB.EQ.-1) THEN
            MonBloc%Last = .TRUE.
          ELSE
            MonBloc%Last   = .FALSE.
          ENDIF
          MonBloc%LastPiv= NPIV
          LAST_CALL=.FALSE.
          CALL DMUMPS_688(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
          IF (IFLAG .LT. 0 ) RETURN
        ENDIF
       CALL DMUMPS_234(IBEG_BLOCK,
     &            NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &            LDA, IOLDPS,POSELT, NBOLKJ, NBTLKJ,KEEP(4),
     &            POSTPONE_COL_UPDATE,
     &            KEEP,KEEP8)
       IF (IFINB.EQ.-1) THEN 
         CALL DMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &         LDA, IOLDPS,POSELT, KEEP,KEEP8,
     &         POSTPONE_COL_UPDATE, ETATASS,
     &         TYPEFile, LAFAC, MonBloc, NextPiv2beWritten,
     &         LIWFAC, MYID, IFLAG)
     &
         GOTO 500
       ENDIF
      GO TO 50
 500  CONTINUE
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+KEEP(IXSZ))
          LAST_CALL=.TRUE.
          CALL DMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           NextPiv2beWritten, IDUMMY,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
          IF (IFLAG < 0 ) RETURN
          CALL DMUMPS_644 (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_140
      SUBROUTINE DMUMPS_222 
     &   (NFRONT,NASS,N,INODE,IW,LIW,
     &    A,LA, INOPV,
     &    NNEG,
     &    IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,PIVSIZ,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, XSIZE,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk,
     &     PP_LastPIVRPTRIndexFilled,MAXFROMM,IS_MAXFROMM_AVAIL)
#if defined (PROFILE_BLAS_ASS_G)
      USE DMUMPS_LOAD
#endif
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER NFRONT,NASS,N,LIW,INODE,IFLAG,INOPV,
     &        IOLDPS, NNEG
      INTEGER PIVSIZ,LPIV, XSIZE
      DOUBLE PRECISION A(LA) 
      DOUBLE PRECISION UU, UULOC, SEUIL
      INTEGER IW(LIW)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk
      INTEGER PP_LastPIVRPTRIndexFilled
      DOUBLE PRECISION, intent(in) :: MAXFROMM
      LOGICAL, intent(inout) :: IS_MAXFROMM_AVAIL
      include 'mpif.h'
      INTEGER (8) :: POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX
      DOUBLE PRECISION RMAX,AMAX,TMAX,TOL
      DOUBLE PRECISION MAXPIV
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      DOUBLE PRECISION PIVOT,DETPIV
      PARAMETER(TOL = 1.0D-20)
      INCLUDE 'mumps_headers.h'
      INTEGER :: J
      INTEGER(8) :: APOS, J1, J2, JJ, NFRONT8, KK, J1_ini, JJ_ini
      INTEGER    :: LDA
      INTEGER(8) :: LDA8
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,K
      INTRINSIC max
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO = 0.0D0 )
      PARAMETER( ONE = 1.0D0 )
      DOUBLE PRECISION RZERO,RONE
      PARAMETER(RZERO=0.0D0, RONE=1.0D0)
      LOGICAL OMP_FLAG
      INTEGER I_PIVRPTR, I_PIVR, NBPANELS_L
      PIVNUL = DKEEP(1)
      FIXA   = DKEEP(2)
      CSEUIL = SEUIL
      LDA     = NFRONT
      LDA8    = int(LDA,8)
      NFRONT8 = int(NFRONT,8)
      IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
             CALL DMUMPS_667(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR, I_PIVR, IOLDPS+2*NFRONT+6+KEEP(IXSZ),
     &       IW, LIW)
      ENDIF
      UULOC = UU
      PIVSIZ = 1
      NPIV    = IW(IOLDPS+1+XSIZE)
      NPIVP1  = NPIV + 1
      NASSW   = iabs(IW(IOLDPS+3+XSIZE))
      IF(INOPV .EQ. -1) THEN
         APOS = POSELT + (LDA8+1_8) * int(NPIV,8)
         POSPV1 = APOS
         IF(abs(A(APOS)).LT.SEUIL) THEN
            IF(dble(A(APOS)) .GE. RZERO) THEN
               A(APOS) = CSEUIL
            ELSE
               A(APOS) = -CSEUIL
               NNEG = NNEG+1
            ENDIF
            KEEP(98) = KEEP(98)+1
         ELSE IF (KEEP(258) .NE. 0) THEN
            CALL DMUMPS_762( A(APOS), DKEEP(6), KEEP(259) )
         ENDIF
              IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
                CALL DMUMPS_680( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, NPIVP1,
     &               PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
         GO TO 420
      ENDIF
      INOPV   = 0
      DO 460 IPIV=NPIVP1,NASSW
         APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIV,8)
         POSPV1 = APOS + int(IPIV - NPIVP1,8)
         PIVOT = A(POSPV1)
         IF (UULOC.EQ.RZERO) THEN 
            IF (abs(A(APOS)).EQ.RZERO) GO TO 630
            IF (A(APOS).LT.RZERO) NNEG = NNEG+1
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_762(A(APOS), DKEEP(6), KEEP(259))
            ENDIF
            GO TO 420
         ENDIF
         AMAX = RZERO
         JMAX = 0
         IF ( IS_MAXFROMM_AVAIL ) THEN
            IF ( MAXFROMM > PIVNUL .AND. abs(PIVOT) > TOL) THEN
               IF ( abs(PIVOT) .GT. max(UULOC*MAXFROMM,SEUIL) ) THEN
                 IF (PIVOT .LT. RZERO) NNEG = NNEG+1
                 IF (KEEP(258) .NE. 0) THEN
                   CALL DMUMPS_762(PIVOT, DKEEP(6), KEEP(259))
                 ENDIF
                 GOTO 415
               ENDIF
            ENDIF
            IS_MAXFROMM_AVAIL = .FALSE.
         ENDIF
         J1 = APOS
         J2 = POSPV1 - 1_8
         DO JJ=J1,J2
            IF(abs(A(JJ)) .GT. AMAX) THEN
               AMAX = abs(A(JJ))
               JMAX = IPIV - int(POSPV1-JJ)
            ENDIF
         ENDDO
         J1 = POSPV1 + LDA8
         DO J=1, NASSW - IPIV
            IF(abs(A(J1)) .GT. AMAX) THEN
               AMAX = abs(A(J1))
               JMAX = IPIV + J
            ENDIF
            J1 = J1 + LDA8
         ENDDO
           RMAX = RZERO
           J1_ini = J1
           IF ( (NFRONT - KEEP(253) - NASSW).GE.300 ) THEN
             OMP_FLAG = .TRUE.
           ELSE
             OMP_FLAG = .FALSE.
           ENDIF
           DO J=1, NFRONT - KEEP(253) - NASSW
              J1 = J1_ini + int(J-1,8) * LDA8
              RMAX = max(abs(A(J1)),RMAX)
           ENDDO
         IF (max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL) THEN
            KEEP(109) = KEEP(109)+1
            PIVNUL_LIST(KEEP(109)) = -1
            IF(dble(FIXA).GT.RZERO) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDA8
               DO J=1, NASSW - IPIV
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               DO J=1,NFRONT - NASSW
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               A(POSPV1) = ONE
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
         IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
           IF (max(AMAX,RMAX,abs(PIVOT)).LE.TOL) THEN
            IF(SEUIL .GT. epsilon(SEUIL)) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = CSEUIL
               ELSE
                  A(POSPV1) = -CSEUIL
                  NNEG = NNEG+1
               ENDIF
               PIVOT = A(POSPV1)
               KEEP(98) = KEEP(98)+1
               GO TO 415
            ENDIF
           ENDIF
         ENDIF
         IF (max(AMAX,abs(PIVOT)).LE.TOL) GO TO 460
         IF (abs(PIVOT).GT.max(UULOC*max(RMAX,AMAX),SEUIL)) THEN
               IF (PIVOT .LT. ZERO) NNEG = NNEG+1
               IF (KEEP(258) .NE.0 ) THEN
                 CALL DMUMPS_762(PIVOT, DKEEP(6), KEEP(259))
               ENDIF
               GO TO 415
         END IF
         IF (AMAX.LE.TOL) GO TO 460
         IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  IF(int(POSPV1-JJ) .NE. IPIV-JMAX) THEN
                     RMAX = max(RMAX,abs(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDA8
               DO J=1,NASS-IPIV
                  IF(IPIV+J .NE. JMAX) THEN
                     RMAX = max(abs(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDA8
               ENDDO
           ENDIF
           APOSJ = POSELT + int(JMAX-1,8)*LDA8 + int(NPIV,8)
           POSPV2 = APOSJ + int(JMAX - NPIVP1,8)
           IF (IPIV.LT.JMAX) THEN
              OFFDAG = APOSJ + int(IPIV - NPIVP1,8)
           ELSE
              OFFDAG = APOS + int(JMAX - NPIVP1,8)
           END IF
           TMAX = RZERO
           IF(JMAX .LT. IPIV) THEN
              JJ_ini = POSPV2
              OMP_FLAG = (NFRONT-JMAX-KEEP(253). GE. 300)
              DO K = 1, NFRONT - JMAX - KEEP(253)
                 JJ = JJ_ini+ int(K,8)*NFRONT8
                 IF (JMAX+K.NE.IPIV) THEN
                    TMAX=max(TMAX,abs(A(JJ)))
                 ENDIF
              ENDDO
              DO KK =  APOSJ, POSPV2-1_8
                 TMAX = max(TMAX,abs(A(KK)))
              ENDDO
           ELSE
              JJ_ini = POSPV2
              OMP_FLAG = (NFRONT-JMAX-KEEP(253). GE. 300)
              DO K = 1, NFRONT-JMAX-KEEP(253)
                 JJ = JJ_ini + int(K,8)*NFRONT8
                 TMAX=max(TMAX,abs(A(JJ)))
              ENDDO
              DO KK =  APOSJ, POSPV2 - 1_8
                 IF (KK.NE.OFFDAG) THEN
                    TMAX = max(TMAX,abs(A(KK)))
                 ENDIF
              ENDDO
           ENDIF
           DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
           IF (SEUIL.GT.RZERO) THEN
                IF (sqrt(abs(DETPIV)) .LE. SEUIL ) GOTO 460
           ENDIF
           MAXPIV = max(abs(A(POSPV1)),abs(A(POSPV2)))
           IF (MAXPIV.EQ.RZERO) MAXPIV = RONE
           IF (abs(DETPIV)/MAXPIV.LE.TOL) GO TO 460
           IF ((abs(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     &          abs(DETPIV)) GO TO 460
           IF ((abs(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     &          abs(DETPIV)) GO TO 460
           IF (KEEP(258) .NE.0 ) THEN
             CALL DMUMPS_762(DETPIV, DKEEP(6), KEEP(259))
           ENDIF
           PIVSIZ = 2
           KEEP(103) = KEEP(103)+1
           IF(DETPIV .LT. RZERO) THEN
             NNEG = NNEG+1
           ELSE IF(A(POSPV2) .LT. RZERO) THEN
             NNEG = NNEG+2
           ENDIF
 415       CONTINUE
           DO K=1,PIVSIZ
              IF (PIVSIZ .EQ. 2) THEN
                IF (K==1) THEN
                  LPIV = min(IPIV,JMAX)
                ELSE
                  LPIV   = max(IPIV,JMAX)
                ENDIF
              ELSE
                LPIV = IPIV
              ENDIF
              IF (LPIV.EQ.NPIVP1) THEN
                 GOTO 416
              ENDIF
              CALL DMUMPS_319( A, LA, IW, LIW,
     &             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     &             LDA, NFRONT, 1, KEEP(219), KEEP(50),
     &             KEEP(IXSZ)) 
 416          CONTINUE
              IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
                CALL DMUMPS_680( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, LPIV, PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
              NPIVP1 = NPIVP1 + 1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+(LDA8+1_8)*int(NPIV,8)+1_8) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      PIVSIZ = 0
      IFLAG = -10
  420 CONTINUE
      IS_MAXFROMM_AVAIL = .FALSE.
      RETURN
      END SUBROUTINE DMUMPS_222
      SUBROUTINE DMUMPS_680( PIVRPTR, NBPANELS, PIVR, NASS,
     &                                  K, P, LastPanelonDisk,
     &                                  LastPIVRPTRIndexFilled)
      IMPLICIT NONE
      INTEGER, intent(in) :: NBPANELS, NASS, K, P
      INTEGER, intent(inout) :: PIVRPTR(NBPANELS), PIVR(NASS)
      INTEGER LastPanelonDisk, LastPIVRPTRIndexFilled
      INTEGER I
      IF ( LastPanelonDisk+1 > NBPANELS ) THEN
           WRITE(*,*) "INTERNAL ERROR IN DMUMPS_680!"
           WRITE(*,*) "NASS=",NASS,"PIVRPTR=",PIVRPTR(1:NBPANELS)
           WRITE(*,*) "K=",K, "P=",P, "LastPanelonDisk=",LastPanelonDisk
           WRITE(*,*) "LastPIVRPTRIndexFilled=", LastPIVRPTRIndexFilled
           CALL MUMPS_ABORT()
      ENDIF
      PIVRPTR(LastPanelonDisk+1) = K + 1
      IF (LastPanelonDisk.NE.0) THEN
        PIVR(K - PIVRPTR(1) + 1) = P
        DO I = LastPIVRPTRIndexFilled + 1, LastPanelonDisk
          PIVRPTR(I)=PIVRPTR(LastPIVRPTRIndexFilled)
        ENDDO
      ENDIF
      LastPIVRPTRIndexFilled = LastPanelonDisk + 1
      RETURN
      END SUBROUTINE DMUMPS_680
      SUBROUTINE DMUMPS_226(IBEG_BLOCK,
     &     NFRONT,NASS,N,INODE,IW,LIW,
     &     A,LA,LDA, POSTPONE_COL_UPDATE,
     &     IOLDPS,POSELT,IFINB,LKJIB,PIVSIZ,XSIZE,
     &     MAXFROMM, IS_MAXFROMM_AVAIL, IS_MAX_USEFUL,
     &     KEEP253)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,IFINB,LKJIB,
     &        NPBEG, IBEG_BLOCK
      INTEGER LDA
      INTEGER(8) :: LA
      INTEGER(8) :: NFRONT8
      DOUBLE PRECISION    A(LA)
      LOGICAL POSTPONE_COL_UPDATE
      INTEGER IW(LIW)
      DOUBLE PRECISION    VALPIV
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, intent(out) :: MAXFROMM
      LOGICAL, intent(out) :: IS_MAXFROMM_AVAIL
      LOGICAL, intent(in) :: IS_MAX_USEFUL
      INTEGER, INTENT(in) :: KEEP253
      DOUBLE PRECISION :: MAXFROMMTMP
      INTEGER IOLDPS, NCB1
      INTEGER(8) :: LDA8
      INTEGER(8) :: K1POS
      INTEGER NPIV,JROW2
      INTEGER NEL2,NEL
      INTEGER XSIZE
      DOUBLE PRECISION ONE, ZERO
      INTEGER(8) :: APOS, LPOS, LPOS1, LPOS2
      INTEGER(8) :: POSPV1, POSPV2
      INTEGER PIVSIZ,NPIV_NEW,J2,I
      INTEGER(8) :: OFFDAG, OFFDAG_OLD, IBEG, IEND
      INTEGER(8) :: JJ, K1, K2, IROW
      DOUBLE PRECISION SWOP,DETPIV,MULT1,MULT2
      INCLUDE 'mumps_headers.h'
      PARAMETER(ONE  = 1.0D0,
     &          ZERO = 0.0D0)
      LDA8   = int(LDA,8)
      NFRONT8= int(NFRONT,8)
      NPIV   = IW(IOLDPS+1+XSIZE)
      NPIV_NEW = NPIV + PIVSIZ
      NEL    = NFRONT - NPIV_NEW
      IFINB  = 0
      IS_MAXFROMM_AVAIL = .FALSE.
      JROW2 = IW(IOLDPS+3+XSIZE)
      NPBEG = IBEG_BLOCK
      NEL2   = JROW2 - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         A(APOS) = VALPIV
         LPOS   = APOS + LDA8
         MAXFROMM = 0.0D00
         IF (NEL2 > 0) THEN
           IF (.NOT. IS_MAX_USEFUL) THEN
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               DO JJ=1_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ELSE
             IS_MAXFROMM_AVAIL = .TRUE.
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               A(K1POS+1_8)=A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
               MAXFROMM=max( MAXFROMM,abs(A(K1POS+1_8)) )
               DO JJ = 2_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ENDIF
         ENDIF
         IF (POSTPONE_COL_UPDATE) THEN
           NCB1 = NASS   - JROW2
         ELSE
           NCB1 = NFRONT - JROW2
         ENDIF
         IF (.NOT. IS_MAX_USEFUL) THEN
           DO I=NEL2+1, NEL2 + NCB1
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             DO JJ = 1_8, int(NEL2,8)
               A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
             ENDDO
           ENDDO
         ELSE
           MAXFROMMTMP=0.0D0
           DO I=NEL2+1, NEL2 + NCB1 - KEEP253
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             IF (NEL2 > 0) THEN
               A(K1POS+1_8) = A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
               MAXFROMMTMP=max(MAXFROMMTMP, abs(A(K1POS+1_8)))
               DO JJ = 2_8, int(NEL2,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDIF
           ENDDO
           DO I = NEL2 + NCB1 - KEEP253 + 1, NEL2 + NCB1
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             DO JJ = 1_8, int(NEL2,8)
               A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
             ENDDO
           ENDDO
           MAXFROMM=max(MAXFROMM, MAXFROMMTMP)
         ENDIF
      ELSE
         POSPV1 = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         POSPV2 = POSPV1 + NFRONT8 + 1_8
         OFFDAG_OLD = POSPV2 - 1_8
         OFFDAG = POSPV1 + 1_8
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
         A(POSPV2) = A(POSPV1)/DETPIV
         A(POSPV1) = SWOP/DETPIV
         A(OFFDAG) = -A(OFFDAG_OLD)/DETPIV
         A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDA8 - 1_8
         LPOS2   = LPOS1 + 1_8
         CALL dcopy(NFRONT-NPIV_NEW, A(LPOS1), LDA, A(POSPV1+2_8), 1)
         CALL dcopy(NFRONT-NPIV_NEW, A(LPOS2), LDA, A(POSPV2+1_8), 1)
         JJ = POSPV2 + NFRONT8-1_8
         IBEG = JJ + 2_8
         IEND = IBEG
         DO J2 = 1,NEL2
            K1 = JJ
            K2 = JJ+1_8
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1 + 2_8
            K2 = POSPV2 + 1_8
            DO IROW = IBEG, IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1_8
               K2 = K2 + 1_8
            ENDDO
            A( JJ       ) = -MULT1
            A( JJ + 1_8 ) = -MULT2
            IBEG = IBEG + NFRONT8
            IEND = IEND + NFRONT8 + 1_8
            JJ = JJ+NFRONT8
         ENDDO
         IEND = IEND-1_8
         DO J2 = JROW2+1,NFRONT
            K1 = JJ
            K2 = JJ+1_8
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1 + 2_8
            K2 = POSPV2 + 1_8
            DO IROW = IBEG, IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1_8
               K2 = K2 + 1_8
            ENDDO
            A( JJ       ) = -MULT1
            A( JJ + 1_8 ) = -MULT2
            IBEG = IBEG + NFRONT8
            IEND = IEND + NFRONT8
            JJ   = JJ   + NFRONT8
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_226
      SUBROUTINE DMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     &           IOLDPS,POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LIW
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    VALPIV
      INTEGER (8) :: APOS, POSELT, LPOS, NFRONT8
      INTEGER IOLDPS,NEL
      INTEGER JROW
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
        APOS   = POSELT 
        VALPIV = ONE/A(APOS)
        A(APOS) = VALPIV
        NEL    = NFRONT - 1
        IF (NEL.EQ.0) GO TO 500
        NFRONT8 = int(NFRONT,8)
        LPOS   = APOS + NFRONT8
        CALL DMUMPS_XSYR('U',NEL, -VALPIV, 
     &             A(LPOS), NFRONT, A(LPOS+1_8), NFRONT)
          DO JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT8
          END DO
  500   CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_230
      SUBROUTINE DMUMPS_234(IBEG_BLOCK,
     &    NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &    LDA,
     &    IOLDPS,POSELT,LKJIB_ORIG,LKJIB,LKJIT,
     &    POSTPONE_COL_UPDATE,
     &    KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,LIW, IBEG_BLOCK
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB_ORIG, LKJIB, INODE, KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: POSELT
      INTEGER LDA
      INTEGER(8) :: LDA8
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1, NEL11
      INTEGER LBP, HF
      INTEGER(8) :: LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER I, Block
      INTEGER BLSIZE
      LOGICAL POSTPONE_COL_UPDATE
      DOUBLE PRECISION ONE, ALPHA
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      LDA8 = int(LDA,8)
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
      JROW2  = iabs(IW(IOLDPS+3+KEEP(IXSZ)))
      NPBEG  = IBEG_BLOCK
      HF     = 6 + IW(IOLDPS+5+KEEP(IXSZ)) + KEEP(IXSZ)
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      NEL11  = NFRONT - NPIV
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB_ORIG + NONEL
          IW(IOLDPS+3+KEEP(IXSZ))= min0(NPIV+LKJIB,NASS)
          LKJIB       = min0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3+KEEP(IXSZ)) = NASS
        ENDIF
        IBEG_BLOCK = NPIV + 1
      ELSEIF (JROW2.LT.NASS) THEN
          IBEG_BLOCK   = NPIV + 1
          IW(IOLDPS+3+KEEP(IXSZ)) = min0(JROW2+LKJIB,NASS)
          LKJIB  = min0(LKJIB,NASS-NPIV)
      ENDIF
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
#if defined(SAK_BYROW)
         DO IROW = JROW2+1, NASS, BLSIZE
           Block = min( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + int(IROW  - 1,8) * LDA8 + int(NPBEG - 1,8)
           UPOS = POSELT + int(NPBEG - 1,8) * LDA8 + int(IROW - 1,8)
           APOS = POSELT + int(IROW  - 1,8) * LDA8 + int(JROW2,8)
           CALL dgemm( 'N','N', IROW + Block - JROW2 - 1, Block, LKJIW,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         ENDDO
#else
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = min( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + int( IROW - 1,8) * LDA8 + int(NPBEG - 1,8)
           UPOS = POSELT + int(NPBEG - 1,8) * LDA8 + int( IROW - 1,8)
           APOS = POSELT + int( IROW - 1,8) * LDA8 + int( IROW - 1,8)
           CALL dgemm( 'N','N', Block, NASS - IROW + 1, LKJIW,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         END DO
#endif
        END IF
       LPOS = POSELT + int(NASS,8)*LDA8 + int(NPBEG - 1,8)
       UPOS = POSELT + int(NPBEG-1,8) * LDA8 + int(JROW2,8)
       APOS = POSELT + int(NASS,8)*LDA8 + int(JROW2,8)
       IF ( .NOT. POSTPONE_COL_UPDATE ) THEN
         CALL dgemm('N', 'N', NEL1, NFRONT-NASS, LKJIW, ALPHA, 
     &              A(UPOS), LDA, A(LPOS), LDA, ONE, 
     &              A(APOS), LDA)
       END IF
      ENDIF
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_234
        SUBROUTINE DMUMPS_319( A, LA, IW, LIW,
     &                       IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     &                       LDA, NFRONT, LEVEL, K219, K50, XSIZE )
        IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER LIW, IOLDPS, NPIVP1, IPIV
      INTEGER LDA, NFRONT, NASS, LEVEL, K219, K50, XSIZE
      DOUBLE PRECISION A( LA )
      INTEGER IW( LIW )
      INCLUDE 'mumps_headers.h'
      INTEGER ISW, ISWPS1, ISWPS2, HF
      INTEGER(8) :: IDIAG, APOS
      INTEGER(8) :: LDA8
      DOUBLE PRECISION SWOP
            LDA8 = int(LDA,8)
            APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIVP1-1,8)
            IDIAG = APOS + int(IPIV - NPIVP1,8)
            HF = 6 + IW( IOLDPS + 5 + XSIZE) + XSIZE
            ISWPS1 = IOLDPS + HF + NPIVP1 - 1
            ISWPS2 = IOLDPS + HF + IPIV - 1
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            ISW = IW(ISWPS1+NFRONT)
            IW(ISWPS1+NFRONT) = IW(ISWPS2+NFRONT)
            IW(ISWPS2+NFRONT) = ISW
            IF ( LEVEL .eq. 2 ) THEN
              CALL dswap( NPIVP1 - 1,
     &            A( POSELT + int(NPIVP1-1,8) ), LDA,
     &            A( POSELT + int(IPIV-1,8)   ), LDA )
            END IF
            CALL dswap( NPIVP1-1,
     &           A( POSELT+int(NPIVP1-1,8) * LDA8 ), 1,
     &           A( POSELT + int(IPIV-1,8) * LDA8 ), 1 )
             CALL dswap( IPIV - NPIVP1 - 1,
     &           A( POSELT+int(NPIVP1,8) * LDA8 + int(NPIVP1-1,8) ),
     &           LDA, A( APOS + 1_8 ), 1 )
            SWOP = A(IDIAG)
            A(IDIAG) = A( POSELT+int(NPIVP1-1,8)*LDA8+int(NPIVP1-1,8) )
            A( POSELT + int(NPIVP1-1,8)*LDA8 + int(NPIVP1-1,8) ) = SWOP
            CALL dswap( NASS - IPIV, A( APOS + LDA8 ), LDA,
     &                  A( IDIAG + LDA8 ), LDA )
            IF ( LEVEL .eq. 1 ) THEN
              CALL dswap( NFRONT - NASS,
     &        A( APOS  + int(NASS-IPIV+1,8) * LDA8 ), LDA,
     &        A( IDIAG + int(NASS-IPIV+1,8) * LDA8 ), LDA )
            END IF
            IF (K219.NE.0 .AND.K50.EQ.2) THEN
             IF ( LEVEL .eq. 2) THEN
              APOS                 = POSELT+LDA8*LDA8-1_8
              SWOP                 = A(APOS+int(NPIVP1,8))
              A(APOS+int(NPIVP1,8))= A(APOS+int(IPIV,8))
              A(APOS+int(IPIV,8))  = SWOP
             ENDIF
            ENDIF
        RETURN
        END SUBROUTINE DMUMPS_319
      SUBROUTINE DMUMPS_237(NFRONT,NASS,N,INODE,
     &    IW,LIW,A,LA,
     &    LDA,
     &    IOLDPS,POSELT,KEEP,KEEP8,
     &    POSTPONE_COL_UPDATE, ETATASS,
     &    TYPEFile, LAFAC, MonBloc, NextPiv2beWritten,
     &    LIWFAC, MYID, IFLAG
     &    )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,INODE,LIW
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW) 
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: POSELT
      INTEGER LDA
      INTEGER IOLDPS, ETATASS
      LOGICAL POSTPONE_COL_UPDATE
      INTEGER(8) :: LAFAC
      INTEGER TYPEFile, NextPiv2beWritten
      INTEGER LIWFAC, MYID, IFLAG
      TYPE(IO_BLOCK):: MonBloc
      INTEGER IDUMMY
      LOGICAL LAST_CALL
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: UPOS, APOS, LPOS
      INTEGER(8) :: LDA8
      INTEGER BLSIZE, BLSIZE2, Block, IROW, NPIV, I, IROWEND
      INTEGER I2, I2END, Block2
      DOUBLE PRECISION  ONE, ALPHA, BETA, ZERO
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      LDA8 = int(LDA,8)
      IF (ETATASS.EQ.1) THEN
        BETA = ZERO
      ELSE
        BETA = ONE
      ENDIF
      IF ( NFRONT - NASS > KEEP(57) ) THEN
        BLSIZE = KEEP(58)
      ELSE
        BLSIZE = NFRONT - NASS
      END IF
      BLSIZE2 = KEEP(218)
      NPIV = IW( IOLDPS + 1 + KEEP(IXSZ))
      IF ( NFRONT - NASS .GT. 0 ) THEN
       IF ( POSTPONE_COL_UPDATE ) THEN
         CALL dtrsm( 'L', 'U', 'T', 'U',
     &               NPIV, NFRONT-NPIV, ONE,
     &               A( POSELT ), LDA,
     &               A( POSELT + LDA8 * int(NPIV,8) ), LDA )
       ENDIF
       DO IROWEND = NFRONT - NASS, 1, -BLSIZE
        Block = min( BLSIZE, IROWEND )
        IROW  = IROWEND - Block + 1
        LPOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8
        APOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8 +
     &                  int(NASS + IROW - 1,8)
        UPOS = POSELT + int(NASS,8)
        IF (.NOT. POSTPONE_COL_UPDATE) THEN
          UPOS = POSELT + int(NASS + IROW - 1,8)
        ENDIF
        IF (POSTPONE_COL_UPDATE) THEN
         DO I = 1, NPIV
          CALL dcopy( Block, A( LPOS+int(I-1,8) ), LDA,
     &                       A( UPOS+int(I-1,8)*LDA8 ), 1 )
          CALL dscal( Block, A(POSELT+(LDA8+1_8)*int(I-1,8)),
     &                A( LPOS + int(I - 1,8) ), LDA )
         ENDDO
        ENDIF
        DO I2END = Block, 1, -BLSIZE2
          Block2 = min(BLSIZE2, I2END)
          I2 = I2END - Block2+1
          CALL dgemm('N', 'N', Block2, Block-I2+1, NPIV, ALPHA,
     &               A(UPOS+int(I2-1,8)), LDA,
     &               A(LPOS+int(I2-1,8)*LDA8), LDA,
     &               BETA,
     &               A(APOS + int(I2-1,8) + int(I2-1,8)*LDA8), LDA)
          IF (KEEP(201).EQ.1) THEN
            IF (NextPiv2beWritten.LE.NPIV) THEN
              LAST_CALL=.FALSE.
              CALL DMUMPS_688(
     &        STRAT_TRY_WRITE, TYPEFile,
     &        A(POSELT), LAFAC, MonBloc,
     &        NextPiv2beWritten, IDUMMY,
     &        IW(IOLDPS), LIWFAC, MYID,
     &        KEEP8(31),
     &        IFLAG,LAST_CALL )
              IF (IFLAG .LT. 0 ) RETURN
            ENDIF
          ENDIF
        ENDDO
        IF ( NFRONT - NASS - IROW + 1 - Block > 0 ) THEN
        CALL dgemm( 'N', 'N', Block, NFRONT-NASS-Block-IROW+1, NPIV,
     &              ALPHA,  A( UPOS ), LDA,
     &              A( LPOS + LDA8 * int(Block,8) ), LDA,
     &              BETA,
     &              A( APOS + LDA8 * int(Block,8) ), LDA )
        ENDIF
       END DO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_237
      SUBROUTINE DMUMPS_320( BUF, BLOCK_SIZE,
     &                           MYROW, MYCOL, NPROW, NPCOL,
     &                           A, LOCAL_M, LOCAL_N, N, MYID, COMM )
      IMPLICIT NONE
      INTEGER BLOCK_SIZE, NPROW, NPCOL, LOCAL_M, LOCAL_N, N, COMM
      INTEGER MYROW, MYCOL, MYID
      DOUBLE PRECISION BUF( BLOCK_SIZE * BLOCK_SIZE )
      DOUBLE PRECISION A( LOCAL_M, LOCAL_N )
      INTEGER NBLOCK, IBLOCK, JBLOCK, IBLOCK_SIZE, JBLOCK_SIZE
      INTEGER ROW_SOURCE, ROW_DEST, COL_SOURCE, COL_DEST
      INTEGER IGLOB, JGLOB
      INTEGER IROW_LOC_SOURCE, JCOL_LOC_SOURCE
      INTEGER IROW_LOC_DEST, JCOL_LOC_DEST
      INTEGER PROC_SOURCE, PROC_DEST
      NBLOCK = ( N - 1 ) / BLOCK_SIZE + 1
      DO IBLOCK = 1, NBLOCK
        IF ( IBLOCK .NE. NBLOCK
     &    ) THEN
          IBLOCK_SIZE = BLOCK_SIZE
        ELSE
          IBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
        END IF
        ROW_SOURCE = mod( IBLOCK - 1, NPROW ) 
        COL_DEST   = mod( IBLOCK - 1, NPCOL )
        IGLOB = ( IBLOCK - 1 ) * BLOCK_SIZE + 1
        IROW_LOC_SOURCE = BLOCK_SIZE *
     &                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     &                  + mod( IGLOB - 1, BLOCK_SIZE ) + 1
        JCOL_LOC_DEST   = BLOCK_SIZE *
     &                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     &                  + mod( IGLOB - 1, BLOCK_SIZE ) + 1
        DO JBLOCK = 1, IBLOCK
          IF ( JBLOCK .NE. NBLOCK
     &      ) THEN
            JBLOCK_SIZE = BLOCK_SIZE
          ELSE
            JBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
          END IF
          COL_SOURCE = mod( JBLOCK - 1, NPCOL )
          ROW_DEST   = mod( JBLOCK - 1, NPROW )
          PROC_SOURCE = ROW_SOURCE * NPCOL + COL_SOURCE
          PROC_DEST   = ROW_DEST   * NPCOL + COL_DEST
          IF ( PROC_SOURCE .eq. PROC_DEST ) THEN
           IF ( MYID .eq. PROC_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     &                  ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     &                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            IF ( IBLOCK .eq. JBLOCK ) THEN
              IF ( IBLOCK_SIZE .ne. JBLOCK_SIZE ) THEN
                WRITE(*,*) MYID,': Error in calling transdiag:unsym'
                CALL MUMPS_ABORT()
              END IF
              CALL DMUMPS_327( A( IROW_LOC_SOURCE,
     &                 JCOL_LOC_SOURCE),
     &                 IBLOCK_SIZE, LOCAL_M )
            ELSE
              CALL DMUMPS_326(
     &           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ),
     &           A( IROW_LOC_DEST, JCOL_LOC_DEST ),
     &           IBLOCK_SIZE, JBLOCK_SIZE, LOCAL_M )
            END IF
           END IF
          ELSE IF (  MYROW .eq. ROW_SOURCE 
     &    .AND. MYCOL .eq. COL_SOURCE ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     &                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL DMUMPS_293( BUF,
     &           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ), LOCAL_M,
     &           IBLOCK_SIZE, JBLOCK_SIZE, COMM, PROC_DEST )
          ELSE IF ( MYROW .eq. ROW_DEST 
     &    .AND.     MYCOL .eq. COL_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     &                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL DMUMPS_281( BUF,
     &           A( IROW_LOC_DEST, JCOL_LOC_DEST ), LOCAL_M,
     &           JBLOCK_SIZE, IBLOCK_SIZE, COMM, PROC_SOURCE )
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_320
      SUBROUTINE DMUMPS_293( BUF, A, LDA, M, N, COMM, DEST )
      IMPLICIT NONE
      INTEGER M, N, LDA, DEST, COMM
      DOUBLE PRECISION BUF(*), A(LDA,*)
      INTEGER I, IBUF, IERR
      INTEGER J
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IBUF = 1
      DO J = 1, N
        BUF( IBUF: IBUF + M - 1 ) = A( 1 : M, J )
        DO I = 1, M
        END DO
        IBUF = IBUF + M
      END DO
      CALL MPI_SEND( BUF, M * N, MPI_DOUBLE_PRECISION,
     &     DEST, SYMMETRIZE, COMM, IERR )
      RETURN
      END SUBROUTINE DMUMPS_293
      SUBROUTINE DMUMPS_281( BUF, A, LDA, M, N, COMM, SOURCE )
      IMPLICIT NONE
      INTEGER LDA, M, N, COMM, SOURCE
      DOUBLE PRECISION BUF(*), A( LDA, *)
      INTEGER I, IBUF, IERR
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      CALL MPI_RECV( BUF(1), M * N, MPI_DOUBLE_PRECISION, SOURCE,
     &               SYMMETRIZE, COMM, STATUS, IERR )
      IBUF = 1
      DO I = 1, M
        CALL dcopy( N, BUF(IBUF), 1, A(I,1), LDA )
        IBUF = IBUF + N
      END DO
      RETURN
      END SUBROUTINE DMUMPS_281
      SUBROUTINE DMUMPS_327( A, N, LDA )
      IMPLICIT NONE
      INTEGER N,LDA
      DOUBLE PRECISION A( LDA, * )
      INTEGER I, J
      DO I = 2, N
        DO J = 1, I - 1
          A( J, I ) = A( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_327
      SUBROUTINE DMUMPS_326( A1, A2, M, N, LD )
      IMPLICIT NONE
      INTEGER M,N,LD
      DOUBLE PRECISION A1( LD,* ), A2( LD, * )
      INTEGER I, J
      DO J = 1, N
        DO I = 1, M
          A2( J, I ) = A1( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_326
      RECURSIVE SUBROUTINE DMUMPS_274( 
     &   COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS,  
     &    PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) PTRAST(KEEP(28)), PTRFAC(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), ND(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      DOUBLE PRECISION DBLARR( max(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER PIVI
      INTEGER (8) POSPV1,POSPV2,OFFDAG,LPOS1
      INTEGER J2
      DOUBLE PRECISION MULT1,MULT2
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LP
      INTEGER INODE, POSITION, NPIV, IERR
      INTEGER NCOL
      INTEGER(8) LAELL, POSBLOCFACTO
      INTEGER(8) POSELT
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11
      INTEGER(8) LPOS, LPOS2, DPOS, UPOS
      INTEGER (8) IPOS, KPOS
      INTEGER I, IPIV, FPERE, NSLAVES_TOT,
     &        NSLAVES_FOLLOW, NB_BLOC_FAC
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok, TO_UPDATE_CPT_END
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: UIP21K
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPivDummy
      LOGICAL LAST_CALL
      TYPE(IO_BLOCK) :: MonBloc
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      FPERE    = -1
      POSITION = 0
      TO_UPDATE_CPT_END = -654321
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     &                 MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LAELL = int(NPIV,8) * int(NCOL,8)
      IF ( NPIV.GT.0 ) THEN
       IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          CALL MUMPS_731(LAELL-LRLUS, IERROR)
          IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE IN DMUMPS_274,
     & REAL WORKSPACE TOO SMALL"
          GOTO 700
        END IF
        CALL DMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     &       LRLU, IPTRLU,
     &       IWPOS, IWPOSCB, PTRIST, PTRAST,
     &       STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &       KEEP(IXSZ))
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     &       ,LRLU,LRLUS
             IFLAG = -9
             CALL MUMPS_731(LAELL-LRLUS,IERROR)
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE IN DMUMPS_274,
     & INTEGER WORKSPACE TOO SMALL"
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
       END IF
       LRLU  = LRLU - LAELL
       LRLUS = LRLUS - LAELL
      ENDIF
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL DMUMPS_471(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8,LAELL,KEEP,KEEP8,LRLU)
      IF ( NPIV.GT.0 ) THEN
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              A(POSBLOCFACTO), NPIV*NCOL, MPI_DOUBLE_PRECISION,
     &              COMM, IERR )
      ENDIF
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
         DO WHILE ( PTRIST(STEP(INODE)) .EQ. 0 )
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &      BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &      MSGSOU, MAITRE_DESC_BANDE,
     &      STATUS,
     &      BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &      IFLAG, IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF,
     &      NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &      LPTRAR, NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
          IF ( IFLAG .LT. 0 ) GOTO 600
         END DO
      ENDIF
      DO WHILE ( NBPROCFILS(STEP(INODE)) .NE. 0 ) 
        BLOCKING = .TRUE.
        SET_IRECV=.FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, CONTRIB_TYPE2,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END  DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC, 
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS + KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 + KEEP(IXSZ))
      NROW1  = IW( IOLDPS + 2 + KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 + KEEP(IXSZ))
      NSLAV1 = IW( IOLDPS + 5 + KEEP(IXSZ))
      NSLAVES_FOLLOW = NSLAV1 - XTRA_SLAVES_SYM
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      IF ( LASTBL ) THEN
        TO_UPDATE_CPT_END = ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     &                       NB_BLOC_FAC
      END IF
      IF (NPIV.GT.0) THEN
        IF ( NPIV1 + NCOL .NE. NASS1 ) THEN
          WRITE(*,*) 'SymBLFC Error: NPIV1 + NCOL .NE. NASS1 :',
     &               NPIV1,NCOL,NASS1
          CALL MUMPS_ABORT()
        END IF
        ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
        DO I = 1, NPIV
          PIVI = abs(IW(IPIV+I-1))
          IF (PIVI.EQ.I) CYCLE
          ISW = IW(ICT11+I)
          IW(ICT11+I) = IW(ICT11+PIVI)
          IW(ICT11+PIVI) = ISW
          IPOS = POSELT + int(NPIV1 + I - 1,8)
          KPOS = POSELT + int(NPIV1 + PIVI - 1,8)
          CALL dswap(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
        ENDDO
        ALLOCATE( UIP21K( NPIV * NROW1 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR UIP21K IN DMUMPS_274"
          IFLAG = -13
          IERROR = NPIV * NROW1
          GOTO 700
        END IF
        IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
          ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     &            stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR LIST_SLAVES_FOLLOW
     & IN DMUMPS_274"
            IFLAG = -13
            IERROR = NSLAVES_FOLLOW
            GOTO 700
          END IF
          LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     &    IW(IOLDPS+6+XTRA_SLAVES_SYM+KEEP(IXSZ):
     &     IOLDPS+5+XTRA_SLAVES_SYM+KEEP(IXSZ)+NSLAVES_FOLLOW)
        END IF
        CALL dtrsm( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     &               A( POSBLOCFACTO ), NCOL,
     &               A(POSELT+int(NPIV1,8)), NCOL1 )
        LPOS = POSELT + int(NPIV1,8)
        UPOS = 1_8
        DO I = 1, NROW1
          UIP21K( UPOS: UPOS + int(NPIV-1,8) ) = 
     &                       A(LPOS: LPOS+int(NPIV-1,8))
          LPOS = LPOS + int(NCOL1,8)
          UPOS = UPOS + int(NPIV,8)
        END DO
        LPOS = POSELT + int(NPIV1,8)
        DPOS = POSBLOCFACTO
        I = 1
        DO
          IF(I .GT. NPIV) EXIT
          IF(IW(IPIV+I-1) .GT. 0) THEN
            CALL dscal( NROW1, A(DPOS), A(LPOS), NCOL1 )
            LPOS = LPOS + 1_8
            DPOS = DPOS + int(NCOL + 1,8)
            I = I+1
          ELSE
            POSPV1 = DPOS
            POSPV2 = DPOS+ int(NCOL + 1,8)
            OFFDAG = POSPV1+1_8
            LPOS1 = LPOS
            DO J2 = 1,NROW1
               MULT1 = A(POSPV1)*A(LPOS1)+A(OFFDAG)*A(LPOS1+1_8)
               MULT2 = A(OFFDAG)*A(LPOS1)+A(POSPV2)*A(LPOS1+1_8)
               A(LPOS1) = MULT1
               A(LPOS1+1_8) = MULT2
               LPOS1 = LPOS1 + int(NCOL1,8)
            ENDDO
            LPOS = LPOS + 2_8
            DPOS = POSPV2 + int(NCOL + 1,8)
            I = I+2
          ENDIF
        ENDDO
      ENDIF
      IF (KEEP(201).eq.1) THEN
        MonBloc%INODE = INODE
        MonBloc%MASTER = .FALSE.
        MonBloc%Typenode = 2
        MonBloc%NROW = NROW1  
        MonBloc%NCOL = NCOL1  
        MonBloc%NFS  = NASS1
        MonBloc%LastPiv = NPIV1 + NPIV 
        NULLIFY(MonBloc%INDICES)
        MonBloc%Last = LASTBL
        STRAT = STRAT_TRY_WRITE 
        NextPivDummy      = -8888 
        LIWFAC = IW(IOLDPS+XXI)
        CALL MUMPS_729(LAFAC, IW(IOLDPS+XXR))
        LAST_CALL=.FALSE.
        CALL DMUMPS_688( STRAT, TYPEF_L, A(POSELT),
     &       LAFAC, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
      ENDIF
      IF (NPIV.GT.0) THEN
        LPOS2 = POSELT + int(NPIV1,8)
        UPOS = POSBLOCFACTO+int(NPIV,8)
        LPOS  = LPOS2 + int(NPIV,8)
        CALL dgemm('N','N', NCOL-NPIV,NROW1,NPIV,ALPHA,A(UPOS),NCOL,
     &           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
        DPOS = POSELT + int(NCOL1 - NROW1,8)
        IF ( NROW1 .GT. KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NROW1
        ENDIF
        IF ( NROW1 .GT. 0 ) THEN
          DO IROW = 1, NROW1, BLSIZE
            Block = min( BLSIZE, NROW1 - IROW + 1 )
            DPOS  = POSELT + int(NCOL1 - NROW1,8)
     &            + int( IROW - 1, 8 ) * int( NCOL1 + 1, 8 )
            LPOS2 = POSELT + int(NPIV1,8)
     &            + int( IROW - 1, 8 ) * int( NCOL1, 8 )
            UPOS  = int( IROW - 1, 8 ) * int(NPIV, 8) + 1_8
            DO I = 1, Block
              CALL dgemv( 'T', NPIV, Block-I+1, ALPHA,
     &                A( LPOS2 + int(I - 1,8) * int(NCOL1,8) ), NCOL1,
     &                UIP21K( UPOS + int(NPIV,8) * int( I - 1, 8 ) ),
     &                1, ONE, A(DPOS+int(NCOL1+1,8)*int(I-1,8)),NCOL1 )
            END DO
           IF ( NROW1-IROW+1-Block .ne. 0 )
     &     CALL dgemm( 'T', 'N', Block, NROW1-IROW+1-Block, NPIV, ALPHA,
     &             UIP21K( UPOS ), NPIV,
     &             A( LPOS2 + int(Block,8) * int(NCOL1,8) ), NCOL1, ONE,
     &             A( DPOS + int(Block,8) * int(NCOL1,8) ), NCOL1 )
          ENDDO
        ENDIF
        FLOP1 = dble(NROW1) * dble(NPIV) *
     &           dble( 2 * NCOL  - NPIV + NROW1 +1 )
        FLOP1 = -FLOP1
        CALL DMUMPS_190( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF 
      IW(IOLDPS+KEEP(IXSZ)) = IW(IOLDPS+KEEP(IXSZ)) - NPIV
      IW(IOLDPS + 3+KEEP(IXSZ)) = IW(IOLDPS+3+KEEP(IXSZ)) + NPIV
      IF (LASTBL) IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS + 3+KEEP(IXSZ))
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      POSFAC = POSFAC - LAELL
      IWPOS = IWPOS - NPIV
      CALL DMUMPS_471(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLU)
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
          IERR = -1
           DO WHILE ( IERR .eq. -1 )
            CALL DMUMPS_64(
     &                    INODE, NPIVSENT, FPERE,
     &                    IPOSK, JPOSK,
     &                    UIP21K, NROW1,
     &                    NSLAVES_FOLLOW,
     &                    LIST_SLAVES_FOLLOW(1),
     &                    COMM, IERR )
            IF (IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV= .FALSE.
              MESSAGE_RECEIVED = .FALSE.
              CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS, 
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &         PTLUST_S, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, COMM,
     &         NBPROCFILS,
     &         IPOOL, LPOOL, LEAF,
     &         NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &         LPTRAR, NELT, FRTPTR, FRTELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
             IF ( IFLAG .LT. 0 ) GOTO 600
            END IF
           END DO
           IF ( IERR .eq. -2 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING
     & DMUMPS_274"
             WRITE(LP,*) "NPIV=", NPIV, "NROW1=",NROW1
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           IF ( IERR .eq. -3 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING
     & DMUMPS_274"
             IFLAG = -20
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           DEALLOCATE(LIST_SLAVES_FOLLOW)
      END IF
      IF ( NPIV .NE. 0 ) DEALLOCATE( UIP21K )
      IOLDPS = PTRIST(STEP(INODE))
      IF (LASTBL) THEN
         IW(IOLDPS+6+KEEP(IXSZ)) =  IW(IOLDPS+6+KEEP(IXSZ)) -
     &                            TO_UPDATE_CPT_END 
         IF ( IW(IOLDPS+6+KEEP(IXSZ) ) .eq. 0
     &        .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 0
     &        .and. NSLAVES_TOT.NE.1)THEN
         DEST = MUMPS_275( PROCNODE_STEPS(STEP(INODE)), SLAVEF )
         CALL DMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     &                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
           IFLAG = -99
           GOTO 700
         END IF
         ENDIF
      END IF
      IF (LASTBL) THEN 
        IF (IW(IOLDPS+6+KEEP(IXSZ)) .eq. 0 ) THEN 
          CALL DMUMPS_759( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root,
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_274
      RECURSIVE SUBROUTINE DMUMPS_759(
     &    COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root,
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        USE DMUMPS_LOAD
        IMPLICIT NONE
        INCLUDE 'dmumps_root.h'
        INCLUDE 'mumps_headers.h'
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
      INTEGER INODE, FPERE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS, LA
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER IWPOS, IWPOSCB
      INTEGER LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP, IFLAG, IERROR
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER NBFIN, SLAVEF
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC( N + KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER ND( KEEP(28) )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( max(1,KEEP(13)) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER ITYPE2
      INTEGER IHDR_REC
      PARAMETER (ITYPE2=2)
      INTEGER IOLDPS, NROW, LDA
      INTEGER NPIV, LCONT, NELIM, NASS, NCOL_TO_SEND,
     &        SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON
      INTEGER(8) :: SHIFT_VAL_SON
      INTEGER(8) MEM_GAIN 
        IF (KEEP(50).EQ.0) THEN
          IHDR_REC=6
        ELSE
          IHDR_REC=8
        ENDIF
        IOLDPS = PTRIST(STEP(INODE))
        IW(IOLDPS+XXS)=S_ALL
         IF (KEEP(214).EQ.1) THEN
          CALL DMUMPS_314( N, INODE,
     &    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     &     )
          IOLDPS = PTRIST(STEP(INODE))
          IF (KEEP(38).NE.FPERE) THEN
            IW(IOLDPS+XXS)=S_NOLCBNOCONTIG
            IF (KEEP(216).NE.3) THEN
             MEM_GAIN=int(IW( IOLDPS + 2 + KEEP(IXSZ) ),8)*
     &                int(IW( IOLDPS + 3 + KEEP(IXSZ) ),8)
             LRLUS = LRLUS+MEM_GAIN
             CALL DMUMPS_471(.FALSE.,.FALSE.,
     &              LA-LRLUS,0_8,-MEM_GAIN,KEEP,KEEP8,LRLU)
            ENDIF
          ENDIF
          IF (KEEP(216).EQ.2) THEN
           IF (FPERE.NE.KEEP(38)) THEN
           CALL DMUMPS_627(A,LA,PTRAST(STEP(INODE)),
     &         IW( IOLDPS + 2 + KEEP(IXSZ) ),
     &         IW( IOLDPS + KEEP(IXSZ) ),
     &         IW( IOLDPS + 3 + KEEP(IXSZ) )+
     &         IW( IOLDPS + KEEP(IXSZ) ), 0,
     &         IW( IOLDPS + XXS ), 0_8 )
           IW(IOLDPS+XXS)=S_NOLCBCONTIG
           IW(IOLDPS+XXS)=S_NOLCBCONTIG
           ENDIF
          ENDIF 
         ENDIF 
      IF ( KEEP(38).EQ.FPERE) THEN
       LCONT  = IW(IOLDPS+KEEP(IXSZ))
       NROW   = IW(IOLDPS+2+KEEP(IXSZ))
       NPIV   = IW(IOLDPS+3+KEEP(IXSZ))
       NASS   = IW(IOLDPS+4+KEEP(IXSZ))
       NELIM  = NASS-NPIV
       NCOL_TO_SEND =  LCONT-NELIM
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5+KEEP(IXSZ)) + KEEP(IXSZ)
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NASS
       SHIFT_VAL_SON      = int(NASS,8)
       LDA                = LCONT + NPIV
      IF (IW(IOLDPS+IHDR_REC+KEEP(IXSZ)).EQ.S_ROOTBAND_INIT) THEN
        IW(IOLDPS+IHDR_REC+KEEP(IXSZ)) = S_REC_CONTSTATIC
      ELSE
      ENDIF
       CALL DMUMPS_80( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    PTRIST, PTRAST, 
     &    root, NROW, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     &    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA, 
     &    ROOT_CONT_STATIC, MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, .FALSE., ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
       IF ( IFLAG < 0 ) GOTO 600
       IF (NELIM.EQ.0) THEN
         IF (KEEP(214).EQ.2) THEN
          CALL DMUMPS_314( N, INODE,  
     &    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     &    )
         ENDIF
         CALL DMUMPS_626( N, INODE, PTRIST, PTRAST, IW, LIW,
     &        A, LA, LRLU, LRLUS, IWPOSCB, IPTRLU, STEP,
     &        MYID, KEEP
     &         )
       ELSE
         IOLDPS = PTRIST(STEP(INODE))
         IF (IW(IOLDPS+IHDR_REC+KEEP(IXSZ)).EQ.S_ROOT2SON_CALLED) THEN
           CALL DMUMPS_626( N, INODE, PTRIST, PTRAST, IW, LIW,
     &        A, LA, LRLU, LRLUS, IWPOSCB, IPTRLU, STEP,
     &        MYID, KEEP
     &         )
         ELSE
          IW(IOLDPS+IHDR_REC+KEEP(IXSZ)) = S_ROOTBAND_INIT
          IF (KEEP(214).EQ.1.AND.KEEP(216).NE.3) THEN
           IW(IOLDPS+XXS)=S_NOLCBNOCONTIG38
           CALL DMUMPS_628( IW(IOLDPS),
     &                     LIW-IOLDPS+1,
     &                     MEM_GAIN, KEEP(IXSZ) )
           LRLUS = LRLUS + MEM_GAIN
              CALL DMUMPS_471(.FALSE.,.FALSE.,
     &                LA-LRLUS,0_8,-MEM_GAIN,KEEP,KEEP8,LRLU)
            IF (KEEP(216).EQ.2) THEN
              CALL DMUMPS_627(A,LA,PTRAST(STEP(INODE)),
     &         IW( IOLDPS + 2 + KEEP(IXSZ) ),
     &         IW( IOLDPS + KEEP(IXSZ) ),
     &         IW( IOLDPS + 3 + KEEP(IXSZ) )+
     &         IW( IOLDPS + KEEP(IXSZ) ),
     &         IW( IOLDPS + 4 + KEEP(IXSZ) ) -
     &         IW( IOLDPS + 3 + KEEP(IXSZ) ),
     &         IW( IOLDPS + XXS ),0_8)
              IW(IOLDPS+XXS)=S_NOLCBCONTIG38
            ENDIF
          ENDIF
         ENDIF 
       ENDIF 
      ENDIF 
 600  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_759
      SUBROUTINE DMUMPS_141( COMM_LOAD, ASS_IRECV, 
     &           N, INODE, FPERE, IW, LIW, A, LA,
     &           UU, NOFFW,
     &           NPVW,
     &             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &            DKEEP,PIVNUL_LIST,LPN_LIST )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, NOFFW, NPVW
      INTEGER(8) :: LA
      DOUBLE PRECISION A( LA )
      DOUBLE PRECISION UU, SEUIL
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER NBFIN, IFLAG, IERROR, LEAF, LPOOL, SLAVEF, 
     &        IWPOS, IWPOSCB, COMP 
      INTEGER NB_BLOC_FAC
      INTEGER ICNTL(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, TARGET :: IW( LIW )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL), ITLOC(N+KEEP(253))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND( KEEP(28) )
      INTEGER FRERE(KEEP(28)), FILS(N) 
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)),
     & PTLUST_S(KEEP(28)),
     &        
     & PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     &        PROCNODE_STEPS(KEEP(28)), STEP(N)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION  DBLARR(max(1,KEEP(13)))
      LOGICAL AVOID_DELAYED 
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER(8) :: POSELT
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEGKJI, NBOLKJ, NBTLKJ
      INTEGER NASS, IEND, IOLDPS, LDAFS,allocok, IBEG_BLOCK
      LOGICAL LASTBL 
      LOGICAL RESET_TO_ONE, TO_UPDATE
      INTEGER K109_ON_ENTRY
      INTEGER I,J,JJ,K,IDEB
      DOUBLE PRECISION UUTEMP
      INCLUDE 'mumps_headers.h'
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, NextPiv2beWritten,
     &        IDUMMY
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER PP_FIRST2SWAP_L, IFLAG_OOC
      INTEGER PP_LastPIVRPTRFilled 
      EXTERNAL DMUMPS_223, DMUMPS_235,
     &         DMUMPS_227, DMUMPS_294,
     &         DMUMPS_44
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      INTEGER PIVSIZ,IWPOSPIV
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0D0)
      INOPV = 0
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
        SEUIL_LOC=SEUIL
        UUTEMP=UU
      ENDIF
      RESET_TO_ONE = ((KEEP(110).GT.0).AND.(DKEEP(2).LE.0.0D0))
      IF (RESET_TO_ONE) THEN
        K109_ON_ENTRY = KEEP(109)
      ENDIF
      IBEG_BLOCK=1
      NB_BLOC_FAC = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST( STEP( INODE ))
      NFRONT = IW(IOLDPS+KEEP(IXSZ))
      NASS   = iabs(IW(IOLDPS+2+KEEP(IXSZ)))
      LDAFS  = NASS
      IF (NASS .GT. KEEP(3)) THEN
        NBOLKJ = min( KEEP(6), NASS )
      ELSE
        NBOLKJ = min( KEEP(5), NASS )
      ENDIF
      NBTLKJ = NBOLKJ
      IW(IOLDPS+3+KEEP(IXSZ)) = min0(NASS,NBTLKJ)
      IF (KEEP(201).EQ.1) THEN 
        IDUMMY    = -9876
        CALL MUMPS_729(LAFAC,IW(IOLDPS+XXR))
        LIWFAC    = IW(IOLDPS+XXI)
        TYPEFile  = TYPEF_L
        NextPiv2beWritten = 1 
        PP_FIRST2SWAP_L = NextPiv2beWritten 
        MonBloc%LastPanelWritten_L = 0 
        MonBloc%INODE    = INODE
        MonBloc%MASTER   = .TRUE.
        MonBloc%Typenode = 2
        MonBloc%NROW     = NASS
        MonBloc%NCOL     = NASS
        MonBloc%NFS      = NASS
        MonBloc%Last     = .FALSE.   
        MonBloc%LastPiv  = -66666    
        MonBloc%INDICES =>
     &  IW(IOLDPS+6+NFRONT+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ))
     &    :IOLDPS+5+2*NFRONT+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ)))
      ENDIF
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID,' : FACTO_NIV2 :failed to allocate ',NASS,
     & ' integers'
        IFLAG=-13
        IERROR=NASS
        GO TO 490
      END IF
 50   CONTINUE
      IBEGKJI = IBEG_BLOCK
      CALL DMUMPS_223(
     &                NFRONT,NASS,IBEGKJI, NASS, IPIV,
     &                N,INODE,IW,LIW,A,LA,NOFFW,INOPV,
     &                IFLAG,IOLDPS,POSELT,UU, SEUIL_LOC,
     &                KEEP,KEEP8,PIVSIZ,
     &           DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &           PP_FIRST2SWAP_L, MonBloc%LastPanelWritten_L,
     &           PP_LastPIVRPTRFilled)
      IF (IFLAG.LT.0) GOTO 490
      IF(KEEP(109).GT. 0) THEN
         IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
            IWPOSPIV = IOLDPS+IW(IOLDPS+1+KEEP(IXSZ))+6
     &              +IW(IOLDPS+5+KEEP(IXSZ))
            PIVNUL_LIST(KEEP(109)) = IW(IWPOSPIV+KEEP(IXSZ))
         ENDIF
      ENDIF
         IF(INOPV.EQ. 1 .AND. STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
      IF (INOPV.GE.1) THEN
          LASTBL = (INOPV.EQ.1)
          IEND = IW(IOLDPS+1+KEEP(IXSZ))
          CALL DMUMPS_294( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, LDAFS,
     &             IBEGKJI, IEND, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     &
     &             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
      ENDIF
      IF (INOPV.EQ.1) GO TO 500
      IF (INOPV.EQ.2) THEN
         CALL DMUMPS_235(IBEG_BLOCK,
     &            NASS,N,INODE,IW,LIW,A,LA,
     &            LDAFS, 
     &            IOLDPS,POSELT,NBOLKJ, NBTLKJ,KEEP(4),KEEP,KEEP8)
         GOTO 50
      ENDIF
      NPVW = NPVW + PIVSIZ
      IF (NASS.LE.1) THEN
        IFINB = -1
        IF (NASS == 1) A(POSELT)=ONE/A(POSELT)
      ELSE
         CALL DMUMPS_227(IBEG_BLOCK,
     &             NASS, N,INODE,IW,LIW,A,LA,
     &             LDAFS, IOLDPS,POSELT,IFINB,
     &             NBTLKJ,KEEP(4),PIVSIZ,KEEP(IXSZ))
         IF(PIVSIZ .EQ. 2) THEN
            IWPOSPIV = IOLDPS+KEEP(IXSZ)+IW(IOLDPS+1+KEEP(IXSZ))+6+
     &                 IW(IOLDPS+5+KEEP(IXSZ))
            IW(IWPOSPIV+NFRONT) = -IW(IWPOSPIV+NFRONT)
         ENDIF
      ENDIF
      IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+1+KEEP(IXSZ)) + PIVSIZ
       IF (IFINB.EQ.0) GOTO 50
       IF ((IFINB.EQ.1).OR.(IFINB.EQ.-1)) THEN
          LASTBL = (IFINB.EQ.-1) 
          IEND = IW(IOLDPS+1+KEEP(IXSZ))
          CALL DMUMPS_294(COMM_LOAD, ASS_IRECV, 
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, LDAFS, 
     &             IBEGKJI, IEND, IPIV, NASS, LASTBL,NB_BLOC_FAC,
     &
     &             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
       ENDIF
       IF (IFINB.EQ.(-1)) GOTO 500
       NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
       IF (KEEP(201).EQ.1) THEN
        IF (.NOT.RESET_TO_ONE.OR.K109_ON_ENTRY.EQ.KEEP(109)) THEN
         MonBloc%Last   = .FALSE.
         MonBloc%LastPiv= NPIV
         LAST_CALL=.FALSE.
         CALL DMUMPS_688(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
         IF (IFLAG_OOC .LT. 0 ) IFLAG = IFLAG_OOC
         IF (IFLAG .LT. 0 ) RETURN
        ENDIF
       ENDIF
      CALL DMUMPS_235(IBEG_BLOCK,
     &            NASS,N,INODE,IW,LIW,A,LA,
     &            LDAFS, 
     &            IOLDPS,POSELT,NBOLKJ,NBTLKJ,KEEP(4),KEEP,KEEP8)
      IF (KEEP(201).EQ.1) THEN
         IF (RESET_TO_ONE.AND.K109_ON_ENTRY.LT.KEEP(109)) THEN
          IDEB =  IOLDPS+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ))+6
          JJ= IDEB
          TO_UPDATE=.FALSE. 
          DO K = K109_ON_ENTRY+1, KEEP(109)
           I = PIVNUL_LIST(K)  
           DO J=JJ,JJ+NASS
            IF (IW(J).EQ.I) THEN
              TO_UPDATE=.TRUE. 
              EXIT
            ENDIF
           ENDDO
           IF (TO_UPDATE) THEN
            JJ= J  
            J =  J-IDEB+1
            A(POSELT+int(J-1,8)+int(LDAFS,8)*int(J-1,8))= ONE
            TO_UPDATE=.FALSE. 
           ELSE
            IF (ICNTL(1).GT.0) THEN
             write(ICNTL(1),*) ' Internal error related ', 
     &                 'to null pivot row detection'
            ENDIF
            EXIT
           ENDIF
          ENDDO
         ENDIF
         K109_ON_ENTRY = KEEP(109)
         MonBloc%Last   = .FALSE.
         MonBloc%LastPiv= NPIV
         LAST_CALL=.FALSE.
         CALL DMUMPS_688(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
         IF (IFLAG_OOC .LT. 0 ) IFLAG = IFLAG_OOC
         IF (IFLAG .LT. 0 ) RETURN
      ENDIF
      GO TO 50
 490  CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
 500  CONTINUE
      IF (RESET_TO_ONE.AND.K109_ON_ENTRY.LT.KEEP(109)) THEN
       IDEB =  IOLDPS+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ))+6
       JJ= IDEB
       TO_UPDATE=.FALSE. 
       DO K = K109_ON_ENTRY+1, KEEP(109)
        I = PIVNUL_LIST(K)  
        DO J=JJ,JJ+NASS
         IF (IW(J).EQ.I) THEN
           TO_UPDATE=.TRUE. 
           EXIT
         ENDIF
        ENDDO
        IF (TO_UPDATE) THEN
            JJ= J  
           J =  J-IDEB+1
           A(POSELT+int(J-1,8)+int(LDAFS,8)*int(J-1,8))= ONE
           TO_UPDATE=.FALSE. 
        ELSE
        IF (ICNTL(1).GT.0) THEN
         write(ICNTL(1),*) ' Internal error related ', 
     &                'to null pivot row detection'
        ENDIF
         EXIT
        ENDIF
       ENDDO
      ENDIF
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+KEEP(IXSZ))
          LAST_CALL = .TRUE.
          CALL DMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           NextPiv2beWritten, IDUMMY,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
         IF (IFLAG_OOC .LT. 0 ) IFLAG = IFLAG_OOC
         IF (IFLAG .LT. 0 ) RETURN
          CALL DMUMPS_644 (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE DMUMPS_141
      SUBROUTINE DMUMPS_223( NFRONT, NASS,
     &                   IBEGKJI, NASS2, TIPIV,
     &                   N, INODE, IW, LIW,
     &                   A, LA, NNEG, 
     &                   INOPV, IFLAG,
     &                   IOLDPS, POSELT, UU, 
     &                   SEUIL,KEEP,KEEP8,PIVSIZ,
     &                   DKEEP,PIVNUL_LIST,LPN_LIST,
     &                   PP_FIRST2SWAP_L, PP_LastPanelonDisk,
     &                   PP_LastPIVRPTRIndexFilled)
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,IFLAG,INOPV
      INTEGER NASS2, IBEGKJI, NNEG
      INTEGER TIPIV( NASS2 )
      INTEGER PIVSIZ,LPIV
      INTEGER(8) :: LA
      DOUBLE PRECISION A(LA) 
      DOUBLE PRECISION UU, UULOC, SEUIL
      DOUBLE PRECISION CSEUIL
      INTEGER IW(LIW) 
      INTEGER   IOLDPS
      INTEGER(8) :: POSELT
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk
      INTEGER PP_LastPIVRPTRIndexFilled
      include 'mpif.h'
      INTEGER(8) :: POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX
      DOUBLE PRECISION RMAX,AMAX,TMAX,TOL
      DOUBLE PRECISION MAXPIV
      DOUBLE PRECISION PIVOT,DETPIV
      PARAMETER(TOL = 1.0D-20)
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: APOSMAX
      INTEGER(8) :: APOS
      INTEGER(8) :: J1, J2, JJ, KK
      INTEGER    :: LDAFS
      INTEGER(8) :: LDAFS8
      DOUBLE PRECISION, PARAMETER :: RZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER :: RONE  = 1.0D0
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO = 0.0D0 )
      PARAMETER( ONE = 1.0D0 )
      DOUBLE PRECISION PIVNUL, VALTMP
      DOUBLE PRECISION FIXA
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,ILOC,K,J
      INTRINSIC max
      INTEGER I_PIVRPTR, I_PIVR, NBPANELS_L
      PIVNUL = DKEEP(1)
      FIXA   = DKEEP(2)
      CSEUIL = SEUIL
      LDAFS  = NASS
      LDAFS8 = int(LDAFS,8)
      IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
             CALL DMUMPS_667(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR, I_PIVR, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+KEEP(IXSZ))
     &              +KEEP(IXSZ),
     &       IW, LIW)
      ENDIF
        UULOC = UU
        PIVSIZ = 1
        NPIV    = IW(IOLDPS+1+KEEP(IXSZ))
        NPIVP1  = NPIV + 1
        ILOC = NPIVP1 - IBEGKJI + 1
        TIPIV( ILOC ) = ILOC
        NASSW   = iabs(IW(IOLDPS+3+KEEP(IXSZ)))
        APOSMAX = POSELT+LDAFS8*LDAFS8-1_8
        IF(INOPV .EQ. -1) THEN
           APOS = POSELT + LDAFS8*int(NPIVP1-1,8) + int(NPIV,8)
           POSPV1 = APOS
           IF(abs(A(APOS)).LT.SEUIL) THEN
              IF(dble(A(APOS)) .GE. RZERO) THEN
                 A(APOS) = CSEUIL
              ELSE
                 A(APOS) = -CSEUIL
              ENDIF
           ELSE IF (KEEP(258) .NE.0 ) THEN
             CALL DMUMPS_762( A(APOS), DKEEP(6), KEEP(259) )
           ENDIF
           IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
             CALL DMUMPS_680( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
        DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + LDAFS8*int(IPIV-1,8) + int(NPIV,8)
            POSPV1 = APOS + int(IPIV - NPIVP1,8)
            PIVOT = A(POSPV1)
            IF (UULOC.EQ.RZERO) THEN 
              IF (abs(A(APOS)).EQ.RZERO) GO TO 630
              IF (A(APOS).LT.RZERO) NNEG = NNEG+1
              IF (KEEP(258) .NE. 0) THEN
                CALL DMUMPS_762(A(APOS), DKEEP(6), KEEP(259))
              ENDIF
              GO TO 420
            ENDIF
            AMAX = RZERO
            JMAX = 0
            J1 = APOS
            J2 = POSPV1 - 1_8
            DO JJ=J1,J2
               IF(abs(A(JJ)) .GT. AMAX) THEN
                  AMAX = abs(A(JJ))
                  JMAX = IPIV - int(POSPV1-JJ)
               ENDIF
            ENDDO
            J1 = POSPV1 + LDAFS8
            DO J=1, NASSW - IPIV
               IF(abs(A(J1)) .GT. AMAX) THEN
                  AMAX = max(abs(A(J1)),AMAX)
                  JMAX = IPIV + J
               ENDIF
               J1 = J1 + LDAFS8
            ENDDO
            IF (KEEP(219).NE.0) THEN
             RMAX = dble(A(APOSMAX+int(IPIV,8)))
            ELSE
             RMAX = RZERO
            ENDIF
            DO J=1,NASS - NASSW
               RMAX = max(abs(A(J1)),RMAX)
               J1 = J1 + LDAFS8
            ENDDO
         IF (max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL) THEN
            KEEP(109) = KEEP(109)+1
            PIVNUL_LIST(KEEP(109)) = -1
            IF (dble(FIXA).GT.RZERO) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDAFS8
               DO J=1, NASSW - IPIV
                  A(J1) = ZERO
                  J1 = J1 + LDAFS8
               ENDDO
               DO J=1,NASS - NASSW
                  A(J1) = ZERO
                  J1 = J1 + LDAFS8
               ENDDO
                VALTMP = max(1.0D10*RMAX, sqrt(huge(RMAX))/1.0D8)
                A(POSPV1) = VALTMP
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
        IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
         IF (max(AMAX,RMAX,abs(PIVOT)).LE.TOL) THEN
            IF(SEUIL .GT. epsilon(SEUIL)) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = CSEUIL
               ELSE
                  A(POSPV1) = -CSEUIL
                  NNEG = NNEG+1
               ENDIF
               PIVOT = A(POSPV1)
               WRITE(*,*) 'WARNING matrix may be singular'
               KEEP(98) = KEEP(98)+1
               GO TO 415
            ENDIF
         ENDIF
        ENDIF
        IF (max(AMAX,abs(PIVOT)) .LE. TOL) GOTO 460
        IF (abs(PIVOT).GT.max(UULOC*max(RMAX,AMAX),SEUIL)) THEN
          IF (A(POSPV1).LT.RZERO) NNEG = NNEG+1
          IF (KEEP(258) .NE.0 ) THEN
            CALL DMUMPS_762(PIVOT, DKEEP(6), KEEP(259))
          ENDIF
          GO TO 415
        END IF
            IF (AMAX.LE.TOL) GO TO 460
            IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  IF(int(POSPV1-JJ) .NE. IPIV-JMAX) THEN
                     RMAX = max(RMAX,abs(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDAFS8
               DO J=1,NASS-IPIV
                  IF(IPIV+J .NE. JMAX) THEN
                     RMAX = max(abs(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDAFS8
               ENDDO
            ENDIF            
            APOSJ = POSELT + int(JMAX-1,8)*LDAFS8 + int(NPIV,8)
            POSPV2 = APOSJ + int(JMAX - NPIVP1,8)
            IF (IPIV.LT.JMAX) THEN
               OFFDAG = APOSJ + int(IPIV - NPIVP1,8)
            ELSE
               OFFDAG = APOS + int(JMAX - NPIVP1,8)
            END IF
            IF (KEEP(219).NE.0) THEN
             TMAX = max(SEUIL/UULOC,dble(A(APOSMAX+int(JMAX,8))))
            ELSE
             TMAX = SEUIL/UULOC
            ENDIF
            IF(JMAX .LT. IPIV) THEN
               JJ = POSPV2
               DO K = 1, NASS-JMAX
                  JJ = JJ+int(NASS,8) 
                  IF (JMAX+K.NE.IPIV) THEN
                     TMAX=max(TMAX,abs(A(JJ)))
                  ENDIF
               ENDDO
               DO KK =  APOSJ, POSPV2-1_8
                  TMAX = max(TMAX,abs(A(KK)))
               ENDDO
            ELSE
               JJ = POSPV2
               DO K = 1, NASS-JMAX
                  JJ = JJ+int(NASS,8) 
                  TMAX=max(TMAX,abs(A(JJ)))
               ENDDO
               DO KK =  APOSJ, POSPV2 - 1_8
                  IF (KK.NE.OFFDAG) THEN
                     TMAX = max(TMAX,abs(A(KK)))
                  ENDIF
               ENDDO
            ENDIF
            DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
            IF (SEUIL.GT.RZERO) THEN
               IF (sqrt(abs(DETPIV)) .LE. SEUIL ) GOTO 460
            ENDIF
            MAXPIV = max(abs(A(POSPV1)),abs(A(POSPV2)))
            IF (MAXPIV.EQ.RZERO) MAXPIV = RONE
            IF (abs(DETPIV)/MAXPIV.LE.TOL) GO TO 460
            IF ((abs(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     &           abs(DETPIV)) GO TO 460
            IF ((abs(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     &           abs(DETPIV)) GO TO 460
           IF (KEEP(258).NE.0) THEN
             CALL DMUMPS_762(DETPIV, DKEEP(6), KEEP(259))
           ENDIF
           PIVSIZ = 2
           KEEP(105) = KEEP(105)+1
           IF(DETPIV .LT. RZERO) THEN
             NNEG = NNEG+1
           ELSE IF(A(POSPV2) .LT. RZERO) THEN
             NNEG = NNEG+2
           ENDIF
 415       CONTINUE
           DO K=1,PIVSIZ
              IF (PIVSIZ .EQ. 2 ) THEN
                IF (K==1) THEN
                  LPIV = min(IPIV, JMAX)
                  TIPIV(ILOC) = -(LPIV - IBEGKJI + 1)
                ELSE
                  LPIV = max(IPIV, JMAX)
                  TIPIV(ILOC+1) = -(LPIV - IBEGKJI + 1)
                ENDIF
              ELSE
                LPIV = IPIV
                TIPIV(ILOC) = IPIV - IBEGKJI + 1
              ENDIF
              IF (LPIV.EQ.NPIVP1) THEN
                 GOTO 416
              ENDIF
              CALL DMUMPS_319( A, LA, IW, LIW,
     &             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     &             LDAFS, NFRONT, 2, KEEP(219), KEEP(50),
     &             KEEP(IXSZ))
 416          CONTINUE
              IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
                CALL DMUMPS_680( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, LPIV, PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
              NPIVP1 = NPIVP1+1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+LDAFS8*int(NPIV,8)+int(NPIV+1,8)) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
  420 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_223
      SUBROUTINE DMUMPS_235(
     &                 IBEG_BLOCK,
     &                 NASS, N, INODE,
     &                 IW, LIW, A, LA,
     &                 LDAFS,
     &                 IOLDPS, POSELT,
     &                 LKJIB_ORIG, LKJIB, LKJIT, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER NASS,N,LIW
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB_ORIG, LKJIB, INODE, KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER (8) :: POSELT
      INTEGER (8) :: LDAFS8
      INTEGER LDAFS, IBEG_BLOCK
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1
      INTEGER HF
      INTEGER(8) :: LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER J, Block
      INTEGER BLSIZE
      DOUBLE PRECISION ONE, ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INCLUDE 'mumps_headers.h'
      LDAFS8 = int(LDAFS,8)
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
      JROW2  = iabs(IW(IOLDPS+3+KEEP(IXSZ)))
      NPBEG  = IBEG_BLOCK
      HF     = 6 + IW(IOLDPS+5+KEEP(IXSZ)) + KEEP(IXSZ)
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB_ORIG + NONEL
          IW(IOLDPS+3+KEEP(IXSZ))= min0(NPIV+LKJIB,NASS)
          LKJIB       = min0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3+KEEP(IXSZ)) = NASS
        ENDIF
      ELSEIF (JROW2.LT.NASS) THEN
          IW(IOLDPS+3+KEEP(IXSZ)) = min0(JROW2+LKJIB,NASS)
      ENDIF
      IBEG_BLOCK = NPIV + 1
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = min( BLSIZE, NASS - IROW + 1 )
          LPOS = POSELT + int(IROW - 1,8) * LDAFS8 + int(NPBEG - 1,8)
          UPOS = POSELT + int(NPBEG - 1,8) * LDAFS8 + int(IROW - 1,8)
          APOS =  POSELT + int(IROW-1,8) * LDAFS8 + int(IROW - 1,8)
          DO J=1, Block
            CALL dgemv( 'T', LKJIW, Block - J + 1, ALPHA,
     &                  A( LPOS ), LDAFS, A( UPOS ), LDAFS,
     &                  ONE, A( APOS ), LDAFS )
            LPOS = LPOS + LDAFS8
            APOS = APOS + LDAFS8 + 1_8
            UPOS = UPOS + 1_8
          END DO
          LPOS = POSELT + int(IROW-1+Block,8) * LDAFS8
     &                  + int(NPBEG-1,8)
          UPOS = POSELT + int(NPBEG-1,8) * LDAFS8 + int(IROW-1,8)
          APOS = POSELT + int( IROW - 1 + Block,8 ) * LDAFS8
     &                  + int(IROW - 1,8)
          CALL dgemm( 'N','N', Block, NASS - IROW + 1 - Block, LKJIW,
     &                ALPHA, A( UPOS ), LDAFS,
     &                A( LPOS ), LDAFS, ONE, A( APOS ), LDAFS )
         END DO
        END IF
      END IF
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_235
      SUBROUTINE DMUMPS_227
     &     ( IBEG_BLOCK, NASS, N, INODE, IW, LIW,
     &     A, LA, LDAFS, 
     &     IOLDPS,POSELT,IFINB,LKJIB,LKJIT,PIVSIZ,
     &     XSIZE)
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER    :: LIW
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    VALPIV
      INTEGER IOLDPS, NCB1
      INTEGER LKJIT, IBEG_BLOCK
      INTEGER NPIV,JROW2
      INTEGER(8) :: APOS
      INTEGER(8) :: LPOS, LPOS1, LPOS2, K1POS 
      INTEGER(8) :: JJ, K1, K2
      INTEGER(8) :: POSPV1, POSPV2, OFFDAG, OFFDAG_OLD
      INTEGER(8) :: LDAFS8
      INTEGER NASS,N,INODE,IFINB,LKJIB,LDAFS,
     &        NPBEG
      INTEGER NEL2
      INTEGER XSIZE
      DOUBLE PRECISION ONE, ALPHA
      DOUBLE PRECISION ZERO
      INTEGER PIVSIZ,NPIV_NEW
      INTEGER(8) :: IBEG, IEND, IROW
      INTEGER    :: J2
      DOUBLE PRECISION SWOP,DETPIV,MULT1,MULT2
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      INCLUDE 'mumps_headers.h'
      LDAFS8 = int(LDAFS,8)
      NPIV   = IW(IOLDPS+1+XSIZE)
      NPIV_NEW = NPIV + PIVSIZ
      IFINB  = 0
      IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
         IW(IOLDPS+3+XSIZE) = min0(NASS,LKJIB)
      ENDIF
      JROW2 = IW(IOLDPS+3+XSIZE)
      NPBEG = IBEG_BLOCK
      NEL2   = JROW2 - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + int(NPIV,8)*(LDAFS8 + 1_8)
         VALPIV = ONE/A(APOS)
         A(APOS) = VALPIV
         LPOS   = APOS + LDAFS8
         CALL dcopy(NASS-NPIV_NEW, A(LPOS), LDAFS, A(APOS+1_8), 1)
         CALL DMUMPS_XSYR('U', NEL2, -VALPIV, A(LPOS), LDAFS,
     &        A(LPOS+1_8), LDAFS)
         CALL dscal(NASS-NPIV_NEW, VALPIV, A(LPOS), LDAFS)
         IF (NEL2.GT.0) THEN
            K1POS = LPOS + int(NEL2,8)*LDAFS8
            NCB1  = NASS - JROW2
            CALL dger(NEL2, NCB1 , ALPHA, A(APOS+1_8), 1, 
     &           A(K1POS), LDAFS, A(K1POS+1_8), LDAFS)
         ENDIF
      ELSE
         POSPV1 = POSELT + int(NPIV,8)*(LDAFS8 + 1_8)
         POSPV2 = POSPV1+LDAFS8+1_8
         OFFDAG_OLD = POSPV2 - 1_8
         OFFDAG = POSPV1+1_8
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
         A(POSPV2) = A(POSPV1)/DETPIV
         A(POSPV1) = SWOP/DETPIV
         A(OFFDAG) = -A(OFFDAG_OLD)/DETPIV
         A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDAFS8 - 1_8
         LPOS2   = LPOS1 + 1_8
         CALL dcopy(NASS-NPIV_NEW, A(LPOS1), LDAFS, A(POSPV1+2_8), 1)
         CALL dcopy(NASS-NPIV_NEW, A(LPOS2), LDAFS, A(POSPV2+1_8), 1)
         JJ = POSPV2 + int(NASS-1,8)
         IBEG = JJ + 2_8
         IEND = IBEG
         DO J2 = 1,NEL2
            K1 = JJ
            K2 = JJ+1_8
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2_8
            K2 = POSPV2+1_8
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1_8
               K2 = K2 + 1_8
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1_8) = -MULT2
            IBEG = IBEG + int(NASS,8) 
            IEND = IEND + int(NASS + 1,8)
            JJ = JJ+int(NASS,8)
         ENDDO
         IEND = IEND-1_8
         DO J2 = JROW2+1,NASS
            K1 = JJ
            K2 = JJ+1_8
            MULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1_8)*A(K2))
            MULT2 = - (A(POSPV1+1_8)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2_8
            K2 = POSPV2+1_8
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1_8
               K2 = K2 + 1_8
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1_8) = -MULT2
            IBEG = IBEG + int(NASS,8) 
            IEND = IEND + int(NASS,8) 
            JJ = JJ+int(NASS,8) 
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_227
      RECURSIVE SUBROUTINE DMUMPS_263(
     &   COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS,  
     &    PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &     )
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS, LA
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER NBPROCFILS(KEEP(28)), STEP(N), PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28))
      INTEGER ITLOC( N + KEEP(253)), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER ND( KEEP(28) ), FRERE_STEPS( KEEP(28) )
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      DOUBLE PRECISION  DBLARR( max(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      INTEGER INODE, IPOSK, JPOSK, NCOLU, NPIV, POSITION, IERR
      INTEGER(8) POSELT, POSBLOCFACTO
      INTEGER(8) LAELL
      INTEGER IOLDPS, LCONT1, NROW1, NCOL1, NPIV1
      INTEGER NSLAVES_TOT, HS, DEST, NSLAVES_FOLLOW
      INTEGER FPERE
      INTEGER(8) CPOS, LPOS
      LOGICAL DYNAMIC
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER allocok
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: UDYNAMIC
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      DYNAMIC = .FALSE.
      POSITION  = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPOSK, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, JPOSK, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF ( NPIV .LE. 0 ) THEN
      NPIV = - NPIV
        WRITE(*,*) MYID,':error, received negative NPIV in BLFAC'
        CALL MUMPS_ABORT()
      END IF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOLU, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LAELL = int(NPIV,8) * int(NCOLU,8)
      IF ( LRLU .LT. LAELL ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          CALL MUMPS_731(LAELL - LRLUS, IERROR)
          GOTO 700
        END IF
        CALL DMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ))
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     &       ,LRLU,LRLUS
             IFLAG = -9
             CALL MUMPS_731(LAELL - LRLU, IERROR)
             GOTO 700
        END IF
      END IF
      LRLU  = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL DMUMPS_471(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8, LAELL,KEEP,KEEP8,LRLU)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*NCOLU,
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) DYNAMIC = .TRUE.
      IF ( (PTRIST(STEP( INODE )).NE.0) .AND.
     &  (IPOSK + NPIV -1 .GT.
     &   IW(PTRIST(STEP(INODE))+3+KEEP(IXSZ))) )THEN
        DYNAMIC = .TRUE.
      ENDIF
      IF (DYNAMIC)  THEN
        ALLOCATE(UDYNAMIC(LAELL), stat=allocok)
        if (allocok .GT. 0) THEN
          write(*,*) MYID, ' : PB allocation U in blfac_slave '
     &     , LAELL
          IFLAG = -13 
          CALL MUMPS_731(LAELL,IERROR)
          GOTO 700
        endif
        UDYNAMIC(1_8:LAELL) = A(POSBLOCFACTO:POSBLOCFACTO+LAELL-1_8)
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL DMUMPS_471(.FALSE.,.FALSE.,
     &          LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLU)
      ENDIF
      DO WHILE ( PTRIST(STEP(INODE)) .EQ. 0 )
        MSGSOU = MUMPS_275( PROCNODE_STEPS(STEP(INODE)),
     &           SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_329( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MSGSOU, MAITRE_DESC_BANDE,
     &    STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
        IF ( IFLAG .LT. 0 ) GOTO 600
      ENDDO
      DO WHILE ( IPOSK + NPIV -1 .GT.
     &            IW( PTRIST(STEP( INODE )) + 3 +KEEP(IXSZ)) )
        MSGSOU = MUMPS_275( PROCNODE_STEPS(STEP(INODE)), SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_329( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MSGSOU, BLOC_FACTO_SYM,
     &    STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL DMUMPS_329( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &    STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
      IOLDPS  = PTRIST(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      LCONT1 = IW( IOLDPS + KEEP(IXSZ) )
      NROW1  = IW( IOLDPS + 2  + KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3  + KEEP(IXSZ))
      NSLAVES_TOT = IW( IOLDPS + 5  + KEEP(IXSZ))
      HS     = 6 + NSLAVES_TOT + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      CPOS = POSELT + int(JPOSK - 1,8)
      LPOS = POSELT + int(IPOSK - 1,8)
      IF ( NPIV .GT. 0 ) THEN
       IF (DYNAMIC) THEN
        CALL dgemm( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     &            UDYNAMIC(1), NPIV,
     &            A( LPOS ), NCOL1, ONE,
     &            A( CPOS ), NCOL1 )
       ELSE
        CALL dgemm( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     &            A( POSBLOCFACTO ), NPIV,
     &            A( LPOS ), NCOL1, ONE,
     &            A( CPOS ), NCOL1 )
       ENDIF
       FLOP1 = dble(NCOLU*NPIV)*dble(2*NROW1)
       FLOP1 = -FLOP1
       CALL DMUMPS_190(1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF
      IW( IOLDPS + 6 + KEEP(IXSZ) ) = IW( IOLDPS + 6+ KEEP(IXSZ) ) + 1
      IF (DYNAMIC) THEN
       DEALLOCATE(UDYNAMIC)
      ELSE
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL DMUMPS_471(.FALSE.,.FALSE.,
     &                      LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLU)
      ENDIF
      NSLAVES_FOLLOW = IW( IOLDPS + 5 +KEEP(IXSZ) ) - XTRA_SLAVES_SYM
      IF ( IW( IOLDPS + 6  +KEEP(IXSZ)) .eq. 0 .and.
     &     KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 0 )
     &     THEN
         DEST = MUMPS_275( PROCNODE_STEPS(STEP(INODE)), SLAVEF )
         CALL DMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     &                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_BLFAC_SLAVE.'
           IFLAG = -99
           GOTO 700
         END IF
      END IF
      IF (IW(PTRIST(STEP(INODE)) + 6+KEEP(IXSZ) ) .eq. 0) THEN
         CALL DMUMPS_759( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root,
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
       ENDIF 
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_263
      SUBROUTINE DMUMPS_38( NROW_SON, NCOL_SON, INDROW_SON, 
     &                     INDCOL_SON, NSUPCOL, VAL_SON, VAL_ROOT,
     &                     LOCAL_M, LOCAL_N,
     &                     RHS_ROOT, NLOC_ROOT, CBP )
      IMPLICIT NONE
      INTEGER NCOL_SON, NROW_SON, NSUPCOL
      INTEGER, intent(in) :: CBP   
      INTEGER INDROW_SON( NROW_SON ), INDCOL_SON( NCOL_SON )
      INTEGER LOCAL_M, LOCAL_N
      DOUBLE PRECISION VAL_SON( NCOL_SON, NROW_SON )
      DOUBLE PRECISION VAL_ROOT( LOCAL_M, LOCAL_N )
      INTEGER NLOC_ROOT
      DOUBLE PRECISION RHS_ROOT( LOCAL_M, NLOC_ROOT )
      INTEGER I, J
      IF (CBP .EQ. 0) THEN
        DO I = 1, NROW_SON
          DO J = 1, NCOL_SON-NSUPCOL
          VAL_ROOT( INDROW_SON( I ), INDCOL_SON( J ) ) =
     &    VAL_ROOT( INDROW_SON( I ), INDCOL_SON( J ) ) + VAL_SON(J,I)
          END DO
          DO J = NCOL_SON-NSUPCOL+1, NCOL_SON
            RHS_ROOT( INDROW_SON( I ), INDCOL_SON(J)) =
     &      RHS_ROOT(INDROW_SON(I),INDCOL_SON(J)) + VAL_SON(J,I)
          ENDDO
        END DO
      ELSE
        DO I=1, NROW_SON  
          DO J = 1, NCOL_SON  
           RHS_ROOT( INDROW_SON( I ), INDCOL_SON(J)) =
     &     RHS_ROOT(INDROW_SON(I),INDCOL_SON(J)) + VAL_SON(J,I)
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_38
      RECURSIVE SUBROUTINE DMUMPS_80
     &  ( COMM_LOAD, ASS_IRECV, N, ISON, IROOT,
     &    PTRI, PTRR,
     &    root,
     &    NBROW, NBCOL, SHIFT_LIST_ROW_SON,
     &    SHIFT_LIST_COL_SON,
     &    SHIFT_VAL_SON, LDA, TAG,
     &    MYID, COMM,
     &
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, INVERT, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE DMUMPS_OOC        
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, ISON, IROOT, TAG
      INTEGER PTRI( KEEP(28) )
      INTEGER(8) :: PTRR( KEEP(28) )
      INTEGER NBROW, NBCOL, LDA
      INTEGER(8) :: SHIFT_VAL_SON
      INTEGER SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON
      INTEGER MYID, COMM
      LOGICAL INVERT
      INCLUDE 'mpif.h'
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS, LA
      INTEGER IWPOS, IWPOSCB
      INTEGER LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) ), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28)), NSTK( N )
      INTEGER COMP, IFLAG, IERROR
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER NBFIN, SLAVEF
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER ITLOC( N + KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( max(1,KEEP(13)) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTRROW,  PTRCOL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NSUPROW, NSUPCOL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ROW_INDEX_LIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: COL_INDEX_LIST
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER I, POS_IN_ROOT, IROW, JCOL, IGLOB, JGLOB
      INTEGER PDEST, IERR
      INTEGER LOCAL_M, LOCAL_N
      INTEGER(8) :: POSROOT
      INTEGER NSUBSET_ROW, NSUBSET_COL
      INTEGER NRLOCAL, NCLOCAL
      LOGICAL SET_IRECV, BLOCKING, MESSAGE_RECEIVED
      INTEGER NBROWS_ALREADY_SENT
      INTEGER SIZE_MSG
      INTEGER LP
      INCLUDE 'mumps_headers.h'
      LOGICAL SKIPLAST_RHS_ROWS, BCP_SYM_NONEMPTY
      INTEGER BBPCBP
      BBPCBP  = 0   
      LP = ICNTL(1)
      IF ( ICNTL(4) .LE. 0 ) LP = -1
      ALLOCATE(PTRROW(root%NPROW + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPROW + 1
      endif
      ALLOCATE(PTRCOL(root%NPCOL + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPCOL + 1
      endif
      ALLOCATE(NSUPROW(root%NPROW + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPROW + 1
      endif
      ALLOCATE(NSUPCOL(root%NPCOL + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPCOL + 1
      endif
      IF (IFLAG.LT.0) THEN
         IF (LP > 0) write(6,*) MYID, ' : MEMORY ALLOCATION ',
     &     'FAILURE in DMUMPS_80'
         CALL DMUMPS_44( MYID, SLAVEF, COMM )
         RETURN
      ENDIF
      SKIPLAST_RHS_ROWS = ((KEEP(253).GT.0).AND.(KEEP(50).EQ.0))
      BCP_SYM_NONEMPTY = .FALSE.
      PTRROW = 0
      PTRCOL = 0
      NSUPROW = 0
      NSUPCOL = 0
      DO I = 1, NBROW                   
        IGLOB  =  IW( PTRI(STEP(ISON)) +
     &                          SHIFT_LIST_ROW_SON + I - 1 )
        IF (SKIPLAST_RHS_ROWS.AND.(IGLOB.GT.N)) CYCLE
        IF ( .NOT. INVERT ) THEN 
          IF (IGLOB.GT.N) THEN
            BCP_SYM_NONEMPTY = .TRUE.
            POS_IN_ROOT = IGLOB - N
            JCOL =  mod((POS_IN_ROOT-1)/root%NBLOCK,root%NPCOL) 
            NSUPCOL(JCOL+1) = NSUPCOL(JCOL+1) + 1
            PTRCOL( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
          ELSE  
            POS_IN_ROOT = root%RG2L_ROW( IGLOB ) 
            IROW  = mod((POS_IN_ROOT-1)/root%MBLOCK,root%NPROW)
            PTRROW ( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
          ENDIF
        ELSE        
          IF (IGLOB .GT. N) THEN 
            POS_IN_ROOT = IGLOB - N
          ELSE  
            POS_IN_ROOT = root%RG2L_COL( IGLOB ) 
          ENDIF
          JCOL =  mod( ( POS_IN_ROOT - 1 ) / root%NBLOCK, root%NPCOL ) 
          IF (IGLOB.GT.N)  
     &               NSUPCOL(JCOL+1) = NSUPCOL(JCOL+1) + 1
          PTRCOL( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
        END IF
      END DO
      IF ((KEEP(50).NE.0).AND.(.NOT.INVERT).AND.BCP_SYM_NONEMPTY) 
     &             BBPCBP = 1
      DO I = 1, NBCOL                   
        JGLOB   =  IW( PTRI(STEP(ISON)) +
     &                SHIFT_LIST_COL_SON + I - 1 ) 
        IF ((KEEP(50).GT.0) .AND. (JGLOB.GT.N)) CYCLE  
        IF ( .NOT. INVERT ) THEN
          IF (KEEP(50).EQ.0) THEN
            IF (JGLOB.LE.N) THEN
              POS_IN_ROOT = root%RG2L_COL(JGLOB)  
            ELSE
              POS_IN_ROOT = JGLOB-N
            ENDIF
            JCOL =  mod((POS_IN_ROOT-1) / root%NBLOCK, root%NPCOL ) 
            IF (JGLOB.GT.N) THEN
             NSUPCOL(JCOL+1) = NSUPCOL(JCOL+1)  + 1  
            ENDIF
            PTRCOL ( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
          ELSE 
            POS_IN_ROOT = root%RG2L_COL(JGLOB) 
            JCOL =  mod((POS_IN_ROOT-1) / root%NBLOCK, root%NPCOL )
            PTRCOL ( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
            IF (BCP_SYM_NONEMPTY) THEN
             POS_IN_ROOT = root%RG2L_ROW(JGLOB) 
             IROW  = mod((POS_IN_ROOT-1)/root%MBLOCK,root%NPROW)
             NSUPROW(IROW+1) = NSUPROW(IROW+1)+1
             PTRROW( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
            ENDIF
          ENDIF
        ELSE  
          IF (JGLOB.LE.N) THEN
           POS_IN_ROOT = root%RG2L_ROW( JGLOB ) 
          ELSE
           POS_IN_ROOT = JGLOB-N
          ENDIF
          IROW        = mod( ( POS_IN_ROOT - 1 ) /
     &                  root%MBLOCK, root%NPROW )
          PTRROW ( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
        END IF
      END DO
      PTRROW( 1 ) = 1
      DO IROW = 2, root%NPROW + 1
        PTRROW( IROW ) = PTRROW( IROW ) + PTRROW( IROW - 1 )
      END DO
      PTRCOL( 1 ) = 1
      DO JCOL = 2, root%NPCOL + 1
        PTRCOL( JCOL ) = PTRCOL( JCOL ) + PTRCOL( JCOL - 1 )
      END DO
      ALLOCATE(ROW_INDEX_LIST(PTRROW(root%NPROW+1)-1+1),
     &         stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = PTRROW(root%NPROW+1)-1+1
      endif
      ALLOCATE(COL_INDEX_LIST(PTRCOL(root%NPCOL+1)-1+1),
     &         stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = PTRCOL(root%NPCOL+1)-1+1
      endif
      DO I = 1, NBROW
        IGLOB  =  IW( PTRI(STEP(ISON)) +
     &                          SHIFT_LIST_ROW_SON + I - 1 )
        IF (SKIPLAST_RHS_ROWS.AND.(IGLOB.GT.N)) CYCLE
        IF ( .NOT. INVERT ) THEN
          IF (IGLOB.GT.N) CYCLE   
          POS_IN_ROOT = root%RG2L_ROW( IGLOB )
          IROW        = mod( ( POS_IN_ROOT - 1 ) / root%MBLOCK,
     &                       root%NPROW )
          ROW_INDEX_LIST( PTRROW( IROW + 1 ) ) = I 
          PTRROW ( IROW + 1 ) = PTRROW( IROW + 1 ) + 1
        ELSE
          IF (IGLOB.LE.N) THEN
           POS_IN_ROOT = root%RG2L_COL( IGLOB )
          ELSE
           POS_IN_ROOT = IGLOB - N  
          ENDIF
          JCOL        = mod( ( POS_IN_ROOT - 1 ) / root%NBLOCK,
     &                       root%NPCOL )
          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I 
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1 
        END IF
      END DO
      DO I = 1, NBCOL 
        JGLOB =  IW( PTRI(STEP(ISON))+SHIFT_LIST_COL_SON+I - 1 ) 
        IF ((KEEP(50).GT.0) .AND. (JGLOB.GT.N)) CYCLE  
        IF ( .NOT. INVERT ) THEN
          IF ( JGLOB.LE.N ) THEN
           POS_IN_ROOT = root%RG2L_COL( JGLOB )
          ELSE
           POS_IN_ROOT = JGLOB - N
          ENDIF
          JCOL        = mod( ( POS_IN_ROOT - 1 ) /
     &               root%NBLOCK, root%NPCOL )
          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I 
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1
        ELSE
          IF ( JGLOB.LE.N ) THEN
           POS_IN_ROOT = root%RG2L_ROW( JGLOB )
          ELSE
           POS_IN_ROOT = JGLOB - N
          ENDIF
          IROW        = mod( ( POS_IN_ROOT - 1 ) /
     &                root%MBLOCK, root%NPROW )
          ROW_INDEX_LIST( PTRROW( IROW + 1 ) ) = I    
          PTRROW( IROW + 1 ) = PTRROW( IROW + 1 ) + 1 
        END IF
      END DO
      IF (BCP_SYM_NONEMPTY) THEN
        DO I = 1, NBROW
          IGLOB  =  IW( PTRI(STEP(ISON)) +
     &                         SHIFT_LIST_ROW_SON + I - 1 )
          IF (IGLOB.LE.N) CYCLE  
          POS_IN_ROOT = IGLOB - N
          JCOL =  mod((POS_IN_ROOT-1)/root%NBLOCK,root%NPCOL)
          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I 
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1
        ENDDO
        DO I=1, NBCOL
         JGLOB =  IW( PTRI(STEP(ISON))+SHIFT_LIST_COL_SON+I - 1 ) 
         IF (JGLOB.GT.N) THEN
           EXIT
         ELSE
           POS_IN_ROOT = root%RG2L_ROW(JGLOB) 
         ENDIF
         IROW  = mod((POS_IN_ROOT-1)/root%MBLOCK,root%NPROW)
         ROW_INDEX_LIST( PTRROW( IROW + 1 ) ) = I    
         PTRROW( IROW + 1 ) = PTRROW( IROW + 1 ) + 1 
        ENDDO
      ENDIF
      DO IROW = root%NPROW, 2, -1
        PTRROW( IROW ) = PTRROW( IROW - 1 )
      END DO
      PTRROW( 1 ) = 1
      DO JCOL = root%NPCOL, 2, -1
        PTRCOL( JCOL ) = PTRCOL( JCOL - 1 )
      END DO
      PTRCOL( 1 ) = 1
      JCOL  = root%MYCOL
      IROW  = root%MYROW
      IF ( root%yes ) THEN
        if (IROW .ne. root%MYROW .or. JCOL.ne.root%MYCOL) then
        write(*,*) ' error in grid position buildandsendcbroot'
        CALL MUMPS_ABORT()
        end if
        IF ( PTRIST(STEP(IROOT)).EQ.0.AND.
     &       PTLUST_S(STEP(IROOT)).EQ.0) THEN
           NBPROCFILS( STEP(IROOT) ) = -1
           CALL DMUMPS_284(root, IROOT, N, IW, LIW,
     &                     A, LA,
     &                     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &                     LRLU, IPTRLU,
     &                     IWPOS, IWPOSCB, PTRIST, PTRAST,
     &                     STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &                     COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
           IF (IFLAG.LT.0) THEN
                CALL DMUMPS_44( MYID, SLAVEF, COMM )
                RETURN
           ENDIF
        ELSE
           NBPROCFILS(STEP(IROOT)) = NBPROCFILS(STEP(IROOT)) - 1
           IF ( NBPROCFILS( STEP(IROOT) ) .eq. 0 ) THEN
              IF (KEEP(201).EQ.1) THEN
                 CALL DMUMPS_681(IERR)
              ELSE IF (KEEP(201).EQ.2) THEN
                 CALL DMUMPS_580(IERR)              
              ENDIF
              CALL DMUMPS_507(N, IPOOL, LPOOL, PROCNODE_STEPS,
     &        SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     &        STEP, IROOT+N )
              IF (KEEP(47) .GE. 3) THEN
                 CALL DMUMPS_500(
     &                IPOOL, LPOOL, 
     &                PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &                MYID, STEP, N, ND, FILS )
              ENDIF
          END IF
        END IF
       IF (KEEP(60) .NE. 0 ) THEN
         LOCAL_M = root%SCHUR_LLD
         LOCAL_N = root%SCHUR_NLOC
            NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
            NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
            CALL DMUMPS_285( N,
     &        root%SCHUR_POINTER(1),
     &        LOCAL_M, LOCAL_N,
     &        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     &        NBCOL, NBROW,
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     &        LDA, A( PTRR( STEP(ISON)) + SHIFT_VAL_SON ),
     &        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     &        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     &        NRLOCAL,
     &        NCLOCAL,
     &        NSUPROW(IROW+1), NSUPCOL(JCOL+1),
     &        root%RG2L_ROW(1), root%RG2L_COL(1), INVERT,
     &        KEEP,
     &        root%RHS_ROOT(1,1), root%RHS_NLOC )
       ELSE
        IF ( PTRIST(STEP( IROOT )) .GE. 0 ) THEN
          IF ( PTRIST(STEP( IROOT )) .EQ. 0 ) THEN
            LOCAL_N = IW( PTLUST_S(STEP(IROOT)) + 1 + KEEP(IXSZ))
            LOCAL_M = IW( PTLUST_S(STEP(IROOT)) + 2 + KEEP(IXSZ))
            POSROOT = PTRFAC(IW( PTLUST_S(STEP(IROOT)) +4+KEEP(IXSZ) ))
          ELSE
            LOCAL_N = - IW( PTRIST(STEP(IROOT)) +KEEP(IXSZ))
            LOCAL_M = IW( PTRIST(STEP(IROOT)) + 1 +KEEP(IXSZ))
            POSROOT = PAMASTER(STEP( IROOT ))
          ENDIF
          NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
          NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
          CALL DMUMPS_285( N, A( POSROOT ),
     &        LOCAL_M, LOCAL_N,
     &        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     &        NBCOL, NBROW,
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     &        LDA, A( PTRR( STEP(ISON)) + SHIFT_VAL_SON ),
     &        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     &        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     &        NRLOCAL,
     &        NCLOCAL,
     &        NSUPROW(IROW+1), NSUPCOL(JCOL+1),
     &        root%RG2L_ROW(1), root%RG2L_COL(1), INVERT,
     &        KEEP,
     &        root%RHS_ROOT(1,1), root%RHS_NLOC )
        END IF
       ENDIF
      END IF
      DO IROW = 0, root%NPROW - 1
        DO JCOL = 0, root%NPCOL - 1
          PDEST = IROW * root%NPCOL + JCOL
          IF ( (root%MYROW.eq.IROW.and.root%MYCOL.eq.JCOL) .and.
     &      MYID.ne.PDEST) THEN
            write(*,*) 'error: myrow,mycol=',root%MYROW,root%MYCOL
            write(*,*) ' MYID,PDEST=',MYID,PDEST
            CALL MUMPS_ABORT()
          END IF
          IF ( root%MYROW .NE. IROW .OR. root%MYCOL .NE. JCOL) THEN
            NBROWS_ALREADY_SENT = 0
            IERR = -1
            DO WHILE ( IERR .EQ. -1 )
              NSUBSET_ROW = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
              NSUBSET_COL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
              IF ( LRLU .LT. int(NSUBSET_ROW,8) * int(NSUBSET_COL,8)
     &        .AND. LRLUS .GT. int(NSUBSET_ROW,8) * int(NSUBSET_COL,8) )
     &        THEN
                CALL DMUMPS_94(N, KEEP(28),
     &          IW, LIW, A, LA,
     &          LRLU, IPTRLU,
     &          IWPOS, IWPOSCB, PTRIST, PTRAST,
     &          STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &          KEEP(IXSZ))
                COMP = COMP + 1
                IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) MYID,': Error in b&scbroot: pb compress'
                  WRITE(*,*) MYID,': LRLU, LRLUS=',LRLU,LRLUS
                  CALL MUMPS_ABORT()
                END IF
              END IF
              CALL DMUMPS_648( N, ISON,
     &        NBCOL, NBROW,
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     &        LDA, A( PTRR(STEP(ISON)) + SHIFT_VAL_SON ),
     &        TAG,
     &        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     &        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     &        NSUBSET_ROW, NSUBSET_COL,
     &        NSUPROW(IROW+1), NSUPCOL(JCOL+1),
     &        root%NPROW, root%NPCOL, root%MBLOCK,
     &        root%RG2L_ROW, root%RG2L_COL,
     &        root%NBLOCK, PDEST,
     &        COMM, IERR, A( POSFAC ), LRLU, INVERT,
     &        SIZE_MSG, NBROWS_ALREADY_SENT, KEEP, BBPCBP )
              IF ( IERR .EQ. -1 ) THEN
                  BLOCKING  = .FALSE.
                  SET_IRECV = .TRUE.
                  MESSAGE_RECEIVED = .FALSE.
                  CALL DMUMPS_329( COMM_LOAD, ASS_IRECV, 
     &            BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &            MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &            STATUS, BUFR, LBUFR,
     &            LBUFR_BYTES, PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB,
     &            IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA,
     &            PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &            PIMASTER, PAMASTER, NSTK,
     &            COMP, IFLAG, IERROR, COMM, NBPROCFILS, IPOOL, LPOOL,
     &            LEAF, NBFIN, MYID, SLAVEF, root,
     &            OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS,
     &            PTRARW, PTRAIW, INTARR, DBLARR, ICNTL, KEEP,KEEP8,
     &            ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     &            ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.  )
                  IF ( IFLAG .LT. 0 ) GOTO 500
              END IF
            END DO
            IF ( IERR == -2 ) THEN
              IFLAG  = -17
              IERROR = SIZE_MSG
              IF (LP > 0) WRITE(LP, *) "FAILURE, SEND BUFFER TOO
     & SMALL DURING DMUMPS_80"
              CALL DMUMPS_44( MYID, SLAVEF, COMM )
              GOTO 500
            ENDIF
            IF ( IERR == -3 ) THEN
              IF (LP > 0) WRITE(LP, *) "FAILURE, RECV BUFFER TOO
     & SMALL DURING DMUMPS_80"
              IFLAG  = -20
              IERROR = SIZE_MSG
              CALL DMUMPS_44( MYID, SLAVEF, COMM )
              GOTO 500
            ENDIF
          END IF
        END DO
      END DO
 500  CONTINUE
      DEALLOCATE(PTRROW)
      DEALLOCATE(PTRCOL)
      DEALLOCATE(ROW_INDEX_LIST)
      DEALLOCATE(COL_INDEX_LIST)
      RETURN
      END SUBROUTINE DMUMPS_80
      SUBROUTINE DMUMPS_285( N, VAL_ROOT,
     &   LOCAL_M, LOCAL_N,
     &   NPCOL, NPROW, MBLOCK, NBLOCK, NBCOL_SON, NBROW_SON, INDCOL_SON,
     &   INDROW_SON, LD_SON, VAL_SON, SUBSET_ROW, SUBSET_COL,
     &   NSUBSET_ROW, NSUBSET_COL, NSUPROW, NSUPCOL,
     &   RG2L_ROW, RG2L_COL, INVERT,
     &   KEEP, RHS_ROOT, NLOC  )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER N, LOCAL_M, LOCAL_N
      DOUBLE PRECISION VAL_ROOT( LOCAL_M, LOCAL_N )
      INTEGER NPCOL, NPROW, MBLOCK, NBLOCK
      INTEGER NBCOL_SON, NBROW_SON
      INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
      INTEGER LD_SON 
      INTEGER NSUPROW, NSUPCOL
      DOUBLE PRECISION VAL_SON( LD_SON, NBROW_SON )
      INTEGER KEEP(500)
      INTEGER NSUBSET_ROW, NSUBSET_COL
      INTEGER SUBSET_ROW( NSUBSET_ROW ), SUBSET_COL( NSUBSET_COL )
      INTEGER RG2L_ROW( N ), RG2L_COL( N )
      LOGICAL INVERT
      INTEGER NLOC
      DOUBLE PRECISION RHS_ROOT( LOCAL_M, NLOC)
      INTEGER ISUB, JSUB, I, J, IPOS_ROOT, JPOS_ROOT
      INTEGER ILOC_ROOT, JLOC_ROOT, IGLOB, JGLOB
      IF (KEEP(50).EQ.0) THEN
        DO ISUB = 1, NSUBSET_ROW
          I         = SUBSET_ROW( ISUB )
          IGLOB     = INDROW_SON( I )
          IPOS_ROOT = RG2L_ROW( IGLOB )   
          ILOC_ROOT = MBLOCK
     &            * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &            + mod( IPOS_ROOT - 1, MBLOCK ) + 1
          DO JSUB = 1, NSUBSET_COL-NSUPCOL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDCOL_SON( J )
            JPOS_ROOT = RG2L_COL( JGLOB )
            JLOC_ROOT = NBLOCK
     &              * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &              + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &           VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
          END DO
          DO JSUB = NSUBSET_COL-NSUPCOL+1, NSUBSET_COL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDCOL_SON( J )
             JPOS_ROOT = JGLOB - N  
             JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              RHS_ROOT(ILOC_ROOT, JLOC_ROOT) =  
     &            RHS_ROOT(ILOC_ROOT, JLOC_ROOT) + VAL_SON( J, I )
          ENDDO
        END DO
      ELSE
        IF ( .NOT. INVERT ) THEN
          DO ISUB = 1, NSUBSET_ROW - NSUPROW 
            I         = SUBSET_ROW( ISUB )
            IGLOB     = INDROW_SON( I )
            IPOS_ROOT = RG2L_ROW( IGLOB )
            ILOC_ROOT = MBLOCK
     &            * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &            + mod( IPOS_ROOT - 1, MBLOCK ) + 1
            DO JSUB = 1, NSUBSET_COL -NSUPCOL
              J         = SUBSET_COL( JSUB )
              JGLOB     = INDCOL_SON( J )
              JPOS_ROOT = RG2L_COL( JGLOB )
              JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
            END DO
          END DO
          DO JSUB = NSUBSET_COL -NSUPCOL+1, NSUBSET_COL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDROW_SON( J )  
            JPOS_ROOT = JGLOB - N  
            JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            DO ISUB = NSUBSET_ROW - NSUPROW +1, NSUBSET_ROW
              I         = SUBSET_ROW( ISUB )
              IGLOB     = INDCOL_SON( I )  
              IPOS_ROOT = RG2L_ROW(IGLOB)
              ILOC_ROOT = MBLOCK
     &            * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &            + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              RHS_ROOT(ILOC_ROOT, JLOC_ROOT) =  
     &            RHS_ROOT(ILOC_ROOT, JLOC_ROOT) + VAL_SON( I, J )
            END DO
          END DO
        ELSE
          DO ISUB = 1, NSUBSET_COL-NSUPCOL 
            I         = SUBSET_COL( ISUB )
            IGLOB     = INDROW_SON( I )
            JPOS_ROOT = RG2L_COL( IGLOB )
            JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            DO JSUB = 1, NSUBSET_ROW
              J         = SUBSET_ROW( JSUB )
              JGLOB     = INDCOL_SON( J )
              IPOS_ROOT = RG2L_ROW( JGLOB )  
              ILOC_ROOT = MBLOCK
     &                * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
            END DO
           ENDDO
           DO ISUB = NSUBSET_COL-NSUPCOL+1, NSUBSET_COL
            I         = SUBSET_COL( ISUB )
            IGLOB     = INDROW_SON( I )
            JPOS_ROOT = IGLOB - N 
            JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            DO JSUB = 1, NSUBSET_ROW
              J         = SUBSET_ROW( JSUB )
              JGLOB     = INDCOL_SON( J )
              IPOS_ROOT = RG2L_ROW( JGLOB )  
              ILOC_ROOT = MBLOCK
     &                * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              RHS_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &            RHS_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
            END DO
           ENDDO
        END IF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_285
      SUBROUTINE DMUMPS_164
     &( MYID, NPROCS, N, root, COMM_ROOT, IROOT, FILS,
     &  K50, K46, K51
     &     , K60, IDNPROW, IDNPCOL, IDMBLOCK, IDNBLOCK
     & )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER MYID, MYID_ROOT
      TYPE (DMUMPS_ROOT_STRUC)::root
      INTEGER COMM_ROOT
      INTEGER N, IROOT, NPROCS, K50, K46, K51
      INTEGER FILS( N )
      INTEGER K60, IDNPROW, IDNPCOL, IDMBLOCK, IDNBLOCK
      INTEGER INODE, NPROWtemp, NPCOLtemp
      LOGICAL SLAVE
      root%ROOT_SIZE     = 0
      root%TOT_ROOT_SIZE = 0
      SLAVE = ( MYID .ne. 0 .or.
     &        ( MYID .eq. 0 .and. K46 .eq. 1 ) )
      INODE = IROOT
      DO WHILE ( INODE .GT. 0 )
        INODE = FILS( INODE )
        root%ROOT_SIZE = root%ROOT_SIZE + 1
      END DO
      IF ( ( K60 .NE. 2 .AND. K60 .NE. 3 ) .OR.
     &       IDNPROW .LE. 0 .OR. IDNPCOL .LE. 0
     &      .OR. IDMBLOCK .LE.0 .OR. IDNBLOCK.LE.0
     &      .OR. IDNPROW * IDNPCOL .GT. NPROCS ) THEN
        root%MBLOCK = K51
        root%NBLOCK = K51
        CALL DMUMPS_99( NPROCS, root%NPROW, root%NPCOL,
     &                         root%ROOT_SIZE, K50 )
        IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
          IDNPROW = root%NPROW
          IDNPCOL = root%NPCOL
          IDMBLOCK = root%MBLOCK
          IDNBLOCK = root%NBLOCK
        ENDIF
      ELSE IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
        root%NPROW = IDNPROW
        root%NPCOL = IDNPCOL
        root%MBLOCK = IDMBLOCK
        root%NBLOCK = IDNBLOCK
      ENDIF
      IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
        IF (SLAVE) THEN
          root%LPIV = 0
          IF (K46.EQ.0) THEN
            MYID_ROOT=MYID-1
          ELSE
            MYID_ROOT=MYID
          ENDIF
          IF (MYID_ROOT < root%NPROW*root%NPCOL) THEN
            root%MYROW = MYID_ROOT / root%NPCOL
            root%MYCOL = mod(MYID_ROOT, root%NPCOL)
            root%yes  = .true.
          ELSE
            root%MYROW = -1
            root%MYCOL = -1
            root%yes  = .FALSE.
          ENDIF
        ELSE
          root%yes  = .FALSE.
        ENDIF
      ELSE IF ( SLAVE ) THEN
        IF ( root%gridinit_done) THEN
           CALL blacs_gridexit( root%CNTXT_BLACS )
           root%gridinit_done = .FALSE.
        END IF
        root%CNTXT_BLACS = COMM_ROOT
        CALL blacs_gridinit( root%CNTXT_BLACS, 'R',
     &                       root%NPROW, root%NPCOL )
        root%gridinit_done = .TRUE.
        CALL blacs_gridinfo( root%CNTXT_BLACS,
     &                       NPROWtemp, NPCOLtemp,
     &                       root%MYROW, root%MYCOL )
        IF ( root%MYROW .NE. -1 ) THEN
          root%yes = .true.
        ELSE
          root%yes = .false.
        END IF
        root%LPIV = 0
      ELSE
        root%yes = .FALSE.
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_164
      SUBROUTINE DMUMPS_165( N, root, FILS, IROOT,
     &                                 KEEP, INFO )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE ( DMUMPS_ROOT_STRUC ):: root
      INTEGER N, IROOT, INFO(40), KEEP(500)
      INTEGER FILS( N )
      INTEGER INODE, I, allocok
      IF ( associated( root%RG2L_ROW ) ) DEALLOCATE( root%RG2L_ROW )
      IF ( associated( root%RG2L_COL ) ) DEALLOCATE( root%RG2L_COL )
      ALLOCATE( root%RG2L_ROW( N ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        RETURN
      ENDIF
      ALLOCATE( root%RG2L_COL( N ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        RETURN
      ENDIF
      INODE = IROOT
      I = 1
      DO WHILE ( INODE .GT. 0 )
        root%RG2L_ROW( INODE ) = I
        root%RG2L_COL( INODE ) = I
        I = I + 1
        INODE = FILS( INODE )
      END DO
      RETURN
      END SUBROUTINE DMUMPS_165
      SUBROUTINE DMUMPS_99( NPROCS, NPROW, NPCOL, SIZE, K50 )
      IMPLICIT NONE
      INTEGER NPROCS, NPROW, NPCOL, SIZE, K50
      INTEGER NPROWtemp, NPCOLtemp, NPROCSused, FLATNESS
      LOGICAL KEEPIT
      IF ( K50 .EQ. 1 ) THEN
        FLATNESS = 2
      ELSE
        FLATNESS = 3
      ENDIF
      NPROW  = int(sqrt(dble(NPROCS)))
      NPROWtemp = NPROW
      NPCOL  = int(NPROCS / NPROW)
      NPCOLtemp = NPCOL
      NPROCSused = NPROWtemp * NPCOLtemp
 10   CONTINUE
      IF ( NPROWtemp >= NPCOLtemp/FLATNESS .AND. NPROWtemp > 1) THEN
        NPROWtemp = NPROWtemp - 1
        NPCOLtemp = int(NPROCS / NPROWtemp)
        KEEPIT=.FALSE.
        IF ( NPROWtemp * NPCOLtemp .GE. NPROCSused ) THEN
          IF ( ( K50 .NE. 1 .AND. NPROWtemp >= NPCOLtemp/FLATNESS)
     &         .OR. NPROWtemp * NPCOLtemp .GT. NPROCSused )
     &         KEEPIT=.TRUE.
        END IF
        IF ( KEEPIT ) THEN
          NPROW = NPROWtemp
          NPCOL = NPCOLtemp
          NPROCSused = NPROW * NPCOL
        END IF
        GO TO 10
      END IF
      RETURN
      END SUBROUTINE DMUMPS_99
      SUBROUTINE DMUMPS_290(MYID, M, N, ASEQ,
     &                    LOCAL_M, LOCAL_N,
     &                    MBLOCK, NBLOCK,
     &                    APAR,
     &                    MASTER_ROOT,
     &                    NPROW, NPCOL,
     &                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      DOUBLE PRECISION APAR( LOCAL_M, LOCAL_N )
      DOUBLE PRECISION ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, IDEST, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      DOUBLE PRECISION WK( MBLOCK * NBLOCK )
      LOGICAL JUPDATE
        IAPAR = 1
        JAPAR = 1
        DO J = 1, N, NBLOCK
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = mod ( IBLOCK, NPROW )
            ICOL = mod ( JBLOCK, NPCOL )
            IDEST = IROW * NPCOL + ICOL
            IF ( IDEST .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  WK(KK)=ASEQ(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SSEND( WK, SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         IDEST, 128, COMM, IERR )
              ELSE IF ( MYID .EQ. IDEST ) THEN
                CALL MPI_RECV( WK(1),
     &                         SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         MASTER_ROOT,128,COMM,STATUS,IERR)
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  APAR(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     &              JAPAR:JAPAR+SIZE_JBLOCK-1 )
     &        = ASEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
      RETURN
      END SUBROUTINE DMUMPS_290
      SUBROUTINE DMUMPS_156(MYID, M, N, ASEQ,
     &                    LOCAL_M, LOCAL_N,
     &                    MBLOCK, NBLOCK,
     &                    APAR,
     &                    MASTER_ROOT,
     &                    NPROW, NPCOL,
     &                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      DOUBLE PRECISION APAR( LOCAL_M, LOCAL_N )
      DOUBLE PRECISION ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, ISOUR, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      DOUBLE PRECISION WK( MBLOCK * NBLOCK )
      LOGICAL JUPDATE
        IAPAR = 1
        JAPAR = 1
        DO J = 1, N, NBLOCK
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = mod ( IBLOCK, NPROW )
            ICOL = mod ( JBLOCK, NPCOL )
            ISOUR = IROW * NPCOL + ICOL
            IF ( ISOUR .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                CALL MPI_RECV( WK(1), SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         ISOUR, 128, COMM, STATUS, IERR )
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  ASEQ(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
              ELSE IF ( MYID .EQ. ISOUR ) THEN
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  WK(KK)=APAR(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SSEND( WK( 1 ),
     &                         SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         MASTER_ROOT,128,COMM,IERR)
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              ASEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
     &        = APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     &                JAPAR:JAPAR+SIZE_JBLOCK-1 )
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
      RETURN
      END SUBROUTINE DMUMPS_156
      SUBROUTINE DMUMPS_284(root, IROOT, N,
     &                  IW, LIW, A, LA,
     &                  FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &                  LRLU, IPTRLU,
     &                  IWPOS, IWPOSCB, PTRIST, PTRAST,
     &                  STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &                  COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER MYID
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      TYPE (DMUMPS_ROOT_STRUC ) :: root
      INTEGER(8) :: LA, LRLU, IPTRLU, LRLUS
      INTEGER IROOT, LIW, N, IWPOS, IWPOSCB
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER PTRIST(KEEP(28)), STEP(N)
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PIMASTER(KEEP(28))
      INTEGER ITLOC( N + KEEP(253) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER COMP, IFLAG, IERROR
      INCLUDE 'mumps_headers.h'
      INTEGER FILS( N ), PTRAIW(N), PTRARW( N )
      INTEGER INTARR(max(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      INTEGER numroc
      EXTERNAL numroc
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER(8) :: LREQA_ROOT
      INTEGER LREQI_ROOT, LOCAL_M, LOCAL_N, allocok
            LOCAL_M = numroc( root%ROOT_SIZE, root%MBLOCK,
     &                root%MYROW, 0, root%NPROW )
            LOCAL_M = max( 1, LOCAL_M )
            LOCAL_N = numroc( root%ROOT_SIZE, root%NBLOCK,
     &                root%MYCOL, 0, root%NPCOL )
            IF (KEEP(253).GT.0) THEN
              root%RHS_NLOC = numroc( KEEP(253), root%NBLOCK,
     &                root%MYCOL, 0, root%NPCOL )
              root%RHS_NLOC = max(1, root%RHS_NLOC)
            ELSE
              root%RHS_NLOC = 1
            ENDIF
            IF (associated( root%RHS_ROOT) ) 
     &               DEALLOCATE (root%RHS_ROOT)
            ALLOCATE(root%RHS_ROOT(LOCAL_M,root%RHS_NLOC), 
     &                stat=allocok)
            IF ( allocok.GT.0) THEN
              IFLAG=-13
              IERROR = LOCAL_M*root%RHS_NLOC
              RETURN
            ENDIF
            IF (KEEP(253).NE.0) THEN
              root%RHS_ROOT = ZERO
              CALL DMUMPS_760 (  N, FILS, 
     &             root, KEEP, RHS_MUMPS, 
     &             IFLAG, IERROR )
              IF ( IFLAG .LT. 0 ) RETURN
            ENDIF
            IF (KEEP(60) .NE. 0) THEN
              PTRIST(STEP(IROOT)) = -6666666
              RETURN
            ENDIF
            LREQI_ROOT = 2 + KEEP(IXSZ)
            LREQA_ROOT = int(LOCAL_M,8) * int(LOCAL_N,8)
            IF (LREQA_ROOT.EQ.0_8) THEN
              PTRIST(STEP(IROOT)) = -9999999
              RETURN
            ENDIF
            CALL DMUMPS_22(.FALSE.,0_8,.FALSE.,.FALSE.,
     &                     MYID,N,KEEP,KEEP8,IW,LIW,A,LA,
     &                     LRLU, IPTRLU,
     &                     IWPOS, IWPOSCB, PTRIST, PTRAST,
     &                     STEP, PIMASTER, PAMASTER, LREQI_ROOT,
     &                     LREQA_ROOT, IROOT, S_NOTFREE, .TRUE., COMP,
     &                     LRLUS, IFLAG, IERROR
     &           )
            IF ( IFLAG .LT. 0 ) RETURN
            PTRIST  ( STEP(IROOT) ) = IWPOSCB + 1
            PAMASTER( STEP(IROOT) ) = IPTRLU  + 1_8
            IW( IWPOSCB + 1 + KEEP(IXSZ)) = - LOCAL_N
            IW( IWPOSCB + 2 + KEEP(IXSZ)) =   LOCAL_M
      RETURN
      END SUBROUTINE DMUMPS_284
      SUBROUTINE DMUMPS_760
     &           ( N, FILS, root, KEEP, RHS_MUMPS, 
     &             IFLAG, IERROR )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER N, KEEP(500), IFLAG, IERROR
      INTEGER FILS(N)
      TYPE (DMUMPS_ROOT_STRUC ) :: root
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER JCOL, IPOS_ROOT, JPOS_ROOT,
     &        IROW_GRID, JCOL_GRID, ILOCRHS, JLOCRHS,
     &        INODE
      INODE = KEEP(38)
      DO WHILE (INODE.GT.0)
        IPOS_ROOT = root%RG2L_ROW( INODE ) 
        IROW_GRID  = mod( ( IPOS_ROOT - 1 ) / root%MBLOCK, root%NPROW )
        IF (  IROW_GRID .NE. root%MYROW ) GOTO 100 
        ILOCRHS = root%MBLOCK * ( ( IPOS_ROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOS_ROOT - 1, root%MBLOCK ) + 1
        DO JCOL = 1, KEEP(253) 
          JPOS_ROOT = JCOL
          JCOL_GRID  = mod((JPOS_ROOT-1)/root%NBLOCK, root%NPCOL)
          IF (JCOL_GRID.NE.root%MYCOL ) CYCLE
           JLOCRHS = root%NBLOCK * ( ( JPOS_ROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOS_ROOT - 1, root%NBLOCK ) + 1
          root%RHS_ROOT(ILOCRHS, JLOCRHS) =
     &                 RHS_MUMPS(INODE+(JCOL-1)*KEEP(254))
        ENDDO
 100    CONTINUE
        INODE=FILS(INODE)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_760
      INTEGER FUNCTION DMUMPS_IXAMAX(n,x,incx)
      DOUBLE PRECISION x(*)
      integer incx,n
      INTEGER idamax
      DMUMPS_IXAMAX = idamax(n,x,incx)
      return
      END FUNCTION DMUMPS_IXAMAX
      SUBROUTINE DMUMPS_XSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
      CHARACTER          UPLO
      INTEGER            INCX, LDA, N
      DOUBLE PRECISION            ALPHA
      DOUBLE PRECISION            A( LDA, * ), X( * )
      CALL dsyr( UPLO, N, ALPHA, X, INCX, A, LDA )
      RETURN
      END SUBROUTINE DMUMPS_XSYR
