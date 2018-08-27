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
      RECURSIVE SUBROUTINE DMUMPS_210( COMM_LOAD, ASS_IRECV,
     &  BUFR, LBUFR, LBUFR_BYTES,
     &
     &  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     &  NFRONT_PERE, NASS_PERE, NFS4FATHER,LMAP, TROW,
     &  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &  LRLUS, N, IW,
     &  LIW, A, LA,
     &  PTRIST, PTLUST_S, PTRFAC,
     &  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     &  NBFIN, ICNTL, KEEP,KEEP8,
     &  root, OPASSW, OPELIW,
     &  ITLOC, RHS_MUMPS,
     &  FILS, PTRARW, PTRAIW, INTARR, DBLARR, ND, FRERE,
     &  LPTRAR, NELT, FRTPTR, FRTELT, 
     &
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &  )
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC ) :: root
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ICNTL( 40 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER SLAVEF, NBFIN
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS, POSFAC
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER COMP
      INTEGER NSTK( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER INODE_PERE, ISON
      INTEGER NFS4FATHER
      INTEGER NBROWS_ALREADY_SENT
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE( * )
      INTEGER LMAP 
      INTEGER TROW( LMAP )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER NOSLA, I, ISTCHK, ISTCHK_LOC
      INTEGER INDICE_PERE
      INTEGER INDICE_PERE_ARRAY_ARG(1)
      INTEGER PDEST, PDEST_MASTER
      INTEGER NFRONT
      INTEGER(8) :: SIZFR
      INTEGER LDA_SON
      INTEGER NSLSON, NBCOLS, NROW, NROWS_TO_SEND, IROW_SON,
     &        NPIV, NROWS_TO_STACK, II, COLLIST
      INTEGER(8) :: POSROW, SHIFTCB_SON
      INTEGER NBCOLS_EFF
      INTEGER PDEST_MASTER_ISON, IPOS_IN_SLAVE
      LOGICAL DESCLU, SLAVE_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED 
      INTEGER MSGSOU, MSGTAG
      INTEGER LP
      INTEGER ITMP
      LOGICAL SAME_PROC, COMPRESSCB
      LOGICAL IS_ERROR_BROADCASTED, IS_ofType5or6
      INTEGER ITYPE, TYPESPLIT
      INTEGER KEEP253_LOC
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_275, MUMPS_330, MUMPS_810
      EXTERNAL MUMPS_275, MUMPS_330, MUMPS_810
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM
      IS_ERROR_BROADCASTED = .FALSE.
      TYPESPLIT = MUMPS_810(PROCNODE_STEPS(STEP(INODE_PERE)),
     &                  SLAVEF)
      IS_ofType5or6 = ((TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6))
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     &  ' : PB allocation SLAVES_PERE in DMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      IF (NSLAVES_PERE.GT.0) 
     &SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = MUMPS_275( PROCNODE_STEPS(STEP(INODE_PERE)),
     &                 SLAVEF )
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP>0) write(LP,*) MYID,
     &  ' : PB allocation NBROW in DMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 670 
      endif
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP>0) THEN
        write(LP,*) MYID, ' : PB allocation LMAP in DMUMPS_210'
        ENDIF
        IFLAG  =-13
        IERROR = LMAP
        GOTO 680
      endif
      MAP( 1 : LMAP ) = TROW( 1 : LMAP )
      PDEST_MASTER_ISON = MUMPS_275(PROCNODE_STEPS(STEP(ISON)),
     &                    SLAVEF)
      SLAVE_ISON = PDEST_MASTER_ISON .NE. MYID
      IF (SLAVE_ISON) THEN
        DO WHILE ( PTRIST(STEP( ISON )) .EQ. 0 )
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    PDEST_MASTER_ISON, MAITRE_DESC_BANDE,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &    )
          IF ( IFLAG .LT. 0 ) THEN
            IS_ERROR_BROADCASTED = .TRUE.  
            GOTO 670                       
          ENDIF
        END DO
        DO WHILE (
     &     ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &       IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) .OR.
     &     ( KEEP(50) .NE. 0 .AND.
     &       IW( PTRIST(STEP(ISON)) + 6 + KEEP(IXSZ) ) .NE. 0 ) )
          IF ( KEEP(50).eq.0) THEN
            MSGSOU = PDEST_MASTER_ISON
            MSGTAG = BLOC_FACTO
          ELSE
            IF ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &           IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) THEN
              MSGSOU = PDEST_MASTER_ISON
              MSGTAG = BLOC_FACTO_SYM
            ELSE
              MSGSOU = MPI_ANY_SOURCE
              MSGTAG = BLOC_FACTO_SYM_SLAVE
            END IF
          END IF
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MSGSOU, MSGTAG,
     &    STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
          IF ( IFLAG .LT. 0 ) THEN
            IS_ERROR_BROADCASTED = .TRUE.  
            GOTO 670                       
          ENDIF
        END DO
      ENDIF
      IF ( NSLAVES_PERE .EQ. 0 ) THEN
        NBROW( 0 ) = LMAP
      ELSE
        DO I = 0, NSLAVES_PERE
          NBROW( I ) = 0
        END DO
        DO I = 1, LMAP_LOC
          INDICE_PERE = MAP( I )
          CALL MUMPS_47(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          NBROW( NOSLA ) = NBROW( NOSLA ) + 1
        END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM(LMAP_LOC), stat=allocok)
      IF (allocok .GT. 0) THEN
          IF (LP.GT.0) THEN
          write(LP,*) MYID,': PB allocation PERM in DMUMPS_210'
          ENDIF
          IFLAG  =-13
          IERROR = LMAP_LOC
          GOTO 670
      ENDIF
      ISTCHK     = PTRIST(STEP(ISON))
      NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
      KEEP253_LOC   = 0
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          IF (INDICE_PERE > NFRONT_PERE - KEEP(253)) THEN
             KEEP253_LOC = KEEP253_LOC + 1
          ENDIF
          CALL MUMPS_47(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          PERM( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      PDEST_MASTER = SLAVES_PERE(0)
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .EQ. MYID ) THEN 
            NBPROCFILS(STEP(INODE_PERE)) =
     &                 NBPROCFILS(STEP(INODE_PERE)) - 1
            IF ( PDEST .EQ. PDEST_MASTER ) THEN
              NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
            ENDIF
            ISTCHK     = PTRIST(STEP(ISON))
            NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
            NROW       = IW(ISTCHK+2+KEEP(IXSZ))
            NPIV       = IW(ISTCHK+3+KEEP(IXSZ))
            NSLSON     = IW(ISTCHK+5+KEEP(IXSZ))
            NFRONT     = NPIV + NBCOLS
            COMPRESSCB = (IW(ISTCHK+XXS).EQ.S_CB1COMP)
            CALL MUMPS_729(SIZFR, IW(ISTCHK+XXR))
            IF (IW(ISTCHK+XXS).EQ.S_NOLCBCONTIG) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = int(NPIV,8)*int(NROW,8)
            ELSE IF (IW(ISTCHK+XXS).EQ.S_NOLCLEANED) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = 0_8
            ELSE
               LDA_SON     = NFRONT
               SHIFTCB_SON = int(NPIV,8)
            ENDIF
            IF (I == NSLAVES_PERE) THEN
              NROWS_TO_STACK=LMAP_LOC-NBROW(I)+1
            ELSE
              NROWS_TO_STACK=NBROW(I+1)-NBROW(I)
            ENDIF
            IF (PDEST .NE. PDEST_MASTER) THEN
               IF ( KEEP(55) .eq. 0 ) THEN
                 CALL DMUMPS_539
     &           (N, INODE_PERE, IW, LIW,
     &           A, LA, NROWS_TO_STACK, NBCOLS,
     &           OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &           ITLOC, RHS_MUMPS,
     &           FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &           KEEP,KEEP8, MYID )
               ELSE
                 CALL DMUMPS_123(NELT, FRTPTR, FRTELT,
     &           N, INODE_PERE, IW, LIW,
     &           A, LA, NROWS_TO_STACK, NBCOLS, 
     &           OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &           ITLOC, RHS_MUMPS,
     &           FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &           KEEP, KEEP8, MYID )
               ENDIF
            ENDIF
            DO II = 1,NROWS_TO_STACK
              IROW_SON = PERM(NBROW(I)+II-1)
              INDICE_PERE=MAP(IROW_SON)
              CALL MUMPS_47(
     &        KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &        NASS_PERE,
     &        NFRONT_PERE - NASS_PERE,
     &        NSLAVES_PERE,
     &        INDICE_PERE,
     &        NOSLA,
     &        IPOS_IN_SLAVE )
              INDICE_PERE = IPOS_IN_SLAVE
              IF ( COMPRESSCB ) THEN
                IF (NBCOLS - NROW .EQ. 0 ) THEN
                  ITMP = IROW_SON 
                  POSROW = PTRAST(STEP(ISON))+
     &                     int(ITMP,8) * int(ITMP-1,8) / 2_8
                ELSE
                  ITMP = IROW_SON + NBCOLS - NROW
                  POSROW = PTRAST(STEP(ISON))
     &               + int(ITMP,8) * int(ITMP-1,8) / 2_8
     &               - int(NBCOLS-NROW,8) * int(NBCOLS-NROW+1,8)/2_8
                ENDIF
              ELSE 
                POSROW = PTRAST(STEP(ISON)) + SHIFTCB_SON
     &               +int(IROW_SON-1,8)*int(LDA_SON,8)
              ENDIF
              IF (PDEST == PDEST_MASTER) THEN
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
                 IF ((IS_ofType5or6).AND.(KEEP(50).EQ.0)) THEN
                   CALL DMUMPS_39(N, INODE_PERE, IW, LIW, 
     &             A, LA, ISON, NROWS_TO_STACK, NBCOLS_EFF, 
     &             INDICE_PERE_ARRAY_ARG,
     &             A(POSROW), PTLUST_S, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &             )
                   EXIT  
                 ELSE IF ( (KEEP(50).NE.0) .AND. 
     &              (.NOT.COMPRESSCB).AND.(IS_ofType5or6) ) THEN
                   CALL DMUMPS_39(N, INODE_PERE, IW, LIW, 
     &             A, LA, ISON, NROWS_TO_STACK,
     &             NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &             A(POSROW), PTLUST_S, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &)
                   EXIT
                 ELSE
                   CALL DMUMPS_39(N, INODE_PERE, IW, LIW, 
     &             A, LA, ISON, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &             A(POSROW), PTLUST_S, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &)
                 ENDIF
              ELSE
                 ISTCHK  = PTRIST(STEP(ISON))
                 COLLIST = ISTCHK + 6 + KEEP(IXSZ) 
     &                   + IW( ISTCHK + 5 +KEEP(IXSZ)) + NROW + NPIV
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
                 IF ( (IS_ofType5or6) .AND.
     &                 (
     &                  ( KEEP(50).EQ.0) 
     &                    .OR. 
     &                  ( (KEEP(50).NE.0).and. (.NOT.COMPRESSCB) )
     &                 )
     &               ) THEN
                   CALL DMUMPS_40(N, INODE_PERE,
     &             IW, LIW,
     &             A, LA, NROWS_TO_STACK, NBCOLS, 
     &             INDICE_PERE_ARRAY_ARG,
     &             IW( COLLIST ), A(POSROW),
     &             OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &             ITLOC, RHS_MUMPS,
     &             FILS, ICNTL, KEEP,KEEP8,
     &             MYID, IS_ofType5or6, LDA_SON)
                   EXIT
                 ELSE
                   CALL DMUMPS_40(N, INODE_PERE,
     &             IW, LIW,
     &             A, LA, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &             IW( COLLIST ), A(POSROW),
     &             OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &             ITLOC, RHS_MUMPS,
     &             FILS, ICNTL, KEEP,KEEP8,
     &             MYID, IS_ofType5or6, LDA_SON)
                 ENDIF
              ENDIF
            ENDDO
            IF (PDEST.EQ.PDEST_MASTER) THEN 
             IF (KEEP(219).NE.0) THEN
               IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
                  IF (COMPRESSCB) THEN
                    WRITE(*,*) "Error 1 in PARPIV/DMUMPS_210"
                    CALL MUMPS_ABORT()
                  ELSE
                    POSROW = PTRAST(STEP(ISON))+SHIFTCB_SON+
     &                       int(NBROW(1)-1,8)*int(LDA_SON,8)
                  ENDIF
                  CALL DMUMPS_617(NFS4FATHER,IERR)
                  IF (IERR .NE.0) THEN
                    IF (LP .GT. 0) THEN
                      WRITE(LP, *) "MAX_ARRAY allocation failed"
                    ENDIF
                    IFLAG=-13
                    IERROR=NFS4FATHER
                    GOTO 600
                  ENDIF
                  ITMP=-9999
                  IF ( LMAP_LOC-NBROW(1)+1-KEEP253_LOC .NE. 0 ) THEN
                  CALL DMUMPS_618(
     &                 A(POSROW),
     & SIZFR-SHIFTCB_SON-int(NBROW(1)-1,8)*int(LDA_SON,8),
     &                 LDA_SON, LMAP_LOC-NBROW(1)+1-KEEP253_LOC,
     &                 BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,ITMP)
                  ELSE
                       CALL DMUMPS_757(
     &                 BUF_MAX_ARRAY, NFS4FATHER)
                  ENDIF
                  CALL DMUMPS_619(N, INODE_PERE, IW, LIW, 
     &                 A, LA, ISON, NFS4FATHER,
     &                 BUF_MAX_ARRAY, PTLUST_S, PTRAST,
     &                 STEP, PIMASTER,
     &                 OPASSW,IWPOSCB,MYID, KEEP,KEEP8)
               ENDIF
             ENDIF 
             IF ( NBPROCFILS(STEP(ISON)) .EQ. 0) THEN
               ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
               IF (SAME_PROC) THEN
                 CALL DMUMPS_530(N, ISON, INODE_PERE,
     &             IWPOSCB, PIMASTER, PTLUST_S, IW, LIW, STEP,
     &             KEEP,KEEP8)
               ENDIF
               IF (SAME_PROC) THEN
                 ISTCHK_LOC = PTRIST(STEP(ISON))
                 PTRIST(STEP( ISON) ) = -99999999
               ELSE
                 PIMASTER(STEP( ISON )) = -99999999
               ENDIF
               CALL DMUMPS_152(.FALSE., MYID, N,
     &            ISTCHK_LOC,
     &            PAMASTER(STEP(ISON)),
     &            IW, LIW, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &            LA, KEEP,KEEP8, .FALSE.
     &            )
             ENDIF
             IF ( NBPROCFILS(STEP(INODE_PERE)) .EQ. 0 ) THEN
               CALL DMUMPS_507( N, IPOOL, LPOOL,
     &           PROCNODE_STEPS,
     &           SLAVEF, KEEP(28), KEEP(76), KEEP(80),
     &           KEEP(47), STEP, INODE_PERE+N )
               IF (KEEP(47) .GE. 3) THEN
                 CALL DMUMPS_500(
     &          IPOOL, LPOOL, 
     &          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &          MYID, STEP, N, ND, FILS )
               ENDIF
             END IF
            ELSE
             CALL DMUMPS_531
     &       (N, INODE_PERE, IW, LIW,
     &       NBROW(I), STEP, PTRIST, ITLOC, RHS_MUMPS,
     &       KEEP,KEEP8)
            END IF
        END IF
      END DO
      DO I = NSLAVES_PERE, 0, -1
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
            DESCLU = .FALSE.
            NBROWS_ALREADY_SENT = 0
            IF (I == NSLAVES_PERE) THEN
              NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
            ELSE
              NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
            ENDIF
            COMPRESSCB=(IW(PTRIST(STEP(ISON))+XXS).EQ.S_CB1COMP)
 95         CONTINUE
            IF ( PTRIST(STEP(ISON)) .lt.0 .or.
     &         IW ( PTRIST(STEP(ISON) )+KEEP(IXSZ) ) .GT. N ) THEN
              WRITE(*,*) MYID,': Internal error in Maplig'
              WRITE(*,*) MYID,': PTRIST(STEP(ISON))/N=',
     &                           PTRIST(STEP(ISON)), N
              WRITE(*,*) MYID,': I, NBROW(I)=',I, NBROW(I)
              WRITE(*,*) MYID,': NSLAVES_PERE=',NSLAVES_PERE
              WRITE(*,*) MYID,': ISON, INODE_PERE=',ISON,INODE_PERE
              WRITE(*,*) MYID,': Son header=',
     &        IW(PTRIST(STEP(ISON)): PTRIST(STEP(ISON))+5+KEEP(IXSZ))
              CALL MUMPS_ABORT()
            END IF
            CALL DMUMPS_67( NBROWS_ALREADY_SENT,
     &      DESCLU, INODE_PERE,
     &      NFRONT_PERE, NASS_PERE, NFS4FATHER,
     &           NSLAVES_PERE, ISON,
     &      NROWS_TO_SEND, LMAP_LOC, MAP,
     &      PERM(min(LMAP_LOC,NBROW(I))),
     &      IW( PTRIST(STEP(ISON))),
     &      A(PTRAST(STEP(ISON))), I, PDEST, PDEST_MASTER, 
     &      COMM, IERR, 
     &
     &      KEEP,KEEP8, STEP, N, SLAVEF,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, COMPRESSCB,
     &      KEEP253_LOC )
            IF ( IERR .EQ. -2 ) THEN
              IFLAG  = -17
              IF (LP .GT. 0) THEN
                WRITE(LP,*)
     &          "FAILURE: SEND BUFFER TOO SMALL IN DMUMPS_210"
              ENDIF
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+KEEP(IXSZ))
     &        * KEEP( 35 )
              GO TO 600
            END IF
            IF ( IERR .EQ. -3 ) THEN
              IF (LP .GT. 0) THEN
                WRITE(LP,*)
     &          "FAILURE: RECV BUFFER TOO SMALL IN DMUMPS_210"
              ENDIF
              IFLAG  = -20
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+KEEP(IXSZ))
     &        * KEEP( 35 )
              GOTO 600
            ENDIF
            IF (KEEP(219).NE.0) THEN
             IF ( IERR .EQ. -4 ) THEN
               IFLAG  = -13
              IERROR = NFS4FATHER
              IF (LP .GT. 0) THEN
                WRITE(LP, *)
     & "FAILURE: MAX_ARRAY allocation failed IN DMUMPS_210"
              ENDIF
              GO TO 600
             END IF
            END IF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL DMUMPS_329( COMM_LOAD,
     &          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &          MPI_ANY_SOURCE, MPI_ANY_TAG,
     &          STATUS,
     &          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &          IWPOS, IWPOSCB, IPTRLU,
     &          LRLU, LRLUS, N, IW, LIW, A, LA,
     &          PTRIST, PTLUST_S, PTRFAC,
     &          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &          IFLAG, IERROR, COMM,
     &          NBPROCFILS,
     &          IPOOL, LPOOL, LEAF,
     &          NBFIN, MYID, SLAVEF,
     &
     &          root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, 
     &          PTRARW, PTRAIW,
     &          INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     &          NELT, FRTPTR, FRTELT, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
              IF ( IFLAG .LT. 0 ) THEN
                IS_ERROR_BROADCASTED=.TRUE.
                GOTO 600
              ENDIF
              GO TO 95
            END IF
        END IF
      END DO
      ITYPE = MUMPS_330(PROCNODE_STEPS(STEP(ISON)), SLAVEF)
      IF (KEEP(214) .EQ. 2) THEN
        CALL DMUMPS_314( N, ISON,
     &    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE
     &     )
        IF (IFLAG .LT. 0) THEN
          IS_ERROR_BROADCASTED = .TRUE.
          GOTO 600
        ENDIF
      ENDIF
      CALL DMUMPS_626( N, ISON, PTRIST, PTRAST, IW, LIW,
     &             A, LA, LRLU, LRLUS, IWPOSCB, IPTRLU,
     &             STEP, MYID, KEEP
     &)
 600  CONTINUE
      DEALLOCATE(PERM)
 670  CONTINUE
      DEALLOCATE(MAP)
 680  CONTINUE
      DEALLOCATE(NBROW)
      DEALLOCATE(SLAVES_PERE)
 700  CONTINUE
      IF (IFLAG .LT. 0 .AND. .NOT. IS_ERROR_BROADCASTED) THEN
        CALL DMUMPS_44( MYID, SLAVEF, COMM )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_210
      SUBROUTINE DMUMPS_211( COMM_LOAD, ASS_IRECV, 
     &  BUFR, LBUFR, LBUFR_BYTES,
     &
     &  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     &  NFRONT_PERE, NASS_PERE, NFS4FATHER, LMAP, TROW,
     &  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &  LRLUS, N, IW,
     &  LIW, A, LA,
     &  PTRIST, PTLUST_S, PTRFAC,
     &  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     &  NBFIN, ICNTL, KEEP,KEEP8, root,
     &  OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &  FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &  ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     &
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &  )
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER ICNTL( 40 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER SLAVEF, NBFIN
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS, POSFAC
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      DOUBLE PRECISION A( LA )
      INTEGER COMP
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER INODE_PERE, ISON
      INTEGER NFS4FATHER
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE(NSLAVES_PERE)
      INTEGER NELIM, LMAP, TROW( LMAP )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      INTEGER LPTRAR, NELT
      INTEGER IW( LIW )
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK( KEEP(28) ), ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER LP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER NOSLA, I, ISTCHK, ISTCHK_LOC
      INTEGER NBROWS_ALREADY_SENT 
      INTEGER INDICE_PERE
      INTEGER INDICE_PERE_ARRAY_ARG(1)
      INTEGER PDEST, PDEST_MASTER, NFRONT
      LOGICAL SAME_PROC, DESCLU
      INTEGER(8) :: APOS, POSROW, ASIZE
      INTEGER NSLSON, NBCOLS, NROW, NROWS_TO_SEND,
     &        NPIV, NROWS_TO_STACK, II, IROW_SON,
     &        IPOS_IN_SLAVE
      INTEGER NBCOLS_EFF
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      LOGICAL COMPRESSCB
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      if (NSLAVES_PERE.le.0) then
       write(6,*) ' error 2 in maplig_fils_niv1 ', NSLAVES_PERE
       CALL MUMPS_ABORT()
      endif
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      IF (allocok .GT. 0) THEN
        IF (LP > 0)
     &  write(LP,*) MYID,
     &  ' : PB allocation NBROW in DMUMPS_211'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat =allocok)
      IF ( allocok .GT. 0 ) THEN
        IF (LP > 0) write(LP,*) MYID,
     &  ' : PB allocation SLAVES_PERE in DMUMPS_211'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = MUMPS_275( 
     &                       PROCNODE_STEPS(STEP(INODE_PERE)),
     &                       SLAVEF )
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     &   ' : PB allocation LMAP in DMUMPS_211'
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
      MAP( 1 : LMAP_LOC ) = TROW( 1 : LMAP_LOC )
      DO I = 0, NSLAVES_PERE
        NBROW( I ) = 0
      END DO
      IF (NSLAVES_PERE == 0) THEN
        NBROW(0) = LMAP_LOC
      ELSE
       DO I = 1, LMAP_LOC
        INDICE_PERE = MAP( I )
        CALL MUMPS_47(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &         NASS_PERE,
     &         NFRONT_PERE - NASS_PERE,
     &         NSLAVES_PERE,
     &         INDICE_PERE,
     &         NOSLA,
     &         IPOS_IN_SLAVE )
        NBROW( NOSLA ) = NBROW( NOSLA ) + 1
       END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     &  ': PB allocation PERM in DMUMPS_211'
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          CALL MUMPS_47(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          PERM( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      PDEST_MASTER = MYID
      IF ( SLAVES_PERE(0) .NE. MYID ) THEN
        WRITE(*,*) 'Error 1 in MAPLIG_FILS_NIV1:',MYID, SLAVES_PERE
        CALL MUMPS_ABORT()
      END IF
      PDEST        = PDEST_MASTER
        I = 0
        NBPROCFILS(STEP(INODE_PERE))=NBPROCFILS(STEP(INODE_PERE))-1
        NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
        NELIM      = IW(ISTCHK+1+KEEP(IXSZ))
        NROW       = IW(ISTCHK+2+KEEP(IXSZ))
        NPIV       = IW(ISTCHK+3+KEEP(IXSZ))
        IF (NPIV.LT.0) THEN
         write(6,*) ' Error 2 in DMUMPS_211 ', NPIV
         CALL MUMPS_ABORT()
        ENDIF
        NSLSON     = IW(ISTCHK+5+KEEP(IXSZ))
        NFRONT     = NPIV + NBCOLS
        COMPRESSCB=(IW(PTRIST(STEP(ISON))+XXS) .eq. S_CB1COMP)
        IF (I == NSLAVES_PERE) THEN
          NROWS_TO_STACK=LMAP_LOC-NBROW(I)+1
        ELSE
          NROWS_TO_STACK=NBROW(I+1)-NBROW(I)
        ENDIF
        DO II = 1,NROWS_TO_STACK
          IROW_SON=PERM(NBROW(I)+II-1)
          INDICE_PERE = MAP(IROW_SON)
          CALL MUMPS_47(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &         NASS_PERE,
     &         NFRONT_PERE - NASS_PERE,
     &         NSLAVES_PERE,
     &         INDICE_PERE,
     &         NOSLA,
     &         IPOS_IN_SLAVE )
          INDICE_PERE = IPOS_IN_SLAVE
          IF (COMPRESSCB) THEN
            IF (NELIM.EQ.0) THEN
            POSROW = PAMASTER(STEP(ISON)) +
     &         int(IROW_SON,8)*int(IROW_SON-1,8)/2_8
            ELSE
            POSROW = PAMASTER(STEP(ISON)) +
     &         int(NELIM+IROW_SON,8)*int(NELIM+IROW_SON-1,8)/2_8
            ENDIF
          ELSE
            POSROW = PAMASTER(STEP(ISON)) +
     &             int(NELIM+IROW_SON-1,8)*int(NBCOLS,8)
          ENDIF
          IF (KEEP(50).NE.0) THEN
            NBCOLS_EFF = NELIM + IROW_SON
          ELSE
            NBCOLS_EFF = NBCOLS
          ENDIF
          INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
          CALL DMUMPS_39(N, INODE_PERE, IW, LIW, 
     &    A, LA, ISON, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &    A(POSROW), PTLUST_S, PTRAST,
     &    STEP, PIMASTER, OPASSW, IWPOSCB, 
     &    MYID, KEEP,KEEP8,.FALSE.,NBCOLS_EFF)
        ENDDO
        IF (KEEP(219).NE.0) THEN
         IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
           IF (COMPRESSCB) THEN
             POSROW = PAMASTER(STEP(ISON))
     &          + int(NELIM+NBROW(1),8)*int(NELIM+NBROW(1)-1,8)/2_8
             ASIZE  = int(LMAP_LOC+NELIM,8)*int(NELIM+LMAP_LOC+1,8)/2_8
     &          - int(NELIM+NBROW(1),8)*int(NELIM+NBROW(1)-1,8)/2_8
           ELSE
             POSROW = PAMASTER(STEP(ISON)) +
     &                 int(NELIM+NBROW(1)-1,8)*int(NBCOLS,8)
             ASIZE  = int(LMAP_LOC-NBROW(1)+1,8) * int(NBCOLS,8)
           ENDIF
           CALL DMUMPS_617(NFS4FATHER,IERR)
           IF (IERR .NE.0) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &    ": PB allocation MAX_ARRAY during DMUMPS_211"
              IFLAG=-13
              IERROR=NFS4FATHER
              GOTO 700
           ENDIF
           IF  ( LMAP_LOC-NBROW(1)+1-KEEP(253).GT. 0 ) THEN
           CALL DMUMPS_618(
     &          A(POSROW),ASIZE,NBCOLS,
     &          LMAP_LOC-NBROW(1)+1-KEEP(253),
     &          BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,
     &          NELIM+NBROW(1))
           ELSE
                CALL DMUMPS_757(BUF_MAX_ARRAY,
     &          NFS4FATHER)
           ENDIF
           CALL DMUMPS_619(N, INODE_PERE, IW, LIW, 
     &          A, LA, ISON, NFS4FATHER,
     &          BUF_MAX_ARRAY, PTLUST_S, PTRAST,
     &          STEP, PIMASTER, OPASSW,
     &          IWPOSCB,MYID, KEEP,KEEP8)
         ENDIF
        ENDIF 
          IF ( NBPROCFILS(STEP(ISON)) .EQ. 0) THEN
               ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
               IF (SAME_PROC) THEN
                 CALL DMUMPS_530(N, ISON, INODE_PERE,
     &            IWPOSCB, PIMASTER, PTLUST_S, IW, LIW, STEP,
     &            KEEP,KEEP8)
               ENDIF
          ENDIF
          IF ( NBPROCFILS(STEP(INODE_PERE)) .EQ. 0 ) THEN
            CALL DMUMPS_507( N, IPOOL, LPOOL,
     &        PROCNODE_STEPS,
     &        SLAVEF, KEEP(28), KEEP(76), KEEP(80),
     &        KEEP(47), STEP, INODE_PERE+N )
            IF (KEEP(47) .GE. 3) THEN
              CALL DMUMPS_500(
     &       IPOOL, LPOOL, 
     &       PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &       MYID, STEP, N, ND, FILS )
            ENDIF
          END IF
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
           NBROWS_ALREADY_SENT = 0
 95        CONTINUE
           NFRONT = IW(PIMASTER(STEP(ISON))+KEEP(IXSZ))
           NELIM  = IW(PIMASTER(STEP(ISON))+1+KEEP(IXSZ))
           APOS = PAMASTER(STEP(ISON))
           DESCLU = .TRUE.
           IF (I == NSLAVES_PERE) THEN
             NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
           ELSE
             NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
           ENDIF
           CALL DMUMPS_67(NBROWS_ALREADY_SENT,
     &      DESCLU, INODE_PERE,
     &      NFRONT_PERE, NASS_PERE, NFS4FATHER, 
     &           NSLAVES_PERE,
     &      ISON, NROWS_TO_SEND, LMAP_LOC,
     &      MAP, PERM(min(LMAP_LOC,NBROW(I))),
     &      IW(PIMASTER(STEP(ISON))),
     &      A(APOS), I, PDEST, PDEST_MASTER, COMM, IERR,
     &
     &      KEEP,KEEP8, STEP, N, SLAVEF,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &      COMPRESSCB, KEEP(253))
            IF ( IERR .EQ. -2 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING DMUMPS_211"
              IFLAG  = -17
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            END IF
            IF ( IERR .EQ. -3 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING DMUMPS_211"
              IFLAG  = -20
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            ENDIF
            IF (KEEP(219).NE.0) THEN
             IF ( IERR .EQ. -4 ) THEN
               IFLAG  = -13
               IERROR = BUF_LMAX_ARRAY
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, MAX_ARRAY ALLOC FAILED DURING DMUMPS_211"
               GO TO 700
             ENDIF
            ENDIF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .FALSE.
              MESSAGE_RECEIVED = .FALSE.
              CALL DMUMPS_329( COMM_LOAD,
     &          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &          MPI_ANY_SOURCE, MPI_ANY_TAG,
     &          STATUS, 
     &          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &          IWPOS, IWPOSCB, IPTRLU,
     &          LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &          PTLUST_S, PTRFAC,
     &          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &          IFLAG, IERROR, COMM,
     &          NBPROCFILS,
     &          IPOOL, LPOOL, LEAF,
     &          NBFIN, MYID, SLAVEF,
     &          root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &          FILS, PTRARW, PTRAIW,
     &          INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &          LPTRAR, NELT, FRTPTR, FRTELT, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
              IF ( IFLAG .LT. 0 ) GOTO 600
              GO TO 95
            END IF
        END IF
      END DO
      ISTCHK = PTRIST(STEP(ISON))
      PTRIST(STEP( ISON )) = -77777777
            IF ( IW(ISTCHK+KEEP(IXSZ)) .GE. 0 ) THEN
              WRITE(*,*) 'error 3 in DMUMPS_211'
              CALL MUMPS_ABORT()
            ENDIF
      CALL DMUMPS_152(.FALSE., MYID, N, ISTCHK,
     &     PAMASTER(STEP(ISON)),
     &     IW, LIW, LRLU, LRLUS, IPTRLU,
     &     IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &     )
 600  CONTINUE
      DEALLOCATE(NBROW)
      DEALLOCATE(MAP)
      DEALLOCATE(PERM)
      DEALLOCATE(SLAVES_PERE)
      RETURN
 700  CONTINUE
      CALL DMUMPS_44(MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_211
      SUBROUTINE DMUMPS_93(SIZE_INPLACE,
     &MYID,N,IOLDPS,TYPE,IW, LIW, A, LA,
     &POSFAC, LRLU, LRLUS, IWPOS, PTRAST, PTRFAC, STEP, KEEP,KEEP8,
     &SSARBR,INODE,IERR)
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER MYID
      INTEGER IOLDPS, TYPE, LIW, N, KEEP(500)
      INTEGER(8) :: SIZE_INPLACE, LA, POSFAC, LRLU, LRLUS
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) KEEP8(150)
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER IWPOS, LDLT
      INTEGER STEP( N )
      INTEGER (8) :: PTRFAC(KEEP(28))
      LOGICAL SSARBR
      INTEGER IOLDSHIFT, IPSSHIFT
      INCLUDE 'mumps_headers.h'
      INTEGER LCONT, NELIM, NROW, NPIV, INTSIZ
      INTEGER NFRONT, NSLAVES
      INTEGER IPS, IPSIZE
      INTEGER(8) :: SIZELU, SIZECB, IAPOS, I
      LOGICAL MOVEPTRAST
      INTEGER INODE
      INTEGER IERR
      IERR=0
      LDLT = KEEP(50)
      IOLDSHIFT = IOLDPS + KEEP(IXSZ)
      IF ( IW( IOLDSHIFT ) < 0 ) THEN
        write(*,*) ' ERROR 1 compressLU:Should not point to a band.'
        CALL MUMPS_ABORT()
      ELSE IF ( IW( IOLDSHIFT + 2 ) < 0 ) THEN
        write(*,*) ' ERROR 2 compressLU:Stack not performed yet',
     &  IW(IOLDSHIFT + 2)
        CALL MUMPS_ABORT()
      ENDIF
      LCONT  = IW( IOLDSHIFT )
      NELIM  = IW( IOLDSHIFT + 1 )
      NROW   = IW( IOLDSHIFT + 2 )
      NPIV   = IW( IOLDSHIFT + 3 )
      IAPOS  = PTRFAC(IW( IOLDSHIFT + 4 ))
      NSLAVES= IW( IOLDSHIFT + 5 )
      NFRONT = LCONT + NPIV
      INTSIZ = IW(IOLDPS+XXI)
      IF ( (NSLAVES > 0  .AND. TYPE .NE. 2) .OR. 
     &   (NSLAVES .eq. 0 .AND. TYPE .EQ. 2 ) ) THEN
          WRITE(*,*) ' ERROR 3 compressLU: problem with level of inode'
          CALL MUMPS_ABORT()
      END IF
      IF (LDLT.EQ.0) THEN
        SIZELU = int(LCONT + NROW, 8) * int(NPIV,8)
      ELSE
        SIZELU =   int(NROW,8) * int(NPIV,8)
      ENDIF
      IF ( TYPE .EQ. 2 ) THEN
        IF (LDLT.EQ.0) THEN
          SIZECB = int(NELIM,8) * int(LCONT,8)
        ELSE
          IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
            SIZECB = int(NELIM+1,8) * int(NELIM + NPIV,8)
          ELSE
            SIZECB = int(NELIM,8) * int(NELIM + NPIV,8)
          ENDIF
        ENDIF
      ELSE
        IF (LDLT.EQ.0) THEN
         SIZECB = int(LCONT,8) * int(LCONT,8)
        ELSE
         SIZECB = int(NROW,8) * int(LCONT,8)
        ENDIF
      END IF
      CALL MUMPS_724( IW(IOLDPS+XXR), SIZECB )
      IF ((SIZECB.EQ.0_8).AND.(KEEP(201).EQ.0)) THEN
         GOTO 500
      ENDIF
      IF (KEEP(201).EQ.2) THEN
         KEEP8(31)=KEEP8(31)+SIZELU
         CALL DMUMPS_576(INODE,PTRFAC,KEEP,KEEP8,
     &        A,LA,SIZELU, IERR)
         IF(IERR.LT.0)THEN
            WRITE(*,*)MYID,': Internal error in DMUMPS_576'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      IF ( IOLDPS + INTSIZ .NE. IWPOS ) THEN
         IPS = IOLDPS + INTSIZ
         MOVEPTRAST = .FALSE.
         DO WHILE ( IPS .NE. IWPOS )
           IPSIZE = IW(IPS+XXI)
           IPSSHIFT = IPS + KEEP(IXSZ)
           IF ( IW( IPSSHIFT + 2 ) < 0 ) THEN
             NFRONT = IW( IPSSHIFT )
             IF(KEEP(201).EQ.0)THEN
               PTRFAC(IW( IPSSHIFT + 4 )) = 
     &                      PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
             ELSE
               PTRFAC(IW(IPSSHIFT+4))=PTRFAC(IW(IPSSHIFT+4)) -
     &               SIZECB - SIZELU
             ENDIF
             MOVEPTRAST = .TRUE.
             IF(KEEP(201).EQ.0)THEN
               PTRAST(IW(IPSSHIFT+4))=PTRAST(IW(IPSSHIFT+4))-SIZECB
             ELSE
               PTRAST(IW(IPSSHIFT+4))=PTRAST(IW(IPSSHIFT+4))-SIZECB
     &               - SIZELU
             ENDIF
           ELSE IF ( IW( IPSSHIFT ) < 0 ) THEN
             IF(KEEP(201).EQ.0)THEN
               PTRFAC(IW(IPSSHIFT+3)) = PTRFAC(IW(IPSSHIFT+3))-SIZECB
             ELSE
               PTRFAC(IW(IPSSHIFT+3)) = PTRFAC(IW(IPSSHIFT+3))
     &                                  -SIZECB-SIZELU
             ENDIF
           ELSE
             NFRONT = IW( IPSSHIFT ) + IW( IPSSHIFT + 3 )
             IF(KEEP(201).EQ.0)THEN
                PTRFAC(IW( IPSSHIFT + 4 )) = 
     &                    PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
             ELSE
                PTRFAC(IW( IPSSHIFT + 4 )) = 
     &               PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
     &               - SIZELU
             ENDIF
           END IF
           IPS = IPS + IPSIZE
         END DO
         IF ((SIZECB .NE. 0_8).OR.(KEEP(201).NE.0)) THEN
            IF (KEEP(201).NE.0) THEN
               DO I=IAPOS, POSFAC - SIZECB - SIZELU - 1_8
                  A( I ) = A( I + SIZECB + SIZELU)
               END DO
            ELSE
               DO I=IAPOS + SIZELU, POSFAC - SIZECB - 1_8
                  A( I ) = A( I + SIZECB )
               END DO
            ENDIF
         END IF
      ENDIF
      IF (KEEP(201).NE.0) THEN
        POSFAC = POSFAC  - (SIZECB+SIZELU)
        LRLU   = LRLU    + (SIZECB+SIZELU)
        LRLUS  = LRLUS   + (SIZECB+SIZELU) - SIZE_INPLACE
      ELSE
        POSFAC = POSFAC - SIZECB
        LRLU   = LRLU   + SIZECB
        LRLUS  = LRLUS  + SIZECB - SIZE_INPLACE
      ENDIF
 500  CONTINUE
      CALL DMUMPS_471(SSARBR,.FALSE.,
     &         LA-LRLUS,SIZELU,-SIZECB+SIZE_INPLACE,KEEP,KEEP8,LRLU)
      RETURN
      END SUBROUTINE DMUMPS_93
      SUBROUTINE DMUMPS_314( N, ISON, 
     &    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA, 
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, TYPE_SON
     &     )
      USE DMUMPS_OOC
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER(8) :: LA, LRLU, LRLUS, POSFAC, IPTRLU
      INTEGER N, ISON, LIW, IWPOS, IWPOSCB,
     &        COMP, IFLAG, IERROR, SLAVEF, MYID, COMM,
     &        TYPE_SON
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), STEP(N), 
     & PIMASTER(KEEP(28)), IW(LIW)
      INTEGER PTLUST_S(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      DOUBLE PRECISION OPELIW
      DOUBLE PRECISION FLOP1, FLOP1_EFFECTIVE
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: LREQA, POSA, POSALOC, OLDPOS, JJ
      INTEGER  NFRONT, NCOL_L, NROW_L, LREQI, NSLAVES_L,
     &         POSI, I, IROW_L, ICOL_L, LDA_BAND, NASS
      LOGICAL NONEED_TO_COPY_FACTORS
      INTEGER(8) :: LAFAC, LREQA_HEADER
      INTEGER LIWFAC, STRAT, TYPEFile, NextPivDummy,
     &        IOLDPS_CB
      LOGICAL LAST_CALL
      TYPE(IO_BLOCK) :: MonBloc 
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0d0)
      FLOP1 = ZERO
      NCOL_L = IW( PTRIST(STEP( ISON )) + 3 + KEEP(IXSZ) )
      NROW_L = IW( PTRIST(STEP( ISON )) + 2 + KEEP(IXSZ) )
      NSLAVES_L = IW( PTRIST(STEP( ISON )) + 5 + KEEP(IXSZ) )
      LDA_BAND = NCOL_L + IW( PTRIST(STEP( ISON )) + KEEP(IXSZ) )
      IF  ( KEEP(50) .eq. 0 ) THEN
        NFRONT = LDA_BAND
      ELSE
        NFRONT = IW( PTRIST(STEP( ISON )) + 7 + KEEP(IXSZ) )
      END IF
      IF (KEEP(201).EQ.1) THEN 
          IOLDPS_CB = PTRIST(STEP( ISON ))
          CALL MUMPS_729(LAFAC, IW(IOLDPS_CB+XXR))
          LIWFAC    = IW(IOLDPS_CB+XXI)
          TYPEFile  = TYPEF_L
          NextPivDummy      = -8888 
          MonBloc%INODE    = ISON
          MonBloc%MASTER   = .FALSE.   
          MonBloc%Typenode =  2        
          MonBloc%NROW     = NROW_L
          MonBloc%NCOL     = LDA_BAND
          MonBloc%NFS      = IW(IOLDPS_CB+1+KEEP(IXSZ))
          MonBloc%LastPiv  = NCOL_L    
          NULLIFY(MonBloc%INDICES)
          STRAT        = STRAT_WRITE_MAX
          LAST_CALL    = .TRUE.
          MonBloc%Last = .TRUE.
          CALL DMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A(PTRAST(STEP(ISON))), LAFAC, MonBloc,
     &           NextPivDummy, NextPivDummy,
     &           IW(IOLDPS_CB), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG,LAST_CALL )
          IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) THEN 
          ENDIF
      ENDIF  
      NONEED_TO_COPY_FACTORS = (KEEP(201).EQ.1) .OR. (KEEP(201).EQ.-1)
      IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) THEN 
        GOTO 80
      ENDIF
      LREQI   = 4 + NCOL_L + NROW_L + KEEP(IXSZ)
      LREQA_HEADER =  int(NCOL_L,8) * int(NROW_L,8)
      IF (NONEED_TO_COPY_FACTORS) THEN 
        LREQA = 0_8
      ELSE
        LREQA   = LREQA_HEADER
      ENDIF
      IF ( LRLU .LT. LREQA .OR.
     &  IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LREQA ) THEN
          IFLAG  = -9
          CALL MUMPS_731(LREQA - LRLUS, IERROR)
          GO TO 700
        END IF
        CALL DMUMPS_94( N,KEEP(28), IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ))
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB compress Stack_band:LRLU,LRLUS=',
     &            LRLU, LRLUS
                  IFLAG = -9
                  CALL MUMPS_731(LREQA - LRLUS, IERROR)
                  GOTO 700
        END IF
        IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
          IFLAG  = -8
          IERROR = IWPOS + LREQI - 1 - IWPOSCB
          GOTO 700
        END IF
      END IF
      IF (.NOT. NONEED_TO_COPY_FACTORS) THEN
        POSA = POSFAC
        POSFAC = POSFAC + LREQA
        LRLU = LRLU - LREQA
        LRLUS = LRLUS - LREQA
        KEEP8(67) = min(LRLUS, KEEP8(67))
        IF(KEEP(201).NE.2)THEN
           CALL DMUMPS_471(.FALSE.,.FALSE.,
     &          LA-LRLUS,LREQA,LREQA,KEEP,KEEP8,LRLU)
        ELSE
           CALL DMUMPS_471(.FALSE.,.FALSE.,
     &          LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLU)
        ENDIF
      ENDIF
      POSI = IWPOS
      IWPOS = IWPOS + LREQI
      PTLUST_S(STEP( ISON )) = POSI
      IW(POSI+XXI)=LREQI
      CALL MUMPS_730(LREQA, IW(POSI+XXR))
      CALL MUMPS_730(LREQA_HEADER, IW(POSI+XXR))
      IW(POSI+XXS)=-9999
      POSI=POSI+KEEP(IXSZ)
      IW( POSI     ) = - NCOL_L
      IW( POSI + 1 ) =   NROW_L
      IW( POSI + 2 ) =   NFRONT - NCOL_L
      IW( POSI + 3 ) =   STEP(ISON)
      IF (.NOT. NONEED_TO_COPY_FACTORS) THEN
        PTRFAC(STEP(ISON)) = POSA
      ELSE
        PTRFAC(STEP(ISON)) = -77777_8
      ENDIF
      IROW_L = PTRIST(STEP(ISON)) + 6 + NSLAVES_L + KEEP(IXSZ)
      ICOL_L = PTRIST(STEP(ISON)) + 6 + NROW_L + NSLAVES_L + KEEP(IXSZ)
      DO I = 1, NROW_L
        IW( POSI+3+I ) = IW( IROW_L+I-1 )
      ENDDO
      DO I = 1, NCOL_L
        IW( POSI+NROW_L+3+I) = IW( ICOL_L+I-1 )
      ENDDO
      IF (.NOT.NONEED_TO_COPY_FACTORS) THEN
        POSALOC = POSA
        DO I = 1, NROW_L
          OLDPOS =  PTRAST( STEP(ISON)) + int(I-1,8)*int(LDA_BAND,8)
          DO JJ = 0_8, int(NCOL_L-1,8)
            A( POSALOC+JJ ) = A( OLDPOS+JJ )
          ENDDO
          POSALOC = POSALOC + int(NCOL_L,8)
        END DO
      ENDIF
      IF (KEEP(201).EQ.2) THEN
       KEEP8(31)=KEEP8(31)+LREQA
      ENDIF
      KEEP8(10) = KEEP8(10) + int(NCOL_L,8) * int(NROW_L,8)
      IF (KEEP(201).EQ.2) THEN 
       CALL DMUMPS_576(ISON,PTRFAC,KEEP,KEEP8,A,LA,LREQA,IFLAG)
       IF(IFLAG.LT.0)THEN
         WRITE(*,*)MYID,': Internal error in DMUMPS_576'
         IERROR=0
         GOTO 700
       ENDIF
      ENDIF
      IF (KEEP(201).EQ.2) THEN
        POSFAC = POSFAC - LREQA
        LRLU = LRLU + LREQA
        LRLUS = LRLUS + LREQA
        CALL DMUMPS_471(.FALSE.,.FALSE.,
     &            LA-LRLUS,LREQA,0_8,KEEP,KEEP8,LRLU)
      ENDIF
  80  CONTINUE
      IF (TYPE_SON == 1) THEN
         GOTO 90
      ENDIF
      IF ( KEEP(50) .eq. 0 ) THEN
         FLOP1 = dble( NCOL_L * NROW_L) +
     &     dble(NROW_L*NCOL_L)*dble(2*NFRONT-NCOL_L-1)
      ELSE
         FLOP1 = dble( NCOL_L ) * dble( NROW_L )
     &         * dble( 2 * LDA_BAND - NROW_L - NCOL_L + 1)
      END IF
      OPELIW = OPELIW + FLOP1
      FLOP1_EFFECTIVE = FLOP1
      NASS = IW( PTRIST(STEP( ISON )) + 4 + KEEP(IXSZ) )
      IF ( NCOL_L .NE. NASS ) THEN
        IF ( KEEP(50).eq.0 ) THEN
           FLOP1 = dble( NASS * NROW_L) +
     &     dble(NROW_L*NASS)*dble(2*NFRONT-NASS-1)
        ELSE
           FLOP1 = dble( NASS ) * dble( NROW_L ) *
     &     dble( 2 * LDA_BAND - NROW_L - NASS + 1)
        END IF
      END IF
      CALL DMUMPS_190(1,.FALSE.,FLOP1_EFFECTIVE-FLOP1,
     &                        KEEP,KEEP8)
      CALL DMUMPS_190(2,.FALSE.,-FLOP1,KEEP,KEEP8)
 90   CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_314
      SUBROUTINE DMUMPS_626( N, ISON, 
     &    PTRIST, PTRAST, IW, LIW, A, LA, 
     &    LRLU, LRLUS, IWPOSCB,
     &    IPTRLU, STEP, MYID, KEEP
     &     )
      IMPLICIT NONE
      include 'mumps_headers.h'
      INTEGER(8) :: LRLU, LRLUS, IPTRLU, LA
      INTEGER ISON, MYID, N, IWPOSCB
      INTEGER KEEP(500), STEP(N)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER PTRIST(KEEP(28))
      INTEGER LIW
      INTEGER IW(LIW)
      DOUBLE PRECISION A(LA)
      INTEGER ISTCHK
      ISTCHK = PTRIST(STEP(ISON))
      CALL DMUMPS_152(.FALSE.,MYID, N, ISTCHK,
     &     PTRAST(STEP(ISON)),
     &     IW, LIW, LRLU, LRLUS, IPTRLU,
     &     IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &     )
      PTRIST(STEP( ISON )) = -9999888
      PTRAST(STEP( ISON )) = -9999888_8
      RETURN
      END SUBROUTINE DMUMPS_626
      SUBROUTINE DMUMPS_214( KEEP,KEEP8,
     &           MYID, N, NELT, LNA, NZ, NA_ELT, NSLAVES,
     &           MEMORY_MBYTES, EFF, OOC_STRAT, PERLU_ON,
     &           MEMORY_BYTES )
      IMPLICIT NONE
      LOGICAL,   INTENT(IN)  :: EFF, PERLU_ON
      INTEGER,   INTENT(IN)  :: OOC_STRAT
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER MYID, N, NELT, NSLAVES, LNA, NZ, NA_ELT
      INTEGER(8), INTENT(OUT) :: MEMORY_BYTES
      INTEGER,   INTENT(OUT) :: MEMORY_MBYTES
      LOGICAL    :: I_AM_SLAVE, I_AM_MASTER
      INTEGER    :: PERLU, NBRECORDS
      INTEGER(8) :: NB_REAL, MAXS_MIN
      INTEGER(8) :: TEMP, NB_BYTES, NB_INT
      INTEGER    :: DMUMPS_LBUF_INT, DMUMPS_LBUFR_BYTES, DMUMPS_LBUF
      INTEGER    :: NBUFS
      INTEGER(8) :: TEMPI
      INTEGER(8) :: TEMPR
      INTEGER    :: MIN_PERLU
      INTEGER(8) :: BUF_OOC, BUF_OOC_PANEL, BUF_OOC_NOPANEL
      INTEGER(8) :: OOC_NB_FILE_TYPE
      INTEGER(8) :: NSTEPS8, N8, NELT8
      INTEGER(8) :: I8OVERI
      I8OVERI   = int(KEEP(10),8)
      PERLU     = KEEP(12)
      NSTEPS8   = int(KEEP(28),8)
      N8        = int(N,8)
      NELT8     = int(NELT,8)
      IF (.NOT.PERLU_ON) PERLU = 0
      I_AM_MASTER = ( MYID .eq. 0 )
      I_AM_SLAVE  = ( KEEP(46).eq. 1 .or. MYID .ne. 0 )
      TEMP    = 0_8
      NB_REAL = 0_8
      NB_BYTES = 0_8
      NB_INT  = 0_8
      NB_INT = NB_INT + 5_8 * NSTEPS8
      NB_INT = NB_INT + NSTEPS8 + int(KEEP(56),8)*int(NSLAVES+2,8)
      NB_INT = NB_INT + 3_8 * N8
      IF (KEEP(23).ne.0 .and. I_AM_MASTER) NB_INT=NB_INT + N8
      IF (KEEP(55).eq.0) THEN
        NB_INT = NB_INT + 2_8 * N8
      ELSE
        NB_INT = NB_INT + 2_8 * ( NELT8 + 1_8 )
      ENDIF
      IF (KEEP(55) .ne. 0 ) THEN
        NB_INT = NB_INT + N8 + 1_8 + NELT8
      END IF
      NB_INT = NB_INT + int(LNA,8)
      IF ( OOC_STRAT .GT. 0 .OR. OOC_STRAT .EQ. -1 ) THEN
        MAXS_MIN = KEEP8(14)
      ELSE
        MAXS_MIN = KEEP8(12)
      ENDIF
      IF ( .NOT. EFF ) THEN
        IF ( KEEP8(24).EQ.0_8 ) THEN
         NB_REAL = NB_REAL + MAXS_MIN +
     &             int(PERLU,8)*(MAXS_MIN / 100_8 + 1_8 )
        ENDIF
      ELSE
        NB_REAL = NB_REAL + KEEP8(67)
      ENDIF
      IF ( OOC_STRAT .GT. 0 .AND. I_AM_SLAVE ) THEN
        BUF_OOC_NOPANEL = 2_8 * KEEP8(119)
        IF (KEEP(50).EQ.0)THEN
          BUF_OOC_PANEL = 8_8 * int(KEEP(226),8)
        ELSE
          BUF_OOC_PANEL = 4_8 * int(KEEP(226),8)
        ENDIF
        IF (OOC_STRAT .EQ. 2) THEN
          BUF_OOC = BUF_OOC_NOPANEL
        ELSE
          BUF_OOC = BUF_OOC_PANEL
        ENDIF
        NB_REAL = NB_REAL + min(BUF_OOC + int(max(PERLU,0),8) *
     &          (BUF_OOC/100_8+1_8),12000000_8)
        IF (OOC_STRAT .EQ. 2) THEN
          OOC_NB_FILE_TYPE = 1_8
        ELSE
          IF (KEEP(50).EQ.0) THEN
            OOC_NB_FILE_TYPE = 2_8
          ELSE
            OOC_NB_FILE_TYPE = 1_8
          ENDIF
        ENDIF
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8 * I8OVERI
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8 * I8OVERI
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8
      ENDIF
      NB_REAL = NB_REAL + int(KEEP(13),8)
      IF (KEEP(252).EQ.1 .AND. .NOT. I_AM_MASTER) THEN
        NB_REAL = NB_REAL + N8
      ENDIF
      IF ( .not. ( I_AM_SLAVE .and. I_AM_MASTER .and. KEEP(52) .eq. 0
     &         .and. KEEP(55) .ne. 0 ) ) THEN
        NB_INT  = NB_INT  + int(KEEP(14),8)
      END IF
      IF ( I_AM_SLAVE .and. KEEP(38) .ne. 0 ) THEN
        NB_INT = NB_INT + 2_8 * N8
      END IF
      TEMPI= 0_8
      TEMPR = 0_8
      NBRECORDS = KEEP(39)
      IF (KEEP(55).eq.0) THEN
        NBRECORDS = min(KEEP(39), NZ)
      ELSE
        NBRECORDS = min(KEEP(39), NA_ELT)
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
        IF ( I_AM_MASTER ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            NBUFS = NSLAVES 
          ELSE
            NBUFS = NSLAVES - 1
            IF (KEEP(55) .eq. 0 )
     &      TEMPI = TEMPI + 2_8 * N8
          END IF
          TEMPI = TEMPI + 2_8 * int(NBRECORDS,8) * int(NBUFS,8)
          TEMPR = TEMPR + int(NBRECORDS,8) * int(NBUFS,8)
        ELSE
          IF ( KEEP(55) .eq. 0 )THEN
            TEMPI = TEMPI + 2_8 * int(NBRECORDS,8)
            TEMPR = TEMPR + int(NBRECORDS,8)
          END IF
        END IF
      ELSE
        IF ( I_AM_SLAVE ) THEN
          TEMPI = TEMPI + int(1+4*NSLAVES,8) * int(NBRECORDS,8)
          TEMPR = TEMPR + int(1+2*NSLAVES,8) * int(NBRECORDS,8)
        END IF
      END IF
      TEMP = max( NB_BYTES + (NB_INT + TEMPI) * int(KEEP(34),8)
     &           + (NB_REAL+TEMPR) * int(KEEP(35),8)
     &            , TEMP )
      IF ( I_AM_SLAVE ) THEN
        DMUMPS_LBUFR_BYTES = KEEP( 44 ) * KEEP( 35 )
        DMUMPS_LBUFR_BYTES = max( DMUMPS_LBUFR_BYTES,
     &                      100000 )
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU=2
        ELSE
          MIN_PERLU=0
        ENDIF
        DMUMPS_LBUFR_BYTES = DMUMPS_LBUFR_BYTES
     &        + int( 2.0D0 * dble(max(PERLU,MIN_PERLU))*
     &        dble(DMUMPS_LBUFR_BYTES)/100D0)
        NB_BYTES = NB_BYTES + int(DMUMPS_LBUFR_BYTES,8)
        DMUMPS_LBUF = int( dble(KEEP(213)) / 100.0D0
     &                     * dble(KEEP( 43 ) * KEEP( 35 )) )
        DMUMPS_LBUF = max( DMUMPS_LBUF, 100000 )
        DMUMPS_LBUF = DMUMPS_LBUF
     &                 + int( 2.0D0 * dble(max(PERLU,0))*
     &                   dble(DMUMPS_LBUF)/100D0)
        DMUMPS_LBUF = max(DMUMPS_LBUF, DMUMPS_LBUFR_BYTES)
        NB_BYTES = NB_BYTES + int(DMUMPS_LBUF,8)
        DMUMPS_LBUF_INT = ( KEEP(56) + 
     &         NSLAVES * NSLAVES ) * 5
     &               * KEEP(34)
        NB_BYTES = NB_BYTES + int(DMUMPS_LBUF_INT,8)
        IF ( EFF ) THEN
          IF (OOC_STRAT .GT. 0) THEN
            NB_INT = NB_INT + int(KEEP(225),8)
          ELSE
            NB_INT = NB_INT + int(KEEP(15),8)
          ENDIF
        ELSE
          IF (OOC_STRAT .GT. 0) THEN
            NB_INT = NB_INT +  int(
     &           KEEP(225) + 2 * max(PERLU,10) *
     &           ( KEEP(225) / 100 + 1 )
     &                              ,8)
          ELSE
            NB_INT = NB_INT +  int(
     &           KEEP(15) + 2 * max(PERLU,10) *
     &           ( KEEP(15) / 100 + 1 )
     &                              ,8)
          ENDIF
        ENDIF
        NB_INT = NB_INT + NSTEPS8
        NB_INT = NB_INT + NSTEPS8 * I8OVERI
        NB_INT = NB_INT + N8 + 5_8 * NSTEPS8 + 3_8
        NB_INT = NB_INT + 2_8 * NSTEPS8 * I8OVERI
      END IF
      MEMORY_BYTES = NB_BYTES + NB_INT * int(KEEP(34),8) +
     &               NB_REAL * int(KEEP(35),8)
      MEMORY_BYTES = max( MEMORY_BYTES, TEMP )
      MEMORY_MBYTES = int( MEMORY_BYTES / 1000000_8 )  + 1
      RETURN
      END SUBROUTINE DMUMPS_214
      SUBROUTINE DMUMPS_757(M_ARRAY, M_SIZE)
      IMPLICIT NONE
      INTEGER M_SIZE
      DOUBLE PRECISION M_ARRAY(M_SIZE)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      M_ARRAY=ZERO
      RETURN
      END SUBROUTINE DMUMPS_757
      SUBROUTINE DMUMPS_618(
     &     A,ASIZE,NCOL,NROW,
     &     M_ARRAY,NMAX,COMPRESSCB,LROW1)
      IMPLICIT NONE
      INTEGER(8) :: ASIZE
      INTEGER NROW,NCOL,NMAX,LROW1
      LOGICAL COMPRESSCB
      DOUBLE PRECISION A(ASIZE)
      DOUBLE PRECISION M_ARRAY(NMAX)
      INTEGER I
      INTEGER(8):: APOS, J, LROW
      DOUBLE PRECISION ZERO,TMP
      PARAMETER (ZERO=0.0D0)
      M_ARRAY(1:NMAX) = ZERO
      APOS = 0_8
      IF (COMPRESSCB) THEN
        LROW=int(LROW1,8)
      ELSE
        LROW=int(NCOL,8)
      ENDIF
      DO I=1,NROW
         DO J=1_8,int(NMAX,8)
            TMP = abs(A(APOS+J))
            IF(TMP.GT.M_ARRAY(J)) M_ARRAY(J) = TMP
         ENDDO
         APOS = APOS + LROW
         IF (COMPRESSCB) LROW=LROW+1_8
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_618
      SUBROUTINE DMUMPS_710 (id, NB_INT,NB_CMPLX )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC) :: id
      INTEGER(8) NB_INT, NB_CMPLX 
      INTEGER(8) NB_REAL
      NB_INT        = 0_8
      NB_CMPLX      = 0_8
      NB_REAL       = 0_8
      IF (associated(id%IS))          NB_INT=NB_INT+size(id%IS)
      IF (associated(id%IS1))         NB_INT=NB_INT+size(id%IS1)
      NB_INT=NB_INT+size(id%KEEP)
      NB_INT=NB_INT+size(id%ICNTL)
      NB_INT=NB_INT+size(id%INFO)
      NB_INT=NB_INT+size(id%INFOG)
      IF (associated(id%MAPPING))     NB_INT=NB_INT+size(id%MAPPING)
      IF (associated(id%POIDS))       NB_INT=NB_INT+size(id%POIDS)
      IF (associated(id%BUFR))        NB_INT=NB_INT+size(id%BUFR)
      IF (associated(id%STEP))        NB_INT=NB_INT+size(id%STEP)
      IF (associated(id%NE_STEPS  ))  NB_INT=NB_INT+size(id%NE_STEPS  )
      IF (associated(id%ND_STEPS))    NB_INT=NB_INT+size(id%ND_STEPS)
      IF (associated(id%Step2node))   NB_INT=NB_INT+size(id%Step2node)
      IF (associated(id%FRERE_STEPS)) NB_INT=NB_INT+size(id%FRERE_STEPS)
      IF (associated(id%DAD_STEPS))   NB_INT=NB_INT+size(id%DAD_STEPS)
      IF (associated(id%FILS))        NB_INT=NB_INT+size(id%FILS)
      IF (associated(id%PTRAR))       NB_INT=NB_INT+size(id%PTRAR)
      IF (associated(id%FRTPTR))      NB_INT=NB_INT+size(id%FRTPTR)
      NB_INT=NB_INT+size(id%KEEP8) * id%KEEP(10)
      IF (associated(id%PTRFAC)) NB_INT=NB_INT+size(id%PTRFAC) *
     &                                         id%KEEP(10)
      IF (associated(id%FRTELT))      NB_INT=NB_INT+size(id%FRTELT)
      IF (associated(id%NA))          NB_INT=NB_INT+size(id%NA)
      IF       (associated(id%PROCNODE_STEPS))
     &  NB_INT=NB_INT+size(id%PROCNODE_STEPS)
      IF (associated(id%PTLUST_S)) NB_INT=NB_INT+size(id%PTLUST_S)
      IF (associated(id%PROCNODE)) NB_INT=NB_INT+size(id%PROCNODE)
      IF (associated(id%INTARR)) NB_INT=NB_INT+size(id%INTARR)
      IF (associated(id%ELTPROC))  NB_INT=NB_INT+size(id%ELTPROC)
      IF (associated(id%CANDIDATES))
     &  NB_INT=NB_INT+size(id%CANDIDATES)
      IF       (associated(id%ISTEP_TO_INIV2))
     &  NB_INT=NB_INT+size(id%ISTEP_TO_INIV2)
      IF       (associated(id%FUTURE_NIV2))
     &  NB_INT=NB_INT+size(id%FUTURE_NIV2)
      IF (associated(id%TAB_POS_IN_PERE))
     &  NB_INT=NB_INT+size(id%TAB_POS_IN_PERE)
      IF (associated(id%I_AM_CAND))
     &  NB_INT=NB_INT+size(id%I_AM_CAND)
      IF (associated(id%MEM_DIST)) 
     &  NB_INT=NB_INT+size(id%MEM_DIST)
      IF (associated(id%POSINRHSCOMP))
     &  NB_INT=NB_INT+size(id%POSINRHSCOMP)
      IF       (associated(id%MEM_SUBTREE))
     &  NB_INT=NB_INT+size(id%MEM_SUBTREE)
      IF       (associated(id%MY_ROOT_SBTR))
     &  NB_INT=NB_INT+size(id%MY_ROOT_SBTR)
      IF       (associated(id%MY_FIRST_LEAF))
     &  NB_INT=NB_INT+size(id%MY_FIRST_LEAF)
      IF (associated(id%MY_NB_LEAF)) NB_INT=NB_INT+size(id%MY_NB_LEAF)
      IF (associated(id%DEPTH_FIRST)) NB_INT=NB_INT+size(id%DEPTH_FIRST)
      IF (associated(id%COST_TRAV)) NB_INT=NB_INT+size(id%COST_TRAV)
      IF (associated(id%CB_SON_SIZE)) NB_INT=NB_INT+size(id%CB_SON_SIZE)
      IF       (associated(id%OOC_INODE_SEQUENCE))
     &  NB_INT=NB_INT+size(id%OOC_INODE_SEQUENCE)
      IF       (associated(id%OOC_SIZE_OF_BLOCK))
     &  NB_INT=NB_INT+size(id%OOC_SIZE_OF_BLOCK)
      IF       (associated(id%OOC_VADDR)) 
     &  NB_INT=NB_INT+size(id%OOC_VADDR)
      IF       (associated(id%OOC_TOTAL_NB_NODES))
     &  NB_INT=NB_INT+size(id%OOC_TOTAL_NB_NODES)
      IF       (associated(id%OOC_NB_FILES))
     &  NB_INT=NB_INT+size(id%OOC_NB_FILES)
      IF       (associated(id%OOC_FILE_NAME_LENGTH))
     &  NB_INT=NB_INT+size(id%OOC_FILE_NAME_LENGTH)
      IF (associated(id%PIVNUL_LIST)) NB_INT=NB_INT+size(id%PIVNUL_LIST)
      IF (associated(id%SUP_PROC))    NB_INT=NB_INT+size(id%SUP_PROC)
      IF (associated(id%DBLARR))  NB_CMPLX=NB_CMPLX+size(id%DBLARR)
      IF (associated(id%RHSCOMP)) NB_CMPLX=NB_CMPLX+size(id%RHSCOMP)
      IF (associated(id%S))       NB_CMPLX=NB_CMPLX+id%KEEP8(23)
      IF (associated(id%COLSCA))  NB_REAL=NB_REAL+size(id%COLSCA)
      IF (associated(id%ROWSCA))  NB_REAL=NB_REAL+size(id%ROWSCA)
      NB_REAL=NB_REAL+size(id%CNTL)
      NB_REAL=NB_REAL+size(id%RINFO)
      NB_REAL=NB_REAL+size(id%RINFOG)
      NB_REAL=NB_REAL+size(id%DKEEP)
      NB_CMPLX = NB_CMPLX + NB_REAL
      RETURN
      END SUBROUTINE DMUMPS_710 
      SUBROUTINE DMUMPS_756(N8,SRC,DEST)
      IMPLICIT NONE
      INTEGER(8) :: N8
      DOUBLE PRECISION, intent(in)  :: SRC(N8)
      DOUBLE PRECISION, intent(out) :: DEST(N8)
      INTEGER(8) :: SHIFT8, HUG8
      INTEGER    :: I, I4SIZE
      HUG8=int(huge(I4SIZE),8)
      DO I = 1, int(( N8 + HUG8 - 1_8 ) / HUG8)
        SHIFT8 = 1_8 + int(I-1,8) * HUG8
        I4SIZE = int(min(HUG8, N8-SHIFT8+1_8))
        CALL dcopy(I4SIZE, SRC(SHIFT8), 1, DEST(SHIFT8), 1)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_756
      SUBROUTINE DMUMPS_22( INPLACE, MIN_SPACE_IN_PLACE,
     &   SSARBR, PROCESS_BANDE,
     &   MYID,N, KEEP,KEEP8,
     &   IW, LIW, A, LA,
     &   LRLU, IPTRLU,IWPOS,IWPOSCB,
     &   PTRIST,PTRAST,STEP,PIMASTER,PAMASTER,
     &   LREQ, LREQCB, NODE_ARG, STATE_ARG, SET_HEADER,
     &   COMP, LRLUS, IFLAG, IERROR )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N,LIW, KEEP(500)
      INTEGER(8) LA, LRLU, IPTRLU, LRLUS, LREQCB
      INTEGER(8) PAMASTER(KEEP(28)), PTRAST(KEEP(28))
      INTEGER IWPOS,IWPOSCB
      INTEGER(8) :: MIN_SPACE_IN_PLACE
      INTEGER NODE_ARG, STATE_ARG
      INTEGER(8) KEEP8(150)
      INTEGER IW(LIW),PTRIST(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER MYID, IXXP
      DOUBLE PRECISION A(LA)
      LOGICAL INPLACE, PROCESS_BANDE, SSARBR, SET_HEADER
      INTEGER COMP, LREQ, IFLAG, IERROR
      INCLUDE 'mumps_headers.h'
      INTEGER INODE_LOC,NPIV,NASS,NROW,NCB
      INTEGER ISIZEHOLE
      INTEGER(8) :: MEM_GAIN, RSIZEHOLE, LREQCB_EFF, LREQCB_WISHED
      LOGICAL DONE
      IF ( INPLACE ) THEN
        LREQCB_EFF = MIN_SPACE_IN_PLACE
        IF ( MIN_SPACE_IN_PLACE > 0_8 ) THEN
          LREQCB_WISHED = LREQCB
        ELSE
          LREQCB_WISHED = 0_8
        ENDIF
      ELSE
        LREQCB_EFF = LREQCB
        LREQCB_WISHED = LREQCB
      ENDIF
      IF (IWPOSCB.EQ.LIW) THEN
        IF (LREQ.NE.KEEP(IXSZ).OR.LREQCB.NE.0_8
     &      .OR. .NOT. SET_HEADER) THEN
          WRITE(*,*) "Internal error in DMUMPS_22",
     &      SET_HEADER, LREQ, LREQCB
          CALL MUMPS_ABORT()
        ENDIF
        IF (IWPOSCB-IWPOS+1 .LT. KEEP(IXSZ)) THEN
          WRITE(*,*) "Problem with integer stack size",IWPOSCB,
     &               IWPOS, KEEP(IXSZ)
          IFLAG  = -8
          IERROR = LREQ
          RETURN
        ENDIF
        IWPOSCB=IWPOSCB-KEEP(IXSZ)
        IW(IWPOSCB+1+XXI)=KEEP(IXSZ)
        CALL MUMPS_730(0_8,IW(IWPOSCB+1+XXR))
        IW(IWPOSCB+1+XXN)=-919191
        IW(IWPOSCB+1+XXS)=S_NOTFREE
        IW(IWPOSCB+1+XXP)=TOP_OF_STACK
        RETURN
      ENDIF
      IF (KEEP(214).EQ.1.AND.
     &    KEEP(216).EQ.1.AND.
     &    IWPOSCB.NE.LIW) THEN
       IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG.OR.
     &     IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG38) THEN
        NCB  = IW( IWPOSCB+1 + KEEP(IXSZ) )
        NROW = IW( IWPOSCB+1 + KEEP(IXSZ) + 2)
        NPIV = IW( IWPOSCB+1 + KEEP(IXSZ) + 3)
        INODE_LOC= IW( IWPOSCB+1 + XXN)
        CALL DMUMPS_632(IWPOSCB+1,IW,LIW,
     &                          ISIZEHOLE,RSIZEHOLE)
        IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG) THEN
          CALL DMUMPS_627(A,LA,IPTRLU+1_8,
     &                           NROW,NCB,NPIV+NCB,0,
     &                           IW(IWPOSCB+1 + XXS),RSIZEHOLE)
          IW(IWPOSCB+1 + XXS) =S_NOLCLEANED
          MEM_GAIN            = int(NROW,8)*int(NPIV,8)
        ENDIF
        IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG38) THEN
          NASS = IW( IWPOSCB+1 + KEEP(IXSZ) + 4)
          CALL DMUMPS_627(A,LA,IPTRLU+1_8,
     &                           NROW,NCB,NPIV+NCB,NASS-NPIV,
     &                           IW(IWPOSCB+1 + XXS),RSIZEHOLE)
          IW(IWPOSCB+1 + XXS) =S_NOLCLEANED38
          MEM_GAIN = int(NROW,8)*int(NPIV+NCB-(NASS-NPIV),8)
        ENDIF
        IF (ISIZEHOLE.NE.0) THEN
          CALL DMUMPS_630( IW,LIW,IWPOSCB+1,
     &                       IWPOSCB+IW(IWPOSCB+1+XXI),
     &                       ISIZEHOLE )
          IWPOSCB=IWPOSCB+ISIZEHOLE
          IW(IWPOSCB+1+XXP+IW(IWPOSCB+1+XXI))=IWPOSCB+1
          PTRIST(STEP(INODE_LOC))=PTRIST(STEP(INODE_LOC))+
     &    ISIZEHOLE
        ENDIF
        CALL MUMPS_724(IW(IWPOSCB+1+XXR), MEM_GAIN)
        IPTRLU              = IPTRLU+MEM_GAIN+RSIZEHOLE
        LRLU                = LRLU+MEM_GAIN+RSIZEHOLE
        PTRAST(STEP(INODE_LOC))=
     &  PTRAST(STEP(INODE_LOC))+MEM_GAIN+RSIZEHOLE
       ENDIF
      ENDIF
      DONE =.FALSE.
      IF ((IPTRLU.LT.LREQCB_WISHED).OR.(LRLU.LT.LREQCB_WISHED)) THEN
        IF (LRLUS.LT.LREQCB_EFF) THEN
          GOTO 620
        ELSE
          CALL DMUMPS_94(N,KEEP(28),IW,LIW,A,LA,
     &                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     &                    PTRIST,PTRAST,
     &                    STEP, PIMASTER,PAMASTER,KEEP(216),LRLUS,
     &                    KEEP(IXSZ))
          IF ( LRLU .NE. LRLUS ) THEN
            WRITE(*,*) 'PB compress... alloc_cb',
     &      'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 620
          END IF
          DONE = .TRUE.
          COMP = COMP + 1
        ENDIF
      ENDIF
      IF (IWPOSCB-IWPOS+1 .LT. LREQ) THEN
       IF (DONE) GOTO 600
                 CALL DMUMPS_94(N,KEEP(28),IW,LIW,A,LA,
     &                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     &                    PTRIST,PTRAST,
     &                    STEP, PIMASTER,PAMASTER,KEEP(216),LRLUS,
     &                    KEEP(IXSZ))
                 IF ( LRLU .NE. LRLUS ) THEN
                   WRITE(*,*) 'PB compress... alloc_cb',
     &             'LRLU,LRLUS=',LRLU,LRLUS
                   GOTO 620
                 END IF
          COMP = COMP + 1
          IF (IWPOSCB-IWPOS+1 .LT. LREQ) GOTO 600
      ENDIF
      IXXP=IWPOSCB+XXP+1
      IF (IXXP.GT.LIW) THEN
        WRITE(*,*) "Internal error 3 in DMUMPS_22",IXXP
      ENDIF
      IF (IW(IXXP).GT.0) THEN
        WRITE(*,*) "Internal error 2 in DMUMPS_22",IW(IXXP),IXXP
      ENDIF
      IWPOSCB = IWPOSCB - LREQ
      IF (SET_HEADER) THEN
        IW(IXXP)= IWPOSCB + 1
        IW(IWPOSCB+1+XXI)=LREQ
        CALL MUMPS_730(LREQCB, IW(IWPOSCB+1+XXR))
        IW(IWPOSCB+1+XXS)=STATE_ARG
        IW(IWPOSCB+1+XXN)=NODE_ARG
        IW(IWPOSCB+1+XXP)=TOP_OF_STACK
      ENDIF
      IPTRLU = IPTRLU - LREQCB
      LRLU   = LRLU - LREQCB
      LRLUS  = LRLUS - LREQCB_EFF
      KEEP8(67) = min(LRLUS, KEEP8(67))
#if ! defined(OLD_LOAD_MECHANISM)
      CALL DMUMPS_471(SSARBR,PROCESS_BANDE,
     &              LA-LRLUS,0_8,LREQCB_EFF,KEEP,KEEP8,LRLU)
#else
#if defined (CHECK_COHERENCE)
      CALL DMUMPS_471(SSARBR,PROCESS_BANDE,
     &              LA-LRLUS,0_8,LREQCB_EFF,KEEP,KEEP8,LRLU)
#else
      CALL DMUMPS_471(SSARBR,.FALSE.,
     &              LA-LRLUS,0_8,LREQCB_EFF,KEEP,KEEP8,LRLU)
#endif
#endif
      RETURN
 600  IFLAG  = -8
      IERROR = LREQ
      RETURN
 620  IFLAG  = -9
      CALL MUMPS_731(LREQCB_EFF - LRLUS, IERROR)
      RETURN
      END SUBROUTINE DMUMPS_22
      SUBROUTINE DMUMPS_244(N, NSTEPS,
     & A, LA, IW, LIW, SYM_PERM, NA, LNA,
     & NE_STEPS, NFSIZ, FILS,
     & STEP, FRERE, DAD, CAND, 
     & ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     & PTRAR, LDPTRAR,
     & PTRIST, PTLUST_S, PTRFAC, IW1, IW2, ITLOC, RHS_MUMPS,
     & POOL, LPOOL,  
     & CNTL1, ICNTL, INFO, RINFO, KEEP,KEEP8,PROCNODE_STEPS,
     & SLAVEF,
     & COMM_NODES, MYID, MYID_NODES,
     & BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,
     & root, NELT, FRTPTR, FRTELT, COMM_LOAD,
     & ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     & MEM_DISTRIB,
     & DKEEP,PIVNUL_LIST,LPN_LIST)
      USE DMUMPS_LOAD 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER(8) :: LA
      INTEGER N,NSTEPS,LIW,LPOOL,SLAVEF,COMM_NODES
      INTEGER MYID, MYID_NODES,LNA
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION RINFO(40)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER BUFR( LBUFR )
      INTEGER NELT, LDPTRAR
      INTEGER FRTPTR(*), FRTELT(*)
      DOUBLE PRECISION CNTL1
      INTEGER   ICNTL(40)
      INTEGER   INFO(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER   IW(LIW), SYM_PERM(N), NA(LNA),
     &          NE_STEPS(KEEP(28)), FILS(N),
     &          FRERE(KEEP(28)), NFSIZ(KEEP(28)), 
     &          DAD(KEEP(28))
      INTEGER   CAND(SLAVEF+1, max(1,KEEP(56)))
      INTEGER   STEP(N)
      INTEGER   PTRAR(LDPTRAR,2)
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER   PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER   IW1(3*KEEP(28)), ITLOC(N+KEEP(253)), POOL(LPOOL)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: IW2(2*KEEP(28))
      INTEGER   PROCNODE_STEPS(KEEP(28))
      INTEGER   COMM_LOAD, ASS_IRECV
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER   INTARR(max(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      DOUBLE PRECISION SEUIL, SEUIL_LDLT_NIV2
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
       INTEGER MUMPS_275
       EXTERNAL MUMPS_275
      DOUBLE PRECISION UULOC 
      INTEGER LP, MPRINT
      INTEGER NSTK,PTRAST, NBPROCFILS
      INTEGER PIMASTER, PAMASTER
      LOGICAL PROK
      DOUBLE PRECISION ZERO, ONE
      DATA ZERO /0.0D0/
      DATA ONE /1.0D0/
      INTRINSIC int,real,log
      INTEGER IERR
      INTEGER NTOTPV, NTOTPVTOT, NMAXNPIV
      INTEGER(8) :: POSFAC, LRLU, IPTRLU, LRLUS
      INTEGER IWPOS, LEAF, NBROOT, NROOT
      KEEP(41)=0
      KEEP(42)=0
      NSTEPS   = 0
      LP     = ICNTL(1)
      MPRINT = ICNTL(2)
      PROK   = (MPRINT.GT.0)
      UULOC = CNTL1
      IF (UULOC.GT.ONE)   UULOC=ONE
      IF (UULOC.LT.ZERO)  UULOC=ZERO
      IF (KEEP(50).NE.0.AND.UULOC.GT.0.5D0) THEN
        UULOC = 0.5D0
      ENDIF
      PIMASTER   = 1
      NSTK       = PIMASTER + KEEP(28)
      NBPROCFILS = NSTK + KEEP(28)
      PTRAST = 1
      PAMASTER = 1 + KEEP(28)
      IF (KEEP(4).LE.0) KEEP(4)=32
      IF (KEEP(5).LE.0) KEEP(5)=16
      IF (KEEP(5).GT.KEEP(4)) KEEP(5) = KEEP(4)
      IF (KEEP(6).LE.0) KEEP(6)=24
      IF (KEEP(3).LE.KEEP(4)) KEEP(3)=KEEP(4)*2
      IF (KEEP(6).GT.KEEP(3)) KEEP(6) = KEEP(3)
      POSFAC = 1_8
      IWPOS  = 1
      LRLU = LA
      LRLUS = LRLU
      KEEP8(67) = LRLUS
      IPTRLU = LRLU
      NTOTPV   = 0
      NMAXNPIV = 0
      IW1(NSTK:NSTK+KEEP(28)-1) = NE_STEPS(1:KEEP(28))
      CALL MUMPS_362(N, LEAF, NBROOT, NROOT,
     &                     MYID_NODES,
     &                     SLAVEF, NA, LNA,
     &                     KEEP,KEEP8, STEP,
     &                     PROCNODE_STEPS,
     &                     POOL, LPOOL)
      CALL DMUMPS_506(POOL, LPOOL, LEAF)
      CALL DMUMPS_555(POOL, LPOOL,KEEP,KEEP8)
      IF ( KEEP( 38 ) .NE. 0 ) THEN
        NBROOT = NBROOT + root%NPROW * root%NPCOL - 1
      END IF
      IF ( root%yes )  THEN 
         IF ( MUMPS_275( PROCNODE_STEPS(STEP(KEEP(38))), SLAVEF )
     &         .NE. MYID_NODES ) THEN
             NROOT = NROOT + 1
         END IF
      END IF
      CALL DMUMPS_251(N,IW,LIW,A,LA,IW1(NSTK),IW1(NBPROCFILS),
     &         INFO(1),NFSIZ,FILS,STEP,FRERE, DAD, CAND,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &         INFO(11), NTOTPV, NMAXNPIV, PTRIST,IW2(PTRAST),
     &         IW1(PIMASTER), IW2(PAMASTER), PTRAR(1,2), 
     &         PTRAR(1,1), 
     &         ITLOC, RHS_MUMPS, INFO(2), POOL, LPOOL, 
     &         RINFO, POSFAC,IWPOS,LRLU,IPTRLU, 
     &         LRLUS, LEAF, NROOT, NBROOT,
     &         UULOC,ICNTL,PTLUST_S,PTRFAC,NSTEPS,INFO,
     &         KEEP,KEEP8,
     &         PROCNODE_STEPS,SLAVEF,MYID,COMM_NODES,
     &         MYID_NODES, BUFR,LBUFR, LBUFR_BYTES,
     &         INTARR, DBLARR, root, SYM_PERM,
     &         NELT, FRTPTR, FRTELT, LDPTRAR, 
     &         COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     &         MEM_DISTRIB,NE_STEPS,
     &     DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
      POSFAC = POSFAC -1_8
      IWPOS = IWPOS -1
      IF (KEEP(201).LE.0) THEN
        KEEP8(31) = POSFAC
      ENDIF
      KEEP(32) = IWPOS
      CALL MUMPS_735(KEEP8(31), INFO(9))
      INFO(10) = KEEP(32)
      KEEP8(67) = LA - KEEP8(67)
      KEEP(89)  = NTOTPV
      KEEP(246) = NMAXNPIV
      INFO(23) = KEEP(89)
      CALL MPI_ALLREDUCE(NTOTPV, NTOTPVTOT, 1, MPI_INTEGER, MPI_SUM, 
     &                COMM_NODES, IERR)
      IF ( ( (INFO(1).EQ.-10 .OR. INFO(1).EQ.-40)
     &       .AND. (NTOTPVTOT.EQ.N) )
     &              .OR. ( NTOTPVTOT.GT.N ) ) THEN
       write(*,*) ' Error 1 in mc51d NTOTPVTOT=', NTOTPVTOT 
       CALL MUMPS_ABORT()
      ENDIF
      IF ( (KEEP(19).NE.0 ) .AND. (NTOTPVTOT.NE.N) .AND. 
     & (INFO(1).GE.0) )  THEN
       write(*,*) ' Error 2 in mc51d NTOTPVTOT=', NTOTPVTOT 
       CALL MUMPS_ABORT()
      ENDIF
      IF ( (INFO(1) .GE. 0 ) 
     &      .AND. (NTOTPVTOT.NE.N) ) THEN
         INFO(1) = -10
         INFO(2) = NTOTPVTOT
      ENDIF
      IF (PROK) THEN
        WRITE (MPRINT,99980) INFO(1), INFO(2),
     &       KEEP(28), KEEP8(31), INFO(10), INFO(11), INFO(12),
     &       INFO(13), INFO(14), INFO(25), RINFO(2), RINFO(3)
      ENDIF
      RETURN
99980 FORMAT (/' LEAVING FACTORIZATION PHASE WITH ...'/
     &      ' INFO (1)                                      =',I15/
     &      '  --- (2)                                      =',I15/
     &      '           NUMBER OF NODES IN THE TREE         =',I15/
     &      ' INFO (9)  REAL SPACE FOR FACTORS              =',I15/
     &      '  --- (10) INTEGER SPACE FOR FACTORS           =',I15/
     &      '  --- (11) MAXIMUM SIZE OF FRONTAL MATRICES    =',I15/
     &      '  --- (12) NUMBER OF OFF DIAGONAL PIVOTS       =',I15/
     &      '  --- (13) NUMBER OF DELAYED PIVOTS            =',I15/
     &      '  --- (14) NUMBER OF MEMORY COMPRESSES         =',I15/
     &      '  --- (25) NUMBER OF ENTRIES IN FACTORS        =',I15/
     &  ' RINFO(2)  OPERATIONS DURING NODE ASSEMBLY     =',1PD10.3/
     &  ' -----(3)  OPERATIONS DURING NODE ELIMINATION  =',1PD10.3)
99990 FORMAT(/ ' NUMBER OF MSGS RECVD FOR DYNAMIC LOAD  =',I15)
      END SUBROUTINE DMUMPS_244
      SUBROUTINE DMUMPS_269( MYID,KEEP,KEEP8,
     &           BUFR, LBUFR, LBUFR_BYTES,
     &           IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS,
     &           N, IW, LIW, A, LA,
     &           PTRIST, PTRAST, STEP, PIMASTER, PAMASTER,
     &           NSTK_S, COMP,
     &           FPERE, FLAG, IFLAG, IERROR, COMM,
     &           ITLOC, RHS_MUMPS )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER MYID
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER KEEP(500), BUFR( LBUFR )
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: PTRAST  (KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) )
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP, FPERE
      LOGICAL FLAG
      INTEGER NSTK_S( KEEP(28) ), ITLOC( N + KEEP(253) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER IFLAG, IERROR, COMM
      INTEGER POSITION, FINODE, FLCONT, LREQ
      INTEGER(8) :: LREQCB
      INTEGER(8) :: IPOS_NODE, ISHIFT_PACKET
      INTEGER SIZE_PACKET
      INTEGER NBROWS_ALREADY_SENT, NBROWS_PACKET
      INCLUDE 'mumps_headers.h'
      LOGICAL COMPRESSCB
      FLAG = .FALSE.
      POSITION = 0
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                FINODE, 1, MPI_INTEGER,
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                FPERE, 1, MPI_INTEGER, 
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                FLCONT, 1, MPI_INTEGER,
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                NBROWS_PACKET, 1, MPI_INTEGER,
     &                COMM, IERR)
      COMPRESSCB = (FLCONT.LT.0) 
      IF (COMPRESSCB) THEN
        FLCONT   = -FLCONT
        LREQCB  = (int(FLCONT,8) * int(FLCONT+1,8)) / 2_8
      ELSE
        LREQCB  = int(FLCONT,8) * int(FLCONT,8)
      ENDIF
      IF (NBROWS_ALREADY_SENT == 0) THEN
        LREQ    = 2 * FLCONT + 6 + KEEP(IXSZ)
        IF (IPTRLU.LT.0_8) WRITE(*,*) 'before alloc_cb:IPTRLU = ',IPTRLU
        CALL DMUMPS_22( .FALSE., 0_8, .FALSE.,.FALSE.,
     &  MYID,N, KEEP,KEEP8, IW, LIW, A, LA,
     &  LRLU, IPTRLU,IWPOS,IWPOSCB,
     &  PTRIST,PTRAST,STEP, PIMASTER, PAMASTER,
     &  LREQ, LREQCB, FINODE, S_NOTFREE, .TRUE.,
     &  COMP, LRLUS, IFLAG, IERROR
     &     )
        IF (IPTRLU.LT.0_8) WRITE(*,*) 'after alloc_cb:IPTRLU = ',IPTRLU
        IF ( IFLAG .LT. 0 ) RETURN
        PIMASTER(STEP( FINODE )) = IWPOSCB + 1
        PAMASTER(STEP( FINODE )) = IPTRLU  + 1_8
        IF (COMPRESSCB)  IW(IWPOSCB + 1 + XXS ) = S_CB1COMP
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        IW(IWPOSCB + 1+KEEP(IXSZ)), LREQ-KEEP(IXSZ),
     &        MPI_INTEGER, COMM, IERR)
      ENDIF
      IF (COMPRESSCB) THEN
        ISHIFT_PACKET = int(NBROWS_ALREADY_SENT,8) *
     &                  int(NBROWS_ALREADY_SENT+1,8) / 2_8
        SIZE_PACKET = (NBROWS_PACKET * (NBROWS_PACKET+1))/2 +
     &                 NBROWS_ALREADY_SENT * NBROWS_PACKET
      ELSE
        ISHIFT_PACKET = int(NBROWS_ALREADY_SENT,8) * int(FLCONT,8)
        SIZE_PACKET = NBROWS_PACKET * FLCONT
      ENDIF
      IF (NBROWS_PACKET.NE.0) THEN
        IF ( LREQCB .ne. 0_8 ) THEN
        IPOS_NODE = PAMASTER(STEP(FINODE))-1_8
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        A(IPOS_NODE + 1_8 + ISHIFT_PACKET),
     &        SIZE_PACKET, MPI_DOUBLE_PRECISION, COMM, IERR)
        END IF
      ENDIF
      IF (NBROWS_ALREADY_SENT+NBROWS_PACKET == FLCONT) THEN
        NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
        IF ( NSTK_S(STEP(FPERE)).EQ.0 ) THEN
          FLAG = . TRUE.
        END IF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_269
      SUBROUTINE DMUMPS_270( TOT_ROOT_SIZE,
     &    TOT_CONT_TO_RECV, root,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND)
      USE DMUMPS_LOAD
      USE DMUMPS_OOC        
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) ), ND( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, COMM_LOAD
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC(N+KEEP(253)), FILS(N), PTRARW(N), PTRAIW(N)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      INTEGER ::  allocok
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: TMP
      INTEGER NEW_LOCAL_M, NEW_LOCAL_N
      INTEGER OLD_LOCAL_M, OLD_LOCAL_N
      INTEGER I, J
      INTEGER LREQI, IROOT
      INTEGER(8) :: LREQA
      INTEGER POSHEAD, IPOS_SON,IERR
      LOGICAL MASTER_OF_ROOT
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INCLUDE 'mumps_headers.h'
      INTEGER numroc, MUMPS_275
      EXTERNAL numroc, MUMPS_275
      IROOT = KEEP( 38 )
      root%TOT_ROOT_SIZE = TOT_ROOT_SIZE
      MASTER_OF_ROOT = ( MYID .EQ. 
     &                   MUMPS_275( PROCNODE_STEPS(STEP(IROOT)),
     &                   SLAVEF ) )
      NEW_LOCAL_M  = numroc( TOT_ROOT_SIZE, root%MBLOCK,
     &               root%MYROW, 0, root%NPROW )
      NEW_LOCAL_M  = max( 1, NEW_LOCAL_M )
      NEW_LOCAL_N  = numroc( TOT_ROOT_SIZE, root%NBLOCK,
     &               root%MYCOL, 0, root%NPCOL )
      IF ( PTRIST(STEP( IROOT )).GT.0) THEN
        OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) + KEEP(IXSZ) )
        OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1  + KEEP(IXSZ))
      ELSE
        OLD_LOCAL_N = 0
        OLD_LOCAL_M = NEW_LOCAL_M
      ENDIF
      IF (KEEP(60) .NE. 0) THEN
        IF (root%yes) THEN
        IF ( NEW_LOCAL_M .NE. root%SCHUR_MLOC .OR.
     &       NEW_LOCAL_N .NE. root%SCHUR_NLOC ) THEN
          WRITE(*,*) "Internal error 1 in DMUMPS_270"
          CALL MUMPS_ABORT()
        ENDIF
        ENDIF
        PTLUST_S(STEP(IROOT)) = -4444
        PTRFAC(STEP(IROOT)) = -4445_8
        PTRIST(STEP(IROOT)) = 0
        IF ( MASTER_OF_ROOT ) THEN
          LREQI=6+2*TOT_ROOT_SIZE+KEEP(IXSZ)
          LREQA=0_8
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
           CALL DMUMPS_94( N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &           KEEP(IXSZ))
           COMP = COMP + 1
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB1 compress root2slave:LRLU,LRLUS=',
     &            LRLU, LRLUS
                  IFLAG = -9
                  CALL MUMPS_731(LREQA-LRLUS, IERROR)
                  GOTO 700
           END IF
          ENDIF
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
            IFLAG = -8
            IERROR = IWPOS + LREQI - 1 - IWPOSCB
            GOTO 700
          ENDIF
          PTLUST_S(STEP(IROOT))= IWPOS
          IWPOS = IWPOS + LREQI
          POSHEAD = PTLUST_S( STEP(IROOT))
          IW( POSHEAD + XXI )=LREQI
          CALL MUMPS_730( LREQA, IW(POSHEAD + XXR))
          IW( POSHEAD + XXS )=-9999
          IW( POSHEAD +KEEP(IXSZ)) = 0
          IW( POSHEAD + 1 +KEEP(IXSZ)) = -1
          IW( POSHEAD + 2 +KEEP(IXSZ)) = -1
          IW( POSHEAD + 4 +KEEP(IXSZ)) = STEP(IROOT)
          IW( POSHEAD + 5 +KEEP(IXSZ)) = 0
          IW( POSHEAD + 3 +KEEP(IXSZ)) = TOT_ROOT_SIZE
        ENDIF
        GOTO 100
      ENDIF
      IF ( MASTER_OF_ROOT ) THEN
        LREQI = 6 + 2 * TOT_ROOT_SIZE+KEEP(IXSZ)
      ELSE
        LREQI = 6+KEEP(IXSZ)
      END IF
      LREQA = int(NEW_LOCAL_M, 8) * int(NEW_LOCAL_N, 8)
      IF ( LRLU . LT. LREQA .OR.
     &     IWPOS + LREQI - 1. GT. IWPOSCB )THEN
           IF ( LRLUS .LT. LREQA ) THEN
             IFLAG  = -9
             CALL MUMPS_731(LREQA - LRLUS, IERROR)
             GOTO 700
           END IF
           CALL DMUMPS_94( N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &           KEEP(IXSZ))
           COMP = COMP + 1
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB2 compress root2slave:LRLU,LRLUS=',
     &            LRLU, LRLUS
                  IFLAG = -9
                  CALL MUMPS_731(LREQA - LRLUS, IERROR)
                  GOTO 700
           END IF
           IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
              IFLAG  = -8
              IERROR = IWPOS + LREQI - 1 - IWPOSCB
              GOTO 700
           END IF
      END IF
      PTLUST_S(STEP( IROOT )) = IWPOS
      IWPOS           = IWPOS + LREQI
      IF (LREQA.EQ.0_8) THEN
        PTRAST (STEP(IROOT)) = max(POSFAC-1_8,1_8)
        PTRFAC (STEP(IROOT)) = max(POSFAC-1_8,1_8)
      ELSE
        PTRAST (STEP(IROOT)) = POSFAC
        PTRFAC (STEP(IROOT)) = POSFAC
      ENDIF
      POSFAC           = POSFAC + LREQA
      LRLU   = LRLU  - LREQA
      LRLUS  = LRLUS - LREQA
      KEEP8(67) = min(KEEP8(67), LRLUS)
      CALL DMUMPS_471(.FALSE.,.FALSE.,
     &          LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLU)
      POSHEAD = PTLUST_S( STEP(IROOT))
      IW( POSHEAD + XXI )     = LREQI
      CALL MUMPS_730( LREQA, IW(POSHEAD + XXR))
      IW( POSHEAD + XXS ) = S_NOTFREE
      IW( POSHEAD + KEEP(IXSZ) ) = 0
      IW( POSHEAD + 1 + KEEP(IXSZ) ) = NEW_LOCAL_N
      IW( POSHEAD + 2 + KEEP(IXSZ) ) = NEW_LOCAL_M
      IW( POSHEAD + 4 + KEEP(IXSZ) ) = STEP(IROOT)
      IW( POSHEAD + 5 + KEEP(IXSZ) ) = 0
      IF ( MASTER_OF_ROOT ) THEN
        IW( POSHEAD + 3 + KEEP(IXSZ) ) = TOT_ROOT_SIZE
      ELSE
        IW( POSHEAD + 3 + KEEP(IXSZ) ) = 0
      ENDIF
      IF ( KEEP( 50 ) .eq. 0 .OR. KEEP(50) .eq. 2 ) THEN
      OPELIW = OPELIW + ( dble(2*TOT_ROOT_SIZE) *
     &         dble(TOT_ROOT_SIZE) * dble(TOT_ROOT_SIZE ) / dble(3)
     &       - 0.5d0 * dble( TOT_ROOT_SIZE ) * dble( TOT_ROOT_SIZE )
     &       - dble( TOT_ROOT_SIZE ) / dble( 6 ) )
     &       / dble( root%NPROW * root%NPCOL )
      ELSE
      OPELIW = OPELIW + ( dble(TOT_ROOT_SIZE) *
     &         dble( TOT_ROOT_SIZE) * 
     &         dble( TOT_ROOT_SIZE + 1 ) )
     &          / dble( 3 * root%NPROW * root%NPCOL )
      END IF
      IF ( PTRIST(STEP( IROOT )) .LE. 0 ) THEN
        PTRIST(STEP( IROOT ))            = 0
        PAMASTER(STEP( IROOT ))          = 0_8
        IF (LREQA.GT.0_8) A(PTRAST(STEP(IROOT)):
     &    PTRAST(STEP(IROOT))+LREQA-1_8) = ZERO
      ELSE
        OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) + KEEP(IXSZ) )
        OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1  + KEEP(IXSZ))
        IF ( TOT_ROOT_SIZE .eq. root%ROOT_SIZE ) THEN
          IF ( LREQA .NE. int(OLD_LOCAL_M,8) * int(OLD_LOCAL_N,8) )
     &    THEN
             write(*,*) 'error 1 in PROCESS_ROOT2SLAVE',
     &       OLD_LOCAL_M, OLD_LOCAL_N
             CALL MUMPS_ABORT()
          END IF
          CALL DMUMPS_756(LREQA,
     &                          A( PAMASTER(STEP(IROOT)) ),
     &                          A( PTRAST  (STEP(IROOT)) ) )
        ELSE
          CALL DMUMPS_96( A( PTRAST(STEP(IROOT))), 
     &        NEW_LOCAL_M,
     &        NEW_LOCAL_N, A( PAMASTER( STEP(IROOT)) ), OLD_LOCAL_M,
     &        OLD_LOCAL_N )
        END IF
        IF ( PTRIST( STEP( IROOT ) ) .GT. 0 ) THEN
           IPOS_SON= PTRIST( STEP(IROOT))
           CALL DMUMPS_152(.FALSE., MYID, N, IPOS_SON,
     &          PAMASTER(STEP(IROOT)),
     &          IW, LIW, LRLU, LRLUS, IPTRLU,
     &          IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &         )
           PTRIST(STEP( IROOT ))   = 0
           PAMASTER(STEP( IROOT )) = 0_8
        END IF
      END IF
       IF (NEW_LOCAL_M.GT.OLD_LOCAL_M) THEN
          TMP => root%RHS_ROOT
          NULLIFY(root%RHS_ROOT)
          ALLOCATE (root%RHS_ROOT(NEW_LOCAL_M, root%RHS_NLOC), 
     &                stat=allocok)
          IF ( allocok.GT.0) THEN
              IFLAG=-13
              IERROR = NEW_LOCAL_M*root%RHS_NLOC
              GOTO 700
          ENDIF
          DO J = 1, root%RHS_NLOC
            DO I = 1, OLD_LOCAL_M
              root%RHS_ROOT(I,J)=TMP(I,J)
            ENDDO
            DO I = OLD_LOCAL_M+1, NEW_LOCAL_M
              root%RHS_ROOT(I,J) = ZERO
            ENDDO
          ENDDO
          DEALLOCATE(TMP)
          NULLIFY(TMP) 
       ENDIF
 100  CONTINUE
      NBPROCFILS(STEP(IROOT))=NBPROCFILS(STEP(IROOT)) + TOT_CONT_TO_RECV
      IF ( NBPROCFILS(STEP(IROOT)) .eq. 0 ) THEN
         IF (KEEP(201).EQ.1) THEN 
            CALL DMUMPS_681(IERR)
         ELSE IF (KEEP(201).EQ.2) THEN 
            CALL DMUMPS_580(IERR)              
         ENDIF
        CALL DMUMPS_507( N, IPOOL, LPOOL, PROCNODE_STEPS,
     &       SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     &       STEP, IROOT + N )
        IF (KEEP(47) .GE. 3) THEN
           CALL DMUMPS_500(
     &          IPOOL, LPOOL, 
     &          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &          MYID, STEP, N, ND, FILS )
        ENDIF
      END IF
      RETURN
 700  CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_270
      SUBROUTINE DMUMPS_96
     &( NEW, M_NEW, N_NEW,OLD, M_OLD, N_OLD )
      INTEGER M_NEW, N_NEW, M_OLD, N_OLD
      DOUBLE PRECISION NEW( M_NEW, N_NEW ), OLD( M_OLD, N_OLD )
      INTEGER J
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      DO J = 1, N_OLD
        NEW( 1: M_OLD, J ) = OLD( 1: M_OLD, J )
        NEW( M_OLD + 1: M_NEW, J ) = ZERO
      END DO
      NEW( 1: M_NEW,N_OLD + 1: N_NEW ) = ZERO
      RETURN
      END SUBROUTINE DMUMPS_96
      INTEGER FUNCTION DMUMPS_505(KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DMUMPS_505 = KEEP(28) + 1 + 3
      RETURN
      END FUNCTION DMUMPS_505
      SUBROUTINE DMUMPS_506(IPOOL, LPOOL, LEAF)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER LPOOL, LEAF
      INTEGER IPOOL(LPOOL)
      IPOOL(LPOOL-2) = 0
      IPOOL(LPOOL-1) = 0
      IPOOL(LPOOL)   = LEAF-1
      RETURN
      END SUBROUTINE DMUMPS_506
      SUBROUTINE DMUMPS_507
     &           (N, POOL, LPOOL, PROCNODE, SLAVEF,
     &           K28, K76, K80, K47, STEP, INODE)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N, INODE, LPOOL, K28, SLAVEF, K76, K80, K47
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(K28)
      EXTERNAL MUMPS_170
      LOGICAL MUMPS_170, ATM_CURRENT_NODE
      INTEGER NBINSUBTREE, NBTOP, INODE_EFF,POS_TO_INSERT
      INTEGER IPOS1, IPOS2, ISWAP
      INTEGER NODE,J,I
      ATM_CURRENT_NODE = ( K76 == 2 .OR. K76 ==3 .OR.
     &     K76==4 .OR. K76==5)
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      IF (INODE > N ) THEN
        INODE_EFF = INODE - N
      ELSE IF (INODE < 0) THEN
        INODE_EFF = - INODE
      ELSE
        INODE_EFF = INODE
      ENDIF
      IF(((INODE.GT.0).AND.(INODE.LE.N)).AND.(.NOT.
     &     MUMPS_170(PROCNODE(STEP(INODE_EFF)),
     &               SLAVEF))
     &  ) THEN
         IF ((K80 == 1 .AND. K47 .GE. 1) .OR.
     &     (( K80 == 2 .OR. K80==3 ) .AND.
     &          ( K47 == 4 ))) THEN
            CALL DMUMPS_514(INODE,1)
         ENDIF
      ENDIF
      IF ( MUMPS_170(PROCNODE(STEP(INODE_EFF)),
     &                             SLAVEF) ) THEN
        POOL(NBINSUBTREE + 1 ) = INODE
        NBINSUBTREE = NBINSUBTREE + 1
      ELSE
         POS_TO_INSERT=NBTOP+1
         IF((K76.EQ.4).OR.(K76.EQ.5))THEN
#if defined(NOT_ATM_POOL_SPECIAL)
            J=NBTOP
#else
            IF((INODE.GT.N).OR.(INODE.LE.0))THEN
               DO J=NBTOP,1,-1
                  IF((POOL(LPOOL-2-J).GT.0)
     &                 .AND.(POOL(LPOOL-2-J).LE.N))THEN
                     GOTO 333
                  ENDIF
                  IF ( POOL(LPOOL-2-J) < 0 ) THEN
                     NODE=-POOL(LPOOL-2-J)
                  ELSE IF ( POOL(LPOOL-2-J) > N ) THEN
                     NODE = POOL(LPOOL-2-J) - N
                  ELSE
                     NODE = POOL(LPOOL-2-J)
                  ENDIF
                  IF(K76.EQ.4)THEN
                     IF(DEPTH_FIRST_LOAD(STEP(NODE)).GE.
     &                    DEPTH_FIRST_LOAD(STEP(INODE_EFF)))THEN
                        GOTO 333
                     ENDIF
                  ENDIF
                  IF(K76.EQ.5)THEN
                     IF(COST_TRAV(STEP(NODE)).LE.
     &                    COST_TRAV(STEP(INODE_EFF)))THEN
                        GOTO 333
                     ENDIF
                  ENDIF
                  POS_TO_INSERT=POS_TO_INSERT-1
               ENDDO
               IF(J.EQ.0) J=1
 333           CONTINUE
               DO I=NBTOP,POS_TO_INSERT,-1
                  POOL(LPOOL-2-I-1)=POOL(LPOOL-2-I)
               ENDDO
               POOL(LPOOL-2-POS_TO_INSERT)=INODE
               NBTOP = NBTOP + 1
               GOTO 20
            ENDIF
            DO J=NBTOP,1,-1
               IF((POOL(LPOOL-2-J).GT.0).AND.(POOL(LPOOL-2-J).LE.N))THEN
                  GOTO 888
               ENDIF
               POS_TO_INSERT=POS_TO_INSERT-1
            ENDDO
 888        CONTINUE
#endif
            DO I=J,1,-1
#if defined(NOT_ATM_POOL_SPECIAL)
               IF ( POOL(LPOOL-2-I) < 0 ) THEN
                  NODE=-POOL(LPOOL-2-I)
               ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                  NODE = POOL(LPOOL-2-I) - N
               ELSE
                  NODE = POOL(LPOOL-2-I)
               ENDIF
#else
               NODE=POOL(LPOOL-2-I)
#endif
               IF(K76.EQ.4)THEN
                  IF(DEPTH_FIRST_LOAD(STEP(NODE)).GE.
     &                 DEPTH_FIRST_LOAD(STEP(INODE_EFF)))THEN
                     GOTO 999
                  ENDIF
               ENDIF
               IF(K76.EQ.5)THEN
                  IF(COST_TRAV(STEP(NODE)).LE.
     &                 COST_TRAV(STEP(INODE_EFF)))THEN
                     GOTO 999
                  ENDIF
               ENDIF
               POS_TO_INSERT=POS_TO_INSERT-1
            ENDDO
            IF(I.EQ.0) I=1
 999        CONTINUE
            DO J=NBTOP,POS_TO_INSERT,-1
               POOL(LPOOL-2-J-1)=POOL(LPOOL-2-J)
            ENDDO
            POOL(LPOOL-2-POS_TO_INSERT)=INODE
            NBTOP = NBTOP + 1
            GOTO 20
         ENDIF
         POOL( LPOOL - 2 - ( NBTOP + 1 ) ) = INODE
         NBTOP = NBTOP + 1
        IPOS1 = LPOOL - 2 - NBTOP
        IPOS2 = LPOOL - 2 - NBTOP + 1
 10     CONTINUE
        IF ( IPOS2 == LPOOL - 2 ) GOTO 20
        IF ( POOL(IPOS1) < 0 ) GOTO 20
        IF ( POOL(IPOS2) < 0 ) GOTO 30
        IF ( ATM_CURRENT_NODE ) THEN
          IF ( POOL(IPOS1) > N ) GOTO 20
          IF ( POOL(IPOS2) > N ) GOTO 30
        END IF
        GOTO 20
 30     CONTINUE
        ISWAP = POOL(IPOS1)
        POOL(IPOS1) = POOL(IPOS2)
        POOL(IPOS2) = ISWAP
        IPOS1 = IPOS1 + 1
        IPOS2 = IPOS2 + 1
        GOTO 10
 20     CONTINUE
      ENDIF
      POOL(LPOOL) = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      RETURN
      END SUBROUTINE DMUMPS_507
      LOGICAL FUNCTION DMUMPS_508(POOL, LPOOL)
      IMPLICIT NONE
      INTEGER LPOOL
      INTEGER POOL(LPOOL)
      INTEGER NBINSUBTREE, NBTOP
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      DMUMPS_508 = (NBINSUBTREE + NBTOP == 0)
      RETURN
      END FUNCTION DMUMPS_508
      SUBROUTINE DMUMPS_509( N, POOL, LPOOL, PROCNODE, SLAVEF,
     &           STEP, INODE, KEEP,KEEP8, MYID, ND,
     &           FORCE_EXTRACT_TOP_SBTR )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE, LPOOL, SLAVEF, N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(KEEP(28)),
     &        ND(KEEP(28))
      EXTERNAL MUMPS_167, MUMPS_283, DMUMPS_508
      LOGICAL MUMPS_167, MUMPS_283, DMUMPS_508
      EXTERNAL MUMPS_275
      INTEGER MUMPS_275
      INTEGER NBINSUBTREE, NBTOP, INSUBTREE, INODE_EFF, MYID
      LOGICAL LEFT, ATOMIC_SUBTREE,UPPER,FLAG_MEM,SBTR_FLAG,PROC_FLAG
      LOGICAL FORCE_EXTRACT_TOP_SBTR
      INTEGER NODE_TO_EXTRACT,I,J,MIN_PROC
#if defined(POOL_EXTRACT_MNG)
      INTEGER POS_TO_EXTRACT
#endif
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF ( KEEP(76) > 6 .OR. KEEP(76) < 0 ) THEN
         WRITE(*,*) "Error 2 in DMUMPS_509: unknown strategy"
         CALL MUMPS_ABORT()
      ENDIF
      ATOMIC_SUBTREE =  ( KEEP(76) == 1 .OR. KEEP(76) == 3)
      IF ( DMUMPS_508(POOL, LPOOL) ) THEN
         WRITE(*,*) "Error 1 in DMUMPS_509"
         CALL MUMPS_ABORT()
      ENDIF
      IF ( .NOT. ATOMIC_SUBTREE ) THEN
         LEFT = (NBTOP == 0)
         IF(.NOT.LEFT)THEN
            IF((KEEP(76).EQ.4).OR.(KEEP(76).EQ.5))THEN
               IF(NBINSUBTREE.EQ.0)THEN
                  LEFT=.FALSE.
               ELSE
                  IF ( POOL(NBINSUBTREE) < 0 ) THEN
                     I = -POOL(NBINSUBTREE)
                  ELSE IF ( POOL(NBINSUBTREE) > N ) THEN
                     I = POOL(NBINSUBTREE) - N
                  ELSE
                     I = POOL(NBINSUBTREE)
                  ENDIF
                  IF ( POOL(LPOOL-2-NBTOP) < 0 ) THEN
                     J = -POOL(LPOOL-2-NBTOP)
                  ELSE IF ( POOL(LPOOL-2-NBTOP) > N ) THEN
                     J = POOL(LPOOL-2-NBTOP) - N
                  ELSE
                     J = POOL(LPOOL-2-NBTOP)
                  ENDIF
                  IF(KEEP(76).EQ.4)THEN
                     IF(DEPTH_FIRST_LOAD(STEP(J)).GE.
     &                    DEPTH_FIRST_LOAD(STEP(I)))THEN
                        LEFT=.TRUE.
                     ELSE
                        LEFT=.FALSE.
                     ENDIF
                  ENDIF
                  IF(KEEP(76).EQ.5)THEN
                     IF(COST_TRAV(STEP(J)).LE.
     &                    COST_TRAV(STEP(I)))THEN
                        LEFT=.TRUE.
                     ELSE
                        LEFT=.FALSE.
                     ENDIF
                  ENDIF
               ENDIF           
            ENDIF
         ENDIF
      ELSE
         IF ( INSUBTREE == 1 ) THEN
            IF (NBINSUBTREE == 0) THEN
               WRITE(*,*) "Error 3 in DMUMPS_509"
               CALL MUMPS_ABORT()
            ENDIF
            LEFT = .TRUE.
         ELSE
            LEFT = ( NBTOP == 0)
         ENDIF
      ENDIF
 222  CONTINUE
      IF ( LEFT ) THEN
         INODE = POOL( NBINSUBTREE )
         IF(KEEP(81).EQ.2)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
            IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
               CALL DMUMPS_561(INODE,POOL,LPOOL,N,
     &              STEP,KEEP,KEEP8,PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     &              PROC_FLAG,MIN_PROC)
               IF(.NOT.SBTR_FLAG)THEN
                  WRITE(*,*)MYID,': ca a change pour moi'
                  LEFT=.FALSE.
                  GOTO 222
               ENDIF
#if ! defined(NOT_ATM_POOL_SPECIAL)
            ENDIF
#endif
         ELSEIF(KEEP(81).EQ.3)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
            IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
               NODE_TO_EXTRACT=INODE
               FLAG_MEM=.FALSE.
               CALL DMUMPS_820(FLAG_MEM)
               IF(FLAG_MEM)THEN
                  CALL DMUMPS_561(INODE,POOL,LPOOL,N,
     &                 STEP,KEEP,KEEP8,
     &                 PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     &                 PROC_FLAG,MIN_PROC)
                  IF(.NOT.SBTR_FLAG)THEN
                     LEFT=.FALSE.
                     WRITE(*,*)MYID,': ca a change pour moi (2)'
                     GOTO 222
                  ENDIF
               ENDIF
#if ! defined(NOT_ATM_POOL_SPECIAL)
            ENDIF
#endif
         ENDIF
         NBINSUBTREE = NBINSUBTREE - 1
         IF ( INODE < 0 ) THEN
            INODE_EFF = -INODE
         ELSE IF ( INODE > N ) THEN
            INODE_EFF = INODE - N
         ELSE
            INODE_EFF = INODE
         ENDIF
         IF ( MUMPS_167( PROCNODE(STEP(INODE_EFF)), SLAVEF) ) THEN
            IF((KEEP(47).GE.2.AND.KEEP(81).EQ.1).AND.
     &           (INSUBTREE.EQ.0))THEN
               CALL DMUMPS_513(.TRUE.)
            ENDIF
            INSUBTREE = 1
         ELSE IF ( MUMPS_283( PROCNODE(STEP(INODE_EFF)), 
     &           SLAVEF)) THEN
            IF((KEEP(47).GE.2.AND.KEEP(81).EQ.1).AND.
     &           (INSUBTREE.EQ.1))THEN
               CALL DMUMPS_513(.FALSE.)
            ENDIF
            INSUBTREE = 0
         END IF
      ELSE
         IF (NBTOP < 1 ) THEN
            WRITE(*,*) "Error 5 in DMUMPS_509", NBTOP
            CALL MUMPS_ABORT()
         ENDIF
         INODE = POOL( LPOOL - 2 - NBTOP )
         IF(KEEP(81).EQ.1)THEN
            CALL DMUMPS_520
     &           (INODE,UPPER,SLAVEF,KEEP,KEEP8,
     &            STEP,POOL,LPOOL,PROCNODE,N)
            IF(UPPER)THEN
               GOTO 666
            ELSE
               NBINSUBTREE=NBINSUBTREE-1
               IF ( MUMPS_167( PROCNODE(STEP(INODE)), 
     &              SLAVEF) ) THEN
                  INSUBTREE = 1
               ELSE IF ( MUMPS_283( PROCNODE(STEP(INODE)), 
     &                 SLAVEF)) THEN
                  INSUBTREE = 0
               ENDIF
               GOTO 777
            ENDIF
         ENDIF
         IF(KEEP(81).EQ.2)THEN
            CALL DMUMPS_561(INODE,POOL,LPOOL,N,STEP,
     &           KEEP,KEEP8,
     &           PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
            IF(SBTR_FLAG)THEN
               LEFT=.TRUE. 
               WRITE(*,*)MYID,': ca a change pour moi (3)'              
               GOTO 222
            ENDIF
         ELSE
#if defined(POOL_EXTRACT_MNG)
            IF(KEEP(76).EQ.4)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
               IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif                 
                  POS_TO_EXTRACT=-1
                  NODE_TO_EXTRACT=-1
                  DO I=NBTOP,1,-1
                     IF(NODE_TO_EXTRACT.LT.0)THEN
                        POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                        INODE_EFF = POOL(LPOOL-2-I)
                        IF ( POOL(LPOOL-2-I) < 0 ) THEN
                           NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                        ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                        ELSE
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                        ENDIF
#else
                        NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                     ELSE
                        IF(DEPTH_FIRST_LOAD(STEP(POOL(LPOOL-2-I))).LT.
     &                       DEPTH_FIRST_LOAD(STEP(NODE_TO_EXTRACT)))
     &                       THEN
                           POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                           INODE_EFF = POOL(LPOOL-2-I)
                           IF ( POOL(LPOOL-2-I) < 0 ) THEN
                              NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                           ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                           ELSE
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                           ENDIF
#else
                           NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                        ENDIF
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)
                  INODE = NODE_TO_EXTRACT
#else
                  INODE = INODE_EFF
#endif
                  DO I=POS_TO_EXTRACT,NBTOP
                     IF(I.NE.NBTOP)THEN
                        POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)
               ENDIF
#endif
            ENDIF
            IF(KEEP(76).EQ.5)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
               IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
                  POS_TO_EXTRACT=-1
                  NODE_TO_EXTRACT=-1
                  DO I=NBTOP,1,-1
                     IF(NODE_TO_EXTRACT.LT.0)THEN
                        POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                        INODE_EFF = POOL(LPOOL-2-I)
                        IF ( POOL(LPOOL-2-I) < 0 ) THEN
                           NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                        ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                        ELSE
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                        ENDIF
#else
                        NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                     ELSE
                        IF(COST_TRAV(STEP(POOL(LPOOL-2-I))).GT.
     &                       COST_TRAV(STEP(NODE_TO_EXTRACT)))THEN
                           POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                           INODE_EFF = POOL(LPOOL-2-I)
                           IF ( POOL(LPOOL-2-I) < 0 ) THEN
                              NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                           ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                           ELSE
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                           ENDIF
#else
                           NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                        ENDIF
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)              
                  INODE = NODE_TO_EXTRACT
#else
                  INODE = INODE_EFF
#endif
                  DO I=POS_TO_EXTRACT,NBTOP
                     IF(I.NE.NBTOP)THEN
                        POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)
               ENDIF
#endif
            ENDIF
#endif
            IF(KEEP(81).EQ.3)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
               IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
                  NODE_TO_EXTRACT=INODE
                  FLAG_MEM=.FALSE.
                  CALL DMUMPS_820(FLAG_MEM)
                  IF(FLAG_MEM)THEN
                     CALL DMUMPS_561(INODE,POOL,LPOOL,N,
     &                    STEP,KEEP,KEEP8,
     &                    PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     &                    PROC_FLAG,MIN_PROC)
                     IF(SBTR_FLAG)THEN
                        LEFT=.TRUE.
                        WRITE(*,*)MYID,': ca a change pour moi (4)'
                        GOTO 222
                     ENDIF
                  ELSE
                     CALL DMUMPS_819(INODE)
                  ENDIF
#if ! defined(NOT_ATM_POOL_SPECIAL)
               ENDIF
#endif
            ENDIF
         ENDIF
 666     CONTINUE
         NBTOP = NBTOP - 1
         IF((INODE.GT.0).AND.(INODE.LE.N))THEN
            IF ((( KEEP(80) == 2 .OR. KEEP(80)==3 ) .AND.
     &           ( KEEP(47) == 4 ))) THEN
               CALL DMUMPS_514(INODE,2)
            ENDIF
         ENDIF
         IF ( INODE < 0 ) THEN
            INODE_EFF = -INODE
         ELSE IF ( INODE > N ) THEN
            INODE_EFF = INODE - N
         ELSE
            INODE_EFF = INODE
         ENDIF
      END IF
 777  CONTINUE
      POOL(LPOOL)     = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      POOL(LPOOL - 2) = INSUBTREE
      RETURN
      END SUBROUTINE DMUMPS_509
      SUBROUTINE DMUMPS_552(INODE,POOL,LPOOL,N,STEP,
     &     KEEP,KEEP8,
     &     PROCNODE,SLAVEF,MYID,SBTR,FLAG_SAME_PROC,MIN_PROC)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE,LPOOL,N,MYID,SLAVEF,PROC,MIN_PROC
      INTEGER POOL(LPOOL),KEEP(500),STEP(N),PROCNODE(KEEP(28))
      INTEGER(8) KEEP8(150)
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      LOGICAL SBTR,FLAG_SAME_PROC
      INTEGER POS_TO_EXTRACT,NODE_TO_EXTRACT,NBTOP,I,INSUBTREE,
     &     NBINSUBTREE
      DOUBLE PRECISION MIN_COST, TMP_COST
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      MIN_COST=huge(MIN_COST) 
      TMP_COST=huge(TMP_COST)
      FLAG_SAME_PROC=.FALSE.
      SBTR=.FALSE.
      MIN_PROC=-9999
#if ! defined(NOT_ATM_POOL_SPECIAL)
      IF((INODE.GT.0).AND.(INODE.LE.N))THEN
#endif
         POS_TO_EXTRACT=-1
         NODE_TO_EXTRACT=-1
         DO I=NBTOP,1,-1
            IF(NODE_TO_EXTRACT.LT.0)THEN
               POS_TO_EXTRACT=I
               NODE_TO_EXTRACT=POOL(LPOOL-2-I)
               CALL DMUMPS_818(NODE_TO_EXTRACT,
     &                                       TMP_COST,PROC)
               MIN_COST=TMP_COST
               MIN_PROC=PROC
            ELSE
               CALL DMUMPS_818(POOL(LPOOL-2-I),
     &                                       TMP_COST,PROC)
               IF((PROC.NE.MIN_PROC).OR.(TMP_COST.NE.MIN_COST))THEN
                  FLAG_SAME_PROC=.TRUE.
               ENDIF
               IF(TMP_COST.GT.MIN_COST)THEN
                  POS_TO_EXTRACT=I
                  NODE_TO_EXTRACT=POOL(LPOOL-2-I)
                  MIN_COST=TMP_COST
                  MIN_PROC=PROC
               ENDIF
            ENDIF
         ENDDO
         IF((KEEP(47).EQ.4).AND.(NBINSUBTREE.NE.0))THEN
            CALL DMUMPS_554(NBINSUBTREE,INSUBTREE,NBTOP,
     &           MIN_COST,SBTR)
            IF(SBTR)THEN
               WRITE(*,*)MYID,': selecting from subtree'
               RETURN
            ENDIF
         ENDIF
         IF((.NOT.SBTR).AND.(.NOT.FLAG_SAME_PROC))THEN
            WRITE(*,*)MYID,': I must search for a task
     &           to save My friend'
            RETURN
         ENDIF
         INODE = NODE_TO_EXTRACT
         DO I=POS_TO_EXTRACT,NBTOP
            IF(I.NE.NBTOP)THEN
               POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
            ENDIF
         ENDDO
         POOL(LPOOL-2-NBTOP)=INODE
         CALL DMUMPS_819(INODE)
#if ! defined(NOT_ATM_POOL_SPECIAL)
      ELSE
      ENDIF
#endif
      END SUBROUTINE DMUMPS_552
      SUBROUTINE DMUMPS_561(INODE,POOL,LPOOL,N,STEP,
     &     KEEP,KEEP8,
     &     PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE,LPOOL,N,SLAVEF,MYID,MIN_PROC
      INTEGER POOL(LPOOL),KEEP(500),PROCNODE(KEEP(28)),STEP(N)
      INTEGER(8) KEEP8(150)
      LOGICAL SBTR_FLAG,PROC_FLAG
      EXTERNAL MUMPS_167
      LOGICAL MUMPS_167
      INTEGER NODE_TO_EXTRACT,I,POS_TO_EXTRACT,NBTOP,NBINSUBTREE
      NBTOP= POOL(LPOOL - 1)
      NBINSUBTREE = POOL(LPOOL)
      IF(NBTOP.GT.0)THEN
         WRITE(*,*)MYID,': NBTOP=',NBTOP
      ENDIF
      SBTR_FLAG=.FALSE.
      PROC_FLAG=.FALSE.
      CALL DMUMPS_552(INODE,POOL,LPOOL,N,STEP,KEEP,KEEP8,
     &     PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
      IF(SBTR_FLAG)THEN
         RETURN
      ENDIF
      IF(MIN_PROC.EQ.-9999)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
         IF((INODE.GT.0).AND.(INODE.LT.N))THEN
#endif
            SBTR_FLAG=(NBINSUBTREE.NE.0)
#if ! defined(NOT_ATM_POOL_SPECIAL)
         ENDIF
#endif
         RETURN
      ENDIF
      IF(.NOT.PROC_FLAG)THEN
         NODE_TO_EXTRACT=INODE
         IF((INODE.GE.0).AND.(INODE.LE.N))THEN
            CALL DMUMPS_553(MIN_PROC,POOL,
     &           LPOOL,INODE)
            IF(MUMPS_167(PROCNODE(STEP(INODE)),
     &           SLAVEF))THEN
               WRITE(*,*)MYID,': Extracting from a subtree
     &              for helping',MIN_PROC
               SBTR_FLAG=.TRUE.
               RETURN
            ELSE
               IF(NODE_TO_EXTRACT.NE.INODE)THEN
                  WRITE(*,*)MYID,': Extracting from top
     &                 inode=',INODE,'for helping',MIN_PROC
               ENDIF
               CALL DMUMPS_819(INODE)
            ENDIF
         ENDIF
         DO I=1,NBTOP
            IF (POOL(LPOOL-2-I).EQ.INODE)THEN
               GOTO 452
            ENDIF
         ENDDO
 452     CONTINUE
         POS_TO_EXTRACT=I
         DO I=POS_TO_EXTRACT,NBTOP-1
            POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
         ENDDO
         POOL(LPOOL-2-NBTOP)=INODE
      ENDIF
      END SUBROUTINE DMUMPS_561
      SUBROUTINE DMUMPS_574
     &           ( IPOOL, LPOOL, III, LEAF, 
     &             INODE, STRATEGIE )
            IMPLICIT NONE
      INTEGER, INTENT(IN) :: STRATEGIE, LPOOL
      INTEGER IPOOL (LPOOL)
      INTEGER III,LEAF
      INTEGER, INTENT(OUT) :: INODE
         LEAF  = LEAF - 1
         INODE = IPOOL( LEAF )
      RETURN
      END SUBROUTINE DMUMPS_574
      SUBROUTINE DMUMPS_128(N, NELT, ELTPTR, ELTVAR, LIW,
     &            IKEEP, PTRAR,
     &            IORD, NFSIZ, FILS, FRERE, 
     &            LISTVAR_SCHUR, SIZE_SCHUR,
     &            ICNTL, INFO, KEEP,KEEP8,
     &            ELTNOD, NSLAVES, 
     &            XNODEL, NODEL)
      IMPLICIT NONE
      INTEGER N,NELT,LIW,IORD, SIZE_SCHUR, NSLAVES
      INTEGER PTRAR(N,3), NFSIZ(N), FILS(N), FRERE(N)
      INTEGER ELTPTR(NELT+1) 
      INTEGER XNODEL(N+1), NODEL(ELTPTR(NELT+1)-1)
      INTEGER ELTVAR(ELTPTR(NELT+1)-1)
      INTEGER IKEEP(N,3)
      INTEGER LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER INFO(40), ICNTL(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER ELTNOD(NELT)
      INTEGER K,I,L1,L2,IWFR,NCMPA,LLIW,IFSON,IN
      INTEGER NEMIN, MPRINT, LP, MP, LDIAG
      INTEGER NZ, allocok, ITEMP
      LOGICAL PROK, NOSUPERVAR
      INTEGER(8) :: K79REF
      PARAMETER(K79REF=12000000_8)
      LOGICAL SPLITROOT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW2
          INTEGER OPT_METIS_SIZE, NUMFLAG
          PARAMETER(OPT_METIS_SIZE = 8, NUMFLAG = 1)
          INTEGER OPTIONS_METIS(OPT_METIS_SIZE)
      INTEGER IDUM 
      EXTERNAL MUMPS_197, DMUMPS_130, DMUMPS_131,
     &         DMUMPS_129, DMUMPS_132, 
     &         DMUMPS_133, DMUMPS_134,
     &         DMUMPS_199,
     &         DMUMPS_557, DMUMPS_201
#if defined(OLDDFS)
      EXTERNAL DMUMPS_200
#endif
        ALLOCATE( IW ( LIW ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = LIW
          RETURN
        ENDIF
      MPRINT= ICNTL(3)
      PROK  = (MPRINT.GT.0)
      LP    = ICNTL(1)
      MP    = ICNTL(3)
      LDIAG = ICNTL(4)
      IF (KEEP(60).NE.0) THEN
       NOSUPERVAR=.TRUE.
       IF (IORD.GT.1) IORD = 0
      ELSE
       NOSUPERVAR=.FALSE.
      ENDIF
      IF (IORD == 7) THEN
         IF ( N < 10000 ) THEN
           IORD = 0
         ELSE
#if defined(metis) || defined(parmetis)
           IORD = 5
#else
           IORD = 0
#endif
         ENDIF
      END IF
#if ! defined(metis) && ! defined(parmetis)
      IF (IORD == 5) IORD = 0
#endif
      IF (KEEP(1).LT.1) KEEP(1) = 1
      NEMIN = KEEP(1)
      IF (LDIAG.LE.2 .OR. MP.LE.0) GO TO 10
      WRITE (MP,99999) N, NELT, LIW, INFO(1)
      K = min0(10,NELT+1)
      IF (LDIAG.EQ.4) K = NELT+1
      IF (K.GT.0) WRITE (MP,99998) (ELTPTR(I),I=1,K)
      K = min0(10,ELTPTR(NELT+1)-1)
      IF (LDIAG.EQ.4) K = ELTPTR(NELT+1)-1
      IF (K.GT.0) WRITE (MP,99995) (ELTVAR(I),I=1,K)
      K = min0(10,N)
      IF (LDIAG.EQ.4) K = N
      IF (IORD.EQ.1 .AND. K.GT.0) THEN
        WRITE (MP,99997) (IKEEP(I,1),I=1,K)
      ENDIF
   10 L1 = 1
      L2 = L1 + N
      IF (LIW .LT. 3*N) THEN
          INFO(1)= -2002
          INFO(2) = LIW
      ENDIF
#if defined(metis) || defined(parmetis)
      IF ( IORD == 5 ) THEN
        IF (LIW .LT. N+N+1) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        ENDIF
      ELSE
#endif
      IF (NOSUPERVAR) THEN
        IF ( LIW .LT. 2*N ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        END IF
      ELSE
        IF ( LIW .LT.  4*N+4 ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        END IF
      ENDIF
#if defined(metis) || defined(parmetis)
      ENDIF
#endif
      IDUM=0
      CALL DMUMPS_258(NELT, N, ELTPTR(NELT+1)-1, ELTPTR, ELTVAR,
     &           XNODEL, NODEL, IW(L1), IDUM, ICNTL)
      IF (IORD.NE.1 .AND. IORD .NE. 5) THEN
        IORD = 0
        IF (NOSUPERVAR) THEN
          CALL DMUMPS_129(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              PTRAR(1,2), IW(L1))
        ELSE
         CALL DMUMPS_130(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              PTRAR(1,2), 4*N+4, IW(L1))
        ENDIF
        LLIW = max(NZ,N)
        ALLOCATE( IW2(LLIW), stat = allocok )
        IF (allocok.GT.0) THEN
          INFO(1) = -7
          INFO(2) = LLIW
          RETURN
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL DMUMPS_132(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IW2, LLIW, PTRAR, PTRAR(1,2),
     &              IW(L1), IWFR)
        ELSE
         CALL DMUMPS_131(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IW2, LLIW, PTRAR, PTRAR(1,2),
     &              IW(L1), IWFR)
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL MUMPS_162(N, LLIW, PTRAR, IWFR, PTRAR(1,2), IW2,
     &   IW(L1), IKEEP,
     &   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), PTRAR(1,3),
     &   LISTVAR_SCHUR, SIZE_SCHUR)
         IF (KEEP(60) == 1) THEN
           KEEP(20) = LISTVAR_SCHUR(1)
         ELSEIF (KEEP(60) == 2 .OR. KEEP(60) == 3 ) THEN
           KEEP(38) = LISTVAR_SCHUR(1)
         ELSE
           WRITE(*,*) "Internal error in DMUMPS_128",KEEP(60)
           CALL MUMPS_ABORT()
         ENDIF
        ELSE
         CALL MUMPS_23(N, LLIW, PTRAR, IWFR, PTRAR(1,2), IW2, 
     &   IW(L1), IKEEP, 
     &   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), PTRAR(1,3))
        ENDIF
      ELSE
#if defined(metis) || defined(parmetis)
        IF (IORD.EQ.5) THEN
         IF (PROK) THEN
          WRITE(MPRINT,'(A)') ' Ordering based on METIS '
         ENDIF
         CALL DMUMPS_129(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              PTRAR(1,2), IW(L1))
         LLIW = max(NZ,N)
         ALLOCATE( IW2(LLIW), stat = allocok )
         IF (allocok.GT.0) THEN
           INFO(1) = -7
           INFO(2) = LLIW
           RETURN
         ENDIF
         CALL DMUMPS_538(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IW2, LLIW, IW(L2), PTRAR(1,2),
     &              IW(L1), IWFR)
          OPTIONS_METIS(1) = 0
          CALL METIS_NODEND(N, IW(L2), IW2(1), NUMFLAG, OPTIONS_METIS,
     &       IKEEP(1,2), IKEEP(1,1) )
           DEALLOCATE(IW2)
        ELSE IF (IORD.NE.1) THEN
          WRITE(*,*) IORD
          WRITE(*,*)  'bad option for ordering'
          CALL MUMPS_ABORT()
        ENDIF
#endif
       DO K=1,N
         IW(L1+K) = 0
       ENDDO
       DO K=1,N
         IF ((IKEEP(K,1).LE.0).OR.(IKEEP(K,1).GT.N)) 
     &    GO TO 40
         IF (IW(L1+IKEEP(K,1)).EQ.1) THEN
          GOTO 40
         ELSE
          IW(L1+IKEEP(K,1)) = 1
         ENDIF
       ENDDO
       CALL DMUMPS_133(N, NZ, NELT, ELTPTR(NELT+1)-1,
     &             ELTPTR, ELTVAR, XNODEL, NODEL, 
     &             IKEEP, PTRAR(1,2), IW(L1))
       LLIW = NZ+N
       ALLOCATE( IW2(LLIW), stat = allocok )
       IF (allocok.GT.0) THEN
         INFO(1) = -7
         INFO(2) = LLIW
         RETURN
       ENDIF
       CALL DMUMPS_134(N, NZ, NELT, ELTPTR(NELT+1)-1,
     &             ELTPTR, ELTVAR, XNODEL, NODEL, 
     &             IKEEP, IW2, LLIW, PTRAR, PTRAR(1,2),
     &             IW(L1), IWFR)
       IF (KEEP(60) == 0) THEN
         ITEMP = 0 
       ELSE
         ITEMP = SIZE_SCHUR
         IF (KEEP(60) == 1) THEN
           KEEP(20) = LISTVAR_SCHUR(1)
         ELSEIF (KEEP(60) == 2 .OR. KEEP(60) == 3 ) THEN
           KEEP(38) = LISTVAR_SCHUR(1)
         ELSE
           WRITE(*,*) "Internal error in DMUMPS_128",KEEP(60)
           CALL MUMPS_ABORT()
         ENDIF
       ENDIF
       CALL DMUMPS_199(N, PTRAR, IW2, LLIW, IWFR, IKEEP, 
     &    IKEEP(1,2), IW(L1),
     &    IW(L2), NCMPA, ITEMP)
      ENDIF
#if defined(OLDDFS)
      CALL DMUMPS_200(N, PTRAR, IW(L1), IKEEP, IKEEP(1,2),
     &     IKEEP(1,3),
     &     NFSIZ, INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN, KEEP(60))
#else
      CALL DMUMPS_557(N, PTRAR, IW(L1), IKEEP, IKEEP(1,2),
     &     IKEEP(1,3),
     &     NFSIZ, PTRAR(1,2), 
     &     INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN, 
     &     IW(L2),  KEEP(60), KEEP(20), KEEP(38),
     &     IW2,KEEP(104),IW(L2+N),KEEP(50),
     &     ICNTL(13), KEEP(37), NSLAVES, KEEP(250).EQ.1)
#endif
      DEALLOCATE(IW2)
      IF (KEEP(60).NE.0) THEN
         IF (KEEP(60)==1) THEN
           IN = KEEP(20)
         ELSE
           IN = KEEP(38)
         ENDIF
         DO WHILE (IN.GT.0)
          IN = FILS (IN)
         END DO
         IFSON = -IN
         IF (KEEP(60)==1) THEN
           IN = KEEP(20)
         ELSE
           IN = KEEP(38)
         ENDIF
         DO I=2,SIZE_SCHUR
          FILS(IN) = LISTVAR_SCHUR (I)
          IN       = FILS(IN)
          FRERE (IN) = N+1
         ENDDO
         FILS(IN) = -IFSON
      ENDIF
      CALL DMUMPS_201(IKEEP(1,2),
     &  PTRAR(1,3), INFO(6),
     &  INFO(5), KEEP(2),KEEP(50),
     &  KEEP(101), KEEP(108),KEEP(5),
     &  KEEP(6), KEEP(226), KEEP(253))
      IF ( KEEP(53) .NE. 0 ) THEN
        CALL MUMPS_209( N, FRERE, FILS, NFSIZ, KEEP(20) )
      END IF
      IF ( KEEP(48) == 4 .OR.
     &   ( (KEEP(24).NE.0).AND.(KEEP8(21).GT.0_8) ) ) THEN
          CALL DMUMPS_510(KEEP8(21), KEEP(2),
     &    KEEP(48), KEEP(50), NSLAVES)
      END IF
      IF (KEEP(210).LT.0.OR.KEEP(210).GT.2) KEEP(210)=0
      IF (KEEP(210).EQ.0.AND.KEEP(201).GT.0) KEEP(210)=1 
      IF (KEEP(210).EQ.0.AND.KEEP(201).EQ.0) KEEP(210)=2 
      IF (KEEP(210).EQ.2) KEEP8(79)=huge(KEEP8(79))
      IF (KEEP(210).EQ.1.AND.KEEP8(79).LE.0_8) THEN
        IF ( huge(KEEP8(79)) / K79REF + 1_8 .GE. int(NSLAVES,8) ) THEN
        KEEP8(79)=huge(KEEP8(79))
        ELSE
        KEEP8(79)=K79REF * int(NSLAVES,8)
        ENDIF
      ENDIF
      IF (KEEP(79).EQ.0) THEN
       IF (KEEP(210).EQ.1) THEN
        SPLITROOT = .FALSE. 
        IF ( KEEP(62).GE.1) THEN
         CALL DMUMPS_97(N, FRERE, FILS, NFSIZ,INFO(6),
     &                       NSLAVES, KEEP,KEEP8, SPLITROOT,
     &                       MP, LDIAG, INFO(1), INFO(2))
         IF (INFO(1).LT.0) RETURN
        ENDIF
       ENDIF
      ENDIF
      SPLITROOT = ((ICNTL(13).GT.0) .AND. (NSLAVES.GE.ICNTL(13)))
      IF (SPLITROOT) THEN
         CALL DMUMPS_97(N, FRERE, FILS, NFSIZ,INFO(6),
     &                    NSLAVES, KEEP,KEEP8, SPLITROOT,
     &                    MP, LDIAG, INFO(1), INFO(2))
         IF (INFO(1).LT.0) RETURN
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
       K = min0(10,N)
       IF (LDIAG.EQ.4) K = N
       IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
       IF (K.GT.0) WRITE (MP,99991) (IKEEP(I,2),I=1,K)
       IF (K.GT.0) WRITE (MP,99990) (IKEEP(I,3),I=1,K)
       IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
   40 INFO(1) = -4
      INFO(2) = K
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996) INFO(1)
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
   90 CONTINUE
      DEALLOCATE(IW)
      RETURN
99999 FORMAT (/'Entering analysis phase with ...'/
     & '                N         NELT       LIW       INFO(1)'/,
     & 9X, I8, I11, I12, I14)
99998 FORMAT ('Element pointers:  ELTPTR()   '/(9X, 7I10))
99995 FORMAT ('Element variables: ELTVAR()   '/(9X, 7I10))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT (/'** Error return ** from Analysis   *  INFO(1)=', I3)
99991 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99990 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99989 FORMAT ('FILS (.)  =', 10I6/(12X, 10I6))
99988 FORMAT ('FRERE(.)  =', 10I6/(12X, 10I6))
99987 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99982 FORMAT ('Error in permutation array KEEP   INFO(2)=', I3)
      END SUBROUTINE DMUMPS_128
      SUBROUTINE DMUMPS_258( NELT, N, NELNOD, XELNOD, ELNOD,
     &                        XNODEL, NODEL, FLAG, IERROR, ICNTL ) 
      IMPLICIT NONE
      INTEGER NELT, N, NELNOD, IERROR, ICNTL(40)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     &        FLAG(N)
      INTEGER I, J, K, MP, NBERR
      MP = ICNTL(2)
      FLAG(1:N) = 0
      XNODEL(1:N) = 0
      IERROR = 0
      DO I = 1, NELT
        DO K = XELNOD(I), XELNOD(I+1)-1
          J = ELNOD(K)
          IF ( J.LT.1 .OR. J.GT.N ) THEN
            IERROR = IERROR + 1
          ELSE
            IF ( FLAG(J).NE.I ) THEN
              XNODEL(J) = XNODEL(J) + 1
              FLAG(J) = I
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF ( IERROR.GT.0 .AND. MP.GT.0 .AND. ICNTL(4).GE.2 ) THEN
        NBERR = 0
        WRITE(MP,99999)
        DO I = 1, NELT
          DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF ( J.LT.1 .OR. J.GT.N ) THEN
              NBERR = NBERR + 1
              IF (NBERR.LE.10) THEN
                WRITE(MP,'(A,I8,A,I8,A)')
     &          'Element ',I,' variable ',J,' ignored.' 
              ELSE
                GO TO 100
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
  100 CONTINUE
      K = 1
      DO I = 1, N
         K = K + XNODEL(I)
         XNODEL(I) = K
      ENDDO
      XNODEL(N+1) = XNODEL(N)
      FLAG(1:N) = 0
      DO I = 1, NELT
         DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF (FLAG(J).NE.I) THEN   
              XNODEL(J) = XNODEL(J) - 1
              NODEL(XNODEL(J)) = I
              FLAG(J) = I
            ENDIF
         ENDDO
      ENDDO
      RETURN
99999 FORMAT (/'*** Warning message from subroutine DMUMPS_258 ***')
      END SUBROUTINE DMUMPS_258
      SUBROUTINE DMUMPS_129(N, NZ, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  LEN, FLAG)
      IMPLICIT NONE
      INTEGER N, NELT, NELNOD, NZ
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER  XNODEL(N+1), NODEL(NELNOD),
     &        FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN 
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                LEN(I) = LEN(I) + 1
                LEN(J) = LEN(J) + 1
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0
      DO I = 1,N
        NZ = NZ + LEN(I)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_129
      SUBROUTINE DMUMPS_538(N, NZ, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER LEN(N)
      INTEGER IPE(N+1)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER  XNODEL(N+1), NODEL(NELNOD), 
     &          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IWFR = IWFR + LEN(I)
          IPE(I) = IWFR
      ENDDO 
      IPE(N+1)=IPE(N)
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                IPE(I) = IPE(I) - 1
                IW(IPE(I)) = J
                IPE(J) = IPE(J) - 1
                IW(IPE(J)) = I
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_538
      SUBROUTINE DMUMPS_132(N, NZ, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER LEN(N)
      INTEGER IPE(N)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER  XNODEL(N+1), NODEL(NELNOD), 
     &          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IWFR = IWFR + LEN(I)
        IF (LEN(I).GT.0) THEN
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                IPE(I) = IPE(I) - 1
                IW(IPE(I)) = J
                IPE(J) = IPE(J) - 1
                IW(IPE(J)) = I
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_132
      SUBROUTINE DMUMPS_133(N, NZ, NELT, NELNOD,
     & XELNOD, ELNOD, XNODEL, NODEL, 
     & PERM, LEN, FLAG)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  LEN(I) = LEN(I) + 1
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0  
      DO I = 1,N 
        NZ = NZ + LEN(I)
      ENDDO   
      RETURN  
      END SUBROUTINE DMUMPS_133
      SUBROUTINE DMUMPS_134(N, NZ, NELT, NELNOD,
     & XELNOD, ELNOD, XNODEL, NODEL, 
     & PERM, IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER  XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER IPE(N), LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), IW(LW), 
     &          FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 0
      DO I = 1,N
        IWFR = IWFR + LEN(I) + 1
        IPE(I) = IWFR 
      ENDDO
      IWFR = IWFR + 1
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  IW(IPE(I)) = J
                  IPE(I) = IPE(I) - 1
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DO I = 1,N
        J = IPE(I)
        IW(J) = LEN(I)
        IF (LEN(I).EQ.0) IPE(I) = 0
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_134
      SUBROUTINE DMUMPS_25( MYID, SLAVEF, N, 
     &           PROCNODE, STEP, PTRAIW, PTRARW, 
     &           NELT, FRTPTR, FRTELT, 
     &           KEEP,KEEP8, ICNTL, SYM )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, N, NELT, SYM
      INTEGER KEEP( 500 ), ICNTL( 40 )
      INTEGER(8) KEEP8(150)
      INTEGER PTRAIW( NELT+1 ), PTRARW( NELT+1 )
      INTEGER STEP( N )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PROCNODE( KEEP(28) )
      INTEGER MUMPS_330, MUMPS_275
      EXTERNAL MUMPS_330, MUMPS_275
      INTEGER ELT, I, K, IPTRI, IPTRR, NVAR
      INTEGER TYPE_PARALL, ITYPE, IRANK
      TYPE_PARALL = KEEP(46)
      PTRAIW( 1:NELT ) = 0
      DO I = 1, N
        IF (STEP(I).LT.0) CYCLE
        ITYPE = MUMPS_330( PROCNODE(abs(STEP(I))), SLAVEF )
        IRANK = MUMPS_275( PROCNODE(abs(STEP(I))), SLAVEF )
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK = IRANK + 1
        END IF
        IF ( (ITYPE .EQ. 2) .OR.
     &       (ITYPE .EQ. 1 .AND. IRANK .EQ. MYID) ) THEN
          DO K = FRTPTR(I),FRTPTR(I+1)-1
            ELT = FRTELT(K)
            PTRAIW( ELT ) = PTRARW(ELT+1) - PTRARW(ELT)
          ENDDO
        ELSE 
        END IF
      END DO
      IPTRI = 1
      DO ELT = 1,NELT
        NVAR = PTRAIW( ELT )
        PTRAIW( ELT ) = IPTRI
        IPTRI = IPTRI + NVAR
      ENDDO
      PTRAIW( NELT+1 ) = IPTRI
      KEEP( 14 ) = IPTRI - 1
      IF ( .TRUE. ) THEN  
        IF (SYM .EQ. 0) THEN
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + NVAR*NVAR
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ELSE
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + (NVAR*(NVAR+1))/2
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ENDIF
      ELSE
        IF (SYM .EQ. 0) THEN
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + NVAR*NVAR 
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ELSE
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + (NVAR*(NVAR+1))/2  
          ENDDO 
          PTRARW( NELT+1 ) = IPTRR
        ENDIF
      ENDIF 
      KEEP( 13 ) = IPTRR - 1
      RETURN
      END SUBROUTINE DMUMPS_25
      SUBROUTINE DMUMPS_120( N, NELT, ELTPROC, SLAVEF, PROCNODE )
      IMPLICIT NONE
      INTEGER N, NELT, SLAVEF
      INTEGER PROCNODE( N ), ELTPROC( NELT )
      INTEGER ELT, I, ITYPE, MUMPS_330, MUMPS_275
      EXTERNAL MUMPS_330, MUMPS_275
      DO ELT = 1, NELT
          I = ELTPROC(ELT)
          IF ( I .NE. 0) THEN
           ITYPE = MUMPS_330(PROCNODE(I),SLAVEF)
           IF (ITYPE.EQ.1) THEN
             ELTPROC(ELT) = MUMPS_275(PROCNODE(I),SLAVEF)
           ELSE IF (ITYPE.EQ.2) THEN
             ELTPROC(ELT) = -1
           ELSE
            ELTPROC(ELT) = -2
           ENDIF
          ELSE
           ELTPROC(ELT) = -3
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_120
      SUBROUTINE DMUMPS_153(N, NELT, NELNOD, FRERE, FILS, NA, NE,
     &           XNODEL, NODEL, FRTPTR, FRTELT, ELTNOD) 
      IMPLICIT NONE
      INTEGER  N, NELT, NELNOD
      INTEGER  FRERE(N), FILS(N), NA(N), NE(N)
      INTEGER  FRTPTR(N+1), FRTELT(NELT), ELTNOD(NELT)
      INTEGER  XNODEL(N+1), NODEL(NELNOD) 
      INTEGER TNSTK( N ), IPOOL( N )
      INTEGER I, K, IFATH
      INTEGER INODE, LEAF, NBLEAF, NBROOT, III, IN
      TNSTK = NE
      LEAF = 1
      IF (N.EQ.1) THEN
        NBROOT = 1
        NBLEAF = 1
        IPOOL(1) = 1
        LEAF = LEAF + 1
      ELSEIF (NA(N).LT.0) THEN
        NBLEAF = N
        NBROOT = N
        DO 20 I=1,NBLEAF-1
           INODE = NA(I)
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
 20     CONTINUE
        INODE = -NA(N)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSEIF (NA(N-1).LT.0) THEN
        NBLEAF = N-1
        NBROOT = NA(N)
        IF (NBLEAF-1.GT.0) THEN
         DO 30 I=1,NBLEAF-1
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 30      CONTINUE
        ENDIF
        INODE = -NA(N-1)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSE
        NBLEAF = NA(N-1)
        NBROOT = NA(N)
        DO 40 I = 1,NBLEAF
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 40     CONTINUE
      ENDIF
      ELTNOD(1:NELT) = 0
      III = 1
 90   CONTINUE
        IF (III.NE.LEAF) THEN
           INODE=IPOOL(III)
           III = III + 1
        ELSE 
           WRITE(6,*) ' ERROR 1 in file DMUMPS_153 '
           CALL MUMPS_ABORT()
        ENDIF
 95     CONTINUE 
        IN = INODE
 100    CONTINUE
        DO K = XNODEL(IN),XNODEL(IN+1)-1
          I = NODEL(K)
          IF (ELTNOD(I).EQ.0) ELTNOD(I) = INODE
        ENDDO
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
        IN = INODE
 110    IN = FRERE(IN)
        IF (IN.GT.0) GO TO 110
        IF (IN.EQ.0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
        ELSE
         IFATH = -IN
        ENDIF
        TNSTK(IFATH) = TNSTK(IFATH) - 1
        IF ( TNSTK(IFATH) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
        ELSE
            GOTO 90
        ENDIF
  115 CONTINUE
      FRTPTR(1:N) = 0
      DO I = 1,NELT
        IF (ELTNOD(I) .NE. 0) THEN
         FRTPTR(ELTNOD(I)) = FRTPTR(ELTNOD(I)) + 1
        ENDIF
      ENDDO
      K = 1
      DO I = 1,N
        K = K + FRTPTR(I)
        FRTPTR(I) = K
      ENDDO
      FRTPTR(N+1) = FRTPTR(N)
      DO K = 1,NELT
        INODE = ELTNOD(K)
        IF (INODE .NE. 0) THEN
         FRTPTR(INODE) = FRTPTR(INODE) - 1
         FRTELT(FRTPTR(INODE)) = K
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_153
      SUBROUTINE DMUMPS_130(N, NZ, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  LEN, LW, IW)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     &        IW(LW)
      INTEGER I,J,K1,K2,K3,LP,NSUP,SUPVAR
      INTEGER INFO44(6)
      EXTERNAL DMUMPS_315
      LP = 6
      CALL DMUMPS_315(N,NELT,XELNOD(NELT+1)-1,ELNOD,XELNOD,
     &           NSUP,IW(3*N+3+1),3*N+3,IW,LP,INFO44)
      IF (INFO44(1) .LT. 0) THEN
        IF (LP.GE.0) WRITE(LP,*) 
     &     'Error return from DMUMPS_315. INFO(1) = ',INFO44(1) 
      ENDIF
      IW(1:NSUP) = 0
      LEN(1:N) = 0
      DO I = 1,N
        SUPVAR = IW(3*N+3+1+I)
        IF (SUPVAR .EQ. 0) CYCLE
        IF (IW(SUPVAR).NE.0) THEN
          LEN(I) = -IW(SUPVAR)
        ELSE
          IW(SUPVAR) = I
        ENDIF
      ENDDO
      IW(N+1:2*N) = 0
      NZ = 0
      DO SUPVAR = 1,NSUP
        I = IW(SUPVAR)
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J).GE.0) THEN
                IF ((I.NE.J) .AND. (IW(N+J).NE.I)) THEN
                  IW(N+J) = I
                  LEN(I) = LEN(I) + 1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        NZ = NZ + LEN(I)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_130
      SUBROUTINE DMUMPS_131(N, NZ, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER IPE(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     &          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IF (LEN(I).GT.0) THEN
          IWFR = IWFR + LEN(I)
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        IF (LEN(I).LE.0) CYCLE
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J) .GT. 0) THEN
                IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                  IPE(I) = IPE(I) - 1
                  IW(IPE(I)) = J
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_131
      SUBROUTINE DMUMPS_315(N,NELT,NZ,ELTVAR,ELTPTR,NSUP,SVAR,
     &                 LIW,IW,LP,INFO)
      INTEGER LIW,LP,N,NELT,NSUP,NZ
      INTEGER INFO(6)
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER IW(LIW),SVAR(0:N)
      INTEGER FLAG,NEW,VARS
      EXTERNAL DMUMPS_316
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      INFO(4) = 0
      IF (N.LT.1) GO TO 10
      IF (NELT.LT.1) GO TO 20
      IF (NZ.LT.ELTPTR(NELT+1)-1) GO TO 30
      IF (LIW.LT.6) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      END IF
      NEW = 1
      VARS = NEW + LIW/3
      FLAG = VARS + LIW/3
      CALL DMUMPS_316(N,NELT,ELTPTR,NZ,ELTVAR,SVAR,NSUP,LIW/3-1,
     &           IW(NEW),IW(VARS),IW(FLAG),INFO)
      IF (INFO(1).EQ.-4) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      ELSE
         INFO(4) = 3*NSUP + 3
      END IF
      GO TO 50
   10 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   20 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   30 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   40 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) INFO(4)
      END IF
   50 RETURN
 9000 FORMAT (/3X,'Error message from DMUMPS_315: INFO(1) = ',I2)
 9010 FORMAT (3X,'LIW is insufficient. Upper bound on required work',
     &       'space is ',I8)
      END SUBROUTINE DMUMPS_315
      SUBROUTINE DMUMPS_316( N, NELT, ELTPTR, NZ, ELTVAR,
     &           SVAR, NSUP, MAXSUP, NEW, VARS, FLAG, INFO )
      INTEGER MAXSUP,N,NELT,NSUP,NZ
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER INFO(6)
      INTEGER FLAG(0:MAXSUP), NEW(0:MAXSUP),SVAR(0:N),
     &          VARS(0:MAXSUP)
      INTEGER I,IS,J,JS,K,K1,K2
      DO 10 I = 0,N
         SVAR(I) = 0
   10 CONTINUE
      VARS(0) = N + 1
      NEW(0) = -1
      FLAG(0) = 0
      NSUP = 0
      DO 40 J = 1,NELT
         K1 = ELTPTR(J)
         K2 = ELTPTR(J+1) - 1
         DO 20 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) THEN
               INFO(2) = INFO(2) + 1
               GO TO 20
            END IF
            IS = SVAR(I)
            IF (IS.LT.0) THEN
               ELTVAR(K) = 0
               INFO(3) = INFO(3) + 1
               GO TO 20
            END IF
            SVAR(I) = SVAR(I) - N - 2
            VARS(IS) = VARS(IS) - 1
   20    CONTINUE
         DO 30 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) GO TO 30
            IS = SVAR(I) + N + 2
            IF (FLAG(IS).LT.J) THEN
               FLAG(IS) = J
               IF (VARS(IS).GT.0) THEN
                  NSUP = NSUP + 1
                  IF (NSUP.GT.MAXSUP) THEN
                     INFO(1) = -4
                     RETURN
                  END IF
                  VARS(NSUP) = 1
                  FLAG(NSUP) = J
                  NEW(IS) = NSUP
                  SVAR(I) = NSUP
               ELSE
                  VARS(IS) = 1
                  NEW(IS) = IS
                  SVAR(I) = IS
               END IF
            ELSE
               JS = NEW(IS)
               VARS(JS) = VARS(JS) + 1
               SVAR(I) = JS
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_316
      SUBROUTINE DMUMPS_36( COMM_LOAD, ASS_IRECV,
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, IFLAG,
     &    IERROR, ND, 
     &    FILS, FRERE, DAD, MAXFRW, root,
     &    OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC, PTRAST, 
     &    STEP, PIMASTER, PAMASTER,PTRARW, 
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, 
     &    ICNTL, KEEP,KEEP8,INTARR,DBLARR, 
     &
     &    NSTK_S,NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER IZERO 
      PARAMETER (IZERO=0)
      INTEGER NELT,N,LIW,NSTEPS
      INTEGER(8) LA, LRLU, LRLUS, IPTRLU, POSFAC
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     &        IWPOS, IWPOSCB, COMP
      INTEGER IDUMMY(1)
      INTEGER IW(LIW), ITLOC(N+KEEP(253)),
     &        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)), PERM(N), 
     &        FILS(N), FRERE(KEEP(28)), DAD(KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28)),
     &              PAMASTER(KEEP(28))
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      LOGICAL SON_LEVEL2
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(max(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER IPOOL( LPOOL )
      INTEGER NBPROCFILS(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INTEGER ETATASS
      INCLUDE 'mumps_headers.h'
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      INTEGER  MUMPS_330
      EXTERNAL MUMPS_330
      INTEGER LP, HS, HF
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFS4FATHER
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER(8) NFRONT8
      INTEGER(8) LAELL8, APOS, APOS2, LAPOS2
      INTEGER(8) POSELT, POSEL1, ICT12, ICT21
      INTEGER(8) IACHK
      INTEGER(8) JJ2
      INTEGER(8) LSTK8, SIZFR8
#if defined(ALLOW_NON_INIT)
      INTEGER(8) :: JJ8
#endif
      INTEGER NBPANELS_L, NBPANELS_U, LREQ_OOC
      INTEGER SIZFI, NCB
      INTEGER JJ,J1,J2
      INTEGER NCOLS, NROWS, LDA_SON
      INTEGER NELIM,JJ1,J3,
     &        IORG, IBROT
      INTEGER JPOS,ICT11, IJROW
      INTEGER Pos_First_NUMORG,NUMORG,IOLDPS,
     &        NUMELT, ELBEG
      INTEGER AINPUT, 
     &        AII, J
      INTEGER NSLAVES, NSLSON, NPIVS, NPIV_ANA, NPIV
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, I
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER NCBSON
      LOGICAL SAME_PROC
      INTRINSIC real
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      LOGICAL  MUMPS_167, SSARBR
      EXTERNAL MUMPS_167
      DOUBLE PRECISION FLOP1,FLOP1_EFF
      EXTERNAL MUMPS_170
      LOGICAL MUMPS_170
      NFS4FATHER = -1
      ETATASS    = 0  
      COMPRESSCB=.FALSE.
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      LEVEL = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF (LEVEL.NE.1) THEN 
       write(6,*) 'Error1 in mpi51f_niv1 '
       CALL MUMPS_ABORT()
      END IF
      NSLAVES = 0
      HF = 6 + NSLAVES + KEEP(IXSZ)
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .ne. 0 ) THEN
        ELBEG  = FRT_PTR(INODE)
      ELSE
        ELBEG  = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NPIV_ANA=NUMORG
      NSTEPS = NSTEPS + 1
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 + KEEP(IXSZ))
         ISON = FRERE(STEP(ISON))
      END DO
      NFRONT = ND(STEP(INODE)) + NASS  + KEEP(253)
      NASS1 = NASS + NUMORG
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1) THEN 
        CALL DMUMPS_684(KEEP(50), NFRONT, NFRONT, NASS1,
     &       NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC   
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL DMUMPS_94(N, KEEP(28),
     &        IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ))
          COMP = COMP+1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          END IF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      END IF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .TRUE.
      IF (KEEP(50).EQ.0) THEN
        CALL  MUMPS_124(
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF, NFRONT, NFRONT_EFF,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     &        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     &        INTARR, KEEP(14), ITLOC, RHS_MUMPS, FILS, FRERE,
     &        KEEP,
     &        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG,
     &        DAD,PROCNODE_STEPS, SLAVEF, 
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG )
      ELSE
        CALL MUMPS_125( 
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     &        INTARR, KEEP(14), ITLOC, RHS_MUMPS, FILS, FRERE,
     &        KEEP,
     &        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG,
     &        DAD,PROCNODE_STEPS, SLAVEF,
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG)
        IF (IFLAG.LT.0) GOTO 300
      END IF
      IF (NFRONT_EFF.NE.NFRONT) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
           IF(MUMPS_170(PROCNODE_STEPS(STEP(INODE)),
     &          SLAVEF))THEN
              NPIV=NASS1-(NFRONT_EFF-ND(STEP(INODE)))
              CALL MUMPS_511(ND(STEP(INODE))+KEEP(253),
     &                                 NPIV,NPIV,
     &                                 KEEP(50),1,FLOP1)             
              NPIV=NPIV_ANA
              CALL MUMPS_511(ND(STEP(INODE))+KEEP(253),
     &                                 NPIV,NPIV,
     &                                 KEEP(50),1,FLOP1_EFF)
              CALL DMUMPS_190(0,.FALSE.,FLOP1-FLOP1_EFF,
     &             KEEP,KEEP8)
           ENDIF
        IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF))
        NFRONT = NFRONT_EFF
        LREQ = HF + 2 * NFRONT + LREQ_OOC   
        ELSE
         Write(*,*) ' ERROR 1 during ass_niv1_ELT'
         GOTO 270
        ENDIF
      ENDIF
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
        CALL DMUMPS_691(KEEP(50),
     &       NBPANELS_L, NBPANELS_U, NASS1, 
     &       IOLDPS + HF + 2 * NFRONT, IW, LIW)
      ENDIF
      NCB   = NFRONT - NASS1
      MAXFRW = max0(MAXFRW, NFRONT)
      ICT11 = IOLDPS + HF - 1 + NFRONT
      NFRONT8=int(NFRONT,8)
      LAELL8 = NFRONT8*NFRONT8
      IF (LRLU .LT. LAELL8) THEN
        IF (LRLUS .LT. LAELL8) THEN
          GOTO 280
        ELSE
          CALL DMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     &         LRLU, IPTRLU,
     &         IWPOS, IWPOSCB, PTRIST, PTRAST,
     &         STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &         KEEP(IXSZ))
          COMP = COMP + 1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          END IF
        END IF
      END IF
      LRLU = LRLU - LAELL8
      LRLUS = LRLUS - LAELL8
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL8
      SSARBR=MUMPS_167(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      CALL DMUMPS_471(SSARBR,.FALSE.,LA-LRLUS,0_8,LAELL8,
     &    KEEP,KEEP8,
     &     LRLU)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
        A(POSELT:LAPOS2) = ZERO
      ELSE
        APOS = POSELT
        DO JJ8 = 0_8, int(NFRONT -1,8)
          A(APOS:APOS+JJ8) = ZERO
          APOS = APOS + NFRONT8
        END DO
      END IF
#endif
      NASS = NASS1
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      PTLUST_S(STEP(INODE)) = IOLDPS
      IW(IOLDPS+XXI) = LREQ  
      CALL MUMPS_730(LAELL8,IW(IOLDPS+XXR))
      IW(IOLDPS+XXS) =-9999
      IW(IOLDPS+XXS+1:IOLDPS+KEEP(IXSZ)-1)=-99999
      IW(IOLDPS+KEEP(IXSZ))   = NFRONT
      IW(IOLDPS+KEEP(IXSZ)+ 1) = 0
      IW(IOLDPS+KEEP(IXSZ) + 2) = -NASS1
      IW(IOLDPS+KEEP(IXSZ) + 3) = -NASS1
      IW(IOLDPS+KEEP(IXSZ) + 4) = STEP(INODE)
      IW(IOLDPS+KEEP(IXSZ)+5)   = NSLAVES
      IF (NUMSTK.NE.0) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK    = PIMASTER(STEP(ISON))
          LSTK      = IW(ISTCHK+KEEP(IXSZ))
          LSTK8     = int(LSTK,8)
          NELIM     = IW(ISTCHK + 1+KEEP(IXSZ))
          NPIVS     = IW(ISTCHK + 3+KEEP(IXSZ))
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON    = IW(ISTCHK + 5+KEEP(IXSZ))
          HS        = 6 + NSLSON + KEEP(IXSZ)
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          IF ( SAME_PROC ) THEN
             COMPRESSCB =
     &           ( IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
             COMPRESSCB =  ( IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          LEVEL1    = NSLSON.EQ.0
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW( ISTCHK + 2+KEEP(IXSZ))
          ELSE
           NROWS = NCOLS
          ENDIF
          SIZFI   = HS + NROWS + NCOLS 
          J1 = ISTCHK + HS + NROWS + NPIVS
          IF ( .NOT. LEVEL1 .AND. NELIM.EQ.0 ) GOTO 205
          IF (LEVEL1) THEN
           J2 = J1 + LSTK - 1
           IF (COMPRESSCB) THEN
             SIZFR8 = (LSTK8*(LSTK8+1_8)/2_8)
           ELSE
             SIZFR8 = LSTK8*LSTK8
           ENDIF
          ELSE
           IF ( KEEP(50).eq.0 ) THEN
             SIZFR8 = int(NELIM,8) * LSTK8
           ELSE
             SIZFR8 = int(NELIM,8) * int(NELIM,8)
           END IF
           J2 = J1 + NELIM - 1
          ENDIF
          OPASSW = OPASSW + dble(SIZFR8)
          IACHK = PAMASTER(STEP(ISON))
          IF ( KEEP(50) .eq. 0 ) THEN
            POSEL1 = PTRAST(STEP(INODE)) - NFRONT8
            IF (J2.GE.J1) THEN
              DO 170 JJ = J1, J2
                APOS = POSEL1 + int(IW(JJ),8) * NFRONT8
                DO 160 JJ1 = 1, LSTK
                  JJ2 = APOS + int(IW(J1 + JJ1 - 1) - 1,8)
                  A(JJ2) = A(JJ2) + A(IACHK + int(JJ1 - 1,8))
  160           CONTINUE
                IACHK = IACHK + LSTK8
  170         CONTINUE
            END IF
          ELSE
            IF (LEVEL1) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = SIZFR8
            ELSE
              LCB = int(LDA_SON,8)*int(J2 - J1 + 1,8)
            ENDIF
            CALL DMUMPS_178(A, LA,
     &           PTRAST(STEP( INODE )), NFRONT, NASS1,
     &           IACHK, LDA_SON, LCB,
     &           IW( J1 ), J2 - J1 + 1, NELIM, ETATASS, 
     &           COMPRESSCB,
     &           .FALSE. 
     &          )
          ENDIF
  205     IF (LEVEL1) THEN 
           IF (SAME_PROC) ISTCHK = PTRIST(STEP(ISON))
           IF (SAME_PROC) THEN
             IF (KEEP(50).NE.0) THEN
              J2 = J1 + LSTK - 1
              DO JJ = J1, J2
               IW(JJ) = IW(JJ - NROWS)
              END DO
             ELSE
              J2 = J1 + LSTK - 1
              J3 = J1 + NELIM
              DO JJ = J3, J2
               IW(JJ) = IW(JJ - NROWS)
              END DO
              IF (NELIM .NE. 0) THEN
                J3 = J3 - 1
                DO JJ = J1, J3
                 JPOS = IW(JJ) + ICT11
                 IW(JJ) = IW(JPOS)
                END DO
              ENDIF
             ENDIF
           ENDIF
           IF ( SAME_PROC ) THEN
               PTRIST(STEP( ISON )) = -99999999
           ELSE
               PIMASTER(STEP( ISON )) = -99999999
           ENDIF
           CALL DMUMPS_152(SSARBR, MYID, N, ISTCHK,
     &     IACHK,
     &     IW, LIW, LRLU, LRLUS, IPTRLU,
     &     IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &     )
          ELSE
           PDEST = ISTCHK + 6 + KEEP(IXSZ)
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL MUMPS_49( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
               CALL DMUMPS_210( COMM_LOAD, ASS_IRECV, 
     &                           BUFR, LBUFR, LBUFR_BYTES,
     &                           INODE, ISON, NSLAVES, IDUMMY,
     &                           NFRONT, NASS1,NFS4FATHER,
     &         TROW_SIZE, IW( INDX ),
     &         PROCNODE_STEPS,
     &         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &         LRLUS, N, IW,
     &         LIW, A, LA, 
     &         PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &         PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL,
     &         LEAF, NBFIN, ICNTL, KEEP,KEEP8, root,
     &         OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ND, FRERE,
     &         NELT+1, NELT, FRT_PTR, FRT_ELT,
     &   
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &         )
               IF ( IFLAG .LT. 0 ) GOTO 500
               EXIT
             ENDIF
           END DO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  DMUMPS_71( INODE, NFRONT, 
     &       NASS1, NFS4FATHER,ISON, MYID,
     &       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     &       COMM, IERR, IW(PDEST), NSLSON, 
     &       SLAVEF, 
     &       KEEP,KEEP8, STEP, N, 
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
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
     &         INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &         NELT+1, NELT, FRT_PTR, FRT_ELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
               IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
           ENDIF
          ENDIF
      ISON = FRERE(STEP(ISON))
  220 CONTINUE
      END IF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          AINPUT    = AII + II - J1
          ICT12 = POSELT + int(I-1,8) * NFRONT8
          DO JJ=J1,J2
           APOS2 = ICT12 + int(INTARR(JJ) - 1,8)
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
           AINPUT = AINPUT + SIZE_ELTI
          END DO
         ELSE
          ICT12 = POSELT + int(- NFRONT + I - 1,8)
          ICT21 = POSELT + int(I-1,8)*NFRONT8 - 1_8
          DO JJ=II,J2
           J =  INTARR(JJ)
           IF (I.LT.J) THEN
              APOS2 = ICT12 + int(J,8)*NFRONT8
           ELSE
              APOS2 = ICT21 + int(J,8)
           ENDIF
           A(APOS2) = A(APOS2) + DBLARR(AII)
           AII = AII + 1
          END DO
         END IF
        END DO
      END DO
      IF (KEEP(253).GT.0) THEN
       POSELT = PTRAST(STEP(INODE))
       IBROT = INODE
       IJROW = Pos_First_NUMORG  
       DO IORG = 1, NUMORG
        IF (KEEP(50).EQ.0) THEN
          DO JJ=1, KEEP(253)
            APOS = POSELT+
     &           int(IJROW-1,8)               * NFRONT8 +
     &           int(NFRONT-KEEP(253)+JJ-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ELSE
          DO JJ=1, KEEP(253)
            APOS = POSELT+
     &           int(NFRONT-KEEP(253)+JJ-1,8)  * NFRONT8 +
     &           int(IJROW-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ENDIF
       IBROT = FILS(IBROT)
       IJROW = IJROW+1
       ENDDO
      ENDIF
      GOTO 500
  270 CONTINUE
      IFLAG = -8
      IERROR = LREQ
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE IN INTEGER ALLOCATION DURING DMUMPS_36'
      ENDIF
      GOTO 490
  280 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, WORKSPACE TOO SMALL DURING DMUMPS_36'
      ENDIF
      IFLAG = -9
      CALL MUMPS_731(LAELL8-LRLUS, IERROR)
      GOTO 500
  290 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &  ' FAILURE, SEND BUFFER TOO SMALL DURING DMUMPS_36'
      ENDIF
      IFLAG = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &  ' FAILURE, RECV BUFFER TOO SMALL DURING DMUMPS_36'
      ENDIF
      IFLAG = -20
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION DURING DMUMPS_36'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK 
  490 CALL  DMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_36
      SUBROUTINE DMUMPS_37( COMM_LOAD, ASS_IRECV,
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, IFLAG,
     &    IERROR, ND, FILS, FRERE, DAD,
     &    CAND, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW, root,
     &    OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_S,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS, 
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, 
     &    ICNTL, KEEP,KEEP8,INTARR,DBLARR, 
     &    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, LEAF, IPOOL, LPOOL,
     &    PERM,
     &    MEM_DISTRIB)
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER NELT, N,LIW,NSTEPS, NBFIN
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER(8) LRLU, IPTRLU, LRLUS, POSFAC, LA
      INTEGER(8) LAELL8
      INTEGER JJ
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     &        LPOOL, LEAF, 
     &        IWPOS, 
     &        IWPOSCB, COMP, SLAVEF
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER IW(LIW), ITLOC(N+KEEP(253)),
     &        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)),
     &        FILS(N), FRERE(KEEP(28)), STEP(N), DAD (KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     & PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), PERM(N)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRFAC(KEEP(28)), PAMASTER(KEEP(28)),
     &              PTRAST(KEEP(28))
      INTEGER   CAND(SLAVEF+1,max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(max(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      INTEGER MYID, COMM
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER NBPROCFILS(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INCLUDE 'mumps_headers.h'
      INTEGER LP, HS, HF, HF_OLD, NSLAVES_OLD,NCBSON
      INTEGER NCBSON_MAX
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      INTEGER NFS4FATHER
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER LREQ_OOC, NBPANELS_L, NBPANELS_U
      INTEGER NCB
      INTEGER J1,J2
      INTEGER(8) NFRONT8, LAPOS2, POSELT, POSEL1, LDAFS8,
     &           JJ2, IACHK, ICT12, ICT21
#if defined(ALLOW_NON_INIT)
      INTEGER(8) :: JJ8
#endif
      INTEGER(8) APOS, APOS2
      INTEGER NELIM,LDAFS,JJ1,NPIVS,NCOLS,NROWS,
     &        IORG
      INTEGER LDA_SON, IJROW, IBROT
      INTEGER Pos_First_NUMORG,NBCOL,NUMORG,IOLDPS
      INTEGER AINPUT
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, ISLAVE, PDEST
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, ELBEG, NUMELT, I, J, AII
      LOGICAL SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      logical :: force_cand
      INTEGER(8) APOSMAX
      DOUBLE PRECISION  MAXARR
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok
      INTEGER NUMORG_SPLIT, TYPESPLIT,
     &        NCB_SPLIT, SIZE_LIST_SPLIT, NBSPLIT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST, COPY_CAND
      INTEGER IZERO 
      INTEGER IDUMMY(1)
      INTEGER PDEST1(1)
      INTEGER ETATASS
      PARAMETER( IZERO = 0 )
      INTEGER MUMPS_275, MUMPS_330, MUMPS_810
      EXTERNAL MUMPS_275, MUMPS_330, MUMPS_810
      INTRINSIC real
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION RZERO
      PARAMETER( RZERO = 0.0D0 )
      PARAMETER( ZERO = 0.0D0 )
      COMPRESSCB=.FALSE.
      ETATASS = 0  
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      NSTEPS = NSTEPS + 1
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .NE. 0 ) THEN
        ELBEG = FRT_PTR(INODE)
      ELSE
        ELBEG = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      NCBSON_MAX  = 0
      DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         IF ( KEEP(48)==5 .AND. 
     &        MUMPS_330(PROCNODE_STEPS(STEP(ISON)),
     &        SLAVEF) .EQ. 1) THEN
            NCBSON_MAX =
     &      max(NCBSON_MAX,IW(PIMASTER(STEP(ISON))+KEEP(IXSZ)))
         END IF
         NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 + KEEP(IXSZ))
         ISON = FRERE(STEP(ISON))
      END DO
      NFRONT = ND(STEP(INODE)) + NASS + KEEP(253)
      MAXFRW = max0(MAXFRW, NFRONT)
      NASS1 = NASS + NUMORG
      NCB   = NFRONT - NASS1
      IF((KEEP(24).eq.0).or.(KEEP(24).eq.1)) then
         force_cand=.FALSE.
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      end if
      IF (force_cand) THEN
         INIV2 = ISTEP_TO_INIV2( STEP( INODE ))
         SIZE_TMP_SLAVES_LIST = CAND( SLAVEF+1, INIV2 )
      ELSE
         INIV2 = 1
         SIZE_TMP_SLAVES_LIST = SLAVEF - 1
      ENDIF
      ALLOCATE(TMP_SLAVES_LIST(SIZE_TMP_SLAVES_LIST),stat=allocok)
      IF (allocok > 0 ) THEN
        GOTO 265
      ENDIF
       TYPESPLIT =  MUMPS_810 (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
       IF  ( (TYPESPLIT.EQ.4) 
     &               .OR.(TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) 
     &     )  THEN
        IF (TYPESPLIT.EQ.4) THEN
         ALLOCATE(COPY_CAND(SLAVEF+1),stat=allocok)
         IF (allocok > 0 ) THEN
           GOTO 265
         ENDIF
         CALL DMUMPS_791 (
     &      INODE, STEP, N, SLAVEF, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      CAND(1,INIV2), ICNTL, COPY_CAND,
     &      NBSPLIT, NUMORG_SPLIT, TMP_SLAVES_LIST(1),
     &      SIZE_TMP_SLAVES_LIST 
     &                                    )
         NCB_SPLIT = NCB-NUMORG_SPLIT
         SIZE_LIST_SPLIT = SIZE_TMP_SLAVES_LIST - NBSPLIT
         CALL DMUMPS_472( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, COPY_CAND,
     &     MEM_DISTRIB(0), NCB_SPLIT, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST(NBSPLIT+1),
     &     SIZE_LIST_SPLIT,INODE )
         DEALLOCATE (COPY_CAND)
         CALL DMUMPS_790 (
     &      INODE, STEP, N, SLAVEF, NBSPLIT, NCB,
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      ICNTL, 
     &      TAB_POS_IN_PERE(1,INIV2),
     &      NSLAVES
     &                                    )
        ELSE
         ISTCHK    = PIMASTER(STEP(IFSON))
         PDEST     = ISTCHK + 6 + KEEP(IXSZ)
         NSLSON    = IW(ISTCHK + KEEP(IXSZ) + 5)
         TMP_SLAVES_LIST (1:SIZE_TMP_SLAVES_LIST) 
     &             = CAND (1:SIZE_TMP_SLAVES_LIST, INIV2)
         CALL DMUMPS_792 (
     &      INODE, TYPESPLIT, IFSON, 
     &      IW(PDEST), NSLSON,
     &      STEP, N, SLAVEF, NBSPLIT, NCB,
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      ICNTL, ISTEP_TO_INIV2, INIV2,
     &      TAB_POS_IN_PERE, NSLAVES, 
     &      TMP_SLAVES_LIST,
     &      SIZE_TMP_SLAVES_LIST
     &                                    )
        ENDIF
       ELSE
        CALL DMUMPS_472( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, CAND(1,INIV2),
     &     MEM_DISTRIB(0), NCB, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST,
     &     SIZE_TMP_SLAVES_LIST,INODE )
       ENDIF
      HF   = NSLAVES + 6 + KEEP(IXSZ)
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1) THEN 
        CALL DMUMPS_684(KEEP(50), NASS1, NFRONT, NASS1,
     &                               NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL DMUMPS_94(N, KEEP(28),
     &        IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ))
          COMP = COMP+1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv2'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          ENDIF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      ENDIF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .FALSE.
      IF (KEEP(50).EQ.0) THEN
        CALL  MUMPS_124(
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF, NFRONT,NFRONT_EFF,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     &        INTARR, KEEP(14), ITLOC, RHS_MUMPS, FILS, FRERE,
     &        KEEP, SON_LEVEL2, NIV1, NBPROCFILS, IFLAG, 
     &        DAD,PROCNODE_STEPS, SLAVEF,
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG)
      ELSE
        CALL MUMPS_125(
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     &        INTARR, KEEP(14), ITLOC, RHS_MUMPS, FILS, FRERE,
     &        KEEP, SON_LEVEL2, NIV1, NBPROCFILS, IFLAG,
     &        DAD,PROCNODE_STEPS, SLAVEF,
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG)
        IF (IFLAG.LT.0) GOTO 250
      ENDIF
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (
     &        (TYPESPLIT.EQ.5) .OR. (TYPESPLIT.EQ.6)) THEN
          WRITE(6,*) ' SPLITTING NOT YET READY FOR THAT'
          CALL MUMPS_ABORT()
        ENDIF
        IF (NFRONT.GT.NFRONT_EFF) THEN
            NCB    = NFRONT_EFF - NASS1
            NSLAVES_OLD = NSLAVES
            HF_OLD      = HF
            IF (TYPESPLIT.EQ.4) THEN
             WRITE(6,*) ' Internal error 2 in fac_ass_elt due',
     &     ' to plitting ',
     &     ' INODE, NFRONT, NFRONT_EFF =', INODE, NFRONT, NFRONT_EFF 
             CALL MUMPS_ABORT()
            ELSE
             CALL DMUMPS_472( NCBSON_MAX, 
     &        SLAVEF, KEEP,KEEP8,ICNTL,
     &        CAND(1,INIV2),
     &        MEM_DISTRIB(0), NCB, NFRONT_EFF, NSLAVES,
     &        TAB_POS_IN_PERE(1,INIV2),
     &        TMP_SLAVES_LIST, SIZE_TMP_SLAVES_LIST,INODE )
            ENDIF
            HF = NSLAVES + 6 + KEEP(IXSZ)
            IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF)) -
     &                   (NSLAVES_OLD - NSLAVES)
            IF (NSLAVES_OLD .NE. NSLAVES) THEN
              IF (NSLAVES_OLD > NSLAVES) THEN
               DO JJ=0,2*NFRONT_EFF-1
                 IW(IOLDPS+HF+JJ)=IW(IOLDPS+HF_OLD+JJ)
               ENDDO
              ELSE
               IF (IWPOS - 1 > IWPOSCB ) GOTO 270
               DO JJ=2*NFRONT_EFF-1, 0, -1
                 IW(IOLDPS+HF+JJ) = IW(IOLDPS+HF_OLD+JJ)
               ENDDO
              END IF
            END IF
            NFRONT = NFRONT_EFF
            LREQ = HF + 2 * NFRONT + LREQ_OOC
        ELSE
          Write(*,*) ' ERROR 2 during ass_niv2'
          GOTO 270
        ENDIF
      ENDIF
      NFRONT8=int(NFRONT,8)
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
        CALL DMUMPS_691(KEEP(50),
     &       NBPANELS_L, NBPANELS_U, NASS1, 
     &       IOLDPS + HF + 2 * NFRONT, IW, LIW)
      ENDIF
      MAXFRW = max0(MAXFRW, NFRONT)
      PTLUST_S(STEP(INODE)) = IOLDPS
      IW(IOLDPS + 1+KEEP(IXSZ)) = 0
      IW(IOLDPS + 2+KEEP(IXSZ)) = -NASS1
      IW(IOLDPS + 3+KEEP(IXSZ)) = -NASS1
      IW(IOLDPS + 4+KEEP(IXSZ)) = STEP(INODE)
      IW(IOLDPS+KEEP(IXSZ))   = NFRONT
      IW(IOLDPS+5+KEEP(IXSZ)) = NSLAVES
      IW(IOLDPS+6+KEEP(IXSZ):IOLDPS+5+NSLAVES+KEEP(IXSZ))=
     &                     TMP_SLAVES_LIST(1:NSLAVES)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE)
      IF (KEEP(73) .EQ. 0) THEN
#endif
#endif
        CALL DMUMPS_461(MYID, SLAVEF, COMM_LOAD,
     &     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &     NASS1, KEEP, KEEP8, IW(IOLDPS+6+KEEP(IXSZ)), NSLAVES,INODE)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE) 
      ENDIF
#endif
#endif
      IF(KEEP(86).EQ.1)THEN
         IF(mod(KEEP(24),2).eq.0)THEN
            CALL DMUMPS_533(SLAVEF,CAND(SLAVEF+1,INIV2),
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, NSLAVES,INODE)
         ELSEIF((KEEP(24).EQ.0).OR.(KEEP(24).EQ.1))THEN
            CALL DMUMPS_533(SLAVEF,SLAVEF-1,
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8,TMP_SLAVES_LIST, NSLAVES,INODE)
         ENDIF
      ENDIF         
      DEALLOCATE(TMP_SLAVES_LIST)
      IF (KEEP(50).EQ.0) THEN
        LAELL8 = int(NASS1,8) * NFRONT8
        LDAFS = NFRONT
        LDAFS8 = NFRONT8
      ELSE
        LAELL8 = int(NASS1,8)*int(NASS1,8)
        IF (KEEP(219).NE.0) THEN
          IF(KEEP(50) .EQ. 2) LAELL8 = LAELL8+int(NASS1,8)
        ENDIF
        LDAFS = NASS1
        LDAFS8 = int(NASS1,8)
      ENDIF
      IF (LRLU .LT. LAELL8) THEN
        IF (LRLUS .LT. LAELL8) THEN
          GOTO 280
        ELSE
          CALL DMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ))
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv2'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          ENDIF
        ENDIF
      ENDIF
      LRLU = LRLU - LAELL8
      LRLUS = LRLUS - LAELL8
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSELT = POSFAC
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      POSFAC = POSFAC + LAELL8
      IW(IOLDPS+XXI)   = LREQ  
      CALL MUMPS_730(LAELL8,IW(IOLDPS+XXR)) 
      IW(IOLDPS+XXS) =-9999
      IW(IOLDPS+XXS+1:IOLDPS+KEEP(IXSZ)-1)=-99999
      CALL DMUMPS_471(.FALSE.,.FALSE.,LA-LRLUS,0_8,LAELL8,
     & KEEP,KEEP8,
     &LRLU)
      POSEL1 = POSELT - LDAFS8
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
        A(POSELT:LAPOS2) = ZERO
      ELSE
        APOS = POSELT
        DO JJ8 = 0_8, LDAFS8 - 1_8
          A(APOS:APOS+JJ8) = ZERO
          APOS = APOS + LDAFS8
        END DO
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
          A(APOS:APOS+LDAFS8-1_8)=ZERO
        ENDIF
      END IF
#endif
      IF ((NUMSTK.NE.0).AND.(NASS.NE.0)) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK = PIMASTER(STEP(ISON))
          NELIM = IW(ISTCHK + KEEP(IXSZ) + 1)
          IF (NELIM.EQ.0) GOTO 210
          LSTK    = IW(ISTCHK + KEEP(IXSZ))
          NPIVS   = IW(ISTCHK + KEEP(IXSZ) + 3)
          IF (NPIVS.LT.0) NPIVS=0
          NSLSON  = IW(ISTCHK + KEEP(IXSZ) + 5)
          HS      = 6 + KEEP(IXSZ) + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          IF ( SAME_PROC ) THEN
            COMPRESSCB=( IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
            COMPRESSCB=( IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW(ISTCHK + KEEP(IXSZ) + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          OPASSW = OPASSW + dble(NELIM*LSTK)
          J1 = ISTCHK + HS + NROWS + NPIVS
          J2 = J1 + NELIM - 1
          IACHK = PAMASTER(STEP(ISON))
          IF (KEEP(50).eq.0) THEN
           DO 170 JJ = J1, J2
            APOS = POSEL1 + int(IW(JJ),8) * LDAFS8
            DO 160 JJ1 = 1, LSTK
              JJ2 = APOS + int(IW(J1 + JJ1 - 1),8) - 1_8
              A(JJ2) = A(JJ2) + A(IACHK + int(JJ1,8) - 1_8)
  160       CONTINUE
            IACHK = IACHK + int(LSTK,8)
  170      CONTINUE
          ELSE
            IF (NSLSON.EQ.0) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = (int(NELIM,8)*int(NELIM+1,8))/2_8
            ELSE
              LCB = int(LDA_SON,8)*int(NELIM,8)
            ENDIF
            CALL DMUMPS_178(A, LA,
     &           POSELT, LDAFS, NASS1,
     &           IACHK, LDA_SON, LCB,
     &           IW( J1 ), NELIM, NELIM, ETATASS,
     &           COMPRESSCB,
     &           .FALSE. 
     &          )
          ENDIF
  210     ISON = FRERE(STEP(ISON))
  220   CONTINUE
      ENDIF
      APOSMAX = POSELT + int(NASS1,8)*int(NASS1,8)
      IF (KEEP(219).NE.0) THEN
        IF (KEEP(50).EQ.2) THEN
          A( APOSMAX: APOSMAX+int(NASS1-1,8))=ZERO
        ENDIF
      ENDIF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          IF (I.LE.NASS1) THEN
           AINPUT    = AII + II - J1
           ICT12 = POSELT + int(I-1,8) * LDAFS8
           DO JJ=J1,J2
            APOS2 = ICT12 + int(INTARR(JJ) - 1,8)
            A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            AINPUT = AINPUT + SIZE_ELTI
           END DO
          ENDIF
         ELSE
          ICT12 = POSELT - LDAFS8 + int(I,8) - 1_8
          ICT21 = POSELT + int(I-1,8)*LDAFS8 - 1_8
          IF ( I .GT. NASS1 ) THEN
           IF (KEEP(219).NE.0) THEN
            IF (KEEP(50).EQ.2) THEN
              AINPUT=AII
              DO JJ=II,J2
               J=INTARR(JJ)
               IF (J.LE.NASS1) THEN
                A(APOSMAX+int(J-1,8))=
     &              max(dble(A(APOSMAX+int(J-1,8))),
     &                  abs(DBLARR(AINPUT)))
               ENDIF
               AINPUT=AINPUT+1
              ENDDO
            ELSE
              AII = AII + J2 - II + 1
              CYCLE
            ENDIF
           ELSE
             AII = AII + J2 - II + 1
             CYCLE
           ENDIF
          ELSE
            IF (KEEP(219).NE.0) THEN
              MAXARR = RZERO
            ENDIF
            DO JJ=II,J2
              J =  INTARR(JJ)
              IF ( J .LE. NASS1) THEN
                IF (I.LT.J) THEN
                  APOS2 = ICT12 + int(J,8)*LDAFS8
                ELSE
                  APOS2 = ICT21 + int(J,8)
                ENDIF
                A(APOS2) = A(APOS2) + DBLARR(AII)
              ELSE IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                MAXARR = max(MAXARR,abs(DBLARR(AII)))
              ENDIF
              AII = AII + 1
            END DO
            IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
                A(APOSMAX+int(I-1,8)) = 
     &             max( MAXARR, dble(A(APOSMAX+int(I-1,8))))
            ENDIF
          ENDIF 
         END IF 
        END DO
      END DO
      IF (KEEP(253).GT.0) THEN
       POSELT = PTRAST(STEP(INODE))
       IBROT = INODE
       IJROW = Pos_First_NUMORG  
       DO IORG = 1, NUMORG
        IF (KEEP(50).EQ.0) THEN
          DO JJ = 1, KEEP(253)
            APOS = POSELT +
     &             int(IJROW-1,8) * int(LDAFS,8) +
     &             int(LDAFS-KEEP(253)+JJ-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ENDIF
        IBROT = FILS(IBROT)
        IJROW = IJROW+1
       ENDDO
      ENDIF
      PTRCOL = IOLDPS + HF + NFRONT 
      PTRROW = IOLDPS + HF + NASS1 
      PDEST  = IOLDPS + 6 + KEEP(IXSZ)
      DO ISLAVE = 1, NSLAVES
              CALL MUMPS_49( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB,
     &                NSLAVES, 
     &                NBLIG, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
        IERR = -1
        DO WHILE (IERR .EQ.-1)
         IF ( KEEP(50) .eq. 0 ) THEN
           NBCOL =  NFRONT
           CALL DMUMPS_68( INODE,
     &      NBPROCFILS(STEP(INODE)),
     &      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     &      IZERO, IDUMMY,
     &      IW(PDEST), NFRONT, COMM, IERR)
         ELSE
           NBCOL = NASS1+SHIFT_INDEX+NBLIG
           CALL DMUMPS_68( INODE,
     &      NBPROCFILS(STEP(INODE)),
     &      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     &      NSLAVES-ISLAVE,
     &      IW( PTLUST_S(STEP(INODE))+6+KEEP(IXSZ)+ISLAVE),
     &      IW(PDEST), NFRONT, COMM, IERR)
         ENDIF
         IF (IERR.EQ.-1) THEN
          BLOCKING  = .FALSE.
          SET_IRECV = .TRUE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &     MPI_ANY_SOURCE, MPI_ANY_TAG,
     &     STATUS, BUFR, LBUFR,
     &     LBUFR_BYTES,
     &     PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &     LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &     PTLUST_S, PTRFAC,
     &     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &     IERROR, COMM,
     &     NBPROCFILS,
     &     IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &     root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &     FILS, PTRARW, PTRAIW,
     &     INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &     NELT+1, NELT, FRT_PTR, FRT_ELT,
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.)
           IF ( IFLAG .LT. 0 ) GOTO 500
          IF (MESSAGE_RECEIVED) THEN
           IOLDPS = PTLUST_S(STEP(INODE))
           PTRCOL = IOLDPS + HF + NFRONT
           PTRROW = IOLDPS + HF + NASS1 + SHIFT_INDEX
          ENDIF
         ENDIF
        END DO
        IF (IERR .EQ. -2) GOTO 300
        IF (IERR .EQ. -3) GOTO 305
        PTRROW = PTRROW + NBLIG
        PDEST  = PDEST + 1
      END DO
      IF (NUMSTK.EQ.0) GOTO 500
      ISON = IFSON
      DO IELL = 1, NUMSTK
        ISTCHK = PIMASTER(STEP(ISON))
        NELIM = IW(ISTCHK + 1 + KEEP(IXSZ))
        LSTK    = IW(ISTCHK + KEEP(IXSZ))
        NPIVS   = IW(ISTCHK + 3 + KEEP(IXSZ))
        IF ( NPIVS .LT. 0 ) NPIVS = 0
        NSLSON  = IW(ISTCHK + 5 + KEEP(IXSZ))
        HS      = 6 + NSLSON + KEEP(IXSZ)
        NCOLS     = NPIVS + LSTK
        SAME_PROC     = (ISTCHK.LE.IWPOS)
        IF (.NOT.SAME_PROC) THEN
         NROWS = IW(ISTCHK + 2 + KEEP(IXSZ) )
        ELSE
         NROWS = NCOLS
        ENDIF
        PDEST   = ISTCHK + 6 + KEEP(IXSZ)
        NCBSON  = LSTK - NELIM
        PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
           NFS4FATHER = NCBSON
           DO I=0,NCBSON-1
              IF(IW(PTRCOL+I) .GT. NASS1) THEN
                 NFS4FATHER = I
                 EXIT
              ENDIF
           ENDDO
           NFS4FATHER=NFS4FATHER + NELIM
        ELSE
          NFS4FATHER = 0
        ENDIF
        IF (NSLSON.EQ.0) THEN
          NSLSON = 1
          PDEST1(1)  = MUMPS_275(PROCNODE_STEPS(STEP(ISON)),
     &                                SLAVEF)
          IF (PDEST1(1).EQ.MYID) THEN
            CALL DMUMPS_211( COMM_LOAD, ASS_IRECV, 
     &      BUFR, LBUFR, LBUFR_BYTES,
     &      INODE, ISON, NSLAVES, 
     &      IW( PTLUST_S(STEP(INODE)) + 6 +KEEP(IXSZ)),
     &      NFRONT, NASS1, NFS4FATHER,NCBSON, IW( PTRCOL ),
     &      PROCNODE_STEPS,
     &      SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &      LRLUS, N, IW, LIW, A, LA, 
     &      PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &      PIMASTER, PAMASTER, NSTK_S, COMP,
     &      IFLAG, IERROR, MYID, COMM, NBPROCFILS, 
     &      IPOOL, LPOOL, LEAF, NBFIN, ICNTL, KEEP,KEEP8, root,
     &      OPASSW, OPELIW,
     &      ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &      ND, FRERE, NELT+1, NELT, 
     &      FRT_PTR, FRT_ELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE )
           IF ( IFLAG .LT. 0 ) GOTO 500
          ELSE
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            CALL  DMUMPS_71( 
     &           INODE, NFRONT,NASS1,NFS4FATHER,
     &           ISON, MYID,
     &      NSLAVES, IW( PTLUST_S(STEP(INODE)) + 6 +KEEP(IXSZ)),
     &      IW(PTRCOL), NCBSON,
     &      COMM, IERR, PDEST1, NSLSON, SLAVEF, 
     &      KEEP,KEEP8, STEP, N, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &        MPI_ANY_SOURCE, MPI_ANY_TAG,
     &        STATUS, BUFR, LBUFR, LBUFR_BYTES,
     &        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &        LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &        PTLUST_S, PTRFAC,
     &        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &        IERROR, COMM,
     &        NBPROCFILS,
     &        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &        root,OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS,
     &        PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.)
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
          ENDIF
        ELSE
          DO ISLAVE = 0, NSLSON-1
            IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
               CALL MUMPS_49( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON,
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX        = PTRCOL + SHIFT_INDEX
              CALL DMUMPS_210( COMM_LOAD, ASS_IRECV, 
     &        BUFR, LBUFR, LBUFR_BYTES,
     &        INODE, ISON, NSLAVES, 
     &        IW( PTLUST_S(STEP(INODE))+6+KEEP(IXSZ)),
     &        NFRONT, NASS1,NFS4FATHER,
     &        TROW_SIZE, IW( INDX ),
     &        PROCNODE_STEPS,
     &        SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &        LRLUS, N, IW,
     &        LIW, A, LA, 
     &        PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &        PIMASTER, PAMASTER, NSTK_S, COMP,
     &        IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     &        NBFIN, ICNTL, KEEP,KEEP8, root,
     &        OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &        INTARR, DBLARR, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     & 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &        )
              IF ( IFLAG .LT. 0 ) GOTO 500
              EXIT
            ENDIF
          END DO
          IF (PIMASTER(STEP(ISON)).GT.0) THEN
          IERR = -1
          DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  DMUMPS_71( 
     &           INODE, NFRONT, NASS1, NFS4FATHER,
     &           ISON, MYID,
     &      NSLAVES, IW(PTLUST_S(STEP(INODE))+6+KEEP(IXSZ)),
     &      IW(PTRCOL), NCBSON,
     &      COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     &      KEEP,KEEP8, STEP, N, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &        MPI_ANY_SOURCE, MPI_ANY_TAG,
     &        STATUS, BUFR, LBUFR,
     &        LBUFR_BYTES,
     &        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &        LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &        PTLUST_S, PTRFAC,
     &        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &        IERROR, COMM,
     &        NBPROCFILS,
     &        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &        root,OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &        FILS, PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
             IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
          END DO
          IF (IERR .EQ. -2) GOTO 290
          IF (IERR .EQ. -3) GOTO 295
          ENDIF
        ENDIF
       ISON = FRERE(STEP(ISON))
      END DO
      GOTO 500
  250 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION during assembly'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK + 1
      GOTO 490
  265 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE ALLOCATING TMP_SLAVES_LIST',
     &                 ' DYNAMIC ALLOCATION during assembly'
      ENDIF
      IFLAG  = -13
      IERROR = SIZE_TMP_SLAVES_LIST
      GOTO 490
  270 CONTINUE
      IFLAG = -8
      IERROR = LREQ
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &  ' FAILURE IN INTEGER ALLOCATION DURING DMUMPS_37'
      ENDIF
      GOTO 490
  280 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &  ' FAILURE, WORKSPACE TOO SMALL DURING DMUMPS_37'
      ENDIF
      IFLAG = -9
      CALL MUMPS_731(LAELL8 - LRLUS, IERROR)
      GOTO 490
  290 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (1) DURING DMUMPS_37'
      ENDIF
      IFLAG = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (1) DURING DMUMPS_37'
      ENDIF
      IFLAG = -20
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, SENDBUFFER TOO SMALL (2) DURING DMUMPS_37'
      ENDIF
      IFLAG = -17
      LREQ = NBLIG + NBCOL + 4+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  305 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, RECVBUFFER TOO SMALL (2) DURING DMUMPS_37'
      ENDIF
      IFLAG = -17
      LREQ = NBLIG + NBCOL + 4+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  490 CALL DMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_37
      SUBROUTINE DMUMPS_123(
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, 
     &    NBROWS, NBCOLS,
     &    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     &    RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR, 
     &    ICNTL, KEEP, KEEP8, MYID)
      IMPLICIT NONE
      INTEGER NELT, N,LIW
      INTEGER(8) :: LA
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER INODE, MYID
      INTEGER NBROWS, NBCOLS 
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER IW(LIW), ITLOC(N + KEEP(253)), STEP(N),
     &        PTRIST(KEEP(28)),
     &        FILS(N), PTRARW(NELT+1), 
     &        PTRAIW(NELT+1)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      DOUBLE PRECISION A(LA),
     &        DBLARR(max(1,KEEP(13)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER(8) :: POSELT, APOS2, ICT12, APOS
      INTEGER IOLDPS, NBCOLF, NBROWF, NSLAVES, HF,
     &        K1,K2,K,I,J,JPOS,NASS,JJ,
     &        IN,AINPUT,J1,J2,IJROW,ILOC, 
     &        ELBEG, NUMELT, ELTI, SIZE_ELTI, IPOS, 
     &        IPOS1, IPOS2, AII, II, IELL
      INTEGER    :: K1RHS, K2RHS, JFirstRHS
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INCLUDE 'mumps_headers.h'
      IOLDPS  = PTRIST(STEP(INODE))
      POSELT  = PTRAST(STEP(INODE))
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NASS    = IW(IOLDPS+1+KEEP(IXSZ))
      NSLAVES = IW(IOLDPS+5+KEEP(IXSZ))
      HF      = 6 + NSLAVES+KEEP(IXSZ)
      IF (NASS.LT.0) THEN
          NASS         = -NASS
          IW(IOLDPS+1+KEEP(IXSZ)) = NASS
          A(POSELT:POSELT+int(NBROWF,8)*int(NBCOLF,8)-1_8) =
     &    ZERO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -JPOS
           JPOS     = JPOS + 1
          END DO
          K1 = IOLDPS + HF 
          K2 = K1 + NBROWF - 1
          JPOS = 1
          IF ((KEEP(253).GT.0).AND.(KEEP(50).NE.0)) THEN
           K1RHS = 0
           K2RHS = -1
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
            IF ((K1RHS.EQ.0).AND.(J.GT.N)) THEN
             K1RHS = K
             JFirstRHS=J-N 
            ENDIF
            JPOS     = JPOS + 1
           ENDDO
           IF (K1RHS.GT.0) K2RHS=K2
           IF ( K2RHS.GE.K1RHS ) THEN
             IN = INODE
             DO WHILE (IN.GT.0) 
               IJROW = -ITLOC(IN)  
               DO K = K1RHS, K2RHS
                J    = IW(K)       
                I    = ITLOC(J)    
                ILOC = mod(I,NBCOLF) 
              APOS = POSELT+int(ILOC-1,8)*int(NBCOLF,8) + 
     &               int(IJROW-1,8) 
              A(APOS) = A(APOS) + RHS_MUMPS(
     &                 (JFirstRHS+(K-K1RHS)-1)*KEEP(254)+ IN)
             ENDDO
             IN = FILS(IN)
            ENDDO
           ENDIF
          ELSE  
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
            JPOS     = JPOS + 1
           END DO
          ENDIF
          ELBEG  = FRT_PTR(INODE)
          NUMELT = FRT_PTR(INODE+1) - ELBEG
          DO IELL=ELBEG,ELBEG+NUMELT-1
           ELTI = FRT_ELT(IELL)
           J1= PTRAIW(ELTI)
           J2= PTRAIW(ELTI+1)-1
           AII = PTRARW(ELTI)
           SIZE_ELTI = J2 - J1 + 1
           DO II=J1,J2
            I = ITLOC(INTARR(II))
            IF (KEEP(50).EQ.0) THEN
             IF (I.LE.0) CYCLE
             AINPUT    = AII + II - J1
             IPOS = mod(I,NBCOLF)
             ICT12 = POSELT + int(IPOS-1,8) * int(NBCOLF,8)
             DO JJ = J1, J2
              JPOS = ITLOC(INTARR(JJ))
              IF (JPOS.LE.0) THEN 
                   JPOS = -JPOS
              ELSE
                   JPOS = JPOS/NBCOLF
              END IF
              APOS2    = ICT12 + int(JPOS - 1,8)
              A(APOS2) = A(APOS2) +  DBLARR(AINPUT)
              AINPUT   = AINPUT + SIZE_ELTI
             END DO
            ELSE
              IF ( I .EQ. 0 ) THEN 
               AII = AII + J2 - II + 1
               CYCLE
              ENDIF
              IF ( I .LE. 0 ) THEN 
               IPOS1 = -I
               IPOS2 = 0
              ELSE 
               IPOS1 = I/NBCOLF
               IPOS2 = mod(I,NBCOLF)
              END IF
              ICT12 =  POSELT + int(IPOS2-1,8)*int(NBCOLF,8)
              DO JJ=II,J2
               AII = AII + 1
               J = ITLOC(INTARR(JJ))
               IF ( J .EQ. 0 ) CYCLE
               IF ( IPOS2.EQ.0 .AND. J.LE.0) CYCLE
               IF ( J .LE. 0 ) THEN
                JPOS = -J
               ELSE
                JPOS = J/NBCOLF
               END IF
               IF ( (IPOS1.GE.JPOS) .AND. (IPOS2.GT.0) ) THEN
                 APOS2 = ICT12  + int(JPOS - 1,8)
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
               IF ( (IPOS1.LT.JPOS) .AND. (J.GT.0) ) THEN
                 IPOS = mod(J,NBCOLF)
                 JPOS = IPOS1
                 APOS2 = POSELT + int(IPOS-1,8)*int(NBCOLF,8)
     &                          + int(JPOS - 1,8)
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
              END DO
            END IF
           END DO
          END DO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          END DO
      END IF
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          END DO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_123
      SUBROUTINE DMUMPS_126(
     &            N, NELT, NA_ELT,
     &            COMM, MYID, SLAVEF,
     &            IELPTR_LOC, RELPTR_LOC,
     &            ELTVAR_LOC, ELTVAL_LOC,
     &            KEEP,KEEP8, MAXELT_SIZE,
     &            FRTPTR, FRTELT, A, LA, FILS,
     &            id, root )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, NELT, NA_ELT
      INTEGER COMM, MYID, SLAVEF, MAXELT_SIZE, MSGLEN
      INTEGER(8), intent(IN) :: LA
      INTEGER FRTPTR( N+1 )
      INTEGER FRTELT( NELT ), FILS ( N )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IELPTR_LOC( NELT + 1 ), RELPTR_LOC( NELT + 1 )
      INTEGER ELTVAR_LOC( max(1,KEEP(14)) )
      DOUBLE PRECISION ELTVAL_LOC( max(1,KEEP(13)) )
      DOUBLE PRECISION A( LA )
      TYPE(DMUMPS_STRUC)     :: id
      TYPE(DMUMPS_ROOT_STRUC) :: root
      INTEGER numroc
      EXTERNAL numroc
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR_MPI
      INTEGER MSGTAG
      INTEGER allocok
      INTEGER I, DEST, MAXELT_REAL_SIZE, MPG, IEL, SIZEI, SIZER
      INTEGER NBRECORDS, NBUF
      INTEGER RECV_IELTPTR, RECV_RELTPTR
      INTEGER IELTPTR, RELTPTR, INODE
      LOGICAL FINI, PROKG, I_AM_SLAVE
      INTEGER(8) :: PTR_ROOT
      INTEGER LOCAL_M, LOCAL_N, LP, IBEG, IGLOB, JGLOB
      INTEGER ARROW_ROOT
      INTEGER IELT, J, K, NB_REC, IREC
      INTEGER ILOCROOT, JLOCROOT, IPOSROOT, JPOSROOT, IPTR
      INTEGER JCOL_GRID, IROW_GRID
      INTEGER IVALPTR
      INTEGER NBELROOT
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      DOUBLE PRECISION  VAL
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, DIMENSION( :, : ), ALLOCATABLE :: BUFI
      DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE :: BUFR
      DOUBLE PRECISION, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_R
      INTEGER, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_I
      INTEGER, DIMENSION( : ), ALLOCATABLE :: ELROOTPOS
      INTEGER, DIMENSION( : ), ALLOCATABLE, TARGET :: RG2LALLOC
      INTEGER, DIMENSION( : ), POINTER     :: RG2L
      MPG = id%ICNTL(3)
      LP  = id%ICNTL(1)
      I_AM_SLAVE = ( KEEP(46) .eq. 1 .or. MYID .ne.MASTER )
      PROKG = ( MPG > 0 .and. MYID .eq. MASTER )
      KEEP(49) = 0
      ARROW_ROOT = 0
      IF ( MYID .eq. MASTER ) THEN
        IF ( KEEP(46) .eq. 0 ) THEN
          NBUF = SLAVEF
        ELSE
          NBUF = SLAVEF - 1
        END IF
        NBRECORDS = min(KEEP(39),NA_ELT)
        IF ( KEEP(50) .eq. 0 ) THEN
          MAXELT_REAL_SIZE = MAXELT_SIZE * MAXELT_SIZE
        ELSE
          MAXELT_REAL_SIZE = MAXELT_SIZE * (MAXELT_SIZE+1)/2
        END IF
        IF ( MAXELT_REAL_SIZE .GT. KEEP(39) ) THEN
          NBRECORDS = MAXELT_REAL_SIZE
          IF ( MPG .GT. 0 ) THEN
            WRITE(MPG,*)
     & ' ** Warning : For element distrib NBRECORDS set to ',
     & MAXELT_REAL_SIZE,' because one element is large'
          END IF
        END IF
        ALLOCATE( BUFI( 2*NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = 2*NBRECORDS + 1
          GOTO 100
        END IF
        ALLOCATE( BUFR( NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = NBRECORDS + 1
          GOTO 100
        END IF
        IF ( KEEP(52) .ne. 0 ) THEN
          ALLOCATE( TEMP_ELT_R( MAXELT_REAL_SIZE ), stat =allocok )
          IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_REAL_SIZE
            GOTO 100
          END IF
        END IF
        ALLOCATE( TEMP_ELT_I( MAXELT_SIZE ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_SIZE
            GOTO 100
        END IF
        IF ( KEEP(38) .ne. 0 ) THEN
          NBELROOT = FRTPTR(KEEP(38)+1)-FRTPTR(KEEP(38))
          ALLOCATE( ELROOTPOS( max(NBELROOT,1) ),
     &              stat = allocok )
          IF ( allocok .gt. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = NBELROOT
              GOTO 100
          END IF
          IF (KEEP(46) .eq. 0 ) THEN
           ALLOCATE( RG2LALLOC( N ), stat = allocok )
           IF ( allocok .gt. 0 ) THEN
               id%INFO(1) = -13
               id%INFO(2) = N
               GOTO 100
           END IF
           INODE = KEEP(38)
           I     = 1
           DO WHILE ( INODE .GT. 0 )
             RG2LALLOC( INODE ) = I
             INODE = FILS( INODE )
             I = I + 1
           END DO
           RG2L => RG2LALLOC
          ELSE 
           RG2L => root%RG2L_ROW
          END IF
        END IF
        DO I = 1, NBUF
          BUFI( 1, I ) = 0
          BUFR( 1, I ) = ZERO
        END DO
      END IF
 100  CONTINUE
      CALL MUMPS_276( id%ICNTL(1), id%INFO(1), COMM, MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( NBRECORDS, 1, MPI_INTEGER, MASTER,
     &                COMM, IERR_MPI )
      RECV_IELTPTR = 1
      RECV_RELTPTR = 1
      IF ( MYID .eq. MASTER ) THEN
        NBELROOT = 0
        RELTPTR = 1
        RELPTR_LOC(1) = 1
        DO IEL = 1, NELT
          IELTPTR = id%ELTPTR( IEL )
          SIZEI   = id%ELTPTR( IEL + 1 ) - IELTPTR
          IF ( KEEP( 50 ) .eq. 0 ) THEN
            SIZER = SIZEI * SIZEI
          ELSE
            SIZER = SIZEI * ( SIZEI + 1 ) / 2
          END IF
          DEST = id%ELTPROC( IEL )
          IF ( DEST .eq. -2 ) THEN
            NBELROOT = NBELROOT + 1
            FRTELT( FRTPTR(KEEP(38)) + NBELROOT - 1 ) = IEL
            ELROOTPOS( NBELROOT ) = RELTPTR
            GOTO 200
          END IF
          IF ( DEST .ge. 0 .and. KEEP(46) .eq. 0 ) DEST = DEST + 1
          IF ( KEEP(52) .ne. 0 ) THEN
            CALL DMUMPS_288( N, SIZEI, SIZER,
     &               id%ELTVAR( IELTPTR ), id%A_ELT( RELTPTR ),
     &               TEMP_ELT_R(1), MAXELT_REAL_SIZE,
     &               id%ROWSCA(1), id%COLSCA(1), KEEP(50) )
          END IF
          IF ( DEST .eq. 0 .or. ( DEST .eq. -1 .and. KEEP(46) .ne. 0 ) )
     &      THEN
            ELTVAR_LOC( RECV_IELTPTR: RECV_IELTPTR + SIZEI - 1 )
     &      = id%ELTVAR( IELTPTR: IELTPTR + SIZEI - 1 )
            RECV_IELTPTR = RECV_IELTPTR + SIZEI
            IF ( KEEP(52) .ne. 0 ) THEN
              ELTVAL_LOC( RECV_RELTPTR: RECV_RELTPTR + SIZER - 1)
     &        = TEMP_ELT_R( 1: SIZER )
              RECV_RELTPTR = RECV_RELTPTR + SIZER
            END IF
          END IF
          IF ( DEST .NE. 0 .AND. DEST. NE. -3 ) THEN
            IF ( KEEP(52) .eq. 0 ) THEN
              CALL DMUMPS_127(
     &           id%ELTVAR(IELTPTR),
     &           id%A_ELT (RELTPTR),
     &           SIZEI, SIZER,
     &
     &           DEST, NBUF, NBRECORDS,
     &           BUFI, BUFR, COMM )
            ELSE
              CALL DMUMPS_127(
     &           id%ELTVAR(IELTPTR),
     &           TEMP_ELT_R( 1 ),
     &           SIZEI, SIZER,
     &
     &           DEST, NBUF, NBRECORDS,
     &           BUFI, BUFR, COMM )
            END IF
          END IF
 200      CONTINUE
          RELTPTR = RELTPTR + SIZER
          IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
            RELPTR_LOC( IEL + 1 ) = RELTPTR
          ELSE
            RELPTR_LOC( IEL + 1 ) = RECV_RELTPTR
          ENDIF
        END DO
        IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
          KEEP(13) = RELTPTR - 1
        ELSE
          KEEP(13) = RECV_RELTPTR - 1
        ENDIF
        IF ( RELTPTR - 1 .ne. id%NA_ELT ) THEN
          WRITE(*,*) ' ** ERROR ELT DIST: RELPTR - 1 / id%NA_ELT=',
     &               RELTPTR - 1,id%NA_ELT
          CALL MUMPS_ABORT()
        END IF
        DEST = -2
        IELTPTR = 1
        RELTPTR = 1
        SIZEI   = 1
        SIZER   = 1
        CALL DMUMPS_127(
     &           id%ELTVAR(IELTPTR),
     &           id%A_ELT (RELTPTR),
     &           SIZEI, SIZER,
     &
     &           DEST, NBUF, NBRECORDS,
     &           BUFI, BUFR, COMM )
        IF ( KEEP(52) .NE. 0 ) DEALLOCATE( TEMP_ELT_R )
      ELSE
        FINI = ( RECV_IELTPTR .eq. IELPTR_LOC( NELT+1 )
     &     .and. RECV_RELTPTR .eq. RELPTR_LOC( NELT+1 ) )
        DO WHILE ( .not. FINI )
          CALL MPI_PROBE( MASTER, MPI_ANY_TAG,
     &                    COMM, STATUS, IERR_MPI )
          MSGTAG = STATUS( MPI_TAG    )
          SELECT CASE ( MSGTAG )
             CASE( ELT_INT )
               CALL MPI_GET_COUNT( STATUS, MPI_INTEGER,
     &                             MSGLEN, IERR_MPI )
               CALL MPI_RECV( ELTVAR_LOC( RECV_IELTPTR ), MSGLEN,
     &            MPI_INTEGER, MASTER, ELT_INT,
     &            COMM, STATUS, IERR_MPI )
               RECV_IELTPTR = RECV_IELTPTR + MSGLEN
             CASE( ELT_REAL )
                CALL MPI_GET_COUNT( STATUS, MPI_DOUBLE_PRECISION,
     &                              MSGLEN, IERR_MPI )
                CALL MPI_RECV( ELTVAL_LOC( RECV_RELTPTR ), MSGLEN,
     &            MPI_DOUBLE_PRECISION, MASTER, ELT_REAL,
     &            COMM, STATUS, IERR_MPI )
                RECV_RELTPTR = RECV_RELTPTR + MSGLEN
          END SELECT
          FINI = ( RECV_IELTPTR .eq. IELPTR_LOC( NELT+1 )
     &       .and. RECV_RELTPTR .eq. RELPTR_LOC( NELT+1 ) )
        END DO
      END IF
      IF ( KEEP(38) .NE. 0 ) THEN
         IF ( I_AM_SLAVE .and. root%yes ) THEN
          IF (KEEP(60)==0) THEN
           LOCAL_M = numroc( root%ROOT_SIZE, root%MBLOCK,
     &             root%MYROW, 0, root%NPROW )
           LOCAL_M = max( 1, LOCAL_M )
           LOCAL_N = numroc( root%ROOT_SIZE, root%NBLOCK,
     &               root%MYCOL, 0, root%NPCOL )
           PTR_ROOT = LA - int(LOCAL_M,8) * int(LOCAL_N,8) + 1_8
           IF ( PTR_ROOT .LE. LA ) THEN
             A( PTR_ROOT:LA ) = ZERO
           END IF
          ELSE
           DO I = 1, root%SCHUR_NLOC
            root%SCHUR_POINTER((I-1)*root%SCHUR_LLD+1:
     &      (I-1)*root%SCHUR_LLD+root%SCHUR_MLOC)=ZERO
           ENDDO
          ENDIF
         END IF
        IF ( MYID .NE. MASTER ) THEN
          ALLOCATE( BUFI( NBRECORDS * 2 + 1, 1 ), stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS * 2 + 1
            GOTO 250
          END IF
          ALLOCATE( BUFR( NBRECORDS, 1 )        , stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS
          END IF
        END IF
 250    CONTINUE
        CALL MUMPS_276( id%ICNTL(1), id%INFO(1), COMM, MYID )
        IF ( id%INFO(1) .LT. 0 ) RETURN
        IF ( MYID .eq. MASTER ) THEN
        DO IPTR = FRTPTR(KEEP(38)), FRTPTR(KEEP(38)+1) - 1
          IELT = FRTELT( IPTR )
          SIZEI = id%ELTPTR( IELT + 1 ) - id%ELTPTR( IELT )
          DO I = 1, SIZEI
            TEMP_ELT_I( I ) = RG2L
     &              ( id%ELTVAR( id%ELTPTR(IELT) + I - 1 ) )
          END DO
          IVALPTR = ELROOTPOS( IPTR - FRTPTR(KEEP(38)) + 1 ) - 1
          K = 1
          DO J = 1, SIZEI
            JGLOB = id%ELTVAR( id%ELTPTR( IELT ) + J - 1 )
            IF ( KEEP(50).eq. 0 ) THEN
              IBEG = 1
            ELSE
              IBEG = J
            END IF
            DO I = IBEG, SIZEI
              IGLOB = id%ELTVAR( id%ELTPTR( IELT ) + I - 1 )
              IF ( KEEP(52) .eq. 0 ) THEN
                VAL = id%A_ELT( IVALPTR + K )
              ELSE
                VAL = id%A_ELT( IVALPTR + K ) *
     &                id%ROWSCA( IGLOB ) * id%COLSCA( JGLOB )
              END IF
              IF ( KEEP(50).eq.0 ) THEN
                IPOSROOT = TEMP_ELT_I( I )
                JPOSROOT = TEMP_ELT_I( J )
              ELSE
                IF ( TEMP_ELT_I(I) .GT. TEMP_ELT_I(J) ) THEN
                  IPOSROOT = TEMP_ELT_I(I)
                  JPOSROOT = TEMP_ELT_I(J)
                ELSE
                  IPOSROOT = TEMP_ELT_I(J)
                  JPOSROOT = TEMP_ELT_I(I)
                END IF
              END IF
              IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK,
     &                           root%NPROW )
              JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK,
     &                           root%NPCOL )
              IF ( KEEP(46) .eq. 0 ) THEN
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID + 1
              ELSE
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID
              END IF
              IF ( DEST .eq. MASTER ) THEN
                ARROW_ROOT = ARROW_ROOT + 1
                ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
                JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
              IF (KEEP(60)==0) THEN
                A( PTR_ROOT
     &             + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &             + int(ILOCROOT - 1,8) )
     &          =  A( PTR_ROOT
     &             + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &             + int(ILOCROOT - 1,8) )
     &          + VAL
              ELSE
                root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                            * int(root%SCHUR_LLD,8)
     &                            + int(ILOCROOT,8) )
     &          = root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                            * int(root%SCHUR_LLD,8)
     &                            + int(ILOCROOT,8) )
     &          + VAL
              ENDIF
              ELSE
                CALL DMUMPS_34(
     &          IPOSROOT, JPOSROOT, VAL, DEST, BUFI, BUFR, NBRECORDS,
     &          NBUF, LP, COMM, KEEP(46) )
              END IF
              K = K + 1
            END DO
          END DO
        END DO
        CALL DMUMPS_18(
     &          BUFI, BUFR, NBRECORDS,
     &          NBUF, LP, COMM, KEEP(46) )
        ELSE
          FINI = .FALSE.
          DO WHILE ( .not. FINI )
            CALL MPI_RECV( BUFI(1,1), 2*NBRECORDS+1,
     &                MPI_INTEGER, MASTER,
     &                ARROWHEAD,
     &                COMM, STATUS, IERR_MPI )
            NB_REC = BUFI(1,1)
            IF (NB_REC.LE.0) THEN
              FINI = .TRUE.
              NB_REC = -NB_REC
            ENDIF
            IF (NB_REC.EQ.0) EXIT
            CALL MPI_RECV( BUFR(1,1), NBRECORDS, MPI_DOUBLE_PRECISION,
     &                     MASTER, ARROWHEAD,
     &                     COMM, STATUS, IERR_MPI )
            ARROW_ROOT = ARROW_ROOT + NB_REC
            DO IREC = 1, NB_REC
              IPOSROOT = BUFI( IREC * 2, 1 )
              JPOSROOT = BUFI( IREC * 2 + 1, 1 )
              VAL      = BUFR( IREC, 1 )
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
              IF (KEEP(60).eq.0) THEN
                 A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &                       + int(ILOCROOT-1,8))
     &        =  A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &                       + int(ILOCROOT-1,8))
     &           + VAL
              ELSE
                root%SCHUR_POINTER(int(JLOCROOT-1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8) )
     &        = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8))
     &          + VAL
              ENDIF
            END DO
          END DO
          DEALLOCATE( BUFI )
          DEALLOCATE( BUFR )
        END IF
      END IF
      IF ( MYID .eq. MASTER ) THEN
        DEALLOCATE( BUFI )
        DEALLOCATE( BUFR )
        IF (KEEP(38).ne.0) THEN 
          DEALLOCATE(ELROOTPOS)
          IF (KEEP(46) .eq. 0 ) THEN
             DEALLOCATE(RG2LALLOC)
          ENDIF
        ENDIF
        DEALLOCATE( TEMP_ELT_I )
      END IF
      KEEP(49) = ARROW_ROOT
      RETURN
      END SUBROUTINE DMUMPS_126
      SUBROUTINE DMUMPS_127(
     &       ELNODES, ELVAL, SIZEI, SIZER,
     &       DEST, NBUF, NBRECORDS, BUFI, BUFR, COMM )
      IMPLICIT NONE
      INTEGER SIZEI, SIZER, DEST, NBUF, NBRECORDS, COMM
      INTEGER ELNODES( SIZEI ), BUFI( 2*NBRECORDS + 1, NBUF )
      DOUBLE PRECISION ELVAL( SIZER ), BUFR( NBRECORDS + 1, NBUF )
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER I, IBEG, IEND, IERR_MPI, NBRECR
      INTEGER NBRECI
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      IF ( DEST .lt. 0 ) THEN
        IBEG = 1
        IEND = NBUF
      ELSE
        IBEG = DEST
        IEND = DEST
      END IF
      DO I = IBEG, IEND
        NBRECI = BUFI(1,I)
        IF ( NBRECI .ne.0  .and.
     &       ( DEST.eq.-2 .or.
     &         NBRECI + SIZEI .GT. 2*NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFI(2, I), NBRECI, MPI_INTEGER,
     &                    I, ELT_INT, COMM, IERR_MPI )
           BUFI(1,I) = 0
           NBRECI    = 0
        END IF
        NBRECR = int(dble(BUFR(1,I))+0.5D0)
        IF ( NBRECR .ne.0  .and.
     &       ( DEST.eq.-2 .or.
     &         NBRECR + SIZER .GT. NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFR(2, I), NBRECR, MPI_DOUBLE_PRECISION,
     &                    I, ELT_REAL, COMM, IERR_MPI )
           BUFR(1,I) = ZERO
           NBRECR    = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          BUFI( 2 + NBRECI : 2 + NBRECI + SIZEI - 1, I ) =
     &    ELNODES( 1: SIZEI )
          BUFR( 2 + NBRECR : 2 + NBRECR + SIZER - 1, I ) =
     &    ELVAL( 1: SIZER )
          BUFI(1,I) = NBRECI + SIZEI
          BUFR(1,I) = dble( NBRECR + SIZER )
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_127
      SUBROUTINE DMUMPS_213( ELTPTR, NELT, MAXELT_SIZE )
      INTEGER NELT, MAXELT_SIZE
      INTEGER ELTPTR( NELT + 1 )
      INTEGER I, S
      MAXELT_SIZE = 0
      DO I = 1, NELT
        S = ELTPTR( I + 1 ) - ELTPTR( I )
        MAXELT_SIZE = max( S, MAXELT_SIZE )
      END DO
      RETURN
      END SUBROUTINE DMUMPS_213
      SUBROUTINE DMUMPS_288( N, SIZEI, SIZER,
     &               ELTVAR, ELTVAL,
     &               SELTVAL, LSELTVAL,
     &               ROWSCA, COLSCA, K50 )
      INTEGER N, SIZEI, SIZER, LSELTVAL, K50
      INTEGER ELTVAR( SIZEI )
      DOUBLE PRECISION ELTVAL( SIZER )
      DOUBLE PRECISION SELTVAL( LSELTVAL )
      DOUBLE PRECISION ROWSCA( N ), COLSCA( N )
      INTEGER I, J, K
      K = 1
      IF ( K50 .eq. 0 ) THEN
        DO J = 1, SIZEI
          DO I = 1, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     &                   ROWSCA(ELTVAR(I)) *
     &                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      ELSE
        DO J = 1, SIZEI
          DO I = J, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     &                   ROWSCA(ELTVAR(I)) *
     &                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_288
      SUBROUTINE DMUMPS_F77( JOB, SYM, PAR, COMM_F77, N, ICNTL, CNTL,
     &                      NZ, IRN, IRNhere, JCN, JCNhere, A, Ahere,
     &                      NZ_loc, IRN_loc, IRN_lochere,
     &                      JCN_loc, JCN_lochere,
     &                      A_loc, A_lochere,
     &                      NELT, ELTPTR, ELTPTRhere,  ELTVAR,
     &                      ELTVARhere, A_ELT, A_ELThere,
     &                      PERM_IN, PERM_INhere,
     &                      RHS, RHShere, REDRHS, REDRHShere,
     &                      INFO, RINFO, INFOG, RINFOG,
     &                      DEFICIENCY, LWK_USER,
     &                      SIZE_SCHUR, LISTVAR_SCHUR,
     &                      LISTVAR_SCHURhere, SCHUR, SCHURhere,
     &                      WK_USER, WK_USERhere,
     &                      COLSCA, COLSCAhere, ROWSCA, ROWSCAhere,
     &                      INSTANCE_NUMBER, NRHS, LRHS, LREDRHS,
     &
     &                      RHS_SPARSE, RHS_SPARSEhere,
     &                      SOL_loc, SOL_lochere,
     &                      IRHS_SPARSE, IRHS_SPARSEhere,
     &                      IRHS_PTR, IRHS_PTRhere,
     &                      ISOL_loc, ISOL_lochere,
     &                      NZ_RHS, LSOL_loc
     &                      , 
     & SCHUR_MLOC,
     & SCHUR_NLOC,
     & SCHUR_LLD,
     & MBLOCK,
     & NBLOCK,
     & NPROW,
     & NPCOL,
     &
     & OOC_TMPDIR,
     & OOC_PREFIX,
     & WRITE_PROBLEM,
     & TMPDIRLEN,
     & PREFIXLEN,
     & WRITE_PROBLEMLEN
     &
     & )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER OOC_PREFIX_MAX_LENGTH, OOC_TMPDIR_MAX_LENGTH
      INTEGER PB_MAX_LENGTH
      PARAMETER(OOC_PREFIX_MAX_LENGTH=63, OOC_TMPDIR_MAX_LENGTH=255)
      PARAMETER(PB_MAX_LENGTH=255)
      INTEGER JOB, SYM, PAR, COMM_F77, N, NZ, NZ_loc, NELT,
     &        DEFICIENCY, LWK_USER, SIZE_SCHUR, INSTANCE_NUMBER,
     &        NRHS, LRHS,
     &        NZ_RHS, LSOL_loc, LREDRHS
      INTEGER ICNTL(40), INFO(40), INFOG(40)
      INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
      INTEGER MBLOCK, NBLOCK, NPROW, NPCOL
      INTEGER TMPDIRLEN, PREFIXLEN, WRITE_PROBLEMLEN
      DOUBLE PRECISION CNTL(15), RINFO(40), RINFOG(40)
      INTEGER, TARGET :: IRN(*), JCN(*), ELTPTR(*), ELTVAR(*)
      INTEGER, TARGET :: PERM_IN(*), IRN_loc(*), JCN_loc(*)
      INTEGER, TARGET :: LISTVAR_SCHUR(*)
      INTEGER, TARGET :: IRHS_PTR(*), IRHS_SPARSE(*), ISOL_loc(*)
      DOUBLE PRECISION, TARGET :: A(*), A_ELT(*), A_loc(*), RHS(*)
      DOUBLE PRECISION, TARGET :: WK_USER(*)
      DOUBLE PRECISION, TARGET :: REDRHS(*)
      DOUBLE PRECISION, TARGET :: ROWSCA(*), COLSCA(*)
      DOUBLE PRECISION, TARGET :: SCHUR(*)
      DOUBLE PRECISION, TARGET :: RHS_SPARSE(*), SOL_loc(*)
      INTEGER, INTENT(in) :: OOC_TMPDIR(OOC_TMPDIR_MAX_LENGTH)
      INTEGER, INTENT(in) :: OOC_PREFIX(OOC_PREFIX_MAX_LENGTH)
      INTEGER, INTENT(in) :: WRITE_PROBLEM(PB_MAX_LENGTH)
      INTEGER IRNhere, JCNhere, Ahere, ELTPTRhere, ELTVARhere,
     &        A_ELThere, PERM_INhere, WK_USERhere,
     &        RHShere, REDRHShere, IRN_lochere,
     &        JCN_lochere, A_lochere, LISTVAR_SCHURhere,
     &        SCHURhere, COLSCAhere, ROWSCAhere, RHS_SPARSEhere,
     &        SOL_lochere, IRHS_PTRhere, IRHS_SPARSEhere, ISOL_lochere
      INCLUDE 'mpif.h'
      TYPE DMUMPS_STRUC_PTR
          TYPE (DMUMPS_STRUC), POINTER :: PTR
      END TYPE DMUMPS_STRUC_PTR
      TYPE (DMUMPS_STRUC), POINTER :: mumps_par
      TYPE (DMUMPS_STRUC_PTR), DIMENSION (:), POINTER, SAVE ::
     &  mumps_par_array
      TYPE (DMUMPS_STRUC_PTR), DIMENSION (:), POINTER ::
     &  mumps_par_array_bis
      INTEGER, SAVE :: DMUMPS_STRUC_ARRAY_SIZE = 0
      INTEGER, SAVE :: N_INSTANCES = 0
      INTEGER A_ELT_SIZE, I, Np, IERR
      INTEGER DMUMPS_STRUC_ARRAY_SIZE_INIT
      PARAMETER (DMUMPS_STRUC_ARRAY_SIZE_INIT=10)
      EXTERNAL MUMPS_AFFECT_MAPPING,
     &         MUMPS_AFFECT_PIVNUL_LIST,
     &         MUMPS_AFFECT_SYM_PERM,
     &         MUMPS_AFFECT_UNS_PERM
      IF (JOB == -1) THEN
        DO I = 1, DMUMPS_STRUC_ARRAY_SIZE
          IF ( .NOT. associated(mumps_par_array(I)%PTR) ) GOTO 10
        END DO
        ALLOCATE( mumps_par_array_bis(DMUMPS_STRUC_ARRAY_SIZE +
     &  DMUMPS_STRUC_ARRAY_SIZE_INIT), stat=IERR)
        IF (IERR /= 0) THEN
          WRITE(*,*) ' ** Allocation Error 1 in DMUMPS_F77.'
          CALL MUMPS_ABORT()
        END IF
        DO I = 1, DMUMPS_STRUC_ARRAY_SIZE
          mumps_par_array_bis(I)%PTR=>mumps_par_array(I)%PTR
        ENDDO
        IF (associated(mumps_par_array)) DEALLOCATE(mumps_par_array)
        mumps_par_array=>mumps_par_array_bis
        NULLIFY(mumps_par_array_bis)
        DO I = DMUMPS_STRUC_ARRAY_SIZE+1, DMUMPS_STRUC_ARRAY_SIZE +
     &  DMUMPS_STRUC_ARRAY_SIZE_INIT
          NULLIFY(mumps_par_array(I)%PTR)
        ENDDO
        I = DMUMPS_STRUC_ARRAY_SIZE+1
        DMUMPS_STRUC_ARRAY_SIZE = DMUMPS_STRUC_ARRAY_SIZE +
     &  DMUMPS_STRUC_ARRAY_SIZE_INIT
 10     CONTINUE
        INSTANCE_NUMBER = I
        N_INSTANCES = N_INSTANCES+1
        ALLOCATE( mumps_par_array(INSTANCE_NUMBER)%PTR,stat=IERR )
        IF (IERR /= 0) THEN
          WRITE(*,*) '** Allocation Error 2 in DMUMPS_F77.'
          CALL MUMPS_ABORT()
        ENDIF
        mumps_par_array(INSTANCE_NUMBER)%PTR%KEEP(40) = 0
        mumps_par_array(INSTANCE_NUMBER)%PTR%INSTANCE_NUMBER =
     &  INSTANCE_NUMBER
      END IF
      IF ( INSTANCE_NUMBER .LE. 0 .OR. INSTANCE_NUMBER .GT.
     &     DMUMPS_STRUC_ARRAY_SIZE ) THEN
        WRITE(*,*) ' ** Instance Error 1 in DMUMPS_F77',
     &             INSTANCE_NUMBER
        CALL MUMPS_ABORT()
      END IF
      IF ( .NOT. associated ( mumps_par_array(INSTANCE_NUMBER)%PTR ) )
     &  THEN
        WRITE(*,*) ' Instance Error 2 in DMUMPS_F77',
     &             INSTANCE_NUMBER
        CALL MUMPS_ABORT()
      END IF
      mumps_par => mumps_par_array(INSTANCE_NUMBER)%PTR
      mumps_par%SYM = SYM
      mumps_par%PAR = PAR
      mumps_par%JOB = JOB
      mumps_par%N   = N
      mumps_par%NZ  = NZ
      mumps_par%NZ_loc  = NZ_loc
      mumps_par%LWK_USER = LWK_USER
      mumps_par%SIZE_SCHUR  = SIZE_SCHUR
      mumps_par%NELT= NELT
      mumps_par%ICNTL(1:40)=ICNTL(1:40)
      mumps_par%CNTL(1:15)=CNTL(1:15)
      mumps_par%NRHS  = NRHS
      mumps_par%LRHS  = LRHS
      mumps_par%LREDRHS = LREDRHS
      mumps_par%NZ_RHS   = NZ_RHS
      mumps_par%LSOL_loc = LSOL_loc
      mumps_par%SCHUR_MLOC   = SCHUR_MLOC
      mumps_par%SCHUR_NLOC   = SCHUR_NLOC
      mumps_par%SCHUR_LLD    = SCHUR_LLD
      mumps_par%MBLOCK = MBLOCK
      mumps_par%NBLOCK = NBLOCK
      mumps_par%NPROW  = NPROW
      mumps_par%NPCOL  = NPCOL
      IF ( COMM_F77 .NE. -987654 ) THEN
        mumps_par%COMM = COMM_F77
      ELSE
        mumps_par%COMM = MPI_COMM_WORLD
      ENDIF
      CALL MPI_BCAST(NRHS,1,MPI_INTEGER,0,mumps_par%COMM,IERR)
      IF ( IRNhere /= 0 ) mumps_par%IRN => IRN(1:NZ)
      IF ( JCNhere /= 0 ) mumps_par%JCN => JCN(1:NZ)
      IF ( Ahere /= 0 )   mumps_par%A   => A(1:NZ)
      IF ( IRN_lochere /= 0 ) mumps_par%IRN_loc => IRN_loc(1:NZ_loc)
      IF ( JCN_lochere /= 0 ) mumps_par%JCN_loc => JCN_loc(1:NZ_loc)
      IF ( A_lochere /= 0 )   mumps_par%A_loc   => A_loc(1:NZ_loc)
      IF ( ELTPTRhere /= 0 ) mumps_par%ELTPTR => ELTPTR(1:NELT+1)
      IF ( ELTVARhere /= 0 ) mumps_par%ELTVAR =>
     &   ELTVAR(1:ELTPTR(NELT+1)-1)
      IF ( A_ELThere /= 0 ) THEN
        A_ELT_SIZE = 0
        DO I = 1, NELT
          Np = ELTPTR(I+1) -ELTPTR(I)
          IF (SYM == 0) THEN
            A_ELT_SIZE = A_ELT_SIZE + Np * Np
          ELSE
            A_ELT_SIZE = A_ELT_SIZE + Np * ( Np + 1 ) / 2
          END IF
        END DO
        mumps_par%A_ELT => A_ELT(1:A_ELT_SIZE)
      END IF
      IF ( PERM_INhere /= 0) mumps_par%PERM_IN => PERM_IN(1:N)
      IF ( LISTVAR_SCHURhere /= 0)
     &   mumps_par%LISTVAR_SCHUR =>LISTVAR_SCHUR(1:SIZE_SCHUR)
      IF ( SCHURhere /= 0 ) THEN
        mumps_par%SCHUR_CINTERFACE=>SCHUR(1:1)
      ENDIF
      IF (NRHS .NE. 1) THEN
        IF ( RHShere /= 0 ) mumps_par%RHS => RHS(1:NRHS*LRHS)
        IF (REDRHShere /= 0)mumps_par%REDRHS=>REDRHS(1:NRHS*LREDRHS)
      ELSE
        IF ( RHShere /= 0 ) mumps_par%RHS => RHS(1:N)
        IF (REDRHShere /= 0)mumps_par%REDRHS=>REDRHS(1:SIZE_SCHUR)
      ENDIF
      IF ( WK_USERhere /=0 ) THEN
        IF (LWK_USER > 0 ) THEN
          mumps_par%WK_USER => WK_USER(1:LWK_USER)
        ELSE
          mumps_par%WK_USER => WK_USER(1_8:-int(LWK_USER,8)*1000000_8)
        ENDIF
      ENDIF
      IF ( COLSCAhere /= 0) mumps_par%COLSCA => COLSCA(1:N)
      IF ( ROWSCAhere /= 0) mumps_par%ROWSCA => ROWSCA(1:N)
      IF ( RHS_SPARSEhere /=0 ) mumps_par%RHS_SPARSE=>
     &                          RHS_SPARSE(1:NZ_RHS)
      IF ( IRHS_SPARSEhere /=0 ) mumps_par%IRHS_SPARSE=>
     &                          IRHS_SPARSE(1:NZ_RHS)
      IF ( SOL_lochere /=0 ) mumps_par%SOL_loc=>
     &                          SOL_loc(1:LSOL_loc*NRHS)
      IF ( ISOL_lochere /=0 ) mumps_par%ISOL_loc=>
     &                          ISOL_loc(1:LSOL_loc)
      IF ( IRHS_PTRhere /=0 ) mumps_par%IRHS_PTR=>
     &                          IRHS_PTR(1:NRHS+1)
      DO I=1,TMPDIRLEN
        mumps_par%OOC_TMPDIR(I:I)=char(OOC_TMPDIR(I))
      ENDDO
      DO I=TMPDIRLEN+1,OOC_TMPDIR_MAX_LENGTH
        mumps_par%OOC_TMPDIR(I:I)=' '
      ENDDO
      DO I=1,PREFIXLEN
        mumps_par%OOC_PREFIX(I:I)=char(OOC_PREFIX(I))
      ENDDO
      DO I=PREFIXLEN+1,OOC_PREFIX_MAX_LENGTH
        mumps_par%OOC_PREFIX(I:I)=' '
      ENDDO
      DO I=1,WRITE_PROBLEMLEN
        mumps_par%WRITE_PROBLEM(I:I)=char(WRITE_PROBLEM(I))
      ENDDO
      DO I=WRITE_PROBLEMLEN+1,PB_MAX_LENGTH
        mumps_par%WRITE_PROBLEM(I:I)=' '
      ENDDO
      CALL DMUMPS( mumps_par )
      INFO(1:40)=mumps_par%INFO(1:40)
      INFOG(1:40)=mumps_par%INFOG(1:40)
      RINFO(1:40)=mumps_par%RINFO(1:40)
      RINFOG(1:40)=mumps_par%RINFOG(1:40)
      ICNTL(1:40) = mumps_par%ICNTL(1:40)
      CNTL(1:15) = mumps_par%CNTL(1:15)
      SYM = mumps_par%SYM
      PAR = mumps_par%PAR
      JOB = mumps_par%JOB
      N   = mumps_par%N
      NZ  = mumps_par%NZ
      NRHS = mumps_par%NRHS
      LRHS = mumps_par%LRHS
      LREDRHS = mumps_par%LREDRHS
      NZ_loc  = mumps_par%NZ_loc
      NZ_RHS  = mumps_par%NZ_RHS
      LSOL_loc= mumps_par%LSOL_loc
      SIZE_SCHUR  = mumps_par%SIZE_SCHUR
      LWK_USER = mumps_par%LWK_USER
      NELT= mumps_par%NELT
      DEFICIENCY = mumps_par%Deficiency
      SCHUR_MLOC   = mumps_par%SCHUR_MLOC
      SCHUR_NLOC   = mumps_par%SCHUR_NLOC
      SCHUR_LLD    = mumps_par%SCHUR_LLD
      MBLOCK       = mumps_par%MBLOCK
      NBLOCK       = mumps_par%NBLOCK
      NPROW        = mumps_par%NPROW
      NPCOL        = mumps_par%NPCOL
      IF ( associated (mumps_par%MAPPING) ) THEN
         CALL MUMPS_AFFECT_MAPPING(mumps_par%MAPPING(1))
      ELSE
         CALL MUMPS_NULLIFY_C_MAPPING()
      ENDIF
      IF ( associated (mumps_par%PIVNUL_LIST) ) THEN
         CALL MUMPS_AFFECT_PIVNUL_LIST(mumps_par%PIVNUL_LIST(1))
      ELSE
         CALL MUMPS_NULLIFY_C_PIVNUL_LIST()
      ENDIF
      IF ( associated (mumps_par%SYM_PERM) ) THEN
         CALL MUMPS_AFFECT_SYM_PERM(mumps_par%SYM_PERM(1))
      ELSE
         CALL MUMPS_NULLIFY_C_SYM_PERM()
      ENDIF
      IF ( associated (mumps_par%UNS_PERM) ) THEN
         CALL MUMPS_AFFECT_UNS_PERM(mumps_par%UNS_PERM(1))
      ELSE
         CALL MUMPS_NULLIFY_C_UNS_PERM()
      ENDIF
      IF ( JOB == -2 ) THEN
         IF (associated(mumps_par_array(INSTANCE_NUMBER)%PTR))THEN
           DEALLOCATE(mumps_par_array(INSTANCE_NUMBER)%PTR)
           NULLIFY   (mumps_par_array(INSTANCE_NUMBER)%PTR)
           N_INSTANCES = N_INSTANCES - 1
           IF ( N_INSTANCES == 0 ) THEN
             DEALLOCATE(mumps_par_array)
             DMUMPS_STRUC_ARRAY_SIZE = 0
           END IF
         ELSE
           WRITE(*,*) "** Warning: instance already freed"
           WRITE(*,*) "            this should normally not happen."
         ENDIF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_F77
