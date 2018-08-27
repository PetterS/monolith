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
      SUBROUTINE SMUMPS_301( id)
      USE SMUMPS_STRUC_DEF
      USE MUMPS_SOL_ES
      USE SMUMPS_COMM_BUFFER
      USE SMUMPS_OOC
      USE TOOLS_COMMON
      IMPLICIT NONE
      INTERFACE
      SUBROUTINE SMUMPS_710( id, NB_INT,NB_CMPLX )
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC) :: id
      INTEGER(8)        :: NB_INT,NB_CMPLX
      END SUBROUTINE SMUMPS_710
      SUBROUTINE SMUMPS_758 
     &(idRHS, idINFO, idN, idNRHS, idLRHS)
      REAL, DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      END SUBROUTINE SMUMPS_758
      END INTERFACE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_headers.h'
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      TYPE (SMUMPS_STRUC), TARGET :: id
      INTEGER MP,LP, MPG
      LOGICAL PROK, PROKG
      INTEGER MTYPE, ICNTL21
      LOGICAL LSCAL, ERANAL, GIVSOL
      INTEGER ICNTL10, ICNTL11
      INTEGER I,K,JPERM, J, II, IZ2
      INTEGER IZ, NZ_THIS_BLOCK
      INTEGER LIW
      INTEGER(8) :: LA, LA_PASSED
      INTEGER LIW_PASSED
      INTEGER LWCB_MIN, LWCB, LWCB_SOL_C
      INTEGER(8) :: TMP_LWCB8
      INTEGER SMUMPS_LBUF, SMUMPS_LBUF_INT
      INTEGER MSG_MAX_BYTES_SOLVE, MSG_MAX_BYTES_GTHRSOL
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF,
     &        IBEG_GLOB_DEF, IEND_GLOB_DEF,
     &        IROOT_DEF_RHS_COL1
      INTEGER NITREF, NOITER, SOLVET, KASE, JOBIREF
      REAL RSOL(1)
      LOGICAL INTERLEAVE_PAR, DO_PERMUTE_RHS
      INTEGER     :: NRHS_NONEMPTY
      INTEGER     :: STRAT_PERMAM1
      INTEGER     :: K220(0:id%NSLAVES)
      LOGICAL     :: DO_NULL_PIV
      INTEGER, DIMENSION(:), POINTER :: IRHS_PTR_COPY
      INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE_COPY
      REAL, DIMENSION(:), POINTER :: RHS_SPARSE_COPY
      LOGICAL IRHS_SPARSE_COPY_ALLOCATED, IRHS_PTR_COPY_ALLOCATED,
     &        RHS_SPARSE_COPY_ALLOCATED
      INTEGER         :: NBCOL, COLSIZE, JBEG_RHS, JEND_RHS, JBEG_NEW,
     &                   NBCOL_INBLOC, IPOS, NBT
      INTEGER     :: PERM_PIV_LIST(max(id%KEEP(112),1))
      INTEGER     :: MAP_PIVNUL_LIST(max(id%KEEP(112),1))
      INTEGER, DIMENSION(:), ALLOCATABLE :: PERM_RHS
      REAL ONE
      REAL ZERO
      PARAMETER( ONE = 1.0E0 )
      PARAMETER( ZERO = 0.0E0 )
      REAL RZERO, RONE
      PARAMETER( RZERO = 0.0E0, RONE = 1.0E0 )
      REAL, DIMENSION(:), POINTER :: RHS_MUMPS
      REAL, DIMENSION(:), POINTER :: WORK_WCB
      REAL, DIMENSION(:), POINTER :: PTR_RHS_ROOT
      INTEGER :: LPTR_RHS_ROOT
      REAL, ALLOCATABLE :: SAVERHS(:), C_RW1(:),
     &                                 C_RW2(:),
     &                                 SRW3(:), C_Y(:),
     &                                 C_W(:)
      REAL, ALLOCATABLE :: CWORK(:)
      REAL, ALLOCATABLE :: R_RW1(:), R_Y(:), D(:)
      REAL, ALLOCATABLE :: R_W(:)
      REAL,    ALLOCATABLE, DIMENSION(:) :: R_LOCWK54
      REAL, ALLOCATABLE, DIMENSION(:) :: C_LOCWK54
      INTEGER, DIMENSION(:), ALLOCATABLE :: UNS_PERM_INV, 
     &                                      POSINRHSCOMP_N
      INTEGER LIWK_SOLVE, LIWCB
      INTEGER, ALLOCATABLE :: IW1(:), IWK_SOLVE(:), IWCB(:)
      INTEGER(8)       :: MAXS
      REAL, DIMENSION(:), POINTER :: CNTL
      INTEGER, DIMENSION (:), POINTER :: KEEP,ICNTL,INFO
      INTEGER(8), DIMENSION (:), POINTER :: KEEP8
      INTEGER, DIMENSION (:), POINTER :: IS
      REAL, DIMENSION(:),POINTER::   RINFOG
          type scaling_data_t
            SEQUENCE
            REAL, dimension(:), pointer :: SCALING
            REAL, dimension(:), pointer :: SCALING_LOC
          end type scaling_data_t
          type (scaling_data_t) :: scaling_data
          REAL, DIMENSION(:), POINTER :: PT_SCALING
          REAL, TARGET                :: Dummy_SCAL(1)
      REAL ARRET
      REAL C_DUMMY(1)
      REAL R_DUMMY(1)
      INTEGER IDUMMY(1), JDUMMY(1), KDUMMY(1), LDUMMY(1), MDUMMY(1)
      INTEGER, TARGET :: IDUMMY_TARGET(1)
      REAL, TARGET :: CDUMMY_TARGET(1)
      INTEGER JJ, WHAT
      INTEGER allocok
      INTEGER NBRHS, NBRHS_EFF, BEG_RHS, NB_RHSSKIPPED,
     &        IBEG, LD_RHS, KDEC, 
     &        MASTER_ROOT, MASTER_ROOT_IN_COMM
      INTEGER IPT_RHS_ROOT, SIZE_ROOT, LD_REDRHS
      INTEGER IBEG_REDRHS, IBEG_RHSCOMP, LD_RHSCOMP, LENRHSCOMP
      INTEGER NB_K133, IRANK, TSIZE
      INTEGER KMAX_246_247
      LOGICAL WORKSPACE_MINIMAL_PREFERRED, WK_USER_PROVIDED
      INTEGER(8)  NB_BYTES     
      INTEGER(8)  NB_BYTES_MAX 
      INTEGER(8)  NB_BYTES_EXTRA 
      INTEGER(8) NB_INT, NB_CMPLX, K34_8, K35_8, NB_BYTES_ON_ENTRY
      INTEGER(8) K16_8, ITMP8
#if defined(V_T)
      INTEGER soln_drive_class, glob_comm_ini, perm_scal_ini, soln_dist,
     &        soln_assem, perm_scal_post
#endif
      LOGICAL I_AM_SLAVE, BUILD_POSINRHSCOMP
      LOGICAL WORK_WCB_ALLOCATED, IS_INIT_OOC_DONE 
      LOGICAL STOP_AT_NEXT_EMPTY_COL
      INTEGER  MTYPE_LOC
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln driver',soln_drive_class,IERR)
      CALL VTFUNCDEF( 'glob_comm_ini',soln_drive_class,
     &     glob_comm_ini,IERR)
      CALL VTFUNCDEF( 'perm_scal_ini',soln_drive_class,
     &     perm_scal_ini,IERR)
      CALL VTFUNCDEF( 'soln_dist',soln_drive_class,soln_dist,IERR)
      CALL VTFUNCDEF( 'soln_assem',soln_drive_class,soln_assem,IERR)
      CALL VTFUNCDEF( 'perm_scal_post',soln_drive_class,
     &     perm_scal_post,IERR)
#endif
      IRHS_PTR_COPY => IDUMMY_TARGET
      IRHS_PTR_COPY_ALLOCATED = .FALSE.
      IRHS_SPARSE_COPY => IDUMMY_TARGET
      IRHS_SPARSE_COPY_ALLOCATED=.FALSE.
      RHS_SPARSE_COPY => CDUMMY_TARGET
      RHS_SPARSE_COPY_ALLOCATED=.FALSE.
      NULLIFY(RHS_MUMPS)
      NULLIFY(WORK_WCB)
      IS_INIT_OOC_DONE   = .FALSE.
      WK_USER_PROVIDED   = .FALSE.
      WORK_WCB_ALLOCATED = .FALSE.
      CNTL =>id%CNTL
      KEEP =>id%KEEP
      KEEP8=>id%KEEP8
      IS   =>id%IS
      ICNTL=>id%ICNTL
      INFO =>id%INFO
      RINFOG =>id%RINFOG
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LP  = id%ICNTL( 1 )
      PROK  = (MP.GT.0)
      PROKG = (MPG.GT.0 .and. id%MYID.eq.MASTER)
      IF ( PROK  ) WRITE(MP,100)
      IF ( PROKG ) WRITE(MPG,100)
      NB_BYTES       = 0_8
      NB_BYTES_MAX   = 0_8
      NB_BYTES_EXTRA = 0_8
      K34_8    = int(KEEP(34), 8)
      K35_8    = int(KEEP(35), 8)
      K16_8    = int(KEEP(16), 8)
      NB_RHSSKIPPED = 0 
      LSCAL              = .FALSE.
      WORK_WCB_ALLOCATED = .FALSE.
      ICNTL21  = -99998  
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
      CALL SMUMPS_710 (id, NB_INT,NB_CMPLX  )
      NB_BYTES = NB_BYTES + NB_INT * K34_8 + NB_CMPLX * K35_8
      NB_BYTES_ON_ENTRY = NB_BYTES  
      NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
      IF (id%MYID .EQ. MASTER) THEN
          CALL SMUMPS_807(id) 
          id%KEEP(111) = id%ICNTL(25)
          id%KEEP(248) = id%ICNTL(20)
          ICNTL21      = id%ICNTL(21)
          IF (ICNTL21 .ne.0.and.ICNTL21.ne.1) ICNTL21=0
          IF ( id%ICNTL(30) .NE.0 ) THEN 
            id%KEEP(237) = 1
          ELSE
            id%KEEP(237) = 0
          ENDIF
          IF (id%KEEP(248) .eq.0.and. id%KEEP(237).ne.0) THEN 
             id%KEEP(248)=1
          ENDIF
          IF (id%KEEP(248) .ne.0.and.id%KEEP(248).ne.1) id%KEEP(248)=0
          IF ((id%KEEP(221).EQ.2 ).AND.(id%KEEP(248).NE.0) ) THEN
           id%KEEP(248) = 0
          ENDIF
          IF ((id%KEEP(221).EQ.2 ).AND.(id%KEEP(235).NE.0) ) THEN
           id%KEEP(235) = 0
          ENDIF
          IF ( (id%KEEP(248).EQ.0).AND.(id%KEEP(111).EQ.0) ) THEN
            id%KEEP(235) = 0
          ENDIF
          MTYPE = ICNTL(  9 )
          IF (id%KEEP(237).NE.0) MTYPE = 1 
      ENDIF
      CALL MPI_BCAST(MTYPE,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST( id%KEEP(111), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(248), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( ICNTL21, 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(235), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%KEEP(221), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      CALL MPI_BCAST( id%NRHS,1, MPI_INTEGER, MASTER, id%COMM,IERR)
      CALL MPI_BCAST( id%KEEP(237), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
      IF ( id%MYID .EQ. MASTER ) THEN
          IF (KEEP(201) .EQ. -1) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Solve impossible because factors not kept'
             id%INFO(1)=-44
             id%INFO(2)=KEEP(251)
             GOTO 333
          ELSE IF (KEEP(221).EQ.0 .AND. KEEP(251) .EQ. 2
     &            .AND. KEEP(252).EQ.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Solve impossible because factors not kept'
             id%INFO(1)=-44
             id%INFO(2)=KEEP(251)
             GOTO 333
          ENDIF
          IF (KEEP(252).NE.0 .AND. id%NRHS .NE. id%KEEP(253)) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: id%NRHS not allowed to change when ICNTL(32)=1'
             id%INFO(1)=-42
             id%INFO(2)=id%KEEP(253)
             GOTO 333
          ENDIF
          IF (KEEP(252).NE.0 .AND. ICNTL(9).NE.1) THEN
             INFO(1) = -43
             INFO(2) = 9
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Transpose system (ICNTL(9).NE.0) not ',
     &       ' compatible with forward performed during',
     &       ' factorization (ICNTL(32)=1)'
             GOTO 333
          ENDIF
          IF (KEEP(248) .NE. 0.AND.KEEP(252).NE.0) THEN
             INFO(1) = -43
             IF (KEEP(237).NE.0) THEN
              INFO(2) = 30 
              IF (PROKG) WRITE(MPG,'(A)')
     &        ' ERROR: A-1 functionality incompatible with forward',
     &        ' performed during factorization (ICNTL(32)=1)'
             ELSE
              INFO(2) = 20 
              IF (PROKG) WRITE(MPG,'(A)')
     &        ' ERROR: sparse RHS incompatible with forward',
     &        ' performed during factorization (ICNTL(32)=1)'
             ENDIF
             GOTO 333
          ENDIF
          IF (KEEP(237) .NE. 0 .AND. ICNTL21.NE.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality is incompatible',
     &       ' with distributed solution.'
             INFO(1)=-48
             INFO(2)=21
             GOTO 333
          ENDIF
          IF (KEEP(237) .NE. 0 .AND. KEEP(60) .NE.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality is incompatible',
     &       ' with Schur.'
             INFO(1)=-48
             INFO(2)=19
             GOTO 333
          ENDIF
          IF (KEEP(237) .NE. 0 .AND. KEEP(111) .NE.0) THEN
             IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality is incompatible',
     &       ' with null space.'
             INFO(1)=-48
             INFO(2)=25
             GOTO 333
          ENDIF
          IF (id%NRHS .LE. 0) THEN
             id%INFO(1)=-45
             id%INFO(2)=id%NRHS
             GOTO 333
          ENDIF
          IF ( (id%KEEP(237).EQ.0) ) THEN
             IF ((id%KEEP(248) == 0 .AND.KEEP(221).NE.2) 
     &            .OR. ICNTL21==0) THEN
               CALL SMUMPS_758
     &         (id%RHS,id%INFO,id%N,id%NRHS,id%LRHS)
               IF (id%INFO(1) .LT. 0) GOTO 333
             ENDIF
          ELSE
            IF (id%NRHS .NE. id%N) THEN
              id%INFO(1)=-47
              id%INFO(2)=id%NRHS
              GOTO 333
            ENDIF
          ENDIF
          IF (id%KEEP(248) == 1) THEN
            IF ((id%NZ_RHS.LE.0).AND.(KEEP(237).NE.0)) THEN
              id%INFO(1)=-46
              id%INFO(2)=id%NZ_RHS
              GOTO 333
            ENDIF
            IF (( id%NZ_RHS .LE.0 ).AND.(KEEP(221).EQ.1)) THEN
              id%INFO(1)=-46
              id%INFO(2)=id%NZ_RHS
              GOTO 333
            ENDIF
            IF ( .not. associated(id%RHS_SPARSE) )THEN
              id%INFO(1)=-22
              id%INFO(2)=10
              GOTO 333
            ENDIF
            IF ( .not. associated(id%IRHS_SPARSE) )THEN
              id%INFO(1)=-22
              id%INFO(2)=11
              GOTO 333
            ENDIF
            IF ( .not. associated(id%IRHS_PTR) )THEN
              id%INFO(1)=-22
              id%INFO(2)=12
              GOTO 333
            ENDIF
            IF (size(id%IRHS_PTR) < id%NRHS + 1) THEN
              id%INFO(1)=-22
              id%INFO(2)=12
              GOTO 333
            END IF
            IF (id%IRHS_PTR(id%NRHS + 1).ne.id%NZ_RHS+1) THEN
              id%INFO(1)=-27
              id%INFO(2)=id%IRHS_PTR(id%NRHS+1)
              GOTO 333
            END IF
            IF (dble(id%N)*dble(id%NRHS).LT.dble(id%NZ_RHS)) THEN
              IF (PROKG) THEN
                 write(MPG,*)id%MYID,
     &              " Incompatible values for sparse RHS ",
     &              " id%NZ_RHS,id%N,id%NRHS =", 
     &              id%NZ_RHS,id%N,id%NRHS
              ENDIF
              id%INFO(1)=-22
              id%INFO(2)=11
              GOTO 333
            END IF
            IF (id%IRHS_PTR(1).ne.1) THEN
              id%INFO(1)=-28
              id%INFO(2)=id%IRHS_PTR(1)
              GOTO 333
            END IF
            IF (size(id%IRHS_SPARSE) < id%NZ_RHS) THEN
              id%INFO(1)=-22
              id%INFO(2)=11
              GOTO 333
            END IF
            IF (size(id%RHS_SPARSE) < id%NZ_RHS) THEN
              id%INFO(1)=-22
              id%INFO(2)=10
              GOTO 333
            END IF
          ENDIF
          CALL SMUMPS_634(ICNTL(1),KEEP(1),MPG,INFO(1))
          IF (INFO(1) .LT. 0) GOTO 333
          IF (KEEP(111).eq.-1.AND.id%NRHS.NE.KEEP(112)+KEEP(17))THEN
                INFO(1)=-32
                INFO(2)=id%NRHS
                GOTO 333
          ENDIF
          IF (KEEP(111).gt.0 .AND. id%NRHS .NE. 1) THEN
                INFO(1)=-32
                INFO(2)=id%NRHS
                GOTO 333
          ENDIF
          IF (KEEP(111) .NE. 0 .AND. id%KEEP(50) .EQ. 0) THEN
              IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: null space not available for unsymmetric matrices'
              INFO(1) = -37
              INFO(2) = 0  
              GOTO 333
          ENDIF
          IF (KEEP(248) .NE.0.AND.KEEP(111).NE.0) THEN
              IF (PROKG) WRITE(MPG,'(A)')
     &        ' ERROR: ICNTL(20) and ICNTL(30) functionalities ',
     &        ' incompatible with null space'
              INFO(1) = -37
              IF (KEEP(237).NE.0) THEN
                 INFO(2) = 30 
                 IF (PROKG) WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(30) functionality ',
     &           ' incompatible with null space'
              ELSE
                 IF (PROKG) WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(20) functionality ',
     &           ' incompatible with null space'
                 INFO(2) = 20  
              ENDIF
              GOTO 333
          ENDIF
          IF (( KEEP(111) .LT. -1 ) .OR.
     &         (KEEP(111).GT.KEEP(112)+KEEP(17)) .OR.
     &         (KEEP(111) .EQ.-1 .AND. KEEP(112)+KEEP(17).EQ.0))
     &         THEN
                 INFO(1)=-36
                 INFO(2)=KEEP(111)
                 GOTO 333
          ENDIF
      END IF  
      IF (ICNTL21==1) THEN
          IF ( id%MYID .ne. MASTER  .OR.
     &       ( id%MYID .eq. MASTER .AND.
     &               id%KEEP(46) .eq. 1 ) ) THEN
            IF ( id%LSOL_loc < id%KEEP(89) ) THEN
              id%INFO(1)= -29
              id%INFO(2)= id%LSOL_loc
              GOTO 333
            ENDIF
            IF (id%KEEP(89) .NE. 0) THEN
              IF ( .not. associated(id%ISOL_loc) )THEN
                id%INFO(1)=-22
                id%INFO(2)=13
                GOTO 333
              ENDIF
              IF ( .not. associated(id%SOL_loc) )THEN
                id%INFO(1)=-22
                id%INFO(2)=14
                GOTO 333
              ENDIF
              IF (size(id%ISOL_loc) < id%KEEP(89) ) THEN
                id%INFO(1)=-22
                id%INFO(2)=13
                GOTO 333
              END IF
              IF (size(id%SOL_loc) < 
     &              (id%NRHS-1)*id%LSOL_loc+id%KEEP(89)) THEN  
                id%INFO(1)=-22
                id%INFO(2)=14
                GOTO 333
              END IF
            ENDIF
          ENDIF
      ENDIF
      IF (id%MYID .NE. MASTER) THEN
          IF (id%KEEP(248) == 1) THEN
           IF ( associated( id%RHS ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 7
             GOTO 333
           END IF
           IF ( associated( id%RHS_SPARSE ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 10
             GOTO 333
           END IF
           IF ( associated( id%IRHS_SPARSE ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 11
             GOTO 333
           END IF
           IF ( associated( id%IRHS_PTR ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 12
             GOTO 333
           END IF
          END IF
      ENDIF
      IF (id%MYID.EQ.MASTER) THEN
          CALL SMUMPS_769(id)
      END IF 
      IF (id%INFO(1) .LT. 0) GOTO 333
 333  CONTINUE
      CALL MUMPS_276( id%ICNTL(1),
     &                      id%INFO(1),
     &                      id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GO TO 90
      IF ((id%KEEP(248).EQ.1).AND.(id%KEEP(237).EQ.0)) THEN
       CALL MPI_BCAST(id%NZ_RHS,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
       IF (id%NZ_RHS.EQ.0) THEN
        IF ((ICNTL21.EQ.1).AND.(I_AM_SLAVE)) THEN
          LIW_PASSED=max(1,KEEP(32))
          IF (KEEP(89) .GT. 0) THEN
            CALL SMUMPS_535( MTYPE, id%ISOL_loc(1),
     &               id%PTLUST_S(1),
     &               id%KEEP(1),id%KEEP8(1),
     &               id%IS(1), LIW_PASSED,id%MYID_NODES,
     &               id%N, id%STEP(1), id%PROCNODE_STEPS(1),
     &               id%NSLAVES, scaling_data, LSCAL )
            DO J=1, id%NRHS
              DO I=1, KEEP(89)
                id%SOL_loc((J-1)*id%LSOL_loc + I) =ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        IF (ICNTL21.NE.1) THEN  
         IF (id%MYID.EQ.MASTER) THEN
            DO J=1, id%NRHS
              DO I=1, id%N
                id%RHS((J-1)*id%LRHS + I) =ZERO
              ENDDO
            ENDDO
         ENDIF
        ENDIF
        IF ( PROKG )  THEN 
           WRITE( MPG, 150 )
     &        id%NRHS, ICNTL(27), ICNTL(9), ICNTL(10), ICNTL(11),
     &        ICNTL(20), ICNTL(21), ICNTL(30)
           IF (KEEP(221).NE.0) THEN 
            WRITE (MPG, 152) KEEP(221)
           ENDIF
           IF (KEEP(252).GT.0) THEN   
            WRITE (MPG, 153) KEEP(252)
           ENDIF
        ENDIF
        GOTO 90 
       ENDIF  
      ENDIF 
      IF (id%MYID.EQ.MASTER) THEN
         IF ((KEEP(111).NE.0)) THEN
            KEEP(242)   = 0
         ENDIF
      ENDIF
      INTERLEAVE_PAR   =.FALSE.  
      DO_PERMUTE_RHS   =.FALSE.  
      IF ((id%KEEP(235).NE.0).or.(id%KEEP(237).NE.0)) THEN
         IF (id%KEEP(237).NE.0.AND.
     &       id%KEEP(248).EQ.0) THEN
          IF (LP.GT.0) THEN
           WRITE(LP,'(A,I4,I4)') 
     &     ' Internal Error in solution driver (A-1) ', 
     &       id%KEEP(237), id%KEEP(248)
          ENDIF
          CALL MUMPS_ABORT()
         ENDIF
         NBT = 0
         CALL MUMPS_733(id%Step2node, id%KEEP(28), id%INFO, LP, 
     &        FORCE=.TRUE.,
     &        STRING='id%Step2node (Solve)', MEMCNT=NBT, ERRCODE=-13)
         CALL MUMPS_276( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF ( INFO(1).LT.0 ) RETURN
         IF (NBT.NE.0) THEN
          DO I=1, id%N
           IF (id%STEP(I).LE.0) CYCLE  
           id%Step2node(id%STEP(I)) = I
          ENDDO
         ENDIF
         NB_BYTES = NB_BYTES + int(NBT,8)*K34_8
         NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
         NB_BYTES_EXTRA = NB_BYTES_EXTRA + int(NBT,8) * K34_8
      ENDIF 
      IF ( I_AM_SLAVE )
     &  CALL MUMPS_804(id%OOC_SIZE_OF_BLOCK, id%KEEP(201))
      DO_NULL_PIV = .TRUE.
      NBCOL_INBLOC = -9998
      NZ_THIS_BLOCK= -9998
      JBEG_RHS  = -9998
      IF (id%MYID.EQ.MASTER) THEN 
        IF( KEEP(111)==0 .AND. KEEP(248)==1 ) THEN
           NRHS_NONEMPTY = 0
           DO I=1, id%NRHS
              IF (id%IRHS_PTR(I).LT.id%IRHS_PTR(I+1)) 
     &             NRHS_NONEMPTY = NRHS_NONEMPTY+1  
           ENDDO
           IF (NRHS_NONEMPTY.LE.0) THEN
            IF (LP.GT.0) 
     &        WRITE(LP,*) 'Internal error : NRHS_NONEMPTY=',
     &        NRHS_NONEMPTY
              CALL MUMPS_ABORT()
           ENDIF
        ELSE
           NRHS_NONEMPTY = id%NRHS
        ENDIF
      ENDIF
      BUILD_POSINRHSCOMP = .TRUE.
      IF (KEEP(221).EQ.2 .AND. KEEP(252).EQ.0) THEN
         BUILD_POSINRHSCOMP = .FALSE.
      ENDIF
      SIZE_ROOT   = -33333
      IF ( KEEP( 38 ) .ne. 0 ) THEN
            MASTER_ROOT = MUMPS_275(
     &                    id%PROCNODE_STEPS(id%STEP( KEEP(38))),
     &                    id%NSLAVES )
            IF (id%MYID_NODES .eq. MASTER_ROOT) THEN
              SIZE_ROOT = id%root%TOT_ROOT_SIZE
            ELSE IF ((id%MYID.EQ.MASTER).AND.KEEP(60).NE.0) THEN 
              SIZE_ROOT=id%KEEP(116)
            ENDIF
      ELSE IF (KEEP( 20 ) .ne. 0 ) THEN
            MASTER_ROOT = MUMPS_275(
     &                    id%PROCNODE_STEPS(id%STEP(KEEP(20))),
     &                    id%NSLAVES )
            IF (id%MYID_NODES .eq. MASTER_ROOT) THEN
              SIZE_ROOT = id%IS(
     &               id%PTLUST_S(id%STEP(KEEP(20)))+KEEP(IXSZ) + 3)
            ELSE IF ((id%MYID.EQ.MASTER).AND.KEEP(60).NE.0) THEN 
              SIZE_ROOT=id%KEEP(116)
            ENDIF
      ELSE
            MASTER_ROOT = -44444
      END IF
      IF (id%MYID .eq. MASTER) THEN
        KEEP(84) = ICNTL(27)
        IF (KEEP(252).NE.0) THEN
          NBRHS = KEEP(253)
        ELSE
          IF (KEEP(201) .EQ. 0 .OR. KEEP(84) .GT. 0) THEN
            NBRHS = abs(KEEP(84))
          ELSE
            NBRHS = -2*KEEP(84)
          END IF
          IF (NBRHS .GT. NRHS_NONEMPTY ) NBRHS = NRHS_NONEMPTY
        ENDIF
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(glob_comm_ini,IERR)
#endif
      CALL MPI_BCAST(NRHS_NONEMPTY,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(NBRHS,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      IF (KEEP(201).GT.0) THEN
          IF (I_AM_SLAVE) THEN
            IF (KEEP(201).EQ.1 
     &        .AND.KEEP(50).EQ.0 
     &        .AND.KEEP(251).NE.2 
     &        ) THEN
              OOC_NB_FILE_TYPE=2 
            ELSE
              OOC_NB_FILE_TYPE=1 
            ENDIF
          ENDIF
          WORKSPACE_MINIMAL_PREFERRED = .FALSE.
          IF (id%MYID .eq. MASTER) THEN
             KEEP(107) = max(0,KEEP(107))
             IF ((KEEP(107).EQ.0).AND.
     &            (KEEP(204).EQ.0).AND.(KEEP(211).NE.1) ) THEN
              WORKSPACE_MINIMAL_PREFERRED=.TRUE.
             ENDIF
          ENDIF
          CALL MPI_BCAST( KEEP(107), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(204), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(208), 2, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( WORKSPACE_MINIMAL_PREFERRED, 1,
     &                  MPI_LOGICAL,
     &                  MASTER, id%COMM, IERR )
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        NB_K133     = 3
        IF ( KEEP( 38 ) .NE. 0 .OR. KEEP( 20 ) .NE. 0 ) THEN 
          IF ( MASTER_ROOT .eq. id%MYID_NODES ) THEN
            IF (
     &          .NOT. associated(id%root%RHS_CNTR_MASTER_ROOT)
     &         ) THEN
                NB_K133 = NB_K133 + 1
            ENDIF
          END IF
        ENDIF
        LWCB_MIN = NB_K133*KEEP(133)*NBRHS
        WK_USER_PROVIDED = (id%LWK_USER.NE.0)
        IF (id%LWK_USER.EQ.0) THEN
          ITMP8 = 0_8
        ELSE IF (id%LWK_USER.GT.0) THEN
          ITMP8= int(id%LWK_USER,8)
        ELSE
          ITMP8 = -int(id%LWK_USER,8)* 1000000_8 
        ENDIF
        IF (KEEP(201).EQ.0) THEN  
          IF (ITMP8.NE.KEEP8(24)) THEN
            INFO(1) = -41
            INFO(2) = id%LWK_USER
            GOTO 99    
          ENDIF
        ELSE
          KEEP8(24)=ITMP8
        ENDIF
        MAXS = 0_8
        IF (WK_USER_PROVIDED) THEN
           MAXS = KEEP8(24)
           IF (MAXS.LT. KEEP8(20)) THEN 
                  INFO(1)= -11
                  ITMP8  = KEEP8(20)+1_8-MAXS
                  CALL  MUMPS_731(ITMP8, INFO(2))
           ENDIF
           IF (INFO(1) .GE. 0 ) id%S => id%WK_USER(1:KEEP8(24))
        ELSE IF (associated(id%S)) THEN 
           MAXS = KEEP8(23)
        ELSE
          IF (KEEP(201).EQ.0) THEN  
            WRITE(*,*) ' Working array S not allocated ',
     &                ' on entry to solve phase (in core) '
            CALL MUMPS_ABORT()
          ELSE
            IF ( KEEP(209).EQ.-1 .AND. WORKSPACE_MINIMAL_PREFERRED)
     &        THEN 
              MAXS = KEEP8(20) + 1_8
            ELSE IF ( KEEP(209) .GE.0 ) THEN
              MAXS = max(int(KEEP(209),8), KEEP8(20) + 1_8)
            ELSE
              MAXS  = id%KEEP8(14) 
            ENDIF
            ALLOCATE (id%S(MAXS), stat = allocok)
            KEEP8(23)=MAXS
            IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem allocation of S at solve'
              INFO(1) = -13
              CALL MUMPS_731(MAXS, INFO(2))
              NULLIFY(id%S)
              KEEP8(23)=0_8
            ENDIF
            NB_BYTES = NB_BYTES + KEEP8(23) * K35_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ENDIF
        ENDIF
        IF(KEEP(201).EQ.0)THEN
           LA  = KEEP8(31)
        ELSE
           LA = MAXS
           IF(MAXS.GT.KEEP8(31)+KEEP8(20)*int(KEEP(107)+1,8))THEN
             LA=KEEP8(31)+KEEP8(20)*int(KEEP(107)+1,8)
           ENDIF
        ENDIF
        IF ( MAXS-LA .GT. int(LWCB_MIN,8) ) THEN
           TMP_LWCB8 = min( MAXS - LA, int(huge(LWCB),8) )
           LWCB      = int( TMP_LWCB8, kind(LWCB) )
           WORK_WCB => id%S(LA+1_8:LA+TMP_LWCB8)
           WORK_WCB_ALLOCATED=.FALSE.
        ELSE
           LWCB = LWCB_MIN
           ALLOCATE(WORK_WCB(LWCB_MIN), stat = allocok)
           IF (allocok < 0 ) THEN
                   INFO(1)=-13
                   INFO(2)=LWCB_MIN
           ENDIF
           WORK_WCB_ALLOCATED=.TRUE.
           NB_BYTES = NB_BYTES + int(size(WORK_WCB),8)*K35_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        ENDIF
      ENDIF 
  99  CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) < 0) GOTO 90
      IF ( I_AM_SLAVE ) THEN
        IF (KEEP(201).GT.0) THEN
          CALL SMUMPS_590(LA)
          CALL SMUMPS_586(id)
          IS_INIT_OOC_DONE = .TRUE.
        ENDIF
      ENDIF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) < 0) GOTO 90
      IF (id%MYID .eq. MASTER) THEN
         IF ( (KEEP(242).NE.0) .or. (KEEP(243).NE.0) ) THEN
           IF ( (KEEP(237) .EQ.0) .and. (KEEP(111).EQ.0) ) THEN
            KEEP(242) = 0
            KEEP(243) = 0
           ENDIF
         ENDIF
        IF ( PROKG )  THEN 
           WRITE( MPG, 150 )
     &        id%NRHS, NBRHS, ICNTL(9), ICNTL(10), ICNTL(11),
     &        ICNTL(20), ICNTL(21), ICNTL(30)
           IF (KEEP(111).NE.0) THEN    
            WRITE (MPG, 151) KEEP(111)
           ENDIF
           IF (KEEP(221).NE.0) THEN 
            WRITE (MPG, 152) KEEP(221)
           ENDIF
           IF (KEEP(252).GT.0) THEN   
            WRITE (MPG, 153) KEEP(252)
           ENDIF
        ENDIF
        LSCAL = (((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 8)) .OR. (
     &    KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        ERANAL = ((ICNTL(11) .GT. 0) .OR. (ICNTL(10) .GT. 0))
        IF ( (KEEP(55).eq.0) .AND. KEEP(54).eq.0 .AND. 
     &      .NOT.associated(id%A) ) THEN
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        ELSE
          ICNTL10 = ICNTL(10)
          ICNTL11 = ICNTL(11)
        ENDIF
        IF ((KEEP(111).NE.0).OR.(KEEP(237).NE.0).OR.
     &       (KEEP(252).NE.0) ) THEN
          IF (ICNTL10 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(10) treated as if set to 0 '
          ENDIF
          IF (ICNTL11 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(11) treated as if set to 0 '
          ENDIF
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        END IF
        IF (KEEP(221).NE.0) THEN
          IF (ICNTL10 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(10) treated as if set to 0 (reduced RHS))'
          ENDIF
          IF (ICNTL11 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: ICNTL(11) treated as if set to 0 (reduced RHS)'
          ENDIF
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        END IF
        IF ((ERANAL .AND. NBRHS > 1) .OR. ICNTL(21) > 0) THEN
          IF (ICNTL11 > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &     ' WARNING: ICNTL(11) treated as if set to zero'
            ICNTL11=0
          ENDIF
          IF (ICNTL10 > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     &     ' WARNING: ICNTL(10) treated as if set to zero'
            ICNTL10=0
          ENDIF
          ERANAL = .FALSE.     
        ENDIF
        IF (ERANAL) THEN
            ALLOCATE(SAVERHS(id%N*NBRHS),stat = allocok)
            IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem in solve: error allocating SAVERHS'
              INFO(1) = -13
              INFO(2) = id%N*NBRHS
              GOTO 111
            END IF
            NB_BYTES = NB_BYTES + int(size(SAVERHS),8)*K35_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        ENDIF
        IF (KEEP(237).NE.0 .AND.KEEP(111).NE.0) THEN
          IF (PROKG) WRITE(MPG,'(A)')
     &    ' WARNING: KEEP(237) treated as if set to 0 (null space)'
          KEEP(237)=0
        ENDIF
        IF (KEEP(242).EQ.0) KEEP(243)=0 
      END IF
      CALL MPI_BCAST(ICNTL10,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL11,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL21,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(ERANAL,1,MPI_LOGICAL,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(LSCAL,1,MPI_LOGICAL,MASTER,
     &               id%COMM,IERR)
       CALL MPI_BCAST(KEEP(111),1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(KEEP(242),1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      CALL MPI_BCAST(KEEP(243),1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      DO_PERMUTE_RHS = (KEEP(242).NE.0)
      IF ( KEEP(242).NE.0) THEN
         IF ((KEEP(237).EQ.0).AND.(KEEP(111).EQ.0)) THEN
          IF (MP.GT.0) THEN
           write(MP,*) ' Warning incompatible options ',
     &       ' permute RHS reset to false '
          ENDIF
          DO_PERMUTE_RHS = .FALSE.
         ENDIF
      ENDIF
      IF ( (id%NSLAVES.GT.1) .AND. (KEEP(243).NE.0)
     &        )  THEN
         IF ((KEEP(237).NE.0).or.(KEEP(111).GT.0)) THEN
            INTERLEAVE_PAR= .TRUE.
         ELSE
          IF (PROKG) THEN
           write(MPG,*) ' Warning incompatible options ',
     &       ' interleave RHS reset to false '
          ENDIF
         ENDIF
      ENDIF
#if defined(check)
      IF ( id%MYID_NODES .EQ. MASTER ) THEN
          WRITE(*,*) " ES  A-1  DO_Perm   Interleave ="
          WRITE(*,144) id%KEEP(235), id%KEEP(237),
     &         id%KEEP(242),id%KEEP(243)
      ENDIF
#endif
        MSG_MAX_BYTES_SOLVE =  ( 4 + KEEP(133) ) * KEEP(34) +
     &                           KEEP(133) * NBRHS * KEEP(35)
     &  + 16 * KEEP(34) 
        IF (KEEP(237).EQ.0) THEN
         KMAX_246_247 = max(KEEP(246),KEEP(247))
         MSG_MAX_BYTES_GTHRSOL =  ( ( 2 + KMAX_246_247 ) * KEEP(34) +
     &                              KMAX_246_247 * NBRHS * KEEP(35) )
        ELSE
         MSG_MAX_BYTES_GTHRSOL =  (  3  * KEEP(34) + KEEP(35) )
        ENDIF
        id%LBUFR_BYTES = max(MSG_MAX_BYTES_SOLVE, MSG_MAX_BYTES_GTHRSOL)
        TSIZE = int(min(100_8*int(MSG_MAX_BYTES_GTHRSOL,8),
     &              10000000_8))
        id%LBUFR_BYTES = max(id%LBUFR_BYTES,TSIZE)
        id%LBUFR = ( id%LBUFR_BYTES + KEEP(34) - 1 ) / KEEP(34)
        IF ( associated (id%BUFR) ) THEN 
          NB_BYTES = NB_BYTES - int(size(id%BUFR),8)*K34_8
          DEALLOCATE(id%BUFR)
          NULLIFY(id%BUFR)
        ENDIF
        ALLOCATE (id%BUFR(id%LBUFR),stat=allocok)
        IF ( allocok .GT. 0 ) THEN
            IF (LP.GT.0) 
     &      WRITE(LP,*) id%MYID, 
     &      ' Problem in solve: error allocating BUFR'
            INFO(1) = -13
            INFO(2) = id%LBUFR
            GOTO 111
        ENDIF
        NB_BYTES = NB_BYTES + int(size(id%BUFR),8)*K34_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
      IF ( I_AM_SLAVE ) THEN
        SMUMPS_LBUF_INT = ( 20 + id%NSLAVES * id%NSLAVES  * 4 )
     &                 * KEEP(34)
        CALL SMUMPS_55( SMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = SMUMPS_LBUF_INT
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     &      ':Error allocating small Send buffer:IERR=',IERR
          END IF
          GOTO 111
        END IF
        SMUMPS_LBUF = (MSG_MAX_BYTES_SOLVE + 2*KEEP(34) )*id%NSLAVES
        SMUMPS_LBUF = min(SMUMPS_LBUF, 100 000 000)
        SMUMPS_LBUF = max(SMUMPS_LBUF,
     &      (MSG_MAX_BYTES_SOLVE+2*KEEP(34)) * min(id%NSLAVES,3))
        SMUMPS_LBUF = SMUMPS_LBUF + KEEP(34)
        CALL SMUMPS_53( SMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = SMUMPS_LBUF/KEEP(34) + 1   
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     &      ':Error allocating Send buffer:IERR=', IERR
          END IF
          GOTO 111
        END IF
      ENDIF
      IF (
     &  ( id%MYID .NE. MASTER ) 
     &     .or.
     &    ( I_AM_SLAVE .AND. id%MYID .EQ. MASTER .AND. 
     &      ICNTL21 .NE.0 .AND.
     &      ( KEEP(248).ne.0 .OR. KEEP(221).EQ.2 
     &          .OR. KEEP(111).NE.0 ) 
     &    )
     &     .or.
     &    ( id%MYID .EQ. MASTER .AND. (KEEP(237).NE.0) )
     &    ) THEN
        ALLOCATE(RHS_MUMPS(id%N*NBRHS),stat=IERR)
        NB_BYTES = NB_BYTES + int(size(RHS_MUMPS),8)*K35_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        IF ( IERR .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=id%N*NBRHS
          IF (LP > 0)
     &      WRITE(LP,*) 'ERROR while allocating RHS on a slave'
          GOTO 111
        END IF
      ELSE
        RHS_MUMPS=>id%RHS
      ENDIF 
      IF ( I_AM_SLAVE ) THEN
        LD_RHSCOMP = max(KEEP(89),1)
        IF (id%MYID.EQ.MASTER) THEN
            LD_RHSCOMP = max(LD_RHSCOMP, KEEP(247))
        ENDIF
        IF (KEEP(221).EQ.2 .AND. KEEP(252).EQ.0) THEN
           IF (.NOT.associated(id%RHSCOMP)) THEN
             INFO(1) = -35
             INFO(2) = 1
             GOTO 111
           ENDIF
           IF (.NOT.associated(id%POSINRHSCOMP)) THEN
             INFO(1) = -35
             INFO(2) = 2
             GOTO 111
           ENDIF
           LENRHSCOMP = size(id%RHSCOMP)
           LD_RHSCOMP = LENRHSCOMP/id%NRHS
        ELSE IF (KEEP(221).EQ.1) THEN
          IF (associated(id%RHSCOMP)) THEN 
            NB_BYTES = NB_BYTES - int(size(id%RHSCOMP),8)*K35_8
            DEALLOCATE(id%RHSCOMP)
          ENDIF
          LENRHSCOMP = LD_RHSCOMP*id%NRHS
          ALLOCATE (id%RHSCOMP(LENRHSCOMP))
          NB_BYTES = NB_BYTES + int(size(id%RHSCOMP),8)*K35_8
          IF (associated(id%POSINRHSCOMP)) THEN 
            NB_BYTES = NB_BYTES - int(size(id%POSINRHSCOMP),8)*K34_8
            DEALLOCATE(id%POSINRHSCOMP)
          ENDIF
          ALLOCATE (id%POSINRHSCOMP(KEEP(28)) )
          NB_BYTES = NB_BYTES + int(size(id%POSINRHSCOMP),8)*K34_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        ELSE
          LENRHSCOMP = LD_RHSCOMP*NBRHS
          IF (associated(id%RHSCOMP)) THEN 
           NB_BYTES = NB_BYTES - int(size(id%RHSCOMP),8)*K35_8
           DEALLOCATE(id%RHSCOMP)
           NULLIFY(id%RHSCOMP)
          ENDIF
          ALLOCATE (id%RHSCOMP(LENRHSCOMP))
          NB_BYTES = NB_BYTES + int(size(id%RHSCOMP),8)*K35_8
          IF (associated(id%POSINRHSCOMP)) THEN 
            NB_BYTES = NB_BYTES - int(size(id%POSINRHSCOMP),8)*K34_8
            DEALLOCATE(id%POSINRHSCOMP)
          ENDIF
          ALLOCATE (id%POSINRHSCOMP(KEEP(28)) )
          NB_BYTES = NB_BYTES + int(size(id%POSINRHSCOMP),8)*K34_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        ENDIF
        LIWK_SOLVE = 4 * KEEP(28) + 1
        IF (KEEP(201).EQ.1) THEN
          LIWK_SOLVE = LIWK_SOLVE + KEEP(228) + 1
        ELSE
          LIWK_SOLVE = LIWK_SOLVE + 1
        ENDIF
        ALLOCATE ( IWK_SOLVE( LIWK_SOLVE), stat = allocok )
        IF (allocok .GT. 0 ) THEN
         INFO(1)=-13
         INFO(2)=LIWK_SOLVE
         GOTO 111
        END IF
        NB_BYTES = NB_BYTES + int(LIWK_SOLVE,8)*K34_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        LIWCB =  20*NB_K133*2 + KEEP(133) 
        ALLOCATE ( IWCB( LIWCB), stat = allocok )
        IF (allocok .GT. 0 ) THEN
         INFO(1)=-13
         INFO(2)=LIWCB
         GOTO 111
        END IF
        NB_BYTES = NB_BYTES + int(LIWCB,8)*K34_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        LIW = KEEP(32)
        ALLOCATE(SRW3(KEEP(133)), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=KEEP(133)
          GOTO 111
        END IF
        NB_BYTES = NB_BYTES + int(size(SRW3),8)*K35_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        IF ( (KEEP(111).NE.0) .OR. (KEEP(248).NE.0) ) THEN
          ALLOCATE(POSINRHSCOMP_N(id%N), stat=allocok)
          IF ( allocok .GT. 0 ) THEN
            IF (LP.GT.0) WRITE(LP,*)
     &     ' ERROR in SMUMPS_301: allocating POSINRHSCOMP_N'
            INFO(1) = -13
            INFO(2) = id%N
            GOTO 111
          END IF
          NB_BYTES = NB_BYTES + int(size(POSINRHSCOMP_N),8)*K34_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        END IF
      ELSE
        LIW=0
      END IF
      IF (allocated(UNS_PERM_INV)) DEALLOCATE(UNS_PERM_INV)
      IF ( ( id%MYID .eq. MASTER.AND.(KEEP(23).GT.0) .AND.
     &         (MTYPE .NE. 1).AND.(KEEP(248).NE.0) 
     &       )
     &      .OR. ( ( KEEP(237).NE.0 ) .AND. (KEEP(23).NE.0)  )
     &      )  THEN
           ALLOCATE(UNS_PERM_INV(id%N),stat=allocok) 
           if (allocok .GT.0 ) THEN
             INFO(1)=-13
             INFO(2)=id%N
             GOTO 111
           endif
           NB_BYTES = NB_BYTES + int(id%N,8)*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
           IF (id%MYID.EQ.MASTER) THEN
            DO I = 1, id%N
              UNS_PERM_INV(id%UNS_PERM(I))=I
            ENDDO
           ENDIF
        ELSE
           ALLOCATE(UNS_PERM_INV(1), stat=allocok)
           if (allocok .GT.0 ) THEN
             INFO(1)=-13
             INFO(2)=1
             GOTO 111
           endif
           NB_BYTES = NB_BYTES + 1_8*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
        ENDIF
 111  CONTINUE
#if defined(V_T)
      CALL VTEND(glob_comm_ini,IERR)
#endif
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
        IF ( ( KEEP(237).NE.0 ) .AND. (KEEP(23).NE.0)  ) THEN
          CALL MPI_BCAST(UNS_PERM_INV,id%N,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      ENDIF
      IF ( ICNTL21==1 ) THEN
        IF (LSCAL) THEN
          IF (id%MYID.NE.MASTER) THEN
            IF (MTYPE == 1) THEN
              ALLOCATE(id%COLSCA(id%N),stat=allocok)
            ELSE
              ALLOCATE(id%ROWSCA(id%N),stat=allocok)
            ENDIF
            IF (allocok > 0) THEN
              IF (LP > 0) THEN
               WRITE(LP,*) 'Error allocating temporary scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
            NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ENDIF
          IF (MTYPE == 1) THEN
              CALL MPI_BCAST(id%COLSCA(1),id%N,
     &                       MPI_REAL,MASTER,
     &                       id%COMM,IERR)
              scaling_data%SCALING=>id%COLSCA
          ELSE
              CALL MPI_BCAST(id%ROWSCA(1),id%N,
     &                       MPI_REAL,MASTER,
     &                       id%COMM,IERR)
              scaling_data%SCALING=>id%ROWSCA
          ENDIF
          IF (I_AM_SLAVE) THEN
            ALLOCATE(scaling_data%SCALING_LOC(id%KEEP(89)),
     &               stat=allocok)
            IF (allocok > 0) THEN
              IF (LP > 0) THEN
                WRITE(LP,*) 'Error allocating local scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%KEEP(89)
              GOTO 40
            ENDIF
            NB_BYTES = NB_BYTES + int(id%KEEP(89),8)*K16_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ENDIF
        ENDIF
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED=max(1,LIW)
          IF (KEEP(89) .GT. 0) THEN
            CALL SMUMPS_535( MTYPE, id%ISOL_loc(1),
     &               id%PTLUST_S(1),
     &               id%KEEP(1),id%KEEP8(1),
     &               id%IS(1), LIW_PASSED,id%MYID_NODES,
     &               id%N, id%STEP(1), id%PROCNODE_STEPS(1),
     &               id%NSLAVES, scaling_data, LSCAL )
          ENDIF
          IF (id%MYID.NE.MASTER .AND. LSCAL) THEN
            IF (MTYPE == 1) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ELSE
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
            NB_BYTES = NB_BYTES - int(id%N,8)*K16_8
          ENDIF
        ENDIF
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
          IF (id%MYID.NE.MASTER) THEN
            ALLOCATE(id%UNS_PERM(id%N),stat=allocok)
            IF (allocok > 0) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
          ENDIF
        ENDIF
 40     CONTINUE
        CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
        IF (INFO(1) .LT.0 ) GOTO 90
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
          CALL MPI_BCAST(id%UNS_PERM(1),id%N,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
          IF (I_AM_SLAVE) THEN
            DO I=1, KEEP(89)
              id%ISOL_loc(I) = id%UNS_PERM(id%ISOL_loc(I))
            ENDDO
          ENDIF
          IF (id%MYID.NE.MASTER) THEN
            DEALLOCATE(id%UNS_PERM)
            NULLIFY(id%UNS_PERM)
          ENDIF
        ENDIF
      ENDIF
      IF ( ( KEEP(221) .EQ. 1 ) .OR.
     &     ( KEEP(221) .EQ. 2 ) 
     &   ) THEN
         IF (KEEP(46).EQ.1) THEN
             MASTER_ROOT_IN_COMM=MASTER_ROOT
         ELSE
             MASTER_ROOT_IN_COMM =MASTER_ROOT+1
         ENDIF
         IF ( id%MYID .EQ. MASTER ) THEN
             IF (id%NRHS.EQ.1) THEN
               LD_REDRHS = id%KEEP(116)
             ELSE
               LD_REDRHS = id%LREDRHS
             ENDIF
         ENDIF
         IF (MASTER.NE.MASTER_ROOT_IN_COMM) THEN
            IF ( id%MYID .EQ. MASTER ) THEN
             CALL MPI_SEND(LD_REDRHS,1,MPI_INTEGER,
     &       MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
            ELSEIF ( id%MYID.EQ.MASTER_ROOT_IN_COMM) THEN
             CALL MPI_RECV(LD_REDRHS,1,MPI_INTEGER,
     &       MASTER, 0, id%COMM,STATUS,IERR)
            ENDIF
         ENDIF
      ENDIF
      IF ( KEEP(248)==1 ) THEN  
        JEND_RHS = 0 
        IF (DO_PERMUTE_RHS) THEN
           ALLOCATE(PERM_RHS(id%NRHS),stat=allocok)
           IF (allocok > 0) THEN
                INFO(1) = -13
                INFO(2) = id%NRHS
                GOTO 109
           ENDIF
           NB_BYTES = NB_BYTES +  int(id%NRHS,8)*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
           IF (id%MYID.EQ.MASTER) THEN
              STRAT_PERMAM1 = KEEP(242)
                 CALL MUMPS_780
     &             (STRAT_PERMAM1, id%SYM_PERM(1),
     &             id%IRHS_PTR(1), id%NRHS+1,
     &             PERM_RHS, id%NRHS,
     &             IERR
     &           )
           ENDIF
        ENDIF
      ENDIF
109   CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
      IF (id%NSLAVES .EQ. 1) THEN
       IF (DO_PERMUTE_RHS .AND. KEEP(111).NE.0 ) THEN
            WRITE(*,*) id%MYID, ':INTERNAL ERROR 1 : ',
     &                 ' PERMUTE RHS during null space computation ',
     &                 ' not available yet '
            CALL MUMPS_ABORT()
       ENDIF                     
      ELSE
         IF (DO_PERMUTE_RHS .AND. KEEP(111).NE.0 ) THEN 
            WRITE(*,*) id%MYID, ':INTERNAL ERROR 2 : ',
     &                 ' PERMUTE RHS during null space computation ',
     &                 ' not available yet '
            CALL MUMPS_ABORT()
         ENDIF                  
         IF (INTERLEAVE_PAR) THEN
          IF ( KEEP(111).NE.0 ) THEN
            WRITE(*,*) id%MYID, ':INTERNAL ERROR 3 : ',
     &        ' INTERLEAVE RHS during null space computation ',
     &        ' not available yet '
            CALL MUMPS_ABORT()
          ELSE
           IF (id%MYID.EQ.MASTER) THEN
             CALL MUMPS_772
     &           (PERM_RHS, id%NRHS, id%N, id%KEEP(28),
     &           id%PROCNODE_STEPS(1), id%STEP(1), id%NSLAVES,
     &           id%Step2node(1),
     &           IERR)
           ENDIF               
          ENDIF                 
         ENDIF                  
      ENDIF                  
      IF (DO_PERMUTE_RHS.AND.(KEEP(111).EQ.0)) THEN
        CALL MPI_BCAST(PERM_RHS(1),
     &            id%NRHS,
     &            MPI_INTEGER, 
     &            MASTER, id%COMM,IERR)
      ENDIF
       BEG_RHS=1
       DO WHILE (BEG_RHS.LE.NRHS_NONEMPTY)
        NBRHS_EFF    = min(NRHS_NONEMPTY-BEG_RHS+1, NBRHS)
         IF (IRHS_SPARSE_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES - 
     &        int(size(IRHS_SPARSE_COPY),8)*K34_8 
             DEALLOCATE(IRHS_SPARSE_COPY)
             IRHS_SPARSE_COPY_ALLOCATED=.FALSE.
             NULLIFY(IRHS_SPARSE_COPY)
         ENDIF 
         IF (IRHS_PTR_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES - 
     &        int(size(IRHS_PTR_COPY),8)*K34_8 
             DEALLOCATE(IRHS_PTR_COPY)
             IRHS_PTR_COPY_ALLOCATED=.FALSE.
             NULLIFY(IRHS_PTR_COPY)
         ENDIF
         IF (RHS_SPARSE_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES - 
     &        int(size(RHS_SPARSE_COPY),8)*K35_8 
             DEALLOCATE(RHS_SPARSE_COPY)
             RHS_SPARSE_COPY_ALLOCATED=.FALSE.
             NULLIFY(RHS_SPARSE_COPY)
         ENDIF 
        IF (
     &      ( id%MYID .NE. MASTER ) 
     &     .or.
     &    ( I_AM_SLAVE .AND. id%MYID .EQ. MASTER .AND. 
     &      ICNTL21 .NE.0 .AND.
     &      ( KEEP(248).ne.0 .OR. KEEP(221).EQ.2 
     &          .OR. KEEP(111).NE.0 ) 
     &    )
     &     .or.
     &    ( id%MYID .EQ. MASTER .AND. (KEEP(237).NE.0) )
     &    ) THEN
          LD_RHS = id%N
          IBEG   = 1
        ELSE
          IF ( associated(id%RHS) ) THEN
              LD_RHS    = max(id%LRHS, id%N)
          ELSE
              LD_RHS    = id%N
          ENDIF
          IBEG      = (BEG_RHS-1) * LD_RHS + 1
        ENDIF
       JBEG_RHS = BEG_RHS        
       IF ( (id%MYID.EQ.MASTER) .AND. 
     &        KEEP(248)==1  ) THEN
        JBEG_RHS = JEND_RHS + 1
        IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
           DO WHILE ( id%IRHS_PTR(PERM_RHS(JBEG_RHS)) .EQ. 
     &          id%IRHS_PTR(PERM_RHS(JBEG_RHS)+1) ) 
              IF ((KEEP(237).EQ.0).AND.(ICNTL21.EQ.0).AND.
     &              (KEEP(221).NE.1) ) THEN
                 DO I=1, id%N
                   RHS_MUMPS((PERM_RHS(JBEG_RHS) -1)*LD_RHS+I) 
     &                     = ZERO
                 ENDDO
              ENDIF
              JBEG_RHS = JBEG_RHS +1
              CYCLE  
           ENDDO
        ELSE
          DO WHILE( id%IRHS_PTR(JBEG_RHS) .EQ. 
     &          id%IRHS_PTR(JBEG_RHS+1) ) 
            IF ((KEEP(237).EQ.0).AND.(ICNTL21.EQ.0).AND.
     &           (KEEP(221).NE.1)) THEN
              DO I=1, id%N
                RHS_MUMPS((JBEG_RHS -1)*LD_RHS + I) = ZERO
              ENDDO
            ENDIF
            IF (KEEP(221).EQ.1) THEN
              DO I = 1, id%SIZE_SCHUR
                id%REDRHS((JBEG_RHS-1)*LD_REDRHS + I) =  ZERO
              ENDDO
            ENDIF
            JBEG_RHS = JBEG_RHS +1
          ENDDO
        ENDIF                   
        NB_RHSSKIPPED = JBEG_RHS - (JEND_RHS + 1) 
        IF ((KEEP(248).EQ.1).AND.(KEEP(237).EQ.0)
     &         .AND. (ICNTL21.EQ.0)) 
     &         THEN
             IBEG      = (JBEG_RHS-1) * LD_RHS + 1
        ENDIF
       ENDIF
         CALL MPI_BCAST( JBEG_RHS,1, MPI_INTEGER,
     &            MASTER, id%COMM,IERR)
        IF (id%MYID.EQ.MASTER .AND. KEEP(221).NE.0) THEN
          IBEG_REDRHS= (JBEG_RHS-1)*LD_REDRHS + 1
        ELSE
          IBEG_REDRHS=-142424  
        ENDIF
        IF ( I_AM_SLAVE ) THEN
          IF ( KEEP(221).EQ.0 ) THEN
             IBEG_RHSCOMP= 1
          ELSE
             IBEG_RHSCOMP= (JBEG_RHS-1)*LD_RHSCOMP + 1
          ENDIF
        ELSE
          IBEG_RHSCOMP=-152525  
        ENDIF
#if defined(V_T)
      CALL VTBEGIN(perm_scal_ini,IERR)
#endif
      IF (id%MYID .eq. MASTER) THEN
        IF (KEEP(248)==1) THEN
          NBCOL        = 0
          NBCOL_INBLOC = 0
          NZ_THIS_BLOCK = 0
          STOP_AT_NEXT_EMPTY_COL = .FALSE.
          DO I=JBEG_RHS, id%NRHS
            NBCOL_INBLOC = NBCOL_INBLOC +1
            IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
               COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1) 
     &                    - id%IRHS_PTR(PERM_RHS(I))
           ELSE
               COLSIZE = id%IRHS_PTR(I+1) - id%IRHS_PTR(I)
           ENDIF
           IF ((.NOT.STOP_AT_NEXT_EMPTY_COL).AND.(COLSIZE.GT.0).AND.
     &          (KEEP(237).EQ.0)) 
     &         STOP_AT_NEXT_EMPTY_COL =.TRUE.
           IF (COLSIZE.GT.0) THEN 
             NBCOL = NBCOL+1
             NZ_THIS_BLOCK = NZ_THIS_BLOCK + COLSIZE
           ELSE IF (STOP_AT_NEXT_EMPTY_COL) THEN
             NBCOL_INBLOC = NBCOL_INBLOC -1
             NBRHS_EFF = NBCOL
             EXIT
           ENDIF
           IF (NBCOL.EQ.NBRHS_EFF) EXIT
          ENDDO
          IF (NBCOL.NE.NBRHS_EFF) THEN
              WRITE(6,*) 'INTERNAL ERROR 1 in SMUMPS_301 ',
     &            NBCOL, NBRHS_EFF
              call MUMPS_ABORT()
          ENDIF
          ALLOCATE(IRHS_PTR_COPY(NBCOL_INBLOC+1),stat=allocok)
          if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NBCOL_INBLOC+1
              GOTO 30
          endif
          IRHS_PTR_COPY_ALLOCATED = .TRUE.
          NB_BYTES =  NB_BYTES +  int(NBCOL_INBLOC+1,8)*K34_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          JEND_RHS =JBEG_RHS + NBCOL_INBLOC - 1 
          IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
            IPOS = 1
            J = 0
            DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
               J = J+1
               IRHS_PTR_COPY(J) = IPOS
               COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1) 
     &                    - id%IRHS_PTR(PERM_RHS(I))
               IPOS = IPOS + COLSIZE
            ENDDO
          ELSE 
            IPOS = 1
            J = 0
            DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
               J = J+1
               IRHS_PTR_COPY(J) = IPOS
               COLSIZE = id%IRHS_PTR(I+1) 
     &                    - id%IRHS_PTR(I)
               IPOS = IPOS + COLSIZE
            ENDDO
          ENDIF                 
          IRHS_PTR_COPY(NBCOL_INBLOC+1)= IPOS
          IF ( IPOS-1 .NE. NZ_THIS_BLOCK ) THEN
                WRITE(*,*) "Error in compressed copy of IRHS_PTR"
                IERR = 99
                call MUMPS_ABORT()
          ENDIF
          IF (KEEP(23) .NE. 0 .and. MTYPE .NE. 1) THEN
            ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),stat=allocok)
            if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 30
            endif
            IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
            NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*K34_8  
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ELSE
            IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &         (KEEP(237).NE.0)) THEN
             ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK), stat=allocok)
             IF (allocok .GT.0 ) THEN
                  IERR = 99
                  GOTO 30
             ENDIF
             IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
             NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*K34_8  
             NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
             IF ( DO_PERMUTE_RHS.OR.INTERLEAVE_PAR ) THEN
              IPOS = 1
              DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
               COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1) 
     &                    - id%IRHS_PTR(PERM_RHS(I))
               IRHS_SPARSE_COPY(IPOS:IPOS+COLSIZE-1) =
     &          id%IRHS_SPARSE(id%IRHS_PTR(PERM_RHS(I)):
     &           id%IRHS_PTR(PERM_RHS(I)+1) -1)
               IPOS = IPOS + COLSIZE
              ENDDO
             ELSE
               IRHS_SPARSE_COPY = id%IRHS_SPARSE(id%IRHS_PTR(JBEG_RHS):
     &             id%IRHS_PTR(JBEG_RHS)+NZ_THIS_BLOCK-1)
             ENDIF
            ELSE
             IRHS_SPARSE_COPY
     &       =>
     &         id%IRHS_SPARSE(id%IRHS_PTR(JBEG_RHS):
     &             id%IRHS_PTR(JBEG_RHS)+NZ_THIS_BLOCK-1)
            ENDIF
          ENDIF
          IF (LSCAL.OR.DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &        (KEEP(237).NE.0)) THEN
            ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK), stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 30
            endif
            RHS_SPARSE_COPY_ALLOCATED = .TRUE.
            NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*K35_8  
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          ELSE 
            IF  ( KEEP(248)==1 ) THEN
              RHS_SPARSE_COPY
     &         => id%RHS_SPARSE(id%IRHS_PTR(JBEG_RHS):
     &             id%IRHS_PTR(JBEG_RHS)+NZ_THIS_BLOCK-1)
            ELSE
              RHS_SPARSE_COPY
     &           => id%RHS_SPARSE(id%IRHS_PTR(BEG_RHS):
     &                       id%IRHS_PTR(BEG_RHS)+NZ_THIS_BLOCK-1)
            ENDIF
          ENDIF
          IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &        (id%KEEP(237).NE.0)) THEN
             IF (id%KEEP(237).NE.0) THEN
              RHS_SPARSE_COPY = ONE
             ELSE IF (.NOT. LSCAL) THEN
              IPOS = 1
              DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
                COLSIZE = id%IRHS_PTR(PERM_RHS(I)+1) 
     &                    - id%IRHS_PTR(PERM_RHS(I))
                IF (COLSIZE .EQ. 0) CYCLE
                RHS_SPARSE_COPY(IPOS:IPOS+COLSIZE-1) =
     &           id%RHS_SPARSE(id%IRHS_PTR(PERM_RHS(I)):
     &            id%IRHS_PTR(PERM_RHS(I)+1) -1)
                IPOS = IPOS + COLSIZE
              ENDDO
             ENDIF
          ENDIF
        ENDIF  
        IF (KEEP(23) .NE. 0) THEN
          IF (MTYPE .NE. 1) THEN
            IF (KEEP(248)==0) THEN
              ALLOCATE( C_RW2( id%N ),stat =allocok )
              IF ( allocok .GT. 0 ) THEN
                INFO(1)=-13
                INFO(2)=id%N
                IF ( LP .GT. 0 ) THEN
                  WRITE(LP,*) id%MYID,
     &            ':Error allocating C_RW2 in SMUMPS_SOLVE_DRIVE'
                END IF
                GOTO 30
              END IF
              DO K = 1, NBRHS_EFF
               KDEC = IBEG+(K-1)*LD_RHS
               DO I = 1, id%N
                C_RW2(I)=RHS_MUMPS(I-1+KDEC)
               END DO
               DO I = 1, id%N
                JPERM = id%UNS_PERM(I)
                RHS_MUMPS(I-1+KDEC) = C_RW2(JPERM)
               END DO
              END DO
              DEALLOCATE(C_RW2)
            ELSE
              IPOS = 1
              DO I=1, NBCOL_INBLOC
                COLSIZE = IRHS_PTR_COPY(I+1) - IRHS_PTR_COPY(I)
                DO K = 1, COLSIZE
                 JPERM = UNS_PERM_INV(IRHS_SPARSE_COPY(IPOS+K-1))
                 IRHS_SPARSE_COPY(IPOS+K-1) = JPERM
                ENDDO
                IPOS = IPOS + COLSIZE
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        IF (ERANAL) THEN
         IF ( KEEP(248) == 0 ) THEN
          DO K = 1, NBRHS_EFF
            KDEC = IBEG+(K-1)*LD_RHS
            DO I = 1, id%N
              SAVERHS(I+(K-1)*id%N) = RHS_MUMPS(KDEC+I-1)
            END DO
          ENDDO
         ENDIF
        ENDIF
        IF (LSCAL) THEN
         IF (KEEP(248)==0) THEN
          IF (MTYPE .EQ. 1) THEN
            DO K =1, NBRHS_EFF 
             KDEC = (K-1) * LD_RHS + IBEG - 1
             DO I = 1, id%N
              RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) *
     &                            id%ROWSCA(I)
             END DO
            ENDDO
          ELSE
            DO K =1, NBRHS_EFF 
             KDEC = (K-1) * LD_RHS + IBEG - 1
             DO I = 1, id%N
              RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) *
     &                            id%COLSCA(I)
             END DO
            ENDDO
          ENDIF
         ELSE  
          KDEC=id%IRHS_PTR(JBEG_RHS)
          IF ((KEEP(248)==1) .AND. 
     &          (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR.OR.
     &           (id%KEEP(237).NE.0))
     &       ) THEN
             IPOS = 1
             J    = 0
             DO I=JBEG_RHS, JBEG_RHS + NBCOL_INBLOC -1
               J = J+1
               COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
               IF (COLSIZE .EQ. 0) CYCLE
               IF (id%KEEP(237).NE.0) THEN
                IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
                  RHS_SPARSE_COPY(IPOS) =   id%ROWSCA(PERM_RHS(I)) *
     &                     ONE
                ELSE
                  RHS_SPARSE_COPY(IPOS) =   id%ROWSCA(I) * ONE
                ENDIF
               ELSE
                DO K = 1, COLSIZE
                 II = id%IRHS_SPARSE(id%IRHS_PTR(PERM_RHS(I))+K-1)  
                 IF (MTYPE.EQ.1) THEN
                  RHS_SPARSE_COPY(IPOS+K-1) =
     &            id%RHS_SPARSE(id%IRHS_PTR(PERM_RHS(I))+K-1)*
     &            id%ROWSCA(II)
                 ELSE
                  RHS_SPARSE_COPY(IPOS+K-1) =
     &            id%RHS_SPARSE(id%IRHS_PTR(PERM_RHS(I))+K-1)*
     &            id%COLSCA(II)
                 ENDIF
                ENDDO
               ENDIF
               IPOS = IPOS + COLSIZE
              ENDDO
          ELSE
            IF (MTYPE .eq. 1) THEN
             DO IZ=1,NZ_THIS_BLOCK
              I=IRHS_SPARSE_COPY(IZ)
              RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*
     &                            id%ROWSCA(I)
             ENDDO
            ELSE
             DO IZ=1,NZ_THIS_BLOCK
              I=IRHS_SPARSE_COPY(IZ)
              RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*
     &                            id%COLSCA(I)
             ENDDO
            ENDIF
          ENDIF
         ENDIF   
        END IF
      ENDIF
#if defined(V_T)
      CALL VTEND(perm_scal_ini,IERR)
#endif
 30   CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
      IF ( I_AM_SLAVE ) THEN
       IF ( (KEEP(111).NE.0) .OR. (KEEP(248).NE.0) .OR. 
     &      (KEEP(252).NE.0) ) THEN
         IF (BUILD_POSINRHSCOMP) THEN
           IF (KEEP(111).NE.0) THEN
             WHAT      = 2
             MTYPE_LOC = 1
           ELSE IF  (KEEP(252).NE.0) THEN
            WHAT = 0
            MTYPE_LOC = 1  
           ELSE
             WHAT      = 1
             MTYPE_LOC = MTYPE
           ENDIF
           LIW_PASSED=max(1,LIW)
           IF (WHAT.EQ.0) THEN
            CALL SMUMPS_639(id%NSLAVES,id%N,
     &           id%MYID_NODES, id%PTLUST_S(1),
     &           id%KEEP(1),id%KEEP8(1), 
     &           id%PROCNODE_STEPS(1), id%IS(1), LIW_PASSED, 
     &           id%STEP(1), 
     &           id%POSINRHSCOMP(1), IDUMMY, 1, MTYPE_LOC,
     &           WHAT )
           ELSE
            CALL SMUMPS_639(id%NSLAVES,id%N,
     &           id%MYID_NODES, id%PTLUST_S(1),
     &           id%KEEP(1),id%KEEP8(1), 
     &           id%PROCNODE_STEPS(1), id%IS(1), LIW_PASSED, 
     &           id%STEP(1), 
     &           id%POSINRHSCOMP(1), POSINRHSCOMP_N(1), 
     &           id%N, MTYPE_LOC,
     &           WHAT )
           ENDIF
           BUILD_POSINRHSCOMP = .FALSE.
         ENDIF
       ENDIF
      ENDIF
      IF (KEEP(248)==1) THEN
            CALL MPI_BCAST( NBCOL_INBLOC,1, MPI_INTEGER,
     &           MASTER, id%COMM,IERR)
      ELSE 
            NBCOL_INBLOC = NBRHS_EFF  
      ENDIF
      JEND_RHS =JBEG_RHS + NBCOL_INBLOC - 1 
      IF ((KEEP(248).EQ.1).AND.(KEEP(237).EQ.0)) THEN
            CALL MPI_BCAST( NBRHS_EFF,1, MPI_INTEGER,
     &           MASTER, id%COMM,IERR)
            CALL MPI_BCAST(NB_RHSSKIPPED,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR)
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(soln_dist,IERR)
#endif
      IF ((KEEP(111).eq.0).AND.(KEEP(252).EQ.0)) THEN
        IF (KEEP(248) == 0) THEN
          IF ( .NOT.I_AM_SLAVE ) THEN
            CALL SMUMPS_638(id%NSLAVES,id%N, id%MYID, id%COMM,
     &          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     &          JDUMMY, id%KEEP(1), id%KEEP8(1), id%PROCNODE_STEPS(1),
     &          IDUMMY, 1,
     &          id%STEP(1), KDUMMY, 1, BUILD_POSINRHSCOMP,
     &          id%ICNTL(1),id%INFO(1))
            BUILD_POSINRHSCOMP=.FALSE.
          ELSE
            LIW_PASSED = max( LIW, 1 )
            CALL SMUMPS_638(id%NSLAVES,id%N, id%MYID, id%COMM,
     &          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     &          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     &          id%PROCNODE_STEPS(1),
     &          IS(1), LIW_PASSED,
     &          id%STEP(1), id%POSINRHSCOMP(1), KEEP(28), 
     &          BUILD_POSINRHSCOMP,
     &          id%ICNTL(1),id%INFO(1))
            BUILD_POSINRHSCOMP=.FALSE.
          ENDIF
          IF (INFO(1).LT.0) GOTO 90
        ELSE   
         CALL MPI_BCAST( NZ_THIS_BLOCK,1, MPI_INTEGER,
     &                   MASTER, id%COMM,IERR)
         IF (id%MYID.NE.MASTER) THEN
           ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),stat=allocok)
           if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 45
           endif
           IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
           ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK),stat=allocok)
           if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 45
           endif
           RHS_SPARSE_COPY_ALLOCATED=.TRUE.
           NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*(K34_8+K35_8)
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
           ALLOCATE(IRHS_PTR_COPY(NBCOL_INBLOC+1),stat=allocok)
           if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NBCOL_INBLOC+1
              GOTO 45
           endif
           IRHS_PTR_COPY_ALLOCATED = .TRUE.
           NB_BYTES = NB_BYTES + int(NBCOL_INBLOC+1,8)*K34_8
           NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
         ENDIF
 45      CONTINUE
         CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
         IF (INFO(1) .LT.0 ) GOTO 90
         CALL MPI_BCAST(IRHS_SPARSE_COPY(1),
     &                NZ_THIS_BLOCK,
     &                MPI_INTEGER,
     &                MASTER, id%COMM,IERR)
         CALL MPI_BCAST(RHS_SPARSE_COPY(1),
     &                NZ_THIS_BLOCK,
     &                MPI_REAL,
     &                MASTER, id%COMM,IERR)
         CALL MPI_BCAST(IRHS_PTR_COPY(1),
     &            NBCOL_INBLOC+1, 
     &            MPI_INTEGER, 
     &            MASTER, id%COMM,IERR)
         IF (IERR.GT.0) THEN
               WRITE (*,*)'NOT OK FOR ALLOC PTR ON SLAVES'
               call MUMPS_ABORT()
         ENDIF
         IF ( I_AM_SLAVE ) THEN
          IF  (KEEP(237).NE.0)  THEN
               K=1              
               RHS_MUMPS(1:NBRHS_EFF*LD_RHS) = ZERO
               IPOS = 1
               DO I = 1, NBCOL_INBLOC
                  COLSIZE = IRHS_PTR_COPY(I+1) - IRHS_PTR_COPY(I)
                  IF (COLSIZE.GT.0) THEN
                   J = I - 1  + JBEG_RHS
                   IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
                          J = PERM_RHS(J)
                   ENDIF
                   IF (POSINRHSCOMP_N(J).NE.0) THEN
                         RHS_MUMPS((K-1) * LD_RHS + J) = 
     &                       RHS_SPARSE_COPY(IPOS)
                   ENDIF               
                   K = K + 1
                   IPOS = IPOS + COLSIZE  
                  ENDIF
               ENDDO
               IF (K.NE.NBRHS_EFF+1) THEN
                  WRITE(6,*) 'INTERNAL ERROR 2 in SMUMPS_301 ',
     &            K, NBRHS_EFF
                  call MUMPS_ABORT()
               ENDIF
          ELSE
           IF ((KEEP(221).EQ.1).AND.(NB_RHSSKIPPED.GT.0)) THEN
             DO K = JBEG_RHS-NB_RHSSKIPPED, JBEG_RHS-1
              DO I = 1, LD_RHSCOMP
                id%RHSCOMP((K-1)*LD_RHSCOMP + I) =  ZERO
              ENDDO
             ENDDO
           ENDIF
           DO K = 1, NBCOL_INBLOC   
            KDEC = (K-1) * LD_RHS + IBEG - 1
            RHS_MUMPS(KDEC+1:KDEC+id%N) = ZERO
            DO IZ=IRHS_PTR_COPY(K), IRHS_PTR_COPY(K+1)-1
              I=IRHS_SPARSE_COPY(IZ)
                IF (POSINRHSCOMP_N(I).NE.0) THEN
                 RHS_MUMPS(KDEC+I)= RHS_SPARSE_COPY(IZ)
                ENDIF
            ENDDO
           ENDDO
          END IF
         ENDIF
        ENDIF  
      ELSE IF (I_AM_SLAVE) THEN
       IF (KEEP(111).NE.0) THEN
        IF (KEEP(111).GT.0) THEN
          IBEG_GLOB_DEF = KEEP(111)
          IEND_GLOB_DEF = KEEP(111)
        ELSE
          IBEG_GLOB_DEF = BEG_RHS
          IEND_GLOB_DEF = BEG_RHS+NBRHS_EFF-1
        ENDIF
        IF ( id%KEEP(112) .GT. 0 .AND. DO_NULL_PIV) THEN
           IF (IBEG_GLOB_DEF .GT.id%KEEP(112)) THEN
              id%KEEP(235) = 0
              DO_NULL_PIV = .FALSE.
           ENDIF
           IF (IBEG_GLOB_DEF .LT.id%KEEP(112)
     &          .AND. IEND_GLOB_DEF .GT.id%KEEP(112)
     &          .AND. DO_NULL_PIV ) THEN
              IEND_GLOB_DEF = id%KEEP(112)
              id%KEEP(235) = 1  
              DO_NULL_PIV = .FALSE.
           ENDIF
        ENDIF
        IF (id%KEEP(235).NE.0) THEN
          NZ_THIS_BLOCK=IEND_GLOB_DEF-IBEG_GLOB_DEF+1
          ALLOCATE(IRHS_PTR_COPY(NZ_THIS_BLOCK+1),stat=allocok)
          if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 50
          endif
          IRHS_PTR_COPY_ALLOCATED = .TRUE.
          ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),stat=allocok)
          if (allocok .GT.0 ) then
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 50
          endif
          IRHS_SPARSE_COPY_ALLOCATED=.TRUE.
          NB_BYTES = NB_BYTES + int(NZ_THIS_BLOCK,8)*(K34_8+K35_8)  
     &                         + K34_8
          NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          IF (id%MYID.eq.MASTER) THEN
             II = 1
             DO I = IBEG_GLOB_DEF, IEND_GLOB_DEF
               IRHS_PTR_COPY(I-IBEG_GLOB_DEF+1)      = I
               IF (INTERLEAVE_PAR .AND.(id%NSLAVES .GT. 1) )THEN
                  IRHS_SPARSE_COPY(II) = PERM_PIV_LIST(I)
               ELSE
                  IRHS_SPARSE_COPY(II) = id%PIVNUL_LIST(I)
               ENDIF
               II = II +1
             ENDDO
             IRHS_PTR_COPY(NZ_THIS_BLOCK+1) = NZ_THIS_BLOCK+1
          ENDIF
 50      CONTINUE
         CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
         IF (INFO(1) .LT.0 ) GOTO 90
          CALL MPI_BCAST(IRHS_SPARSE_COPY(1),
     &                NZ_THIS_BLOCK,
     &                MPI_INTEGER,
     &                MASTER, id%COMM,IERR)
          CALL MPI_BCAST(IRHS_PTR_COPY(1),
     &                NZ_THIS_BLOCK+1,
     &                MPI_INTEGER,
     &                MASTER, id%COMM,IERR)
          RHS_MUMPS( IBEG : 
     &          (IBEG + (NBRHS_EFF)*LD_RHS -1) ) = ZERO
        ENDIF
        DO K=1, NBRHS_EFF
          KDEC = (K-1) *LD_RHSCOMP
          id%RHSCOMP(KDEC+1:KDEC+LD_RHSCOMP)=ZERO
        END DO
        IF ((KEEP(235).NE.0).AND.  INTERLEAVE_PAR) THEN  
         DO I=IBEG_GLOB_DEF,IEND_GLOB_DEF
           IF (id%MYID_NODES .EQ. MAP_PIVNUL_LIST(I) ) THEN
              JJ= POSINRHSCOMP_N(PERM_PIV_LIST(I))
              IF (JJ.GT.LD_RHSCOMP) THEN
                 WRITE(6,*) ' Internal ERROR JJ, LD_RHSCOMP=',
     &                JJ, LD_RHSCOMP
              ENDIF
              IF (JJ.GT.0) THEN
                IF (KEEP(50).EQ.0) THEN
                   id%RHSCOMP(IBEG_RHSCOMP -1+ 
     &             (I-IBEG_GLOB_DEF)*LD_RHSCOMP + JJ) =  
     &             abs(id%DKEEP(2))
                ELSE
                   id%RHSCOMP(IBEG_RHSCOMP -1+ 
     &             (I-IBEG_GLOB_DEF)*LD_RHSCOMP + JJ) = ONE
                ENDIF
              ENDIF
           ENDIF
         ENDDO
        ELSE
         DO I=max(IBEG_GLOB_DEF,KEEP(220)),
     &       min(IEND_GLOB_DEF,KEEP(220)+KEEP(109)-1)
          JJ= POSINRHSCOMP_N(id%PIVNUL_LIST(I-KEEP(220)+1))
          IF (JJ.GT.0) THEN
           IF (KEEP(50).EQ.0) THEN
            id%RHSCOMP(IBEG_RHSCOMP-1+(I-IBEG_GLOB_DEF)*LD_RHSCOMP 
     &          + JJ) =  id%DKEEP(2)
           ELSE
            id%RHSCOMP(IBEG_RHSCOMP-1+(I-IBEG_GLOB_DEF)*LD_RHSCOMP 
     &                + JJ) =  ONE
           ENDIF
          ENDIF
         ENDDO
        ENDIF   
        IF ( KEEP(17).NE.0 .and. id%MYID_NODES.EQ.MASTER_ROOT) THEN
            IBEG_ROOT_DEF  = max(IBEG_GLOB_DEF,KEEP(112)+1)
            IEND_ROOT_DEF  = min(IEND_GLOB_DEF,KEEP(112)+KEEP(17))
            IROOT_DEF_RHS_COL1 = IBEG_ROOT_DEF-IBEG_GLOB_DEF + 1
            IBEG_ROOT_DEF = IBEG_ROOT_DEF-KEEP(112)
            IEND_ROOT_DEF = IEND_ROOT_DEF-KEEP(112)
        ELSE
            IBEG_ROOT_DEF = -90999
            IEND_ROOT_DEF = -90999
        ENDIF
       ELSE  
       ENDIF
      ENDIF  
      IF ( I_AM_SLAVE ) THEN
        LWCB_SOL_C = LWCB
        IF ( id%MYID_NODES .EQ. MASTER_ROOT ) THEN
          IF ( associated(id%root%RHS_CNTR_MASTER_ROOT) ) THEN
            PTR_RHS_ROOT => id%root%RHS_CNTR_MASTER_ROOT
            LPTR_RHS_ROOT = size(id%root%RHS_CNTR_MASTER_ROOT)
          ELSE
            LPTR_RHS_ROOT = NBRHS_EFF * SIZE_ROOT
            IPT_RHS_ROOT  = LWCB - LPTR_RHS_ROOT + 1
            PTR_RHS_ROOT => WORK_WCB(IPT_RHS_ROOT:LWCB)
            LWCB_SOL_C = LWCB_SOL_C - LPTR_RHS_ROOT
          ENDIF
        ELSE
          LPTR_RHS_ROOT = 1
          IPT_RHS_ROOT = LWCB 
          PTR_RHS_ROOT => WORK_WCB(IPT_RHS_ROOT:LWCB)
          LWCB_SOL_C = LWCB_SOL_C - LPTR_RHS_ROOT
        ENDIF
      ENDIF
      IF (KEEP(221) .EQ. 2 ) THEN
         IF ( ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) .AND.
     &        ( id%MYID .EQ. MASTER ) ) THEN
            II = 0
            DO K=1, NBRHS_EFF
             KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS -1
             DO I = 1, SIZE_ROOT
              PTR_RHS_ROOT(II+I) = id%REDRHS(KDEC+I)
             ENDDO
             II = II+SIZE_ROOT
            ENDDO
         ELSE
          IF ( id%MYID .EQ. MASTER) THEN
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               KDEC = IBEG_REDRHS
               CALL MPI_SEND(id%REDRHS(KDEC),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_REAL,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
            ELSE
              DO K=1, NBRHS_EFF
                KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS
                CALL MPI_SEND(id%REDRHS(KDEC),SIZE_ROOT,
     &              MPI_REAL,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
              ENDDO
            ENDIF
          ELSE IF ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) THEN
            II = 1
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               CALL MPI_RECV(PTR_RHS_ROOT(II),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_REAL,
     &              MASTER, 0, id%COMM,STATUS,IERR)
            ELSE
             DO K=1, NBRHS_EFF
              CALL MPI_RECV(PTR_RHS_ROOT(II),SIZE_ROOT,
     &           MPI_REAL,
     &           MASTER, 0, id%COMM,STATUS,IERR)
              II = II + SIZE_ROOT
             ENDDO
            ENDIF
          ENDIF
         ENDIF
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        LIW_PASSED = max( LIW, 1 )
        LA_PASSED  = max( LA, 1_8 )
        IF ((id%KEEP(235).EQ.0).and.(id%KEEP(237).EQ.0) ) THEN
        PRUNED_SIZE_LOADED = 0_8  
        CALL SMUMPS_245(id%root, id%N, id%S(1), LA_PASSED,
     &    IS(1), LIW_PASSED, WORK_WCB(1), LWCB_SOL_C,
     &    IWCB, LIWCB, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF, 
     &    id%NA(1),id%LNA,id%NE_STEPS(1), SRW3, MTYPE,
     &    ICNTL(1), id%STEP(1), id%FRERE_STEPS(1), 
     &    id%DAD_STEPS(1), id%FILS(1), id%PTLUST_S(1), id%PTRFAC(1),
     &    IWK_SOLVE, LIWK_SOLVE, id%PROCNODE_STEPS(1),
     &    id%NSLAVES, INFO(1), KEEP(1), KEEP8(1),
     &    id%COMM_NODES, id%MYID, id%MYID_NODES,
     &    id%BUFR(1), id%LBUFR, id%LBUFR_BYTES, 
     &    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     &    IBEG_ROOT_DEF, IEND_ROOT_DEF,
     &    IROOT_DEF_RHS_COL1,
     &    PTR_RHS_ROOT(1), LPTR_RHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     &    id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, 
     &    id%POSINRHSCOMP(1), BUILD_POSINRHSCOMP
     &    , 1            , 1       ,  1       
     &    , 1
     &    , IDUMMY, 1, JDUMMY, KDUMMY, 1, LDUMMY, 1, MDUMMY
     &    )     
        ELSE
          IF ((KEEP(235).NE.0).AND.KEEP(221).NE.2.AND.
     &          KEEP(111).EQ.0) THEN 
              DO K=1, NBRHS_EFF
               DO I=1, LD_RHSCOMP
                id%RHSCOMP(IBEG_RHSCOMP+(K-1)*LD_RHSCOMP+I-1)= ZERO
               ENDDO
              ENDDO
          ELSEIF  (KEEP(237).NE.0) THEN
              DO K=1, NBRHS_EFF
               DO I=1, LD_RHSCOMP
                id%RHSCOMP(IBEG_RHSCOMP+(K-1)*LD_RHSCOMP+I-1)= ZERO
               ENDDO
              ENDDO
          ENDIF
           IF (.NOT. allocated(PERM_RHS)) THEN
              ALLOCATE(PERM_RHS(1),stat=allocok)
              NB_BYTES = NB_BYTES +  int(size(PERM_RHS),8)*K34_8
              NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
           ENDIF
           CALL SMUMPS_245(id%root, id%N, id%S(1), LA_PASSED,
     &    IS(1), LIW_PASSED, WORK_WCB(1), LWCB_SOL_C, IWCB, LIWCB,
     &    RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF, id%NA(1),id%LNA,
     &    id%NE_STEPS(1), SRW3, MTYPE, ICNTL(1), id%STEP(1),
     &    id%FRERE_STEPS(1), id%DAD_STEPS(1), id%FILS(1),
     &    id%PTLUST_S(1), id%PTRFAC(1), IWK_SOLVE, LIWK_SOLVE,
     &    id%PROCNODE_STEPS(1), id%NSLAVES, INFO(1),KEEP(1),KEEP8(1),
     &    id%COMM_NODES, id%MYID, id%MYID_NODES,
     &    id%BUFR(1), id%LBUFR, id%LBUFR_BYTES, 
     &    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     &    IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1,
     &    PTR_RHS_ROOT(1), LPTR_RHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     &    id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, 
     &    id%POSINRHSCOMP(1), BUILD_POSINRHSCOMP,
     &    NZ_THIS_BLOCK, NBCOL_INBLOC, id%NRHS,
     &    JBEG_RHS , id%Step2node(1), id%KEEP(28), IRHS_SPARSE_COPY(1),
     &    IRHS_PTR_COPY(1),
     &    size(PERM_RHS), PERM_RHS, size(UNS_PERM_INV), UNS_PERM_INV 
     &       )
        ENDIF   
      END IF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF (INFO(1).eq.-2) then
        INFO(1)=-11
        IF (LP.GT.0) 
     &  write(LP,*) 
     &   ' WARNING : -11 error code obtained in solve'
      END IF
      IF (INFO(1).eq.-3) then
        INFO(1)=-14
        IF (LP.GT.0) 
     &  write(LP,*) 
     &    ' WARNING : -14 error code obtained in solve'
      END IF
      IF (INFO(1).LT.0) GO TO 90
      IF ( KEEP(221) .EQ. 1 ) THEN 
         IF ( ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) .AND.
     &        ( id%MYID .EQ. MASTER ) ) THEN
            II = 0
            DO K=1, NBRHS_EFF
             KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS -1
             DO I = 1, SIZE_ROOT
              id%REDRHS(KDEC+I) = PTR_RHS_ROOT(II+I) 
             ENDDO
             II = II+SIZE_ROOT
            ENDDO
         ELSE
          IF ( id%MYID .EQ. MASTER ) THEN
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               KDEC = IBEG_REDRHS
               CALL MPI_RECV(id%REDRHS(KDEC),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_REAL,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,
     &              STATUS,IERR)
            ELSE
             DO K=1, NBRHS_EFF
               KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS
               CALL MPI_RECV(id%REDRHS(KDEC),SIZE_ROOT,
     &              MPI_REAL,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,
     &              STATUS,IERR)
             ENDDO
            ENDIF
          ELSE IF ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) THEN
            II = 1
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               CALL MPI_SEND(PTR_RHS_ROOT(II),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_REAL,
     &              MASTER, 0, id%COMM,IERR)
            ELSE
             DO K=1, NBRHS_EFF
              CALL MPI_SEND(PTR_RHS_ROOT(II),SIZE_ROOT,
     &           MPI_REAL,
     &           MASTER, 0, id%COMM,IERR)
              II = II + SIZE_ROOT
             ENDDO
            ENDIF
          ENDIF
         ENDIF
      ENDIF 
      IF ( KEEP(221) .NE. 1 ) THEN  
       IF (ICNTL21 == 0) THEN 
        IF ((.NOT. I_AM_SLAVE).AND.KEEP(237).EQ.0) THEN
          ALLOCATE( CWORK(KEEP(247)*NBRHS_EFF), stat=allocok)
          IF (allocok > 0) THEN
            ALLOCATE( CWORK(KEEP(247)), stat=allocok)
            IF (allocok > 0) THEN
              INFO(1)=-13
              INFO(2)=KEEP(247)
            ENDIF
          ENDIF
        ENDIF
        CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
        IF (INFO(1).LT.0) GO TO 90
        IF ((id%MYID.NE.MASTER).OR. .NOT.LSCAL) THEN
         PT_SCALING => Dummy_SCAL
        ELSE
          IF (MTYPE.EQ.1) THEN
            PT_SCALING => id%COLSCA
          ELSE
            PT_SCALING => id%ROWSCA
          ENDIF
        ENDIF
        LIW_PASSED = max( LIW, 1 )
        IF ( .NOT.I_AM_SLAVE ) THEN
          IF (KEEP(237).EQ.0) THEN
           CALL SMUMPS_521(id%NSLAVES,id%N,
     &          id%MYID, id%COMM, 
     &          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     &          JDUMMY, id%KEEP(1), id%KEEP8(1), 
     &          id%PROCNODE_STEPS(1), IDUMMY, 1,
     &          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES, 
     &          CWORK(1), size(CWORK), 
     &          LSCAL, PT_SCALING(1), size(PT_SCALING)
     &          )
            DEALLOCATE( CWORK )
          ELSE
           CALL SMUMPS_812(id%NSLAVES,id%N,
     &          id%MYID, id%COMM, 
     &          RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     &          id%KEEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES, 
     &          LSCAL, PT_SCALING(1), size(PT_SCALING)
     &         ,IRHS_PTR_COPY(1), size(IRHS_PTR_COPY), 
     &          IRHS_SPARSE_COPY(1), size(IRHS_SPARSE_COPY),
     &          RHS_SPARSE_COPY(1), size(RHS_SPARSE_COPY),
     &          UNS_PERM_INV, size(UNS_PERM_INV), IDUMMY, 1
     &          )
          ENDIF
        ELSE
          IF (KEEP(237).EQ.0) THEN
           CALL SMUMPS_521(id%NSLAVES,id%N, 
     &          id%MYID, id%COMM,
     &          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     &          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     &          id%PROCNODE_STEPS(1), IS(1), LIW_PASSED,
     &          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &          id%RHSCOMP(1), LENRHSCOMP,
     &          LSCAL, PT_SCALING(1), size(PT_SCALING)
     &         )
           ELSE 
            CALL SMUMPS_812(id%NSLAVES,id%N, 
     &          id%MYID, id%COMM,
     &          RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     &          id%KEEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &          LSCAL, PT_SCALING(1), size(PT_SCALING)
     &         , IRHS_PTR_COPY(1), size(IRHS_PTR_COPY), 
     &          IRHS_SPARSE_COPY(1), size(IRHS_SPARSE_COPY),
     &          RHS_SPARSE_COPY(1), size(RHS_SPARSE_COPY),
     &          UNS_PERM_INV, size(UNS_PERM_INV), POSINRHSCOMP_N, 
     &          id%N
     &          )
           ENDIF
        ENDIF
        IF ( (id%MYID.EQ.MASTER).AND. (KEEP(237).NE.0) 
     &      ) THEN
          IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
            DO J = JBEG_RHS, JBEG_RHS+NBCOL_INBLOC-1
              COLSIZE = id%IRHS_PTR(PERM_RHS(J)+1) - 
     &                  id%IRHS_PTR(PERM_RHS(J))
              IF (COLSIZE.EQ.0) CYCLE
              JJ = J-JBEG_RHS+1
              DO IZ= id%IRHS_PTR(PERM_RHS(J)), 
     &               id%IRHS_PTR(PERM_RHS(J)+1)-1
               I = id%IRHS_SPARSE (IZ)
               DO IZ2 = IRHS_PTR_COPY(JJ),IRHS_PTR_COPY(JJ+1)-1
                 IF (IRHS_SPARSE_COPY(IZ2).EQ.I) EXIT
                 IF (IZ2.EQ.IRHS_PTR_COPY(JJ+1)-1) THEN
                   WRITE(6,*) " INTERNAL ERROR 1 gather sol_driver"
                   CALL MUMPS_ABORT()
                 ENDIF
               ENDDO
               id%RHS_SPARSE(IZ) = RHS_SPARSE_COPY(IZ2)
              ENDDO
            ENDDO
          ELSE            
            DO J = JBEG_RHS, JBEG_RHS+NBCOL_INBLOC-1
              COLSIZE = id%IRHS_PTR(J+1) - id%IRHS_PTR(J)
              IF (COLSIZE.EQ.0) CYCLE
              JJ = J-JBEG_RHS+1
              DO IZ= id%IRHS_PTR(J), id%IRHS_PTR(J+1) - 1
               I = id%IRHS_SPARSE (IZ)
               DO IZ2 = IRHS_PTR_COPY(JJ),IRHS_PTR_COPY(JJ+1)-1
                 IF (IRHS_SPARSE_COPY(IZ2).EQ.I) EXIT
                 IF (IZ2.EQ.IRHS_PTR_COPY(JJ+1)-1) THEN
                   WRITE(6,*) " INTERNAL ERROR 2 gather sol_driver"
                   CALL MUMPS_ABORT()
                 ENDIF
               ENDDO
               id%RHS_SPARSE(IZ) = RHS_SPARSE_COPY(IZ2)
              ENDDO
            ENDDO
          ENDIF           
        ENDIF             
       ELSE   
        IF ( I_AM_SLAVE ) THEN
         LIW_PASSED = max( LIW, 1 )
         IF ( KEEP(89) .GT. 0 ) THEN
           CALL SMUMPS_532(id%NSLAVES,
     &          id%N, id%MYID_NODES,
     &          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     &          id%ISOL_loc(1),
     &          id%SOL_loc(1), JBEG_RHS-NB_RHSSKIPPED, id%LSOL_loc,
     &          id%PTLUST_S(1), id%PROCNODE_STEPS(1),
     &          id%KEEP(1),id%KEEP8(1),
     &          IS(1), LIW_PASSED,
     &          id%STEP(1), scaling_data, LSCAL, NB_RHSSKIPPED )
         ENDIF
        ENDIF
       ENDIF  
      ENDIF 
      IF ( ICNTL10 > 0 .AND. NBRHS_EFF > 1 ) THEN
        DO I = 1, ICNTL10
          write(*,*) 'FIXME: to be implemented'
        END DO
      END IF
      IF (ERANAL) THEN
        IF ((ICNTL10 .GT. 0) .AND. (ICNTL11 .GT. 0)) THEN
          IF (id%MYID .EQ. MASTER) THEN
            GIVSOL = .FALSE.
            IF (MP .GT. 0) WRITE( MP, 170 )
            ALLOCATE(R_RW1(id%N),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 776
            ENDIF
            ALLOCATE(C_RW2(id%N),stat=allocok)
            IF (allocok .GT.0) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 776
            ENDIF
            NB_BYTES = NB_BYTES + int(id%N,8)*K35_8 + int(id%N,8)*K16_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
          END IF
 776      CONTINUE
          CALL MUMPS_276( ICNTL(1), INFO(1),
     &                  id%COMM,id%MYID)
          IF ( INFO(1) .LT. 0 ) GOTO 90
          IF ( KEEP(54) .eq. 0 ) THEN
            IF (id%MYID .EQ. MASTER) THEN
              IF (KEEP(55).EQ.0) THEN
                CALL SMUMPS_278( ICNTL(9), id%N, id%NZ, id%A(1),
     &             id%IRN(1), id%JCN(1),
     &             RHS_MUMPS(IBEG), SAVERHS, R_RW1, C_RW2,
     &             KEEP(1),KEEP8(1) )
              ELSE
                CALL SMUMPS_121( ICNTL(9), id%N, 
     &          id%NELT, id%ELTPTR(1), 
     &          id%LELTVAR, id%ELTVAR(1),
     &          id%NA_ELT, id%A_ELT(1),
     &          RHS_MUMPS(IBEG), SAVERHS, R_RW1, C_RW2,
     &          KEEP(1),KEEP8(1) )
              ENDIF
            END IF
          ELSE
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     &              MPI_REAL, MASTER,
     &              id%COMM, IERR )
            ALLOCATE( C_LOCWK54( id%N ), stat =allocok )
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=id%N
            endif
            CALL MUMPS_276(ICNTL(1), INFO(1), id%COMM, id%MYID)
            IF ( INFO(1) .LT. 0 ) GOTO 90
            NB_BYTES     = NB_BYTES + int(id%N,8)*K35_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
            IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
              CALL SMUMPS_192( id%N, id%NZ_loc,
     &        id%IRN_loc(1), id%JCN_loc(1), id%A_loc(1),
     &        RHS_MUMPS(IBEG), C_LOCWK54, KEEP(50), MTYPE )
            ELSE
              C_LOCWK54 = ZERO
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( C_LOCWK54, C_RW2,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
              C_RW2 = SAVERHS - C_RW2
            ELSE
              CALL MPI_REDUCE( C_LOCWK54, C_DUMMY,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            NB_BYTES = NB_BYTES - int(size(C_LOCWK54),8)*K35_8
            DEALLOCATE( C_LOCWK54 )
            ALLOCATE( R_LOCWK54( id%N ), stat =allocok )
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=id%N
            endif
            CALL MUMPS_276(ICNTL(1), INFO(1), id%COMM, id%MYID)
            IF ( INFO(1) .LT. 0 ) GOTO 90
            NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
            NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
            IF ( I_AM_SLAVE .AND. id%NZ_loc .NE. 0 ) THEN
              CALL SMUMPS_207(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%IRN_loc(1), id%JCN_loc(1),
     &          R_LOCWK54, id%KEEP(1),id%KEEP8(1))
            ELSE
              R_LOCWK54 = RZERO
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( R_LOCWK54, R_RW1,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( R_LOCWK54, R_DUMMY,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            NB_BYTES = NB_BYTES - int(size(R_LOCWK54),8)*K16_8
            DEALLOCATE( R_LOCWK54 )
          END IF
          IF ( id%MYID .EQ. MASTER )  THEN
            CALL SMUMPS_205(ICNTL(9),INFO(1),id%N,id%NZ,
     &        RHS_MUMPS(IBEG), SAVERHS,R_RW1,C_RW2,GIVSOL,
     &        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MP,ICNTL(1),
     &        KEEP(1),KEEP8(1))
            NB_BYTES = NB_BYTES - int(size(R_RW1),8)*K16_8 
     &                          - int(size(C_RW2),8)*K35_8
            DEALLOCATE(R_RW1)
            DEALLOCATE(C_RW2)
          END IF
        END IF
      IF ( PROK  .AND. ICNTL10 .GT. 0 ) WRITE( MP, 270 )
      IF ( PROKG .AND. ICNTL10 .GT. 0 ) WRITE( MPG, 270 )
      ALLOCATE(R_Y(id%N), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=id%N
        GOTO 777
      ENDIF
      NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
      ALLOCATE(C_Y(id%N), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=id%N
        GOTO 777
      ENDIF
      NB_BYTES = NB_BYTES + int(id%N,8)*K35_8
      IF ( id%MYID .EQ. MASTER ) THEN
        ALLOCATE( IW1( 2 * id%N ),stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=2 * id%N
          GOTO 777
        ENDIF
        NB_BYTES = NB_BYTES + int(2*id%N,8)*K34_8
        ALLOCATE( D(id%N),stat =allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=id%N
          GOTO 777
        ENDIF
        NB_BYTES = NB_BYTES + int(id%N,8)*K16_8 
        ALLOCATE( C_W(id%N), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=id%N
          GOTO 777
        ENDIF
        NB_BYTES = NB_BYTES + int(id%N,8)*K35_8
        ALLOCATE( R_W(2*id%N), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=id%N
          GOTO 777
        ENDIF
        NB_BYTES = NB_BYTES + int(2*id%N,8)*K16_8
        NITREF = ICNTL10
        JOBIREF= ICNTL11
        IF ( PROKG .AND. ICNTL10 .GT. 0 )
     &    WRITE( MPG, 240) 'MAXIMUM NUMBER OF STEPS =', NITREF
        DO I = 1, id%N
          D( I ) = RONE
        END DO
      END IF
      ALLOCATE(C_LOCWK54(id%N),stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=id%N
        GOTO 777
      ENDIF
      NB_BYTES = NB_BYTES + int(id%N,8)*K35_8
      ALLOCATE(R_LOCWK54(id%N),stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=id%N
        GOTO 777
      ENDIF
      NB_BYTES = NB_BYTES + int(id%N,8)*K16_8
      KASE = 0
 777  CONTINUE
      NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
      IF ( INFO(1) .LT. 0 ) GOTO 90
  22    CONTINUE
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 0 ) THEN
              IF (KEEP(55).NE.0) THEN 
               CALL SMUMPS_119(MTYPE, id%N, 
     &           id%NELT, id%ELTPTR(1), 
     &           id%LELTVAR, id%ELTVAR(1),
     &           id%NA_ELT, id%A_ELT(1),
     &           R_W(id%N+1), KEEP(1),KEEP8(1) )
              ELSE
               IF ( MTYPE .eq. 1 ) THEN
                 CALL SMUMPS_207
     &   ( id%A(1), id%NZ, id%N, id%IRN(1), id%JCN(1),
     &     R_W(id%N+1), KEEP(1),KEEP8(1))
               ELSE
                 CALL SMUMPS_207
     &   ( id%A(1), id%NZ, id%N, id%JCN(1), id%IRN(1),
     &     R_W(id%N+1), KEEP(1),KEEP8(1))
               END IF
              ENDIF
            ENDIF
          END IF
        ELSE
          IF ( KASE .eq. 0 ) THEN
            IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
              IF ( MTYPE .eq. 1 ) THEN
              CALL SMUMPS_207(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%IRN_loc(1), id%JCN_loc(1),
     &          R_LOCWK54, id%KEEP(1),id%KEEP8(1) )
              ELSE
              CALL SMUMPS_207(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%JCN_loc(1), id%IRN_loc(1),
     &          R_LOCWK54, id%KEEP(1),id%KEEP8(1) )
              END IF
            ELSE
              R_LOCWK54 = RZERO
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( R_LOCWK54, R_W( id%N + 1 ),
     &          id%N, MPI_REAL,
     &          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( R_LOCWK54, R_DUMMY, 
     &          id%N, MPI_REAL,
     &          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
          END IF
        END IF
        IF ( id%MYID .eq. MASTER ) THEN
            ARRET = CNTL(2)
            IF (ARRET .LT. 0.0E0) THEN
              ARRET = sqrt(epsilon(0.0E0))
            END IF
            CALL SMUMPS_206(id%NZ,id%N,SAVERHS,RHS_MUMPS(IBEG),
     &      C_Y, D, R_W, C_W,
     &      IW1, KASE,RINFOG(7),
     &      RINFOG(9), JOBIREF, RINFOG(10), NITREF, NOITER, MP,
     &      KEEP(1),KEEP8(1), ARRET )
        END IF
        IF ( KEEP(54) .ne. 0 ) THEN
          CALL MPI_BCAST( KASE, 1, MPI_INTEGER, MASTER,
     &    id%COMM, IERR )
        END IF
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 14 ) THEN
              IF (KEEP(55).NE.0) THEN
               CALL SMUMPS_122( MTYPE, id%N, 
     &            id%NELT, id%ELTPTR(1), id%LELTVAR,
     &            id%ELTVAR(1), id%NA_ELT, id%A_ELT(1),
     &            SAVERHS, RHS_MUMPS(IBEG),
     &            C_Y, R_W, KEEP(50))
              ELSE
                 IF ( MTYPE .eq. 1 ) THEN
                   CALL SMUMPS_208
     &    (id%A(1), id%NZ, id%N, id%IRN(1), id%JCN(1), SAVERHS,
     &    RHS_MUMPS(IBEG), C_Y, R_W, KEEP(1),KEEP8(1))
                 ELSE
                   CALL SMUMPS_208
     &    (id%A(1), id%NZ, id%N, id%JCN(1), id%IRN(1), SAVERHS,
     &    RHS_MUMPS(IBEG), C_Y, R_W, KEEP(1),KEEP8(1))
                 END IF
              ENDIF
              GOTO 22
            END IF
          END IF
        ELSE
          IF ( KASE.eq.14 ) THEN
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     &              MPI_REAL, MASTER,
     &              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
              CALL SMUMPS_192( id%N, id%NZ_loc,
     &        id%IRN_loc(1), id%JCN_loc(1), id%A_loc(1),
     &        RHS_MUMPS(IBEG), C_LOCWK54, KEEP(50), MTYPE )
            ELSE
              C_LOCWK54 = ZERO
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( C_LOCWK54, C_Y,
     &          id%N, MPI_REAL,
     &          MPI_SUM,MASTER,id%COMM, IERR)
              C_Y = SAVERHS - C_Y
            ELSE
              CALL MPI_REDUCE( C_LOCWK54, C_DUMMY, 
     &          id%N, MPI_REAL,
     &          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and. id%NZ_loc .NE. 0 ) THEN
              CALL SMUMPS_193( id%N, id%NZ_loc,
     &        id%IRN_loc(1), id%JCN_loc(1), id%A_loc(1),
     &        RHS_MUMPS(IBEG), R_LOCWK54, KEEP(50), MTYPE )
            ELSE
              R_LOCWK54 = RZERO
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( R_LOCWK54, R_W,
     &          id%N, MPI_REAL,
     &          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( R_LOCWK54, R_DUMMY, 
     &          id%N, MPI_REAL,
     &          MPI_SUM, MASTER, id%COMM, IERR)
            END IF
            GOTO 22
          END IF
        END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF ( KASE .GT. 0 ) THEN
          IF ( MTYPE .EQ. 1 ) THEN
            SOLVET = KASE - 1
          ELSE
            SOLVET = KASE
          END IF
          IF ( LSCAL ) THEN
            IF ( SOLVET .EQ. 1 ) THEN
              DO K = 1, id%N
                C_Y( K ) = C_Y( K ) * id%ROWSCA( K )
              END DO
            ELSE
              DO K = 1, id%N
                C_Y( K ) = C_Y( K ) * id%COLSCA( K )
              END DO
            END IF
          END IF
        END IF
      END IF
      CALL MPI_BCAST( KASE , 1, MPI_INTEGER, MASTER,
     &                id%COMM, IERR)
      CALL MPI_BCAST( SOLVET, 1, MPI_INTEGER, MASTER,
     &                id%COMM, IERR)
      IF ( KASE .GT. 0 ) THEN
          BUILD_POSINRHSCOMP=.FALSE.  
          IF ( .NOT.I_AM_SLAVE ) THEN
            CALL SMUMPS_638(id%NSLAVES,id%N, id%MYID, id%COMM,
     &          MTYPE, C_Y(1), id%N, 1,
     &          JDUMMY, id%KEEP(1), id%KEEP8(1), id%PROCNODE_STEPS(1),
     &          IDUMMY, 1,
     &          id%STEP(1), KDUMMY, 1, BUILD_POSINRHSCOMP,
     &          id%ICNTL(1),id%INFO(1))
          ELSE
            LIW_PASSED = max( LIW, 1 )
            CALL SMUMPS_638(id%NSLAVES,id%N, id%MYID, id%COMM,
     &          MTYPE, C_Y(1), id%N, 1,
     &          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     &          id%PROCNODE_STEPS(1),
     &          IS(1), LIW_PASSED,
     &          id%STEP(1), id%POSINRHSCOMP(1), KEEP(28), 
     &          BUILD_POSINRHSCOMP,
     &          id%ICNTL(1),id%INFO(1))
          ENDIF
          IF  (INFO(1).LT.0) GOTO 89
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED = max( LIW, 1 )
          LA_PASSED = max( LA, 1_8 )
          CALL SMUMPS_245( id%root, id%N, id%S(1), LA_PASSED,
     & id%IS(1), LIW_PASSED, WORK_WCB(1), LWCB_SOL_C, IWCB, LIWCB, C_Y,
     & id%N, NBRHS_EFF, id%NA(1), id%LNA, id%NE_STEPS(1), SRW3, SOLVET,
     & ICNTL(1), id%STEP(1), id%FRERE_STEPS(1), id%DAD_STEPS(1), id%
     & FILS(1), id%PTLUST_S(1), id%PTRFAC(1), IWK_SOLVE(1), LIWK_SOLVE,
     & id%PROCNODE_STEPS(1), id%NSLAVES, INFO(1), KEEP(1), KEEP8(1),
     & id%COMM_NODES, id%MYID, id%MYID_NODES, id%BUFR(1), id%LBUFR,
     & id%LBUFR_BYTES, id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     & IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1,
     &    PTR_RHS_ROOT(1), LPTR_RHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     &    id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, 
     &    id%POSINRHSCOMP(1), BUILD_POSINRHSCOMP 
     &    , 1            , 1       ,  1       
     &    , 1
     &    , IDUMMY, 1, JDUMMY, KDUMMY, 1, LDUMMY, 1, MDUMMY
     &    )
        END IF
        IF (INFO(1).eq.-2) INFO(1)=-12
        IF (INFO(1).eq.-3) INFO(1)=-15
        IF ( .NOT. I_AM_SLAVE .AND. INFO(1) .GE. 0) THEN
          ALLOCATE( CWORK(KEEP(247)*NBRHS_EFF), stat=allocok)
          IF (allocok > 0) THEN
            ALLOCATE( CWORK(KEEP(247)), stat=allocok)
            IF (allocok > 0) THEN
              INFO(1)=-13
              INFO(2)=KEEP(247)
            ENDIF
          ENDIF
        ENDIF
 89     CALL MUMPS_276( ICNTL(1), INFO(1),
     &                   id%COMM,id%MYID)
        IF (INFO(1).LT.0) GO TO 90
        IF ((id%MYID.NE.MASTER).OR. .NOT.LSCAL) THEN
         PT_SCALING => Dummy_SCAL
        ELSE
          IF (SOLVET.EQ.1) THEN
            PT_SCALING => id%COLSCA
          ELSE
            PT_SCALING => id%ROWSCA
          ENDIF
        ENDIF
        LIW_PASSED = max( LIW, 1 )
        IF ( .NOT. I_AM_SLAVE ) THEN
          CALL SMUMPS_521(id%NSLAVES,id%N,
     &          id%MYID, id%COMM,
     &          SOLVET, C_Y, id%N, NBRHS_EFF,
     &          JDUMMY, id%KEEP(1),id%KEEP8(1), id%PROCNODE_STEPS(1),
     &          IDUMMY, 1,
     &          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &          CWORK, size(CWORK),
     &          LSCAL, PT_SCALING(1), size(PT_SCALING))
          DEALLOCATE( CWORK )
        ELSE
          CALL SMUMPS_521(id%NSLAVES,id%N,
     &          id%MYID, id%COMM,
     &          SOLVET, C_Y, id%N, NBRHS_EFF,
     &          id%PTLUST_S(1), id%KEEP(1),id%KEEP8(1),
     &          id%PROCNODE_STEPS(1),
     &          IS(1), LIW_PASSED,
     &          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES,
     &          id%RHSCOMP(1), LENRHSCOMP,
     &          LSCAL, PT_SCALING(1), size(PT_SCALING))
        ENDIF
        GO TO 22
      ELSEIF ( KASE .LT. 0 ) THEN
           INFO( 1 ) = INFO( 1 ) + 8
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
         NB_BYTES = NB_BYTES - int(size(R_W),8)*K16_8
     &                       - int(size(D  ),8)*K16_8
     &                       - int(size(IW1),8)*K34_8
        DEALLOCATE(R_W,D)
        DEALLOCATE(IW1)
      ENDIF
      IF ( PROKG ) THEN
        IF (NITREF.GT.0) THEN
        WRITE( MPG, 81 ) 
        WRITE( MPG, 141 ) 'NUMBER OF STEPS OF ITERATIVE REFINEMENTS  
     &=', NOITER
       ENDIF
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
       IF ( NITREF .GT. 0 ) THEN
        id%INFOG(15) = NOITER
       END IF
      END IF
      IF ( PROK .AND. ICNTL10 .GT.0 ) WRITE( MP, 131 ) 
      IF (ICNTL11 .GT. 0) THEN
        IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
              CALL SMUMPS_278( MTYPE, id%N, id%NZ, id%A(1),
     &          id%IRN(1), id%JCN(1),
     &          RHS_MUMPS(IBEG), SAVERHS, R_Y, C_W, KEEP(1),KEEP8(1))
            ELSE
              CALL SMUMPS_121( MTYPE, id%N, 
     &          id%NELT, id%ELTPTR(1),
     &          id%LELTVAR, id%ELTVAR(1),
     &          id%NA_ELT, id%A_ELT(1),
     &          RHS_MUMPS(IBEG), SAVERHS, R_Y, C_W, KEEP(1),KEEP8(1))
            ENDIF
          END IF
        ELSE
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     &              MPI_REAL, MASTER, 
     &              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
              CALL SMUMPS_192( id%N, id%NZ_loc,
     &        id%IRN_loc(1), id%JCN_loc(1), id%A_loc(1),
     &        RHS_MUMPS(IBEG), C_LOCWK54, KEEP(50), MTYPE ) 
            ELSE
              C_LOCWK54 = ZERO
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( C_LOCWK54, C_W,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
              C_W = SAVERHS - C_W
            ELSE
              CALL MPI_REDUCE( C_LOCWK54, C_DUMMY,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
              CALL SMUMPS_207(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%IRN_loc(1), id%JCN_loc(1),
     &          R_LOCWK54, id%KEEP(1),id%KEEP8(1) )
            ELSE
              R_LOCWK54 = RZERO
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( R_LOCWK54, R_Y,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( R_LOCWK54, R_DUMMY,
     &        id%N, MPI_REAL,
     &        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
        END IF
        IF (id%MYID .EQ. MASTER) THEN
         IF ((MPG .GT. 0) .AND. (NITREF .GT. 0)) WRITE( MPG, 65 )
         IF ((MPG .GT. 0) .AND. (NITREF .LE. 0)) WRITE( MPG, 170 )
         GIVSOL = .FALSE.
         CALL SMUMPS_205(MTYPE,INFO(1),id%N,id%NZ,RHS_MUMPS(IBEG),
     &        SAVERHS,R_Y,C_W,GIVSOL,
     &        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MPG,ICNTL(1),
     &        KEEP(1),KEEP8(1))
         IF ( MPG .GT. 0 ) THEN
          WRITE( MPG, 115 )
     &'RINFOG(7):COMPONENTWISE SCALED RESIDUAL(W1)=', RINFOG(7)
          WRITE( MPG, 115 )
     &'------(8):---------------------------- (W2)=', RINFOG(8)
          WRITE( MPG, 115 )
     &'------(9):Upper bound ERROR ...............=', RINFOG(9)
          WRITE( MPG, 115 )
     &'-----(10):CONDITION NUMBER (1) ............=', RINFOG(10)
          WRITE( MPG, 115 )
     &'-----(11):CONDITION NUMBER (2) ............=', RINFOG(11)
         END IF
        END IF 
      END IF 
      IF (id%MYID == MASTER) THEN 
         NB_BYTES = NB_BYTES - int(size(C_W),8)*K35_8
         DEALLOCATE(C_W)
      ENDIF
      NB_BYTES = NB_BYTES - 
     &   (int(size(R_Y),8)+int(size(R_LOCWK54),8))*K16_8
      NB_BYTES = NB_BYTES - 
     &   (int(size(C_Y),8)+int(size(C_LOCWK54),8))*K35_8
      DEALLOCATE(R_Y)
      DEALLOCATE(C_Y)
      DEALLOCATE(R_LOCWK54)
      DEALLOCATE(C_LOCWK54)
      END IF
      IF ( id%MYID .EQ. MASTER .AND. ICNTL21==0
     &     .AND. KEEP(23) .NE. 0 .AND. KEEP(237).EQ.0 ) THEN
        IF ((KEEP(221).NE.1 .AND. ICNTL(9) .EQ. 1)
     &     .OR. KEEP(111) .NE.0 .OR. KEEP(252).NE.0 ) THEN
          ALLOCATE( C_RW1( id%N ),stat =allocok ) 
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
            WRITE(*,*) 'could not allocate ', id%N, 'integers.'
            CALL MUMPS_ABORT()
          END IF
          DO K = 1, NBRHS_EFF
           KDEC = (K-1)*LD_RHS+IBEG-1
           DO 70 I = 1, id%N
            C_RW1(I) = RHS_MUMPS(KDEC+I)
 70        CONTINUE
           DO 80 I = 1, id%N
            JPERM = id%UNS_PERM(I)
            RHS_MUMPS( KDEC+JPERM ) = C_RW1( I )
 80        CONTINUE
          END DO
          DEALLOCATE( C_RW1 ) 
        END IF
      END IF
      IF (id%MYID.EQ.MASTER .and.ICNTL21==0.and.KEEP(221).NE.1
     &     .and. KEEP(237).EQ.0 ) THEN
        IF ( INFO(1) .GE. 0 .AND. ICNTL(4).GE.3 .AND. ICNTL(3).GT.0)
     &    THEN
          K = min0(10, id%N)
          IF (ICNTL(4) .eq. 4 ) K = id%N
          J = min0(10,NBRHS_EFF)
          IF (ICNTL(4) .eq. 4 ) J = NBRHS_EFF
          DO II=1, J
            WRITE(ICNTL(3),110) BEG_RHS+II-1
            WRITE(ICNTL(3),160)
     &    (RHS_MUMPS(IBEG+(II-1)*LD_RHS+I-1),I=1,K)
          ENDDO
        END IF
      END IF
        IF ((KEEP(248).EQ.1).AND.(KEEP(237).EQ.0)) THEN
          BEG_RHS = BEG_RHS + NBRHS_EFF 
        ELSE 
          BEG_RHS = BEG_RHS + NBRHS
        ENDIF
      ENDDO  
      IF (   (id%MYID.EQ.MASTER) 
     &       .AND. ( KEEP(248).NE.0 )  
     &       .AND. ( KEEP(237).EQ.0 )  
     &       .AND. ( ICNTL21.EQ.0 )    
     &       .AND. ( KEEP(221) .NE.1 ) 
     &       .AND. ( JEND_RHS .LT. id%NRHS )
     &   )
     &         THEN
        JBEG_NEW = JEND_RHS + 1
        IF (DO_PERMUTE_RHS.OR.INTERLEAVE_PAR) THEN
           DO WHILE ( JBEG_NEW.LE. id%NRHS) 
              DO I=1, id%N
                 RHS_MUMPS((PERM_RHS(JBEG_NEW) -1)*LD_RHS+I) 
     &                     = ZERO
              ENDDO
              JBEG_NEW = JBEG_NEW +1
              CYCLE
           ENDDO
        ELSE
          DO WHILE ( JBEG_NEW.LE. id%NRHS) 
            DO I=1, id%N
                RHS_MUMPS((JBEG_NEW -1)*LD_RHS + I) = ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
          ENDDO
        ENDIF                   
      ENDIF
      IF ( I_AM_SLAVE .AND. (ICNTL21.NE.0) .AND. 
     &      ( JEND_RHS .LT. id%NRHS ) .AND. KEEP(221).NE.1 ) THEN
       JBEG_NEW = JEND_RHS + 1
           DO WHILE ( JBEG_NEW.LE. id%NRHS) 
            DO I=1, KEEP(89)
                id%SOL_loc((JBEG_NEW -1)*id%LSOL_loc + I) = ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
           ENDDO
      ENDIF
      IF ((KEEP(221).EQ.1) .AND. 
     &        ( JEND_RHS .LT. id%NRHS ) ) THEN
       IF (id%MYID .EQ. MASTER) THEN
           JBEG_NEW = JEND_RHS + 1
           DO WHILE ( JBEG_NEW.LE. id%NRHS) 
            DO I=1,  id%SIZE_SCHUR 
              id%REDRHS((JBEG_NEW -1)*LD_REDRHS + I) =  ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
           ENDDO
       ENDIF
       IF (I_AM_SLAVE) THEN
           JBEG_NEW = JEND_RHS + 1
           DO WHILE ( JBEG_NEW.LE. id%NRHS) 
            DO I=1,LD_RHSCOMP
              id%RHSCOMP((JBEG_NEW -1)*LD_RHSCOMP + I) =  ZERO
            ENDDO
            JBEG_NEW = JBEG_NEW +1
           ENDDO
       ENDIF
      ENDIF
      id%INFO(26) = int(NB_BYTES_MAX / 1000000_8, 4)
      CALL MUMPS_243( id%MYID, id%COMM,
     &                           id%INFO(26), id%INFOG(30), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Rank of processor needing largest memory in solve     :',
     &  IRANK
        WRITE( MPG,'(A,I10) ')
     &  ' ** Space in MBYTES used by this processor for solve      :',
     &  id%INFOG(30)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during solve    :',
     &  ( id%INFOG(31)-id%INFO(26) ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during solve    :',
     &  id%INFOG(31) / id%NSLAVES
        END IF
      END IF
 90   CONTINUE
      IF (INFO(1) .LT.0 ) THEN 
      ENDIF
      IF (KEEP(201).GT.0)THEN
        IF (IS_INIT_OOC_DONE) THEN
          CALL SMUMPS_582(IERR)
          IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
        ENDIF
        CALL MUMPS_276( ICNTL(1), INFO(1),
     &         id%COMM,id%MYID)
      ENDIF
        IF (IRHS_SPARSE_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES - 
     &        int(size(IRHS_SPARSE_COPY),8)*K34_8 
             DEALLOCATE(IRHS_SPARSE_COPY)
             IRHS_SPARSE_COPY_ALLOCATED=.FALSE.
             NULLIFY(IRHS_SPARSE_COPY)
       ENDIF 
       IF (IRHS_PTR_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES - 
     &        int(size(IRHS_PTR_COPY),8)*K34_8 
             DEALLOCATE(IRHS_PTR_COPY)
             IRHS_PTR_COPY_ALLOCATED=.FALSE.
             NULLIFY(IRHS_PTR_COPY)
       ENDIF
       IF (RHS_SPARSE_COPY_ALLOCATED) THEN
             NB_BYTES =  NB_BYTES - 
     &        int(size(RHS_SPARSE_COPY),8)*K35_8 
             DEALLOCATE(RHS_SPARSE_COPY)
             RHS_SPARSE_COPY_ALLOCATED=.FALSE.
             NULLIFY(RHS_SPARSE_COPY)
       ENDIF 
      IF (allocated(PERM_RHS)) THEN
        NB_BYTES = NB_BYTES - int(size(PERM_RHS),8)*K34_8
        DEALLOCATE(PERM_RHS)
      ENDIF
      IF (allocated(UNS_PERM_INV)) THEN
        NB_BYTES = NB_BYTES - int(size(UNS_PERM_INV),8)*K34_8
        DEALLOCATE(UNS_PERM_INV)
      ENDIF
      IF (associated(id%BUFR)) THEN
          NB_BYTES = NB_BYTES - int(size(id%BUFR),8)*K34_8
          DEALLOCATE(id%BUFR)
          NULLIFY(id%BUFR)
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        IF (allocated(IWK_SOLVE)) THEN 
          NB_BYTES = NB_BYTES - int(size(IWK_SOLVE),8)*K34_8
          DEALLOCATE( IWK_SOLVE )
        ENDIF
        IF (allocated(IWCB)) THEN 
          NB_BYTES = NB_BYTES - int(size(IWCB),8)*K34_8
          DEALLOCATE( IWCB )
        ENDIF
        CALL SMUMPS_57( IERR )
        CALL SMUMPS_59( IERR )
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF (allocated(SAVERHS)) THEN 
         NB_BYTES = NB_BYTES - int(size(SAVERHS),8)*K35_8
         DEALLOCATE( SAVERHS)
        ENDIF
        IF (
     &       ( 
     &        ( KEEP(248).ne.0 .OR. KEEP(221).EQ.2 
     &          .OR. KEEP(111).NE.0 ) 
     &        .and. ICNTL21.ne.0 ) 
     &     .or.
     &    ( KEEP(237).NE.0 )
     &     )
     &    THEN
          IF ( I_AM_SLAVE ) THEN
           IF (associated(RHS_MUMPS) ) THEN
            NB_BYTES = NB_BYTES - int(size(RHS_MUMPS),8)*K35_8
            DEALLOCATE(RHS_MUMPS)
           ENDIF
          ENDIF
        ENDIF
        NULLIFY(RHS_MUMPS)
      ELSE
        IF (associated(RHS_MUMPS)) THEN
          NB_BYTES = NB_BYTES - int(size(RHS_MUMPS),8)*K35_8
          DEALLOCATE(RHS_MUMPS)
          NULLIFY(RHS_MUMPS)
        END IF
      END IF
      IF (I_AM_SLAVE) THEN
        IF (allocated(SRW3)) THEN 
          NB_BYTES = NB_BYTES - int(size(SRW3),8)*K35_8
          DEALLOCATE(SRW3)
        ENDIF
        IF (allocated(POSINRHSCOMP_N)) THEN 
          NB_BYTES = NB_BYTES - int(size(POSINRHSCOMP_N),8)*K34_8
          DEALLOCATE(POSINRHSCOMP_N)
        ENDIF
        IF (LSCAL .AND. ICNTL21==1) THEN
          NB_BYTES = NB_BYTES - 
     &              int(size(scaling_data%SCALING_LOC),8)*K16_8
          DEALLOCATE(scaling_data%SCALING_LOC)
          NULLIFY(scaling_data%SCALING_LOC)
        ENDIF
        IF (WK_USER_PROVIDED) THEN
          NULLIFY(id%S)
        ELSE IF (associated(id%S).AND.KEEP(201).GT.0) THEN
          NB_BYTES = NB_BYTES - KEEP8(23)*K35_8
          id%KEEP8(23)=0_8
          DEALLOCATE(id%S)
          NULLIFY(id%S)
        ENDIF
        IF (KEEP(221).NE.1) THEN
         IF (associated(id%RHSCOMP)) THEN 
            NB_BYTES = NB_BYTES - int(size(id%RHSCOMP),8)*K35_8
            DEALLOCATE(id%RHSCOMP)
            NULLIFY(id%RHSCOMP)
         ENDIF
         IF (associated(id%POSINRHSCOMP)) THEN
            NB_BYTES = NB_BYTES - int(size(id%POSINRHSCOMP),8)*K34_8
            DEALLOCATE(id%POSINRHSCOMP)
            NULLIFY(id%POSINRHSCOMP)
         ENDIF
        ENDIF
        IF ( WORK_WCB_ALLOCATED ) THEN
          NB_BYTES = NB_BYTES - int(size(WORK_WCB),8)*K35_8
          DEALLOCATE( WORK_WCB )
        ENDIF
        NULLIFY( WORK_WCB )
      ENDIF
      RETURN
 65   FORMAT (//' ERROR ANALYSIS AFTER ITERATIVE REFINEMENT')
 100  FORMAT(//' ****** SOLVE & CHECK STEP ********'/)
 110  FORMAT (//' VECTOR SOLUTION FOR COLUMN ',I12)
 115  FORMAT(1X, A44,1P,D9.2)
 150  FORMAT (/' STATISTICS PRIOR SOLVE PHASE     ...........'/
     &        ' NUMBER OF RIGHT-HAND-SIDES                    =',I12/
     &        ' BLOCKING FACTOR FOR MULTIPLE RHS              =',I12/
     &        ' ICNTL (9)                                     =',I12/
     &        '  --- (10)                                     =',I12/
     &        '  --- (11)                                     =',I12/
     &        '  --- (20)                                     =',I12/
     &        '  --- (21)                                     =',I12/
     &        '  --- (30)                                     =',I12)
 151  FORMAT ('  --- (25)                                     =',I12)
 152  FORMAT ('  --- (26)                                     =',I12)
 153  FORMAT ('  --- (32)                                     =',I12)
 160  FORMAT (' RHS'/(1X,1P,5E14.6))
 170  FORMAT (//' ERROR ANALYSIS' )
 240  FORMAT (1X, A42,I4)
 270  FORMAT (//' BEGIN ITERATIVE REFINEMENT' )
  81  FORMAT (/' STATISTICS AFTER ITERATIVE REFINEMENT ')
 131  FORMAT (/' END   ITERATIVE REFINEMENT ')
 141  FORMAT(1X, A42,I4)
      END SUBROUTINE SMUMPS_301
      SUBROUTINE SMUMPS_245(root, N, A, LA, IW, LIW, W, LWC, 
     & IWCB,LIWW,RHS,LRHS,NRHS,NA,LNA,NE_STEPS, W2,
     & MTYPE, ICNTL,
     & STEP, FRERE, DAD, FILS, PTRIST, PTRFAC, IW1,LIW1,
     & PROCNODE_STEPS, SLAVEF,
     & INFO, KEEP,KEEP8, COMM_NODES, MYID,
     & MYID_NODES,
     & BUFR, LBUFR, LBUFR_BYTES,
     & 
     & ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     & IBEG_ROOT_DEF, IEND_ROOT_DEF,
     & IROOT_DEF_RHS_COL1, PTR_RHS_ROOT, LPTR_RHS_ROOT,
     & SIZE_ROOT, MASTER_ROOT,
     & RHSCOMP, LRHSCOMP, POSINRHSCOMP, BUILD_POSINRHSCOMP
     & , NZ_RHS, NBCOL_INBLOC, NRHS_ORIG 
     & , JBEG_RHS
     & , Step2node, LStep2node 
     & , IRHS_SPARSE           
     & , IRHS_PTR               
     & , SIZE_PERM_RHS, PERM_RHS     
     & , SIZE_UNS_PERM_INV, UNS_PERM_INV    
     & )
      USE SMUMPS_OOC
      USE MUMPS_SOL_ES
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      TYPE ( SMUMPS_ROOT_STRUC ) :: root
      INTEGER(8) :: LA
      INTEGER LWC,N,LIW,MTYPE,LIW1,LIWW,LNA
      INTEGER ICNTL(40),INFO(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IW(LIW),IW1(LIW1),NA(LNA),NE_STEPS(KEEP(28)),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)), FILS(N), PTRIST(KEEP(28)),
     &        DAD(KEEP(28))
      INTEGER(8) ::  PTRFAC(KEEP(28))
      INTEGER LRHS, NRHS, LRHSCOMP
      REAL    A(LA), W(LWC), RHS(LRHS,NRHS),
     &        W2(KEEP(133)), 
     &        RHSCOMP(LRHSCOMP,NRHS)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER PROCNODE_STEPS(KEEP(28)), POSINRHSCOMP(KEEP(28))
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1
      INTEGER SIZE_ROOT, MASTER_ROOT
      INTEGER LPTR_RHS_ROOT
      REAL PTR_RHS_ROOT(LPTR_RHS_ROOT)
      LOGICAL BUILD_POSINRHSCOMP
      INTEGER MP, LP, LDIAG
      INTEGER K,I,II
      INTEGER allocok
      INTEGER LPOOL,MYLEAF,LPANEL_POS
      INTEGER NSTK_S,IPOOL,IPANEL_POS,PTRICB,PTRACB
      INTEGER MTYPE_LOC
      INTEGER IERR
      INTEGER(8) :: IAPOS
      INTEGER       IOLDPS,
     &              LOCAL_M,
     &              LOCAL_N
#if defined(V_T)
      INTEGER soln_c_class, forw_soln, back_soln, root_soln
#endif
      INTEGER IZERO
      LOGICAL DOFORWARD, DOROOT, DOBACKWARD
      LOGICAL I_WORKED_ON_ROOT, SPECIAL_ROOT_REACHED
      INTEGER IROOT
      LOGICAL DOROOT_FWD_OOC, DOROOT_BWD_PANEL
      LOGICAL SWITCH_OFF_ES
      LOGICAL DUMMY_BOOL
      PARAMETER (IZERO = 0 )
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INCLUDE 'mumps_headers.h'
      EXTERNAL SMUMPS_248, SMUMPS_249
      INTEGER, intent(in)      :: NZ_RHS, NBCOL_INBLOC, NRHS_ORIG
      INTEGER, intent(in)      :: SIZE_UNS_PERM_INV 
      INTEGER, intent(in)      :: SIZE_PERM_RHS 
      INTEGER, intent(in) :: JBEG_RHS
      INTEGER, intent(in) :: IRHS_SPARSE(NZ_RHS)
      INTEGER, intent(in) :: IRHS_PTR(NBCOL_INBLOC+1)
      INTEGER, intent(in) :: PERM_RHS(SIZE_PERM_RHS)
      INTEGER, intent(in) :: UNS_PERM_INV(SIZE_UNS_PERM_INV)
      INTEGER, intent(in) :: LStep2node
      INTEGER, intent(in) :: Step2node(LStep2node)
      INTEGER, DIMENSION(:), ALLOCATABLE ::  nodes_RHS
      INTEGER nb_nodes_RHS
      INTEGER nb_prun_leaves
      INTEGER, DIMENSION(:), ALLOCATABLE :: Pruned_Leaves
      INTEGER, DIMENSION(:), ALLOCATABLE ::  Pruned_List
      INTEGER  nb_prun_nodes
      INTEGER nb_prun_roots, JAM1
      INTEGER, DIMENSION(:), ALLOCATABLE ::  Pruned_SONS, Pruned_Roots
      INTEGER, DIMENSION(:), ALLOCATABLE ::  prun_NA
      INTEGER                            :: SIZE_TO_PROCESS
      LOGICAL, DIMENSION(:), ALLOCATABLE :: TO_PROCESS
      INTEGER ISTEP, INODE_PRINC
      LOGICAL AM1, DO_PRUN
      LOGICAL Exploit_Sparsity
      INTEGER :: OOC_FCT_TYPE_TMP
      INTEGER :: MUMPS_808
      EXTERNAL :: MUMPS_808
      MYLEAF = -1
      LP      = ICNTL(1)
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln_c',soln_c_class,ierr)
      CALL VTFUNCDEF( 'forw_soln',soln_c_class,forw_soln,ierr)
      CALL VTFUNCDEF( 'back_soln',soln_c_class,back_soln,ierr)
      CALL VTFUNCDEF( 'root_soln',soln_c_class,root_soln,ierr)
#endif
      NSTK_S   = 1
      PTRICB = NSTK_S + KEEP(28)
      PTRACB = PTRICB + KEEP(28)
      IPOOL  = PTRACB + KEEP(28)
      LPOOL  = KEEP(28)+1
      IPANEL_POS = IPOOL + LPOOL
      IF (KEEP(201).EQ.1) THEN
        LPANEL_POS = KEEP(228)+1
      ELSE
        LPANEL_POS = 1
      ENDIF
      IF (IPANEL_POS + LPANEL_POS -1 .ne. LIW1 )  THEN
         WRITE(*,*)  MYID, ": Internal Error in SMUMPS_245",
     &   IPANEL_POS, LPANEL_POS, LIW1
         CALL MUMPS_ABORT()
      ENDIF
      DOFORWARD = .TRUE.
      DOBACKWARD= .TRUE.
      SPECIAL_ROOT_REACHED = .TRUE.
      SWITCH_OFF_ES    = .FALSE.
      IF ( KEEP(111).NE.0 .OR. KEEP(252).NE.0 ) THEN
        DOFORWARD = .FALSE.
      ENDIF
      IF (KEEP(221).eq.1) DOBACKWARD = .FALSE.
      IF (KEEP(221).eq.2) DOFORWARD  = .FALSE.
      IF ( KEEP(60).EQ.0 .AND.
     &    ( 
     &      (KEEP(38).NE.0 .AND.  root%yes) 
     &  .OR.
     &      (KEEP(20).NE.0 .AND. MYID_NODES.EQ.MASTER_ROOT)
     &    ) 
     &  .AND. KEEP(252).EQ.0
     &   )
     &THEN
        DOROOT = .TRUE.
      ELSE
        DOROOT = .FALSE.
      ENDIF
      DOROOT_BWD_PANEL = DOROOT .AND. MTYPE.NE.1 .AND. KEEP(50).EQ.0
     &                     .AND. KEEP(201).EQ.1
      DOROOT_FWD_OOC = DOROOT .AND. .NOT.DOROOT_BWD_PANEL
      AM1              = (KEEP(237) .NE. 0)
      Exploit_Sparsity = (KEEP(235) .NE. 0) .AND. (.NOT. AM1)
      DO_PRUN          = (Exploit_Sparsity.OR.AM1)
      IF ( DO_PRUN ) THEN
         IF (.not. allocated(Pruned_SONS)) THEN
            ALLOCATE (Pruned_SONS(KEEP(28)), stat=I)
            IF(I.GT.0) THEN
               INFO(1)=-13
               INFO(2)=KEEP(28)
            END IF
            CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
            IF(INFO(1).LT.0) GOTO 500
         END IF
         IF ( allocated(TO_PROCESS)) DEALLOCATE (TO_PROCESS)
         IF (.not. allocated(TO_PROCESS)) THEN
             SIZE_TO_PROCESS = KEEP(28)
             ALLOCATE (TO_PROCESS(SIZE_TO_PROCESS), stat=I)
             IF(I.GT.0) THEN
                INFO(1)=-13
                INFO(2)=KEEP(28)
             END IF
             CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
             IF(INFO(1).LT.0) GOTO 500
         END IF
         TO_PROCESS(:) = .TRUE.
      ENDIF
      IF ( DOFORWARD .AND. DO_PRUN ) THEN
         nb_prun_nodes = 0
         nb_prun_roots = 0
         Pruned_SONS(:) = -1
         IF ( Exploit_Sparsity ) THEN
            nb_nodes_RHS = 0
            DO I = 1, NZ_RHS
               ISTEP       = abs( STEP(IRHS_SPARSE(I)) )
               INODE_PRINC = Step2node( ISTEP )
               IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                  nb_nodes_RHS = nb_nodes_RHS +1
                  Pruned_SONS(ISTEP) = 0 
               ENDIF
            ENDDO
            ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)
            IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_nodes_RHS
            END IF
            CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
            IF(INFO(1).LT.0) GOTO 500
            nb_nodes_RHS = 0
            Pruned_SONS = -1
            DO I = 1, NZ_RHS
               ISTEP       = abs( STEP(IRHS_SPARSE(I)) )
               INODE_PRINC = Step2node( ISTEP )
               IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                  nb_nodes_RHS = nb_nodes_RHS +1
                  nodes_RHS(nb_nodes_RHS)  = INODE_PRINC
                  Pruned_SONS(ISTEP) = 0 
               ENDIF
            ENDDO
         ELSE IF ( AM1 ) THEN  
#if defined(NOT_USED)  
           IF ( KEEP(201).GT.0) THEN
             CALL SMUMPS_789(KEEP(28),
     &           KEEP(38), KEEP(20) )
           ENDIF
#endif
            nb_nodes_RHS = 0
#if defined(check)
            WRITE(*,*) "NBCOL_INBLOC=",NBCOL_INBLOC
            WRITE(*,*) "JBEG SIZE=",JBEG_RHS, SIZE(IRHS_PTR)
#endif
            DO I = 1, NBCOL_INBLOC
              IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
              IF ( (KEEP(242) .NE. 0 ).OR. (KEEP(243).NE.0) ) THEN
                   JAM1 = PERM_RHS(JBEG_RHS+I-1)
              ELSE
                   JAM1 = JBEG_RHS+I-1
              ENDIF       
              ISTEP = abs(STEP(JAM1))
              INODE_PRINC = Step2node(ISTEP)
              IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                 nb_nodes_RHS = nb_nodes_RHS +1
                 Pruned_SONS(ISTEP) = 0                 
              ENDIF
            ENDDO
            ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)
            IF(allocok.GT.0) THEN
               INFO(1)=-13
               INFO(2)=nb_nodes_RHS
            END IF
            CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
            IF(INFO(1).LT.0) GOTO 500
            nb_nodes_RHS = 0
            Pruned_SONS = -1
            DO I = 1, NBCOL_INBLOC
              IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
              IF ( (KEEP(242) .NE. 0 ).OR. (KEEP(243).NE.0) ) THEN
                   JAM1 = PERM_RHS(JBEG_RHS+I-1)
              ELSE
                   JAM1 = JBEG_RHS+I-1
              ENDIF
              ISTEP = abs(STEP(JAM1))
              INODE_PRINC = Step2node(ISTEP)
              IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                 nb_nodes_RHS = nb_nodes_RHS +1
                 nodes_RHS(nb_nodes_RHS)  = INODE_PRINC
                 Pruned_SONS(ISTEP) = 0
              ENDIF
            ENDDO            
         ENDIF                  
         CALL MUMPS_797( 
     &        .FALSE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves ) 
         ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_nodes
         END IF
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_roots
         END IF
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_leaves
         END IF
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         CALL MUMPS_797( 
     &        .TRUE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &        Pruned_List, Pruned_Roots, Pruned_Leaves )
         IF(allocated(nodes_RHS)) DEALLOCATE(nodes_RHS)
         CALL SMUMPS_809(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
         IF ( KEEP(201) .GT. 0) THEN
         OOC_FCT_TYPE_TMP=MUMPS_808
     &                    ('F',MTYPE,KEEP(201),KEEP(50))
         ELSE
         OOC_FCT_TYPE_TMP = -5959 
         ENDIF
         CALL MUMPS_802(
     &        MYID_NODES, N, KEEP(28), KEEP(201), KEEP8(31),
     &        STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_TMP
     &   )
         SPECIAL_ROOT_REACHED = .FALSE.
         DO I= 1, nb_prun_roots
          IF ( (Pruned_Roots(I).EQ.KEEP(38)).OR.
     &         (Pruned_Roots(I).EQ.KEEP(20)) ) THEN
            SPECIAL_ROOT_REACHED = .TRUE.
            EXIT
          ENDIF
         ENDDO
      ENDIF  
      IF (KEEP(201).GT.0) THEN
        IF (DOFORWARD .OR. DOROOT_FWD_OOC) THEN
           CALL SMUMPS_583(PTRFAC,KEEP(28),MTYPE,
     &                                A,LA,DOFORWARD,IERR)
          IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            CALL MUMPS_ABORT()
          ENDIF
        ENDIF
      ENDIF
      IF (DOFORWARD) THEN
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = 1
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(forw_soln,ierr)
#endif
      IF (.NOT.DO_PRUN) THEN
        CALL SMUMPS_248(N, A(1), LA, IW(1), LIW, W(1),
     &           LWC, RHS, LRHS, NRHS,
     &           IW1(PTRICB), IWCB, LIWW,
     &           RHSCOMP,LRHSCOMP,POSINRHSCOMP,BUILD_POSINRHSCOMP,
     &           NE_STEPS, NA, LNA, STEP, FRERE,DAD,FILS,
     &           IW1(NSTK_S),IW1(IPOOL),LPOOL,PTRIST,PTRFAC,
     &           MYLEAF,INFO,
     &           KEEP,KEEP8,
     &           PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     &           BUFR, LBUFR, LBUFR_BYTES,
     &           PTR_RHS_ROOT, LPTR_RHS_ROOT, MTYPE_LOC, 
     & 
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &           )
      ELSE
         ALLOCATE(prun_NA(nb_prun_leaves+nb_prun_roots+2),
     &     STAT=allocok)
         IF(allocok.GT.0) THEN
           INFO(1)=-13
           INFO(2)=nb_prun_leaves+nb_prun_roots+2
         END IF
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF(I.LT.0) GOTO 500
         prun_NA(1) = nb_prun_leaves
         prun_NA(2) = nb_prun_roots
         DO I = 1, nb_prun_leaves
            prun_NA(I+2) = Pruned_Leaves(I)
         ENDDO
         DO I = 1, nb_prun_roots
            prun_NA(I+2+nb_prun_leaves) = Pruned_Roots(I)
         ENDDO
         DEALLOCATE(Pruned_List)
         DEALLOCATE(Pruned_Leaves)
         IF (AM1) THEN 
           DEALLOCATE(Pruned_Roots)
         END IF
         IF ((Exploit_Sparsity).AND.(nb_prun_roots.EQ.NA(2))) THEN
           DEALLOCATE(Pruned_Roots)
           IF ( allocated(TO_PROCESS)) DEALLOCATE (TO_PROCESS)
           SWITCH_OFF_ES = .TRUE.   
         ENDIF
         CALL SMUMPS_248(N, A(1), LA, IW(1), LIW, W(1),
     &        LWC, RHS, LRHS, NRHS,
     &        IW1(PTRICB), IWCB, LIWW,
     &        RHSCOMP,LRHSCOMP,POSINRHSCOMP,BUILD_POSINRHSCOMP,
     &        Pruned_SONS, prun_NA, LNA, STEP, FRERE,DAD,FILS,
     &        IW1(NSTK_S),IW1(IPOOL),LPOOL,PTRIST,PTRFAC,
     &        MYLEAF,INFO,
     &        KEEP,KEEP8,
     &        PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     &        BUFR, LBUFR, LBUFR_BYTES,
     &        PTR_RHS_ROOT, LPTR_RHS_ROOT, MTYPE_LOC, 
     &        
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &        )
         DEALLOCATE(prun_NA)
      ENDIF
      BUILD_POSINRHSCOMP = .FALSE.
#if defined(V_T)
        CALL VTEND(forw_soln,ierr)
#endif
      ENDIF
      CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
      IF ( INFO(1) .LT. 0 ) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP,*) MYID,
     &    ': ** ERROR RETURN FROM SMUMPS_248,INFO(1:2)=',
     &    INFO(1:2)
        END IF
        GOTO 500   
      END IF
      CALL MPI_BARRIER( COMM_NODES, IERR )
      IF (DO_PRUN.AND.SWITCH_OFF_ES) THEN
         DO_PRUN          = .FALSE.
         Exploit_Sparsity = .FALSE.
      ENDIF 
      IF ( DOBACKWARD .AND. DO_PRUN )  THEN
         nb_prun_leaves = 0
        IF ( Exploit_Sparsity .AND. (KEEP(111).EQ.0) ) THEN
          nb_nodes_RHS = nb_prun_roots
          ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)
          IF(allocok.GT.0) THEN
            WRITE(*,*)'Problem with allocation of nodes_RHS'
            INFO(1) = -13
            INFO(2) = nb_nodes_RHS
            CALL MUMPS_ABORT()
          END IF
          nodes_RHS(1:nb_prun_roots)=Pruned_Roots(1:nb_prun_roots)
          DEALLOCATE(Pruned_Roots)
        ELSE
         nb_nodes_RHS = 0
         Pruned_SONS(:) = -1  
         DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)          
            ISTEP = abs(STEP(I))
            INODE_PRINC = Step2node(ISTEP)
            IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
               nb_nodes_RHS = nb_nodes_RHS +1
               Pruned_SONS(ISTEP) = 0
            ENDIF
         ENDDO
         ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)           
         IF(allocok.GT.0) THEN
           WRITE(*,*)'Problem with allocation of nodes_RHS'
           INFO(1) = -13
           INFO(2) = nb_nodes_RHS
           CALL MUMPS_ABORT()
         END IF
         nb_nodes_RHS = 0         
         Pruned_SONS(:) = -1  
         DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
            ISTEP = abs(STEP(I))
            INODE_PRINC = Step2node(ISTEP)
            IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
               nb_nodes_RHS = nb_nodes_RHS +1
               nodes_RHS(nb_nodes_RHS)  = INODE_PRINC
               Pruned_SONS(ISTEP) = 0
            ENDIF
         ENDDO
        ENDIF
        IF ( Exploit_Sparsity ) THEN
           CALL MUMPS_798( 
     &     .FALSE.,
     &     DAD, NE_STEPS, FRERE, KEEP(28),
     &     FILS, STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves
     &     )
           ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_nodes
           END IF
           CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_roots
           END IF
           CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_leaves
           END IF
           CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           CALL MUMPS_798( 
     &     .TRUE.,
     &     DAD, NE_STEPS, FRERE, KEEP(28),
     &     FILS, STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
           CALL SMUMPS_809(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
           IF(allocated(nodes_RHS)) DEALLOCATE(nodes_RHS)
           IF (KEEP(201).GT.0) THEN
             OOC_FCT_TYPE_TMP=MUMPS_808
     &                    ('B',MTYPE,KEEP(201),KEEP(50))
           ELSE
             OOC_FCT_TYPE_TMP = -5959 
           ENDIF
           CALL MUMPS_803(
     &          MYID_NODES, N, KEEP(28), KEEP(201),
     &          KEEP8(31), STEP,
     &          Pruned_List,
     &          nb_prun_nodes, OOC_FCT_TYPE_TMP)
        ENDIF
      ENDIF
      IF(KEEP(201).EQ.1.AND.DOROOT_BWD_PANEL) THEN
         I_WORKED_ON_ROOT = .FALSE. 
         CALL SMUMPS_584(PTRFAC,KEEP(28),MTYPE,
     &   I_WORKED_ON_ROOT, IROOT, A, LA, IERR)
         IF (IERR .LT. 0) THEN
           INFO(1) = -90
           INFO(2) = IERR
         ENDIF 
      ENDIF
      IF (KEEP(201).EQ.1) THEN
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF ( INFO(1) .LT. 0 ) GOTO 500  
      ENDIF
      IF (KEEP(60).NE.0 .AND. KEEP(221).EQ.0
     &   .AND. MYID_NODES .EQ. MASTER_ROOT) THEN
        PTR_RHS_ROOT(1:NRHS*SIZE_ROOT) = ZERO
      ENDIF
      IF ( ( KEEP( 38 ) .NE. 0 ).AND. SPECIAL_ROOT_REACHED ) THEN
        IF ( KEEP(60) .EQ. 0 .AND. KEEP(252) .EQ. 0 ) THEN
          IF ( root%yes ) THEN
            IF (KEEP(201).GT.0) THEN
              IF ( (Exploit_Sparsity.AND.(KEEP(111).NE.0)) .and.
     &            (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) ) THEN
                  write(6,*) " CPA to be double checked "
                  GOTO 1010
              ENDIF
            ENDIF
            IOLDPS = PTRIST(STEP(KEEP(38)))
            LOCAL_M = IW( IOLDPS + 2 + KEEP(IXSZ))
            LOCAL_N = IW( IOLDPS + 1 + KEEP(IXSZ))
            IF (KEEP(201).GT.0) THEN
              CALL SMUMPS_643(
     &           KEEP(38),PTRFAC,KEEP,A,LA,
     &           STEP,KEEP8,N,DUMMY_BOOL,IERR)
              IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                WRITE(*,*) '** ERROR after SMUMPS_643',
     &          INFO(1)
                call MUMPS_ABORT()
              ENDIF
            ENDIF
            IAPOS   = PTRFAC(IW( IOLDPS + 4 + KEEP(IXSZ)))
#if defined(V_T)
            CALL VTBEGIN(root_soln,ierr)
#endif
             CALL SMUMPS_286( NRHS, root%DESCRIPTOR(1), 
     &       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     &       root%MBLOCK, root%NBLOCK,
     &       root%IPIV(1), root%LPIV, MASTER_ROOT, MYID_NODES,
     &       COMM_NODES,
     &       PTR_RHS_ROOT(1),
     &       root%TOT_ROOT_SIZE, A( IAPOS ),
     &       INFO(1), MTYPE, KEEP(50))
            IF(KEEP(201).GT.0)THEN
              CALL SMUMPS_598(KEEP(38),
     &             PTRFAC,KEEP(28),A,LA,.FALSE.,IERR)
              IF(IERR.LT.0)THEN
                 INFO(1)=IERR
                 INFO(2)=0
                 WRITE(*,*)
     &           '** ERROR after SMUMPS_598 ',
     &           INFO(1)
                 call MUMPS_ABORT()
              ENDIF
            ENDIF
          ENDIF  
        ENDIF
      ELSE IF ( ( KEEP(20) .NE. 0) .AND. SPECIAL_ROOT_REACHED ) THEN
        IF ( MYID_NODES .eq.  MASTER_ROOT ) THEN
        END IF 
      END IF 
#if defined(V_T)
      CALL VTEND(root_soln,ierr)
#endif
 1010 CONTINUE
      CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      IF (DOBACKWARD) THEN
        IF ( KEEP(201).GT.0 .AND.  .NOT. DOROOT_BWD_PANEL )
     &    THEN
          I_WORKED_ON_ROOT = DOROOT
          IF (KEEP(111).NE.0) 
     &         I_WORKED_ON_ROOT = .FALSE.
          IF (KEEP(38).gt.0 ) THEN
             IF ( ( Exploit_Sparsity.AND.(KEEP(111).EQ.0) )
     &            .OR. AM1 ) THEN
                IF (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) THEN
                   OOC_STATE_NODE(STEP(KEEP(38)))=-4
                ENDIF
             ENDIF
             IF (Exploit_Sparsity.AND.(KEEP(111).NE.0)) THEN
                IF (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) THEN
                   I_WORKED_ON_ROOT = .FALSE.
                ENDIF
             ENDIF
          ENDIF
        ENDIF                    
        IF ( AM1 ) THEN
         CALL MUMPS_797( 
     &        .FALSE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves) 
         ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_nodes
         END IF
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_roots
         END IF
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_leaves
         END IF
         CALL MUMPS_276(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         CALL MUMPS_797( 
     &        .TRUE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &        Pruned_List, Pruned_Roots, Pruned_Leaves )
          CALL SMUMPS_809(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
          IF (KEEP(201).GT.0) THEN
           OOC_FCT_TYPE_TMP=MUMPS_808
     &                    ('B',MTYPE,KEEP(201),KEEP(50))
          ELSE
           OOC_FCT_TYPE_TMP = -5959 
          ENDIF
          CALL MUMPS_802(
     &    MYID_NODES, N, KEEP(28), KEEP(201), KEEP8(31),
     &    STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_TMP
     &    )
        ENDIF
        IF ( KEEP(201).GT.0 ) THEN
          IROOT = max(KEEP(20),KEEP(38)) 
          CALL SMUMPS_584(PTRFAC,KEEP(28),MTYPE,
     &         I_WORKED_ON_ROOT, IROOT, A, LA, IERR)
        ENDIF
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = IZERO
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(back_soln,ierr)
#endif
        IF ( .NOT.SPECIAL_ROOT_REACHED ) THEN
          PTR_RHS_ROOT(1:NRHS*SIZE_ROOT) = ZERO
        ENDIF
        IF ( .NOT. DO_PRUN ) THEN
           SIZE_TO_PROCESS = 1
           IF ( allocated(TO_PROCESS)) DEALLOCATE (TO_PROCESS)  
           ALLOCATE (TO_PROCESS(SIZE_TO_PROCESS), stat=I)
           TO_PROCESS(:) = .TRUE.
           CALL SMUMPS_249( N, A, LA, IW, LIW, W(1), LWC,
     &          RHS, LRHS, NRHS,
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     &          IW1(PTRICB),IW1(PTRACB),IWCB,LIWW,
     &          W2, NE_STEPS, NA, LNA, STEP, FRERE,DAD,FILS,
     &          IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,INFO,
     &          PROCNODE_STEPS, SLAVEF, COMM_NODES,MYID_NODES,
     &          BUFR, LBUFR, LBUFR_BYTES, KEEP,KEEP8,
     &          PTR_RHS_ROOT, LPTR_RHS_ROOT,
     &          MTYPE_LOC, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, IW1(IPANEL_POS),
     &          LPANEL_POS, TO_PROCESS, SIZE_TO_PROCESS)
        ELSE       
          ALLOCATE(prun_NA(nb_prun_leaves+nb_prun_roots+2),
     &      STAT=allocok)
          IF(allocok.GT.0) THEN
            WRITE(*,*)'Problem with allocation of prun_na'
            CALL MUMPS_ABORT()
          END IF
          prun_NA(1) = nb_prun_leaves
          prun_NA(2) = nb_prun_roots
          DO I = 1, nb_prun_leaves
              prun_NA(I+2) = Pruned_Leaves(I)
          ENDDO
          DO I = 1, nb_prun_roots
              prun_NA(I+2+nb_prun_leaves) = Pruned_Roots(I)
          ENDDO
          CALL SMUMPS_249( N, A, LA, IW, LIW, W(1), LWC,
     &          RHS, LRHS, NRHS,
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     &          IW1(PTRICB),IW1(PTRACB),IWCB,LIWW,       
     &          W2, NE_STEPS, prun_NA, LNA, STEP, FRERE,DAD,FILS,
     &          IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,INFO,
     &          PROCNODE_STEPS, SLAVEF, COMM_NODES,MYID_NODES,
     &          BUFR, LBUFR, LBUFR_BYTES, KEEP,KEEP8,
     &          PTR_RHS_ROOT, LPTR_RHS_ROOT,
     &          MTYPE_LOC, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, IW1(IPANEL_POS),
     &          LPANEL_POS, TO_PROCESS, SIZE_TO_PROCESS)
        ENDIF
#if defined(V_T)
        CALL VTEND(back_soln,ierr)
#endif
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        IF (DOFORWARD) THEN
        K = min0(10,N)
        IF (LDIAG.EQ.4) K = N
        WRITE (MP,99992)
        IF (N.GT.0) WRITE (MP,99993) (RHS(I,1),I=1,K)
        IF (N.GT.0.and.NRHS>1) 
     &              WRITE (MP,99994) (RHS(I,2),I=1,K)
        ENDIF
      ENDIF
500   CONTINUE
      IF ( allocated(TO_PROCESS)) DEALLOCATE (TO_PROCESS)
      IF (Exploit_Sparsity.OR.AM1.OR.SWITCH_OFF_ES) THEN
         IF ( allocated(nodes_RHS)) DEALLOCATE (nodes_RHS)
         IF ( allocated(Pruned_SONS)) DEALLOCATE (Pruned_SONS)
         IF ( allocated(Pruned_Roots)) DEALLOCATE (Pruned_Roots)
         IF ( allocated(prun_NA)) DEALLOCATE (prun_NA)
         IF ( allocated(Pruned_List)) DEALLOCATE (Pruned_List)
         IF ( allocated(Pruned_Leaves)) DEALLOCATE (Pruned_Leaves)
      ENDIF
      RETURN
99993 FORMAT (' RHS    (first column)'/(1X,1P,5E14.6))
99994 FORMAT (' RHS    (2 nd  column)'/(1X,1P,5E14.6))
99992 FORMAT (//' LEAVING SOLVE (MPI41C) WITH')
      END SUBROUTINE SMUMPS_245
      SUBROUTINE SMUMPS_521(NSLAVES, N, MYID, COMM,
     &           MTYPE, RHS, LRHS, NRHS, PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, BUFFER,
     &           SIZE_BUF, SIZE_BUF_BYTES, CWORK, LCWORK,
     &           LSCAL, SCALING, LSCALING)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE
      INTEGER NRHS, LRHS, LCWORK
      REAL RHS   (LRHS, NRHS)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      REAL ::  CWORK(LCWORK)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES
      INTEGER BUFFER(SIZE_BUF)
      LOGICAL, intent(in) :: LSCAL
      INTEGER, intent(in) :: LSCALING
      REAL, intent(in)    :: SCALING(LSCALING)
      INTEGER I, II, J, J1, ISTEP, MASTER,
     &        MYID_NODES, TYPE_PARAL, N2RECV
      INTEGER LIELL, IPOS, NPIV, MAXNPIV_estim, MAXSurf
      INTEGER STATUS(MPI_STATUS_SIZE), IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND
      INTEGER SK38, SK20
      INTEGER, PARAMETER :: FIN = -1
      INTEGER, PARAMETER :: yes =  1
      INTEGER, PARAMETER :: no  = 0
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IROWlist(:) 
      INTEGER :: ONE_PACK
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      TYPE_PARAL = KEEP(46)  
      I_AM_SLAVE = MYID .ne. MASTER .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      IF (NSLAVES.EQ.1 .AND. TYPE_PARAL.EQ.1) THEN 
        IF (LSCAL) THEN
          DO J=1, NRHS
              DO I=1,N
                RHS(I,J) = RHS(I,J)*SCALING(I)
              ENDDO
          ENDDO
        ENDIF
        RETURN
      ENDIF
      IF (NSLAVES.EQ.1 .AND. TYPE_PARAL.EQ.0) THEN
        DO J=1, NRHS
           IF ( I_AM_SLAVE ) THEN
             CALL MPI_SEND(RHS(1, J), N, MPI_REAL, MASTER,
     &                 GatherSol, COMM, IERR)
     & 
           ELSE
             CALL MPI_RECV(RHS(1, J), N, MPI_REAL,
     &                 1,
     &                 GatherSol, COMM, STATUS, IERR )
             IF (LSCAL) THEN
              DO I=1,N
                RHS(I,J) = RHS(I,J)*SCALING(I)
              ENDDO
             ENDIF
           ENDIF
        ENDDO
        RETURN
      ENDIF
      MAXNPIV_estim = max(KEEP(246), KEEP(247))
      MAXSurf       = MAXNPIV_estim*NRHS
      IF (LCWORK .GE. MAXSurf) THEN
        ONE_PACK = yes 
      ELSE IF (LCWORK .GE. MAXNPIV_estim) THEN
        ONE_PACK = no 
      ELSE
        WRITE(*,*)
     &  "Internal error 2 in SMUMPS_521:",
     &  TYPE_PARAL, LCWORK, KEEP(247), NRHS
        CALL MUMPS_ABORT()
      ENDIF
      IF (ONE_PACK .EQ. no .AND. I_AM_SLAVE) THEN
          WRITE(*,*)
     &    "Internal error 1 in SMUMPS_521:",
     &    TYPE_PARAL, LCWORK, KEEP(246),KEEP(247), NRHS
          CALL MUMPS_ABORT()
      ENDIF
      IF (TYPE_PARAL .EQ. 0)
     &CALL MPI_BCAST(ONE_PACK, 1, MPI_INTEGER,
     &     MASTER, COMM, IERR)
      IF (MYID.EQ.MASTER) THEN
        ALLOCATE(IROWlist(KEEP(247)))
      ENDIF
      IF (NSLAVES .EQ. 1 .AND. TYPE_PARAL .EQ. 1) THEN
        CALL MUMPS_ABORT()
      ENDIF
      SIZE1=0
      CALL MPI_PACK_SIZE(MAXNPIV_estim+2,MPI_INTEGER, COMM, 
     &          SIZE1, IERR)
      SIZE2=0
      CALL MPI_PACK_SIZE(MAXSurf,MPI_REAL, COMM,
     &                   SIZE2, IERR)
      RECORD_SIZE_P_1= SIZE1+SIZE2
      IF (RECORD_SIZE_P_1.GT.SIZE_BUF_BYTES) THEN
         write(6,*) MYID, 
     &    ' Internal error 3 in  SMUMPS_521 '
         write(6,*) MYID, ' RECORD_SIZE_P_1, SIZE_BUF_BYTES=', 
     &                 RECORD_SIZE_P_1, SIZE_BUF_BYTES
         CALL MUMPS_ABORT()
      ENDIF
      N2SEND   =0
      N2RECV   =N
      POS_BUF  =0
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      IF (I_AM_SLAVE) THEN
        POS_BUF = 0
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == MUMPS_275(PROCNODE_STEPS(ISTEP),
     &          NSLAVES)) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                    IPOS = PTRIST(ISTEP)+KEEP(IXSZ) 
                    NPIV  = IW(IPOS+3)
                    LIELL = IW(IPOS) + NPIV
                    IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
              ELSE
                  IPOS = PTRIST(ISTEP) + 2+ KEEP(IXSZ)
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              IF (MYID .EQ. MASTER) THEN
                   N2RECV=N2RECV-NPIV
                   IF (NPIV.GT.0.AND.LSCAL) 
     &             CALL SMUMPS_522 ( ONE_PACK, .TRUE. )
              ELSE
                   IF (NPIV.GT.0) 
     &             CALL SMUMPS_522 ( ONE_PACK, .FALSE.)
              ENDIF
          ENDIF
        ENDDO
        CALL SMUMPS_523()   
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     &                 MPI_ANY_SOURCE,
     &                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (NPIV.NE.FIN)
          CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &             IROWlist, NPIV, MPI_INTEGER, COMM, IERR)
          IF (ONE_PACK.EQ.yes) THEN
            CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &             CWORK, NPIV*NRHS, MPI_REAL,
     &             COMM, IERR)
            IF (LSCAL) THEN
              DO J=1, NRHS
                DO I=1,NPIV
                  RHS(IROWlist(I),J)=
     &              CWORK(I+(J-1)*NPIV)*SCALING(IROWlist(I))
                ENDDO
              END DO
            ELSE
              DO J=1, NRHS
                DO I=1,NPIV
                  RHS(IROWlist(I),J)= CWORK(I+(J-1)*NPIV)
                ENDDO
              END DO
            ENDIF
          ELSE 
            DO J=1,NRHS
              CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                   CWORK, NPIV, MPI_REAL,
     &                   COMM, IERR)
              IF (LSCAL) THEN
               DO I=1,NPIV
                RHS(IROWlist(I),J)=CWORK(I)*SCALING(IROWlist(I))
               ENDDO
              ELSE
               DO I=1,NPIV
                RHS(IROWlist(I),J)=CWORK(I)
               ENDDO
              ENDIF
            ENDDO
          ENDIF
          N2RECV=N2RECV-NPIV
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR)
        ENDDO
       ENDDO
       DEALLOCATE(IROWlist)
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE SMUMPS_522 ( ONE_PACK, SCALE_ONLY )
        INTEGER, intent(in) ::  ONE_PACK      
        LOGICAL, intent(in) ::  SCALE_ONLY    
        INTEGER III
        IF (SCALE_ONLY) THEN
         DO II=1,NPIV
              I=IW(J1+II-1)
              DO J=1, NRHS
                RHS(I,J) = RHS(I,J)*SCALING(I)
              ENDDO
         ENDDO
         RETURN
        ENDIF
        DO II=1,NPIV
              I=IW(J1+II-1)
              DO J=1, NRHS
                CWORK(II+(J-1)*NPIV) = RHS(I,J)
              ENDDO
        ENDDO
        CALL MPI_PACK(NPIV, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(IW(J1), NPIV, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        IF (ONE_PACK.EQ.yes) THEN
           CALL MPI_PACK(CWORK(1), NPIV*NRHS, MPI_REAL,
     &                BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     &                IERR)
        ELSE
         III = 1
         DO J=1,NRHS
           CALL MPI_PACK(CWORK(III), NPIV, MPI_REAL,
     &                BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     &                IERR)
           III =III+NPIV
         ENDDO
        ENDIF
        N2SEND=N2SEND+NPIV  
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL SMUMPS_523()
        END IF
        RETURN
        END SUBROUTINE SMUMPS_522
        SUBROUTINE SMUMPS_523()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(FIN, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     &                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE SMUMPS_523
      END SUBROUTINE SMUMPS_521
      SUBROUTINE SMUMPS_812(NSLAVES, N, MYID, COMM,
     &           RHS, LRHS, NRHS, KEEP, BUFFER,
     &           SIZE_BUF, SIZE_BUF_BYTES, 
     &           LSCAL, SCALING, LSCALING,
     &          IRHS_PTR_COPY, LIRHS_PTR_COPY, 
     &          IRHS_SPARSE_COPY, LIRHS_SPARSE_COPY,
     &          RHS_SPARSE_COPY, LRHS_SPARSE_COPY,
     &          UNS_PERM_INV, LUNS_PERM_INV,
     &          POSINRHSCOMP_N, LPOS_N )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM
      INTEGER NRHS, LRHS, LPOS_N
      REAL RHS   (LRHS, NRHS)
      INTEGER KEEP(500)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES
      INTEGER BUFFER(SIZE_BUF)
      INTEGER, intent(in) :: LIRHS_PTR_COPY, LIRHS_SPARSE_COPY, 
     &                       LRHS_SPARSE_COPY, LUNS_PERM_INV
      INTEGER :: IRHS_SPARSE_COPY(LIRHS_SPARSE_COPY), 
     &           IRHS_PTR_COPY(LIRHS_PTR_COPY), 
     &           UNS_PERM_INV(LUNS_PERM_INV), 
     &           POSINRHSCOMP_N(LPOS_N)
      REAL :: RHS_SPARSE_COPY(LRHS_SPARSE_COPY)
      LOGICAL, intent(in) :: LSCAL
      INTEGER, intent(in) :: LSCALING
      REAL, intent(in)    :: SCALING(LSCALING)
      INTEGER COLSIZE, K, IZ, IPREV, NBCOL_INBLOC
      INTEGER I, II, J, MASTER,
     &         TYPE_PARAL, N2RECV
      INTEGER STATUS(MPI_STATUS_SIZE), IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND
      INTEGER, PARAMETER :: FIN = -1
      INCLUDE 'mumps_headers.h'
      TYPE_PARAL = KEEP(46)  
      I_AM_SLAVE = MYID .ne. MASTER .OR. TYPE_PARAL .eq. 1
      NBCOL_INBLOC = size(IRHS_PTR_COPY)-1
      IF (NSLAVES.EQ.1 .AND. TYPE_PARAL.EQ.1) THEN 
        K=1              
        DO J = 1,  NBCOL_INBLOC
           COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
           IF (COLSIZE.EQ.0) CYCLE
           DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
             IF (POSINRHSCOMP_N(I).NE.0) THEN
                IF (LSCAL) THEN
                 RHS_SPARSE_COPY(IZ)=RHS(I,K)*SCALING(I)
                ELSE
                 RHS_SPARSE_COPY(IZ)=RHS(I,K)
                ENDIF
             ENDIF          
           ENDDO
           K = K + 1
        ENDDO
        RETURN
      ENDIF
      IF (I_AM_SLAVE) THEN
        K=1              
        DO J = 1, NBCOL_INBLOC
           COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
           IF (COLSIZE.EQ.0) CYCLE
           DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
             IF (POSINRHSCOMP_N(I).NE.0) THEN
               RHS_SPARSE_COPY(IZ)=RHS(I,K)
             ENDIF          
           ENDDO
           K = K + 1
        ENDDO
      ENDIF
      SIZE1=0
      CALL MPI_PACK_SIZE(3,MPI_INTEGER, COMM,  
     &          SIZE1, IERR)
      SIZE2=0
      CALL MPI_PACK_SIZE(1,MPI_REAL, COMM,
     &                   SIZE2, IERR)
      RECORD_SIZE_P_1= SIZE1+SIZE2
      IF (RECORD_SIZE_P_1.GT.SIZE_BUF_BYTES) THEN
         write(6,*) MYID, 
     &    ' Internal error 3 in  SMUMPS_812 '
         write(6,*) MYID, ' RECORD_SIZE_P_1, SIZE_BUF_BYTES=', 
     &                 RECORD_SIZE_P_1, SIZE_BUF_BYTES
         CALL MUMPS_ABORT()
      ENDIF
      N2SEND   =0
      N2RECV   =size(IRHS_SPARSE_COPY)
      POS_BUF  =0
      IF (I_AM_SLAVE) THEN
        DO J = 1,  NBCOL_INBLOC
            COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
            IF (COLSIZE.LE.0) CYCLE
            K = 0 
            DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
              I = IRHS_SPARSE_COPY(IZ)
              II = I
              IF  (KEEP(23).NE.0) II = UNS_PERM_INV(I)
              IF (POSINRHSCOMP_N(II).NE.0) THEN
               IF (MYID .EQ. MASTER) THEN
                  N2RECV=N2RECV-1
                  IF (LSCAL) 
     &            CALL SMUMPS_813 ( .TRUE. )
                  IRHS_SPARSE_COPY( IRHS_PTR_COPY(J) + K) =
     &               I
                  RHS_SPARSE_COPY( IRHS_PTR_COPY(J) + K) =
     &                RHS_SPARSE_COPY(IZ)
                  K = K+1 
               ELSE
                  CALL SMUMPS_813 (  .FALSE. )
               ENDIF
              ENDIF          
            ENDDO
            IF (MYID.EQ.MASTER) 
     &             IRHS_PTR_COPY(J) = IRHS_PTR_COPY(J) + K
        ENDDO
        CALL SMUMPS_814()   
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     &                 MPI_ANY_SOURCE,
     &                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &                   J, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (J.NE.FIN)
          IZ = IRHS_PTR_COPY(J)
          CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &             I, 1, MPI_INTEGER, COMM, IERR)
          IRHS_SPARSE_COPY(IZ) = I
          CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &             RHS_SPARSE_COPY(IZ), 1, MPI_REAL,
     &             COMM, IERR)
          IF (LSCAL) THEN
              IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
              RHS_SPARSE_COPY(IZ) = RHS_SPARSE_COPY(IZ)*SCALING(I)    
          ENDIF
          N2RECV=N2RECV-1
          IRHS_PTR_COPY(J) = IRHS_PTR_COPY(J) + 1
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                   J, 1, MPI_INTEGER, COMM, IERR)
        ENDDO
       ENDDO
       IPREV = 1
       DO J=1, size(IRHS_PTR_COPY)-1
         I= IRHS_PTR_COPY(J) 
         IRHS_PTR_COPY(J) = IPREV
         IPREV = I
       ENDDO
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE SMUMPS_813 ( SCALE_ONLY )
        LOGICAL, intent(in) ::  SCALE_ONLY    
        INTEGER III
        IF (SCALE_ONLY) THEN
         III = I
         IF (KEEP(23).NE.0) III = UNS_PERM_INV(I)
         IF (LSCAL) THEN
            RHS_SPARSE_COPY(IZ)=RHS_SPARSE_COPY(IZ)*SCALING(III)
         ENDIF
         RETURN
        ENDIF
        CALL MPI_PACK(J, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(I, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(RHS_SPARSE_COPY(IZ), 1, MPI_REAL,
     &                BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     &                IERR)
        N2SEND=N2SEND+1  
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL SMUMPS_814()
        END IF
        RETURN
        END SUBROUTINE SMUMPS_813
        SUBROUTINE SMUMPS_814()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(FIN, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     &                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE SMUMPS_814
      END SUBROUTINE SMUMPS_812
      SUBROUTINE SMUMPS_535(MTYPE, ISOL_LOC,
     &             PTRIST, KEEP,KEEP8,
     &             IW, LIW_PASSED, MYID_NODES, N, STEP,
     &             PROCNODE, NSLAVES, scaling_data, LSCAL)
      IMPLICIT NONE
      INTEGER MTYPE, MYID_NODES, N, NSLAVES
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE(KEEP(28))
      INTEGER ISOL_LOC(KEEP(89)) 
      INTEGER LIW_PASSED
      INTEGER IW(LIW_PASSED)
      INTEGER STEP(N)
      LOGICAL LSCAL
      type scaling_data_t
        SEQUENCE
        REAL, dimension(:), pointer :: SCALING
        REAL, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type (scaling_data_t) :: scaling_data
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      INTEGER ISTEP, K
      INTEGER J1, IPOS, LIELL, NPIV, JJ
      INTEGER SK38,SK20
      INCLUDE 'mumps_headers.h'
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      K=0
      DO ISTEP=1, KEEP(28)
          IF ( MYID_NODES == MUMPS_275( PROCNODE(ISTEP),
     &                   NSLAVES)) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                    IPOS = PTRIST(ISTEP)+KEEP(IXSZ)
                    LIELL = IW(IPOS+3)
                    NPIV = LIELL
                    IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
              ELSE
                  IPOS = PTRIST(ISTEP) + 2 + KEEP(IXSZ)
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 + KEEP(IXSZ))
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                  K=K+1
                  ISOL_LOC(K)=IW(JJ)
                  IF (LSCAL) THEN
                    scaling_data%SCALING_LOC(K)=
     &              scaling_data%SCALING(IW(JJ))
                  ENDIF
              ENDDO
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_535
      SUBROUTINE SMUMPS_532(
     &           SLAVEF, N, MYID_NODES,
     &           MTYPE, RHS, LD_RHS, NRHS,
     &           ISOL_LOC, SOL_LOC, BEG_RHS, LSOL_LOC,
     &           PTRIST,
     &           PROCNODE_STEPS, KEEP,KEEP8, IW, LIW, STEP,
     &           scaling_data, LSCAL, NB_RHSSKIPPED)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      type scaling_data_t
        SEQUENCE
        REAL, dimension(:), pointer :: SCALING
        REAL, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      TYPE (scaling_data_t) :: scaling_data
      LOGICAL LSCAL
      INTEGER SLAVEF, N, MYID_NODES, LIW, MTYPE, NRHS, LD_RHS
      INTEGER LSOL_LOC, BEG_RHS, NB_RHSSKIPPED
      INTEGER ISOL_LOC(LSOL_LOC)
      REAL SOL_LOC( LSOL_LOC, BEG_RHS+NRHS+NB_RHSSKIPPED-1)
      REAL RHS(  LD_RHS , NRHS)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER JJ, J1, ISTEP, K, JEMPTY, JEND
      INTEGER IPOS, LIELL, NPIV
      LOGICAL ROOT
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      K=0
      JEMPTY = BEG_RHS+NB_RHSSKIPPED-1
      JEND   = BEG_RHS+NB_RHSSKIPPED+NRHS-1
        DO ISTEP = 1, KEEP(28)
            IF (MYID_NODES == MUMPS_275(PROCNODE_STEPS(ISTEP),
     &          SLAVEF)) THEN
              ROOT=.false.
              IF (KEEP(38).ne.0) ROOT = STEP(KEEP(38))==ISTEP
              IF (KEEP(20).ne.0) ROOT = STEP(KEEP(20))==ISTEP
              IF ( ROOT ) THEN
                    IPOS = PTRIST(ISTEP) + KEEP(IXSZ)
                    LIELL = IW(IPOS+3)
                    NPIV = LIELL
                    IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
              ELSE
                  IPOS = PTRIST(ISTEP) + 2 +KEEP(IXSZ)
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                K=K+1
                IF (NB_RHSSKIPPED.GT.0) 
     &            SOL_LOC(K, BEG_RHS:JEMPTY) = ZERO                
                IF (LSCAL) THEN
                  SOL_LOC(K,JEMPTY+1:JEND) =
     &            scaling_data%SCALING_LOC(K)*RHS(IW(JJ),1:NRHS)
                ELSE
                  SOL_LOC(K,JEMPTY+1:JEND) =
     &            RHS(IW(JJ),1:NRHS)
                ENDIF
              ENDDO
            ENDIF
        ENDDO
      RETURN
      END SUBROUTINE SMUMPS_532
      SUBROUTINE SMUMPS_638
     &           (NSLAVES, N, MYID, COMM,
     &           MTYPE, RHS, LRHS, NRHS, PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     &           POSINRHSCOMP, LENPOSINRHSCOMP,
     &           BUILD_POSINRHSCOMP, ICNTL, INFO)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE
      INTEGER NRHS, LRHS, LENPOSINRHSCOMP
      INTEGER ICNTL(40), INFO(40)
      REAL RHS   (LRHS, NRHS)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N), POSINRHSCOMP(LENPOSINRHSCOMP)
      LOGICAL BUILD_POSINRHSCOMP
      INTEGER BUF_MAXSIZE, BUF_MAXREF
      PARAMETER (BUF_MAXREF=200000)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF_INDX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: BUF_RHS
      INTEGER ENTRIES_2_PROCESS, PROC_WHO_ASKS, BUF_EFFSIZE
      INTEGER INDX 
      INTEGER allocok
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER I, K, JJ, J1, ISTEP, MASTER,
     &        MYID_NODES, TYPE_PARAL
      INTEGER LIELL, IPOS, NPIV
      INTEGER STATUS(MPI_STATUS_SIZE), IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER SK38, SK20, IPOSINRHSCOMP
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      TYPE_PARAL = KEEP(46)
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      I_AM_SLAVE = MYID .ne. 0 .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      BUF_EFFSIZE = 0
      BUF_MAXSIZE = max(min(BUF_MAXREF,int(2000000/NRHS)),2000)
      ALLOCATE (BUF_INDX(BUF_MAXSIZE),
     &          BUF_RHS(NRHS,BUF_MAXSIZE),
     &          stat=allocok)
      IF (allocok .GT. 0) THEN
        INFO(1)=-13
        INFO(2)=BUF_MAXSIZE*(NRHS+1)
      ENDIF
      CALL MUMPS_276(ICNTL, INFO, COMM, MYID )
      IF (INFO(1).LT.0) RETURN
      IF (MYID.EQ.MASTER) THEN
        ENTRIES_2_PROCESS = N - KEEP(89)
        DO WHILE ( ENTRIES_2_PROCESS .NE. 0)
          CALL MPI_RECV( BUF_INDX, BUF_MAXSIZE, MPI_INTEGER,
     &                 MPI_ANY_SOURCE,
     &                 ScatterRhsI, COMM, STATUS, IERR )
          CALL MPI_GET_COUNT( STATUS, MPI_INTEGER, BUF_EFFSIZE, IERR )
          PROC_WHO_ASKS = STATUS(MPI_SOURCE)
          DO I = 1, BUF_EFFSIZE
              INDX = BUF_INDX( I )
            DO K = 1, NRHS
              BUF_RHS( K, I ) = RHS( INDX, K )
              RHS( BUF_INDX(I), K ) = ZERO
            ENDDO
          ENDDO
          CALL MPI_SEND( BUF_RHS, NRHS*BUF_EFFSIZE,
     &                   MPI_REAL, PROC_WHO_ASKS,
     &                   ScatterRhsR, COMM, IERR)
          ENTRIES_2_PROCESS = ENTRIES_2_PROCESS - BUF_EFFSIZE
        ENDDO
        BUF_EFFSIZE= 0  
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (BUILD_POSINRHSCOMP) THEN
           IPOSINRHSCOMP = 1     
           POSINRHSCOMP = -9678  
        ENDIF
        IF (MYID.NE.MASTER) RHS = ZERO
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == MUMPS_275(PROCNODE_STEPS(ISTEP),
     &          NSLAVES)) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                  IPOS = PTRIST(ISTEP) 
                  LIELL = IW(IPOS+3+KEEP(IXSZ))
                  NPIV = LIELL
                  IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
              ELSE
                  IPOS = PTRIST(ISTEP) + 2+ KEEP(IXSZ)
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
              END IF
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
              ELSE
                   J1=IPOS+1+LIELL
              END IF
              IF (BUILD_POSINRHSCOMP) THEN
                 POSINRHSCOMP(ISTEP) = IPOSINRHSCOMP
                 IPOSINRHSCOMP       = IPOSINRHSCOMP + NPIV
              ENDIF
              IF (MYID.NE.MASTER) THEN
                DO JJ=J1,J1+NPIV-1
                  BUF_EFFSIZE = BUF_EFFSIZE + 1
                  BUF_INDX(BUF_EFFSIZE) = IW(JJ)
                  IF (BUF_EFFSIZE + 1 .GT. BUF_MAXSIZE) THEN
                   CALL SMUMPS_640()
                  ENDIF
                ENDDO
              ENDIF
          ENDIF
        ENDDO
        IF ( BUF_EFFSIZE .NE. 0 .AND. MYID.NE.MASTER ) 
     &              CALL SMUMPS_640()
      ENDIF
      DEALLOCATE (BUF_INDX, BUF_RHS)
      RETURN
      CONTAINS
                  SUBROUTINE SMUMPS_640()
                  CALL MPI_SEND(BUF_INDX, BUF_EFFSIZE, MPI_INTEGER,
     &            MASTER, ScatterRhsI, COMM, IERR )
                  CALL MPI_RECV(BUF_RHS, BUF_EFFSIZE*NRHS,
     &                 MPI_REAL,
     &                 MASTER,
     &                 ScatterRhsR, COMM, STATUS, IERR )
                  DO I = 1, BUF_EFFSIZE
                    INDX = BUF_INDX(I)
                    DO K = 1, NRHS
                      RHS( INDX, K ) = BUF_RHS( K, I )
                    ENDDO
                  ENDDO
                  BUF_EFFSIZE = 0
                  RETURN
                  END SUBROUTINE SMUMPS_640
      END SUBROUTINE SMUMPS_638
      SUBROUTINE SMUMPS_639
     &           (NSLAVES, N, MYID_NODES,
     &           PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     &           POSINRHSCOMP, POSINRHSCOMP_N, LPIRC_N, MTYPE,
     &           WHAT )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID_NODES, LIW
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N), POSINRHSCOMP(KEEP(28))
      INTEGER LPIRC_N, WHAT, MTYPE
      INTEGER POSINRHSCOMP_N(LPIRC_N)
      INTEGER ISTEP
      INTEGER NPIV
      INTEGER SK38, SK20, IPOS, LIELL
      INTEGER JJ, J1
      INTEGER IPOSINRHSCOMP
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      IF (WHAT .NE. 0.AND. WHAT.NE.1.AND.WHAT.NE.2) THEN
        WRITE(*,*) "Internal error in SMUMPS_639"
        CALL MUMPS_ABORT()
      ENDIF
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN 
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      IPOSINRHSCOMP   = 1     
      POSINRHSCOMP = -9678  
      IF (WHAT .NE. 0) THEN
        POSINRHSCOMP_N = 0 
      ENDIF
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_275(PROCNODE_STEPS(ISTEP),
     &     NSLAVES)) THEN
           IPOS = PTRIST(ISTEP)
           NPIV = IW(IPOS+3+KEEP(IXSZ))
           POSINRHSCOMP(ISTEP) = IPOSINRHSCOMP
           IF (WHAT .NE. 0) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                  IPOS = PTRIST(ISTEP)
                  LIELL = IW(IPOS+3+KEEP(IXSZ))
                  NPIV = LIELL
                  IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
              ELSE
                  IPOS = PTRIST(ISTEP) + 2+ KEEP(IXSZ)
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
              ENDIF
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
              ELSE
                   J1=IPOS+1+LIELL
              END IF
              DO JJ = J1, J1+NPIV-1
                POSINRHSCOMP_N(IW(JJ)) = IPOSINRHSCOMP+JJ-J1
              END DO
           ENDIF
           IPOSINRHSCOMP       = IPOSINRHSCOMP + NPIV
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_639
      SUBROUTINE SMUMPS_248(N, A, LA, IW, LIW, WCB, LWCB,
     &    RHS, LRHS, NRHS,
     &    PTRICB, IWCB, LIWCB, 
     &    RHSCOMP, LRHSCOMP, POSINRHSCOMP, BUILD_POSINRHSCOMP,
     &    NE_STEPS, NA, LNA, STEP,
     &    FRERE, DAD, FILS,
     &    NSTK_S, IPOOL, LPOOL, PTRIST, PTRFAC, MYLEAF, INFO,
     &    KEEP,KEEP8,
     &    PROCNODE_STEPS,
     &    SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,
     &    RHS_ROOT, LRHS_ROOT, MTYPE, 
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
      USE SMUMPS_OOC
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER(8) :: LA
      INTEGER N, LIW, LWCB, LPOOL, LIWCB, LNA
      INTEGER SLAVEF, MYLEAF, COMM, MYID
      INTEGER INFO( 40 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER LRHS, NRHS
      REAL A( LA ), RHS( LRHS, NRHS ), WCB( LWCB )
      INTEGER LRHS_ROOT
      REAL RHS_ROOT( LRHS_ROOT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NA( LNA ), NE_STEPS( KEEP(28) )
      INTEGER STEP( N ), FRERE( KEEP(28) ), FILS( N ),
     &        DAD( KEEP(28) )
      INTEGER NSTK_S(KEEP(28)), IPOOL( LPOOL )
      INTEGER PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PTRICB( KEEP(28) ) 
      INTEGER IW( LIW ), IWCB( LIWCB )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER POSINRHSCOMP(KEEP(28)), LRHSCOMP 
      LOGICAL BUILD_POSINRHSCOMP
      REAL RHSCOMP( LRHSCOMP, NRHS )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGTAG, MSGSOU, DUMMY(1)
      LOGICAL FLAG
      INTEGER NBFIN, MYROOT
      INTEGER POSIWCB,POSWCB,PLEFTWCB
      INTEGER INODE
      INTEGER RHSCOMPFREEPOS
      INTEGER I
      INTEGER III, NBROOT,LEAF
      LOGICAL BLOQ
      EXTERNAL MUMPS_275
      INTEGER MUMPS_275
      POSIWCB = LIWCB
      POSWCB  = LWCB
      PLEFTWCB= 1
      IF (BUILD_POSINRHSCOMP) RHSCOMPFREEPOS= 1
      DO I = 1, KEEP(28)
        NSTK_S(I)   = NE_STEPS(I)
      ENDDO
      PTRICB = 0
      CALL MUMPS_362(N, LEAF, NBROOT, MYROOT, MYID,
     &     SLAVEF, NA, LNA, KEEP,KEEP8, STEP,
     &     PROCNODE_STEPS, IPOOL, LPOOL)
      NBFIN = SLAVEF
      IF ( MYROOT .EQ. 0 ) THEN
        NBFIN = NBFIN - 1
        DUMMY(1) = 1
        CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &       RACINE_SOLVE, SLAVEF)
      END IF
      MYLEAF = LEAF - 1
      III    = 1
   50 CONTINUE
      IF (SLAVEF .EQ. 1) THEN
         CALL SMUMPS_574
     &        ( IPOOL(1), LPOOL, III, LEAF, INODE,
     &          KEEP(208) )
        GOTO 60
      ENDIF
      BLOQ = ( ( III .EQ. LEAF )
     &     )
      CALL SMUMPS_303( BLOQ, FLAG,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, NRHS, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &     IWCB, LIWCB,
     &     WCB, LWCB, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB, INFO, KEEP,KEEP8, STEP,
     &     PROCNODE_STEPS,
     &     RHS, LRHS
     &     )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      IF (.not. FLAG) THEN
         IF (III .NE. LEAF) THEN
            CALL SMUMPS_574
     &           (IPOOL(1), LPOOL, III, LEAF, INODE,
     &           KEEP(208) )
            GOTO 60
         ENDIF                  
      ENDIF                     
      GOTO 50
 60   CONTINUE
      CALL SMUMPS_302( INODE, BUFR, LBUFR, LBUFR_BYTES,
     &        MSGTAG, MSGSOU, MYID, SLAVEF, COMM,  N,
     &        IPOOL, LPOOL, III, LEAF, NBFIN, NSTK_S,
     &        IWCB, LIWCB, WCB, LWCB, A, LA,
     &        IW, LIW, RHS, LRHS, NRHS, 
     &        POSWCB, PLEFTWCB, POSIWCB,
     &        PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &        FILS, STEP, FRERE, DAD,
     &        MYROOT, INFO, KEEP,KEEP8, RHS_ROOT, MTYPE, 
     &        RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     &        RHSCOMPFREEPOS, BUILD_POSINRHSCOMP,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &     )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      GOTO 50
  260 CONTINUE
      CALL SMUMPS_150( MYID,COMM,BUFR,
     &                            LBUFR,LBUFR_BYTES )
      RETURN
      END SUBROUTINE SMUMPS_248
      RECURSIVE SUBROUTINE SMUMPS_323
     &     ( BUFR, LBUFR, LBUFR_BYTES,
     &     MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &     N, NRHS, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     &     PTRFAC, IWCB, LIWCB,
     &     WCB, LWCB, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB,
     &     INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS, 
     &     RHS, LRHS 
     &     )
      USE SMUMPS_OOC 
      USE SMUMPS_COMM_BUFFER 
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ),  NSTK_S( N )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      INTEGER PTRICB(KEEP(28)),PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      REAL WCB( LWCB ), A( LA )
      INTEGER LRHS
      REAL RHS(LRHS, NRHS)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, K, JJ
      INTEGER FINODE, FPERE, LONG, NCB, POSITION, NCV, NPIV
      INTEGER PTRX, PTRY, PDEST, I
      INTEGER(8) :: APOS
      LOGICAL DUMMY
      LOGICAL FLAG
      EXTERNAL MUMPS_275
      INTEGER  MUMPS_275
      REAL ALPHA, ONE
      PARAMETER (ONE = 1.0E0, ALPHA=-1.0E0)
      INCLUDE 'mumps_headers.h'
      IF ( MSGTAG .EQ. RACINE_SOLVE ) THEN
         NBFIN = NBFIN - 1
         IF ( NBFIN .eq. 0 ) GOTO 270
      ELSE  IF (MSGTAG .EQ. ContVec ) THEN
         POSITION = 0
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FINODE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FPERE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        NCB, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        LONG, 1, MPI_INTEGER, COMM, IERR )
          IF ( NCB .eq. 0 ) THEN
             PTRICB(STEP(FINODE)) = -1
             NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             END IF
          ELSE
             IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
                PTRICB(STEP(FINODE)) = NCB + 1
             END IF
             IF ( ( POSIWCB - LONG ) .LT. 0 ) THEN
                INFO( 1 ) = -14
                INFO( 2 ) = LONG
                GOTO 260
             END IF
             IF ( POSWCB - PLEFTWCB + 1 .LT. LONG * NRHS) THEN
                INFO( 1 ) = -11
                INFO( 2 ) = PLEFTWCB - POSWCB - 1 + LONG * NRHS
                GOTO 260
             END IF
             IF (LONG .GT. 0) THEN
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               IWCB( 1 ),
     &               LONG, MPI_INTEGER, COMM, IERR )
                DO K = 1, NRHS
                   CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                  WCB( PLEFTWCB ),
     &                  LONG, MPI_REAL, COMM, IERR )
                   DO I = 1, LONG
                      RHS(IWCB(I),K) = RHS(IWCB(I),K) +WCB(PLEFTWCB+I-1)
                   ENDDO
                END DO
                PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - LONG
             ENDIF
             IF ( PTRICB(STEP(FINODE)) == 1 ) THEN
                NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             END IF
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             ENDIF
          END IF
       ELSEIF ( MSGTAG .EQ. Master2Slave ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NCV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NPIV, 1, MPI_INTEGER, COMM, IERR )
          PTRY = PLEFTWCB
          PTRX = PLEFTWCB + NCV * NRHS
          PLEFTWCB = PLEFTWCB + (NPIV + NCV) * NRHS
          IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
             INFO(1) = -11
             INFO(2) = -POSWCB + PLEFTWCB -1
             GO TO 260
          END IF
          DO K=1, NRHS
             CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &            WCB( PTRY + (K-1) * NCV ), NCV,
     &            MPI_REAL, COMM, IERR )
          ENDDO
          IF ( NPIV .GT. 0 ) THEN
             DO K=1, NRHS
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               WCB( PTRX + (K-1)*NPIV ), NPIV,
     &               MPI_REAL, COMM, IERR )
             END DO
          END IF
          IF (KEEP(201).GT.0) THEN
             CALL SMUMPS_643(
     &            FINODE,PTRFAC,KEEP,A,LA,STEP,
     &            KEEP8,N,DUMMY,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          APOS = PTRFAC(STEP(FINODE))
          IF (KEEP(201).EQ.1) THEN
             IF ( NRHS == 1 ) THEN
                CALL sgemv( 'N', NCV, NPIV, ALPHA, A(APOS), NCV,
     &               WCB( PTRX ), 1, ONE,
     &               WCB( PTRY ), 1 )
             ELSE
                CALL sgemm( 'N', 'N', NCV, NRHS, NPIV, ALPHA,
     &               A(APOS), NCV,
     &               WCB( PTRX), NPIV, ONE,
     &               WCB( PTRY), NCV )
             ENDIF
          ELSE                  
             IF ( NRHS == 1 ) THEN
                CALL sgemv( 'T', NPIV, NCV, ALPHA, A(APOS), NPIV,
     &               WCB( PTRX ), 1, ONE,
     &               WCB( PTRY ), 1 )
             ELSE
                CALL sgemm( 'T', 'N', NCV, NRHS, NPIV, ALPHA,
     &               A(APOS), NPIV,
     &               WCB( PTRX), NPIV, ONE,
     &               WCB( PTRY), NCV )
             ENDIF
          ENDIF
          IF (KEEP(201).GT.0) THEN
             CALL SMUMPS_598(FINODE,PTRFAC,
     &            KEEP(28),A,LA,.TRUE.,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          PLEFTWCB = PLEFTWCB - NPIV * NRHS
          PDEST = MUMPS_275( PROCNODE_STEPS(STEP(FPERE)),
     &                            SLAVEF )
          IF ( PDEST .EQ. MYID ) THEN
             IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
                NCB = IW( PTRIST(STEP(FINODE)) + 2 + KEEP(IXSZ) )
                PTRICB(STEP(FINODE)) = NCB + 1
             END IF
             DO I = 1, NCV
                JJ=IW(PTRIST(STEP(FINODE))+3+I+ KEEP(IXSZ) )
                DO K=1, NRHS
                   RHS(JJ,K)= RHS(JJ,K) + WCB(PTRY+I-1+(K-1)*NCV)
                ENDDO
             END DO
             PTRICB(STEP(FINODE)) =
     &            PTRICB(STEP(FINODE)) - NCV
             IF ( PTRICB( STEP( FINODE ) ) == 1 ) THEN
                NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             END IF
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'INTERNAL Error 41r: Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             ENDIF
          ELSE
 210         CONTINUE
             CALL SMUMPS_78( NRHS, FINODE, FPERE,
     &            IW(PTRIST(STEP( FINODE )) + 2 + KEEP(IXSZ) ), NCV,NCV,
     &            IW(PTRIST(STEP(FINODE))+4+ KEEP(IXSZ) ),
     &            WCB( PTRY ), PDEST, ContVec, COMM, IERR )
             IF ( IERR .EQ. -1 ) THEN
                CALL SMUMPS_303( .FALSE., FLAG,
     &               BUFR, LBUFR, LBUFR_BYTES,
     &               MYID, SLAVEF, COMM,
     &               N, NRHS, IPOOL, LPOOL, III, LEAF,
     &               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &               IWCB, LIWCB,
     &               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &               PTRICB, INFO, KEEP,KEEP8, STEP,
     &               PROCNODE_STEPS, 
     &               RHS, LRHS
     &               )
                IF ( INFO( 1 )  .LT. 0 )  GOTO 270
                GOTO 210
             ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &               NCV * KEEP( 35 )
                GOTO 260
             ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &               NCV * KEEP( 35 )
             END IF
          END IF
          PLEFTWCB = PLEFTWCB - NCV * NRHS
       ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GOTO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     &         (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
       ELSE
          INFO(1)=-100
          INFO(2)=MSGTAG
          GO TO 260
       ENDIF
       GO TO 270
 260   CONTINUE
       CALL SMUMPS_44( MYID, SLAVEF, COMM )
 270   CONTINUE
       RETURN
       END SUBROUTINE SMUMPS_323
      SUBROUTINE SMUMPS_302( INODE,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &     N, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S,
     &     IWCB, LIWCB,
     &     WCB, LWCB, A, LA, IW, LIW,
     &     RHS, LRHS, NRHS, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &     FILS, STEP, FRERE, DAD,
     &     MYROOT,
     &     INFO, KEEP,KEEP8, RHS_ROOT, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     &     RHSCOMPFREEPOS, BUILD_POSINRHSCOMP,
     &     
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &     
     &            )
      USE SMUMPS_OOC
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER INODE, LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIWCB, LWCB, LIW, POSWCB, PLEFTWCB, POSIWCB
      INTEGER(8) :: LA
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER MYROOT
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ), NSTK_S(KEEP(28))
      INTEGER IWCB( LIWCB ), IW( LIW )
      INTEGER LRHS, NRHS
      REAL WCB( LWCB ), A( LA )
      REAL RHS(LRHS, NRHS ), RHS_ROOT( * )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER FILS( N ), STEP( N ), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &     TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER POSINRHSCOMP(KEEP(28)), LRHSCOMP, RHSCOMPFREEPOS
      REAL RHSCOMP(LRHSCOMP, NRHS)
      LOGICAL BUILD_POSINRHSCOMP
      EXTERNAL sgemv, strsv, sgemm, strsm, MUMPS_275
      INTEGER MUMPS_275
      REAL ALPHA,ONE,ZERO
      PARAMETER (ZERO=0.0E0, ONE = 1.0E0, ALPHA=-1.0E0)
      INTEGER(8) :: APOS, APOS1, APOS2, APOSOFF
      INTEGER I, J, K, IPOS, NSLAVES, J1, J2, J3, FPERE, NPIV, NCB,
     &     IERR, IFR_ini,
     &     IFR, LIELL, JJ,
     &     NELIM, PLEFT, PCB_COURANT, PPIV_COURANT
      INTEGER IPOSINRHSCOMP
      INTEGER Effective_CB_Size, NUPDATE, ISLAVE, PDEST, FirstIndex
      LOGICAL FLAG, OMP_FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER POSWCB1,POSWCB2
      INTEGER(8) :: APOSDEB
      INTEGER TempNROW, TempNCOL, PANEL_SIZE, LIWFAC, 
     &     JFIN, NBJ, NUPDATE_PANEL,
     &     PPIV_PANEL, PCB_PANEL, NBK, TYPEF
      INTEGER LD_WCBPIV         
      INTEGER LD_WCBCB          
      INTEGER LDAJ, LDAJ_FIRST_PANEL
      INTEGER TMP_NBPANELS,
     &     I_PIVRPTR, I_PIVR, IPANEL
      LOGICAL MUST_BE_PERMUTED
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY( 1 )
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq.KEEP( 20 ) ) THEN
         LIELL = IW( PTRIST( STEP(INODE)) + 3 + KEEP(IXSZ))
         NPIV  = LIELL
         NELIM = 0
         NSLAVES = 0
         IPOS = PTRIST( STEP(INODE)) + 5 + KEEP(IXSZ)
      ELSE
        IPOS = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        NELIM = IW(IPOS-1)
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ) )
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        IF (KEEP(201).GT.0) THEN
           CALL SMUMPS_643(
     &          INODE,PTRFAC,KEEP,A,LA,STEP,
     &          KEEP8,N,MUST_BE_PERMUTED,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
           CALL SMUMPS_755(
     &                 IW(IPOS+1+2*LIELL+1+NSLAVES),
     &                 MUST_BE_PERMUTED )
           ENDIF
        ENDIF                     
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
        IPOS = IPOS + 1 + NSLAVES
      END IF
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
         J1 = IPOS + 1
         J2 = IPOS + LIELL
         J3 = IPOS + NPIV
      ELSE
         J1 = IPOS + LIELL + 1
         J2 = IPOS + 2 * LIELL
         J3 = IPOS + LIELL + NPIV
      END IF
      NCB = LIELL-NPIV
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq. KEEP(20) ) THEN
         IFR = 0
         DO JJ = J1, J3
            J = IW( JJ )
            IFR = IFR + 1
            DO K=1,NRHS
               RHS_ROOT(IFR+NPIV*(K-1)) = RHS(J,K) 
            END DO
         END DO
         IF ( NPIV .LT. LIELL ) THEN
            WRITE(*,*) ' Internal error in SOLVE_NODE for Root node'
            CALL MUMPS_ABORT()
         END IF
         MYROOT = MYROOT - 1
         IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
               DUMMY (1) = 1
               CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     &              COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
         END IF
         GO TO 270
      END IF
      APOS = PTRFAC(STEP(INODE))
      IF (KEEP(201).EQ.1) THEN  
        IF (MTYPE.EQ.1) THEN
            IF ((MTYPE.EQ.1).AND.NSLAVES.NE.0) THEN
              TempNROW= NPIV+NELIM
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ELSE
              TempNROW= LIELL
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ENDIF
            TYPEF=TYPEF_L
        ELSE 
            TempNCOL= LIELL
            TempNROW= NPIV
            LDAJ_FIRST_PANEL=TempNCOL
            TYPEF= TYPEF_U
        ENDIF
        LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
        PANEL_SIZE = SMUMPS_690( LDAJ_FIRST_PANEL )
      ENDIF                     
      PLEFT    = PLEFTWCB
      PPIV_COURANT = PLEFTWCB
      PLEFTWCB = PLEFTWCB + LIELL * NRHS
      IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
         INFO(1) = -11
         INFO(2) = PLEFTWCB - POSWCB - 1
         GO TO 260
      END IF
      IF (KEEP(201).EQ.1) THEN  
         LD_WCBPIV = LIELL
         LD_WCBCB  = LIELL
         PCB_COURANT = PPIV_COURANT + NPIV
         DO K=1, NRHS
            IFR = PPIV_COURANT + (K-1)*LIELL - 1
            DO JJ = J1, J3
               J = IW(JJ)
               IFR = IFR + 1
               WCB(IFR) = RHS(J,K) 
            ENDDO
            IF (NCB.GT.0) THEN
               DO JJ = J3+1, J2
                  J = IW(JJ)
                  IFR = IFR + 1
                  WCB(IFR) = RHS(J,K) 
                  RHS (J,K) = ZERO
               ENDDO
            ENDIF
         END DO
      ELSE                      
         LD_WCBPIV = NPIV
         LD_WCBCB  = NCB
         PCB_COURANT = PPIV_COURANT + NPIV*NRHS
         IFR = PPIV_COURANT - 1
         OMP_FLAG = NRHS.GT.4
         IFR_ini = IFR
         DO 130 JJ = J1, J3
            J = IW(JJ)
            IFR = IFR_ini + (JJ-J1) + 1
            DO K=1, NRHS
               WCB(IFR+(K-1)*NPIV) = RHS(J,K) 
            END DO
 130     CONTINUE
         IFR = PCB_COURANT - 1
         IF (NPIV .LT. LIELL) THEN
            IFR_ini = IFR
            DO 140 JJ = J3 + 1, J2
               J = IW(JJ)
               IFR = IFR_ini + (JJ-J3)
               DO K=1, NRHS
                  WCB(IFR+(K-1)*NCB) = RHS(J,K)
                  RHS(J,K)=ZERO
               ENDDO
 140        CONTINUE
         ENDIF
      ENDIF                     
      IF ( NPIV .NE. 0 ) THEN
         IF (KEEP(201).EQ.1) THEN 
        APOSDEB = APOS
        J = 1
        IPANEL = 0
  10    CONTINUE
          IPANEL = IPANEL + 1
          JFIN    = min(J+PANEL_SIZE-1, NPIV)
          IF (IW(IPOS+ LIELL + JFIN) < 0) THEN
            JFIN=JFIN+1
          ENDIF
          NBJ     = JFIN-J+1
          LDAJ    = LDAJ_FIRST_PANEL-J+1 
          IF ( (KEEP(50).NE.1).AND. MUST_BE_PERMUTED ) THEN
           CALL SMUMPS_667(TYPEF, TMP_NBPANELS,
     &            I_PIVRPTR, I_PIVR, IPOS+1+2*LIELL, IW, LIW)
               IF (NPIV.EQ.(IW(I_PIVRPTR+IPANEL-1)-1)) THEN
                  MUST_BE_PERMUTED=.FALSE. 
               ELSE
                  CALL SMUMPS_698(
     &                 IW( I_PIVR+ IW(I_PIVRPTR+IPANEL-1)-
     &                 IW(I_PIVRPTR)), 
     &                 NPIV-IW(I_PIVRPTR+IPANEL-1)+1, 
     &                 IW(I_PIVRPTR+IPANEL-1)-1, 
     &                          
     &                 A(APOSDEB),
     &                 LDAJ, NBJ, J-1 ) 
               ENDIF
            ENDIF 
            NUPDATE_PANEL = LDAJ - NBJ
            PPIV_PANEL = PPIV_COURANT+J-1
            PCB_PANEL  = PPIV_PANEL+NBJ
            APOS1 = APOSDEB+int(NBJ,8)
            IF  (MTYPE.EQ.1) THEN
               IF ( NRHS == 1 ) THEN
                  CALL strsv( 'L', 'N', 'U', NBJ, A(APOSDEB), LDAJ, 
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL sgemv('N', NUPDATE_PANEL,NBJ,ALPHA, A(APOS1),
     &                    LDAJ,  WCB(PPIV_PANEL), 1, ONE,
     &                    WCB(PCB_PANEL), 1)
                  ENDIF
               ELSE
                  CALL strsm( 'L','L','N','U', NBJ, NRHS, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL sgemm('N', 'N', NUPDATE_PANEL, NRHS, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &                    WCB(PCB_PANEL), LIELL)
                  ENDIF
               ENDIF
            ELSE
               IF (NRHS == 1) THEN
                  CALL strsv( 'L', 'N', 'N', NBJ, A(APOSDEB), LDAJ,
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL sgemv('N',NUPDATE_PANEL, NBJ, ALPHA, A(APOS1),
     &                    LDAJ, WCB(PPIV_PANEL), 1,
     &                    ONE, WCB(PCB_PANEL), 1 )
                  ENDIF
               ELSE
                  CALL strsm('L','L','N','N',NBJ, NRHS, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL)
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL sgemm('N', 'N', NUPDATE_PANEL, NRHS, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &             WCB(PCB_PANEL), LIELL)
                  ENDIF
               ENDIF
            ENDIF
            APOSDEB = APOSDEB+int(LDAJ,8)*int(NBJ,8)
            J=JFIN+1
            IF ( J .LE. NPIV ) GOTO 10
         ELSE                   
            IF (KEEP(50).NE.0) THEN
               IF ( NRHS == 1 ) THEN
                  CALL strsv( 'U', 'T', 'U', NPIV, A(APOS), NPIV,
     &                   WCB(PPIV_COURANT), 1 )
               ELSE
                  CALL strsm( 'L','U','T','U', NPIV, NRHS, ONE,
     &                   A(APOS), NPIV, WCB(PPIV_COURANT),
     &                   NPIV )
               ENDIF
            ELSE
               IF ( MTYPE .eq. 1 ) THEN
                  IF ( NRHS == 1)  THEN
                     CALL strsv( 'U', 'T', 'U', NPIV, A(APOS), LIELL, 
     &                    WCB(PPIV_COURANT), 1 )
                  ELSE
                     CALL strsm( 'L','U','T','U', NPIV, NRHS, ONE,
     &                    A(APOS), LIELL, WCB(PPIV_COURANT),
     &                    NPIV )
                  ENDIF
               ELSE
                  IF (NRHS == 1) THEN
                     CALL strsv( 'L', 'N', 'N', NPIV, A(APOS), LIELL,
     &                    WCB(PPIV_COURANT), 1 )
                  ELSE
                     CALL strsm('L','L','N','N',NPIV, NRHS, ONE,
     &                    A(APOS), LIELL, WCB(PPIV_COURANT),
     &                    NPIV)
                  ENDIF
               END IF
            END IF              
         END IF                 
      END IF                    
      NCB   = LIELL - NPIV
      IF ( MTYPE .EQ. 1 ) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            APOS1 = APOS  + int(NPIV,8) * int(LIELL,8)
         ELSE
            APOS1 = APOS + int(NPIV,8) * int(NPIV,8)
         END IF
         IF ( NSLAVES .EQ. 0 .OR. NPIV .eq. 0 ) THEN
            NUPDATE = NCB
         ELSE
            NUPDATE = NELIM
         END IF
      ELSE
         APOS1 = APOS + int(NPIV,8)
         NUPDATE = NCB
      END IF
      IF (KEEP(201).NE.1) THEN  
         IF ( NPIV .NE. 0 .AND. NUPDATE.NE.0 ) THEN
            IF ( MTYPE .eq. 1 ) THEN
               IF ( NRHS == 1 ) THEN
                  CALL sgemv('T', NPIV, NUPDATE, ALPHA, A(APOS1),
     &            NPIV,  WCB(PPIV_COURANT), 1, ONE,
     &            WCB(PCB_COURANT), 1)
               ELSE
                  CALL sgemm('T', 'N', NUPDATE, NRHS, NPIV, ALPHA,
     &            A(APOS1), NPIV, WCB(PPIV_COURANT), NPIV, ONE,
     &            WCB(PCB_COURANT), NCB)
               END IF
            ELSE                
               IF ( NRHS == 1 ) THEN
                  CALL sgemv('N',NUPDATE, NPIV, ALPHA, A(APOS1),
     &                 LIELL, WCB(PPIV_COURANT), 1,
     &                 ONE, WCB(PCB_COURANT), 1 )
               ELSE
                  CALL sgemm('N', 'N', NUPDATE, NRHS, NPIV, ALPHA,
     &                 A(APOS1), LIELL, WCB(PPIV_COURANT), NPIV, ONE,
     &                 WCB(PCB_COURANT), NCB)
               END IF
            END IF
         END IF
      END IF                    
      IF (BUILD_POSINRHSCOMP) THEN
         POSINRHSCOMP(STEP(INODE)) =  RHSCOMPFREEPOS
         RHSCOMPFREEPOS            = RHSCOMPFREEPOS + NPIV
      ENDIF
      IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
      IF ( KEEP(50) .eq. 0 ) THEN
         DO K=1,NRHS
            IFR =  PPIV_COURANT + (K-1)*LD_WCBPIV
            RHSCOMP(IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1, K) =
     &           WCB(IFR:IFR+NPIV-1)
         ENDDO
      ELSE
         IFR = PPIV_COURANT - 1
         IF (KEEP(201).EQ.1) THEN 
            LDAJ = TempNROW  
         ELSE                
            LDAJ = NPIV 
         ENDIF
         APOS1 = APOS
         JJ    = J1
         IF (KEEP(201).EQ.1) THEN
            NBK   = 0           
         ENDIF
         DO 
            IF(JJ .GT. J3) EXIT
            IFR = IFR + 1
            IF(IW(JJ+LIELL) .GT. 0) THEN
               DO K=1, NRHS
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1 , K ) = 
     &                 WCB( IFR+(K-1)*LD_WCBPIV ) * A( APOS1 )
               END DO
            IF (KEEP(201).EQ.1) THEN
              NBK = NBK+1
              IF (NBK.EQ.PANEL_SIZE) THEN
                NBK = 0
                LDAJ = LDAJ - PANEL_SIZE
              ENDIF
            ENDIF
            APOS1 = APOS1 + int(LDAJ + 1,8)
            JJ = JJ+1
         ELSE
            IF (KEEP(201).EQ.1) THEN
              NBK = NBK+1
            ENDIF
            APOS2 = APOS1+int(LDAJ+1,8)
            IF (KEEP(201).EQ.1) THEN
              APOSOFF = APOS1+int(LDAJ,8)
            ELSE
              APOSOFF=APOS1+1_8
            ENDIF
               DO K=1, NRHS
                  POSWCB1 = IFR+(K-1)*LD_WCBPIV
                  POSWCB2 = POSWCB1+1
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = WCB(POSWCB1)*A(APOS1)
     &                 + WCB(POSWCB2)*A(APOSOFF)
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1+1,K) = 
     &                 WCB(POSWCB1)*A(APOSOFF)
     &                 + WCB(POSWCB2)*A(APOS2)
               END DO
               IF (KEEP(201).EQ.1) THEN
                  NBK = NBK+1
                  IF (NBK.GE.PANEL_SIZE) THEN
                     LDAJ = LDAJ - NBK
                     NBK = 0
                  ENDIF
               ENDIF
               APOS1 = APOS2 + int(LDAJ + 1,8)
               JJ = JJ+2
               IFR = IFR+1
            ENDIF
         ENDDO
      END IF
      IF (KEEP(201).GT.0) THEN
         CALL SMUMPS_598(INODE,PTRFAC,KEEP(28),
     &        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 260
         ENDIF
      END IF
      FPERE = DAD(STEP(INODE))
      IF ( FPERE .EQ. 0 ) THEN
         MYROOT = MYROOT - 1
         PLEFTWCB = PLEFTWCB - LIELL *NRHS
         IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
               DUMMY (1) = 1
               CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     &             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
         END IF
         GO TO 270
      ENDIF
      IF ( NUPDATE .NE. 0 .OR. NCB.eq.0 ) THEN
         IF (MUMPS_275(PROCNODE_STEPS(STEP(FPERE)),
     &        SLAVEF) .EQ. MYID) THEN
            IF ( NCB .ne. 0 ) THEN
               PTRICB(STEP(INODE)) = NCB + 1
               DO 190 I = 1, NUPDATE
                  DO K=1, NRHS
                     RHS( IW(J3 + I), K ) = RHS( IW(J3 + I), K )
     &                    + WCB(PCB_COURANT + I-1 +(K-1)*LD_WCBCB)
                  ENDDO
 190           CONTINUE
               PTRICB(STEP( INODE )) = PTRICB(STEP( INODE )) - NUPDATE
               IF ( PTRICB(STEP(INODE)) == 1 ) THEN
                  NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
                  IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
                        IPOOL( LEAF ) = FPERE
                     LEAF = LEAF + 1
                  ENDIF
               END IF
            ELSE
               PTRICB(STEP( INODE )) = -1
               NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
               IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN 
                     IPOOL( LEAF ) = FPERE 
                  LEAF = LEAF + 1
               ENDIF            
            ENDIF
         ELSE
 210        CONTINUE
            CALL SMUMPS_78( NRHS, INODE, FPERE, NCB, LD_WCBCB,
     &           NUPDATE,
     &           IW( J3 + 1 ), WCB( PCB_COURANT ),
     &           MUMPS_275(PROCNODE_STEPS(STEP(FPERE)), SLAVEF),
     &           ContVec,
     &           COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
               CALL SMUMPS_303( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, III, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, STEP,
     &              PROCNODE_STEPS, 
     &              RHS, LRHS 
     &              )
               IF ( INFO( 1 )  .LT. 0 )  GOTO 270
               GOTO 210
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -17
               INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &              ( NUPDATE + 3 ) * KEEP( 34 )
               GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &              ( NUPDATE + 3 ) * KEEP( 34 )
               GOTO 260
            END IF
         ENDIF
      END IF
      IF ( NSLAVES .NE. 0 .AND. MTYPE .eq. 1
     &     .and. NPIV .NE. 0 ) THEN
         DO ISLAVE = 1, NSLAVES
            PDEST = IW( PTRIST(STEP(INODE)) + 5 + ISLAVE +KEEP(IXSZ))
            CALL MUMPS_49( 
     &           KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &           ISLAVE, NCB - NELIM, 
     &           NSLAVES, 
     &           Effective_CB_Size, FirstIndex )
 222        CALL SMUMPS_72( NRHS,
     &           INODE, FPERE,
     &           Effective_CB_Size, LD_WCBCB, LD_WCBPIV, NPIV,
     &           WCB( PCB_COURANT + NELIM + FirstIndex - 1 ),
     &           WCB( PPIV_COURANT ),
     &           PDEST, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
               CALL SMUMPS_303( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, III, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, STEP,
     &              PROCNODE_STEPS, 
     &              RHS, LRHS 
     &              )
               IF ( INFO( 1 )  .LT. 0 )  GOTO 270
               GOTO 222
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -17
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS*KEEP(35) +
     &              ( Effective_CB_Size + 4 ) * KEEP( 34 )
               GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS*KEEP(35) +
     &              ( Effective_CB_Size + 4 ) * KEEP( 34 )
               GOTO 260
            END IF
         END DO
      END IF
      PLEFTWCB = PLEFTWCB - LIELL*NRHS
 270  CONTINUE
      RETURN
 260  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_302
      RECURSIVE SUBROUTINE SMUMPS_303( BLOQ, FLAG,
     &           BUFR, LBUFR, LBUFR_BYTES,
     &           MYID, SLAVEF, COMM,
     &           N, NRHS, IPOOL, LPOOL, III, LEAF,
     &           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &           IWCB, LIWCB,
     &           WCB, LWCB, POSWCB,
     &           PLEFTWCB, POSIWCB,
     &           PTRICB, INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS,
     &           RHS, LRHS
     &            )
      IMPLICIT NONE
      LOGICAL BLOQ
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL)
      INTEGER NSTK_S( KEEP(28) )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      REAL WCB( LWCB ), A( LA )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LRHS
      REAL RHS(LRHS, NRHS)
      LOGICAL FLAG
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      INTEGER MSGSOU, MSGTAG, MSGLEN
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &                   FLAG, STATUS, IERR )
      END IF
      IF ( FLAG ) THEN
         MSGSOU = STATUS( MPI_SOURCE )
         MSGTAG = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL SMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED,
     &                  MSGSOU, MSGTAG, COMM, STATUS, IERR )
           CALL SMUMPS_323( BUFR, LBUFR, LBUFR_BYTES,
     &          MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &          N, NRHS, IPOOL, LPOOL, III, LEAF,
     &          NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &          IWCB, LIWCB,
     &          WCB, LWCB, POSWCB,
     &          PLEFTWCB, POSIWCB,
     &          PTRICB, INFO, KEEP,KEEP8, STEP,
     &          PROCNODE_STEPS, 
     &          RHS, LRHS 
     &          )
         END IF
      END IF
      RETURN
      END SUBROUTINE SMUMPS_303
      SUBROUTINE SMUMPS_249(N, A, LA, IW, LIW, W, LWC,
     &    RHS, LRHS, NRHS, 
     &    RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     &    PTRICB, PTRACB, IWCB, LIWW, W2, 
     &    NE_STEPS, NA, LNA, STEP,
     &    FRERE, DAD, FILS, IPOOL, LPOOL, PTRIST, PTRFAC, 
     &    MYLEAF, INFO, 
     &    PROCNODE_STEPS,
     &    SLAVEF, COMM,MYID, BUFR, LBUFR, LBUFR_BYTES,
     &    KEEP,KEEP8, RHS_ROOT, LRHS_ROOT, MTYPE, 
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, PANEL_POS, LPANEL_POS
     &    , TO_PROCESS, SIZE_TO_PROCESS
     &    )
      USE SMUMPS_OOC
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER(8) :: LA
      INTEGER N,LIW,LIWW,LWC,LPOOL,LNA
      INTEGER SLAVEF,MYLEAF,COMM,MYID
      INTEGER LPANEL_POS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER NA(LNA),NE_STEPS(KEEP(28))
      INTEGER IPOOL(LPOOL)
      INTEGER PANEL_POS(LPANEL_POS)
      INTEGER INFO(40)
      INTEGER PTRIST(KEEP(28)),
     &        PTRICB(KEEP(28)),PTRACB(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER LRHS, NRHS
      REAL A(LA), RHS(LRHS,NRHS), W(LWC)
      REAL W2(KEEP(133))
      INTEGER IW(LIW),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)),DAD(KEEP(28)),FILS(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER LRHSCOMP, POSINRHSCOMP(KEEP(28))
      REAL RHSCOMP(LRHSCOMP,NRHS)
      INTEGER LRHS_ROOT
      REAL RHS_ROOT( LRHS_ROOT )
      INTEGER, intent(in)           :: SIZE_TO_PROCESS
      LOGICAL, intent(in)           :: TO_PROCESS(SIZE_TO_PROCESS)
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      LOGICAL FLAG
      INTEGER POSIWCB,POSWCB,K
      INTEGER(8) :: APOS, IST
      INTEGER NPIV
      INTEGER IPOS,LIELL,NELIM,IFR,JJ,I
      INTEGER J1,J2,J,NCB,NBFINF
      INTEGER NBLEAF,INODE,NBROOT,NROOT,NBFILS
      INTEGER IN,IF,LONG,POOL_FIRST_POS,TMP
      INTEGER III,IIPOOL,MYLEAFE
      INTEGER NSLAVES
      REAL ALPHA,ONE,ZERO
      PARAMETER (ZERO=0.0E0, ONE = 1.0E0, ALPHA=-1.0E0)
      LOGICAL BLOQ,DEBUT
      INTEGER PROCDEST, DEST
      INTEGER POSINDICES, IPOSINRHSCOMP
      INTEGER DUMMY(1)
      INTEGER PLEFTW, PTWCB
      INTEGER Offset, EffectiveSize, ISLAVE, FirstIndex
      LOGICAL LTLEVEL2, IN_SUBTREE
      INTEGER TYPENODE
      INCLUDE 'mumps_headers.h'
      LOGICAL BLOCK_SEQUENCE
      INTEGER TMP_NBPANELS, I_PIVRPTR, I_PIVR
      LOGICAL MUST_BE_PERMUTED
      LOGICAL NO_CHILDREN
      LOGICAL Exploit_Sparsity, AM1             
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      INTEGER(8) :: APOSDEB, NBENTRIES_ALLPANELS
      INTEGER LDAJ, NBJ, LIWFAC,
     &        NBJLAST, NPIV_LAST, PANEL_SIZE,
     &        PTWCB_PANEL, NCB_PANEL, TYPEF
      INTEGER BEG_PANEL
      LOGICAL TWOBYTWO
      INTEGER NPANELS, IPANEL
      LOGICAL MUMPS_170
      INTEGER MUMPS_330
      EXTERNAL sgemv, strsv, strsm, sgemm,
     &         MUMPS_330, 
     &         MUMPS_170
      PLEFTW = 1
      POSIWCB = LIWW
      POSWCB = LWC
      NROOT = 0
      NBLEAF = NA(1)
      NBROOT = NA(2)
      DO I = NBROOT, 1, -1
        INODE = NA(NBLEAF+I+2)
        IF (MUMPS_275(PROCNODE_STEPS(STEP(INODE)),
     &      SLAVEF) .EQ. MYID) THEN
          NROOT = NROOT + 1
          IPOOL(NROOT) = INODE
        ENDIF
      END DO
      III = 1
      IIPOOL = NROOT + 1                 
      BLOCK_SEQUENCE = .FALSE.
      Exploit_Sparsity = .FALSE.
      AM1 = .FALSE.
      IF (KEEP(235).NE.0) Exploit_Sparsity = .TRUE.
      IF (KEEP(237).NE.0) AM1 = .TRUE.
      NO_CHILDREN = .FALSE.
      IF (Exploit_Sparsity .OR. AM1) MYLEAF = -1
      IF (MYLEAF .EQ. -1) THEN
        MYLEAF = 0
        DO I=1, NBLEAF
          INODE=NA(I+2)
          IF (MUMPS_275(PROCNODE_STEPS(STEP(INODE)),
     &         SLAVEF) .EQ. MYID) THEN
            MYLEAF = MYLEAF + 1
          ENDIF
        ENDDO
      ENDIF
      MYLEAFE=MYLEAF
      NBFINF = SLAVEF
      IF (MYLEAFE .EQ. 0) THEN
        CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, FEUILLE,
     &                  SLAVEF)
        NBFINF = NBFINF - 1
        IF (NBFINF .EQ. 0) THEN
          GOTO 340
        ENDIF
      ENDIF
 50   CONTINUE
      BLOQ = ( (  III .EQ. IIPOOL  )
     &     )
      CALL SMUMPS_41( BLOQ, FLAG, BUFR, LBUFR,
     &     LBUFR_BYTES, MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &     STEP,  FRERE, FILS, PROCNODE_STEPS,
     &     PLEFTW, KEEP,KEEP8,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &     RHS, LRHS, NRHS, MTYPE, 
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP
     &     , TO_PROCESS, SIZE_TO_PROCESS
     &     )
      IF ( INFO(1) .LT. 0 ) GOTO 340
      IF ( .NOT. FLAG ) THEN
        IF (III .NE. IIPOOL) THEN
          INODE = IPOOL(IIPOOL-1)
          IIPOOL = IIPOOL - 1
          GO TO 60
        ENDIF
      END IF                    
      IF ( NBFINF .eq. 0 ) GOTO 340
      GOTO 50
   60 CONTINUE
      IF ( INODE .EQ. KEEP( 38 ) .OR. INODE .EQ. KEEP( 20 ) ) THEN
         IPOS = PTRIST(STEP(INODE))+KEEP(IXSZ)
          NPIV  = IW(IPOS+3)
          LIELL = IW(IPOS) + NPIV  
         IPOS =  PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ)
         IF ( MTYPE .EQ. 1 .AND. KEEP(50) .EQ. 0) THEN
            J1   = IPOS + LIELL + 1
            J2   = IPOS + LIELL + NPIV
         ELSE
            J1   = IPOS + 1
            J2   = IPOS + NPIV
         END IF
         IFR  = 0
         DO JJ = J1, J2
            J  = IW( JJ )
            IFR = IFR + 1
            DO K=1,NRHS
               RHS(J,K) = RHS_ROOT(IFR+NPIV*(K-1))
            END DO
         END DO 
         IN = INODE
 270     IN = FILS(IN)
         IF (IN .GT. 0) GOTO 270
         IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
            ENDIF
            GOTO 50
         ENDIF
         IF   = -IN
         LONG = NPIV
         NBFILS = NE_STEPS(STEP(INODE))
         IF ( AM1 ) THEN
            I = NBFILS
            NBFILS = 0
            DO WHILE (I.GT.0)
               IF ( TO_PROCESS(STEP(IF)) ) NBFILS = NBFILS+1
               IF = FRERE(STEP(IF))
               I = I -1
            ENDDO
            IF (NBFILS.EQ.0) THEN
               NO_CHILDREN = .TRUE.
            ELSE
               NO_CHILDREN = .FALSE.
            ENDIF
            IF = -IN
         ENDIF
         DEBUT = .TRUE.
         DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
         END DO
         POOL_FIRST_POS=IIPOOL
         DO I = 1, NBFILS
            IF ( AM1 ) THEN
 1030          IF ( .NOT.TO_PROCESS(STEP(IF)) ) THEN
                  IF = FRERE(STEP(IF))
                  GOTO 1030
               ENDIF
               NO_CHILDREN = .FALSE.
            ENDIF
            IF (MUMPS_275(PROCNODE_STEPS(STEP(IF)),SLAVEF)
     &           .EQ. MYID) THEN
                  IPOOL(IIPOOL) = IF
                  IIPOOL = IIPOOL + 1
            ELSE
               PROCDEST = MUMPS_275(PROCNODE_STEPS(STEP(IF)),
     &              SLAVEF)
               IF (.NOT. DEJA_SEND( PROCDEST ))  THEN
 600              CALL SMUMPS_78( NRHS, IF, 0, 0,
     &                 LONG, LONG, IW( J1 ),
     &                 RHS_ROOT( 1 ), PROCDEST,
     &                 NOEUD, COMM, IERR )
                  IF ( IERR .EQ. -1 ) THEN
                     CALL SMUMPS_41(
     &                    .FALSE., FLAG,
     &                    BUFR, LBUFR, LBUFR_BYTES,
     &                    MYID, SLAVEF, COMM,
     &                    N, IWCB, LIWW, POSIWCB,
     &                    W, LWC, POSWCB,
     &                    IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                    IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &                    STEP, FRERE, FILS, PROCNODE_STEPS,
     &                    PLEFTW, KEEP,KEEP8,
     &                    PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                    RHS, LRHS, NRHS, MTYPE,
     &                    RHSCOMP, LRHSCOMP, POSINRHSCOMP
     &                    , TO_PROCESS, SIZE_TO_PROCESS
     &                    )
                     IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                     GOTO 600
                  ELSE IF ( IERR .EQ. -2 ) THEN
                     INFO( 1 ) = -17
                     INFO( 2 ) = LONG * KEEP(35) +
     &                    ( LONG + 2 ) * KEEP(34)
                     GOTO 330
                  ELSE IF ( IERR .EQ. -3 ) THEN
                     INFO( 1 ) = -20
                     INFO( 2 ) = LONG * KEEP(35) +
     &                    ( LONG + 2 ) * KEEP(34)
                     GOTO 330
                  END IF
                  DEJA_SEND( PROCDEST ) = .TRUE.
               END IF
               IF ( IERR .NE. 0 ) CALL MUMPS_ABORT()
            ENDIF
            IF = FRERE(STEP(IF))
         ENDDO
         IF (AM1 .AND.NO_CHILDREN) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
               GOTO 50
            ENDIF
         ENDIF
            IF (IIPOOL.NE.POOL_FIRST_POS) THEN
               DO I=1,(IIPOOL-POOL_FIRST_POS)/2
                  TMP=IPOOL(POOL_FIRST_POS+I-1)
                  IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
                  IPOOL(IIPOOL-I)=TMP
               ENDDO
            ENDIF
         GOTO 50
      END IF
      IN_SUBTREE = MUMPS_170( 
     &               PROCNODE_STEPS(STEP(INODE)), SLAVEF ) 
      TYPENODE = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),
     &         SLAVEF)
      LTLEVEL2= ( 
     &   (TYPENODE .eq.2 ) .AND.
     &   (MTYPE.NE.1)   )
      NPIV = IW(PTRIST(STEP(INODE))+2+KEEP(IXSZ)+1)
      IF ((NPIV.NE.0).AND.(LTLEVEL2)) THEN
            IPOS  = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
            LIELL = IW(IPOS-2)+IW(IPOS+1)
            NELIM = IW(IPOS-1)
            IPOS  = IPOS + 1
            NPIV  = IW(IPOS)
            NCB   = LIELL - NPIV - NELIM
            IPOS  = IPOS + 2
            NSLAVES = IW( IPOS )
            Offset = 0  
            IPOS = IPOS + NSLAVES   
            IW(PTRIST(STEP(INODE))+XXS)= C_FINI+NSLAVES
           IF ( POSIWCB - 2 .LT. 0 .or.
     &          POSWCB - NCB*NRHS .LT. PLEFTW - 1 ) THEN
             CALL SMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
             IF ( POSWCB - NCB*NRHS .LT. PLEFTW - 1 ) THEN
               INFO( 1 ) = -11
               INFO( 2 ) = NCB * NRHS - POSWCB - PLEFTW + 1
               GOTO 330
             END IF
             IF ( POSIWCB - 2 .LT. 0 ) THEN
               INFO( 1 ) = -14
               INFO( 2 ) = 2 - POSIWCB
               GO TO 330
             END IF
           END IF
           POSIWCB = POSIWCB - 2
           POSWCB  = POSWCB - NCB*NRHS
           PTRICB(STEP( INODE )) = POSIWCB + 1
           PTRACB(STEP( INODE )) = POSWCB  + 1
           IWCB( PTRICB(STEP( INODE ))     ) = NCB
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
           IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
              POSINDICES = IPOS + LIELL + 1
           ELSE
              POSINDICES = IPOS + 1
           END IF
           IF ( NCB.EQ.0 ) THEN
             write(6,*) ' Internal Error type 2 node with no CB '
             CALL MUMPS_ABORT()
           ENDIF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
               J1 = IPOS + LIELL + NPIV + NELIM +1
               J2 = IPOS + 2 * LIELL
           ELSE
               J1 = IPOS + NPIV + NELIM +1
               J2 = IPOS + LIELL
           END IF
           IFR = PTRACB(STEP( INODE )) - 1
           DO JJ = J1, J2 - KEEP(253)
               J = IW(JJ)
               IFR = IFR + 1
               DO K=1, NRHS
                 W(IFR+(K-1)*NCB) = RHS(J,K)
               ENDDO
           ENDDO
           IF (KEEP(252).NE.0) THEN
             DO JJ = J2-KEEP(253)+1, J2
              IFR = IFR + 1
              DO K=1, NRHS
               IF (K.EQ.JJ-J2+KEEP(253)) THEN
                 W(IFR+(K-1)*NCB) = ALPHA   
               ELSE
                 W(IFR+(K-1)*NCB) = ZERO
               ENDIF
              ENDDO
             ENDDO
           ENDIF
           DO ISLAVE = 1, NSLAVES
              CALL MUMPS_49( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB, 
     &                NSLAVES, 
     &                EffectiveSize,
     &                FirstIndex )
 500         DEST = IW( PTRIST(STEP(INODE))+5+ISLAVE+KEEP(IXSZ))
             CALL SMUMPS_63(NRHS, INODE,
     &             W(Offset+PTRACB(STEP(INODE))), EffectiveSize, 
     &             NCB, DEST,
     &             BACKSLV_MASTER2SLAVE,
     &             COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                 CALL SMUMPS_41(
     &                .FALSE., FLAG,
     &                BUFR, LBUFR, LBUFR_BYTES,
     &                MYID, SLAVEF, COMM,
     &                N, IWCB, LIWW, POSIWCB,
     &                W, LWC, POSWCB,
     &                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &                STEP, FRERE, FILS,
     &                PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     &                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                RHS, LRHS, NRHS, MTYPE,
     &                RHSCOMP, LRHSCOMP, POSINRHSCOMP
     &                , TO_PROCESS, SIZE_TO_PROCESS
     &                )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 500
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = EffectiveSize * KEEP(35) +
     &                            2 * KEEP(34)
                GOTO 330
              ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = EffectiveSize * KEEP(35) +
     &                            2 * KEEP(34)
                GOTO 330
              END IF
              Offset = Offset + EffectiveSize
           END DO
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 0
           CALL SMUMPS_151(NRHS,N, KEEP(28), IWCB, LIWW, W, LWC,
     &             POSWCB,POSIWCB,PTRICB,PTRACB)
           GOTO 50
      ENDIF   
      IPOS = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      NELIM = IW(IPOS-1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      IPOS = IPOS + 1
      IF (KEEP(201).GT.0) THEN
         CALL SMUMPS_643(
     &        INODE,PTRFAC,KEEP,A,LA,STEP,
     &        KEEP8,N,MUST_BE_PERMUTED,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF                     
      APOS = PTRFAC(IW(IPOS))
      NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ) )
      IPOS = IPOS + 1 + NSLAVES
      IF (KEEP(201).EQ.1) THEN 
           LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
           IF (MTYPE.NE.1) THEN
            TYPEF = TYPEF_L
           ELSE
            TYPEF = TYPEF_U
           ENDIF
           PANEL_SIZE =  SMUMPS_690( LIELL )
           IF (KEEP(50).NE.1) THEN
             CALL SMUMPS_755(
     &                   IW(IPOS+1+2*LIELL),
     &                   MUST_BE_PERMUTED )
           ENDIF
      ENDIF  
      LONG = 0
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ELSE
        J1 = IPOS + LIELL + 1
        J2 = IPOS + NPIV + LIELL
      END IF
      IF (IN_SUBTREE) THEN
        PTWCB = PLEFTW
        IF ( POSWCB .LT. LIELL*NRHS ) THEN
          CALL SMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     &                 POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB .LT. LIELL*NRHS ) THEN
            INFO(1) = -11
            INFO(2) = LIELL*NRHS - POSWCB
            GOTO 330
          END IF
        END IF
      ELSE
        IF ( POSIWCB - 2 .LT. 0 .or.
     &     POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
          CALL SMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
            INFO( 1 ) = -11
            INFO( 2 ) = LIELL * NRHS - POSWCB - PLEFTW + 1
            GOTO 330
          END IF
          IF ( POSIWCB - 2 .LT. 0 ) THEN
            INFO( 1 ) = -14
            INFO( 2 ) = 2 - POSIWCB
            GO TO 330
          END IF
        END IF
        POSIWCB = POSIWCB - 2
        POSWCB  = POSWCB - LIELL*NRHS
        PTRICB(STEP( INODE )) = POSIWCB + 1
        PTRACB(STEP( INODE )) = POSWCB  + 1
        IWCB( PTRICB(STEP( INODE ))     ) = LIELL
        IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
        IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
           POSINDICES = IPOS + LIELL + 1
        ELSE
           POSINDICES = IPOS + 1
        END IF
        PTWCB = PTRACB(STEP( INODE )) 
      ENDIF
      IF (KEEP(252).EQ.0) IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
      DO K=1, NRHS
        IF (KEEP(252).NE.0) THEN
         DO JJ = J1, J2
          W(PTWCB+JJ-J1+(K-1)*LIELL) = ZERO
         ENDDO
        ELSE
         DO JJ = J1, J2
          W(PTWCB+JJ-J1+(K-1)*LIELL) = RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
         ENDDO
        ENDIF
      END DO
      IFR   = PTWCB + NPIV - 1
      IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
          J1 = IPOS + LIELL + NPIV + 1
          J2 = IPOS + 2 * LIELL
        ELSE
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
        END IF
        DO JJ = J1, J2-KEEP(253)
          J = IW(JJ)
          IFR = IFR + 1
          DO K=1, NRHS
            W(IFR+(K-1)*LIELL) = RHS(J,K)
          ENDDO
        ENDDO
        IF (KEEP(252).NE.0) THEN
          DO JJ = J2-KEEP(253)+1, J2
           IFR = IFR + 1
           DO K=1, NRHS
            IF (K.EQ.JJ-J2+KEEP(253)) THEN
              W(IFR+(K-1)*LIELL) = ALPHA   
            ELSE
              W(IFR+(K-1)*LIELL) = ZERO
            ENDIF
           ENDDO
          ENDDO
        ENDIF
        NCB = LIELL - NPIV
        IF (NPIV .EQ. 0) GOTO 160
      ENDIF
      IF (KEEP(201).EQ.1) THEN 
       J = NPIV / PANEL_SIZE 
       TWOBYTWO = KEEP(50).EQ.2 .AND.
     & ((TYPENODE.EQ.1.AND.KEEP(103).GT.0) .OR.
     &  (TYPENODE.EQ.2.AND.KEEP(105).GT.0))
       IF (TWOBYTWO) THEN 
         CALL SMUMPS_641(PANEL_SIZE, PANEL_POS, LPANEL_POS,
     &        IW(IPOS+1+LIELL), NPIV, NPANELS, LIELL,
     &        NBENTRIES_ALLPANELS)
       ELSE
         IF (NPIV.EQ.J*PANEL_SIZE) THEN
           NPIV_LAST = NPIV
           NBJLAST   = PANEL_SIZE
           NPANELS   = J
         ELSE
           NPIV_LAST = (J+1)* PANEL_SIZE
           NBJLAST   = NPIV-J*PANEL_SIZE
           NPANELS   = J+1
         ENDIF
            NBENTRIES_ALLPANELS =
     &  int(LIELL,8) * int(NPIV,8) 
     &  - int( ( J * ( J - 1 ) ) / 2,8 ) 
     &    * int(PANEL_SIZE,8) * int(PANEL_SIZE,8) 
     &  - int(J,8)                       
     &    * int(mod(NPIV, PANEL_SIZE),8) 
     &    * int(PANEL_SIZE,8)    
         JJ=NPIV_LAST
       ENDIF
       APOSDEB = APOS + NBENTRIES_ALLPANELS 
       DO IPANEL = NPANELS, 1, -1
            IF (TWOBYTWO) THEN
              NBJ = PANEL_POS(IPANEL+1)-PANEL_POS(IPANEL)
              BEG_PANEL = PANEL_POS(IPANEL)
            ELSE
              IF (JJ.EQ.NPIV_LAST) THEN
                NBJ = NBJLAST
              ELSE
                NBJ = PANEL_SIZE
              ENDIF
              BEG_PANEL = JJ- PANEL_SIZE+1
            ENDIF
            LDAJ    = LIELL-BEG_PANEL+1 
            APOSDEB = APOSDEB - int(NBJ,8)*int(LDAJ,8)
            PTWCB_PANEL = PTWCB + BEG_PANEL - 1
            NCB_PANEL   = LDAJ - NBJ
            IF (KEEP(50).NE.1 .AND. MUST_BE_PERMUTED) THEN
              CALL SMUMPS_667(TYPEF, TMP_NBPANELS,
     &        I_PIVRPTR, I_PIVR, IPOS + 1 + 2 * LIELL, IW, LIW)
              IF (NPIV.EQ.(IW(I_PIVRPTR)-1)) THEN
                MUST_BE_PERMUTED=.FALSE. 
              ELSE
               CALL SMUMPS_698(
     &         IW(I_PIVR + IW(I_PIVRPTR+IPANEL-1)-IW(I_PIVRPTR)),
     &         NPIV-IW(I_PIVRPTR+IPANEL-1)+1,
     &         IW(I_PIVRPTR+IPANEL-1)-1,
     &         A(APOSDEB),
     &         LDAJ, NBJ, BEG_PANEL-1)
              ENDIF
            ENDIF
            IF ( NRHS == 1 ) THEN
              IF (NCB_PANEL.NE.0) THEN
                CALL sgemv( 'T', NCB_PANEL, NBJ, ALPHA, 
     &                A( APOSDEB + int(NBJ,8) ), LDAJ,
     &                W( NBJ + PTWCB_PANEL ),
     &                1, ONE,
     &                W(PTWCB_PANEL), 1 )
              ENDIF
              IF (MTYPE.NE.1) THEN
               CALL strsv('L','T','U', NBJ, A(APOSDEB), LDAJ,
     &              W(PTWCB_PANEL), 1)
              ELSE
               CALL strsv('L','T','N', NBJ, A(APOSDEB), LDAJ,
     &              W(PTWCB_PANEL), 1)
              ENDIF
            ELSE
              IF (NCB_PANEL.NE.0) THEN
                 CALL sgemm( 'T', 'N', NBJ, NRHS, NCB_PANEL, ALPHA,
     &              A(APOSDEB +int(NBJ,8)), LDAJ,
     &              W(NBJ+PTWCB_PANEL),LIELL,
     &              ONE, W(PTWCB_PANEL),LIELL)
              ENDIF
              IF (MTYPE.NE.1) THEN
               CALL strsm('L','L','T','U',NBJ, NRHS, ONE, 
     &           A(APOSDEB), 
     &           LDAJ, W(PTWCB_PANEL), LIELL)
              ELSE
               CALL strsm('L','L','T','N',NBJ, NRHS, ONE, 
     &           A(APOSDEB), 
     &           LDAJ, W(PTWCB_PANEL), LIELL)
              ENDIF
            ENDIF
            IF (.NOT. TWOBYTWO) JJ=BEG_PANEL-1 
       ENDDO 
      ENDIF 
      IF (KEEP(201).EQ.0.OR.KEEP(201).EQ.2)THEN 
       IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          IST = APOS + int(NPIV,8)
          IF (NRHS == 1) THEN
            CALL sgemv( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     &              W(NPIV + PTWCB), 1,
     &              ONE,
     &              W(PTWCB), 1 )
          ELSE
            CALL sgemm('T','N', NPIV, NRHS, NCB, ALPHA, A(IST), LIELL,
     &              W(NPIV+PTWCB), LIELL, ONE,
     &              W(PTWCB), LIELL)
          ENDIF
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
            IST = APOS + int(NPIV,8) * int(LIELL,8)
          ELSE
            IST = APOS + int(NPIV,8) * int(NPIV,8)
          END IF
            IF ( NRHS == 1 ) THEN
              CALL sgemv( 'N', NPIV, NCB, ALPHA, A( IST ), NPIV,
     &                W( NPIV + PTWCB ),
     &                1, ONE,
     &                W(PTWCB), 1 )
            ELSE
                CALL sgemm( 'N', 'N', NPIV, NRHS, NCB, ALPHA,
     &                A(IST), NPIV, W(NPIV+PTWCB),LIELL,
     &                ONE, W(PTWCB),LIELL)
            END IF
        END IF 
       ENDIF  
       IF ( MTYPE .eq. 1 ) THEN
        IF ( NRHS == 1 ) THEN
          CALL strsv('L', 'T', 'N', NPIV, A(APOS), LIELL,
     &              W(PTWCB), 1)
        ELSE
          CALL strsm('L','L','T','N', NPIV, NRHS, ONE, A(APOS),
     &              LIELL, W(PTWCB), LIELL)
        ENDIF
       ELSE
        IF ( KEEP(50) .EQ. 0 ) THEN
          IF ( NRHS == 1 ) THEN
            CALL strsv('U','N','U', NPIV, A(APOS), LIELL,
     &              W(PTWCB), 1)
          ELSE
            CALL strsm('L','U','N','U', NPIV, NRHS, ONE, A(APOS),
     &                 LIELL,W(PTWCB),LIELL)
          END IF
        ELSE
          IF ( NRHS == 1 ) THEN
            CALL strsv('U','N','U', NPIV, A(APOS), NPIV,
     &              W(PTWCB), 1)
          ELSE
            CALL strsm('L','U','N','U',NPIV, NRHS, ONE, A(APOS),
     &           NPIV, W(PTWCB), LIELL)
          END IF
        END IF
       END IF 
      ENDIF 
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      DO 150 I = 1, NPIV
        JJ = IW(J1 + I - 1)
        DO K=1, NRHS
          RHS(JJ, K) = W(PTWCB+I-1+(K-1)*LIELL)
        ENDDO
  150 CONTINUE
  160 CONTINUE
      IF (KEEP(201).GT.0) THEN
         CALL SMUMPS_598(INODE,PTRFAC,KEEP(28),
     &        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF
      IN = INODE
  170 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 170
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (MYLEAFE .EQ. 0) THEN
          CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &                     FEUILLE, SLAVEF )
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        GOTO 50
      ENDIF
      IF = -IN
      NBFILS = NE_STEPS(STEP(INODE))
      IF (AM1) THEN
         I = NBFILS
         NBFILS = 0
         DO WHILE (I.GT.0)
            IF ( TO_PROCESS(STEP(IF)) ) NBFILS = NBFILS+1
            IF = FRERE(STEP(IF))
            I = I -1
         ENDDO
         IF (NBFILS.EQ.0) THEN
            NO_CHILDREN = .TRUE.
         ELSE
            NO_CHILDREN = .FALSE.
         ENDIF
         IF = -IN
      ENDIF
      IF (IN_SUBTREE) THEN
         DO I = 1, NBFILS
            IF ( AM1 ) THEN
 1010          IF ( .NOT.TO_PROCESS(STEP(IF)) )  THEN
                  IF = FRERE(STEP(IF))
                  GOTO 1010
               ENDIF
               NO_CHILDREN = .FALSE.
            ENDIF
               IPOOL((IIPOOL-I+1)+NBFILS-I) = IF
               IIPOOL = IIPOOL + 1
            IF = FRERE(STEP(IF))
         ENDDO
         IF (AM1 .AND. NO_CHILDREN) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
               GOTO 50
            ENDIF
         ENDIF
      ELSE
        DEBUT = .TRUE.
        DO I = 0, SLAVEF - 1
          DEJA_SEND( I ) = .FALSE.
        END DO
        POOL_FIRST_POS=IIPOOL
        DO 190 I = 1, NBFILS
           IF ( AM1 ) THEN
1020      IF ( .NOT.TO_PROCESS(STEP(IF)) ) THEN
                 IF = FRERE(STEP(IF))
                 GOTO 1020
              ENDIF
              NO_CHILDREN = .FALSE.
           ENDIF
          IF (MUMPS_275(PROCNODE_STEPS(STEP(IF)),
     &      SLAVEF) .EQ. MYID) THEN
                IPOOL(IIPOOL) = IF
                IIPOOL = IIPOOL + 1
            IF = FRERE(STEP(IF))
          ELSE
            PROCDEST = MUMPS_275(PROCNODE_STEPS(STEP(IF)),SLAVEF)
            IF (.not. DEJA_SEND( PROCDEST ))  THEN
 400          CONTINUE
              CALL SMUMPS_78( NRHS, IF, 0, 0, LIELL,
     &         LIELL - KEEP(253),
     &         IW( POSINDICES ), 
     &         W   ( PTRACB(STEP( INODE ))), PROCDEST,
     &         NOEUD, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL SMUMPS_41(
     &          .FALSE., FLAG,
     &          BUFR, LBUFR, LBUFR_BYTES,
     &          MYID, SLAVEF, COMM,
     &          N, IWCB, LIWW, POSIWCB,
     &          W, LWC, POSWCB,
     &          IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &          IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &          STEP, FRERE, FILS, PROCNODE_STEPS,
     &          PLEFTW, KEEP,KEEP8,
     &          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &          RHS, LRHS, NRHS, MTYPE, 
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP
     &                , TO_PROCESS, SIZE_TO_PROCESS
     &                )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 400
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
                GOTO 330
              ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
                GOTO 330
              END IF
              DEJA_SEND( PROCDEST ) = .TRUE.
            END IF
            IF = FRERE(STEP(IF))
          ENDIF
  190   CONTINUE
        IF (AM1 .AND. NO_CHILDREN) THEN
           MYLEAFE = MYLEAFE - 1
           IF (MYLEAFE .EQ. 0) THEN
              CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &             FEUILLE, SLAVEF )
              NBFINF = NBFINF - 1
              IF (NBFINF .EQ. 0) GOTO 340
              GOTO 50
           ENDIF
        ENDIF
           DO I=1,(IIPOOL-POOL_FIRST_POS)/2
              TMP=IPOOL(POOL_FIRST_POS+I-1)
              IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
              IPOOL(IIPOOL-I)=TMP
           ENDDO 
        IWCB(PTRICB(STEP(INODE))+1) = IWCB(PTRICB(STEP(INODE))+1)-1
        CALL SMUMPS_151(NRHS,N, KEEP(28), IWCB, LIWW, 
     &     W, LWC,
     &     POSWCB,POSIWCB,PTRICB,PTRACB)
      ENDIF
      GOTO 50
  330 CONTINUE
      CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, TERREUR,
     & SLAVEF)
  340 CONTINUE
      CALL SMUMPS_150( MYID,COMM,BUFR,
     &                            LBUFR,LBUFR_BYTES )
      RETURN
      END SUBROUTINE SMUMPS_249
      RECURSIVE SUBROUTINE SMUMPS_41(
     &     BLOQ, FLAG,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &     STEP, FRERE, FILS, PROCNODE_STEPS,
     &     PLEFTW, KEEP,KEEP8,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, RHS,
     &     LRHS, NRHS, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP
     &     , TO_PROCESS, SIZE_TO_PROCESS
     &     )
      IMPLICIT NONE
      LOGICAL BLOQ, FLAG
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      REAL W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, LPOOL
      INTEGER IPOOL( LPOOL )
      INTEGER LPANEL_POS
      INTEGER PANEL_POS( LPANEL_POS )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER PTRIST(KEEP(28)), IW( LIW )
      INTEGER (8) :: PTRFAC(KEEP(28))
      REAL A( LA ), W2( KEEP(133) )
      INTEGER LRHS, NRHS
      REAL RHS(LRHS, NRHS)
      INTEGER MYLEAFE, MTYPE
      INTEGER LRHSCOMP, POSINRHSCOMP(KEEP(28))
      REAL RHSCOMP(LRHSCOMP,NRHS)
      INTEGER SIZE_TO_PROCESS
      LOGICAL TO_PROCESS(SIZE_TO_PROCESS)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGSOU, MSGTAG, MSGLEN
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &                   FLAG, STATUS, IERR )
      END IF
      IF (FLAG) THEN
         MSGSOU=STATUS(MPI_SOURCE)
         MSGTAG=STATUS(MPI_TAG)
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL SMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV(BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     &                   MSGTAG, COMM, STATUS, IERR)
           CALL SMUMPS_42( MSGTAG, MSGSOU,
     &                BUFR, LBUFR, LBUFR_BYTES,
     &                MYID, SLAVEF, COMM,
     &                N, IWCB, LIWW, POSIWCB,
     &                W, LWC, POSWCB,
     &                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &                FRERE, FILS, PROCNODE_STEPS, PLEFTW,
     &                KEEP,KEEP8,
     &                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                RHS, LRHS, NRHS, MTYPE, 
     &                RHSCOMP, LRHSCOMP, POSINRHSCOMP 
     &               , TO_PROCESS, SIZE_TO_PROCESS
     &          )
         END IF
      END IF
      RETURN
      END SUBROUTINE SMUMPS_41
      RECURSIVE SUBROUTINE SMUMPS_42(
     &                MSGTAG, MSGSOU,
     &                BUFR, LBUFR, LBUFR_BYTES,
     &                MYID, SLAVEF, COMM,
     &                N, IWCB, LIWW, POSIWCB,
     &                W, LWC, POSWCB,
     &                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &                FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     &                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                RHS, LRHS, NRHS, MTYPE, 
     &                RHSCOMP, LRHSCOMP, POSINRHSCOMP 
     &               , TO_PROCESS, SIZE_TO_PROCESS
     &           )
      USE SMUMPS_OOC
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER MSGTAG, MSGSOU
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      REAL W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, LPOOL, LPANEL_POS
      INTEGER IPOOL( LPOOL )
      INTEGER PANEL_POS( LPANEL_POS )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER FRERE(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER IW( LIW ), PTRIST( KEEP(28) )
      INTEGER(8) :: PTRFAC(KEEP(28))
      REAL A( LA ), W2( KEEP(133) )
      INTEGER LRHS, NRHS
      REAL  RHS(LRHS, NRHS)
      INTEGER MYLEAFE, MTYPE
      INTEGER LRHSCOMP, POSINRHSCOMP(KEEP(28))
      REAL RHSCOMP(LRHSCOMP,NRHS)
      INTEGER TMP_NBPANELS, I_PIVRPTR, I_PIVR
      LOGICAL MUST_BE_PERMUTED
      INTEGER  SIZE_TO_PROCESS
      LOGICAL TO_PROCESS(SIZE_TO_PROCESS), NO_CHILDREN
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER POSITION, IF, INODE, IERR, LONG, DUMMY(1)
      INTEGER P_UPDATE, P_SOL_MAS, LIELL, K
      INTEGER(8) :: APOS, IST
      INTEGER NPIV, NROW_L, IPOS, NROW_RECU
      INTEGER I, JJ, IN, PROCDEST, J1, J2, IFR, LDA
      INTEGER NSLAVES, NELIM, J, POSINDICES, INODEPOS,
     &        IPOSINRHSCOMP
      LOGICAL FLAG
      REAL ZERO, ALPHA, ONE
      PARAMETER (ZERO=0.0E0, ONE = 1.0E0, ALPHA=-1.0E0)
      INCLUDE 'mumps_headers.h'
      INTEGER POOL_FIRST_POS, TMP
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275, strsv, strsm, sgemv, sgemm
      INTEGER(8) :: APOSDEB, NBENTRIES_ALLPANELS
      INTEGER LDAJ, NBJ, LIWFAC,
     &        NBJLAST, NPIV_LAST, PANEL_SIZE,
     &        PTWCB_PANEL, NCB_PANEL, TYPEF
      LOGICAL TWOBYTWO
      INTEGER BEG_PANEL
      INTEGER IPANEL, NPANELS
      IF (MSGTAG .EQ. FEUILLE) THEN
          NBFINF = NBFINF - 1
      ELSE IF (MSGTAG .EQ. NOEUD) THEN
          POSITION = 0
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        INODE, 1, MPI_INTEGER,
     &        COMM, IERR)
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        LONG, 1, MPI_INTEGER,
     &        COMM, IERR)
          IF (   POSIWCB - LONG - 2 .LT. 0
     &      .OR. POSWCB - PLEFTW + 1 .LT. LONG ) THEN
            CALL SMUMPS_95(NRHS, N, KEEP(28), IWCB,
     &      LIWW, W, LWC,
     &      POSWCB, POSIWCB, PTRICB, PTRACB)
            IF ((POSIWCB - LONG - 2 ) .LT. 0) THEN
              INFO(1)=-14
              INFO(2)=-POSIWCB + LONG + 2
              WRITE(6,*) MYID,' Internal error in bwd solve COMPSO'
              GOTO 260
            END IF
            IF ( POSWCB - PLEFTW + 1 .LT. LONG ) THEN
              INFO(1) = -11
              INFO(2) = LONG + PLEFTW - POSWCB - 1
              WRITE(6,*) MYID,' Internal error in bwd solve COMPSO'
              GOTO 260
            END IF
          ENDIF
          POSIWCB = POSIWCB - LONG
          POSWCB = POSWCB - LONG
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &          IWCB(POSIWCB + 1), 
     &          LONG, MPI_INTEGER, COMM, IERR)
            DO K=1,NRHS
             CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &          W(POSWCB + 1), LONG, 
     &          MPI_REAL, COMM, IERR)
             DO JJ=0, LONG-1
               RHS(IWCB(POSIWCB+1+JJ),K) = W(POSWCB+1+JJ)
             ENDDO
            ENDDO
            POSIWCB = POSIWCB + LONG
            POSWCB = POSWCB + LONG
          ENDIF
          POOL_FIRST_POS = IIPOOL
          IF ( KEEP(237).GT. 0 ) THEN
             IF (.NOT.TO_PROCESS(STEP(INODE))) 
     &            GOTO 1010
          ENDIF
             IPOOL( IIPOOL ) = INODE
             IIPOOL = IIPOOL + 1
 1010     CONTINUE
          IF = FRERE( STEP(INODE) )
          DO WHILE ( IF .GT. 0 )
             IF ( MUMPS_275(PROCNODE_STEPS(STEP(IF)),
     &            SLAVEF) .eq. MYID ) THEN
                IF ( KEEP(237).GT. 0 ) THEN
                   IF (.NOT.TO_PROCESS(STEP(IF))) THEN
                      IF = FRERE(STEP(IF))
                      CYCLE
                   ENDIF
                ENDIF
                   IPOOL( IIPOOL ) = IF
                   IIPOOL = IIPOOL + 1
             END IF
             IF = FRERE( STEP( IF ) )
          END DO
             DO I=1,(IIPOOL-POOL_FIRST_POS)/2
                TMP=IPOOL(POOL_FIRST_POS+I-1)
                IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
                IPOOL(IIPOOL-I)=TMP
             ENDDO      
      ELSE IF ( MSGTAG .EQ. BACKSLV_MASTER2SLAVE ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   INODE, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NROW_RECU, 1, MPI_INTEGER, COMM, IERR )
        IPOS   = PTRIST( STEP(INODE) ) + KEEP(IXSZ)
        NPIV   = - IW( IPOS     )
        NROW_L =   IW( IPOS + 1 )
        IF (KEEP(201).GT.0) THEN
           CALL SMUMPS_643(
     &     INODE,PTRFAC,KEEP,A,LA,STEP,
     &     KEEP8,N,MUST_BE_PERMUTED,IERR)           
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
        ENDIF                     
        APOS   =   PTRFAC(IW( IPOS + 3 ))
        IF ( NROW_L .NE. NROW_RECU ) THEN
          WRITE(*,*) 'Error1 in 41S : NROW L/RECU=',NROW_L, NROW_RECU
          CALL MUMPS_ABORT()
        END IF
        LONG = NROW_L + NPIV
        IF ( POSWCB - LONG*NRHS .LT. PLEFTW - 1 ) THEN
           CALL SMUMPS_95(NRHS, N, KEEP(28), IWCB,
     &          LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
           IF ( POSWCB - LONG*NRHS .LT. PLEFTW - 1 ) THEN
             INFO(1) = -11
             INFO(2) = LONG * NRHS- POSWCB
             WRITE(6,*) MYID,' Internal error in bwd solve COMPSO'
             GOTO 260
           END IF
        END IF
        P_UPDATE  = PLEFTW
        P_SOL_MAS = PLEFTW + NPIV * NRHS
        PLEFTW    = P_SOL_MAS + NROW_L * NRHS
        DO K=1, NRHS
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   W( P_SOL_MAS+(K-1)*NROW_L),NROW_L,
     &                   MPI_REAL,
     &                   COMM, IERR )
        ENDDO
        IF (KEEP(201).EQ.1) THEN 
          IF ( NRHS == 1 ) THEN
           CALL sgemv( 'T', NROW_L, NPIV, ALPHA, A( APOS ), NROW_L,
     &              W( P_SOL_MAS ), 1, ZERO,
     &              W( P_UPDATE ), 1 )
          ELSE
           CALL sgemm( 'T', 'N', NPIV, NRHS, NROW_L, ALPHA, A(APOS),
     &           NROW_L, W( P_SOL_MAS ), NROW_L, ZERO, W( P_UPDATE ),
     &           NPIV )
          ENDIF
        ELSE
          IF ( NRHS == 1 ) THEN
           CALL sgemv( 'N', NPIV, NROW_L, ALPHA, A( APOS ), NPIV,
     &              W( P_SOL_MAS ), 1, ZERO,
     &              W( P_UPDATE ), 1 )
          ELSE
           CALL sgemm( 'N', 'N', NPIV, NRHS, NROW_L, ALPHA, A(APOS),
     &            NPIV, W( P_SOL_MAS ), NROW_L, ZERO, W( P_UPDATE ),
     &            NPIV )
          END IF
        ENDIF 
        IF (KEEP(201).GT.0) THEN
         CALL SMUMPS_598(INODE,PTRFAC,KEEP(28),
     &          A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 260
         ENDIF
        ENDIF
        PLEFTW = PLEFTW - NROW_L * NRHS
 100    CONTINUE
        CALL SMUMPS_63( NRHS, INODE, W(P_UPDATE),
     &                               NPIV, NPIV,
     &                                MSGSOU, 
     &                                BACKSLV_UPDATERHS,
     &                                COMM, IERR )
        IF ( IERR .EQ. -1 ) THEN
          CALL SMUMPS_41(
     &     .FALSE., FLAG,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &     FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE,
     &     RHS, LRHS, NRHS, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP
     &      , TO_PROCESS, SIZE_TO_PROCESS
     &          )
          IF ( INFO( 1 ) .LT. 0 ) GOTO 270
          GOTO 100
        ELSE IF ( IERR .EQ. -2 ) THEN
          INFO( 1 ) = -17
          INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
          GOTO 260
        ELSE IF ( IERR .EQ. -3 ) THEN
          INFO( 1 ) = -20
          INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
          GOTO 260
        END IF
        PLEFTW = PLEFTW - NPIV * NRHS
      ELSE IF ( MSGTAG .EQ. BACKSLV_UPDATERHS ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   INODE, 1, MPI_INTEGER, COMM, IERR )
        IPOS  = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR )
          NELIM = IW(IPOS-1)
          IPOS = IPOS + 1
          NPIV = IW(IPOS)
          IPOS = IPOS + 1
          NSLAVES = IW( IPOS + 1 )
          IPOS = IPOS + 1 + NSLAVES
          INODEPOS = PTRIST(STEP(INODE)) + KEEP(IXSZ) + 4
          IF ( KEEP(50) .eq. 0 ) THEN
           LDA = LIELL
          ELSE
           LDA = NPIV
          ENDIF
          IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
             J1 = IPOS + 1
             J2 = IPOS + NPIV
          ELSE
             J1 = IPOS + LIELL + 1
             J2 = IPOS + NPIV + LIELL
          END IF
        DO K=1, NRHS
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   W2, NPIV, MPI_REAL,
     &                   COMM, IERR )
         IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
         I = 1
         IF ( (KEEP(253).NE.0) .AND.
     &        (IW(PTRIST(STEP(INODE))+XXS).EQ.C_FINI+NSLAVES) 
     &      )  THEN
          DO JJ = J1,J2   
            RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = W2(I)
            I = I+1
          ENDDO
         ELSE
          DO JJ = J1,J2   
            RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = 
     &      RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) + W2(I)
            I = I+1
          ENDDO
         ENDIF
        ENDDO  
        IW(PTRIST(STEP(INODE))+XXS) = 
     &      IW(PTRIST(STEP(INODE))+XXS) - 1
        IF ( IW(PTRIST(STEP(INODE))+XXS).EQ.C_FINI ) THEN
          IF (KEEP(201).GT.0) THEN
             CALL SMUMPS_643(
     &            INODE,PTRFAC,KEEP,A,LA,STEP,
     &            KEEP8,N,MUST_BE_PERMUTED,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
               CALL SMUMPS_755(
     &              IW(IPOS+1+2*LIELL),
     &              MUST_BE_PERMUTED )
             ENDIF
          ENDIF  
          APOS = PTRFAC(IW(INODEPOS))
          IF (KEEP(201).EQ.1) THEN 
             LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
             TYPEF = TYPEF_L
             NROW_L   = NPIV+NELIM  
             PANEL_SIZE = SMUMPS_690(NROW_L)
             IF (PANEL_SIZE.LT.0) THEN
               WRITE(6,*) ' Internal error in bwd solve PANEL_SIZE=',
     &         PANEL_SIZE
               CALL MUMPS_ABORT()
             ENDIF
          ENDIF 
           IF ( POSIWCB - 2 .LT. 0 .or.
     &         POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
            CALL SMUMPS_95( NRHS, N, KEEP(28), IWCB, 
     &          LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
            IF ( POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
              INFO( 1 ) = -11
              INFO( 2 ) = LIELL * NRHS - POSWCB - PLEFTW + 1
              GOTO 260
            END IF
            IF ( POSIWCB - 2 .LT. 0 ) THEN
              INFO( 1 ) = -14
              INFO( 2 ) = 2 - POSIWCB
              GO TO 260
            END IF
           END IF
           POSIWCB = POSIWCB - 2
           POSWCB  = POSWCB - LIELL*NRHS
           PTRICB(STEP( INODE )) = POSIWCB + 1
           PTRACB(STEP( INODE )) = POSWCB  + 1
           IWCB( PTRICB(STEP( INODE ))     ) = LIELL
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
           IPOS = PTRIST(STEP(INODE)) + KEEP(IXSZ) + 5 + NSLAVES
           IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
             POSINDICES = IPOS + LIELL + 1
           ELSE
             POSINDICES = IPOS + 1
           END IF
           IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
           IFR = PTRACB(STEP( INODE ))
           DO K=1, NRHS
             DO JJ = J1, J2
               W(IFR+JJ-J1+(K-1)*LIELL) = 
     &           RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
             ENDDO
           END DO
           IFR = PTRACB(STEP(INODE))-1+NPIV
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
             J1 = IPOS + LIELL + NPIV + 1
             J2 = IPOS + 2 * LIELL
           ELSE
             J1 = IPOS + NPIV + 1
             J2 = IPOS + LIELL
           END IF
           DO JJ = J1, J2-KEEP(253)   
              J = IW(JJ)
              IFR = IFR + 1
              DO K=1, NRHS
                W(IFR+(K-1)*LIELL) = RHS(J,K)
              ENDDO
           ENDDO
       IF ( KEEP(201).EQ.1 .AND.
     &    (( NELIM .GT. 0 ).OR. (MTYPE.NE.1 )))  THEN
          J = NPIV / PANEL_SIZE  
          TWOBYTWO = KEEP(50).EQ.2 .AND. KEEP(105).GT.0
          IF (TWOBYTWO) THEN
            CALL SMUMPS_641(PANEL_SIZE, PANEL_POS,
     &           LPANEL_POS, IW(IPOS+1+LIELL), NPIV, NPANELS,
     &           NROW_L, NBENTRIES_ALLPANELS)
          ELSE
            IF (NPIV.EQ.J*PANEL_SIZE) THEN
              NPIV_LAST = NPIV
              NBJLAST   = PANEL_SIZE
              NPANELS   = J
            ELSE
              NPIV_LAST = (J+1)* PANEL_SIZE
              NBJLAST   = NPIV-J*PANEL_SIZE
              NPANELS   = J+1
            ENDIF
            NBENTRIES_ALLPANELS =
     &  int(NROW_L,8) * int(NPIV,8) 
     &  - int( ( J * ( J - 1 ) ) / 2,8 ) 
     &    * int(PANEL_SIZE,8) * int(PANEL_SIZE,8) 
     &  - int(J,8)                       
     &    * int(mod(NPIV, PANEL_SIZE),8) 
     &    * int(PANEL_SIZE,8)    
            JJ=NPIV_LAST
          ENDIF
          APOSDEB = APOS + NBENTRIES_ALLPANELS 
          DO IPANEL=NPANELS,1,-1
            IF (TWOBYTWO) THEN
              NBJ = PANEL_POS(IPANEL+1)-PANEL_POS(IPANEL)
              BEG_PANEL = PANEL_POS(IPANEL)
            ELSE
              IF (JJ.EQ.NPIV_LAST) THEN
                NBJ = NBJLAST
              ELSE
                NBJ = PANEL_SIZE
              ENDIF
              BEG_PANEL = JJ- PANEL_SIZE+1
            ENDIF
            LDAJ    = NROW_L-BEG_PANEL+1 
            APOSDEB = APOSDEB - int(NBJ,8)*int(LDAJ,8)
            PTWCB_PANEL =  PTRACB(STEP(INODE)) + BEG_PANEL - 1
            NCB_PANEL   = LDAJ - NBJ
            IF (KEEP(50).NE.1 .AND.MUST_BE_PERMUTED) THEN
              CALL SMUMPS_667(TYPEF, TMP_NBPANELS,
     &        I_PIVRPTR, I_PIVR, IPOS + 1 + 2 * LIELL, IW, LIW)
              CALL SMUMPS_698(
     &        IW(I_PIVR + IW(I_PIVRPTR+IPANEL-1)-IW(I_PIVRPTR)),
     &        NPIV-IW(I_PIVRPTR+IPANEL-1)+1,
     &        IW(I_PIVRPTR+IPANEL-1)-1,
     &        A(APOSDEB),
     &        LDAJ, NBJ, BEG_PANEL-1)
            ENDIF
            IF ( NRHS == 1 ) THEN
              IF (NCB_PANEL.NE.0) THEN
                CALL sgemv( 'T', NCB_PANEL, NBJ, ALPHA, 
     &                A( APOSDEB + int(NBJ,8) ), LDAJ,
     &                W( NBJ + PTWCB_PANEL ),
     &                1, ONE,
     &                W(PTWCB_PANEL), 1 )
              ENDIF
              CALL strsv('L','T','U', NBJ, A(APOSDEB), LDAJ,
     &              W(PTWCB_PANEL), 1)
            ELSE
              IF (NCB_PANEL.NE.0) THEN
                CALL sgemm( 'T', 'N', NBJ, NRHS, NCB_PANEL, ALPHA,
     &              A(APOSDEB + int(NBJ,8)), LDAJ,
     &              W(NBJ+PTWCB_PANEL),LIELL,
     &              ONE, W(PTWCB_PANEL),LIELL)
              ENDIF
              CALL strsm('L','L','T','U',NBJ, NRHS, ONE, 
     &           A(APOSDEB), 
     &           LDAJ, W(PTWCB_PANEL), LIELL)
            ENDIF
            IF (.NOT. TWOBYTWO) JJ=BEG_PANEL-1
          ENDDO 
        GOTO 1234  
       ENDIF 
          IF (NELIM .GT.0) THEN
            IF ( KEEP(50) .eq. 0 ) THEN
                IST = APOS + int(NPIV,8) * int(LIELL,8)
            ELSE
                IST = APOS + int(NPIV,8) * int(NPIV,8)
            END IF
            IF ( NRHS == 1 ) THEN
                CALL sgemv( 'N', NPIV, NELIM, ALPHA,
     &                A( IST ), NPIV,
     &                W( NPIV + PTRACB(STEP(INODE)) ),
     &                1, ONE,
     &                W(PTRACB(STEP(INODE))), 1 )
             ELSE
                CALL sgemm( 'N', 'N', NPIV, NRHS, NELIM, ALPHA,
     &                A(IST), NPIV, W(NPIV+PTRACB(STEP(INODE))),LIELL,
     &                ONE, W(PTRACB(STEP(INODE))),LIELL)
             END IF
          ENDIF 
          IF ( NRHS == 1 ) THEN
              CALL strsv( 'U', 'N', 'U', NPIV, A(APOS), LDA,
     &                  W(PTRACB(STEP(INODE))),1)
          ELSE
             CALL strsm( 'L','U', 'N', 'U', NPIV, NRHS, ONE,
     &                   A(APOS), LDA,
     &                   W(PTRACB(STEP(INODE))),LIELL)
          END IF
 1234     CONTINUE   
          IF (KEEP(201).GT.0) THEN
           CALL SMUMPS_598(INODE,PTRFAC,KEEP(28),
     &          A,LA,.TRUE.,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
          ENDIF
          IPOS =   PTRIST(STEP(INODE)) +  KEEP(IXSZ) + 6 + NSLAVES   
          DO I = 1, NPIV
            JJ = IW( IPOS + I - 1 )
            DO K=1,NRHS
              RHS( JJ, K ) = W( PTRACB(STEP(INODE))+I-1
     &         + (K-1)*LIELL )
            ENDDO
          END DO
          IN = INODE
  200     IN = FILS(IN)
          IF (IN .GT. 0) GOTO 200
          IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
              CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &                       FEUILLE, SLAVEF )
              NBFINF = NBFINF - 1
            ENDIF
            IWCB( PTRICB(STEP(INODE)) + 1 ) = 0
            CALL SMUMPS_151(NRHS, N, KEEP(28),
     &          IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
            GOTO 270
          ENDIF  
          DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
          END DO
          IN = -IN
          IF ( KEEP(237).GT.0 ) THEN 
            NO_CHILDREN = .TRUE. 
          ELSE
            NO_CHILDREN = .FALSE.
          ENDIF
          DO WHILE (IN.GT.0) 
            IF ( KEEP(237).GT.0 ) THEN
               IF (.NOT.TO_PROCESS(STEP(IN))) THEN
                  IN = FRERE(STEP(IN))
                  CYCLE
               ELSE
                 NO_CHILDREN = .FALSE.
               ENDIF
            ENDIF
           POOL_FIRST_POS  = IIPOOL
            IF (MUMPS_275(PROCNODE_STEPS(STEP(IN)),
     &          SLAVEF) .EQ. MYID) THEN
                  IPOOL(IIPOOL ) = IN
                  IIPOOL = IIPOOL + 1
            ELSE
              PROCDEST = MUMPS_275( PROCNODE_STEPS(STEP(IN)), 
     &                   SLAVEF )
              IF ( .NOT. DEJA_SEND( PROCDEST ) ) THEN
 110            CALL SMUMPS_78( NRHS, IN, 0, 0,
     &          LIELL, LIELL-KEEP(253),
     &          IW( POSINDICES ) ,
     &          W( PTRACB(STEP(INODE))),
     &          PROCDEST, NOEUD, COMM, IERR )
                IF ( IERR .EQ. -1 ) THEN
                  CALL SMUMPS_41(
     &            .FALSE., FLAG,
     &            BUFR, LBUFR, LBUFR_BYTES,
     &            MYID, SLAVEF, COMM,
     &            N, IWCB, LIWW, POSIWCB,
     &            W, LWC, POSWCB,
     &            IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &            IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &            FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     &            PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &            RHS, LRHS, NRHS, MTYPE, 
     &            RHSCOMP, LRHSCOMP, POSINRHSCOMP
     &            , TO_PROCESS, SIZE_TO_PROCESS
     &            )
                  IF ( INFO( 1 ) .LT. 0 ) GOTO 270
                  GOTO 110
                ELSE IF ( IERR .eq. -2 ) THEN
                  INFO(1) = -17
                  INFO(2) = LIELL * NRHS * KEEP(35) +
     &                    ( LIELL + 2 ) * KEEP(34)
                  GOTO 260
                ELSE IF ( IERR .eq. -3 ) THEN
                  INFO(1) = -20
                  INFO(2) = LIELL * NRHS * KEEP(35) +
     &                    ( LIELL + 2 ) * KEEP(34)
                  GOTO 260
                END IF
                DEJA_SEND( PROCDEST ) = .TRUE.
              END IF
            END IF
            IN = FRERE( STEP( IN ) )
          END DO
          IF (NO_CHILDREN) THEN
                   MYLEAFE = MYLEAFE - 1
                   IF (MYLEAFE .EQ. 0) THEN
                      CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, 
     &                     COMM, FEUILLE, SLAVEF )
                      NBFINF = NBFINF - 1
                   ENDIF
                   IWCB( PTRICB(STEP(INODE)) + 1 ) = 0
                   CALL SMUMPS_151(NRHS, N, KEEP(28),
     &                  IWCB, LIWW, W, LWC,
     &                  POSWCB, POSIWCB, PTRICB, PTRACB)
                   GOTO 270
           ENDIF
          DO I=1,(IIPOOL-POOL_FIRST_POS)/2
           TMP=IPOOL(POOL_FIRST_POS+I-1)
           IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
           IPOOL(IIPOOL-I)=TMP
          ENDDO 
          IWCB( PTRICB(STEP( INODE )) + 1 ) = 0
          CALL SMUMPS_151(NRHS, N, KEEP(28),
     &          IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
        END IF   
      ELSE IF (MSGTAG.EQ.TERREUR) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GO TO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     &      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
      ELSE
          INFO(1) = -100
          INFO(2) = MSGTAG
          GOTO 260
      ENDIF
      GO TO 270
 260  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
 270  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_42
      SUBROUTINE SMUMPS_641(PANEL_SIZE, PANEL_POS,
     &                           LEN_PANEL_POS, INDICES, NPIV,
     &                           NPANELS, NFRONT_OR_NASS,
     &                           NBENTRIES_ALLPANELS)
      IMPLICIT NONE
      INTEGER, intent (in)   :: PANEL_SIZE, NPIV
      INTEGER, intent (in)   :: INDICES(NPIV)
      INTEGER, intent (in)   :: LEN_PANEL_POS
      INTEGER, intent (out)  :: NPANELS
      INTEGER, intent (out)  :: PANEL_POS(LEN_PANEL_POS)
      INTEGER, intent (in)   :: NFRONT_OR_NASS
      INTEGER(8), intent(out):: NBENTRIES_ALLPANELS
      INTEGER NPANELS_MAX, I, NBeff
      INTEGER(8) :: NBENTRIES_THISPANEL
      NBENTRIES_ALLPANELS = 0_8
      NPANELS_MAX = (NPIV+PANEL_SIZE-1)/PANEL_SIZE
      IF (LEN_PANEL_POS .LT. NPANELS_MAX + 1) THEN
        WRITE(*,*) "Error 1 in SMUMPS_641",
     &              LEN_PANEL_POS,NPANELS_MAX
        CALL MUMPS_ABORT()
      ENDIF
      I = 1
      NPANELS = 0
      IF (I .GT. NPIV) RETURN 
 10   CONTINUE
      NPANELS = NPANELS + 1
      PANEL_POS(NPANELS) = I
      NBeff = min(PANEL_SIZE, NPIV-I+1)
      IF ( INDICES(I+NBeff-1) < 0) THEN
        NBeff=NBeff+1
      ENDIF
      NBENTRIES_THISPANEL = int(NFRONT_OR_NASS-I+1,8) * int(NBeff,8)
      NBENTRIES_ALLPANELS = NBENTRIES_ALLPANELS + NBENTRIES_THISPANEL
      I=I+NBeff
      IF ( I .LE. NPIV ) GOTO 10
      PANEL_POS(NPANELS+1)=NPIV+1
      RETURN
      END SUBROUTINE SMUMPS_641
      SUBROUTINE SMUMPS_286( NRHS, DESCA_PAR, 
     &  CNTXT_PAR,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     &  IPIV,LPIV,MASTER_ROOT,MYID,COMM,
     &  RHS_SEQ,SIZE_ROOT,A,INFO,MTYPE,LDLT )
      IMPLICIT NONE
      INTEGER NRHS, MTYPE
      INTEGER DESCA_PAR( 9 )
      INTEGER LOCAL_M, LOCAL_N, MBLOCK, NBLOCK
      INTEGER CNTXT_PAR, MASTER_ROOT, SIZE_ROOT
      INTEGER MYID, COMM
      INTEGER LPIV, IPIV( LPIV )
      INTEGER INFO(40), LDLT
      REAL RHS_SEQ( SIZE_ROOT *NRHS)
      REAL A( LOCAL_M, LOCAL_N )
      INTEGER IERR, NPROW, NPCOL, MYROW, MYCOL
      INTEGER LOCAL_N_RHS
      REAL, ALLOCATABLE, DIMENSION( :,: ) ::RHS_PAR
      EXTERNAL numroc
      INTEGER  numroc
      INTEGER allocok
      CALL blacs_gridinfo( CNTXT_PAR, NPROW, NPCOL, MYROW, MYCOL )
      LOCAL_N_RHS = numroc(NRHS, NBLOCK, MYCOL, 0, NPCOL)
      LOCAL_N_RHS = max(1,LOCAL_N_RHS)
      ALLOCATE(RHS_PAR(LOCAL_M, LOCAL_N_RHS),stat=allocok)
      IF (allocok > 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root.'
        WRITE(*,*) ' Reduce number of right hand sides.'
        CALL MUMPS_ABORT()
      ENDIF
      CALL SMUMPS_290( MYID, SIZE_ROOT, NRHS, RHS_SEQ,
     &      LOCAL_M, LOCAL_N_RHS,
     &      MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     &      NPROW, NPCOL, COMM )
      CALL SMUMPS_768 (SIZE_ROOT, NRHS, MTYPE,
     &     A, DESCA_PAR, LOCAL_M, LOCAL_N, LOCAL_N_RHS,
     &     IPIV, LPIV, RHS_PAR, LDLT, 
     &     MBLOCK, NBLOCK, CNTXT_PAR,
     &     IERR)
      CALL SMUMPS_156( MYID, SIZE_ROOT, NRHS,
     &    RHS_SEQ, LOCAL_M, LOCAL_N_RHS,
     &    MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     &    NPROW, NPCOL, COMM )
      DEALLOCATE(RHS_PAR)
      RETURN
      END SUBROUTINE SMUMPS_286
      SUBROUTINE SMUMPS_768 (SIZE_ROOT, NRHS, MTYPE,
     &     A, DESCA_PAR, LOCAL_M, LOCAL_N, LOCAL_N_RHS,
     &     IPIV, LPIV, RHS_PAR, LDLT, 
     &     MBLOCK, NBLOCK, CNTXT_PAR,
     &     IERR)
      IMPLICIT NONE
      INTEGER, intent (in) :: SIZE_ROOT, NRHS, LDLT, LOCAL_M, 
     &                        LOCAL_N, LOCAL_N_RHS, 
     &                        MBLOCK, NBLOCK, CNTXT_PAR, MTYPE
      INTEGER, intent (in) :: DESCA_PAR( 9 ) 
      INTEGER, intent (in) :: LPIV, IPIV( LPIV )
      REAL, intent (in) :: A( LOCAL_M, LOCAL_N )
      REAL, intent (inout) :: RHS_PAR(LOCAL_M, LOCAL_N_RHS)
      INTEGER, intent (out) :: IERR
      INTEGER              :: DESCB_PAR( 9 )
      IERR = 0
      CALL DESCINIT( DESCB_PAR, SIZE_ROOT, 
     &      NRHS, MBLOCK, NBLOCK, 0, 0,
     &      CNTXT_PAR, LOCAL_M, IERR )
            IF (IERR.NE.0) THEN
              WRITE(*,*) 'After DESCINIT, IERR = ', IERR
              CALL MUMPS_ABORT()
            END IF
      IF ( LDLT .eq. 0 .OR. LDLT .eq. 2 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          CALL psgetrs('N',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     &      RHS_PAR,1,1,DESCB_PAR,IERR)
        ELSE
          CALL psgetrs('T',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     &      RHS_PAR, 1, 1, DESCB_PAR,IERR)
        END IF
      ELSE
        CALL pspotrs( 'L', SIZE_ROOT, NRHS, A, 1, 1, DESCA_PAR,
     &    RHS_PAR, 1, 1, DESCB_PAR, IERR )
      END IF
      IF ( IERR .LT. 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root'
        CALL MUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE SMUMPS_768
