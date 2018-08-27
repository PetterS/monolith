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
      MODULE ZMUMPS_LOAD
      implicit none
      PUBLIC :: ZMUMPS_188, ZMUMPS_185,
     &  ZMUMPS_189, ZMUMPS_190,
     &  ZMUMPS_183, ZMUMPS_187, 
     &  ZMUMPS_186, ZMUMPS_409,
     &  ZMUMPS_384, ZMUMPS_461,
     &  ZMUMPS_467, ZMUMPS_471,
     &  ZMUMPS_472,
     &  ZMUMPS_791, ZMUMPS_790,
     &  ZMUMPS_792, ZMUMPS_500,
     &  ZMUMPS_501, ZMUMPS_520,
     &  ZMUMPS_513,
     &  ZMUMPS_514, ZMUMPS_512
     &  ,ZMUMPS_533,
     &  ZMUMPS_819, ZMUMPS_818,
     &  ZMUMPS_820, ZMUMPS_554,
     &  ZMUMPS_553,
     &  ZMUMPS_555
      DOUBLE PRECISION, DIMENSION(:),
     &       ALLOCATABLE, SAVE, PRIVATE :: LOAD_FLOPS
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PRIVATE :: BUF_LOAD_RECV
      INTEGER, SAVE, PRIVATE :: LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES
      INTEGER, SAVE, PRIVATE :: K50, K69, K35
      INTEGER(8), SAVE, PRIVATE :: MAX_SURF_MASTER
      LOGICAL, SAVE, PRIVATE :: BDC_MEM, BDC_POOL, BDC_SBTR, 
     &     BDC_POOL_MNG,
     &     BDC_M2_MEM,BDC_M2_FLOPS,BDC_MD,REMOVE_NODE_FLAG,
     &     REMOVE_NODE_FLAG_MEM
      DOUBLE PRECISION, SAVE, PRIVATE :: REMOVE_NODE_COST,
     &     REMOVE_NODE_COST_MEM
      INTEGER, SAVE, PRIVATE :: SBTR_WHICH_M
      DOUBLE PRECISION, DIMENSION(:),
     &       ALLOCATABLE, TARGET, SAVE, PRIVATE :: WLOAD
#if defined(OLD_LOAD_MECHANISM)
#if defined(CHECK_COHERENCE) 
      INTEGER, SAVE, PRIVATE :: NB_LEVEL2
      LOGICAL, PRIVATE :: AMI_CHOSEN,IS_DISPLAYED
#endif
#endif
#if ! defined(OLD_LOAD_MECHANISM)
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PRIVATE :: FUTURE_NIV2
      DOUBLE PRECISION, SAVE, PRIVATE :: DELTA_LOAD, DELTA_MEM
#else
      DOUBLE PRECISION, SAVE, PRIVATE :: LAST_LOAD_SENT, 
     &           DM_LAST_MEM_SENT
#endif
      INTEGER(8), SAVE, PRIVATE :: CHECK_MEM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, TARGET, PRIVATE :: 
     &          IDWLOAD
      DOUBLE PRECISION, SAVE, PRIVATE :: COST_SUBTREE
      DOUBLE PRECISION, SAVE, PRIVATE :: ALPHA
      DOUBLE PRECISION, SAVE, PRIVATE :: BETA
      INTEGER, SAVE, PRIVATE :: MYID, NPROCS, COMM_LD
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE, 
     &           PRIVATE :: POOL_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE, 
     &           SAVE :: SBTR_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, 
     &           PRIVATE, SAVE :: SBTR_CUR
      INTEGER, DIMENSION(:), ALLOCATABLE, 
     &           PRIVATE, SAVE :: NB_SON
      DOUBLE PRECISION, 
     &           PRIVATE, SAVE :: SBTR_CUR_LOCAL
      DOUBLE PRECISION, 
     &           PRIVATE, SAVE :: PEAK_SBTR_CUR_LOCAL
      DOUBLE PRECISION, 
     &           PRIVATE, SAVE :: MAX_PEAK_STK
      DOUBLE PRECISION, SAVE, 
     &           PRIVATE :: POOL_LAST_COST_SENT
      DOUBLE PRECISION, SAVE, 
     &           PRIVATE :: MIN_DIFF
      INTEGER, SAVE :: POS_ID,POS_MEM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: CB_COST_ID
      INTEGER(8), DIMENSION(:), ALLOCATABLE, SAVE
     &           :: CB_COST_MEM
      PUBLIC :: CB_COST_ID, CB_COST_MEM,POS_MEM,POS_ID
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: LU_USAGE
      INTEGER(8), DIMENSION(:), ALLOCATABLE, SAVE,
     &        PRIVATE::MD_MEM, TAB_MAXS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE ::MEM_SUBTREE
      INTEGER  :: NB_SUBTREES,NIV1_FLAG
      INTEGER, PRIVATE  :: INDICE_SBTR,INDICE_SBTR_ARRAY
      INTEGER,SAVE :: INSIDE_SUBTREE
      PUBLIC :: NB_SUBTREES,MEM_SUBTREE,INSIDE_SUBTREE,NIV1_FLAG
      DOUBLE PRECISION, SAVE, PRIVATE :: DM_SUMLU,
     &                   DM_THRES_MEM
      DOUBLE PRECISION, DIMENSION(:),
     &   ALLOCATABLE, SAVE , PRIVATE:: DM_MEM
      INTEGER, SAVE, PRIVATE :: POOL_SIZE,ID_MAX_M2
      DOUBLE PRECISION, SAVE, PRIVATE :: MAX_M2,TMP_M2
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE, PRIVATE:: POOL_NIV2
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,SAVE, 
     &      PRIVATE :: POOL_NIV2_COST, NIV2
      DOUBLE PRECISION, SAVE, PRIVATE  ::      CHK_LD
      INTEGER, DIMENSION(:),POINTER, SAVE, PRIVATE  :: 
     &         PROCNODE_LOAD, STEP_TO_NIV2_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE, PRIVATE  :: KEEP_LOAD
      INTEGER, SAVE, PRIVATE :: N_LOAD
      INTEGER(8), DIMENSION(:), POINTER, SAVE, PRIVATE:: KEEP8_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE :: 
     &         FILS_LOAD, STEP_LOAD,
     &         FRERE_LOAD, ND_LOAD,
     &         NE_LOAD,DAD_LOAD
      INTEGER, DIMENSION(:,:),POINTER, SAVE, PRIVATE :: CAND_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE, 
     &         PRIVATE :: MY_FIRST_LEAF,MY_NB_LEAF, MY_ROOT_SBTR
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE, 
     &         PRIVATE ::SBTR_FIRST_POS_IN_POOL
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,SAVE, 
     &         PRIVATE ::SBTR_PEAK_ARRAY,
     &     SBTR_CUR_ARRAY
      DOUBLE PRECISION,DIMENSION(:),POINTER, SAVE :: COST_TRAV
      INTEGER, DIMENSION(:),POINTER, SAVE :: DEPTH_FIRST_LOAD,
     &     DEPTH_FIRST_SEQ_LOAD,SBTR_ID_LOAD
      PUBLIC :: DEPTH_FIRST_LOAD,COST_TRAV, FILS_LOAD,STEP_LOAD,
     &     FRERE_LOAD, ND_LOAD,NE_LOAD,DAD_LOAD,
     &     DEPTH_FIRST_SEQ_LOAD,SBTR_ID_LOAD
      INTEGER, SAVE     :: ROOT_CURRENT_SUBTREE,CURRENT_BEST,
     &     SECOND_CURRENT_BEST
      PUBLIC :: ROOT_CURRENT_SUBTREE,CURRENT_BEST,
     &     SECOND_CURRENT_BEST
      CONTAINS
      SUBROUTINE ZMUMPS_188( COST_SUBTREE_ARG, K64, K66,
     &     MAXS )
      IMPLICIT NONE
      DOUBLE PRECISION COST_SUBTREE_ARG
      INTEGER K64, K66
      INTEGER(8)::MAXS
      DOUBLE PRECISION T64, T66
      T64 = max ( dble(K64), dble(1) )
      T64 = min ( T64, dble(1000)  )
      T66 = max (dble(K66), dble(100))
      MIN_DIFF     =  ( T64 / dble(1000)  )* 
     &                  T66 * dble(1000000)
      DM_THRES_MEM = dble(MAXS/1000_8)
      COST_SUBTREE = COST_SUBTREE_ARG
      END SUBROUTINE ZMUMPS_188
      SUBROUTINE ZMUMPS_791 ( 
     &      INODE, STEP, N, SLAVEF,
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      CAND, ICNTL, COPY_CAND,
     &      NBSPLIT, NUMORG_SPLIT, SLAVES_LIST,
     &      SIZE_SLAVES_LIST 
     &                                    )
      IMPLICIT NONE
       INTEGER, intent(in) :: INODE, N, SIZE_SLAVES_LIST, SLAVEF, 
     &                        KEEP(500)
       INTEGER, intent(in) :: STEP(N), DAD (KEEP(28)), ICNTL(40),
     &                        PROCNODE_STEPS(KEEP(28)), CAND(SLAVEF+1),
     &                        FILS(N)
       INTEGER, intent(out)   :: NBSPLIT, NUMORG_SPLIT
       INTEGER, intent(inout) :: SLAVES_LIST(SIZE_SLAVES_LIST),
     &                           COPY_CAND(SLAVEF+1)
       INTEGER :: IN, LP, II
       INTEGER  MUMPS_810
       EXTERNAL MUMPS_810
       LP = ICNTL(1) 
       IN = INODE
       NBSPLIT = 0
       NUMORG_SPLIT = 0    
       DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )  
           NBSPLIT = NBSPLIT + 1
           IN = DAD(STEP(IN))
           II = IN
           DO WHILE (II.GT.0)
             NUMORG_SPLIT = NUMORG_SPLIT + 1
             II = FILS(II)
           ENDDO
       END DO
      SLAVES_LIST(1:NBSPLIT) = CAND(1:NBSPLIT)
      COPY_CAND(1:SIZE_SLAVES_LIST-NBSPLIT) = 
     &                   CAND(1+NBSPLIT:SIZE_SLAVES_LIST)
      COPY_CAND(SIZE_SLAVES_LIST-NBSPLIT+1:SLAVEF) = -1
      COPY_CAND(SLAVEF+1) = SIZE_SLAVES_LIST-NBSPLIT
      RETURN
      END SUBROUTINE ZMUMPS_791
      SUBROUTINE ZMUMPS_790 ( 
     &      INODE, STEP, N, SLAVEF, NBSPLIT, NCB, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS, ICNTL,
     &      TAB_POS, NSLAVES_NODE
     &                                    )
      IMPLICIT NONE
       INTEGER, intent(in) :: INODE, N, SLAVEF, NCB, 
     &                        KEEP(500), NBSPLIT
       INTEGER, intent(in) :: STEP(N), DAD (KEEP(28)), ICNTL(40),
     &                        PROCNODE_STEPS(KEEP(28)), 
     &                        FILS(N)
       INTEGER, intent(inout) :: TAB_POS ( SLAVEF+2 ), NSLAVES_NODE
       INTEGER :: IN, LP, II, NUMORG, NBSPLIT_LOC, I
       INTEGER  MUMPS_810
       EXTERNAL MUMPS_810
       DO I= NSLAVES_NODE+1, 1, -1
          TAB_POS(I+NBSPLIT) = TAB_POS(I) 
       END DO
       LP = ICNTL(1) 
       IN = INODE
       NBSPLIT_LOC = 0
       NUMORG = 0
       TAB_POS(1) = 1
       DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )  
           NBSPLIT_LOC = NBSPLIT_LOC + 1
           IN = DAD(STEP(IN))
           II = IN
           DO WHILE (II.GT.0)
             NUMORG = NUMORG + 1
             II = FILS(II)
           ENDDO
           TAB_POS(NBSPLIT_LOC+1) = NUMORG + 1
       END DO
       DO I = NBSPLIT+2, NBSPLIT+NSLAVES_NODE+1
         TAB_POS(I) = TAB_POS(I) + NUMORG
       ENDDO
      NSLAVES_NODE = NSLAVES_NODE + NBSPLIT
      TAB_POS (NSLAVES_NODE+2:SLAVEF+1) = -9999
      TAB_POS ( SLAVEF+2 ) =  NSLAVES_NODE
      RETURN
      END SUBROUTINE ZMUMPS_790
      SUBROUTINE ZMUMPS_792 (
     &      INODE, TYPESPLIT, IFSON, 
     &      SON_SLAVE_LIST, NSLSON,
     &      STEP, N, SLAVEF, NBSPLIT, NCB, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS, ICNTL,
     &      ISTEP_TO_INIV2, INIV2,
     &      TAB_POS_IN_PERE, NSLAVES_NODE,
     &      SLAVES_LIST, SIZE_SLAVES_LIST
     &                                    )
      IMPLICIT NONE
       INTEGER, intent(in) :: INODE, TYPESPLIT, IFSON, N, SLAVEF, 
     &                        NCB, KEEP(500), NBSPLIT,
     &                        NSLSON, SIZE_SLAVES_LIST
       INTEGER, intent(in) :: STEP(N), DAD (KEEP(28)), ICNTL(40),
     &                        PROCNODE_STEPS(KEEP(28)), 
     &                        FILS(N), INIV2,
     &                        SON_SLAVE_LIST (NSLSON),
     &                        ISTEP_TO_INIV2(KEEP(71))
       INTEGER, intent(out)   ::  NSLAVES_NODE
       INTEGER, intent(inout) ::
     &                   TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
       INTEGER, intent(out)   :: SLAVES_LIST (SIZE_SLAVES_LIST)
       INTEGER :: IN, LP, I, NSLAVES_SONS,
     &            INIV2_FILS, ISHIFT
       LP = ICNTL(1) 
       IN = INODE
      INIV2_FILS = ISTEP_TO_INIV2( STEP( IFSON ))
      NSLAVES_SONS = TAB_POS_IN_PERE (SLAVEF+2, INIV2_FILS)
      TAB_POS_IN_PERE (1,INIV2) = 1
      ISHIFT  = TAB_POS_IN_PERE (2, INIV2_FILS) -1 
      DO I = 2, NSLAVES_SONS
         TAB_POS_IN_PERE (I,INIV2) = 
     &            TAB_POS_IN_PERE (I+1,INIV2_FILS) - ISHIFT
         SLAVES_LIST(I-1) =  SON_SLAVE_LIST (I)
      END DO
      TAB_POS_IN_PERE(NSLAVES_SONS+1:SLAVEF+1,INIV2) = -9999
      NSLAVES_NODE = NSLAVES_SONS - 1
      TAB_POS_IN_PERE (SLAVEF+2, INIV2) = NSLAVES_NODE
      RETURN
      END SUBROUTINE ZMUMPS_792 
      SUBROUTINE ZMUMPS_472(
     &  NCBSON_MAX, SLAVEF,
     &  KEEP,KEEP8,ICNTL,
     &  CAND_OF_NODE,
     &  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,INODE)
       IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: ICNTL(40)
      INTEGER, intent(in) :: SLAVEF, NFRONT
      INTEGER, intent (inout) ::NCB   
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(in) :: NCBSON_MAX
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER i
      INTEGER LP,MP
      LP=ICNTL(4)
      MP=ICNTL(2)
      IF ( KEEP(48) == 0 .OR. KEEP(48) .EQ. 3 ) THEN
         CALL ZMUMPS_499(
     &        SLAVEF,
     &        KEEP,KEEP8,
     &        CAND_OF_NODE,
     &        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      ELSE IF ( KEEP(48) == 4 ) THEN
         CALL ZMUMPS_504(
     &        SLAVEF,
     &        KEEP,KEEP8,
     &        CAND_OF_NODE,
     &        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID)
         DO i=1,NSLAVES_NODE
            IF(TAB_POS(i+1)-TAB_POS(i).LE.0)THEN
               WRITE(*,*)'probleme de partition dans 
     &ZMUMPS_545'
               CALL MUMPS_ABORT()
            ENDIF
         ENDDO
      ELSE IF ( KEEP(48) == 5 ) THEN
         CALL ZMUMPS_518(
     &        NCBSON_MAX,
     &        SLAVEF,
     &        KEEP,KEEP8,
     &        CAND_OF_NODE,
     &        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID,INODE,
     &        MP,LP)
         DO i=1,NSLAVES_NODE
            IF(TAB_POS(i+1)-TAB_POS(i).LE.0)THEN
               WRITE(*,*)'problem with partition in
     &ZMUMPS_518'
               CALL MUMPS_ABORT()
            ENDIF
         ENDDO
      ELSE
        WRITE(*,*) "Strategy 6 not implemented"
        CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_472
      SUBROUTINE ZMUMPS_499(
     &  SLAVEF,
     &  KEEP,KEEP8,
     &  CAND_OF_NODE,
     &  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER ITEMP, NMB_OF_CAND, NSLAVES_LESS
      DOUBLE PRECISION MSG_SIZE
      LOGICAL FORCE_CAND
      INTEGER  MUMPS_12
      EXTERNAL MUMPS_12
      IF ( KEEP(48) == 0 .AND. KEEP(50) .NE. 0) THEN
      write(*,*) "Internal error 2 in ZMUMPS_499."
      CALL MUMPS_ABORT()
      END IF
      IF ( KEEP(48) == 3 .AND. KEEP(50) .EQ. 0) THEN
      write(*,*) "Internal error 3 in ZMUMPS_499."
      CALL MUMPS_ABORT()
      END IF
      MSG_SIZE = dble( NFRONT - NCB ) * dble(NCB)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
      END IF
      IF (FORCE_CAND) THEN
        ITEMP=ZMUMPS_409
     &       (MEM_DISTRIB,
     &        CAND_OF_NODE,
     &
     &        KEEP(69), SLAVEF, MSG_SIZE,
     &        NMB_OF_CAND )
      ELSE
        ITEMP=ZMUMPS_186(KEEP(69),MEM_DISTRIB,MSG_SIZE)
        NMB_OF_CAND = SLAVEF - 1
      END IF
      NSLAVES_LESS = max(ITEMP,1)
      NSLAVES_NODE = MUMPS_12(KEEP8(21), KEEP(48),
     &          KEEP(50),SLAVEF,
     &          NCB, NFRONT, NSLAVES_LESS, NMB_OF_CAND)
      CALL MUMPS_441(
     &            KEEP,KEEP8, SLAVEF,
     &            TAB_POS,
     &            NSLAVES_NODE, NFRONT, NCB
     &             )
      IF (FORCE_CAND) THEN
        CALL ZMUMPS_384(MEM_DISTRIB(0),
     &       CAND_OF_NODE, SLAVEF, NSLAVES_NODE,
     &       SLAVES_LIST)
      ELSE
        CALL ZMUMPS_189(MEM_DISTRIB(0),
     &       MSG_SIZE, SLAVES_LIST, NSLAVES_NODE)
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_499
      SUBROUTINE ZMUMPS_185( id, MEMORY_MD_ARG, MAXS )
      USE ZMUMPS_COMM_BUFFER
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(ZMUMPS_STRUC), TARGET :: id
      INTEGER(8), intent(in) :: MEMORY_MD_ARG
      INTEGER(8), intent(in) :: MAXS
      INTEGER K34_LOC,K35_LOC
      INTEGER allocok, IERR, i, BUF_LOAD_SIZE
      DOUBLE PRECISION :: MAX_SBTR
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION MEMORY_SENT
      INTEGER,DIMENSION(:),POINTER:: KEEP
      PARAMETER( ZERO=0.0d0 )
      INTEGER WHAT
      INTEGER(8) MEMORY_MD, LA
      STEP_TO_NIV2_LOAD=>id%ISTEP_TO_INIV2
      CAND_LOAD=>id%CANDIDATES
      ND_LOAD=>id%ND_STEPS
      KEEP_LOAD=>id%KEEP
      KEEP     =>id%KEEP
      KEEP8_LOAD=>id%KEEP8
      FILS_LOAD=>id%FILS
      FRERE_LOAD=>id%FRERE_STEPS
      DAD_LOAD=>id%DAD_STEPS
      PROCNODE_LOAD=>id%PROCNODE_STEPS
      STEP_LOAD=>id%STEP
      NE_LOAD=>id%NE_STEPS
      N_LOAD=id%N
      ROOT_CURRENT_SUBTREE=-9999
      MEMORY_MD=MEMORY_MD_ARG
      LA=MAXS
      MAX_SURF_MASTER=id%MAX_SURF_MASTER+
     & (int(id%KEEP(12),8)*int(id%MAX_SURF_MASTER,8)/int(100,8))
      COMM_LD = id%COMM_LOAD
      MAX_PEAK_STK = 0.0D0
      K69  = KEEP(69)
      IF ( KEEP(47) .le. 0 .OR. KEEP(47) .gt. 4 ) THEN
        write(*,*) "Internal error 1 in ZMUMPS_185"
        CALL MUMPS_ABORT()
      END IF
      CHK_LD=dble(0)
      BDC_MEM      = ( KEEP(47) >= 2 )
      BDC_POOL     = ( KEEP(47) >= 3 )
      BDC_SBTR     = ( KEEP(47) >= 4 )
      BDC_M2_MEM   = ( ( KEEP(80) == 2 .OR. KEEP(80) == 3 )
     &             .AND. KEEP(47) == 4 )
      BDC_M2_FLOPS   = ( KEEP(80) == 1 
     &             .AND. KEEP(47) .GE. 1 )
      BDC_MD       = (KEEP(86)==1)
      SBTR_WHICH_M       = KEEP(90)
      REMOVE_NODE_FLAG=.FALSE.
      REMOVE_NODE_FLAG_MEM=.FALSE.
      REMOVE_NODE_COST_MEM=dble(0)
      REMOVE_NODE_COST=dble(0)
      IF (KEEP(80) .LT. 0 .OR. KEEP(80)>3) THEN
        WRITE(*,*) "Unimplemented KEEP(80) Strategy"
        CALL MUMPS_ABORT()
      ENDIF
      IF ((KEEP(80) == 2 .OR. KEEP(80)==3 ).AND. KEEP(47).NE.4)
     &  THEN
        WRITE(*,*) "Internal error 3 in ZMUMPS_185"
         CALL MUMPS_ABORT()
      END IF
      IF (KEEP(81) == 1 .AND. KEEP(47) < 2) THEN
        WRITE(*,*) "Internal error 2 in ZMUMPS_185"
        CALL MUMPS_ABORT()
      ENDIF
      BDC_POOL_MNG = ((KEEP(81) == 1).AND.(KEEP(47) >= 2))
      IF(KEEP(76).EQ.4)THEN
         DEPTH_FIRST_LOAD=>id%DEPTH_FIRST
      ENDIF
      IF(KEEP(76).EQ.5)THEN
         COST_TRAV=>id%COST_TRAV
      ENDIF
      IF(KEEP(76).EQ.6)THEN
         DEPTH_FIRST_LOAD=>id%DEPTH_FIRST
         DEPTH_FIRST_SEQ_LOAD=>id%DEPTH_FIRST_SEQ
         SBTR_ID_LOAD=>id%SBTR_ID
      ENDIF
      IF (BDC_M2_MEM.OR.BDC_M2_FLOPS) THEN
         ALLOCATE(NIV2(id%NSLAVES), NB_SON(KEEP(28)),
     &            POOL_NIV2(100),POOL_NIV2_COST(100),
     &            stat=allocok)
         NB_SON=id%NE_STEPS
         NIV2=dble(0)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in ZMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES + KEEP(28) + 200
           RETURN
         ENDIF
      ENDIF
      K50      = id%KEEP(50)
      CALL MPI_COMM_RANK( COMM_LD, MYID, IERR )
      NPROCS = id%NSLAVES
      DM_SUMLU=ZERO
      POOL_SIZE=0
      IF(BDC_MD)THEN
         IF ( allocated(MD_MEM) ) DEALLOCATE(MD_MEM)
         ALLOCATE( MD_MEM( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in ZMUMPS_185'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         IF ( allocated(TAB_MAXS) ) DEALLOCATE(TAB_MAXS)
         ALLOCATE( TAB_MAXS( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in ZMUMPS_185'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         TAB_MAXS=0_8
         IF ( allocated(LU_USAGE) ) DEALLOCATE(LU_USAGE)
         ALLOCATE( LU_USAGE( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in ZMUMPS_185'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         LU_USAGE=dble(0)
         MD_MEM=int(0,8)
      ENDIF
      IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
         ALLOCATE(CB_COST_MEM(2*2000*id%NSLAVES), 
     &            stat=allocok)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in ZMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES 
           RETURN
         ENDIF
         CB_COST_MEM=int(0,8)
         ALLOCATE(CB_COST_ID(2000*3), 
     &            stat=allocok)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in ZMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES 
           RETURN
         ENDIF
         CB_COST_ID=0
         POS_MEM=1
         POS_ID=1
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      ALLOCATE(FUTURE_NIV2(NPROCS), stat=allocok)
      IF (allocok > 0 ) THEN
         WRITE(*,*) 'PB allocation in ZMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      ENDIF
      DO i = 1, NPROCS
        FUTURE_NIV2(i) = id%FUTURE_NIV2(i)
        IF(BDC_MD)THEN
           IF(FUTURE_NIV2(i).EQ.0)THEN
              MD_MEM(i-1)=999999999_8
           ENDIF
        ENDIF
      ENDDO
      DELTA_MEM=ZERO
      DELTA_LOAD=ZERO
#endif
      CHECK_MEM=0_8
#if defined(OLD_LOAD_MECHANISM)
#if defined(CHECK_COHERENCE) 
      NB_LEVEL2=0
      AMI_CHOSEN=.FALSE.
      IS_DISPLAYED=.FALSE.
#endif
#endif
      IF(BDC_SBTR.OR.BDC_POOL_MNG)THEN
         NB_SUBTREES=id%NBSA_LOCAL
         IF (allocated(MEM_SUBTREE)) DEALLOCATE(MEM_SUBTREE)
         ALLOCATE(MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
         DO i=1,id%NBSA_LOCAL
            MEM_SUBTREE(i)=id%MEM_SUBTREE(i)
         ENDDO
         MY_FIRST_LEAF=>id%MY_FIRST_LEAF
         MY_NB_LEAF=>id%MY_NB_LEAF
         MY_ROOT_SBTR=>id%MY_ROOT_SBTR
         IF (allocated(SBTR_FIRST_POS_IN_POOL))
     &        DEALLOCATE(SBTR_FIRST_POS_IN_POOL)
         ALLOCATE(SBTR_FIRST_POS_IN_POOL(id%NBSA_LOCAL),stat=allocok)
         INSIDE_SUBTREE=0
         PEAK_SBTR_CUR_LOCAL = dble(0)
         SBTR_CUR_LOCAL      = dble(0)
         IF (allocated(SBTR_PEAK_ARRAY)) DEALLOCATE(SBTR_PEAK_ARRAY)
         ALLOCATE(SBTR_PEAK_ARRAY(id%NBSA_LOCAL),stat=allocok)
         SBTR_PEAK_ARRAY=dble(0)
         IF (allocated(SBTR_CUR_ARRAY)) DEALLOCATE(SBTR_CUR_ARRAY)
         ALLOCATE(SBTR_CUR_ARRAY(id%NBSA_LOCAL),stat=allocok)
         SBTR_CUR_ARRAY=dble(0)
         INDICE_SBTR_ARRAY=1
         NIV1_FLAG=0
         INDICE_SBTR=1
      ENDIF
      IF ( allocated(LOAD_FLOPS) ) DEALLOCATE( LOAD_FLOPS )
      ALLOCATE( LOAD_FLOPS( 0: NPROCS - 1 ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in ZMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( allocated(WLOAD) ) DEALLOCATE( WLOAD )
      ALLOCATE( WLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in ZMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( allocated(IDWLOAD) ) DEALLOCATE( IDWLOAD )
      ALLOCATE( IDWLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in ZMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( BDC_MEM ) THEN
        IF ( allocated(DM_MEM) ) DEALLOCATE( DM_MEM )
        ALLOCATE( DM_MEM( 0:NPROCS-1 ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in ZMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
      END IF
      IF ( BDC_POOL ) THEN
        IF ( allocated(POOL_MEM) ) DEALLOCATE(POOL_MEM)
        ALLOCATE( POOL_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in ZMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        POOL_MEM = dble(0) 
        POOL_LAST_COST_SENT = dble(0)
      END IF
      IF ( BDC_SBTR ) THEN
        IF ( allocated(SBTR_MEM) ) DEALLOCATE(SBTR_MEM)
        ALLOCATE( SBTR_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in ZMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        IF ( allocated(SBTR_CUR) ) DEALLOCATE(SBTR_CUR)
        ALLOCATE( SBTR_CUR(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in ZMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        SBTR_CUR = dble(0) 
        SBTR_MEM = dble(0) 
      END IF
      CALL MUMPS_546(K34_LOC,K35_LOC)
      K35  = K35_LOC
      BUF_LOAD_SIZE = K34_LOC * 2 * ( NPROCS - 1 ) +
     &                NPROCS * ( K35_LOC + K34_LOC )
      IF (BDC_MEM) THEN
        BUF_LOAD_SIZE = BUF_LOAD_SIZE + NPROCS * K35_LOC
      END IF
      IF (BDC_SBTR)THEN
        BUF_LOAD_SIZE = BUF_LOAD_SIZE + NPROCS * K35_LOC
      ENDIF
      LBUF_LOAD_RECV = (BUF_LOAD_SIZE+K34_LOC)/K34_LOC
      LBUF_LOAD_RECV_BYTES = LBUF_LOAD_RECV * K34_LOC
      IF ( allocated(BUF_LOAD_RECV) ) DEALLOCATE(BUF_LOAD_RECV)
      ALLOCATE( BUF_LOAD_RECV( LBUF_LOAD_RECV), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in ZMUMPS_185'
        id%INFO(1) = -13
        id%INFO(2) = LBUF_LOAD_RECV
        RETURN
      ENDIF 
      BUF_LOAD_SIZE = BUF_LOAD_SIZE * 20
      CALL ZMUMPS_54( BUF_LOAD_SIZE, IERR )
      IF ( IERR .LT. 0 ) THEN
         id%INFO(1) = -13
         id%INFO(2) = BUF_LOAD_SIZE
         RETURN
      END IF
      DO i = 0, NPROCS - 1
         LOAD_FLOPS( i ) = ZERO
      END DO
#if defined(OLD_LOAD_MECHANISM)
      LOAD_FLOPS( MYID ) = COST_SUBTREE
      LAST_LOAD_SENT = ZERO
#endif
      IF ( BDC_MEM ) THEN
        DO i = 0, NPROCS - 1
          DM_MEM( i )=ZERO
        END DO
#if defined(OLD_LOAD_MECHANISM)
        DM_LAST_MEM_SENT=ZERO
#endif
      ENDIF
      CALL ZMUMPS_425(KEEP(69))
      IF(BDC_MD)THEN
         MAX_SBTR=0.0D0
         IF(BDC_SBTR)THEN
            DO i=1,id%NBSA_LOCAL
               MAX_SBTR=max(id%MEM_SUBTREE(i),MAX_SBTR)
            ENDDO
         ENDIF
         MD_MEM(MYID)=MEMORY_MD
         WHAT=8
         CALL ZMUMPS_460( WHAT,
     &        COMM_LD, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     &        FUTURE_NIV2,
#endif
     &        dble(MEMORY_MD),dble(0) ,MYID, IERR  )
         WHAT=9
         MEMORY_SENT = dble(LA-MAX_SURF_MASTER)-MAX_SBTR
     &      - max( dble(LA) * dble(3) / dble(100),
     &      dble(2) * dble(max(KEEP(5),KEEP(6))) * dble(KEEP(127)))
         IF (KEEP(12) > 25) THEN
           MEMORY_SENT = MEMORY_SENT -
     &                   dble(KEEP(12))*0.2d0*dble(LA)/100.0d0
         ENDIF
         TAB_MAXS(MYID)=int(MEMORY_SENT,8)
         CALL ZMUMPS_460( WHAT,
     &        COMM_LD, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     &        FUTURE_NIV2,
#endif
     &        MEMORY_SENT,
     &        dble(0),MYID, IERR  )
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_185
      SUBROUTINE ZMUMPS_190( CHECK_FLOPS,PROCESS_BANDE,
     &     INC_LOAD, KEEP,KEEP8 )
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      DOUBLE PRECISION INC_LOAD
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      LOGICAL PROCESS_BANDE
      INTEGER CHECK_FLOPS
      INTEGER IERR
      DOUBLE PRECISION ZERO, SEND_MEM, SEND_LOAD,SBTR_TMP
      PARAMETER( ZERO=0.0d0 )
      INTRINSIC max, abs
      IF (INC_LOAD == 0.0D0) THEN
         IF(REMOVE_NODE_FLAG)THEN
            REMOVE_NODE_FLAG=.FALSE.
         ENDIF
         RETURN
      ENDIF
      IF((CHECK_FLOPS.NE.0).AND.
     &     (CHECK_FLOPS.NE.1).AND.(CHECK_FLOPS.NE.2))THEN
         WRITE(*,*)MYID,': Bad value for CHECK_FLOPS'
         CALL MUMPS_ABORT()
      ENDIF
      IF(CHECK_FLOPS.EQ.1)THEN
         CHK_LD=CHK_LD+INC_LOAD
      ELSE 
         IF(CHECK_FLOPS.EQ.2)THEN
            RETURN
         ENDIF
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( PROCESS_BANDE ) THEN
         RETURN                 
      ENDIF
#endif
      LOAD_FLOPS( MYID ) = max( LOAD_FLOPS( MYID ) + INC_LOAD, ZERO)
      IF(BDC_M2_FLOPS.AND.REMOVE_NODE_FLAG)THEN
         IF(INC_LOAD.NE.REMOVE_NODE_COST)THEN
            IF(INC_LOAD.GT.REMOVE_NODE_COST)THEN
#if ! defined(OLD_LOAD_MECHANISM)
               DELTA_LOAD = DELTA_LOAD +
     &              (INC_LOAD-REMOVE_NODE_COST)
               GOTO 888
#else
               GOTO 888
#endif               
            ELSE
#if ! defined(OLD_LOAD_MECHANISM)
               DELTA_LOAD = DELTA_LOAD -
     &              (REMOVE_NODE_COST-INC_LOAD)
               GOTO 888
#else
               GOTO 888
#endif               
            ENDIF
         ENDIF
         GOTO 333
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      DELTA_LOAD = DELTA_LOAD + INC_LOAD
 888  CONTINUE
      IF ( DELTA_LOAD > MIN_DIFF .OR. DELTA_LOAD < -MIN_DIFF) THEN
         SEND_LOAD = DELTA_LOAD
         IF (BDC_MEM) THEN
           SEND_MEM = DELTA_MEM
         ELSE
           SEND_MEM = ZERO
         END IF
#else
 888   CONTINUE
       IF ( abs( LOAD_FLOPS ( MYID ) -
     &      LAST_LOAD_SENT ).GT.MIN_DIFF)THEN
          IERR = 0
          SEND_LOAD  = LOAD_FLOPS( MYID )
          IF ( BDC_MEM ) THEN
             SEND_MEM = DM_MEM(MYID)
          ELSE
             SEND_MEM = ZERO
          END IF
#endif
         IF(BDC_SBTR)THEN
           SBTR_TMP=SBTR_CUR(MYID)
         ELSE
           SBTR_TMP=dble(0)
         ENDIF
 111     CONTINUE
         CALL ZMUMPS_77( BDC_SBTR,BDC_MEM,
     &        BDC_MD,COMM_LD, NPROCS,
     &        SEND_LOAD,
     &        SEND_MEM,SBTR_TMP,
     &        DM_SUMLU,
#if ! defined(OLD_LOAD_MECHANISM)
     &        FUTURE_NIV2,
#endif
     &        MYID, IERR )
           IF ( IERR == -1 )THEN
             CALL ZMUMPS_467(COMM_LD, KEEP)
             GOTO 111
           ELSE IF ( IERR .NE.0 ) THEN
             WRITE(*,*) "Internal Error in ZMUMPS_190",IERR
             CALL MUMPS_ABORT()
           ENDIF
         IF ( IERR .EQ. 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
           DELTA_LOAD = ZERO
           IF (BDC_MEM) DELTA_MEM  = ZERO
#else          
           LAST_LOAD_SENT = LOAD_FLOPS( MYID )
           IF ( BDC_MEM ) DM_LAST_MEM_SENT = DM_MEM( MYID )
#endif
         END IF
      ENDIF
 333  CONTINUE
      IF(REMOVE_NODE_FLAG)THEN
         REMOVE_NODE_FLAG=.FALSE.
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_190
      SUBROUTINE ZMUMPS_471( SSARBR,
     &           PROCESS_BANDE_ARG, MEM_VALUE, NEW_LU, INC_MEM_ARG,
     &           KEEP,KEEP8,LRLU)
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER(8), INTENT(IN) :: MEM_VALUE, INC_MEM_ARG, NEW_LU,LRLU
      LOGICAL, INTENT(IN) :: PROCESS_BANDE_ARG, SSARBR
      INTEGER IERR, KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION ZERO, SEND_MEM, SBTR_TMP
      PARAMETER( ZERO=0.0d0 )
      INTRINSIC max, abs
      INTEGER(8) :: INC_MEM
      LOGICAL PROCESS_BANDE
#if defined(OLD_LOAD_MECHANISM)
      DOUBLE PRECISION TMP_MEM
#endif
      PROCESS_BANDE=PROCESS_BANDE_ARG
      INC_MEM = INC_MEM_ARG
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( PROCESS_BANDE .AND. NEW_LU .NE. 0_8) THEN
        WRITE(*,*) " Internal Error in ZMUMPS_471."
        WRITE(*,*) " NEW_LU must be zero if called from PROCESS_BANDE"
        CALL MUMPS_ABORT()
      ENDIF
#endif
#if defined(OLD_LOAD_MECHANISM)
#if defined(CHECK_COHERENCE) 
      IF(PROCESS_BANDE)THEN
         PROCESS_BANDE=.FALSE.
         NB_LEVEL2=NB_LEVEL2-1
         IF(NB_LEVEL2.LT.0)THEN
            WRITE(*,*)MYID,': problem with NB_LEVEL2'
         ELSEIF(NB_LEVEL2.EQ.0)THEN
            IF(IS_DISPLAYED)THEN
#if defined(STATS_DYNAMIC_MEMORY)
               WRITE(*,*)MYID,': end of Incoherent state at time=',
     &              MPI_WTIME()-TIME_REF
#endif
               IS_DISPLAYED=.FALSE.
            ENDIF
            AMI_CHOSEN=.FALSE.
         ENDIF
      ENDIF
      IF((KEEP(73).EQ.0).AND.(NB_LEVEL2.NE.0)
     &     .AND.(.NOT.IS_DISPLAYED))THEN
         IS_DISPLAYED=.TRUE.
#if defined(STATS_DYNAMIC_MEMORY)
         WRITE(*,*)MYID,': Begin of Incoherent state (1) at time=',
     &        MPI_WTIME()-TIME_REF
#endif
      ENDIF
#endif
#endif
      DM_SUMLU = DM_SUMLU + dble(NEW_LU)
      IF(KEEP_LOAD(201).EQ.0)THEN
         CHECK_MEM = CHECK_MEM + INC_MEM
      ELSE 
         CHECK_MEM = CHECK_MEM + INC_MEM - NEW_LU
      ENDIF
      IF ( MEM_VALUE .NE. CHECK_MEM ) THEN
         WRITE(*,*)MYID,
     &   ':Problem with increments in ZMUMPS_471',
     &   CHECK_MEM, MEM_VALUE, INC_MEM,NEW_LU
         CALL MUMPS_ABORT()
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (PROCESS_BANDE) THEN
         RETURN
      ENDIF
#endif
      IF(BDC_POOL_MNG) THEN
         IF(SBTR_WHICH_M.EQ.0)THEN
            IF (SSARBR) SBTR_CUR_LOCAL = SBTR_CUR_LOCAL+
     &                                   dble(INC_MEM-NEW_LU)
         ELSE
            IF (SSARBR) SBTR_CUR_LOCAL = SBTR_CUR_LOCAL+
     &                                   dble(INC_MEM)
         ENDIF
      ENDIF
      IF ( .NOT. BDC_MEM ) THEN
         RETURN
      ENDIF
#if defined(OLD_LOAD_MECHANISM)
      IF(KEEP_LOAD(201).EQ.0)THEN
         DM_MEM( MYID ) = dble(CHECK_MEM) - DM_SUMLU
      ELSE
         DM_MEM( MYID ) = dble(CHECK_MEM)
      ENDIF
      TMP_MEM = DM_MEM(MYID)
#endif
      IF (BDC_SBTR .AND. SSARBR) THEN
         IF((SBTR_WHICH_M.EQ.0).AND.(KEEP(201).NE.0))THEN
            SBTR_CUR(MYID) = SBTR_CUR(MYID)+dble(INC_MEM-NEW_LU)
         ELSE
            SBTR_CUR(MYID) = SBTR_CUR(MYID)+dble(INC_MEM)
         ENDIF
         SBTR_TMP = SBTR_CUR(MYID)
      ELSE
        SBTR_TMP=dble(0)
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( NEW_LU > 0_8 ) THEN
        INC_MEM = INC_MEM - NEW_LU
      ENDIF
      DM_MEM( MYID ) = DM_MEM(MYID) + dble(INC_MEM)
      MAX_PEAK_STK=max(MAX_PEAK_STK,DM_MEM(MYID))
      IF(BDC_M2_MEM.AND.REMOVE_NODE_FLAG_MEM)THEN
         IF(dble(INC_MEM).NE.REMOVE_NODE_COST_MEM)THEN
            IF(dble(INC_MEM).GT.REMOVE_NODE_COST_MEM)THEN
               DELTA_MEM = DELTA_MEM +
     &              (dble(INC_MEM)-REMOVE_NODE_COST_MEM)
               GOTO 888               
            ELSE
               DELTA_MEM = DELTA_MEM -
     &              (REMOVE_NODE_COST_MEM-dble(INC_MEM))
               GOTO 888
            ENDIF
         ENDIF
         GOTO 333
      ENDIF
      DELTA_MEM = DELTA_MEM + dble(INC_MEM)
 888  CONTINUE
      IF ((KEEP(48).NE.5).OR.
     &     ((KEEP(48).EQ.5).AND.(abs(DELTA_MEM)
     &      .GE.0.1d0*dble(LRLU))))THEN
         IF ( abs(DELTA_MEM) > DM_THRES_MEM ) THEN
            SEND_MEM = DELTA_MEM
#else
      IF(BDC_M2_MEM.AND.REMOVE_NODE_FLAG_MEM)THEN
         IF(INC_MEM.NE.REMOVE_NODE_COST_MEM)THEN
            IF(INC_MEM.EQ.REMOVE_NODE_COST_MEM)THEN               
               GOTO 333
            ELSEIF(INC_MEM.LT.REMOVE_NODE_COST_MEM)THEN
               GOTO 333
            ENDIF
         ENDIF
      ENDIF
      IF ((KEEP(48).NE.5).OR.
     &     ((KEEP(48).EQ.5).AND.
     &      (abs(TMP_MEM-DM_LAST_MEM_SENT).GE.
     &     0.1d0*dble(LRLU))))THEN
         IF ( abs( TMP_MEM-DM_LAST_MEM_SENT) >
     &     DM_THRES_MEM )  THEN
            IERR = 0
            SEND_MEM = TMP_MEM
#endif
 111        CONTINUE
            CALL ZMUMPS_77( 
     &           BDC_SBTR,
     &           BDC_MEM,BDC_MD, COMM_LD,
     &           NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     &           DELTA_LOAD,
#else
     &           LOAD_FLOPS( MYID ), 
#endif
     &           SEND_MEM,SBTR_TMP,
     &           DM_SUMLU,
#if ! defined(OLD_LOAD_MECHANISM)
     &           FUTURE_NIV2,
#endif
     &           MYID,IERR )
            IF ( IERR == -1 )THEN
               CALL ZMUMPS_467(COMM_LD, KEEP)
               GOTO 111
           ELSE IF ( IERR .NE. 0 ) THEN
             WRITE(*,*) "Internal Error in ZMUMPS_471",IERR
             CALL MUMPS_ABORT()
           ENDIF
           IF ( IERR .EQ. 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
             DELTA_LOAD = ZERO
             DELTA_MEM  = ZERO
#else
             LAST_LOAD_SENT   = LOAD_FLOPS  ( MYID )
             DM_LAST_MEM_SENT = TMP_MEM
#endif
           END IF
         ENDIF
      ENDIF
 333  CONTINUE
      IF(REMOVE_NODE_FLAG_MEM)THEN
         REMOVE_NODE_FLAG_MEM=.FALSE.
      ENDIF
      END SUBROUTINE ZMUMPS_471
      INTEGER FUNCTION ZMUMPS_186( K69, MEM_DISTRIB,MSG_SIZE )
      IMPLICIT NONE
      INTEGER i, NLESS, K69 
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      DOUBLE PRECISION LREF
      DOUBLE PRECISION MSG_SIZE
      NLESS = 0
      DO i=1,NPROCS
         IDWLOAD(i) = i - 1
      ENDDO
      WLOAD(1:NPROCS) = LOAD_FLOPS(0:NPROCS-1)
      IF(BDC_M2_FLOPS)THEN
         DO i=1,NPROCS
            WLOAD(i)=WLOAD(i)+NIV2(i)
         ENDDO
      ENDIF
      IF(K69 .gt. 1) THEN
         CALL ZMUMPS_426(MEM_DISTRIB,MSG_SIZE,IDWLOAD,NPROCS)
      ENDIF
      LREF = LOAD_FLOPS(MYID)
      DO i=1, NPROCS
         IF (WLOAD(i).LT.LREF) NLESS=NLESS+1
      ENDDO
      ZMUMPS_186 = NLESS
      RETURN
      END FUNCTION ZMUMPS_186
      SUBROUTINE ZMUMPS_189(MEM_DISTRIB,MSG_SIZE,DEST,
     &     NSLAVES)
      IMPLICIT NONE
      INTEGER NSLAVES
      INTEGER DEST(NSLAVES)
      INTEGER, DIMENSION(0:NPROCS - 1) :: MEM_DISTRIB
      INTEGER i,J,NBDEST
      DOUBLE PRECISION MSG_SIZE
      IF ( NSLAVES.eq.NPROCS-1 ) THEN
        J = MYID+1
        DO i=1,NSLAVES
           J=J+1
           IF (J.GT.NPROCS) J=1
           DEST(i) = J - 1
        ENDDO
      ELSE
        DO i=1,NPROCS
           IDWLOAD(i) = i - 1
        ENDDO
        CALL MUMPS_558(NPROCS, WLOAD, IDWLOAD)
         NBDEST = 0
         DO i=1, NSLAVES
            J = IDWLOAD(i)
            IF (J.NE.MYID) THEN
               NBDEST = NBDEST+1
               DEST(NBDEST) = J
            ENDIF
         ENDDO
         IF (NBDEST.NE.NSLAVES) THEN
            DEST(NSLAVES) = IDWLOAD(NSLAVES+1)
         ENDIF
         IF(BDC_MD)THEN
            J=NSLAVES+1
            do i=NSLAVES+1,NPROCS
               IF(IDWLOAD(i).NE.MYID)THEN
                  DEST(J)= IDWLOAD(i)
                  J=J+1
               ENDIF
            end do
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_189
      SUBROUTINE ZMUMPS_183( INFO1, IERR )
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER, intent(in) :: INFO1
      INTEGER, intent(out) :: IERR
      IERR=0
      DEALLOCATE( LOAD_FLOPS )
      DEALLOCATE( WLOAD )
      DEALLOCATE( IDWLOAD )
#if ! defined(OLD_LOAD_MECHANISM)
      DEALLOCATE(FUTURE_NIV2)
#endif
      IF(BDC_MD)THEN
         DEALLOCATE(MD_MEM)
         DEALLOCATE(LU_USAGE)
         DEALLOCATE(TAB_MAXS)
      ENDIF
      IF ( BDC_MEM ) DEALLOCATE( DM_MEM )
      IF ( BDC_POOL) DEALLOCATE( POOL_MEM )
      IF ( BDC_SBTR) THEN
         DEALLOCATE( SBTR_MEM )
         DEALLOCATE( SBTR_CUR )
         DEALLOCATE(SBTR_FIRST_POS_IN_POOL)
         NULLIFY(MY_FIRST_LEAF)
         NULLIFY(MY_NB_LEAF)
         NULLIFY(MY_ROOT_SBTR)
      ENDIF
      IF(KEEP_LOAD(76).EQ.4)THEN
         NULLIFY(DEPTH_FIRST_LOAD)
      ENDIF
      IF(KEEP_LOAD(76).EQ.5)THEN
         NULLIFY(COST_TRAV)
      ENDIF
      IF((KEEP_LOAD(76).EQ.4).OR.(KEEP_LOAD(76).EQ.6))THEN
         NULLIFY(DEPTH_FIRST_LOAD)
         NULLIFY(DEPTH_FIRST_SEQ_LOAD)
         NULLIFY(SBTR_ID_LOAD)
      ENDIF
      IF (BDC_M2_MEM.OR.BDC_M2_FLOPS) THEN
        DEALLOCATE(NB_SON,POOL_NIV2,POOL_NIV2_COST, NIV2)
      END IF
      IF((KEEP_LOAD(81).EQ.2).OR.(KEEP_LOAD(81).EQ.3))THEN
         DEALLOCATE(CB_COST_MEM)
         DEALLOCATE(CB_COST_ID)
      ENDIF
      NULLIFY(ND_LOAD)
      NULLIFY(KEEP_LOAD)
      NULLIFY(KEEP8_LOAD)
      NULLIFY(FILS_LOAD)
      NULLIFY(FRERE_LOAD)
      NULLIFY(PROCNODE_LOAD)
      NULLIFY(STEP_LOAD)
      NULLIFY(NE_LOAD)
      NULLIFY(CAND_LOAD)
      NULLIFY(STEP_TO_NIV2_LOAD)
      NULLIFY(DAD_LOAD)
      IF (BDC_SBTR.OR.BDC_POOL_MNG) THEN
         DEALLOCATE(MEM_SUBTREE)
         DEALLOCATE(SBTR_PEAK_ARRAY)
         DEALLOCATE(SBTR_CUR_ARRAY)
      ENDIF
      CALL ZMUMPS_58( IERR )
      CALL ZMUMPS_150( MYID, COMM_LD,
     &     BUF_LOAD_RECV, LBUF_LOAD_RECV,
     &     LBUF_LOAD_RECV_BYTES )
      DEALLOCATE(BUF_LOAD_RECV)
      END SUBROUTINE ZMUMPS_183
#if defined (LAMPORT_)
      RECURSIVE SUBROUTINE ZMUMPS_467(COMM, KEEP)
#else
      SUBROUTINE ZMUMPS_467(COMM, KEEP)
#endif
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MSGTAG, MSGLEN, MSGSOU,COMM
      INTEGER KEEP(500)
      INTEGER STATUS(MPI_STATUS_SIZE)
      LOGICAL FLAG
 10   CONTINUE
      CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &     FLAG, STATUS, IERR )
      IF (FLAG) THEN
        KEEP(65)=KEEP(65)+1
        MSGTAG = STATUS( MPI_TAG )
        MSGSOU = STATUS( MPI_SOURCE )
        IF ( MSGTAG .NE. UPDATE_LOAD) THEN
          write(*,*) "Internal error 1 in ZMUMPS_467",
     &    MSGTAG
          CALL MUMPS_ABORT()
        ENDIF
        CALL MPI_GET_COUNT(STATUS, MPI_PACKED, MSGLEN, IERR)
        IF ( MSGLEN > LBUF_LOAD_RECV_BYTES ) THEN
          write(*,*) "Internal error 2 in ZMUMPS_467",
     &    MSGLEN, LBUF_LOAD_RECV_BYTES
          CALL MUMPS_ABORT()
        ENDIF
        CALL MPI_RECV( BUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES,
     &    MPI_PACKED, MSGSOU, MSGTAG, COMM_LD, STATUS, IERR)
        CALL ZMUMPS_187( MSGSOU, BUF_LOAD_RECV,
     &  LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES )
        GOTO 10
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_467
      RECURSIVE SUBROUTINE ZMUMPS_187
     &   ( MSGSOU, BUFR, LBUFR, LBUFR_BYTES )
      IMPLICIT NONE
      INTEGER MSGSOU, LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INCLUDE 'mpif.h'
      INTEGER POSITION, IERR, WHAT, NSLAVES, i
      DOUBLE PRECISION LOAD_RECEIVED
      INTEGER INODE_RECEIVED,NCB_RECEIVED
      DOUBLE PRECISION SURF
      INTEGER, POINTER, DIMENSION (:) :: LIST_SLAVES
      DOUBLE PRECISION, POINTER, DIMENSION (:) :: LOAD_INCR
      EXTERNAL MUMPS_330
      INTEGER MUMPS_330
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     WHAT, 1, MPI_INTEGER,
     &     COMM_LD, IERR )
      IF ( WHAT == 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
#else
#endif
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
      LOAD_FLOPS( MSGSOU ) = LOAD_FLOPS(MSGSOU) + LOAD_RECEIVED
#else
      LOAD_FLOPS( MSGSOU ) = LOAD_RECEIVED
#endif
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &       LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     &       COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
          DM_MEM(MSGSOU)  = DM_MEM(MSGSOU) + LOAD_RECEIVED
#else
          DM_MEM(MSGSOU)  = LOAD_RECEIVED
#endif
          MAX_PEAK_STK=max(MAX_PEAK_STK,DM_MEM(MSGSOU))
        END IF
        IF(BDC_SBTR)THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &       LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     &       COMM_LD, IERR )
          SBTR_CUR(MSGSOU)=LOAD_RECEIVED
        ENDIF
        IF(BDC_MD)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &          LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     &          COMM_LD, IERR )
           IF(KEEP_LOAD(201).EQ.0)THEN
              LU_USAGE(MSGSOU)=LOAD_RECEIVED
           ENDIF
        ENDIF
      ELSEIF (( WHAT == 1).OR.(WHAT.EQ.19)) THEN
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     NSLAVES, 1, MPI_INTEGER,
     &     COMM_LD, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     INODE_RECEIVED, 1, MPI_INTEGER,
     &     COMM_LD, IERR )
        LIST_SLAVES => IDWLOAD
        LOAD_INCR => WLOAD
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LIST_SLAVES(1), NSLAVES, MPI_INTEGER,
     &     COMM_LD, IERR)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR)
#if defined(OLD_LOAD_MECHANISM)
#if defined(CHECK_COHERENCE) 
        WRITE(*,*)MYID,':Receiving M2A from',MSGSOU
        i=1
        DO WHILE ((i.LE.NSLAVES).AND.(LIST_SLAVES(i).NE.MYID))
           i=i+1
        ENDDO
        IF(i.LT.(NSLAVES+1))THEN
           NB_LEVEL2=NB_LEVEL2+1
           WRITE(*,*)MYID,':NB_LEVEL2=',NB_LEVEL2
           AMI_CHOSEN=.TRUE.
           IF(KEEP_LOAD(73).EQ.1)THEN
              IF(.NOT.IS_DISPLAYED)THEN
              WRITE(*,*)MYID,': Begin of Incoherent state (2) at time=',
     &                MPI_WTIME()-TIME_REF
              IS_DISPLAYED=.TRUE.
              ENDIF
            ENDIF
        ENDIF
        IF(KEEP_LOAD(73).EQ.1) GOTO 344
#endif
#endif
        DO i = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
          IF ( LIST_SLAVES(i) /= MYID ) THEN
#endif
            LOAD_FLOPS(LIST_SLAVES(i)) =
     &      LOAD_FLOPS(LIST_SLAVES(i)) + 
     &      LOAD_INCR(i)
#if defined(OLD_LOAD_MECHANISM)
          END IF
#endif
        END DO
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR)
          DO i = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
            IF ( LIST_SLAVES(i) /= MYID ) THEN
#endif
              DM_MEM(LIST_SLAVES(i)) = DM_MEM(LIST_SLAVES(i)) + 
     &        LOAD_INCR(i)
              MAX_PEAK_STK=max(MAX_PEAK_STK,DM_MEM(LIST_SLAVES(i)))
#if defined(OLD_LOAD_MECHANISM)
            END IF
#endif
          END DO
        END IF
        IF(WHAT.EQ.19)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &          LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &          COMM_LD, IERR)
           CALL ZMUMPS_819(INODE_RECEIVED)
           CB_COST_ID(POS_ID)=INODE_RECEIVED
           CB_COST_ID(POS_ID+1)=NSLAVES
           CB_COST_ID(POS_ID+2)=POS_MEM
           POS_ID=POS_ID+3
           DO i=1,NSLAVES
              WRITE(*,*)MYID,':',LIST_SLAVES(i),'->',LOAD_INCR(i)
              CB_COST_MEM(POS_MEM)=int(LIST_SLAVES(i),8)
              POS_MEM=POS_MEM+1
              CB_COST_MEM(POS_MEM)=int(LOAD_INCR(i),8)
              POS_MEM=POS_MEM+1
           ENDDO
        ENDIF
#if defined(OLD_LOAD_MECHANISM)
#if defined(CHECK_COHERENCE) 
 344    CONTINUE
#endif
#endif
        NULLIFY( LIST_SLAVES )
        NULLIFY( LOAD_INCR )
      ELSE IF (WHAT == 2 ) THEN
        IF ( .not. BDC_POOL ) THEN
          WRITE(*,*) "Internal error 2 in ZMUMPS_187"
          CALL MUMPS_ABORT()
        END IF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR )
        POOL_MEM(MSGSOU)=LOAD_RECEIVED
      ELSE IF ( WHAT == 3 ) THEN
        IF ( .NOT. BDC_SBTR) THEN
          WRITE(*,*) "Internal error 3 in ZMUMPS_187"
          CALL MUMPS_ABORT()
        ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR )
        SBTR_MEM(MSGSOU)=SBTR_MEM(MSGSOU)+LOAD_RECEIVED
#if ! defined(OLD_LOAD_MECHANISM)
      ELSE IF (WHAT == 4) THEN
        FUTURE_NIV2(MSGSOU+1)=0
        IF(BDC_MD)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &          SURF, 1, MPI_DOUBLE_PRECISION,
     &          COMM_LD, IERR )
          MD_MEM(MSGSOU)=999999999_8
          TAB_MAXS(MSGSOU)=TAB_MAXS(MSGSOU)+int(SURF,8)
        ENDIF
#endif
        IF(BDC_M2_MEM.OR.BDC_M2_FLOPS)THEN
        ENDIF
      ELSE IF (WHAT == 5) THEN
         IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
            WRITE(*,*) "Internal error 7 in ZMUMPS_187"
            CALL MUMPS_ABORT()
         ENDIF
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     INODE_RECEIVED, 1,
     &     MPI_INTEGER,
     &     COMM_LD, IERR )
         IF(BDC_M2_MEM) THEN
            CALL ZMUMPS_816(INODE_RECEIVED)
         ELSEIF(BDC_M2_FLOPS) THEN
            CALL ZMUMPS_817(INODE_RECEIVED)
         ENDIF
         IF((KEEP_LOAD(81).EQ.2).OR.(KEEP_LOAD(81).EQ.3))THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           INODE_RECEIVED, 1,
     &           MPI_INTEGER,
     &           COMM_LD, IERR )   
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              NCB_RECEIVED, 1,
     &              MPI_INTEGER,
     &              COMM_LD, IERR )   
            IF(
     &          MUMPS_330(PROCNODE_LOAD(STEP_LOAD(INODE_RECEIVED)),
     &                         NPROCS).EQ.1
     &        )THEN
               CB_COST_ID(POS_ID)=INODE_RECEIVED
               CB_COST_ID(POS_ID+1)=1
               CB_COST_ID(POS_ID+2)=POS_MEM
               POS_ID=POS_ID+3
               CB_COST_MEM(POS_MEM)=int(MSGSOU,8)
               POS_MEM=POS_MEM+1
               CB_COST_MEM(POS_MEM)=int(NCB_RECEIVED,8)*
     &              int(NCB_RECEIVED,8)
               POS_MEM=POS_MEM+1
            ENDIF
         ENDIF
      ELSE IF ( WHAT == 6 ) THEN
         IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
            WRITE(*,*) "Internal error 8 in ZMUMPS_187"
            CALL MUMPS_ABORT()
         ENDIF
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR )
         IF(BDC_M2_MEM) THEN
            NIV2(MSGSOU+1) = LOAD_RECEIVED
         ELSEIF(BDC_M2_FLOPS) THEN
            NIV2(MSGSOU+1) = NIV2(MSGSOU+1) + LOAD_RECEIVED
            IF(NIV2(MSGSOU+1).LT.0.0D0)THEN
               IF(abs(NIV2(MSGSOU+1)).LE.
     &              sqrt(epsilon(LOAD_RECEIVED)))THEN
                  NIV2(MSGSOU+1)=0.0D0
               ELSE
                  WRITE(*,*)'problem with NIV2_FLOPS message',
     &                 NIV2(MSGSOU+1),MSGSOU,LOAD_RECEIVED
                  CALL MUMPS_ABORT()
               ENDIF
            ENDIF
         ENDIF
      ELSEIF(WHAT == 17)THEN
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR )
         IF(BDC_M2_MEM) THEN
            NIV2(MSGSOU+1) = LOAD_RECEIVED
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           LOAD_RECEIVED, 1,
     &           MPI_DOUBLE_PRECISION,
     &           COMM_LD, IERR )
            IF(BDC_MD)THEN
#if ! defined(OLD_LOAD_MECHANISM)
               DM_MEM(MYID)=DM_MEM(MYID)+LOAD_RECEIVED
#else
               DM_MEM(MYID)=LOAD_RECEIVED
#endif
            ELSEIF(BDC_POOL)THEN
               POOL_MEM(MSGSOU)=LOAD_RECEIVED
            ENDIF
         ELSEIF(BDC_M2_FLOPS) THEN
            NIV2(MSGSOU+1) = NIV2(MSGSOU+1) + LOAD_RECEIVED            
            IF(NIV2(MSGSOU+1).LT.0.0D0)THEN
               WRITE(*,*)'problem with NIV2_FLOPS message',
     &              NIV2(MSGSOU+1),MSGSOU,LOAD_RECEIVED
               CALL MUMPS_ABORT()
            ENDIF
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           LOAD_RECEIVED, 1,
     &           MPI_DOUBLE_PRECISION,
     &           COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
            LOAD_FLOPS( MSGSOU ) = LOAD_FLOPS(MSGSOU) + LOAD_RECEIVED
#else
            LOAD_FLOPS( MSGSOU ) = LOAD_RECEIVED
#endif
         ENDIF
      ELSEIF ( WHAT == 7 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 4
     &in ZMUMPS_187'
            CALL MUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     NSLAVES, 1, MPI_INTEGER,
     &     COMM_LD, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     INODE_RECEIVED, 1, MPI_INTEGER,
     &     COMM_LD, IERR )
        LIST_SLAVES => IDWLOAD
        LOAD_INCR => WLOAD
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LIST_SLAVES(1), NSLAVES, MPI_INTEGER,
     &     COMM_LD, IERR)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR)
        DO i = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
          IF ( LIST_SLAVES(i) /= MYID ) THEN
#endif
            MD_MEM(LIST_SLAVES(i)) =
     &      MD_MEM(LIST_SLAVES(i)) + 
     &      int(LOAD_INCR(i),8)
#if ! defined(OLD_LOAD_MECHANISM)
            IF(FUTURE_NIV2(LIST_SLAVES(i)+1).EQ.0)THEN
               MD_MEM(LIST_SLAVES(i))=999999999_8
            ENDIF
#endif
#if defined(OLD_LOAD_MECHANISM)
          END IF
#endif
        END DO
      ELSEIF ( WHAT == 8 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 5
     &in ZMUMPS_187'
            CALL MUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR )
        MD_MEM(MSGSOU)=MD_MEM(MSGSOU)+int(LOAD_RECEIVED,8)
#if ! defined(OLD_LOAD_MECHANISM)
        IF(FUTURE_NIV2(MSGSOU+1).EQ.0)THEN
           MD_MEM(MSGSOU)=999999999_8
        ENDIF
#endif
      ELSEIF ( WHAT == 9 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 6
     &in ZMUMPS_187'
            CALL MUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR )
        TAB_MAXS(MSGSOU)=int(LOAD_RECEIVED,8)
      ELSE
          WRITE(*,*) "Internal error 1 in ZMUMPS_187"
          CALL MUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_187
      integer function ZMUMPS_409
     &                 (MEM_DISTRIB,CAND,
     &                  K69,
     &                  SLAVEF,MSG_SIZE,
     &                  NMB_OF_CAND )
      implicit none      
      integer, intent(in) :: K69, SLAVEF
      INTEGER, intent(in) :: CAND(SLAVEF+1)
      INTEGER, DIMENSION(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      INTEGER, intent(out) :: NMB_OF_CAND
      integer i,nless
      DOUBLE PRECISION lref
      DOUBLE PRECISION MSG_SIZE
      nless = 0 
      NMB_OF_CAND=CAND(SLAVEF+1)
      do i=1,NMB_OF_CAND
         WLOAD(i)=LOAD_FLOPS(CAND(i))
         IF(BDC_M2_FLOPS)THEN
            WLOAD(i)=WLOAD(i)+NIV2(CAND(i)+1)
         ENDIF
      end do
      IF(K69 .gt. 1) THEN
         CALL ZMUMPS_426(MEM_DISTRIB,MSG_SIZE,
     &        CAND,NMB_OF_CAND)
      ENDIF
      lref = LOAD_FLOPS(MYID)
      do i=1, NMB_OF_CAND
         if (WLOAD(i).lt.lref) nless=nless+1
      end do 
      ZMUMPS_409 = nless
      return
      end function ZMUMPS_409
      subroutine ZMUMPS_384
     &           (MEM_DISTRIB,CAND, 
     &
     &            SLAVEF,
     &            nslaves_inode, DEST)
      implicit none
      integer, intent(in) :: nslaves_inode, SLAVEF
      integer, intent(in) :: CAND(SLAVEF+1)
      integer, dimension(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      integer, intent(out) :: DEST(CAND(SLAVEF+1))
      integer i,j,NMB_OF_CAND
      external MUMPS_558
      NMB_OF_CAND = CAND(SLAVEF+1)
      if(nslaves_inode.ge.NPROCS .or.
     &   nslaves_inode.gt.NMB_OF_CAND) then
         write(*,*)'Internal error in ZMUMPS_384',
     &   nslaves_inode, NPROCS, NMB_OF_CAND
         CALL MUMPS_ABORT()
      end if
      if (nslaves_inode.eq.NPROCS-1) then
         j=MYID+1
         do i=1,nslaves_inode
            if(j.ge.NPROCS) j=0
            DEST(i)=j
            j=j+1
         end do
      else
        do i=1,NMB_OF_CAND
               IDWLOAD(i)=i
        end do
        call MUMPS_558(NMB_OF_CAND,
     &       WLOAD(1),IDWLOAD(1) )
        do i=1,nslaves_inode
           DEST(i)= CAND(IDWLOAD(i))
        end do
        IF(BDC_MD)THEN
           do i=nslaves_inode+1,NMB_OF_CAND
              DEST(i)= CAND(IDWLOAD(i))
           end do
        ENDIF
      end if   
      return
      end subroutine ZMUMPS_384
      SUBROUTINE ZMUMPS_425(K69)
      IMPLICIT NONE
      INTEGER K69
      IF (K69 .LE. 4) THEN
         ALPHA = 0.0d0
         BETA = 0.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 5) THEN
         ALPHA = 0.5d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 6) THEN
         ALPHA = 0.5d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 7) THEN
         ALPHA = 0.5d0
         BETA = 150000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 8) THEN
         ALPHA = 1.0d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 9) THEN
         ALPHA = 1.0d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 10) THEN
         ALPHA = 1.0d0
         BETA = 150000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 11) THEN
         ALPHA = 1.5d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 12) THEN
         ALPHA = 1.5d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      ALPHA = 1.5d0
      BETA = 150000.0d0
      RETURN
      END SUBROUTINE ZMUMPS_425
      SUBROUTINE ZMUMPS_426(MEM_DISTRIB,MSG_SIZE,ARRAY_ADM,LEN)
      IMPLICIT NONE
      INTEGER i,LEN
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      DOUBLE PRECISION MSG_SIZE,FORBIGMSG
      INTEGER ARRAY_ADM(LEN)
      DOUBLE PRECISION MY_LOAD
      FORBIGMSG = 1.0d0
      IF (K69 .lt.2) THEN
         RETURN
      ENDIF
      IF(BDC_M2_FLOPS)THEN
         MY_LOAD=LOAD_FLOPS(MYID)+NIV2(MYID+1)
      ELSE
         MY_LOAD=LOAD_FLOPS(MYID)
      ENDIF
      IF((MSG_SIZE * dble(K35) ) .gt. 3200000.0d0) THEN
         FORBIGMSG = 2.0d0
      ENDIF
      IF (K69 .le. 4) THEN
         DO i = 1,LEN
            IF ((MEM_DISTRIB(ARRAY_ADM(i)) .EQ. 1) .AND.
     &      WLOAD(i) .LT. MY_LOAD ) THEN
               WLOAD(i) = WLOAD(i)/MY_LOAD
            ELSE
              IF ( MEM_DISTRIB(ARRAY_ADM(i)) .NE. 1 ) THEN
                WLOAD(i) = WLOAD(i) *
     &              dble(MEM_DISTRIB(ARRAY_ADM(i)))
     &              * FORBIGMSG
     &              + dble(2)
              ENDIF
            ENDIF
         ENDDO
         RETURN
      ENDIF
      DO i = 1,LEN
         IF ((MEM_DISTRIB(ARRAY_ADM(i)) .EQ. 1) .AND.
     &        WLOAD(i) .LT. MY_LOAD ) THEN
            WLOAD(i) = WLOAD(i) /  MY_LOAD
         ELSE
            IF(MEM_DISTRIB(ARRAY_ADM(i)) .NE. 1) THEN     
               WLOAD(i) = (WLOAD(i) +
     &              ALPHA * MSG_SIZE * dble(K35)  +
     &              BETA) * FORBIGMSG
            ENDIF
         ENDIF
      ENDDO  
      RETURN
      END SUBROUTINE ZMUMPS_426
      SUBROUTINE ZMUMPS_461(MYID, SLAVEF, COMM,
     &     TAB_POS, NASS, KEEP,KEEP8, LIST_SLAVES, NSLAVES,INODE)
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: MYID, SLAVEF, COMM, NASS, NSLAVES
      INTEGER, INTENT (IN) :: TAB_POS(SLAVEF+2)
      INTEGER, INTENT (IN) :: LIST_SLAVES( NSLAVES )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NCB, NFRONT, NBROWS_SLAVE
      INTEGER i, IERR,WHAT,INODE
      DOUBLE PRECISION MEM_INCREMENT ( NSLAVES )
      DOUBLE PRECISION FLOPS_INCREMENT( NSLAVES )
      DOUBLE PRECISION CB_BAND( NSLAVES )
      IF((KEEP(81).NE.2).AND.(KEEP(81).NE.3))THEN
         WHAT=1
      ELSE
         WHAT=19
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      FUTURE_NIV2(MYID+1) = FUTURE_NIV2(MYID+1) - 1
      IF ( FUTURE_NIV2(MYID+1) < 0 ) THEN
        WRITE(*,*) "Internal error in ZMUMPS_461"
        CALL MUMPS_ABORT()
      ENDIF
      IF ( FUTURE_NIV2(MYID + 1) == 0 ) THEN
 112    CONTINUE
        CALL ZMUMPS_502(COMM,MYID,SLAVEF,
     &       dble(MAX_SURF_MASTER),IERR)
        IF (IERR == -1 ) THEN
          CALL ZMUMPS_467(COMM_LD, KEEP)
          GOTO 112
        ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in ZMUMPS_461",
     &    IERR
          CALL MUMPS_ABORT()
        ENDIF
      TAB_MAXS(MYID) = TAB_MAXS(MYID) + int(MAX_SURF_MASTER,8)
      ENDIF
#endif
      IF ( NSLAVES /= TAB_POS(SLAVEF + 2) ) THEN
        write(*,*) "Error 1 in ZMUMPS_461",
     &             NSLAVES, TAB_POS(SLAVEF+2)
        CALL MUMPS_ABORT()
      ENDIF
      NCB = TAB_POS(NSLAVES+1) - 1
      NFRONT = NCB + NASS
      DO i = 1, NSLAVES
         NBROWS_SLAVE = TAB_POS(i+1) - TAB_POS(i)
         IF ( KEEP(50) == 0 ) THEN
            FLOPS_INCREMENT( i ) = (dble(NBROWS_SLAVE)*dble( NASS ))+
     &           dble(NBROWS_SLAVE) * dble(NASS) *
     &           dble(2*NFRONT-NASS-1)
         ELSE
            FLOPS_INCREMENT( i ) = dble(NBROWS_SLAVE) * dble(NASS ) *
     &           dble( 2 * ( NASS + TAB_POS(i+1) - 1 ) 
     &           - NBROWS_SLAVE - NASS + 1 )
         ENDIF
         IF ( BDC_MEM ) THEN
            IF ( KEEP(50) == 0 ) THEN
               MEM_INCREMENT( i ) = dble(NBROWS_SLAVE) *
     &              dble(NFRONT)
            ELSE
               MEM_INCREMENT( i ) = dble(NBROWS_SLAVE) *
     &              dble( NASS + TAB_POS(i+1) - 1 )
            END IF
         ENDIF
         IF((KEEP(81).NE.2).AND.(KEEP(81).NE.3))THEN
            CB_BAND(i)=dble(-999999)
         ELSE
            IF ( KEEP(50) == 0 ) THEN
               CB_BAND( i ) = dble(NBROWS_SLAVE) *
     &              dble(NFRONT-NASS)
            ELSE
               CB_BAND( i ) = dble(NBROWS_SLAVE) *
     &              dble(TAB_POS(i+1)-1)
            END IF
         ENDIF
      END DO
      IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
         CB_COST_ID(POS_ID)=INODE
         CB_COST_ID(POS_ID+1)=NSLAVES
         CB_COST_ID(POS_ID+2)=POS_MEM
         POS_ID=POS_ID+3
         DO i=1,NSLAVES
            CB_COST_MEM(POS_MEM)=int(LIST_SLAVES(i),8)
            POS_MEM=POS_MEM+1
            CB_COST_MEM(POS_MEM)=int(CB_BAND(i),8)
            POS_MEM=POS_MEM+1
         ENDDO
      ENDIF
 111  CONTINUE
      CALL ZMUMPS_524(BDC_MEM, COMM, MYID, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     &     FUTURE_NIV2,
#endif
     &     NSLAVES, LIST_SLAVES,INODE,
     &     MEM_INCREMENT,
     &     FLOPS_INCREMENT,CB_BAND, WHAT,IERR)
        IF ( IERR == -1 ) THEN
          CALL ZMUMPS_467(COMM_LD, KEEP)
          GOTO 111
        ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in ZMUMPS_461",
     &    IERR
          CALL MUMPS_ABORT()
        ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (FUTURE_NIV2(MYID+1) .NE. 0) THEN
#endif
        DO i = 1, NSLAVES
          LOAD_FLOPS(LIST_SLAVES(i)) = LOAD_FLOPS(LIST_SLAVES(i))
     &       +  FLOPS_INCREMENT(i)
          IF ( BDC_MEM ) THEN
            DM_MEM(LIST_SLAVES(i)) = DM_MEM(LIST_SLAVES(i))
     &       +  MEM_INCREMENT(i)
          END IF
        ENDDO
#if ! defined(OLD_LOAD_MECHANISM)
      ENDIF
#endif
      RETURN
      END SUBROUTINE ZMUMPS_461
      SUBROUTINE ZMUMPS_500(
     &     POOL, LPOOL,
     &     PROCNODE, KEEP,KEEP8, SLAVEF, COMM, MYID, STEP, N,
     &     ND, FILS )
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER LPOOL, SLAVEF, COMM, MYID 
      INTEGER N, KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER POOL( LPOOL ), PROCNODE( KEEP(28) ), STEP( N )
      INTEGER ND( KEEP(28) ), FILS( N )
      INTEGER i, INODE, NELIM, NFR, LEVEL, IERR, WHAT
      DOUBLE PRECISION COST
      INTEGER NBINSUBTREE,NBTOP,INSUBTREE
      INTEGER MUMPS_330
      EXTERNAL MUMPS_330
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF(BDC_MD)THEN
         RETURN
      ENDIF
      IF((KEEP(76).EQ.0).OR.(KEEP(76).EQ.2))THEN
         IF(NBTOP.NE.0)THEN
            DO i = LPOOL-NBTOP-2, min(LPOOL-3,LPOOL-NBTOP-2+3)
               INODE = POOL( i )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=dble(0) 
            GOTO 30
         ELSE
            DO i = NBINSUBTREE, max(1,NBINSUBTREE-3), -1
               INODE = POOL( i )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=dble(0) 
            GOTO 30
         ENDIF
      ELSE
         IF(KEEP(76).EQ.1)THEN
            IF(INSUBTREE.EQ.1)THEN
               DO i = NBINSUBTREE, max(1,NBINSUBTREE-3), -1
                  INODE = POOL( i )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=dble(0) 
               GOTO 30
            ELSE
               DO i = LPOOL-NBTOP-2, min(LPOOL-3,LPOOL-NBTOP-2+3)
                  INODE = POOL( i )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=dble(0) 
               GOTO 30
            ENDIF
         ELSE
            WRITE(*,*)
     &      'Internal error: Unknown pool management strategy'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
 20   CONTINUE
        i = INODE
        NELIM = 0
 10     CONTINUE
        IF ( i > 0 ) THEN
          NELIM = NELIM + 1
          i = FILS(i)
          GOTO 10
        ENDIF
        NFR = ND( STEP(INODE) )
        LEVEL = MUMPS_330( PROCNODE(STEP(INODE)), SLAVEF )
        IF (LEVEL .EQ. 1) THEN
          COST = dble( NFR ) * dble( NFR )
        ELSE
          IF ( KEEP(50) == 0 ) THEN
            COST = dble( NFR ) * dble( NELIM )
          ELSE
            COST = dble( NELIM ) * dble( NELIM )
          ENDIF
        ENDIF
 30   CONTINUE
      IF ( abs(POOL_LAST_COST_SENT-COST).GT.DM_THRES_MEM ) THEN
        WHAT = 2
 111    CONTINUE
        CALL ZMUMPS_460( WHAT,
     &         COMM, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     &               FUTURE_NIV2,
#endif
     &         COST, dble(0),MYID, IERR  )
        POOL_LAST_COST_SENT = COST
        POOL_MEM(MYID)=COST
        IF ( IERR == -1 )THEN
          CALL ZMUMPS_467(COMM_LD, KEEP)
          GOTO 111
        ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in ZMUMPS_500",
     &    IERR
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_500
      SUBROUTINE ZMUMPS_501(
     &     OK,INODE,POOL,LPOOL,MYID,SLAVEF,COMM,KEEP,KEEP8)
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER LPOOL,MYID,SLAVEF,COMM,INODE
      INTEGER POOL(LPOOL),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER WHAT,IERR
      LOGICAL OK
      DOUBLE PRECISION COST
      LOGICAL FLAG
      EXTERNAL MUMPS_283,MUMPS_170
      LOGICAL MUMPS_283,MUMPS_170
      IF((INODE.LE.0).OR.(INODE.GT.N_LOAD)) THEN
         RETURN
      ENDIF
      IF (.NOT.MUMPS_170(
     &     PROCNODE_LOAD(STEP_LOAD(INODE)), NPROCS)
     &   ) THEN
         RETURN         
      ENDIF
      IF(MUMPS_283(PROCNODE_LOAD(STEP_LOAD(INODE)),NPROCS))THEN
         IF(NE_LOAD(STEP_LOAD(INODE)).EQ.0)THEN
            RETURN
         ENDIF
      ENDIF
      FLAG=.FALSE.
      IF(INDICE_SBTR.LE.NB_SUBTREES)THEN
         IF(INODE.EQ.MY_FIRST_LEAF(INDICE_SBTR))THEN
            FLAG=.TRUE.
         ENDIF
      ENDIF
      IF(FLAG)THEN
         SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY)=MEM_SUBTREE(INDICE_SBTR)
         SBTR_CUR_ARRAY(INDICE_SBTR_ARRAY)=SBTR_CUR(MYID)
         INDICE_SBTR_ARRAY=INDICE_SBTR_ARRAY+1
         WHAT = 3
         IF(dble(MEM_SUBTREE(INDICE_SBTR)).GE.DM_THRES_MEM)THEN
 111        CONTINUE
            CALL ZMUMPS_460(
     &           WHAT, COMM, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     &           FUTURE_NIV2,
#endif
     &           dble(MEM_SUBTREE(INDICE_SBTR)), dble(0),MYID, IERR  )
            IF ( IERR == -1 )THEN
               CALL ZMUMPS_467(COMM_LD, KEEP)
               GOTO 111
            ELSE IF ( IERR .NE. 0 ) THEN
               WRITE(*,*)
     &         "Internal Error 1 in ZMUMPS_501",
     &         IERR
               CALL MUMPS_ABORT()
            ENDIF
         ENDIF
         SBTR_MEM(MYID)=SBTR_MEM(MYID)+
     &        dble(MEM_SUBTREE(INDICE_SBTR))
         INDICE_SBTR=INDICE_SBTR+1
         IF(INSIDE_SUBTREE.EQ.0)THEN
            INSIDE_SUBTREE=1
         ENDIF
      ELSE 
         IF(INODE.EQ.MY_ROOT_SBTR(INDICE_SBTR-1))THEN
            WHAT = 3
            COST=-SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY-1)
            IF(abs(COST).GE.DM_THRES_MEM)THEN
 112           CONTINUE
               CALL ZMUMPS_460(
     &              WHAT, COMM, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     &              FUTURE_NIV2,
#endif
     &              COST, dble(0) ,MYID,IERR  )
               IF ( IERR == -1 )THEN
                  CALL ZMUMPS_467(COMM_LD, KEEP)
                  GOTO 112
               ELSE IF ( IERR .NE. 0 ) THEN
                  WRITE(*,*)
     &        "Internal Error 3 in ZMUMPS_501",
     &        IERR
                  CALL MUMPS_ABORT()
               ENDIF
            ENDIF
            INDICE_SBTR_ARRAY=INDICE_SBTR_ARRAY-1            
            SBTR_MEM(MYID)=SBTR_MEM(MYID)-
     &           SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY)
            SBTR_CUR(MYID)=SBTR_CUR_ARRAY(INDICE_SBTR_ARRAY)
            IF(INDICE_SBTR_ARRAY.EQ.1)THEN
               SBTR_CUR(MYID)=dble(0)
               INSIDE_SUBTREE=0
            ENDIF
         ENDIF
         ENDIF 
         CONTINUE
      END SUBROUTINE ZMUMPS_501
      SUBROUTINE ZMUMPS_504
     &      (SLAVEF,KEEP,KEEP8,PROCS,MEM_DISTRIB,NCB,NFRONT,
     &       NSLAVES_NODE,TAB_POS,
     &       SLAVES_LIST,SIZE_SLAVES_LIST,MYID)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER NUMBER_OF_PROCS,K47, K48, K50
      INTEGER(8) :: K821
      DOUBLE PRECISION DK821
      INTEGER J
      INTEGER KMIN, KMAX
      INTEGER OTHERS,CHOSEN,SMALL_SET,ACC
      DOUBLE PRECISION SOMME,TMP_SUM
      INTEGER AFFECTED
      INTEGER ADDITIONNAL_ROWS,i,X,REF,POS
      INTEGER(8)::TOTAL_MEM
      LOGICAL FORCE_CAND
      DOUBLE PRECISION TEMP(SLAVEF),PEAK
      INTEGER TEMP_ID(SLAVEF),NB_ROWS(SLAVEF)
      EXTERNAL MPI_WTIME
      DOUBLE PRECISION MPI_WTIME
      IF (KEEP8(21) .GT. 0_8) THEN
      write(*,*)MYID,
     & ": Internal Error 1 in ZMUMPS_504"
      CALL MUMPS_ABORT()
      ENDIF
      K821=abs(KEEP8(21))
      DK821=dble(K821)
      K50=KEEP(50)
      K48=KEEP(48)
      K47=KEEP(47)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
      END IF
      IF(K48.NE.4)THEN
         WRITE(*,*)'ZMUMPS_COMPUTE_PARTI_ACTV_MEM_K821
     &      should be called with KEEP(48) different from 4'
         CALL MUMPS_ABORT()
      ENDIF
         KMIN=1
         KMAX=int(K821/int(NFRONT,8))
         IF(FORCE_CAND)THEN
            DO i=1,PROCS(SLAVEF+1)
               WLOAD(i)=DM_MEM(PROCS(i))
               IDWLOAD(i)=PROCS(i)
            ENDDO
            NUMBER_OF_PROCS=PROCS(SLAVEF+1)
            OTHERS=NUMBER_OF_PROCS
         ELSE
            NUMBER_OF_PROCS=SLAVEF
            WLOAD(1:SLAVEF) = DM_MEM(0:NUMBER_OF_PROCS-1)
            DO i=1,NUMBER_OF_PROCS
               IDWLOAD(i) = i - 1
            ENDDO
            OTHERS=NUMBER_OF_PROCS-1
         ENDIF
         NB_ROWS=0
         CALL MUMPS_558(NUMBER_OF_PROCS, WLOAD, IDWLOAD)
         TOTAL_MEM=int(NCB,8)*int(NFRONT,8)
         SOMME=dble(0)
         J=1
         PEAK=dble(0)
         DO i=1,NUMBER_OF_PROCS
            IF((IDWLOAD(i).NE.MYID))THEN
               PEAK=max(PEAK,WLOAD(i))
               TEMP_ID(J)=IDWLOAD(i)
               TEMP(J)=WLOAD(i)
                IF(BDC_SBTR)THEN
                   TEMP(J)=TEMP(J)+SBTR_MEM(IDWLOAD(i))-
     &                  SBTR_CUR(IDWLOAD(i))
                ENDIF
                IF(BDC_POOL)THEN
                   TEMP(J)=TEMP(J)+POOL_MEM(TEMP_ID(J))
                ENDIF
                IF(BDC_M2_MEM)THEN
                   TEMP(J)=TEMP(J)+NIV2(TEMP_ID(J)+1)
                ENDIF
                J=J+1
            ENDIF
         ENDDO
         NUMBER_OF_PROCS=J-1
         CALL MUMPS_558(NUMBER_OF_PROCS, TEMP, TEMP_ID)
         IF(K50.EQ.0)THEN
           PEAK=max(PEAK,
     &       DM_MEM(MYID)+dble(NFRONT)*dble(NFRONT-NCB))
         ELSE
           PEAK=max(PEAK,
     &       DM_MEM(MYID)+dble(NFRONT-NCB)*dble(NFRONT-NCB))
         ENDIF
         PEAK=max(PEAK,TEMP(OTHERS))
         SOMME=dble(0)
         DO i=1,NUMBER_OF_PROCS
           SOMME=SOMME+TEMP(OTHERS)-TEMP(i)
         ENDDO               
         IF(SOMME.LE.dble(TOTAL_MEM)) THEN
            GOTO 096
         ENDIF
 096     CONTINUE
         SOMME=dble(0)
         DO i=1,OTHERS
            SOMME=SOMME+TEMP(OTHERS)-TEMP(i)
         ENDDO
         IF(dble(TOTAL_MEM).GE.SOMME) THEN
#if defined (OLD_PART)
 887        CONTINUE
#endif
            AFFECTED=0
            CHOSEN=0
            ACC=0
            DO i=1,OTHERS
               IF(K50.EQ.0)THEN
                  IF((TEMP(OTHERS)-TEMP(i)).GT.DK821)THEN
                     TMP_SUM=DK821
                  ELSE
                     TMP_SUM=TEMP(OTHERS)-TEMP(i)
                  ENDIF
                  X=int(TMP_SUM/dble(NFRONT))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
               ENDIF
               IF(K50.NE.0)THEN
                  IF((TEMP(OTHERS)-TEMP(i)).GT.DK821)THEN
                     TMP_SUM=DK821
                  ELSE
                     TMP_SUM=TEMP(OTHERS)-TEMP(i)
                  ENDIF
                  X=int((-dble(NFRONT-NCB+ACC)
     &                 +sqrt(((dble(NFRONT-NCB+ACC)*
     &                 dble(NFRONT-NCB+ACC))+dble(4)*
     &                 (TMP_SUM))))/
     &                 dble(2))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
                  IF(X.LE.0) THEN
                     WRITE(*,*)"Internal Error 2 in
     &                    ZMUMPS_504"
                     CALL MUMPS_ABORT()
                  ENDIF
               ENDIF
               NB_ROWS(i)=X
               CHOSEN=CHOSEN+1
               ACC=ACC+X
               IF(NCB-ACC.LT.KMIN) GOTO 111
               IF(NCB.EQ.ACC) GOTO 111
               ENDDO
 111           CONTINUE
               IF((ACC.GT.NCB))THEN
                  X=0
                  DO i=1,OTHERS
                     X=X+NB_ROWS(i)
                  ENDDO
                  WRITE(*,*)'NCB=',NCB,',SOMME=',X
                  WRITE(*,*)MYID,
     &               ": Internal Error 3 in ZMUMPS_504"
                  CALL MUMPS_ABORT()
               ENDIF
               IF((NCB.NE.ACC))THEN
                  IF(K50.NE.0)THEN
                     IF(CHOSEN.NE.0)THEN
                        ADDITIONNAL_ROWS=NCB-ACC
                        NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+ADDITIONNAL_ROWS
                     ELSE
                        TMP_SUM=dble(TOTAL_MEM)/dble(NUMBER_OF_PROCS)
                        CHOSEN=0
                        ACC=0
                        DO i=1,OTHERS
                           X=int((-dble(NFRONT-NCB+ACC)
     &                          +sqrt(((dble(NFRONT-NCB+ACC)*
     &                          dble(NFRONT-NCB+ACC))+dble(4)*
     &                          (TMP_SUM))))/
     &                          dble(2))
                           IF((ACC+X).GT.NCB) X=NCB-ACC
                           NB_ROWS(i)=X
                           CHOSEN=CHOSEN+1
                           ACC=ACC+X
                           IF(NCB-ACC.LT.KMIN) GOTO 002
                           IF(NCB.EQ.ACC) GOTO 002
                        ENDDO
 002                    CONTINUE
                        IF(ACC.LT.NCB)THEN
                           NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+(NCB-ACC)
                        ENDIF
                     ENDIF
                     GOTO 333
                  ENDIF
                  ADDITIONNAL_ROWS=NCB-ACC
                  DO i=CHOSEN,1,-1
                     IF(int(dble(ADDITIONNAL_ROWS)/
     &                    dble(i)).NE.0)THEN
                        GOTO 222 
                     ENDIF
                  ENDDO
 222              CONTINUE
                  X=int(dble(ADDITIONNAL_ROWS)/dble(i))
                  DO J=1,i
                     NB_ROWS(J)=NB_ROWS(J)+X
                     ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                  ENDDO 
                  IF(ADDITIONNAL_ROWS.NE.0) THEN
                     NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                  ENDIF    
               ENDIF
 333           CONTINUE
               IF(NB_ROWS(CHOSEN).EQ.0) CHOSEN=CHOSEN-1
               GOTO 889
            ELSE
               DO i=OTHERS,1,-1
                  SOMME=dble(0)
                  DO J=1,i
                     SOMME=SOMME+TEMP(J)
                  ENDDO
                  SOMME=(dble(i)*TEMP(i))-SOMME
                  IF(dble(TOTAL_MEM).GE.SOMME) GOTO 444
               ENDDO
 444           CONTINUE
               REF=i
               DO J=1,i
                  IF(TEMP(J).EQ.TEMP(i)) THEN
                     SMALL_SET=J
                     GOTO 123
                  ENDIF
               ENDDO
 123           CONTINUE
               IF(i.EQ.1)THEN 
                  NB_ROWS(i)=NCB
                  CHOSEN=1
                  GOTO 666
               ENDIF
 323           CONTINUE
               AFFECTED=0
               CHOSEN=0
               ACC=0
               DO i=1,SMALL_SET
                  IF(K50.EQ.0)THEN
                     IF((TEMP(SMALL_SET)-TEMP(i)).GT.DK821)THEN
                        TMP_SUM=DK821
                     ELSE
                        TMP_SUM=TEMP(SMALL_SET)-TEMP(i)
                     ENDIF
                     X=int(TMP_SUM/dble(NFRONT))
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  IF(K50.NE.0)THEN
                     IF((TEMP(SMALL_SET)-TEMP(i)).GT.DK821)THEN
                        TMP_SUM=DK821
                     ELSE
                        TMP_SUM=TEMP(SMALL_SET)-TEMP(i)
                     ENDIF
                      X=int((-dble(NFRONT-NCB+ACC)
     &                  +sqrt(((dble(NFRONT-NCB+ACC)*
     &                  dble(NFRONT-NCB+ACC))+dble(4)*
     &                  (TMP_SUM))))/
     &                  dble(2))
                     IF(X.LT.0)THEN
                        WRITE(*,*)MYID,
     &             ': Internal error 4 in ZMUMPS_504'
                        CALL MUMPS_ABORT()
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  NB_ROWS(i)=X
                  ACC=ACC+X
                  CHOSEN=CHOSEN+1
                  IF(NCB-ACC.LT.KMIN) GOTO 888
                  IF(NCB.EQ.ACC) GOTO 888
                  IF(ACC.GT.NCB) THEN
                    WRITE(*,*)MYID,
     &            ': Internal error 5 in ZMUMPS_504'
                    CALL MUMPS_ABORT()
                  ENDIF
               ENDDO
 888           CONTINUE
               SOMME=dble(0)
               X=NFRONT-NCB
               IF((ACC.GT.NCB))THEN
                  WRITE(*,*)MYID,
     &           ':Internal error 6 in ZMUMPS_504'
                  CALL MUMPS_ABORT()
               ENDIF
               IF((ACC.LT.NCB))THEN
                  IF(K50.NE.0)THEN
                     IF(SMALL_SET.LT.OTHERS)THEN
                       SMALL_SET=REF+1
                       REF=SMALL_SET
                       GOTO 323
                     ELSE
                       NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+NCB-ACC
                       GOTO 666
                     ENDIF
                 ENDIF
                 ADDITIONNAL_ROWS=NCB-ACC
#if ! defined (OLD_PART)
                 i=CHOSEN+1
                 DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                .AND.(i.LE.NUMBER_OF_PROCS))
                    J=1
#if ! defined (PART1_)
                    X=int(ADDITIONNAL_ROWS/(i-1))
                    IF((X.EQ.0).AND.(ADDITIONNAL_ROWS.NE.0))THEN
                       DO WHILE ((J.LT.i).AND.(ADDITIONNAL_ROWS.GT.0))
                         NB_ROWS(J)=NB_ROWS(J)+1
                         ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                         J=J+1
                       ENDDO
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          WRITE(*,*)MYID,
     &             ':Internal error 7 in ZMUMPS_504'
                         CALL MUMPS_ABORT()
                      ENDIF
                      GOTO 047
                    ENDIF
                    IF((TEMP(1)+dble((NB_ROWS(1)+X)*NFRONT)).LE.
     &                   TEMP(i))THEN
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(J.LT.i))
                          AFFECTED=X
                          IF((AFFECTED+NB_ROWS(J)).GT.
     &                         KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(J)
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                         AFFECTED
                          J=J+1
                       ENDDO
                    ELSE
#endif
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(J.LE.i))
                          AFFECTED=int((TEMP(i)-(TEMP(J)+
     &                         (dble(NB_ROWS(J))*dble(NFRONT))))
     &                         /dble(NFRONT))
                          IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(J)
                          ENDIF
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          J=J+1
                       ENDDO
#if ! defined (PART1_)
                    ENDIF
#endif
                    i=i+1
                 ENDDO
 047             CONTINUE
                 IF((ADDITIONNAL_ROWS.EQ.0).AND.
     &                (i.LT.NUMBER_OF_PROCS))THEN
                    CHOSEN=i-1
                 ELSE
                    CHOSEN=i-2
                 ENDIF
#if ! defined (PART1_)
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     &                 (ADDITIONNAL_ROWS.NE.0))THEN
                    DO i=1,CHOSEN
                       NB_ROWS(i)=NB_ROWS(i)+1
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 048
                    ENDDO
 048                CONTINUE
                 ENDIF
#endif
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     &                (ADDITIONNAL_ROWS.NE.0))THEN
                    i=CHOSEN+1
                    DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                   .AND.(i.LE.NUMBER_OF_PROCS))
                       J=1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(J.LE.i))
                          AFFECTED=int((TEMP(i)-(TEMP(J)+
     &                         (dble(NB_ROWS(J))*
     &                         dble(NFRONT))))/dble(NFRONT))
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          J=J+1
                       ENDDO
                       i=i+1
                    ENDDO
                    CHOSEN=i-2
                 ENDIF
                 CONTINUE
#else
                 DO i=CHOSEN,1,-1
                    IF(int(dble(ADDITIONNAL_ROWS)/
     &                   dble(i)).NE.0)THEN
                       GOTO 555
                    ENDIF
                 ENDDO
 555             CONTINUE
                 X=int(dble(ADDITIONNAL_ROWS)/dble(i))
                 DO J=1,i
                    IF(NB_ROWS(J)+X.GT.K821/NCB)THEN
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                      ((K821/NCB)-NB_ROWS(J))
                       NB_ROWS(J)=(K821/NFRONT)
                    ELSE IF ((TEMP(J)+dble((NB_ROWS(J)+X))*
     &                      dble(NFRONT)).GT.
     &                 PEAK)THEN
                       AFFECTED=int((PEAK-(TEMP(J)+dble(NB_ROWS(J))*
     &                 dble(NFRONT))/dble(NFRONT))
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                       NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                    ELSE
                       NB_ROWS(J)=NB_ROWS(J)+X
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                    ENDIF
                    IF(((TEMP(J)+dble(NB_ROWS(J))*dble(NFRONT))
     &               .GT. PEAK)
     &               .AND.(SMALL_SET.LT.OTHERS))THEN
                       WRITE(*,*)MYID,
     &           ':Internal error 8 in ZMUMPS_504'
                       SMALL_SET=SMALL_SET+1
                       CALL MUMPS_ABORT()
                    ENDIF
                 ENDDO
                 SOMME=dble(0)
                 DO J=1,CHOSEN
                   SOMME=SOMME+NB_ROWS(J)
                 ENDDO
                 IF(ADDITIONNAL_ROWS.NE.0) THEN
                       DO J=1,CHOSEN
                          IF(NB_ROWS(J).LT.0)THEN
                            WRITE(*,*)MYID,
     &           ':Internal error 9 in ZMUMPS_504'
                            CALL MUMPS_ABORT()
                          ENDIF
                          IF ((TEMP(J)+dble(NB_ROWS(J))
     &                 *dble(NFRONT)).GT.
     &                 PEAK)THEN
                             WRITE(*,*)MYID,
     &         ':Internal error 10 in ZMUMPS_504'
                          CALL MUMPS_ABORT()
                       ENDIF
                          IF ((TEMP(J)+dble(NB_ROWS(J)+
     &                 ADDITIONNAL_ROWS)*dble(NFRONT)).GE.
     &                 PEAK)THEN
                            AFFECTED=int((PEAK-(TEMP(J)+
     &                      dble(NB_ROWS(J))*
     &                      dble(NFRONT))/dble(NFRONT))
                            ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                            NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                            IF((TEMP(J)+dble(NFRONT)*
     &                           dble(NB_ROWS(J))).GT.
     &                         PEAK)THEN
                               WRITE(*,*)MYID,
     &         ':Internal error 11 in ZMUMPS_504'
                              CALL MUMPS_ABORT()
                            ENDIF  
                          ELSE
                            NB_ROWS(J)=NB_ROWS(J)+ADDITIONNAL_ROWS
                            ADDITIONNAL_ROWS=0
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0) GOTO 666
                       ENDDO
                       IF(SMALL_SET.EQ.(NUMBER_OF_PROCS)) THEN
                         NB_ROWS=0
                         GOTO 887
                       ENDIF
                       i=CHOSEN+1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                        .AND.(i.LE.NUMBER_OF_PROCS))
                          IF(NB_ROWS(i)+ADDITIONNAL_ROWS.LT.K821/NFRONT)
     &                         THEN
                             NB_ROWS(i)=NB_ROWS(i)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ELSE
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-(K821/
     &                            NFRONT
     &                            -NB_ROWS(i))
                             NB_ROWS(i)=K821/NFRONT
                          ENDIF
                          i=i+1
                       ENDDO
                       CHOSEN=i-1
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                         DO i=CHOSEN,1,-1
                          IF(int(dble(ADDITIONNAL_ROWS)/dble(i))
     &                          .NE.0)THEN
                             GOTO 372 
                          ENDIF
                         ENDDO
 372                     CONTINUE
                         X=int(dble(ADDITIONNAL_ROWS)/dble(i))
                         DO J=1,i
                           NB_ROWS(J)=NB_ROWS(J)+X
                           ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                         ENDDO 
                         IF(ADDITIONNAL_ROWS.NE.0) THEN
                           NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                         ENDIF    
                    ENDIF
                    ENDIF
#endif
              ENDIF
 666          CONTINUE
              SOMME=dble(0)
              X=0
              POS=0
              DO i=1,CHOSEN
                 IF(K50.NE.0) THEN
                    IF((TEMP(i)+dble(NB_ROWS(i))
     &                   *dble(X+NB_ROWS(i)+NFRONT-NCB))
     &                   .GT.PEAK)THEN
                       SMALL_SET=SMALL_SET+1
                    ENDIF
                 ENDIF
                 IF(K50.EQ.0) THEN
                    IF((TEMP(i)+dble(NB_ROWS(i))*dble(NFRONT))
     &                   .GT.PEAK)THEN
                       SMALL_SET=SMALL_SET+1
                    ENDIF
                 ENDIF
                 X=X+NB_ROWS(i)
                 SOMME=SOMME+ dble(NB_ROWS(i))
              ENDDO
           ENDIF
 889       CONTINUE
           J=CHOSEN
           X=0
           DO i=J,1,-1
             IF(NB_ROWS(i).EQ.0)THEN
                IF(X.EQ.1)THEN
                  WRITE(*,*)MYID,
     &         ':Internal error 12 in ZMUMPS_504'
                  CALL MUMPS_ABORT()
                ENDIF
                CHOSEN=CHOSEN-1
             ELSE 
                  IF(NB_ROWS(i).GT.0)THEN
                    X=1
                  ELSE
                    WRITE(*,*)
     &            'Internal error 13 in ZMUMPS_504'
                    CALL MUMPS_ABORT()
                  ENDIF
             ENDIF
          ENDDO
           NSLAVES_NODE=CHOSEN
           TAB_POS(NSLAVES_NODE+1)= NCB+1
           TAB_POS(SLAVEF+2) = CHOSEN
           POS=1
           DO i=1,CHOSEN
              SLAVES_LIST(i)=TEMP_ID(i)
              TAB_POS(i)=POS
              POS=POS+NB_ROWS(i) 
              IF(NB_ROWS(i).LE.0)THEN
                WRITE(*,*)
     &          'Internal error 14 in ZMUMPS_504'
                 CALL MUMPS_ABORT()
              ENDIF
           ENDDO
           DO i=CHOSEN+1,NUMBER_OF_PROCS
              SLAVES_LIST(i)=TEMP_ID(i)
           ENDDO
           IF(POS.NE.(NCB+1))THEN
              WRITE(*,*)
     &        'Internal error 15 in ZMUMPS_504'
             CALL MUMPS_ABORT()
           ENDIF
      END SUBROUTINE ZMUMPS_504
      SUBROUTINE ZMUMPS_518
     &      (NCBSON_MAX,SLAVEF,KEEP,KEEP8,
     &       PROCS,MEM_DISTRIB,NCB,NFRONT,
     &       NSLAVES_NODE,TAB_POS,
     &       SLAVES_LIST,SIZE_SLAVES_LIST,MYID,INODE,MP,LP)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: NCBSON_MAX
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(in) :: MP,LP
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER NUMBER_OF_PROCS,K47,K48, K50,K83,K69
      INTEGER(8) :: K821
      INTEGER J
      INTEGER KMIN, KMAX
      INTEGER OTHERS,CHOSEN,SMALL_SET,ACC
      DOUBLE PRECISION SOMME,TMP_SUM,DELTA,A,B,C,MASTER_WORK
      INTEGER AFFECTED
      INTEGER ADDITIONNAL_ROWS,i,X,REF,POS,NELIM
      INTEGER(8) X8
      LOGICAL FORCE_CAND,SMP
      DOUBLE PRECISION BANDE_K821
      INTEGER NB_SAT,NB_ZERO
      DOUBLE PRECISION TEMP(SLAVEF),TOTAL_COST, MAX_MEM_ALLOW
      INTEGER TEMP_ID(SLAVEF),NB_ROWS(SLAVEF)
      INTEGER NSLAVES_REF,NCB_FILS
      EXTERNAL MPI_WTIME,MUMPS_442
      INTEGER MUMPS_442
      INTEGER POS_MIN_LOAD,SIZE_MY_SMP,WHAT,ADDITIONNAL_ROWS_SPECIAL
      LOGICAL HAVE_TYPE1_SON
      DOUBLE PRECISION MIN_LOAD,MAX_LOAD,TEMP_MAX_LOAD
      DOUBLE PRECISION MPI_WTIME
      DOUBLE PRECISION BUF_SIZE,NELIM_MEM_SIZE
      DOUBLE PRECISION MEM_SIZE_STRONG(SLAVEF),MEM_SIZE_WEAK(SLAVEF)
      K821=abs(KEEP8(21))
      TEMP_MAX_LOAD=dble(0)
      K50=KEEP(50)
      K48=KEEP(48)
      K47=KEEP(47)
      K83=KEEP(83)
      K69=0
      NCB_FILS=NCBSON_MAX
      IF(int(NCB_FILS,8)*int(min(NCB,NCB_FILS),8).GT.K821)THEN
         HAVE_TYPE1_SON=.TRUE.
      ELSE
         HAVE_TYPE1_SON=.FALSE.
      ENDIF
      SMP=(K69.NE.0)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
      END IF
      NELIM=NFRONT-NCB
         KMAX=int(K821/int(NCB,8))
         IF(FORCE_CAND)THEN
           DO i=1,PROCS(SLAVEF+1)
              WLOAD(i)=LOAD_FLOPS(PROCS(i))
              IDWLOAD(i)=PROCS(i)
              IF (WLOAD(i) < -0.5d0 ) THEN
                 IF((MP.GT.0).AND.(LP.GE.2))THEN
                    WRITE(MP,*)MYID,': Warning: negative load ',
     &                   WLOAD(i)
                 ENDIF
              ENDIF 
              WLOAD(i)=max(WLOAD(i),0.0d0)
           ENDDO
           NUMBER_OF_PROCS=PROCS(SLAVEF+1)
            OTHERS=NUMBER_OF_PROCS
         ELSE
            NUMBER_OF_PROCS=SLAVEF
            WLOAD(1:SLAVEF) = LOAD_FLOPS(0:NUMBER_OF_PROCS-1)
            DO i=1,NUMBER_OF_PROCS
               IDWLOAD(i) = i - 1
               IF (WLOAD(i) < -0.5d0 ) THEN
                  IF((MP.GT.0).AND.(LP.GE.2))THEN
                     WRITE(MP,*)MYID,': Negative load ',
     &                    WLOAD(i)
                  ENDIF
               ENDIF 
               WLOAD(i)=max(WLOAD(i),0.0d0)
            ENDDO
            OTHERS=NUMBER_OF_PROCS-1
         ENDIF
         KMAX=int(NCB/OTHERS)
         KMIN=MUMPS_442(int(NCB,8)*int(KMAX,8),K50,KMAX,NCB)
         NB_ROWS=0
         CALL MUMPS_558(NUMBER_OF_PROCS, WLOAD, IDWLOAD)
         IF(K50.EQ.0)THEN
            TOTAL_COST=dble( NELIM ) * dble ( NCB ) +
     &           dble(NCB) * dble(NELIM)*dble(2*NFRONT-NELIM-1)
         ELSE
            TOTAL_COST=dble(NELIM) * dble ( NCB ) *
     &           dble(NFRONT+1)
         ENDIF
         CALL MUMPS_511(NFRONT,NELIM,NELIM,K50,
     &        2,MASTER_WORK)
         SOMME=dble(0)
         J=1
         IF(FORCE_CAND.AND.(NUMBER_OF_PROCS.GT.K83))THEN
            MASTER_WORK=dble(KEEP(88))*MASTER_WORK/dble(100)
         ENDIF
         IF(FORCE_CAND.AND.(NUMBER_OF_PROCS.LE.K83))THEN
            MASTER_WORK=dble(KEEP(87))*MASTER_WORK/dble(100)
         ENDIF
         IF(MASTER_WORK.LT.dble(1))THEN
            MASTER_WORK=dble(1)
         ENDIF
         NSLAVES_REF=int(TOTAL_COST/MASTER_WORK)+1
         IF(FORCE_CAND)THEN
            NSLAVES_REF=min(NSLAVES_REF,NUMBER_OF_PROCS)
         ELSE
            NSLAVES_REF=min(NSLAVES_REF,NUMBER_OF_PROCS-1)
         ENDIF
        DO i=1,NUMBER_OF_PROCS
           IF((IDWLOAD(i).NE.MYID))THEN
              TEMP_ID(J)=IDWLOAD(i)
              TEMP(J)=WLOAD(i)
              IF(BDC_M2_FLOPS)THEN
                 TEMP(J)=TEMP(J)+NIV2(TEMP_ID(J)+1)
              ENDIF
              J=J+1
           ENDIF
        ENDDO               
        NUMBER_OF_PROCS=J-1
        CALL MUMPS_558(NUMBER_OF_PROCS, TEMP, TEMP_ID)
        SOMME=dble(0)
        TMP_SUM=dble(0)
        DO i=1,OTHERS
            SOMME=SOMME+TEMP(OTHERS)-TEMP(i)
            TMP_SUM=TMP_SUM+TEMP(i)
        ENDDO
         TMP_SUM=(TMP_SUM/dble(OTHERS))+
     &        (TOTAL_COST/dble(OTHERS))
         SIZE_MY_SMP=OTHERS
         MIN_LOAD=TEMP(1)
         POS_MIN_LOAD=1
         IF(.NOT.SMP) MAX_LOAD=TEMP(OTHERS)
         IF(SMP)THEN
            J=1
            DO i=1,OTHERS
               IF(MEM_DISTRIB(TEMP_ID(i)).EQ.1)THEN
                  IF(TEMP(i).LE.TMP_SUM)THEN
                     WLOAD(J)=TEMP(i)
                     IDWLOAD(J)=TEMP_ID(i)
                     J=J+1
                  ELSE
                  ENDIF
               ENDIF
            ENDDO
            MAX_LOAD=WLOAD(J-1)
            SIZE_MY_SMP=J-1
            DO i=1,OTHERS
               IF((MEM_DISTRIB(TEMP_ID(i)).NE.1).OR.
     &              ((MEM_DISTRIB(TEMP_ID(i)).EQ.1).AND.
     &              (TEMP(i).GE.TMP_SUM)))THEN
                  WLOAD(J)=TEMP(i)
                  IDWLOAD(J)=TEMP_ID(i)
                  J=J+1
               ENDIF
            ENDDO
            TEMP=WLOAD
            TEMP_ID=IDWLOAD
         ENDIF
        IF(BDC_MD)THEN
           BUF_SIZE=dble(K821)
           IF (KEEP(201).EQ.2) THEN
              A=dble(int((dble(KEEP(100))/dble(2))/dble(NELIM)))
              IF(K50.EQ.0)THEN
                 BUF_SIZE=min(BUF_SIZE,A*dble(NCB))
              ELSE
                 BUF_SIZE=min(BUF_SIZE,A*A)
              ENDIF
           ENDIF
           BUF_SIZE=dble(K821)
           DO i=1,NUMBER_OF_PROCS
              A=dble(MD_MEM(TEMP_ID(i)))/
     &             dble(NELIM)
              A=A*dble(NFRONT)
              IF(K50.EQ.0)THEN
                 B=dble(int(dble(NCB)/dble(NUMBER_OF_PROCS))+1)*
     &                dble(NFRONT)
              ELSE
                 WHAT = 5 
#if ! defined (OOC_PAR_MEM_SLAVE_SELECT) && ! defined (OOC_PAR_MEM2)
                 CALL MUMPS_503(WHAT, KEEP,KEEP8, NCB,
     &                NFRONT, min(NCB,OTHERS), J, X8)
#endif
                 B=dble(X8)+(dble(J)*dble(NELIM))
              ENDIF
              NELIM_MEM_SIZE=A+B
              MEM_SIZE_WEAK(i)=NELIM_MEM_SIZE
            IF((SBTR_WHICH_M.EQ.0).OR.(.NOT.BDC_SBTR))THEN
               IF(BDC_M2_MEM)THEN
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)
               ELSE
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))
               ENDIF
            ELSE
               IF(BDC_SBTR)THEN
                  IF(BDC_M2_MEM)THEN
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ELSE
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ENDIF
               ENDIF
            ENDIF
            IF(min(MEM_SIZE_STRONG(i),MEM_SIZE_WEAK(i)).LT.dble(0))THEN
                IF(MEM_SIZE_STRONG(i).LT.0.0d0)THEN
                   MEM_SIZE_STRONG(i)=dble(0)
                ELSE
                   MEM_SIZE_WEAK(i)=dble(0)
                ENDIF
             ENDIF
          ENDDO
       ELSE
          BUF_SIZE=dble(K821)
          DO i=1,NUMBER_OF_PROCS
            IF((SBTR_WHICH_M.EQ.0).OR.(.NOT.BDC_SBTR))THEN
               IF(BDC_M2_MEM)THEN
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)
               ELSE
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))
               ENDIF
            ELSE
               IF(BDC_SBTR)THEN
                  IF(BDC_M2_MEM)THEN
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ELSE
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ENDIF
               ENDIF
            ENDIF
            MEM_SIZE_STRONG(i)=max(dble(0),MEM_SIZE_STRONG(i))
            MEM_SIZE_WEAK(i)=huge(MEM_SIZE_WEAK(i))
          ENDDO
       ENDIF
       IF((((NUMBER_OF_PROCS.LE.K83).AND.FORCE_CAND).AND.
     &      (TOTAL_COST.GE.SOMME)).OR.
     &      (.NOT.FORCE_CAND).OR.
     &      (((NUMBER_OF_PROCS+1).GT.K83).AND.FORCE_CAND))THEN
               REF=NSLAVES_REF
               SMALL_SET=NSLAVES_REF
               IF(.NOT.SMP)THEN
                  DO i=NSLAVES_REF,1,-1
                     SOMME=dble(0)
                     DO J=1,i
                        SOMME=SOMME+TEMP(J)
                     ENDDO
                     SOMME=(dble(i)*TEMP(i))-SOMME
                     IF(TOTAL_COST.GE.SOMME) GOTO 444
                  ENDDO
 444              CONTINUE
                  REF=i
                  SMALL_SET=REF
                  MAX_LOAD=TEMP(SMALL_SET)
               ELSE
                  X=min(SIZE_MY_SMP,NSLAVES_REF)
 450              CONTINUE
                  SOMME=dble(0)
                  DO J=1,X
                     SOMME=SOMME+(TEMP(X)-TEMP(J))
                  ENDDO
                  IF(SOMME.GT.TOTAL_COST)THEN
                     X=X-1
                     GOTO 450
                  ELSE
                     IF(X.LT.SIZE_MY_SMP) THEN
                        REF=X
                        SMALL_SET=REF
                        MAX_LOAD=TEMP(SMALL_SET)
                     ELSE
                        X=min(SIZE_MY_SMP,NSLAVES_REF)
                        J=X+1
                        MAX_LOAD=TEMP(X)
                        TMP_SUM=MAX_LOAD
                        DO i=X+1,OTHERS
                           IF(TEMP(i).GT.MAX_LOAD)THEN
                              SOMME=SOMME+(dble(i-1)*(TEMP(i)-MAX_LOAD))
                              TMP_SUM=MAX_LOAD
                              MAX_LOAD=TEMP(i)
                           ELSE
                              SOMME=SOMME+(MAX_LOAD-TEMP(i))
                           ENDIF
                           IF(i.EQ.NSLAVES_REF)THEN
                              SMALL_SET=NSLAVES_REF
                              REF=SMALL_SET
                              GOTO 323
                           ENDIF
                           IF(SOMME.GT.TOTAL_COST)THEN
                              REF=i-1
                              SMALL_SET=i-1
                              MAX_LOAD=TMP_SUM
                              GOTO 323
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
 323           CONTINUE
               MAX_LOAD=dble(0)
               DO i=1,SMALL_SET
                  MAX_LOAD=max(MAX_LOAD,TEMP(i))
               ENDDO
               TEMP_MAX_LOAD=MAX_LOAD
               NB_ROWS=0
               TMP_SUM=dble(0)
               CHOSEN=0
               ACC=0
               NB_SAT=0
               NB_ZERO=0
               DO i=1,SMALL_SET
                  IF(K50.EQ.0)THEN
                     X=int(BUF_SIZE/dble(NCB+1))-1
                     BANDE_K821=dble(X)*dble(NFRONT)
                  ELSE
                     A=dble(1)
                     B=dble(ACC+2)
                     C=-BUF_SIZE+dble(ACC+NELIM)
                     DELTA=(B*B)-(dble(4)*A*C)
                     X=int((-B+sqrt(DELTA))/(dble(2)*A))
                     IF(X.GT.NCB-ACC) X=NCB-ACC
                     BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                  ENDIF
                  IF(HAVE_TYPE1_SON)THEN
                     IF(K50.EQ.0)THEN
                        X=int((BUF_SIZE-dble(NFRONT))/dble(NFRONT+1))
                        BANDE_K821=dble(X)*dble(NFRONT)
                     ELSE
                        A=dble(1)
                        B=dble(ACC+2+NELIM)
                        C=-BUF_SIZE+dble(ACC+NELIM)
                        DELTA=(B*B)-(dble(4)*A*C)
                        X=int((-B+sqrt(DELTA))/(dble(2)*A))
                        IF(X.GT.NCB-ACC) X=NCB-ACC
                        BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                     ENDIF
                  ENDIF
                  MAX_MEM_ALLOW=BANDE_K821
                  IF(BDC_MD)THEN
                     MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                    BANDE_K821)
                     MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                  ENDIF
                  IF(K50.EQ.0)THEN
                     KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                     X=int((MAX_LOAD-TEMP(i))/
     &                    (dble(NELIM)*dble(2*NFRONT-NELIM)))
                     IF(X.GE.KMAX)THEN
                        IF(KMAX.GE.KMIN)THEN
                           X=KMAX
                           NB_SAT=NB_SAT+1
                        ELSE
                           X=0
                        ENDIF
                     ELSE
                        IF(X.LT.KMIN)THEN
                           X=0
                        ENDIF                        
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  IF(K50.NE.0)THEN
                        A=dble(1)
                        B=dble(ACC+NELIM)
                        C=dble(-MAX_MEM_ALLOW)
                        DELTA=((B*B)-(dble(4)*A*C))
                        KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                     A=dble(NELIM)
                     B=dble(NELIM)*(dble(NELIM)+dble(2*ACC+1))
                     C=-(MAX_LOAD-TEMP(i))
                     DELTA=(B*B-(dble(4)*A*C))
                     X=int((-B+sqrt(DELTA))/(dble(2)*A))
                     IF(X.LT.0) THEN
                        WRITE(*,*)MYID,
     &    ': Internal error 1 in ZMUMPS_518'
                        CALL MUMPS_ABORT()
                     ENDIF
                     IF(X.GE.KMAX)THEN
                        IF(KMAX.GE.KMIN)THEN
                           X=KMAX
                           NB_SAT=NB_SAT+1
                        ELSE
                           X=0
                        ENDIF
                     ELSE
                        IF(X.LT.KMIN)THEN
                           X=0
                        ENDIF
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF                  
                  NB_ROWS(i)=X
                  ACC=ACC+X
                  CHOSEN=CHOSEN+1
                  IF(SMP)THEN
                     IF(MIN_LOAD.GT.TEMP(i))THEN
                        MIN_LOAD=TEMP(i)
                        POS_MIN_LOAD=i
                     ENDIF
                  ENDIF
                  TMP_SUM=MAX_LOAD
                  IF(K50.EQ.0)THEN
                     MAX_LOAD=max(MAX_LOAD,
     &                    (TEMP(i)+(dble(NELIM) *
     &                    dble(NB_ROWS(i)))+
     &                    (dble(NB_ROWS(i))*dble(NELIM)*
     &                    dble(2*NFRONT-NELIM-1))))
                  ELSE
                     MAX_LOAD=max(MAX_LOAD,
     &               TEMP(i)+(dble(NELIM) * dble(NB_ROWS(i)))*
     &                    dble(2*(NELIM+ACC)-NB_ROWS(i)
     &                    -NELIM+1))
                  ENDIF
                  IF(TMP_SUM.LT.MAX_LOAD)THEN
                  ENDIF
                  IF(NCB-ACC.LT.KMIN) GOTO 888
                  IF(NCB.EQ.ACC) GOTO 888
                  IF(ACC.GT.NCB) THEN
                    WRITE(*,*)MYID,
     &      ': Internal error 2 in ZMUMPS_518'
                    CALL MUMPS_ABORT()
                  ENDIF
               ENDDO
 888           CONTINUE
               SOMME=dble(0)
               X=NFRONT-NCB
               IF((ACC.GT.NCB))THEN
                  WRITE(*,*)MYID,
     &          ': Internal error 3 in ZMUMPS_518'
                  CALL MUMPS_ABORT()
               ENDIF
               IF((ACC.LT.NCB))THEN
                  IF(K50.NE.0)THEN
                     IF(SMALL_SET.LE.OTHERS)THEN
                       IF((NB_SAT.EQ.SMALL_SET).AND.(SMALL_SET.LT.
     &                      NSLAVES_REF))THEN
                          SMALL_SET=REF+1
                          REF=REF+1
                          NB_ROWS=0
                          GOTO 323
                       ENDIF
                       ADDITIONNAL_ROWS_SPECIAL=NCB-ACC
                       DO i=1,SMALL_SET
                          MAX_LOAD=TEMP_MAX_LOAD
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=dble(NELIM)*
     &                         dble(ADDITIONNAL_ROWS)*
     &                         dble(2*NFRONT-ADDITIONNAL_ROWS-NELIM
     &                         +1)
                          SOMME=SOMME/dble(SMALL_SET-NB_SAT)
                          NB_ROWS=0
                          NB_ZERO=0
                          ACC=0
                          CHOSEN=0
                          NB_SAT=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO J=1,SMALL_SET
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                        min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=(dble(NELIM)*dble(NELIM+2*ACC+1))
                             C=-(MAX_LOAD-TEMP(J)+SOMME)
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             X=X+1
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 4 in ZMUMPS_518'
                                CALL MUMPS_ABORT()
                             ENDIF
                             IF(X.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   X=KMAX
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   NB_ZERO=NB_ZERO+1
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.min(KMIN,KMAX))THEN
                                   NB_ZERO=NB_ZERO+1
                                   X=0
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(J)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(J))THEN
                                   MIN_LOAD=TEMP(J)
                                   POS_MIN_LOAD=i
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             TMP_SUM=MAX_LOAD
                             TEMP_MAX_LOAD=max(TEMP_MAX_LOAD,
     &                            TEMP(J)+(dble(NELIM) *
     &                            dble(NB_ROWS(J)))*
     &                            dble(2*(NELIM+
     &                            ACC)-NB_ROWS(J)
     &                            -NELIM+1))
                             IF(REF.LE.NUMBER_OF_PROCS-1)THEN
                                IF(TEMP_MAX_LOAD.GT.TEMP(REF+1))THEN
                                   IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                      SMALL_SET=REF+1
                                      REF=REF+1
                                      NB_ROWS=0
                                      GOTO 323
                                   ENDIF
                                ENDIF
                             ENDIF
                             IF(NCB.EQ.ACC) GOTO 666
                          ENDDO
                          IF(NB_SAT.EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                          IF(NB_ZERO.EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                          IF((NB_SAT+NB_ZERO).EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                       ENDDO
 434                   CONTINUE
                       ADDITIONNAL_ROWS=NCB-ACC
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          IF(ADDITIONNAL_ROWS.LT.KMIN)THEN
                             i=CHOSEN
                             J=ACC
 436                         CONTINUE
                             IF(NB_ROWS(i).NE.0)THEN
                                J=J-NB_ROWS(i)
                                A=dble(1)
                                B=dble(J+2)
                                C=-BUF_SIZE+dble(J+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-J) X=NCB-J
                                BANDE_K821=dble(X)*dble(NELIM+J+X)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=dble(1)
                                   B=dble(J+2+NELIM)
                                   C=-BUF_SIZE+dble(J+NELIM)
                                   DELTA=(B*B)-(dble(4)*A*C)
                                   X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                   IF(X.GT.NCB-J) X=NCB-J
                                   BANDE_K821=dble(X)*dble(NELIM+J+X)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K821
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=min(
     &                         min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                                  BANDE_K821)
                                   MAX_MEM_ALLOW=max(dble(0),
     &                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=dble(1)
                                B=dble(J+NELIM)
                                C=dble(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(dble(4)*A*C))
                                KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(NB_ROWS(i).NE.KMAX)THEN
                                   IF(NCB-J.LE.KMAX)THEN
                                      NB_ROWS(i)=+NCB-J
                                      ADDITIONNAL_ROWS=0
                                   ENDIF
                                ENDIF
                                TEMP_MAX_LOAD=max(TEMP_MAX_LOAD,
     &                               TEMP(i)+
     &                               (dble(NELIM) * dble(NB_ROWS(i)))*
     &                               dble(2*(NELIM+
     &                               ACC)-NB_ROWS(i)
     &                               -NELIM+1))
                                IF(REF.LE.NUMBER_OF_PROCS-1)THEN
                                   IF(TEMP_MAX_LOAD.GT.TEMP(REF+1))THEN
                                      IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                         SMALL_SET=REF+1
                                         REF=REF+1
                                         NB_ROWS=0
                                         GOTO 323
                                      ENDIF
                                   ENDIF
                                ENDIF
                             ELSE
                                i=i-1
                                IF(i.NE.0)GOTO 436
                             ENDIF
                             IF(ADDITIONNAL_ROWS.NE.0)THEN
                                i=CHOSEN
                                IF(i.NE.SMALL_SET)THEN
                                   i=i+1
                                   IF(NB_ROWS(i).NE.0)THEN
                                      WRITE(*,*)MYID,
     &    ': Internal error 5 in ZMUMPS_518'
                                      CALL MUMPS_ABORT()
                                   ENDIF
                                ENDIF
                                NB_ROWS(i)=NB_ROWS(i)+ADDITIONNAL_ROWS
                                ADDITIONNAL_ROWS=0
                             ENDIF
                             CHOSEN=i
                          ENDIF
                       ENDIF
                       i=CHOSEN+1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(i.LE.NUMBER_OF_PROCS))
                          IF((TEMP(i).LE.MAX_LOAD))THEN
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=dble(NELIM)*dble(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(i))
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   X=KMAX
                                ELSE
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.KMIN)THEN
                                   X=0
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(i)=X
                             ACC=ACC+X
                             ADDITIONNAL_ROWS=NCB-ACC
                          ELSE IF((TEMP(i).GT.MAX_LOAD))THEN
                             MAX_LOAD=TEMP(i)
                             NB_SAT=0
                             ACC=0
                             NB_ROWS=0
                             DO J=1,i
                                A=dble(1)
                                B=dble(ACC+2)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=dble(1)
                                   B=dble(ACC+2+NELIM)
                                   C=-BUF_SIZE+dble(ACC+NELIM)
                                   DELTA=(B*B)-(dble(4)*A*C)
                                   X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                   IF(X.GT.NCB-ACC) X=NCB-ACC
                                   BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K821
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                                  BANDE_K821)
                                   MAX_MEM_ALLOW=max(dble(0),
     &                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=dble(1)
                                B=dble(ACC+NELIM)
                                C=dble(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(dble(4)*A*C))
                                KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                                A=dble(NELIM)
                                B=dble(NELIM)*dble(NELIM+2*ACC+1)
                                C=-(MAX_LOAD-TEMP(J))
                                DELTA=(B*B-(dble(4)*A*C))
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.LT.0) THEN
                        WRITE(*,*)MYID,
     &    ': Internal error 6 in ZMUMPS_518'
                                   CALL MUMPS_ABORT()
                                ENDIF
                                IF(X.GE.KMAX)THEN
                                   IF(KMAX.GE.KMIN)THEN
                                      X=KMAX
                                      NB_SAT=NB_SAT+1
                                   ELSE
                                      X=0
                                   ENDIF
                                ELSE
                                   IF(X.LT.min(KMIN,KMAX))THEN
                                      X=0
                                   ENDIF
                                ENDIF
                                IF((ACC+X).GT.NCB) X=NCB-ACC
                                NB_ROWS(J)=X
                                IF(SMP)THEN
                                   IF(MIN_LOAD.GT.TEMP(J))THEN
                                      MIN_LOAD=TEMP(J)
                                      POS_MIN_LOAD=i
                                   ENDIF
                                ENDIF
                                ACC=ACC+X
                                MAX_LOAD=max(MAX_LOAD,
     &                               TEMP(J)+
     &                               (dble(NELIM)*dble(NB_ROWS(J)))*
     &                               dble(2*(NELIM+
     &                               ACC)-NB_ROWS(J)
     &                               -NELIM+1))
                                IF(NCB.EQ.ACC) GOTO 741
                                IF(NCB-ACC.LT.KMIN) GOTO 210
                             ENDDO
 210                         CONTINUE
                          ENDIF
 741                      CONTINUE
                          i=i+1
                          ADDITIONNAL_ROWS=NCB-ACC
                       ENDDO
                       CHOSEN=i-1
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=dble(NELIM)*dble(ADDITIONNAL_ROWS)*
     &                         dble(2*NFRONT-ADDITIONNAL_ROWS-
     &                         NELIM+1)
                          SOMME=SOMME/dble(NUMBER_OF_PROCS)
                          NB_ROWS=0
                          ACC=0
                          CHOSEN=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO i=1,OTHERS
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=dble(NELIM)*dble(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(i)+SOMME)
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 7 in ZMUMPS_518'
                                CALL MUMPS_ABORT()
                             ENDIF
                             IF(X.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   X=KMAX
                                ELSE
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.min(KMIN,KMAX))THEN
                                   X=min(KMAX,KMIN)
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(i)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(i))THEN
                                   MIN_LOAD=TEMP(i)
                                   POS_MIN_LOAD=i
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             IF(NCB.EQ.ACC) GOTO 666
                             IF(NCB-ACC.LT.KMIN) GOTO 488
                          ENDDO
 488                      CONTINUE
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=dble(NELIM)*
     &                         dble(ADDITIONNAL_ROWS)*
     &                         dble(2*NFRONT-ADDITIONNAL_ROWS-
     &                         NELIM+1)
                          SOMME=SOMME/dble(NUMBER_OF_PROCS)
                          NB_ROWS=0
                          ACC=0
                          CHOSEN=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO i=1,OTHERS
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(BANDE_K821,
     &                               MEM_SIZE_STRONG(i))
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=dble(NELIM)*dble(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(i)+SOMME)
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 8 in ZMUMPS_518'
                                CALL MUMPS_ABORT()
                             ENDIF
                             IF(X.GE.KMAX)THEN
                                X=KMAX
                             ELSE
                                IF(X.LT.KMIN)THEN
                                   X=KMIN
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(i)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(i))THEN
                                   MIN_LOAD=TEMP(i)
                                   POS_MIN_LOAD=i
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             IF(NCB.EQ.ACC) GOTO 666
                             IF(NCB-ACC.LT.KMIN) GOTO 477
                          ENDDO
 477                      CONTINUE
                          IF(ACC.NE.NCB)THEN
                             NB_SAT=0
                             ACC=0
                             CHOSEN=0
                             IF(SMP)THEN
                                MIN_LOAD=TEMP(1)
                                POS_MIN_LOAD=1
                             ENDIF
                             DO i=1,OTHERS
                                A=dble(1)
                                B=dble(ACC+2)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=dble(1)
                                   B=dble(ACC+2+NELIM)
                                   C=-BUF_SIZE+dble(ACC+NELIM)
                                   DELTA=(B*B)-(dble(4)*A*C)
                                   X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                   IF(X.GT.NCB-ACC) X=NCB-ACC
                                   BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K821
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=min(BANDE_K821,
     &                                  MEM_SIZE_STRONG(i))
                                   MAX_MEM_ALLOW=max(dble(0),
     &                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=dble(1)
                                B=dble(ACC+NELIM)
                                C=dble(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(dble(4)*A*C))
                                KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                                X=KMAX-NB_ROWS(i)
                                IF((ACC+NB_ROWS(i)+X).GT.NCB)
     &                               X=NCB-(ACC+NB_ROWS(i))
                                NB_ROWS(i)=NB_ROWS(i)+X
                                IF((dble(NB_ROWS(i))*
     &                               dble(NB_ROWS(i)+ACC)).EQ.
     &                               BANDE_K821)THEN
                                   NB_SAT=NB_SAT+1
                                ENDIF
                                ACC=ACC+NB_ROWS(i)
                                IF(SMP)THEN
                                   IF(MIN_LOAD.GT.TEMP(i))THEN
                                      MIN_LOAD=TEMP(i)
                                      POS_MIN_LOAD=i
                                   ENDIF
                                ENDIF
                                CHOSEN=CHOSEN+1
                                IF(NCB.EQ.ACC) GOTO 666
                                IF(NCB-ACC.LT.KMIN) GOTO 834
                             ENDDO
 834                         CONTINUE
                          ENDIF
                          IF(ACC.NE.NCB)THEN
                            ADDITIONNAL_ROWS=NCB-ACC
                            SOMME=dble(NELIM)*
     &                           dble(ADDITIONNAL_ROWS)*
     &                           dble(2*NFRONT-ADDITIONNAL_ROWS-
     &                           NELIM+1)
                            SOMME=SOMME/dble(NUMBER_OF_PROCS-NB_SAT)
                            ACC=0
                            DO i=1,CHOSEN
                               A=dble(1)
                               B=dble(ACC+2)
                               C=-BUF_SIZE+dble(ACC+NELIM)
                               DELTA=(B*B)-(dble(4)*A*C)
                               X=int((-B+sqrt(DELTA))/(dble(2)*A))
                               IF(X.GT.NCB-ACC) X=NCB-ACC
                               BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                               IF(HAVE_TYPE1_SON)THEN
                                  A=dble(1)
                                  B=dble(ACC+2+NELIM)
                                  C=-BUF_SIZE+dble(ACC+NELIM)
                                  DELTA=(B*B)-(dble(4)*A*C)
                                  X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                  IF(X.GT.NCB-ACC) X=NCB-ACC
                                  BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                               ENDIF
                               IF((dble(NB_ROWS(i))*
     &                              dble(NB_ROWS(i)+ACC)).EQ.
     &                              BANDE_K821)THEN
                                  GOTO 102
                               ENDIF
                               A=dble(NELIM)
                               B=dble(NELIM)*
     &                              dble(NELIM+2*(ACC+NB_ROWS(i))+1)
                               C=-(SOMME)
                               DELTA=(B*B-(dble(4)*A*C))
                               X=int((-B+sqrt(DELTA))/(dble(2)*A))
                               A=dble(1)
                               B=dble(ACC+NELIM)
                               C=dble(-BANDE_K821)
                               DELTA=((B*B)-(dble(4)*A*C))
                               KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                               IF(X.LT.0) THEN
                                  WRITE(*,*)MYID,
     &    ': Internal error 9 in ZMUMPS_518'
                                  CALL MUMPS_ABORT()
                               ENDIF
                               IF((ACC+X+NB_ROWS(i)).GT.NCB)THEN
                                  IF((NCB-ACC).GT.KMAX)THEN
                                     NB_ROWS(i)=KMAX
                                  ELSE
                                     NB_ROWS(i)=NCB-ACC
                                  ENDIF
                               ELSE
                                  IF((NB_ROWS(i)+X).GT.KMAX)THEN
                                     NB_ROWS(i)=KMAX
                                  ELSE
                                     NB_ROWS(i)=NB_ROWS(i)+X
                                  ENDIF
                               ENDIF
 102                           CONTINUE
                               ACC=ACC+NB_ROWS(i)
                               IF(NCB.EQ.ACC) THEN
                                  CHOSEN=i
                                  GOTO 666
                               ENDIF
                               IF(NCB-ACC.LT.KMIN) THEN
                                  CHOSEN=i
                                  GOTO 007
                               ENDIF
                            ENDDO
 007                        CONTINUE
                            DO i=1,CHOSEN
                               NB_ROWS(i)=NB_ROWS(i)+1
                               ACC=ACC+1
                               IF(ACC.EQ.NCB)GOTO 666
                            ENDDO
                            IF(ACC.LT.NCB)THEN
                               IF(SMP)THEN
                                  NB_ROWS(1)=NB_ROWS(1)+NCB-ACC
                               ELSE
                                  NB_ROWS(POS_MIN_LOAD)=
     &                                 NB_ROWS(POS_MIN_LOAD)+NCB-ACC
                               ENDIF
                            ENDIF
                         ENDIF
                         GOTO 666
                     ENDIF
                  ENDIF
                  GOTO 666
                 ENDIF
                 ADDITIONNAL_ROWS=NCB-ACC
                 i=CHOSEN+1
                 IF(NB_SAT.EQ.SMALL_SET) GOTO 777
                 DO i=1,SMALL_SET
                    IDWLOAD(i)=i
                    AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                    BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                    IF(HAVE_TYPE1_SON)THEN
                       AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                      (dble(NFRONT+1)))
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                    ENDIF
                    MAX_MEM_ALLOW=BANDE_K821
                    IF(BDC_MD)THEN
                       MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                      BANDE_K821)
                       MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                    ENDIF
                    WLOAD(i)=MAX_MEM_ALLOW
                 ENDDO
                 CALL MUMPS_558(SMALL_SET, WLOAD, IDWLOAD)
                 NB_ZERO=0
                 IF((NB_SAT.EQ.SMALL_SET).AND.
     &                (SMALL_SET.LT.NSLAVES_REF))THEN
                    SMALL_SET=REF+1
                    REF=REF+1
                    NB_ROWS=0
                    GOTO 323
                 ENDIF
                 IF((NB_SAT.EQ.SMALL_SET).AND.
     &                (SMALL_SET.LE.NUMBER_OF_PROCS))GOTO 777
                 AFFECTED=int(ADDITIONNAL_ROWS/(SMALL_SET-NB_SAT))
                 AFFECTED=max(AFFECTED,1)
                 DO i=1,SMALL_SET
                    KMAX=int(WLOAD(i)/dble(NFRONT))
                    IF(NB_ROWS(IDWLOAD(i)).EQ.KMAX)THEN 
                       GOTO 912
                    ENDIF
                    IF((NB_ROWS(IDWLOAD(i))+min(AFFECTED,
     &                   ADDITIONNAL_ROWS)).GT.KMAX)THEN
                       IF(NB_ROWS(IDWLOAD(i)).GT.KMAX)THEN
                       ENDIF
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                      (KMAX-NB_ROWS(IDWLOAD(i)))
                       NB_ROWS(IDWLOAD(i))=KMAX
                       NB_SAT=NB_SAT+1
                       IF(NB_SAT.EQ.SMALL_SET)THEN
                          IF(SMALL_SET.NE.NSLAVES_REF)THEN
                             SMALL_SET=REF+1
                             REF=REF+1
                             NB_ROWS=0
                             GOTO 323
                          ELSE
                             MAX_LOAD=max(MAX_LOAD,
     &                            (TEMP(IDWLOAD(i))+(dble(NELIM) *
     &                            dble(NB_ROWS(IDWLOAD(i))))+
     &                            (dble(NB_ROWS(IDWLOAD(i)))*
     &                            dble(NELIM))*
     &                            dble(2*NFRONT-NELIM-1)))
                             GOTO 777
                          ENDIF
                       ENDIF
                       AFFECTED=int(ADDITIONNAL_ROWS/(SMALL_SET-NB_SAT))
                       AFFECTED=max(AFFECTED,1)
                    ELSE
                       IF((NB_ROWS(IDWLOAD(i))+min(AFFECTED,
     &                      ADDITIONNAL_ROWS)).GE.KMIN)THEN
                          X=min(AFFECTED,ADDITIONNAL_ROWS)
                          NB_ROWS(IDWLOAD(i))=NB_ROWS(IDWLOAD(i))+
     &                         X
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                       ELSE
                          X=int((MAX_LOAD-TEMP(IDWLOAD(i)))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          IF(X+AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ELSE
                             X=AFFECTED+X
                          ENDIF
                          IF(X.GE.KMIN)THEN
                             NB_ROWS(IDWLOAD(i))=NB_ROWS(IDWLOAD(i))+
     &                            X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            X
                          ELSE
                             NB_ZERO=NB_ZERO+1
                          ENDIF
                       ENDIF
                    ENDIF
 912                CONTINUE
                    MAX_LOAD=max(MAX_LOAD,
     &                   (TEMP(IDWLOAD(i))+(dble(NELIM)*
     &                   dble(NB_ROWS(IDWLOAD(i))))+
     &                   (dble(NB_ROWS(IDWLOAD(i)))*dble(NELIM))*
     &                   dble(2*NFRONT-NELIM-1)))
                    IF(SMALL_SET.LT.NUMBER_OF_PROCS)THEN
                       IF(MAX_LOAD.GT.TEMP(SMALL_SET+1))THEN
                          IF(SMALL_SET.LT.NSLAVES_REF)THEN
                             SMALL_SET=REF+1
                             REF=REF+1
                             NB_ROWS=0
                             GOTO 323
                          ENDIF
                       ENDIF
                    ENDIF
                    IF(SMALL_SET.EQ.NB_SAT)GOTO 777
                    IF(ADDITIONNAL_ROWS.EQ.0)THEN
                       CHOSEN=SMALL_SET
                       GOTO 049
                    ENDIF
                 ENDDO
 777             CONTINUE
                 IF((NB_ZERO.NE.0).AND.(ADDITIONNAL_ROWS.GE.KMIN))THEN
                    J=NB_ZERO
 732                CONTINUE
                    X=int(ADDITIONNAL_ROWS/(J))
                    IF(X.LT.KMIN)THEN
                       J=J-1
                       GOTO 732
                    ENDIF
                    IF(X*J.LT.ADDITIONNAL_ROWS)THEN
                       X=X+1
                    ENDIF
                    DO i=1,SMALL_SET
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                         dble(BANDE_K821))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF(NB_ROWS(i).EQ.0)THEN
                          IF(X.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ENDIF
                          IF(X.GT.KMAX)THEN
                             X=KMAX
                          ENDIF
                          IF(X.GT.KMIN)THEN
                             NB_ROWS(i)=X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                             MAX_LOAD=max(MAX_LOAD,
     &                            (TEMP(i)+(dble(NELIM) *
     &                            dble(NB_ROWS(i)))+
     &                            (dble(NB_ROWS(i))*dble(NELIM))*
     &                            dble(2*NFRONT-NELIM-1)))
                          ENDIF
                       ENDIF
                    ENDDO
                 ENDIF
                 i=CHOSEN+1
                 DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                .AND.(i.LE.NUMBER_OF_PROCS))
                    IF((TEMP(i).LE.MAX_LOAD))THEN
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                         BANDE_K821)
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       AFFECTED=int((MAX_LOAD-TEMP(i))/
     &                      (dble(NELIM)*dble(2*NFRONT-NELIM)))
                       IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                          AFFECTED=ADDITIONNAL_ROWS
                       ENDIF
                       IF(NB_ROWS(i).LT.KMAX)THEN
                          IF((AFFECTED+NB_ROWS(i)).GT.KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(i)
                             NB_SAT=NB_SAT+1
                          ELSE
                             IF((AFFECTED+NB_ROWS(i)).LT.
     &                            KMIN)THEN
                                AFFECTED=0
                             ENDIF
                          ENDIF
                          NB_ROWS(i)=NB_ROWS(i)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                       ENDIF
                    ELSE IF((TEMP(i).GT.MAX_LOAD))THEN
                       IF(NB_SAT.EQ.i-1) GOTO 218
                       X=(ADDITIONNAL_ROWS/(i-1-NB_SAT))
                       ACC=1
                       DO J=1,i-1
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(J)+X))
     &                         +(dble(NB_ROWS(J)+X)*dble(NELIM))*
     &                         dble(2*NFRONT-NELIM-1))
                          IF((TEMP(J)+TMP_SUM).GT.MAX_LOAD)THEN
                             ACC=0
                          ENDIF
                       ENDDO
                       IF(ACC.EQ.1)THEN
                          MAX_LOAD=TEMP(i)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                         .AND.(J.LT.i))
                             AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             IF(HAVE_TYPE1_SON)THEN
                                AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             ENDIF
                             AFFECTED=X
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                                AFFECTED=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                                   AFFECTED=KMAX-NB_ROWS(J)
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   IF((AFFECTED+NB_ROWS(J)).LT.
     &                                  KMIN)THEN
                                      AFFECTED=0
                                   ENDIF
                                ENDIF
                                NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               AFFECTED
                             ENDIF
                             J=J+1
                          ENDDO                          
                       ELSE
                          MAX_LOAD=TEMP(i)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                         .AND.(J.LT.i))
                             AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             IF(HAVE_TYPE1_SON)THEN
                                AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             ENDIF
                             TMP_SUM=((dble(NELIM)* dble(NB_ROWS(J)))
     &                            +(dble(NB_ROWS(J))*dble(NELIM))*
     &                            dble(2*NFRONT-NELIM-1))
                             X=int((MAX_LOAD-(TEMP(J)+TMP_SUM))/
     &                            (dble(NELIM)*dble(2*NFRONT-NELIM)))
                             IF(X.LT.0)THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 10 in ZMUMPS_518'
                                CALL MUMPS_ABORT()
                             ENDIF
                             AFFECTED=X
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                                AFFECTED=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                                   AFFECTED=KMAX-NB_ROWS(J)
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   IF((AFFECTED+NB_ROWS(J)).LT.
     &                                  KMIN)THEN
                                      AFFECTED=0
                                   ENDIF
                                ENDIF
                                NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               AFFECTED
                             ENDIF
                             J=J+1
                          ENDDO
                       ENDIF
                    ENDIF
 218                CONTINUE
                    i=i+1
                 ENDDO
                 CHOSEN=i-1
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     &                 (ADDITIONNAL_ROWS.NE.0))THEN
                    DO i=1,CHOSEN
                       IF(NB_ROWS(i)+1.GE.KMIN)THEN
                          NB_ROWS(i)=NB_ROWS(i)+1
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                       ENDIF
                       MAX_LOAD=max(MAX_LOAD,
     &                      (TEMP(i)+(dble(NELIM) *
     &                      dble(NB_ROWS(i)))+
     &                      (dble(NB_ROWS(i))*dble(NELIM))*
     &                      dble(2*NFRONT-NELIM-1)))
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 048
                    ENDDO
 048                CONTINUE
                 ENDIF
                 IF((ADDITIONNAL_ROWS.NE.0))THEN
                    IF(CHOSEN.LT.NUMBER_OF_PROCS)THEN
                       i=CHOSEN+1
                    ELSE
                       IF(CHOSEN.NE.NUMBER_OF_PROCS)THEN
                          WRITE(*,*)MYID,
     &    ': Internal error 11 in ZMUMPS_518'
                          CALL MUMPS_ABORT()
                       ENDIF
                       i=CHOSEN
                    ENDIF
                    DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                   .AND.(i.LE.NUMBER_OF_PROCS))
                       IF(TEMP(i).LE.MAX_LOAD)THEN
                          AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                            BANDE_K821)
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(i)))
     &                         +(dble(NB_ROWS(i))*dble(NELIM))*
     &                         dble(2*NFRONT-NELIM-1))
                          X=int((MAX_LOAD-(TEMP(i)+TMP_SUM))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          AFFECTED=X
                          IF(X.LT.0)THEN
                             WRITE(*,*)MYID,
     &    ': Internal error 12 in ZMUMPS_518'
                             CALL MUMPS_ABORT()
                          ENDIF
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          IF(NB_ROWS(i).LT.KMAX)THEN
                             IF((AFFECTED+NB_ROWS(i)).GT.KMAX)THEN
                                AFFECTED=KMAX-NB_ROWS(i)
                             ELSE
                                IF((AFFECTED+NB_ROWS(i)).LT.
     &                               KMIN)THEN
                                   AFFECTED=0
                                ENDIF
                             ENDIF
                             NB_ROWS(i)=NB_ROWS(i)+AFFECTED
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          ENDIF
                          IF(i.NE.NUMBER_OF_PROCS) GOTO 624
                       ELSE IF((TEMP(i).GT.MAX_LOAD))THEN
                          X=int(ADDITIONNAL_ROWS/i-1)
                          X=max(X,1)
                          IF((MAX_LOAD+((dble(NELIM)*
     &                         dble(X))+(dble(
     &                         X)*dble(NELIM))*dble(
     &                         (2*NFRONT-NELIM-1)))).LE.TEMP(i))THEN
                             AFFECTED=X
                             POS=1
                          ELSE
                             POS=0
                          ENDIF
                          MAX_LOAD=TEMP(i)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                         .AND.(J.LT.i))
                             X=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(X)*dble(NFRONT)
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(HAVE_TYPE1_SON)THEN
                                X=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(X)*dble(NFRONT)
                             ENDIF
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(POS.EQ.0)THEN
                                TMP_SUM=((dble(NELIM) *
     &                               dble(NB_ROWS(J)))
     &                               +(dble(NB_ROWS(J))*dble(NELIM))*
     &                               dble(2*NFRONT-NELIM-1))
                                X=int((TEMP(i)-(TEMP(J)+TMP_SUM))/
     &                               (dble(NELIM)*dble(2*NFRONT-
     &                               NELIM)))
                             ELSE
                                X=int(TMP_SUM)
                             ENDIF
                             IF(X.GT.ADDITIONNAL_ROWS)THEN
                                X=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((X+NB_ROWS(J)).GT.KMAX)THEN
                                   X=KMAX-NB_ROWS(J)
                                ELSE
                                   IF((NB_ROWS(J)+X).LT.
     &                                  KMIN)THEN
                                     X=0
                                  ENDIF
                               ENDIF
                               NB_ROWS(J)=NB_ROWS(J)+X
                               ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                            ENDIF
                            J=J+1
                         ENDDO
                       ENDIF
 624                   CONTINUE
                       i=i+1
                    ENDDO
                    CHOSEN=i-1
                    IF(ADDITIONNAL_ROWS.NE.0)THEN
                       ACC=0
                       DO i=1,CHOSEN
                          X=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(X)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             X=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(X)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                            BANDE_K821)
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(i)))
     &                    +(dble(NB_ROWS(i))*dble(NELIM))*
     &                    dble(2*NFRONT-NELIM-1))
                          X=int((MAX_LOAD-
     &                         (TEMP(i)+TMP_SUM))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          IF(X.LT.0)THEN
                             WRITE(*,*)MYID,
     &    ': Internal error 13 in ZMUMPS_518'
                             CALL MUMPS_ABORT()
                          ENDIF
                          IF(X.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ENDIF
                          IF(NB_ROWS(i).LT.KMAX)THEN
                             IF((X+NB_ROWS(i)).GE.KMAX)THEN
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               (KMAX-NB_ROWS(i))
                                NB_ROWS(i)=KMAX
                             ELSE
                                IF((X+NB_ROWS(i)).GE.
     &                               KMIN)THEN
                                   NB_ROWS(i)=NB_ROWS(i)+X
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                                   ACC=ACC+1
                                ELSE
                                   ACC=ACC+1
                                ENDIF
                             ENDIF
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       IF(CHOSEN.LT.NUMBER_OF_PROCS)THEN
                          CHOSEN=CHOSEN+1
                       ENDIF
                       IF(ACC.EQ.0)THEN
                          ACC=1
                       ENDIF
                       X=int(ADDITIONNAL_ROWS/ACC)
                       X=max(X,1)
                       ACC=0
                       DO i=1,CHOSEN
                          J=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(J)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             J=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(J)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                            BANDE_K821)
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(i)))
     &                         +(dble(NB_ROWS(i))*dble(NELIM))*
     &                         dble(2*NFRONT-NELIM-1))
                          J=int((MAX_LOAD-
     &                         (TEMP(i)+TMP_SUM))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          IF(NB_ROWS(i).LT.KMAX)THEN
                             IF((min(X,J)+NB_ROWS(i)).GE.KMAX)THEN
                                IF((KMAX-NB_ROWS(i)).GT.
     &                               ADDITIONNAL_ROWS)THEN
                                   NB_ROWS(i)=NB_ROWS(i)+
     &                                  ADDITIONNAL_ROWS
                                   ADDITIONNAL_ROWS=0
                                ELSE
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                                  (KMAX-NB_ROWS(i))
                                   NB_ROWS(i)=KMAX
                                ENDIF
                             ELSE
                                IF((min(X,J)+NB_ROWS(i)).GE.
     &                            KMIN)THEN
                                   NB_ROWS(i)=NB_ROWS(i)+min(X,J)
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                                  min(X,J)
                                   ACC=ACC+1
                                ENDIF
                             ENDIF
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       IF(ACC.GT.0)THEN
                          DO i=1,CHOSEN
                             X=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(X)*dble(NFRONT)
                             IF(HAVE_TYPE1_SON)THEN
                                X=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(X)*dble(NFRONT)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(KMAX-NB_ROWS(i).LT.
     &                            ADDITIONNAL_ROWS)THEN
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               (KMAX-NB_ROWS(i))
                                NB_ROWS(i)=KMAX
                             ELSE
                                IF(NB_ROWS(i).EQ.0)THEN
                                   IF(min(KMIN,KMAX).LT.
     &                                  ADDITIONNAL_ROWS)THEN
                                      NB_ROWS(i)=min(KMIN,KMAX)
                                      ADDITIONNAL_ROWS=
     &                                     ADDITIONNAL_ROWS-
     &                                     min(KMIN,KMAX)
                                   ENDIF
                                ELSE
                                   NB_ROWS(i)=NB_ROWS(i)+
     &                                  ADDITIONNAL_ROWS
                                   ADDITIONNAL_ROWS=0
                                ENDIF
                             ENDIF
                             IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                          ENDDO
                       ENDIF
                       DO i=1,CHOSEN
                          IDWLOAD(i)=i
                          AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          ENDIF                          
                          WLOAD(i)=(BANDE_K821-dble(NB_ROWS(i)*NFRONT))
                       ENDDO
                       CALL MUMPS_558(NUMBER_OF_PROCS, WLOAD,
     &                      IDWLOAD)
                       NB_SAT=0
                       DO i=1,CHOSEN
                          X=int(ADDITIONNAL_ROWS/(CHOSEN-NB_SAT))
                          X=max(X,1)
                          AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          ENDIF
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(BANDE_K821,
     &                            MEM_SIZE_STRONG(i))
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          IF(NB_ROWS(IDWLOAD(i)).LT.KMAX)THEN
                             IF((NB_ROWS(IDWLOAD(i))+X).LT.KMAX)THEN
                                NB_ROWS(IDWLOAD(i))=
     &                               NB_ROWS(IDWLOAD(i))+X
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                             ELSE
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               (KMAX-NB_ROWS(IDWLOAD(i)))
                                NB_ROWS(IDWLOAD(i))=KMAX
                             ENDIF
                          ENDIF
                          IF(NB_ROWS(IDWLOAD(i)).EQ.KMAX)THEN
                             NB_SAT=NB_SAT+1
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0) GOTO 049
                       ENDDO
                       DO i=1,CHOSEN
                          X=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(X)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             X=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(X)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(BANDE_K821,
     &                            MEM_SIZE_STRONG(i))
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          IF(KMAX-NB_ROWS(i).LT.ADDITIONNAL_ROWS)THEN
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            (KMAX-NB_ROWS(i))
                             NB_ROWS(i)=KMAX
                          ELSE
                             NB_ROWS(i)=NB_ROWS(i)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       X=int(ADDITIONNAL_ROWS/CHOSEN)
                       X=max(X,1)
                       DO i=1,CHOSEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          NB_ROWS(i)=NB_ROWS(i)+X
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                    ENDIF
                 ENDIF
 049             CONTINUE
              ENDIF
 666          CONTINUE
              SOMME=dble(0)
              X=0
              POS=0
              DO i=1,CHOSEN
                 X=X+NB_ROWS(i)
                 SOMME=SOMME+ dble(NB_ROWS(i))
              ENDDO
              GOTO 890
           ELSE IF((KEEP(83).GE.NUMBER_OF_PROCS).AND.FORCE_CAND)THEN
              MAX_LOAD=dble(0)
              DO i=1,OTHERS
                 MAX_LOAD=max(MAX_LOAD,TEMP(i))
              ENDDO
              ACC=0
              CHOSEN=0
              X=1
              DO i=1,OTHERS
              ENDDO
              DO i=2,OTHERS 
                 IF(TEMP(i).EQ.TEMP(1))THEN
                    X=X+1
                 ELSE 
                    GOTO 329
                 ENDIF
              ENDDO
 329          CONTINUE
              TMP_SUM=TOTAL_COST/dble(X)
              TEMP_MAX_LOAD=dble(0)
              DO i=1,OTHERS
                 IF(K50.EQ.0)THEN
                    X=int(BUF_SIZE/dble(NCB+1))-1
                    BANDE_K821=dble(X)*dble(NFRONT)
                 ELSE
                    A=dble(1)
                    B=dble(ACC+2)
                    C=-BUF_SIZE+dble(ACC+NELIM)
                    DELTA=(B*B)-(dble(4)*A*C)
                    X=int((-B+sqrt(DELTA))/(dble(2)*A))
                    IF(X.GT.NCB-ACC) X=NCB-ACC
                    BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                 ENDIF
                 IF(HAVE_TYPE1_SON)THEN
                    IF(K50.EQ.0)THEN
                       X=int((BUF_SIZE-dble(NFRONT))/dble(NFRONT+1))
                       BANDE_K821=dble(X)*dble(NFRONT)
                    ELSE
                       A=dble(1)
                       B=dble(ACC+2+NELIM)
                       C=-BUF_SIZE+dble(ACC+NELIM)
                       DELTA=(B*B)-(dble(4)*A*C)
                       X=int((-B+sqrt(DELTA))/(dble(2)*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                    ENDIF
                 ENDIF
                 MAX_MEM_ALLOW=BANDE_K821
                 IF(BDC_MD)THEN
                    MAX_MEM_ALLOW=min(BANDE_K821,
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)))
                    MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                 ENDIF
                 IF(K50.EQ.0)THEN
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                    IF(TMP_SUM+TEMP(i).GT.MAX_LOAD)THEN
                       SOMME=MAX_LOAD-TEMP(i)
                    ELSE
                       SOMME=TMP_SUM
                    ENDIF
                    X=int(SOMME/
     &                   (dble(NELIM)*dble(2*NFRONT-NELIM)))
                    IF(X.GT.KMAX)THEN
                       X=KMAX
                    ELSE
                       IF(X.LT.KMIN)THEN
                          X=min(KMIN,KMAX)
                       ENDIF
                    ENDIF
                    IF((ACC+X).GT.NCB) X=NCB-ACC
                 ENDIF
                 IF(K50.NE.0)THEN
                       A=dble(1)
                       B=dble(ACC+NELIM)
                       C=dble(-MAX_MEM_ALLOW)
                       DELTA=((B*B)-(dble(4)*A*C))
                       KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                    A=dble(NELIM)
                    B=dble(NELIM)*dble(NELIM+2*ACC+1)
                    IF(TMP_SUM+TEMP(i).GT.MAX_LOAD)THEN
                       C=-(MAX_LOAD-TEMP(i))
                    ELSE
                       C=-TMP_SUM
                    ENDIF
                    DELTA=(B*B-(dble(4)*A*C))
                    X=int((-B+sqrt(DELTA))/(dble(2)*A))
                    IF(X.LT.0) THEN
                       WRITE(*,*)MYID,
     &    ': Internal error 14 in ZMUMPS_518'
                       CALL MUMPS_ABORT()
                    ENDIF
                    IF(X.GE.KMAX)THEN
                       IF(KMAX.GT.KMIN)THEN
                          X=KMAX
                       ELSE
                          X=0
                       ENDIF
                    ELSE
                       IF(X.LE.min(KMIN,KMAX))THEN
                          IF(KMAX.LT.KMIN)THEN
                             X=0
                          ELSE
                             X=min(KMIN,KMAX)
                          ENDIF
                       ENDIF
                    ENDIF
                    IF((ACC+X).GT.NCB) X=NCB-ACC
                 ENDIF
                 TEMP_MAX_LOAD=max(TEMP_MAX_LOAD,TEMP(i))
                 NB_ROWS(i)=X
                 CHOSEN=CHOSEN+1
                 ACC=ACC+X
                 IF(ACC.EQ.NCB) GOTO 541
              ENDDO
 541          CONTINUE
              IF(ACC.LT.NCB)THEN
                 IF(K50.EQ.0)THEN
                    ADDITIONNAL_ROWS=NCB-ACC
                    DO J=1,CHOSEN
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                         dble(BANDE_K821))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF((NB_ROWS(J)).LT.KMAX)THEN
                          IF(ADDITIONNAL_ROWS.GT.(KMAX-NB_ROWS(J)))THEN
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            (KMAX-NB_ROWS(J))
                             NB_ROWS(J)=KMAX
                          ELSE
                             NB_ROWS(J)=NB_ROWS(J)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ENDIF
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO 
                    X=int(ADDITIONNAL_ROWS/CHOSEN)
                    X=max(X,1)
                    DO J=1,CHOSEN
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(BANDE_K821,
     &                         MEM_SIZE_STRONG(J))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF((NB_ROWS(J)+X).GT.KMAX)THEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                         (KMAX-NB_ROWS(J))
                          NB_ROWS(J)=KMAX
                       ELSE
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          NB_ROWS(J)=NB_ROWS(J)+X
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO 
                    DO i=1,CHOSEN
                       X=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(X)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          X=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(X)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(BANDE_K821,
     &                         MEM_SIZE_STRONG(i))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF(KMAX-NB_ROWS(i).LT.ADDITIONNAL_ROWS)THEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                         (KMAX-NB_ROWS(i))
                          NB_ROWS(i)=KMAX
                       ELSE
                          NB_ROWS(i)=NB_ROWS(i)+ADDITIONNAL_ROWS
                          ADDITIONNAL_ROWS=0
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO
                    DO i=1,NUMBER_OF_PROCS
                       IDWLOAD(i)=i
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF                          
                       WLOAD(i)=(BANDE_K821-(dble(NB_ROWS(i))*
     &                      dble(NFRONT)))
                    ENDDO
                    CALL MUMPS_558(NUMBER_OF_PROCS, WLOAD,
     &                   IDWLOAD)
                    NB_SAT=0
                    DO i=1,CHOSEN
                       X=int(ADDITIONNAL_ROWS/(CHOSEN-NB_SAT))
                       X=max(X,1)
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF(NB_ROWS(IDWLOAD(i)).LT.KMAX)THEN
                          IF((NB_ROWS(IDWLOAD(i))+X).LT.KMAX)THEN
                             NB_ROWS(IDWLOAD(i))=
     &                            NB_ROWS(IDWLOAD(i))+X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          ELSE
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            (KMAX-NB_ROWS(IDWLOAD(i)))
                             NB_ROWS(IDWLOAD(i))=KMAX
                          ENDIF
                       ENDIF
                       IF(NB_ROWS(IDWLOAD(i)).EQ.KMAX)THEN
                          NB_SAT=NB_SAT+1
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 889
                    ENDDO
                    GOTO 994
                 ELSE
                    ACC=0
                    CHOSEN=0
                    DO i=1,OTHERS
                       A=dble(1)
                       B=dble(ACC+2)
                       C=-BUF_SIZE+dble(ACC+NELIM)
                       DELTA=(B*B)-(dble(4)*A*C)
                       X=int((-B+sqrt(DELTA))/(dble(2)*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                       IF(HAVE_TYPE1_SON)THEN
                          A=dble(1)
                          B=dble(ACC+2+NELIM)
                          C=-BUF_SIZE+dble(ACC+NELIM)
                          DELTA=(B*B)-(dble(4)*A*C)
                          X=int((-B+sqrt(DELTA))/(dble(2)*A))
                          IF(X.GT.NCB-ACC) X=NCB-ACC
                          BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(BANDE_K821,
     &                         MEM_SIZE_STRONG(i))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       A=dble(1)
                       B=dble(ACC+NELIM)
                       C=dble(-MAX_MEM_ALLOW)
                       DELTA=((B*B)-(dble(4)*A*C))
                       KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                       X=KMAX-NB_ROWS(i)
                       IF((ACC+NB_ROWS(i)+X).GT.NCB)
     &                            X=NCB-(ACC+NB_ROWS(i))
                       NB_ROWS(i)=NB_ROWS(i)+X
                       ACC=ACC+NB_ROWS(i)
                       CHOSEN=CHOSEN+1
                       IF(NCB.EQ.ACC) GOTO 889
                    ENDDO
                    ADDITIONNAL_ROWS=NCB-ACC
                 ENDIF
                 ACC=0
                 CHOSEN=0
                 DO i=1,OTHERS
                    A=dble(1)
                    B=dble(ACC+2)
                    C=-BUF_SIZE+dble(ACC+NELIM)
                    DELTA=(B*B)-(dble(4)*A*C)
                    X=int((-B+sqrt(DELTA))/(dble(2)*A))
                    IF(X.GT.NCB-ACC) X=NCB-ACC
                    BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                    IF(HAVE_TYPE1_SON)THEN
                       A=dble(1)
                       B=dble(ACC+2+NELIM)
                       C=-BUF_SIZE+dble(ACC+NELIM)
                       DELTA=(B*B)-(dble(4)*A*C)
                       X=int((-B+sqrt(DELTA))/(dble(2)*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                    ENDIF
                    MAX_MEM_ALLOW=BANDE_K821
                    A=dble(1)
                    B=dble(ACC+NELIM)
                    C=dble(-MAX_MEM_ALLOW)
                    DELTA=((B*B)-(dble(4)*A*C))
                    KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                    X=KMAX-NB_ROWS(i)
                    IF((ACC+NB_ROWS(i)+X).GT.NCB)
     &                   X=NCB-(ACC+NB_ROWS(i))
                    NB_ROWS(i)=NB_ROWS(i)+X
                    ACC=ACC+NB_ROWS(i)
                    CHOSEN=CHOSEN+1
                    IF(NCB.EQ.ACC) GOTO 889
                 ENDDO
                 ADDITIONNAL_ROWS=NCB-ACC
 994             CONTINUE
                 X=int(dble(ADDITIONNAL_ROWS)/dble(OTHERS))
                 IF((X*OTHERS).LT.ADDITIONNAL_ROWS)THEN
                    X=X+1
                 ENDIF
                 DO i=1,OTHERS
                    NB_ROWS(i)=NB_ROWS(i)+X
                    ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                    IF(ADDITIONNAL_ROWS.LT.X)X=ADDITIONNAL_ROWS
                 ENDDO
                 CHOSEN=OTHERS
              ENDIF
           ENDIF
 889       CONTINUE
           MAX_LOAD=TEMP_MAX_LOAD
 890       CONTINUE
           J=CHOSEN
           X=0
              DO i=J,1,-1
                 IF(NB_ROWS(i).EQ.0)THEN
                    CHOSEN=CHOSEN-1
                    ELSE 
                       IF(NB_ROWS(i).GT.0)THEN
                          X=1
                       ELSE
                          WRITE(*,*)MYID,
     &    ': Internal error 15 in ZMUMPS_518'
                          CALL MUMPS_ABORT()
                       ENDIF
                    ENDIF
                 ENDDO
           NSLAVES_NODE=CHOSEN
           TAB_POS(NSLAVES_NODE+1)= NCB+1
           TAB_POS(SLAVEF+2) = CHOSEN
           POS=1
           X=1
           DO i=1,J
              IF(NB_ROWS(i).NE.0)THEN
                 SLAVES_LIST(X)=TEMP_ID(i)
                 TAB_POS(X)=POS
                 POS=POS+NB_ROWS(i) 
                 IF(NB_ROWS(i).LE.0)THEN
                    WRITE(*,*)MYID,
     &    ': Internal error 16 in ZMUMPS_518'
                    CALL MUMPS_ABORT()
                 ENDIF
                 X=X+1
               ENDIF
           ENDDO
           DO i=CHOSEN+1,NUMBER_OF_PROCS
              SLAVES_LIST(i)=TEMP_ID(i)
           ENDDO
           IF(POS.NE.(NCB+1))THEN
              WRITE(*,*)MYID,
     &    ': Internal error 17 in ZMUMPS_518',
     &             POS,NCB+1
             CALL MUMPS_ABORT()
           ENDIF
      END SUBROUTINE ZMUMPS_518
      SUBROUTINE ZMUMPS_520
     &      (INODE,UPPER,SLAVEF,KEEP,KEEP8,
     &       STEP,POOL,LPOOL,PROCNODE,N)
      IMPLICIT NONE
      INTEGER INODE, LPOOL, SLAVEF, N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(KEEP(28)), POOL(LPOOL), PROCNODE(KEEP(28))
      LOGICAL UPPER
      INTEGER J
      DOUBLE PRECISION MEM_COST
      INTEGER NBINSUBTREE,i,NBTOP
      EXTERNAL ZMUMPS_508,
     & MUMPS_170
      LOGICAL ZMUMPS_508,
     & MUMPS_170
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      IF(KEEP(47).LT.2)THEN
         WRITE(*,*)'ZMUMPS_520 must
     &        be called with K47>=2'
         CALL MUMPS_ABORT()
      ENDIF        
      IF((INODE.GT.0).AND.(INODE.LE.N))THEN
      MEM_COST=ZMUMPS_543(INODE)
         IF((DM_MEM(MYID)+dble(MEM_COST)+ PEAK_SBTR_CUR_LOCAL-
     &        SBTR_CUR_LOCAL)
     &        .GT.MAX_PEAK_STK)THEN
            DO i=NBTOP-1,1,-1
               INODE = POOL( LPOOL - 2 - i)
               MEM_COST=ZMUMPS_543(INODE)
               IF((INODE.LT.0).OR.(INODE.GT.N)) THEN
                  DO J=i+1,NBTOP,-1
                     POOL(J-1)=POOL(J)
                  ENDDO
                  UPPER=.TRUE.
                  RETURN
               ENDIF
               IF((DM_MEM(MYID)+dble(MEM_COST)+PEAK_SBTR_CUR_LOCAL-
     &              SBTR_CUR_LOCAL).LE.
     &              MAX_PEAK_STK) THEN
                  DO J=i+1,NBTOP,-1
                     POOL(J-1)=POOL(J)
                  ENDDO
                  UPPER=.TRUE.
                  RETURN
               ENDIF
            ENDDO
            IF(NBINSUBTREE.NE.0)THEN
               INODE = POOL( NBINSUBTREE )
               IF(.NOT.MUMPS_170(PROCNODE(STEP(INODE)),
     &              SLAVEF))THEN
                  WRITE(*,*)
     &        'Internal error 1 in ZMUMPS_520'
                  CALL MUMPS_ABORT()
               ENDIF
               UPPER=.FALSE.
               RETURN
            ENDIF
            INODE=POOL(LPOOL-2-NBTOP)
            UPPER=.TRUE.
            RETURN
         ENDIF
      ENDIF
      UPPER=.TRUE.
      END SUBROUTINE ZMUMPS_520
      SUBROUTINE ZMUMPS_513(WHAT)
      IMPLICIT NONE
      LOGICAL WHAT
      IF(.NOT.BDC_POOL_MNG)THEN
         WRITE(*,*)'ZMUMPS_513
     &        should be called when K81>0 and K47>2'
      ENDIF
      IF(WHAT)THEN
         PEAK_SBTR_CUR_LOCAL=PEAK_SBTR_CUR_LOCAL+
     &        dble(MEM_SUBTREE(INDICE_SBTR))
         IF(.NOT.BDC_SBTR) INDICE_SBTR=INDICE_SBTR+1
      ELSE
         PEAK_SBTR_CUR_LOCAL=dble(0)
         SBTR_CUR_LOCAL=dble(0)
      ENDIF
      END SUBROUTINE ZMUMPS_513
      DOUBLE PRECISION FUNCTION ZMUMPS_543( INODE )
      IMPLICIT NONE
      INTEGER INODE,LEVEL,i,NELIM,NFR
      DOUBLE PRECISION COST
      EXTERNAL MUMPS_330
      INTEGER MUMPS_330
      i = INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
        NELIM = NELIM + 1
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      NFR = ND_LOAD( STEP_LOAD(INODE) ) + KEEP_LOAD(253)
      LEVEL = MUMPS_330( PROCNODE_LOAD(STEP_LOAD(INODE)), NPROCS )
      IF (LEVEL .EQ. 1) THEN
        COST =  dble(NFR) * dble(NFR)
      ELSE
        IF ( K50 == 0 ) THEN
           COST =  dble(NFR) * dble(NELIM)
        ELSE
           COST = dble(NELIM) * dble(NELIM)
        ENDIF
      ENDIF
      ZMUMPS_543=COST
      RETURN
      END FUNCTION ZMUMPS_543
      RECURSIVE SUBROUTINE ZMUMPS_515(FLAG,COST,COMM)
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER COMM,WHAT,IERR
      LOGICAL FLAG
      DOUBLE PRECISION COST
      DOUBLE PRECISION TO_BE_SENT
      EXTERNAL MUMPS_330
      INTEGER MUMPS_330
      IF(FLAG)THEN
         WHAT=17 
         IF(BDC_M2_FLOPS)THEN
#if ! defined(OLD_LOAD_MECHANISM)
            TO_BE_SENT=DELTA_LOAD-COST
            DELTA_LOAD=dble(0)
#else
            TO_BE_SENT=LAST_LOAD_SENT-COST
            LAST_LOAD_SENT=LAST_LOAD_SENT-COST
#endif
         ELSE IF(BDC_M2_MEM)THEN
            IF(BDC_POOL.AND.(.NOT.BDC_MD))THEN
               TO_BE_SENT=max(TMP_M2,POOL_LAST_COST_SENT)
               POOL_LAST_COST_SENT=TO_BE_SENT
            ELSE IF(BDC_MD)THEN
#if ! defined(OLD_LOAD_MECHANISM)
               DELTA_MEM=DELTA_MEM+TMP_M2
               TO_BE_SENT=DELTA_MEM
#else
               TO_BE_SENT=DM_LAST_MEM_SENT+TMP_M2
               DM_LAST_MEM_SENT=DM_LAST_MEM_SENT+TMP_M2
#endif
            ELSE
               TO_BE_SENT=dble(0)
            ENDIF
         ENDIF
      ELSE
         WHAT=6
         TO_BE_SENT=dble(0)
      ENDIF
 111  CONTINUE
      CALL ZMUMPS_460( WHAT,
     &         COMM, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     &               FUTURE_NIV2,
#endif
     &         COST, 
     &         TO_BE_SENT, 
     &         MYID, IERR  )
      IF ( IERR == -1 )THEN
         CALL ZMUMPS_467(COMM_LD, KEEP_LOAD)
         GOTO 111
      ELSE IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) "Internal Error in ZMUMPS_500",
     &   IERR
         CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_515
      SUBROUTINE ZMUMPS_512(INODE,STEP,NSTEPS,PROCNODE,FRERE,
     &     NE,COMM,SLAVEF,MYID,KEEP,KEEP8,N)
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER INODE,NSTEPS,MYID,SLAVEF,COMM,N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER FRERE(NSTEPS),NE(NSTEPS),STEP(N),PROCNODE(NSTEPS)
      EXTERNAL MUMPS_170,MUMPS_275
      LOGICAL MUMPS_170
      INTEGER i,NCB,NELIM
      INTEGER MUMPS_275
      INTEGER FATHER_NODE,FATHER,WHAT,IERR
      EXTERNAL MUMPS_330
      INTEGER MUMPS_330
      IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
         WRITE(*,*)MYID,': Problem in ZMUMPS_512'
         CALL MUMPS_ABORT()
      ENDIF
      IF((INODE.LT.0).OR.(INODE.GT.N)) THEN
         RETURN
      ENDIF
      i=INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
         NELIM = NELIM + 1
         i = FILS_LOAD(i)
         GOTO 10
      ENDIF
      NCB=ND_LOAD(STEP_LOAD(INODE))-NELIM + KEEP_LOAD(253)
      WHAT=5
      FATHER_NODE=DAD_LOAD(STEP_LOAD(INODE))
      IF (FATHER_NODE.EQ.0) THEN
         RETURN
      ENDIF
      IF((FRERE(STEP(FATHER_NODE)).EQ.0).AND.
     &     ((FATHER_NODE.EQ.KEEP(38)).OR.
     &     (FATHER_NODE.EQ.KEEP(20))))THEN
         RETURN 
      ENDIF
      IF(MUMPS_170(PROCNODE(STEP(FATHER_NODE)),
     &            SLAVEF)) THEN
         RETURN
      ENDIF
      FATHER=MUMPS_275(PROCNODE(STEP(FATHER_NODE)),SLAVEF)
      IF(FATHER.EQ.MYID)THEN
        IF(BDC_M2_MEM)THEN
           CALL ZMUMPS_816(FATHER_NODE)
        ELSEIF(BDC_M2_FLOPS)THEN
           CALL ZMUMPS_817(FATHER_NODE)
        ENDIF
        IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
           IF(MUMPS_330(PROCNODE_LOAD(STEP_LOAD(INODE)),
     &          NPROCS).EQ.1)THEN
              CB_COST_ID(POS_ID)=INODE
              CB_COST_ID(POS_ID+1)=1
              CB_COST_ID(POS_ID+2)=POS_MEM
              POS_ID=POS_ID+3
              CB_COST_MEM(POS_MEM)=int(MYID,8)
              POS_MEM=POS_MEM+1
              CB_COST_MEM(POS_MEM)=int(NCB,8)*int(NCB,8)
              POS_MEM=POS_MEM+1
           ENDIF
        ENDIF
        GOTO 666
      ENDIF
 111  CONTINUE      
      CALL ZMUMPS_519(WHAT, COMM, NPROCS,
     &     FATHER_NODE,INODE,NCB, KEEP(81),MYID,
     &     FATHER, IERR)
      IF (IERR == -1 ) THEN
        CALL ZMUMPS_467(COMM, KEEP)
        GOTO 111
      ELSE IF ( IERR .NE. 0 ) THEN
        WRITE(*,*) "Internal Error in ZMUMPS_512",
     &  IERR
        CALL MUMPS_ABORT()
      ENDIF 
 666  CONTINUE
      END SUBROUTINE ZMUMPS_512
      SUBROUTINE ZMUMPS_514(INODE,NUM_CALL)
      IMPLICIT NONE
      DOUBLE PRECISION MAXI
      INTEGER i,J,IND_MAXI
      INTEGER INODE,NUM_CALL
      IF(BDC_M2_MEM)THEN
         IF(((NUM_CALL.EQ.1).AND.(BDC_MD)).OR.
     &       ((NUM_CALL.EQ.2).AND.(.NOT.BDC_MD)))THEN
            RETURN
         ENDIF
      ENDIF
      IF((FRERE_LOAD(STEP_LOAD(INODE)).EQ.0).AND.
     &     ((INODE.EQ.KEEP_LOAD(38)).OR.
     &     (INODE.EQ.KEEP_LOAD(20)))) THEN
         RETURN
      ENDIF
      DO i=POOL_SIZE,1,-1
         IF(POOL_NIV2(i).EQ.INODE) GOTO 666
      ENDDO
         NB_SON(STEP_LOAD(INODE))=-1
      RETURN         
 666  CONTINUE
      IF(BDC_M2_MEM)THEN
         IF(POOL_NIV2_COST(i).EQ.MAX_M2)THEN
            TMP_M2=MAX_M2
            MAXI=dble(0)
            IND_MAXI=-9999
            DO J=POOL_SIZE,1,-1
               IF(J.NE.i) THEN
                  IF(POOL_NIV2_COST(J).GT.MAXI)THEN
                     MAXI=POOL_NIV2_COST(J)
                     IND_MAXI=J
                  ENDIF
               ENDIF
            ENDDO
            MAX_M2=MAXI
            J=IND_MAXI
            REMOVE_NODE_FLAG_MEM=.TRUE.
            REMOVE_NODE_COST_MEM=TMP_M2
            CALL ZMUMPS_515(REMOVE_NODE_FLAG,MAX_M2,COMM_LD)
            NIV2(MYID+1)=MAX_M2
         ENDIF
      ELSEIF(BDC_M2_FLOPS)THEN
         REMOVE_NODE_COST=POOL_NIV2_COST(i)
         REMOVE_NODE_FLAG=.TRUE.
         CALL ZMUMPS_515(REMOVE_NODE_FLAG,
     &        -POOL_NIV2_COST(i),COMM_LD)
         NIV2(MYID+1)=NIV2(MYID+1)-POOL_NIV2_COST(i)
      ENDIF
      DO J=i+1,POOL_SIZE
         POOL_NIV2(J-1)=POOL_NIV2(J)
         POOL_NIV2_COST(J-1)=POOL_NIV2_COST(J)
      ENDDO
      POOL_SIZE=POOL_SIZE-1
      END SUBROUTINE ZMUMPS_514
      RECURSIVE SUBROUTINE ZMUMPS_816(INODE)
      IMPLICIT NONE
      INTEGER INODE
      EXTERNAL MUMPS_330
      INTEGER MUMPS_330
      IF((INODE.EQ.KEEP_LOAD(20)).OR.
     &     (INODE.EQ.KEEP_LOAD(38)))THEN
         RETURN
      ENDIF
      IF(NB_SON(STEP_LOAD(INODE)).EQ.-1)THEN
         RETURN
      ELSE
         IF(NB_SON(STEP_LOAD(INODE)).LT.0)THEN
            WRITE(*,*)
     &        'Internal error 1 in ZMUMPS_816'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      NB_SON(STEP_LOAD(INODE))=
     &     NB_SON(STEP_LOAD(INODE))-1
      IF(NB_SON(STEP_LOAD(INODE)).EQ.0)THEN
         POOL_NIV2(POOL_SIZE+1)=INODE
         POOL_NIV2_COST(POOL_SIZE+1)=
     &        ZMUMPS_543(INODE)
         POOL_SIZE=POOL_SIZE+1
         IF(POOL_NIV2_COST(POOL_SIZE).GT.MAX_M2)THEN
            MAX_M2=POOL_NIV2_COST(POOL_SIZE)
            ID_MAX_M2=POOL_NIV2(POOL_SIZE)
            CALL ZMUMPS_515(REMOVE_NODE_FLAG_MEM,MAX_M2,COMM_LD)
            NIV2(1+MYID)=MAX_M2
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_816    
      RECURSIVE SUBROUTINE ZMUMPS_817(INODE)
      IMPLICIT NONE
      INTEGER INODE
      EXTERNAL MUMPS_330
      INTEGER MUMPS_330
      IF((INODE.EQ.KEEP_LOAD(20)).OR.
     &     (INODE.EQ.KEEP_LOAD(38)))THEN
         RETURN
      ENDIF
      IF(NB_SON(STEP_LOAD(INODE)).EQ.-1)THEN
         RETURN
      ELSE
         IF(NB_SON(STEP_LOAD(INODE)).LT.0)THEN
            WRITE(*,*)
     &        'Internal error 1 in ZMUMPS_817'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      NB_SON(STEP_LOAD(INODE))=
     &     NB_SON(STEP_LOAD(INODE))-1
      IF(NB_SON(STEP_LOAD(INODE)).EQ.0)THEN
         POOL_NIV2(POOL_SIZE+1)=INODE
         POOL_NIV2_COST(POOL_SIZE+1)=
     &        ZMUMPS_542(INODE)
         POOL_SIZE=POOL_SIZE+1
         MAX_M2=POOL_NIV2_COST(POOL_SIZE)
         ID_MAX_M2=POOL_NIV2(POOL_SIZE)
         CALL ZMUMPS_515(REMOVE_NODE_FLAG,
     &           POOL_NIV2_COST(POOL_SIZE),
     &        COMM_LD)
         NIV2(MYID+1)=POOL_NIV2_COST(POOL_SIZE)+NIV2(MYID+1)
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_817
      DOUBLE PRECISION FUNCTION ZMUMPS_542(INODE)
      INTEGER INODE
      INTEGER NFRONT,NELIM,i,LEVEL
      EXTERNAL MUMPS_330
      INTEGER MUMPS_330
      DOUBLE PRECISION COST
      i = INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
        NELIM = NELIM + 1
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      NFRONT = ND_LOAD( STEP_LOAD(INODE) ) + KEEP_LOAD(253)
      LEVEL = MUMPS_330( PROCNODE_LOAD(STEP_LOAD(INODE)), NPROCS )
      COST=dble(0)
      CALL MUMPS_511(NFRONT,NELIM,NELIM,
     &                          KEEP_LOAD(50),LEVEL,COST)
      ZMUMPS_542=COST
      RETURN
      END FUNCTION ZMUMPS_542
      INTEGER FUNCTION ZMUMPS_541( INODE )
      IMPLICIT NONE
      INTEGER INODE,NELIM,NFR,SON,IN,i
      INTEGER COST_CB
      COST_CB=0
      i = INODE
 10   CONTINUE
      IF ( i > 0 ) THEN
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      SON=-i
      DO i=1, NE_LOAD(STEP_LOAD(INODE))
         NFR = ND_LOAD( STEP_LOAD(SON) ) + KEEP_LOAD(253)
         IN=SON
         NELIM = 0
 20      CONTINUE
         IF ( IN > 0 ) THEN
            NELIM = NELIM + 1
            IN = FILS_LOAD(IN)
            GOTO 20
         ENDIF
         COST_CB=COST_CB+((NFR-NELIM)*(NFR-NELIM))
         SON=FRERE_LOAD(STEP_LOAD(SON))
      ENDDO
      ZMUMPS_541=COST_CB
      RETURN
      END FUNCTION ZMUMPS_541
      SUBROUTINE ZMUMPS_533(SLAVEF,NMB_OF_CAND,
     &     TAB_POS, NASS, KEEP,KEEP8, LIST_SLAVES,
     &     NSLAVES,INODE)
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, NASS, NSLAVES
      INTEGER, INTENT (IN) :: TAB_POS(SLAVEF+2)
      INTEGER, intent(in) :: NMB_OF_CAND
      INTEGER, INTENT (IN) :: LIST_SLAVES( NMB_OF_CAND )
      INTEGER KEEP(500),INODE
      INTEGER(8) KEEP8(150)
      INTEGER allocok
      DOUBLE PRECISION MEM_COST,FCT_COST
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE ::EMPTY_ARRAY,
     &     DELTA_MD,EMPTY_ARRAY2
      INTEGER NBROWS_SLAVE,i,WHAT,IERR,NPROCS_LOC
      LOGICAL FORCE_CAND
      MEM_COST=dble(0)
      FCT_COST=dble(0)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
        NPROCS_LOC=SLAVEF-1
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
        NPROCS_LOC=NMB_OF_CAND
      END IF
      IF(FORCE_CAND)THEN
         CALL ZMUMPS_540(INODE,FCT_COST,
     &        MEM_COST,NPROCS_LOC,NASS)
      ELSE
         CALL ZMUMPS_540(INODE,FCT_COST,
     &        MEM_COST,SLAVEF-1,NASS)
      ENDIF
      DO i=1,SLAVEF
         IDWLOAD(i)=i-1
      ENDDO
      ALLOCATE(EMPTY_ARRAY(NPROCS_LOC),DELTA_MD(NPROCS_LOC),
     &     EMPTY_ARRAY2(NPROCS_LOC),
     &    stat=allocok)
      DO i = 1, NSLAVES
         NBROWS_SLAVE = TAB_POS(i+1) - TAB_POS(i)
         DELTA_MD( i ) = FCT_COST - dble(NBROWS_SLAVE)*
     &        dble(NASS)
      END DO
      IF(FORCE_CAND)THEN
         DO i=NSLAVES+1,NPROCS_LOC
            DELTA_MD( i ) = FCT_COST
         ENDDO
      ELSE
         DO i=NSLAVES+1,SLAVEF-1
            DELTA_MD( i ) = FCT_COST
         ENDDO
      ENDIF
      WHAT=7
 111  CONTINUE
      CALL ZMUMPS_524(.FALSE., COMM_LD, MYID, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     &     FUTURE_NIV2,
#endif
     &     NPROCS_LOC, LIST_SLAVES,0,
     &     EMPTY_ARRAY,
     &     DELTA_MD,EMPTY_ARRAY2,WHAT, IERR)
      IF ( IERR == -1 ) THEN
          CALL ZMUMPS_467(COMM_LD, KEEP)
          GOTO 111
      ELSE IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) "Internal Error in ZMUMPS_533",
     &   IERR
         CALL MUMPS_ABORT()
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (FUTURE_NIV2(MYID+1) .NE. 0) THEN
#endif
        DO i = 1, NSLAVES
           MD_MEM(LIST_SLAVES(i))=MD_MEM(LIST_SLAVES(i))+
     &          int(DELTA_MD( i ),8)
#if ! defined(OLD_LOAD_MECHANISM)
           IF(FUTURE_NIV2(LIST_SLAVES(i)+1).EQ.0)THEN
              MD_MEM(LIST_SLAVES(i))=999999999_8
           ENDIF
#endif
        ENDDO
#if ! defined(OLD_LOAD_MECHANISM)
      ENDIF
#endif
      DEALLOCATE(EMPTY_ARRAY)
      DEALLOCATE(DELTA_MD)
      END SUBROUTINE ZMUMPS_533
      SUBROUTINE ZMUMPS_540(INODE,FCT_COST,
     &     MEM_COST,NSLAVES,NELIM)
      IMPLICIT NONE
      INTEGER INODE,NSLAVES,NFR,NELIM,IN
      DOUBLE PRECISION MEM_COST,FCT_COST
      NFR=ND_LOAD(STEP_LOAD(INODE)) + KEEP_LOAD(253)
      IN = INODE
      FCT_COST=dble(int(dble(NFR-NELIM)/dble(NSLAVES))+1)*
     &     dble(NELIM)
      MEM_COST=dble(int(dble(NFR-NELIM)/dble(NSLAVES))+1)*
     &     dble(NFR)
      END SUBROUTINE ZMUMPS_540
      SUBROUTINE ZMUMPS_819(INODE)
      IMPLICIT NONE
      INTEGER INODE
      INTEGER i,J,SON,NSLAVES_TEMP,POS_TEMP,K
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      IF((INODE.LT.0).OR.(INODE.GT.N_LOAD))THEN
         RETURN
      ENDIF
      IF(POS_ID.GT.1)THEN
         i=INODE
 10      CONTINUE
         IF ( i > 0 ) THEN
            i = FILS_LOAD(i)
            GOTO 10
         ENDIF
         SON=-i
         IF(POS_ID.LT.NE_LOAD(STEP_LOAD(INODE))*3)THEN
            i=1
         ENDIF
         DO i=1, NE_LOAD(STEP_LOAD(INODE))
            J=1
            DO WHILE (J.LT.POS_ID)
               IF(CB_COST_ID(J).EQ.SON)GOTO 295 
               J=J+3
            ENDDO
 295        CONTINUE
            IF(J.GE.POS_ID)THEN
               IF(MUMPS_275(
     &            PROCNODE_LOAD(STEP_LOAD(INODE)),NPROCS).EQ.MYID)THEN
                  IF(INODE.EQ.KEEP_LOAD(38))THEN
                     GOTO 666
#if ! defined(OLD_LOAD_MECHANISM)
                  ELSE
                     IF(FUTURE_NIV2(MYID+1).NE.0)THEN
                        WRITE(*,*)MYID,': i did not find ',SON
                        CALL MUMPS_ABORT()
                     ENDIF
                     GOTO 666
#endif
                  ENDIF
               ELSE
                  GOTO 666
               ENDIF
            ENDIF
            NSLAVES_TEMP=CB_COST_ID(J+1)       
            POS_TEMP=CB_COST_ID(J+2)
            DO K=J,POS_ID-1
               CB_COST_ID(K)=CB_COST_ID(K+3)
            ENDDO
            K=POS_TEMP
            DO WHILE (K.LE.POS_MEM-1)
               CB_COST_MEM(K)=CB_COST_MEM(K+2*NSLAVES_TEMP)
               K=K+1
            ENDDO
            POS_MEM=POS_MEM-2*NSLAVES_TEMP
            POS_ID=POS_ID-3
            IF((POS_MEM.LT.1).OR.(POS_ID.LT.1))THEN
               WRITE(*,*)MYID,': negative pos_mem or pos_id'
               CALL MUMPS_ABORT()
            ENDIF
 666        CONTINUE
            SON=FRERE_LOAD(STEP_LOAD(SON))
         ENDDO
      ENDIF
      END SUBROUTINE ZMUMPS_819
      SUBROUTINE ZMUMPS_820(FLAG)
      IMPLICIT NONE
      LOGICAL FLAG
      INTEGER i
      DOUBLE PRECISION MEM
      FLAG=.FALSE.
      DO i=0,NPROCS-1
         MEM=DM_MEM(i)+LU_USAGE(i)
         IF(BDC_SBTR)THEN
            MEM=MEM+SBTR_MEM(i)-SBTR_CUR(i)
         ENDIF
         IF((MEM/dble(TAB_MAXS(i))).GT.0.8d0)THEN
            FLAG=.TRUE.
            GOTO 666
         ENDIF
      ENDDO
 666  CONTINUE
      END SUBROUTINE ZMUMPS_820
      SUBROUTINE ZMUMPS_554(NBINSUBTREE,INSUBTREE,NBTOP,
     &           MIN_COST,SBTR)
      IMPLICIT NONE
      INTEGER NBINSUBTREE,INSUBTREE,NBTOP
      DOUBLE PRECISION MIN_COST
      LOGICAL SBTR
      INTEGER i
      DOUBLE PRECISION TMP_COST,TMP_MIN
      TMP_MIN=huge(TMP_MIN)
      DO i=0,NPROCS-1
         IF(i.NE.MYID)THEN
            IF(BDC_SBTR)THEN
               TMP_MIN=min(TMP_MIN,dble(TAB_MAXS(i))-(DM_MEM(i)+
     &              LU_USAGE(i))-(SBTR_MEM(i)-SBTR_CUR(i)))
            ELSE
               TMP_MIN=min(TMP_MIN,dble(TAB_MAXS(i))-
     &              (DM_MEM(i)+LU_USAGE(i)))
            ENDIF
         ENDIF
      ENDDO
      IF(NBINSUBTREE.GT.0)THEN
         IF(INSUBTREE.EQ.1)THEN
            TMP_COST=dble(TAB_MAXS(MYID))-(DM_MEM(MYID)+
     &           LU_USAGE(MYID))
     &           -(SBTR_MEM(MYID)-SBTR_CUR(MYID))
         ELSE
            SBTR=.FALSE.
            GOTO 777
         ENDIF
      ENDIF
      TMP_MIN=min(TMP_COST,TMP_MIN)
      IF(TMP_MIN.GT.MIN_COST) SBTR=.TRUE.
 777  CONTINUE
      END SUBROUTINE ZMUMPS_554
      SUBROUTINE ZMUMPS_818(INODE,MAX_MEM,PROC)
      IMPLICIT NONE
      INTEGER INODE,PROC
      INTEGER i,POS,NSLAVES,SLAVE,NCAND,J,NELIM,NCB,NFRONT,SON,K
      INTEGER allocok
      EXTERNAL MUMPS_330
      INTEGER  MUMPS_330
      DOUBLE PRECISION MAX_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MEM_ON_PROCS,
     &     RECV_BUF
      LOGICAL, DIMENSION(:), ALLOCATABLE :: CONCERNED
      DOUBLE PRECISION MAX_SENT_MSG
#if  defined(NOT_ATM_POOL_SPECIAL)
     DOUBLE PRECISION TMP
#endif
      IF((FRERE_LOAD(STEP_LOAD(INODE)).EQ.0)
     &           .AND.(INODE.EQ.KEEP_LOAD(38)))THEN
         RETURN
      ENDIF
#if  defined(NOT_ATM_POOL_SPECIAL)
      IF((INODE.LT.0).OR.(INODE.GT.N_LOAD))THEN
         MAX_MEM=huge(MAX_MEM)
         DO i=0,NPROCS-1
            TMP=dble(TAB_MAXS(i))-(DM_MEM(i)+LU_USAGE(i))
            IF(BDC_SBTR)THEN
               TMP=TMP-(SBTR_MEM(i)-SBTR_CUR(i))
            ENDIF
            MAX_MEM=min(MAX_MEM,TMP)
         ENDDO
         RETURN
      ENDIF
#endif
      ALLOCATE( MEM_ON_PROCS(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in ZMUMPS_818'
        CALL MUMPS_ABORT()
      ENDIF 
      ALLOCATE( CONCERNED(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in ZMUMPS_818'
        CALL MUMPS_ABORT()
      ENDIF 
      ALLOCATE( RECV_BUF(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in ZMUMPS_818'
        CALL MUMPS_ABORT()
      ENDIF 
      RECV_BUF=dble(0)
      MAX_SENT_MSG=dble(0)
      i = INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
        NELIM = NELIM + 1
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      SON=-i
      NFRONT=ND_LOAD(STEP_LOAD(INODE)) + KEEP_LOAD(253)
      NCB=NFRONT-NELIM
      IF(MUMPS_330(PROCNODE_LOAD(STEP_LOAD(INODE)),
     &     NPROCS).EQ.2)THEN
         NCAND=CAND_LOAD(NPROCS+1, STEP_TO_NIV2_LOAD(STEP_LOAD(INODE)))
      ENDIF
      DO i=0,NPROCS-1
         IF(i.EQ.MYID)THEN
            MEM_ON_PROCS(i)=dble(TAB_MAXS(i))-(DM_MEM(i)+
     &           LU_USAGE(i)+
     &           ZMUMPS_543(INODE))
            IF(BDC_SBTR)THEN
               MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-(SBTR_MEM(i)-SBTR_CUR(i))
            ENDIF
            CONCERNED(i)=.TRUE.
         ELSE
            MEM_ON_PROCS(i)=dble(TAB_MAXS(i))-(DM_MEM(i)+LU_USAGE(i))
            IF(BDC_SBTR)THEN
               MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-(SBTR_MEM(i)-SBTR_CUR(i))
            ENDIF
            IF(BDC_M2_MEM)THEN
               MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-NIV2(i+1)
            ENDIF
         ENDIF
         IF(MUMPS_330(PROCNODE_LOAD(STEP_LOAD(INODE)),
     &        NPROCS).EQ.2)THEN
            IF(BDC_MD.AND.(KEEP_LOAD(48).EQ.5))THEN
               DO J=1,NCAND
                  IF(CAND_LOAD(J, STEP_TO_NIV2_LOAD(STEP_LOAD(INODE)))
     &                 .EQ.i)THEN
                     MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-
     &                 ((dble(NFRONT)*dble(NCB))/dble(NCAND))
                     CONCERNED(i)=.TRUE.
                     GOTO 666
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
 666     CONTINUE
      ENDDO
      DO K=1, NE_LOAD(STEP_LOAD(INODE))
         i=1
         DO WHILE (i.LE.POS_ID)
            IF(CB_COST_ID(i).EQ.SON)GOTO 295
            i=i+3
         ENDDO
 295     CONTINUE
         IF(i.GE.POS_ID)THEN
#if ! defined(OLD_LOAD_MECHANISM)
            IF(FUTURE_NIV2(MYID+1).NE.0)THEN
               WRITE(*,*)MYID,': ',SON,'has not been found
     & in ZMUMPS_818'
               CALL MUMPS_ABORT()
            ENDIF
#endif
            GOTO 777
         ENDIF
         NSLAVES=CB_COST_ID(i+1)
         POS=CB_COST_ID(i+2)
         DO i=1,NSLAVES
            SLAVE=int(CB_COST_MEM(POS))
            IF(.NOT.CONCERNED(SLAVE))THEN
               MEM_ON_PROCS(SLAVE)=MEM_ON_PROCS(SLAVE)+
     &              dble(CB_COST_MEM(POS+1))
            ENDIF
            DO J=0,NPROCS-1
               IF(CONCERNED(J))THEN
                  IF(SLAVE.NE.J)THEN
                     RECV_BUF(J)=max(RECV_BUF(J),
     &                    dble(CB_COST_MEM(POS+1)))
                  ENDIF
               ENDIF
            ENDDO
            POS=POS+2
         ENDDO
 777     CONTINUE
         SON=FRERE_LOAD(STEP_LOAD(SON))
      ENDDO
      MAX_MEM=huge(MAX_MEM)
      WRITE(*,*)'NPROCS=',NPROCS,MAX_MEM
      DO i=0,NPROCS-1
         IF(MAX_MEM.GT.MEM_ON_PROCS(i))THEN
            PROC=i
         ENDIF
         MAX_MEM=min(MEM_ON_PROCS(i),MAX_MEM)
      ENDDO
      DEALLOCATE(MEM_ON_PROCS)
      DEALLOCATE(CONCERNED)
      DEALLOCATE(RECV_BUF)
      END SUBROUTINE ZMUMPS_818
      SUBROUTINE ZMUMPS_553(MIN_PROC,POOL,
     &                      LPOOL,INODE)
      IMPLICIT NONE
      INTEGER INODE,LPOOL,MIN_PROC
      INTEGER POOL(LPOOL)
      EXTERNAL MUMPS_275
      INTEGER MUMPS_275
      INTEGER i,NBTOP,INSUBTREE,NBINSUBTREE,NODE,FATHER,SON,J
      INTEGER SBTR_NB_LEAF,POS,K,allocok,L
      INTEGER, ALLOCATABLE, DIMENSION (:) ::  TMP_SBTR
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF((KEEP_LOAD(47).EQ.4).AND.
     &     ((NBINSUBTREE.NE.0)))THEN
         DO J=INDICE_SBTR,NB_SUBTREES
            NODE=MY_ROOT_SBTR(J) 
            FATHER=DAD_LOAD(STEP_LOAD(NODE))
            i=FATHER
 110        CONTINUE
            IF ( i > 0 ) THEN
               i = FILS_LOAD(i)
               GOTO 110
            ENDIF
            SON=-i
            i=SON
 120        CONTINUE
            IF ( i > 0 ) THEN
               IF(MUMPS_275(PROCNODE_LOAD(STEP_LOAD(i)),NPROCS).EQ.
     &              MIN_PROC)THEN
                  SBTR_NB_LEAF=MY_NB_LEAF(J)
                  POS=SBTR_FIRST_POS_IN_POOL(J)
                  IF(POOL(POS+SBTR_NB_LEAF).NE.MY_FIRST_LEAF(J))THEN
                     WRITE(*,*)MYID,': The first leaf is not ok'
                     CALL MUMPS_ABORT()
                  ENDIF
                  ALLOCATE (TMP_SBTR(SBTR_NB_LEAF), stat=allocok)
                  IF (allocok > 0 ) THEN
                     WRITE(*,*)MYID,': Not enough space
     &                    for allocation'
                     CALL MUMPS_ABORT()
                  ENDIF
                  POS=SBTR_FIRST_POS_IN_POOL(J)
                  DO K=1,SBTR_NB_LEAF
                     TMP_SBTR(K)=POOL(POS+K-1)
                  ENDDO
                  DO K=POS+1,NBINSUBTREE-SBTR_NB_LEAF
                     POOL(K)=POOL(K+SBTR_NB_LEAF)
                  ENDDO
                  POS=1
                  DO K=NBINSUBTREE-SBTR_NB_LEAF+1,NBINSUBTREE
                     POOL(K)=TMP_SBTR(POS)
                     POS=POS+1
                  ENDDO
                  DO K=INDICE_SBTR,J
                     SBTR_FIRST_POS_IN_POOL(K)=SBTR_FIRST_POS_IN_POOL(K)
     &                    -SBTR_FIRST_POS_IN_POOL(J)
                  ENDDO
                  SBTR_FIRST_POS_IN_POOL(J)=NBINSUBTREE-SBTR_NB_LEAF
                  POS=MY_FIRST_LEAF(J)
                  L=MY_NB_LEAF(J)
                  DO K=INDICE_SBTR,J
                     MY_FIRST_LEAF(J)=MY_FIRST_LEAF(J+1)
                     MY_NB_LEAF(J)=MY_NB_LEAF(J+1)
                  ENDDO
                  MY_FIRST_LEAF(INDICE_SBTR)=POS
                  MY_NB_LEAF(INDICE_SBTR)=L
                  INODE=POOL(NBINSUBTREE)
                  DEALLOCATE(TMP_SBTR)
                  RETURN
               ENDIF
               i = FRERE_LOAD(STEP_LOAD(i))
               GOTO 120
            ENDIF           
         ENDDO
      ENDIF
      DO J=NBTOP,1,-1
#if defined(NOT_ATM_POOL_SPECIAL)
         IF ( POOL(LPOOL-2-J) < 0 ) THEN
            NODE=-POOL(LPOOL-2-J)
         ELSE IF ( POOL(LPOOL-2-J) > N_LOAD ) THEN
            NODE = POOL(LPOOL-2-J) - N_LOAD
         ELSE
            NODE = POOL(LPOOL-2-J)
         ENDIF
#else
         NODE=POOL(LPOOL-2-J)
#endif
         FATHER=DAD_LOAD(STEP_LOAD(NODE))
         i=FATHER
 11      CONTINUE
         IF ( i > 0 ) THEN
            i = FILS_LOAD(i)
            GOTO 11
         ENDIF
         SON=-i
         i=SON
 12      CONTINUE
         IF ( i > 0 ) THEN
            IF(MUMPS_275(PROCNODE_LOAD(STEP_LOAD(i)),NPROCS).EQ.
     &           MIN_PROC)THEN
               INODE=NODE
               RETURN
            ENDIF
            i = FRERE_LOAD(STEP_LOAD(i))
            GOTO 12
         ENDIF
      ENDDO
      END SUBROUTINE ZMUMPS_553
      SUBROUTINE ZMUMPS_555(POOL, LPOOL,KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER LPOOL,POOL(LPOOL),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER i,POS
      EXTERNAL MUMPS_283
      LOGICAL MUMPS_283
      IF(.NOT.BDC_SBTR) RETURN
      POS=0
      DO i=NB_SUBTREES,1,-1
         DO WHILE(MUMPS_283(
     &            PROCNODE_LOAD(STEP_LOAD(POOL(POS+1))),
     &            NPROCS))
            POS=POS+1
         ENDDO
         SBTR_FIRST_POS_IN_POOL(i)=POS+1
         POS=POS+MY_NB_LEAF(i)
      ENDDO
      END SUBROUTINE ZMUMPS_555
      END MODULE ZMUMPS_LOAD 
