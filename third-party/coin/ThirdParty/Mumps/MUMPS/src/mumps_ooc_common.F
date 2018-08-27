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
      MODULE MUMPS_OOC_COMMON
      IMPLICIT NONE
      INTEGER, PARAMETER :: FCT = 0 
      INTEGER, PARAMETER, PUBLIC :: TYPEF_INVALID = -999999
      INTEGER, PUBLIC  ::  TYPEF_L, TYPEF_U, TYPEF_CB
      INTEGER OOC_NB_FILE_TYPE, OOC_FCT_TYPE
      INTEGER, DIMENSION(:,:),POINTER :: OOC_INODE_SEQUENCE
      INTEGER(8), DIMENSION(:,:),POINTER :: OOC_VADDR
      INTEGER,DIMENSION(:),POINTER:: KEEP_OOC
      INTEGER ICNTL1 
      INTEGER(8), DIMENSION(:),ALLOCATABLE :: AddVirtLibre
      LOGICAL,SAVE :: STRAT_IO_ASYNC,WITH_BUF,SOLVE
      INTEGER, DIMENSION(:),POINTER :: STEP_OOC,PROCNODE_OOC
      INTEGER, SAVE :: MYID_OOC,SLAVEF_OOC,LOW_LEVEL_STRAT_IO
      INTEGER(8), SAVE :: HBUF_SIZE, DIM_BUF_IO
      INTEGER ERR_STR_OOC_MAX_LEN
      PARAMETER(ERR_STR_OOC_MAX_LEN = 512)
      CHARACTER*1 ERR_STR_OOC(ERR_STR_OOC_MAX_LEN)
      INTEGER DIM_ERR_STR_OOC
      TYPE IO_BLOCK
        INTEGER :: INODE          
        LOGICAL :: MASTER         
        INTEGER :: Typenode       
        INTEGER :: NROW, NCOL, NFS 
        LOGICAL :: Last           
        INTEGER :: LastPiv        
        INTEGER :: LastPanelWritten_L 
        INTEGER :: LastPanelWritten_U 
        INTEGER,POINTER,DIMENSION(:) :: INDICES
      END TYPE
      PUBLIC IO_BLOCK
      INTEGER, PUBLIC :: STRAT_WRITE_MAX, STRAT_TRY_WRITE
      PARAMETER (STRAT_WRITE_MAX=1, STRAT_TRY_WRITE=2)
      END MODULE MUMPS_OOC_COMMON
      SUBROUTINE MUMPS_676(INT1,INT2,BIGINT)
      IMPLICIT NONE
      INTEGER INT1,INT2
      INTEGER(8) BIGINT
      INTEGER(8) TMP1,TMP2,CONV
      PARAMETER (CONV=1073741824_8) 
      TMP1=int(INT1,kind=kind(TMP1))
      TMP2=int(INT2,kind=kind(TMP2))
      BIGINT=(TMP1*CONV)+TMP2
      RETURN
      END SUBROUTINE MUMPS_676
      SUBROUTINE MUMPS_677(INT1,INT2,BIGINT)
      IMPLICIT NONE
      INTEGER INT1,INT2
      INTEGER(8) BIGINT
      INTEGER(8) TMP1,TMP2,CONV
      PARAMETER (CONV=1073741824_8) 
      TMP1=BIGINT/CONV
      TMP2=mod(BIGINT,CONV)
      INT1=int(TMP1)
      INT2=int(TMP2)
      RETURN
      END SUBROUTINE MUMPS_677
      SUBROUTINE MUMPS_796
     &           (TYPEF_L,TYPEF_U,TYPEF_CB,K201, K251, K50,
     &            TYPEF_INVALID)
      IMPLICIT NONE
      INTEGER, intent(out):: TYPEF_L, TYPEF_U, TYPEF_CB
      INTEGER, intent(in) :: K201, K251, K50
      INTEGER, intent(in) :: TYPEF_INVALID
      IF (K201 .EQ. 1 .AND. K50.EQ.0) THEN
        IF ( K251.NE.2 ) THEN
          TYPEF_L  = 1
          TYPEF_U  = 2
          TYPEF_CB = 3 
        ELSE
          TYPEF_U  = 1
          TYPEF_L  = TYPEF_INVALID 
          TYPEF_CB = 2 
        ENDIF
      ELSE 
         TYPEF_L = 1
         TYPEF_U = TYPEF_INVALID 
         TYPEF_CB=2 
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_796
      INTEGER FUNCTION MUMPS_808
     &                 (FWDORBWD, MTYPE, K201, K50)
      USE MUMPS_OOC_COMMON 
      INTEGER,     intent(in) :: MTYPE, K201, K50
      CHARACTER*1, intent(in) :: FWDORBWD
      IF ( (TYPEF_L .NE. 1 .AND. TYPEF_L .NE. TYPEF_INVALID)
     &     .OR. (TYPEF_U .NE. 1 .AND. TYPEF_U .NE. 2 .AND.
     &           TYPEF_U .NE. TYPEF_INVALID) ) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_808",
     &              TYPEF_L, TYPEF_U
        CALL MUMPS_ABORT()
      ENDIF
      IF (FWDORBWD .NE. 'F' .AND. FWDORBWD .NE. 'B') THEN
        WRITE(*,*) "Internal error in MUMPS_808,",FWDORBWD
        CALL MUMPS_ABORT()
      ENDIF
      IF (K201 .EQ. 1) THEN  
        IF (FWDORBWD .EQ. 'F') THEN
         IF((MTYPE.NE.1).AND.(K50.EQ.0))THEN
           MUMPS_808=TYPEF_U
         ELSE
           MUMPS_808=TYPEF_L
         ENDIF
        ELSE 
          IF(K50.EQ.0)THEN
            IF(MTYPE.NE.1)THEN
              MUMPS_808=TYPEF_L
            ELSE
              MUMPS_808=TYPEF_U
            ENDIF
          ELSE
            MUMPS_808=TYPEF_L
          ENDIF
        ENDIF
      ELSE  
        MUMPS_808 = 1
      ENDIF
      RETURN
      END FUNCTION MUMPS_808
