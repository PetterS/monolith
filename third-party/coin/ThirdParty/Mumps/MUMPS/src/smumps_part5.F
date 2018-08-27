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
      SUBROUTINE SMUMPS_26(id)
      USE SMUMPS_LOAD
      USE MUMPS_STATIC_MAPPING
      USE SMUMPS_STRUC_DEF
      USE TOOLS_COMMON
      USE SMUMPS_PARALLEL_ANALYSIS
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      TYPE(SMUMPS_STRUC), TARGET :: id
      INTEGER LIW, IKEEP, FILS, FRERE, PTRAR, NFSIZ
      INTEGER NE, NA
      INTEGER I, allocok
      INTEGER MAXIS1_CHECK
      INTEGER NB_NIV2, IDEST
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LOCAL_M, LOCAL_N
      INTEGER numroc
      EXTERNAL numroc
      INTEGER IRANK
      INTEGER MP, LP, MPG
      LOGICAL PROK, PROKG, LISTVAR_SCHUR_2BE_FREED
      INTEGER SIZE_SCHUR_PASSED
      INTEGER SBUF_SEND, SBUF_REC, TOTAL_MBYTES
      INTEGER(8) SBUF_RECOLD8, MIN_BUF_SIZE8
      INTEGER MIN_BUF_SIZE
      INTEGER(8) MAX_SIZE_FACTOR_TMP
      INTEGER LEAF, INODE, ISTEP, INN, LPTRAR
      INTEGER NBLEAF, NBROOT, MYROW_CHECK, INIV2
      INTEGER(8) K13TMP8, K14TMP8
      REAL PEAK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PAR2_NODES
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWtemp
      INTEGER, DIMENSION(:), ALLOCATABLE :: XNODEL, NODEL
      INTEGER, DIMENSION(:), POINTER :: SSARBR
      INTEGER, POINTER ::  NELT, LELTVAR
      INTEGER, DIMENSION(:), POINTER :: KEEP,INFO, INFOG
      INTEGER(8), DIMENSION(:), POINTER :: KEEP8
      INTEGER(8)                   :: ENTRIES_IN_FACTORS_LOC_MASTERS
      REAL, DIMENSION(:), POINTER :: RINFO
      REAL, DIMENSION(:), POINTER :: RINFOG
      INTEGER, DIMENSION(:), POINTER :: ICNTL
      LOGICAL I_AM_SLAVE, PERLU_ON, COND
      INTEGER :: OOC_STAT
      INTEGER MUMPS_330, MUMPS_275
      EXTERNAL MUMPS_330, MUMPS_275
      INTEGER K,J, IFS
      INTEGER SIZE_TEMP_MEM,SIZE_DEPTH_FIRST,SIZE_COST_TRAV
      LOGICAL IS_BUILD_LOAD_MEM_CALLED
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: TEMP_MEM
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_ROOT
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_LEAF
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_SIZE
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH_FIRST
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH_FIRST_SEQ
      INTEGER, DIMENSION (:), ALLOCATABLE :: SBTR_ID
      REAL, DIMENSION (:), ALLOCATABLE :: COST_TRAV_TMP
      INTEGER(8) :: TOTAL_BYTES
      INTEGER, POINTER, DIMENSION(:) ::  WORK1PTR, WORK2PTR,
     &     NFSIZPTR,
     &     FILSPTR,
     &     FREREPTR
      IS_BUILD_LOAD_MEM_CALLED=.FALSE.
      KEEP   => id%KEEP
      KEEP8  => id%KEEP8
      INFO   => id%INFO
      RINFO  => id%RINFO
      INFOG  => id%INFOG
      RINFOG => id%RINFOG
      ICNTL  => id%ICNTL
      NELT    => id%NELT
      LELTVAR => id%LELTVAR
      KEEP8(24) = 0_8  
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &     ( id%MYID .eq. MASTER .AND.
     &     id%KEEP(46) .eq. 1 ) )
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IF (PROK) WRITE( MP, 220 )
      IF (PROKG.AND.(MPG .NE. MP)) WRITE( MPG, 220 ) id%VERSION_NUMBER
 220  FORMAT( /' SMUMPS ',A ) 
      IF ( PROK ) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            WRITE(MP, '(A)') 'L U Solver for unsymmetric matrices'
         ELSE IF ( KEEP(50) .eq. 1 ) THEN
            WRITE(MP, '(A)') 
     & 'L D L^T Solver for symmetric positive definite matrices'
         ELSE
            WRITE(MP, '(A)') 
     &           'L D L^T Solver for general symmetric matrices'
         END IF
         IF ( KEEP(46) .eq. 1 ) THEN
            WRITE(MP, '(A)') 'Type of parallelism: Working host'
         ELSE
            WRITE(MP, '(A)') 'Type of parallelism: Host not working'
         END IF
      END IF
      IF ( PROKG .AND. (MP.NE.MPG)) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            WRITE(MPG, '(A)') 'L U Solver for unsymmetric matrices'
         ELSE IF ( KEEP(50) .eq. 1 ) THEN
            WRITE(MPG, '(A)') 
     & 'L D L^T Solver for symmetric positive definite matrices'
         ELSE
            WRITE(MPG, '(A)') 
     &           'L D L^T Solver for general symmetric matrices'
         END IF
         IF ( KEEP(46) .eq. 1 ) THEN
            WRITE(MPG, '(A)') 'Type of parallelism: Working host'
         ELSE
            WRITE(MPG, '(A)') 'Type of parallelism: Host not working'
         END IF
      END IF
      IF (PROK) WRITE( MP, 110 )
      IF (PROKG .AND. (MPG.NE.MP)) WRITE( MPG, 110 )
      CALL SMUMPS_647(id)
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( KEEP(60), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      IF (id%KEEP(60) .EQ. 2 .or. id%KEEP(60). EQ. 3) THEN
         CALL MPI_BCAST( id%NPROW, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%NPCOL, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%MBLOCK, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%NBLOCK, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
      ENDIF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( KEEP(54), 2, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(69), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(244), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(251), 3, MPI_INTEGER,MASTER,id%COMM,IERR)
      CALL MPI_BCAST( id%N, 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      IF ( KEEP(55) .EQ. 0) THEN
         IF ( KEEP(54) .eq. 3 ) THEN
          CALL MPI_ALLREDUCE( id%NZ_loc, id%NZ, 1, MPI_INTEGER, 
     &       MPI_SUM, id%COMM, IERR )
         ELSE
            CALL MPI_BCAST( id%NZ, 1, MPI_INTEGER, MASTER,
     &           id%COMM, IERR )
         END IF
      ELSE
         CALL MPI_BCAST( id%NA_ELT, 1, MPI_INTEGER, MASTER,
     &        id%COMM, IERR )
      ENDIF
      IF ( associated(id%MEM_DIST) ) deallocate( id%MEM_DIST )
      allocate( id%MEM_DIST( 0:id%NSLAVES-1 ), STAT=IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO(1) = -7
         INFO(2) = id%NSLAVES
         IF ( LP .GT. 0 ) THEN
            WRITE(LP, 150) 'MEM_DIST'
         END IF
      END IF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      id%MEM_DIST(0:id%NSLAVES-1) = 0
      CALL MUMPS_427(
     &     id%COMM,id%COMM_NODES,KEEP(69),KEEP(46),
     &     id%NSLAVES,id%MEM_DIST,INFO)
      CALL SMUMPS_658(id)
      IF (KEEP(244) .EQ. 1) THEN
         IF ( KEEP(54) .eq. 3 ) THEN
            CALL SMUMPS_664(id)
         END IF
         IF ( id%MYID .eq. MASTER ) THEN
 1234       CONTINUE
            IF ( ( (KEEP(23) .NE. 0) .AND.
     &           ( (KEEP(23).NE.7) .OR. KEEP(50).EQ. 2 ) )
     &           .OR.
     &           ( associated(id%A) .AND. KEEP(52) .EQ. 77 .AND.
     &           (KEEP(50).EQ.2))
     &        .OR.
     &           KEEP(52) .EQ. -2 ) THEN
               IF (.not.associated(id%A)) THEN
                  IF (KEEP(23).GT.2) KEEP(23) = 1
               ENDIF
               CALL SMUMPS_203(id%N, id%NZ, KEEP(23), id%IS1(1), id,
     &              ICNTL(1), INFO(1))
               IF (INFO(1) .LT. 0) THEN
                  KEEP(23) = 0
                  GOTO 10
               END IF
            END IF
            IF (KEEP(55) .EQ. 0) THEN
               IF ( KEEP(256) .EQ. 1 ) THEN 
                  LIW = 2 * id%NZ + 3 * id%N + 2
               ELSE
                  LIW = 2 * id%NZ + 3 * id%N + 2
               ENDIF
               IF (LIW.LT.3*id%N) LIW = 3*id%N
            ELSE
#if defined(metis) || defined(parmetis)
               COND = (KEEP(60) .NE. 0) .OR. (KEEP(256) .EQ. 5)
#else
               COND = (KEEP(60) .NE. 0)
#endif
               IF( COND ) THEN
                  LIW = id%N + id%N + 1
               ELSE
                  LIW =  id%N + id%N + id%N+3 + id%N+1
               ENDIF
            ENDIF
            IF (LIW.LT.3*id%N) LIW = 3*id%N
            IF (KEEP(23) .NE. 0) THEN
               IKEEP = id%N + 1
            ELSE
               IKEEP = 1
            END IF
            NA      = IKEEP +     id%N
            NE      = IKEEP + 2 * id%N
            FILS    = IKEEP + 3 * id%N
            FRERE   = FILS  +     id%N
            PTRAR   = FRERE +     id%N
            IF (KEEP(55) .EQ. 0) THEN
               NFSIZ  = PTRAR + 4 * id%N
               MAXIS1_CHECK = NFSIZ + id%N - 1
            ELSE
               NFSIZ  = PTRAR + 2 * (NELT + 1)
               MAXIS1_CHECK = NFSIZ + id%N -1
            ENDIF
            IF ( id%MAXIS1 .LT. MAXIS1_CHECK ) THEN
               IF (LP.GE.0) THEN
                  WRITE(LP,*) '***********************************'
                  WRITE(LP,*) 'MAXIS1 and MAXIS1_CHECK are different !!'
                  WRITE(LP,*) 'MAXIS1, MAXIS1_CHECK=',id%MAXIS1,
     &                 MAXIS1_CHECK
                  WRITE(LP,*) 'This might cause problems ...'
                  WRITE(LP,*) '***********************************'
               ENDIF
            END IF
            IF ( KEEP(256) .EQ. 1 ) THEN
               DO I = 1, id%N
                  id%IS1( IKEEP + I - 1 ) = id%PERM_IN( I )
               END DO
            END IF
            INFOG(1) = 0
            INFOG(2) = 0
            INFOG(8) = -1
            IF ( .NOT. associated( id%LISTVAR_SCHUR ) ) THEN
               SIZE_SCHUR_PASSED = 1
               LISTVAR_SCHUR_2BE_FREED=.TRUE.
               allocate( id%LISTVAR_SCHUR( 1 ), STAT=allocok )
               IF ( allocok .GT. 0 ) THEN
                  WRITE(*,*)
     &                 'PB allocating an array of size 1 in Schur '
                  CALL MUMPS_ABORT()
               END IF
            ELSE
               SIZE_SCHUR_PASSED=id%SIZE_SCHUR
               LISTVAR_SCHUR_2BE_FREED = .FALSE.
            END IF
            IF (KEEP(55) .EQ. 0) THEN
               CALL SMUMPS_195(id%N, id%NZ, id%IRN(1), id%JCN(1),
     &              LIW, id%IS1(IKEEP),
     &              id%IS1(PTRAR), KEEP(256), id%IS1(NFSIZ),
     &              id%IS1(FILS), id%IS1(FRERE),
     &              id%LISTVAR_SCHUR(1), SIZE_SCHUR_PASSED,
     &              ICNTL(1), INFOG(1), KEEP(1),KEEP8(1),id%NSLAVES, 
     &              id%IS1(1),id)
               IF ( (KEEP(23).LE.-1).AND.(KEEP(23).GE.-6) ) THEN
                  KEEP(23) = -KEEP(23)
                  IF (.NOT. associated(id%A)) KEEP(23) = 1
                  GOTO 1234
               ENDIF
               INFOG(7)     = KEEP(256)
            ELSE
               allocate( IWtemp ( 3*id%N ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = 3*id%N
                  IF ( LP .GT. 0 ) THEN
                     WRITE(LP, 150) 'IWtemp'
                  END IF
                  GOTO 10
               ENDIF
               allocate( XNODEL ( id%N+1 ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = id%N + 1
                  IF ( LP .GT. 0 ) THEN
                     WRITE(LP, 150) 'XNODEL'
                  END IF
                  GOTO 10
               ENDIF
               IF (LELTVAR.ne.id%ELTPTR(NELT+1)-1)  THEN
                  INFO(1) = -2002
                  INFO(2) = id%ELTPTR(NELT+1)-1
                  GOTO 10
               ENDIF
               allocate( NODEL ( LELTVAR ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = LELTVAR
                  IF ( LP .GT. 0 ) THEN
                     WRITE(LP, 150) 'NODEL'
                  END IF
                  GOTO 10
               ENDIF
               CALL SMUMPS_128(id%N, NELT,
     &              id%ELTPTR(1), id%ELTVAR(1), LIW,
     &              id%IS1(IKEEP),
     &              IWtemp(1), KEEP(256), id%IS1(NFSIZ), id%IS1(FILS),
     &              id%IS1(FRERE), id%LISTVAR_SCHUR(1),
     &              SIZE_SCHUR_PASSED, 
     &              ICNTL(1), INFOG(1), KEEP(1),KEEP8(1),
     &              id%ELTPROC(1), id%NSLAVES, 
     &              XNODEL(1), NODEL(1))
               DEALLOCATE(IWtemp)
               INFOG(7)=KEEP(256)
            ENDIF
            IF ( LISTVAR_SCHUR_2BE_FREED ) THEN
               deallocate( id%LISTVAR_SCHUR )
               NULLIFY   ( id%LISTVAR_SCHUR )
            ENDIF
            INFO(1)=INFOG(1)
            INFO(2)=INFOG(2)
            KEEP(28) = INFOG(6)
            IF ( INFO(1) .LT. 0 ) THEN
               GO TO 10
            ENDIF
         ENDIF
      ELSE
         IKEEP   = 1
         NA      = IKEEP +     id%N
         NE      = IKEEP + 2 * id%N
         FILS    = IKEEP + 3 * id%N
         FRERE   = FILS  +     id%N
         PTRAR   = FRERE +     id%N
         NFSIZ   = PTRAR + 4 * id%N
         IF(id%MYID .EQ. MASTER) THEN
            WORK1PTR => id%IS1(IKEEP : IKEEP + 3*id%N-1)
            WORK2PTR => id%IS1(PTRAR : PTRAR + 4*id%N-1)
            NFSIZPTR => id%IS1(NFSIZ : NFSIZ + id%N-1)
            FILSPTR  => id%IS1(FILS  : FILS  + id%N-1)
            FREREPTR => id%IS1(FRERE : FRERE + id%N-1)
         ELSE
            ALLOCATE(WORK1PTR(3*id%N))
            ALLOCATE(WORK2PTR(4*id%N))
         END IF
         CALL SMUMPS_715(id,
     &        WORK1PTR,
     &        WORK2PTR,
     &        NFSIZPTR,
     &        FILSPTR,
     &        FREREPTR)
         IF(id%MYID .EQ. 0) THEN
            NULLIFY(WORK1PTR, WORK2PTR, NFSIZPTR)
            NULLIFY(FILSPTR, FREREPTR)
         ELSE
            DEALLOCATE(WORK1PTR, WORK2PTR)
         END IF
         KEEP(28) = INFOG(6)
      END IF
 10   CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
      IF(id%MYID .EQ. MASTER) THEN
         CALL MUMPS_633(KEEP(12),ICNTL(14),
     &        KEEP(50),KEEP(54),ICNTL(6),KEEP(52))
         CALL SMUMPS_348(id%N, id%IS1(FILS), id%IS1(FRERE),
     &        id%IS1(IKEEP+2*id%N), id%IS1(IKEEP+id%N))
         IF (id%NSLAVES .EQ. 1) THEN
            id%NBSA = 0
            IF ( (id%KEEP(60).EQ.0).
     &           AND.(id%KEEP(53).EQ.0))  THEN 
               id%KEEP(20)=0
               id%KEEP(38)=0
            ENDIF
            id%KEEP(56)=0
            id%PROCNODE = 0
            IF (id%KEEP(60) .EQ. 2 .OR. id%KEEP(60).EQ.3) THEN
               CALL SMUMPS_564(id%KEEP(38), id%PROCNODE(1),
     &              1+2*id%NSLAVES, id%IS1(FILS),id%N)
            ENDIF
         ELSE
       PEAK = real(id%INFOG(5))*real(id%INFOG(5)) + 
     &        real(id%KEEP(2))*real(id%KEEP(2))     
            SSARBR => id%IS1(IKEEP:IKEEP+id%N-1)
            CALL SMUMPS_537(id%N,id%NSLAVES,ICNTL(1),
     &           INFOG(1),
     &           id%IS1(NE),
     &           id%IS1(NFSIZ),
     &           id%IS1(FRERE),
     &           id%IS1(FILS),
     &           KEEP(1),KEEP8(1),id%PROCNODE(1),
     &           SSARBR(1),id%NBSA,PEAK,IERR
     &           )
            NULLIFY(SSARBR)
            if(IERR.eq.-999) then 
               write(6,*) ' Internal error in MUMPS_369'
               INFO(1) = IERR
               GOTO 11
            ENDIF
            IF(IERR.NE.0) THEN 
               INFO(1) = -135
               INFO(2) = IERR
               GOTO 11
            ENDIF
            CALL SMUMPS_348(id%N, id%IS1(FILS),
     &           id%IS1(FRERE), id%IS1(IKEEP+2*id%N),
     &           id%IS1(IKEEP+id%N))
         ENDIF
 11      CONTINUE
      ENDIF
      CALL MUMPS_276( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
      CALL MPI_BCAST( id%NELT, 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      IF (KEEP(55) .EQ. 0) THEN
         if (associated(id%FRTPTR)) DEALLOCATE(id%FRTPTR)
         if (associated(id%FRTELT)) DEALLOCATE(id%FRTELT)
         allocate( id%FRTPTR(1), id%FRTELT(1) )
      ELSE
         LPTRAR = id%NELT+id%NELT+2
         CALL MUMPS_733(id%PTRAR, LPTRAR, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%PTRAR (Analysis)', ERRCODE=-7)
         CALL MUMPS_733(id%FRTPTR, id%N+1, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%FRTPTR (Analysis)', ERRCODE=-7)
         CALL MUMPS_733(id%FRTELT, id%NELT, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%FRTELT (Analysis)', ERRCODE=-7)
         CALL MUMPS_276( ICNTL(1), INFO(1), id%COMM, id%MYID )
         IF ( INFO(1) < 0 ) RETURN
         IF(id%MYID .EQ. MASTER) THEN
            CALL SMUMPS_153(
     &           id%N, NELT, id%ELTPTR(NELT+1)-1, id%IS1(FRERE),
     &           id%IS1(FILS),
     &           id%IS1(IKEEP+id%N), id%IS1(IKEEP+2*id%N), XNODEL, 
     &           NODEL, id%FRTPTR(1), id%FRTELT(1), id%ELTPROC(1))
            DO I=1, id%NELT+1
               id%PTRAR(id%NELT+I+1)=id%ELTPTR(I)
            ENDDO
            deallocate(XNODEL)
            deallocate(NODEL)
         END IF
         CALL MPI_BCAST( id%PTRAR(id%NELT+2), id%NELT+1, MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%FRTPTR(1), id%N+1, MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%FRTELT(1), id%NELT,  MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
      ENDIF
      IF(id%MYID .EQ. MASTER) THEN
         IF ( INFO( 1 ) .LT. 0 ) GOTO 12
         IF ( KEEP(55) .ne. 0 ) THEN
            CALL SMUMPS_120(id%N, NELT, id%ELTPROC(1),id%NSLAVES,
     &           id%PROCNODE(1))
         END IF
         NB_NIV2 = KEEP(56)
         IF ( NB_NIV2.GT.0 ) THEN
            allocate(PAR2_NODES(NB_NIV2),
     &           STAT=allocok)
            IF (allocok .GT.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'PAR2_NODES'
               END IF
               GOTO 12
            END IF
         ENDIF
         IF ((NB_NIV2.GT.0) .AND. (KEEP(24).EQ.0)) THEN
            INIV2 = 0
            DO 777 INODE = 1, id%N
               IF ( ( id%IS1(FRERE+INODE-1) .NE. id%N+1 ) .AND.
     &              ( MUMPS_330(id%PROCNODE(INODE),id%NSLAVES)
     &              .eq. 2) ) THEN
                  INIV2 = INIV2 + 1
                  PAR2_NODES(INIV2) = INODE
               END IF
 777        CONTINUE
            IF ( INIV2 .NE. NB_NIV2 ) THEN
               WRITE(*,*) "Internal Error 2 in SMUMPS_26",
     &              INIV2, NB_NIV2
               CALL MUMPS_ABORT()
            ENDIF
         ENDIF
         IF ( (KEEP(24) .NE. 0) .AND. (NB_NIV2.GT.0) ) THEN
            IF ( associated(id%CANDIDATES)) deallocate(id%CANDIDATES)
            allocate( id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     &           stat=allocok)
            if (allocok .gt.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2*(id%NSLAVES+1)
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'CANDIDATES'
               END IF
               GOTO 12
            END IF
            CALL MUMPS_393
     &           (PAR2_NODES,id%CANDIDATES,IERR)
            IF(IERR.NE.0)  THEN
               INFO(1) = -2002
               GOTO 12
            ENDIF
            CALL MUMPS_494()
            IF(IERR.NE.0)  THEN
               INFO(1) = -2002
               GOTO 12
            ENDIF
         ELSE
            IF (associated(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
            allocate(id%CANDIDATES(1,1), stat=allocok)
            IF (allocok .NE. 0) THEN
               INFO(1)= -7
               INFO(2)= 1
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'CANDIDATES'
               END IF
               GOTO 12
            ENDIF
         ENDIF
 12      CONTINUE
         KEEP(84) = ICNTL(27)   
      END IF
      CALL MUMPS_276( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
      CALL MPI_BCAST( id%KEEP(1), 110, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MUMPS_749( id%KEEP8(21), MASTER,
     &                     id%MYID, id%COMM, IERR)
      CALL MPI_BCAST( id%KEEP(205), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MPI_BCAST( id%NBSA, 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      IF (id%MYID==MASTER) KEEP(127)=INFOG(5)
      CALL MPI_BCAST( id%KEEP(127), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(226), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MUMPS_733(id%STEP, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%STEP (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_733(id%PROCNODE_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE.,
     &     STRING='id%PROCNODE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_733(id%NE_STEPS, id%KEEP(28), id%INFO, LP, 
     &     FORCE=.TRUE., 
     &     STRING='id%NE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_733(id%ND_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE., 
     &     STRING='id%ND_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_733(id%FRERE_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE., 
     &     STRING='id%FRERE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_733(id%DAD_STEPS, id%KEEP(28), id%INFO, LP, 
     &     FORCE=.TRUE., 
     &     STRING='id%DAD_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_733(id%FILS, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%FILS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_733(id%SYM_PERM, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%SYM_PERM (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      IF (KEEP(55) .EQ. 0) THEN
        LPTRAR = id%N+id%N
        CALL MUMPS_733(id%PTRAR, LPTRAR, id%INFO, LP, FORCE=.TRUE.,
     &       STRING='id%PTRAR (Analysis)', ERRCODE=-7)
        IF(INFO(1).LT.0) GOTO 94
      ENDIF
      IF ( associated( id%UNS_PERM ) ) deallocate(id%UNS_PERM)
      IF ( id%MYID == MASTER .AND. id%KEEP(23) .NE. 0 ) THEN
         allocate(id%UNS_PERM(id%N),stat=allocok)
         IF ( allocok .ne. 0) THEN
            INFO(1) = -7
            INFO(2) = id%N
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%UNS_PERM'
            END IF
            GOTO 94
         ENDIF
         DO I=1,id%N
            id%UNS_PERM(I) = id%IS1(I)
         END DO
      ENDIF
 94   CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%MYID .EQ. MASTER ) THEN
         DO I=1,id%N
            id%FILS(I) = id%IS1(FILS+I-1)
         ENDDO
      END IF
      IF (id%MYID .EQ. MASTER ) THEN
         IF (id%N.eq.1) THEN
            NBROOT = 1
            NBLEAF = 1
         ELSE IF (id%IS1(NA+id%N-1) .LT.0) THEN
            NBLEAF = id%N
            NBROOT = id%N
         ELSE IF (id%IS1(NA+id%N-2) .LT.0) THEN
            NBLEAF = id%N-1
            NBROOT = id%IS1(NA+id%N-1)
         ELSE
            NBLEAF = id%IS1(NA+id%N-2)
            NBROOT = id%IS1(NA+id%N-1)
         ENDIF
         id%LNA = 2+NBLEAF+NBROOT
      ENDIF
      CALL MPI_BCAST( id%LNA, 1, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MUMPS_733(id%NA, id%LNA, id%INFO, LP, FORCE=.TRUE., 
     &     STRING='id%NA (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 96
      IF (id%MYID .EQ.MASTER ) THEN
         id%NA(1) = NBLEAF
         id%NA(2) = NBROOT
         LEAF = 3
         IF ( id%N == 1 ) THEN
            id%NA(LEAF) = 1
            LEAF = LEAF + 1
         ELSE IF (id%IS1(NA+id%N-1) < 0) THEN
            id%NA(LEAF) = - id%IS1(NA+id%N-1)-1
            LEAF = LEAF + 1
            DO I = 1, NBLEAF - 1
               id%NA(LEAF) = id%IS1(NA+I-1)
               LEAF = LEAF + 1
            ENDDO
         ELSE IF (id%IS1(NA+id%N-2) < 0 ) THEN
            INODE = - id%IS1(NA+id%N-2) - 1
            id%NA(LEAF) = INODE
            LEAF =LEAF + 1
            IF ( NBLEAF > 1 ) THEN
               DO I = 1, NBLEAF - 1
                  id%NA(LEAF) = id%IS1(NA+I-1)
                  LEAF = LEAF + 1
               ENDDO
            ENDIF
         ELSE
            DO I = 1, NBLEAF
               id%NA(LEAF) = id%IS1(NA+I-1)
               LEAF = LEAF + 1
            ENDDO
         END IF
      END IF
 96   CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF (associated(id%Step2node))      THEN
        DEALLOCATE(id%Step2node)
        NULLIFY(id%Step2node)
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
         ISTEP = 0
         DO I = 1, id%N
            IF ( id%IS1(FRERE+I-1) .ne. id%N + 1 ) THEN
               ISTEP = ISTEP + 1
               id%STEP(I)=ISTEP
               INN = id%IS1(FILS+I-1)
               DO WHILE ( INN .GT. 0 )
                  id%STEP(INN) = - ISTEP
                  INN = id%IS1(FILS + INN -1)
               END DO
               IF (id%IS1(FRERE+I-1) .eq. 0) THEN
                  id%NA(LEAF) = I
                  LEAF = LEAF + 1
               ENDIF
            ENDIF
         END DO
         IF ( LEAF - 1 .NE. 2+NBROOT + NBLEAF ) THEN
            WRITE(*,*) 'Internal error 2 in SMUMPS_26'
            CALL MUMPS_ABORT()
         ENDIF
         IF ( ISTEP .NE. id%KEEP(28) ) THEN
            write(*,*) 'Internal error 3 in SMUMPS_26'
            CALL MUMPS_ABORT()
         ENDIF
         DO I = 1, id%N
            IF (id%IS1(FRERE+I-1) .NE. id%N+1) THEN
               id%PROCNODE_STEPS(id%STEP(I)) = id%PROCNODE( I )
               id%FRERE_STEPS(id%STEP(I))    = id%IS1(FRERE+I-1)
               id%NE_STEPS(id%STEP(I))    = id%IS1(NE+I-1)
               id%ND_STEPS(id%STEP(I))    = id%IS1(NFSIZ+I-1)
            ENDIF
         ENDDO
         DO I = 1, id%N
            IF ( id%STEP(I) .LE. 0) CYCLE
            IF (id%IS1(FRERE+I-1) .eq. 0) THEN
               id%DAD_STEPS(id%STEP(I)) = 0
            ENDIF
            IFS = id%IS1(FILS+I-1)
            DO WHILE ( IFS .GT. 0 )
               IFS= id%IS1(FILS + IFS -1)
            END DO
            IFS = -IFS
            DO WHILE (IFS.GT.0) 
               id%DAD_STEPS(id%STEP(IFS)) = I
               IFS   = id%IS1(FRERE+IFS-1)
            ENDDO
         END DO
         deallocate(id%PROCNODE)
         NULLIFY(id%PROCNODE)
         deallocate(id%IS1)
         NULLIFY(id%IS1)
               CALL SMUMPS_363(id%N, id%FRERE_STEPS(1),
     &              id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &              id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1), 
     &              id%KEEP(28), .TRUE., id%KEEP(28), id%KEEP(70),
     &              id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(215),
     &              id%KEEP(234), id%KEEP(55),
     &              id%PROCNODE_STEPS(1),id%NSLAVES,PEAK,id%KEEP(90)
     &              )
         IF ((id%KEEP(76).GE.4).OR.(id%KEEP(76).GE.6).OR.
     &              (id%KEEP(47).EQ.4).OR.((id%KEEP(81).GT.0)
     &              .AND.(id%KEEP(47).GE.2)))THEN
            IS_BUILD_LOAD_MEM_CALLED=.TRUE.
            IF ((id%KEEP(47) .EQ. 4).OR.
     &           (( id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
               IF(id%NSLAVES.GT.1) THEN
                  SIZE_TEMP_MEM = id%NBSA
               ELSE
                  SIZE_TEMP_MEM = id%NA(2)
               ENDIF
            ELSE
               SIZE_TEMP_MEM = 1
            ENDIF
            IF((id%KEEP(76).EQ.4).OR.(id%KEEP(76).EQ.6))THEN
               SIZE_DEPTH_FIRST=id%KEEP(28)
            ELSE
               SIZE_DEPTH_FIRST=1
            ENDIF
            allocate(TEMP_MEM(SIZE_TEMP_MEM,id%NSLAVES),STAT=allocok) 
            IF (allocok .NE.0) THEN
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'TEMP_MEM'
               END IF
               GOTO 80
            END IF
            allocate(TEMP_LEAF(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'TEMP_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(TEMP_SIZE(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'TEMP_SIZE'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(TEMP_ROOT(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'TEMP_ROOT'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(DEPTH_FIRST(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'DEPTH_FIRST'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            ALLOCATE(DEPTH_FIRST_SEQ(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'DEPTH_FIRST_SEQ'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            ALLOCATE(SBTR_ID(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'SBTR_ID'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            IF(id%KEEP(76).EQ.5)THEN
               SIZE_COST_TRAV=id%KEEP(28)
            ELSE
               SIZE_COST_TRAV=1
            ENDIF
            allocate(COST_TRAV_TMP(SIZE_COST_TRAV),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'COST_TRAV_TMP'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_COST_TRAV
               GOTO 80
            END IF
            IF(id%KEEP(76).EQ.5)THEN
               IF(id%KEEP(70).EQ.0)THEN
                  id%KEEP(70)=5
               ENDIF
               IF(id%KEEP(70).EQ.1)THEN
                  id%KEEP(70)=6
               ENDIF
            ENDIF
            IF(id%KEEP(76).EQ.4)THEN
               IF(id%KEEP(70).EQ.0)THEN
                  id%KEEP(70)=3
               ENDIF
               IF(id%KEEP(70).EQ.1)THEN
                  id%KEEP(70)=4
               ENDIF
            ENDIF
            CALL SMUMPS_364(id%N, id%FRERE_STEPS(1),
     &           id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &           id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1), 
     &           id%KEEP(28), .TRUE., id%KEEP(28), id%KEEP(70),
     &           id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(47),
     &           id%KEEP(81),id%KEEP(76),id%KEEP(215),
     &           id%KEEP(234), id%KEEP(55),
     &           id%PROCNODE_STEPS(1),TEMP_MEM,id%NSLAVES, 
     &           SIZE_TEMP_MEM, PEAK,id%KEEP(90),SIZE_DEPTH_FIRST,
     &           SIZE_COST_TRAV,DEPTH_FIRST(1),DEPTH_FIRST_SEQ(1),
     &           COST_TRAV_TMP(1),
     &           TEMP_LEAF,TEMP_SIZE,TEMP_ROOT,SBTR_ID(1)
     &              )
         END IF
         CALL SMUMPS_181(id%N, id%NA(1), id%LNA,
     &        id%NE_STEPS(1), id%SYM_PERM(1),
     &        id%FILS(1), id%DAD_STEPS(1),
     &        id%STEP(1), id%KEEP(28), id%INFO(1) )
      ENDIF
 80   CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      CALL MPI_BCAST( id%FILS(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%NA(1), id%LNA, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%STEP(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%PROCNODE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%DAD_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%FRERE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR)
      CALL MPI_BCAST( id%NE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%ND_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%SYM_PERM(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      IF (KEEP(55) .EQ. 0) THEN
         CALL SMUMPS_746(id, id%PTRAR(1))
         IF(id%MYID .EQ. MASTER) THEN
            IF ( (KEEP(244) .EQ. 1) .AND. (KEEP(54) .EQ. 3) ) THEN
               DEALLOCATE( id%IRN )
               DEALLOCATE( id%JCN )
            END IF
         END IF
      ENDIF
      IF((id%KEEP(76).EQ.4).OR.(id%KEEP(76).EQ.6))THEN
         IF(associated(id%DEPTH_FIRST))
     &        deallocate(id%DEPTH_FIRST)
         allocate(id%DEPTH_FIRST(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         IF(associated(id%DEPTH_FIRST_SEQ))
     *        DEALLOCATE(id%DEPTH_FIRST_SEQ)
         ALLOCATE(id%DEPTH_FIRST_SEQ(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(associated(id%SBTR_ID))
     *        DEALLOCATE(id%SBTR_ID)
         ALLOCATE(id%SBTR_ID(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%DEPTH_FIRST(1:id%KEEP(28))=DEPTH_FIRST(1:id%KEEP(28))
            id%DEPTH_FIRST_SEQ(1:id%KEEP(28))=
     &           DEPTH_FIRST_SEQ(1:id%KEEP(28))
            id%SBTR_ID(1:KEEP(28))=SBTR_ID(1:KEEP(28))
         ENDIF
         CALL MPI_BCAST( id%DEPTH_FIRST(1), id%KEEP(28), MPI_INTEGER,
     &           MASTER, id%COMM, IERR )         
         CALL MPI_BCAST( id%DEPTH_FIRST_SEQ(1), id%KEEP(28),
     &           MPI_INTEGER,MASTER, id%COMM, IERR )  
         CALL MPI_BCAST( id%SBTR_ID(1), id%KEEP(28),
     &           MPI_INTEGER,MASTER, id%COMM, IERR )  
      ELSE
         IF(associated(id%DEPTH_FIRST))
     &        deallocate(id%DEPTH_FIRST)
         allocate(id%DEPTH_FIRST(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         IF(associated(id%DEPTH_FIRST_SEQ))
     *        DEALLOCATE(id%DEPTH_FIRST_SEQ)
         ALLOCATE(id%DEPTH_FIRST_SEQ(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(associated(id%SBTR_ID))
     *        DEALLOCATE(id%SBTR_ID)
         ALLOCATE(id%SBTR_ID(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         id%SBTR_ID(1)=0
         id%DEPTH_FIRST(1)=0
         id%DEPTH_FIRST_SEQ(1)=0
      ENDIF
      IF(id%KEEP(76).EQ.5)THEN
         IF(associated(id%COST_TRAV))
     &        deallocate(id%COST_TRAV)
         allocate(id%COST_TRAV(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%COST_TRAV'
            END IF
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%COST_TRAV(1:id%KEEP(28))=
     &      dble(COST_TRAV_TMP(1:id%KEEP(28)))
         ENDIF
         CALL MPI_BCAST( id%COST_TRAV(1), id%KEEP(28),
     &        MPI_DOUBLE_PRECISION,MASTER, id%COMM, IERR )         
      ELSE
         IF(associated(id%COST_TRAV))
     &        deallocate(id%COST_TRAV)
         allocate(id%COST_TRAV(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%COST_TRAV(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         id%COST_TRAV(1)=0.0d0
      ENDIF
      IF (id%KEEP(47) .EQ. 4 .OR.
     &     ((id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
         IF(id%MYID .EQ. MASTER)THEN
            DO K=1,id%NSLAVES
               DO J=1,SIZE_TEMP_MEM
                  IF(TEMP_MEM(J,K) < 0.0D0) GOTO 666 
               ENDDO
 666           CONTINUE
               J=J-1
               IF (id%KEEP(46) == 1) THEN
                  IDEST = K - 1
               ELSE
                  IDEST = K
               ENDIF
               IF (IDEST .NE. MASTER) THEN
                  CALL MPI_SEND(J,1,MPI_INTEGER,IDEST,0,
     &                 id%COMM,IERR)
                  CALL MPI_SEND(TEMP_MEM(1,K),J,MPI_DOUBLE_PRECISION,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_LEAF(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_SIZE(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_ROOT(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)             
               ELSE
                  IF(associated(id%MEM_SUBTREE))
     &                 deallocate(id%MEM_SUBTREE)
                  allocate(id%MEM_SUBTREE(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LP .GT. 0 ) THEN
                        WRITE(LP, 150) 'id%MEM_SUBTREE'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%NBSA_LOCAL = J
                  id%MEM_SUBTREE(1:J)=TEMP_MEM(1:J,1)
                  IF(associated(id%MY_ROOT_SBTR))
     &                 deallocate(id%MY_ROOT_SBTR)
                  allocate(id%MY_ROOT_SBTR(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LP .GT. 0 ) THEN
                        WRITE(LP, 150) 'id%MY_ROOT_SBTR'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_ROOT_SBTR(1:J)=TEMP_ROOT(1:J,1)
                  IF(associated(id%MY_FIRST_LEAF))
     &                 deallocate(id%MY_FIRST_LEAF)
                  allocate(id%MY_FIRST_LEAF(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LP .GT. 0 ) THEN
                        WRITE(LP, 150) 'id%MY_FIRST_LEAF'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_FIRST_LEAF(1:J)=TEMP_LEAF(1:J,1)
                  IF(associated(id%MY_NB_LEAF))
     &                 deallocate(id%MY_NB_LEAF)
                  allocate(id%MY_NB_LEAF(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LP .GT. 0 ) THEN
                        WRITE(LP, 150) 'id%MY_NB_LEAF'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_NB_LEAF(1:J)=TEMP_SIZE(1:J,1)
               ENDIF
            ENDDO
         ELSE
            CALL MPI_RECV(id%NBSA_LOCAL,1,MPI_INTEGER,
     &           MASTER,0,id%COMM,STATUS, IERR)
            IF(associated(id%MEM_SUBTREE))
     &           deallocate(id%MEM_SUBTREE)
            allocate(id%MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'id%MEM_SUBTREE'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_ROOT_SBTR))
     &           deallocate(id%MY_ROOT_SBTR)
            allocate(id%MY_ROOT_SBTR(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'id%MY_ROOT_SBTR'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_FIRST_LEAF))
     &           deallocate(id%MY_FIRST_LEAF)
            allocate(id%MY_FIRST_LEAF(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'MY_FIRST_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_NB_LEAF))
     &           deallocate(id%MY_NB_LEAF)
            allocate(id%MY_NB_LEAF(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'MY_NB_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            CALL MPI_RECV(id%MEM_SUBTREE(1),id%NBSA_LOCAL,
     &           MPI_DOUBLE_PRECISION,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_FIRST_LEAF(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_NB_LEAF(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_ROOT_SBTR(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
         ENDIF
      ELSE
         id%NBSA_LOCAL = -999999
         IF(associated(id%MEM_SUBTREE))
     &        deallocate(id%MEM_SUBTREE)
         allocate(id%MEM_SUBTREE(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%MEM_SUBTREE(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_ROOT_SBTR))
     &        deallocate(id%MY_ROOT_SBTR)
         allocate(id%MY_ROOT_SBTR(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%MY_ROOT_SBTR(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_FIRST_LEAF))
     &        deallocate(id%MY_FIRST_LEAF)
         allocate(id%MY_FIRST_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%MY_FIRST_LEAF(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_NB_LEAF))
     &        deallocate(id%MY_NB_LEAF)
         allocate(id%MY_NB_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%MY_NB_LEAF(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
      ENDIF
      IF(id%MYID.EQ.MASTER)THEN
         IF(IS_BUILD_LOAD_MEM_CALLED)THEN 
            deallocate(TEMP_MEM)
            deallocate(TEMP_SIZE)
            deallocate(TEMP_ROOT)
            deallocate(TEMP_LEAF)
            deallocate(COST_TRAV_TMP)
            deallocate(DEPTH_FIRST)
            deallocate(DEPTH_FIRST_SEQ)
            deallocate(SBTR_ID)
         ENDIF
      ENDIF
 87   CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      NB_NIV2 = KEEP(56)        
      IF (  NB_NIV2.GT.0  ) THEN
         if (id%MYID.ne.MASTER) then
            IF (associated(id%CANDIDATES)) deallocate(id%CANDIDATES)
            allocate(PAR2_NODES(NB_NIV2),
     &           id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     &           STAT=allocok)
            IF (allocok .ne.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2*(id%NSLAVES+1)
               IF ( LP .GT. 0 ) THEN
                  WRITE(LP, 150) 'PAR2_NODES/id%CANDIDATES'
               END IF
            end if
         end if
         CALL MUMPS_276( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF ( INFO(1).LT.0 ) RETURN
         CALL MPI_BCAST(PAR2_NODES(1),NB_NIV2,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         IF (KEEP(24) .NE.0 ) THEN
            CALL MPI_BCAST(id%CANDIDATES(1,1),
     &           (NB_NIV2*(id%NSLAVES+1)),
     &           MPI_INTEGER, MASTER, id%COMM, IERR )
         ENDIF
      ENDIF
      IF ( associated(id%ISTEP_TO_INIV2)) THEN
         deallocate(id%ISTEP_TO_INIV2)
         NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
      IF ( associated(id%I_AM_CAND)) THEN
         deallocate(id%I_AM_CAND)
         NULLIFY(id%I_AM_CAND)
      ENDIF
      IF (NB_NIV2.EQ.0) THEN 
         id%KEEP(71) = 1
      ELSE
         id%KEEP(71) = id%KEEP(28)
      ENDIF
      allocate(id%ISTEP_TO_INIV2(id%KEEP(71)),
     &     id%I_AM_CAND(max(NB_NIV2,1)),
     &     stat=allocok)
      IF (allocok .gt.0) THEN
         IF ( LP .GT. 0 ) THEN
            WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
            WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
         END IF
         INFO(1)= -7
         IF (NB_NIV2.EQ.0) THEN
            INFO(2)= 2
         ELSE
            INFO(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
         END IF
         GOTO 321
      ENDIF
      IF ( NB_NIV2 .GT.0 ) THEN
         DO INIV2 = 1, NB_NIV2
            INN = PAR2_NODES(INIV2)
            id%ISTEP_TO_INIV2(id%STEP(INN)) = INIV2
         END DO 
         CALL SMUMPS_649( id%NSLAVES,
     &        NB_NIV2, id%MYID_NODES,
     &        id%CANDIDATES(1,1), id%I_AM_CAND(1) )
      ENDIF
#if                             ! defined(OLD_LOAD_MECHANISM)
      IF (associated(id%FUTURE_NIV2)) THEN
         deallocate(id%FUTURE_NIV2)
         NULLIFY(id%FUTURE_NIV2)
      ENDIF
      allocate(id%FUTURE_NIV2(id%NSLAVES), stat=allocok)
      IF (allocok .gt.0) THEN
         IF ( LP .GT. 0 ) THEN
            WRITE(LP, 150) 'FUTURE_NIV2'
         END IF
         INFO(1)= -7
         INFO(2)= id%NSLAVES
         GOTO 321
      ENDIF
      id%FUTURE_NIV2=0
      DO INIV2 = 1, NB_NIV2
         IDEST = MUMPS_275(
     &        id%PROCNODE_STEPS(id%STEP(PAR2_NODES(INIV2))),
     &        id%NSLAVES)
         id%FUTURE_NIV2(IDEST+1)=id%FUTURE_NIV2(IDEST+1)+1
      ENDDO
#endif
      IF ( I_AM_SLAVE ) THEN
         IF ( associated(id%TAB_POS_IN_PERE)) THEN
            deallocate(id%TAB_POS_IN_PERE)
            NULLIFY(id%TAB_POS_IN_PERE)
         ENDIF
         allocate(id%TAB_POS_IN_PERE(id%NSLAVES+2,max(NB_NIV2,1)),
     &        stat=allocok)
         IF (allocok .gt.0) THEN
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
               WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
            END IF
            INFO(1)= -7
            IF (NB_NIV2.EQ.0) THEN
               INFO(2)= 2
            ELSE
               INFO(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
            END IF
            GOTO 321
         ENDIF
      END IF
      IF (NB_NIV2.GT.0) deallocate (PAR2_NODES)
 321  CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( KEEP(23).NE.0 .and. id%MYID .EQ. MASTER ) THEN
         IKEEP = id%N + 1
      ELSE
         IKEEP = 1
      END IF
      FILS   = IKEEP + 3 * id%N
      NE     = IKEEP + 2 * id%N
      NA     = IKEEP +     id%N
      FRERE  = FILS  + id%N
      PTRAR  = FRERE + id%N
      IF (KEEP(55) .EQ. 0) THEN
         IF ( id%MYID.EQ.MASTER ) THEN
            NFSIZ   = PTRAR  + 4 * id%N
         ELSE
            NFSIZ   = PTRAR  + 2 * id%N
         ENDIF
      ELSE
         NFSIZ   = PTRAR  + 2 * (NELT + 1)
      END IF
      IF ( KEEP(38) .NE. 0 ) THEN
         CALL SMUMPS_164( id%MYID,
     &        id%NSLAVES, id%N, id%root,
     &        id%COMM_NODES, KEEP( 38 ), id%FILS(1),
     &        id%KEEP(50), id%KEEP(46),
     &        id%KEEP(51)
     &        , id%KEEP(60), id%NPROW, id%NPCOL, id%MBLOCK, id%NBLOCK
     &        )
      ELSE
         id%root%yes = .FALSE.
      END IF
      IF ( KEEP(38) .NE. 0 .and. I_AM_SLAVE ) THEN
         CALL MPI_ALLREDUCE(id%root%MYROW, MYROW_CHECK, 1,
     &        MPI_INTEGER, MPI_MAX, id%COMM_NODES, IERR)
         IF ( MYROW_CHECK .eq. -1) THEN
            INFO(1) = -25
            INFO(2) = 0
         END IF
         IF ( id%root%MYROW .LT. -1 .OR.
     &        id%root%MYCOL .LT. -1 ) THEN
            INFO(1) = -25
            INFO(2) = 0
         END IF
         IF ( LP > 0 .AND. INFO(1) == -25 ) THEN
            WRITE(LP, '(A)')
     &           'Problem with your version of the BLACS.'
            WRITE(LP, '(A)') 'Try using a BLACS version from netlib.'
         ENDIF
      END IF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( I_AM_SLAVE ) THEN
         IF (KEEP(55) .EQ. 0) THEN
            CALL SMUMPS_24( id%MYID,
     &           id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     &           id%STEP(1), id%PTRAR(1),
     &           id%PTRAR(id%N +1),
     &           id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &           KEEP(1),KEEP8(1), ICNTL(1), id )
         ELSE
            CALL SMUMPS_25( id%MYID,
     &           id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     &           id%STEP(1),
     &           id%PTRAR(1),
     &           id%PTRAR(id%NELT+2 ),
     &           id%NELT, 
     &           id%FRTPTR(1), id%FRTELT(1),
     &           KEEP(1), KEEP8(1), ICNTL(1), id%KEEP(50) )
         ENDIF
      ENDIF
      IF ( I_AM_SLAVE ) THEN
            IF ( id%root%yes ) THEN
               LOCAL_M = numroc( id%ND_STEPS(id%STEP(KEEP(38))),
     &              id%root%MBLOCK, id%root%MYROW, 0,
     &              id%root%NPROW )
               LOCAL_M = max(1, LOCAL_M)
               LOCAL_N = numroc( id%ND_STEPS(id%STEP(KEEP(38))),
     &              id%root%NBLOCK, id%root%MYCOL, 0,
     &              id%root%NPCOL )
            ELSE
               LOCAL_M = 0
               LOCAL_N = 0
            END IF
            IF  ( KEEP(60) .EQ. 2 .OR. KEEP(60) .EQ. 3 ) THEN
               id%SCHUR_MLOC=LOCAL_M
               id%SCHUR_NLOC=LOCAL_N
               id%root%SCHUR_MLOC=LOCAL_M
               id%root%SCHUR_NLOC=LOCAL_N
            ENDIF
               IF ( .NOT. associated(id%CANDIDATES)) THEN
                  ALLOCATE(id%CANDIDATES(id%NSLAVES+1,1))
               ENDIF
               CALL SMUMPS_246( id%MYID_NODES, id%N,
     &              id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     &              id%NA(1), id%LNA, id%NE_STEPS(1), id%DAD_STEPS(1),
     &              id%ND_STEPS(1), id%PROCNODE_STEPS(1),
     &              id%NSLAVES,
     &              KEEP8(11), KEEP(26), KEEP(15),
     &              KEEP8(12),  
     &              KEEP8(14),  
     &              KEEP(224), KEEP(225),
     &              KEEP(27), RINFO(1),
     &              KEEP(1), KEEP8(1), LOCAL_M, LOCAL_N, SBUF_RECOLD8,
     &              SBUF_SEND, SBUF_REC, id%COST_SUBTREES, KEEP(28),
     &              id%I_AM_CAND(1), max(KEEP(56),1),
     &              id%ISTEP_TO_INIV2(1), id%CANDIDATES(1,1), 
     &              INFO(1), INFO(2)
     &              ,KEEP8(15)
     &              ,MAX_SIZE_FACTOR_TMP, KEEP8(9) 
     &              ,ENTRIES_IN_FACTORS_LOC_MASTERS
     &           )
            id%MAX_SURF_MASTER = KEEP8(15)
            KEEP8(19)=MAX_SIZE_FACTOR_TMP
            KEEP( 29 ) = KEEP(15) + 2* max(KEEP(12),10)
     &           * ( KEEP(15) / 100 + 1)
            INFO( 19 ) = KEEP(225) + 2* max(KEEP(12),10)
     &           * ( KEEP(225) / 100 + 1)
            KEEP8(13)  = KEEP8(12) + int(KEEP(12),8) *
     &           ( KEEP8(12) / 100_8 + 1_8 )
            KEEP8(17)  = KEEP8(14) + int(KEEP(12),8) *
     &           ( KEEP8(14) /100_8 +1_8)
         CALL MUMPS_736 ( SBUF_RECOLD8, KEEP8(22), MPI_MAX,
     &                            id%COMM_NODES )
         SBUF_SEND = max(SBUF_SEND,KEEP(27))
         SBUF_REC  = max(SBUF_REC ,KEEP(27))
         CALL MPI_ALLREDUCE (SBUF_REC, KEEP(44), 1, 
     &        MPI_INTEGER, MPI_MAX,
     &        id%COMM_NODES, IERR)
         IF (KEEP(48)==5) THEN
            KEEP(43)=KEEP(44)
         ELSE
            KEEP(43)=SBUF_SEND
         ENDIF
         MIN_BUF_SIZE8 = KEEP8(22) / int(KEEP(238),8)
         MIN_BUF_SIZE8 = min( MIN_BUF_SIZE8, int(huge (KEEP(43)),8))
         MIN_BUF_SIZE  = int( MIN_BUF_SIZE8 )
         KEEP(44) = max(KEEP(44), MIN_BUF_SIZE)
         KEEP(43) = max(KEEP(43), MIN_BUF_SIZE)
            IF ( MP .GT. 0 ) THEN
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated INTEGER space for factors         :',
     &              KEEP(26)
               WRITE(MP,'(A,I10) ') 
     &              ' INFO(3), est. real space to store factors   :',
     &              KEEP8(11)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated number of entries in factors      :',
     &              KEEP8(9)
               WRITE(MP,'(A,I10) ') 
     &              ' Current value of space relaxation parameter :',
     &              KEEP(12)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated size of IS (In Core factorization):',
     &              KEEP(29)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated size of S  (In Core factorization):',
     &              KEEP8(13)
               WRITE(MP,'(A,I10) ') 
     &              ' Estimated size of S  (OOC factorization)    :',
     &              KEEP8(17)
            END IF
      ELSE
         ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
         KEEP8(13) = 0_8
         KEEP(29) = 0
         KEEP8(17)= 0_8
         INFO(19) = 0
         KEEP8(11) = 0_8
         KEEP(26) = 0
         KEEP(27) = 0
         RINFO(1) = 0.0E0
      END IF
      CALL MUMPS_736( ENTRIES_IN_FACTORS_LOC_MASTERS, 
     &     KEEP8(109), MPI_SUM, id%COMM)
      CALL MUMPS_736( KEEP8(19), KEEP8(119),
     &     MPI_MAX, id%COMM)
      CALL MPI_ALLREDUCE( KEEP(27), KEEP(127), 1,
     &     MPI_INTEGER, MPI_MAX,
     &     id%COMM, IERR)
      CALL MPI_ALLREDUCE( KEEP(26), KEEP(126), 1,
     &     MPI_INTEGER, MPI_SUM,
     &     id%COMM, IERR)
      CALL MUMPS_646( KEEP8(11), KEEP8(111), MPI_SUM,
     &     MASTER, id%COMM )
      CALL MUMPS_735( KEEP8(111), INFOG(3) )
      CALL MPI_ALLREDUCE( RINFO(1), RINFOG(1), 1,
     &     MPI_REAL, MPI_SUM,
     &     id%COMM, IERR)
      CALL MUMPS_735( KEEP8(11), INFO(3) )
      INFO ( 4 ) = KEEP(  26 )
      INFO ( 5 ) = KEEP(  27 )
      INFO ( 7 ) = KEEP(  29 )
      CALL MUMPS_735( KEEP8(13), INFO(8) )
      CALL MUMPS_735( KEEP8(17), INFO(20) )
      CALL MUMPS_735( KEEP8(9), INFO(24) )
      INFOG( 4 ) = KEEP( 126 )
      INFOG( 5 ) = KEEP( 127 )
      CALL MUMPS_735( KEEP8(109), INFOG(20) )
      CALL SMUMPS_100(id%MYID, id%COMM, KEEP(1), KEEP8(1),
     &     INFO(1), INFOG(1), RINFO(1), RINFOG(1), ICNTL(1))
         OOC_STAT = KEEP(201)
         IF (KEEP(201) .NE. -1) OOC_STAT=0 
         PERLU_ON = .FALSE.     
         CALL SMUMPS_214( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         KEEP8(2) = TOTAL_BYTES    
         PERLU_ON  = .TRUE.
         CALL SMUMPS_214( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         IF ( MP .gt. 0 ) THEN
            WRITE(MP,'(A,I10) ')
     & ' Estimated space in MBYTES for IC factorization            :',
     &           TOTAL_MBYTES
         END IF
         id%INFO(15) = TOTAL_MBYTES
      CALL MUMPS_243( id%MYID, id%COMM,
     &     id%INFO(15), id%INFOG(16), IRANK )
      IF ( PROKG ) THEN
         WRITE( MPG,'(A,I10) ')
     & ' ** Rank of proc needing largest memory in IC facto        :',
     &        IRANK
         WRITE( MPG,'(A,I10) ')
     & ' ** Estimated corresponding MBYTES for IC facto            :',
     &        id%INFOG(16)
         IF ( KEEP(46) .eq. 0 ) THEN
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (IC)     :'
     &           ,(id%INFOG(17)-id%INFO(15))/id%NSLAVES
         ELSE
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (IC)     :'
     &           ,id%INFOG(17)/id%NSLAVES
         END IF
         WRITE(MPG,'(A,I10) ')
     & ' ** TOTAL     space in MBYTES for IC factorization         :'
     &        ,id%INFOG(17)
      END IF
         OOC_STAT = KEEP(201)
#if defined(OLD_OOC_NOPANEL)
         IF (OOC_STAT .NE. -1) OOC_STAT=2
#else
         IF (OOC_STAT .NE. -1) OOC_STAT=1
#endif
         PERLU_ON = .FALSE.     
         CALL SMUMPS_214( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         KEEP8(3) = TOTAL_BYTES
         PERLU_ON  = .TRUE.     
         CALL SMUMPS_214( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     &        id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STAT, PERLU_ON, TOTAL_BYTES)
         id%INFO(17) = TOTAL_MBYTES
      CALL MUMPS_243( id%MYID, id%COMM,
     &     id%INFO(17), id%INFOG(26), IRANK )
      IF ( PROKG  ) THEN
         WRITE( MPG,'(A,I10) ')
     & ' ** Rank of proc needing largest memory for OOC facto      :',
     &        IRANK
         WRITE( MPG,'(A,I10) ')
     & ' ** Estimated corresponding MBYTES for OOC facto           :',
     &        id%INFOG(26)
         IF ( KEEP(46) .eq. 0 ) THEN
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (OOC)    :'
     &           ,(id%INFOG(27)-id%INFO(15))/id%NSLAVES
         ELSE
            WRITE( MPG,'(A,I10) ')
     & ' ** Estimated avg. MBYTES per work. proc at facto (OOC)    :'
     &           ,id%INFOG(27)/id%NSLAVES
         END IF
         WRITE(MPG,'(A,I10) ')
     & ' ** TOTAL     space in MBYTES for OOC factorization        :'
     &        ,id%INFOG(27)
      END IF
      IF ( id%MYID. eq. MASTER .AND. KEEP(54) .eq. 1 ) THEN
         IF (associated( id%MAPPING))
     &        deallocate( id%MAPPING)
         allocate( id%MAPPING(id%NZ), stat=allocok)
         IF ( allocok .GT. 0 ) THEN
            INFO(1) = -7
            INFO(2) = id%NZ
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%MAPPING'
            END IF
            GOTO 92
         END IF
         allocate(IWtemp( id%N ), stat=allocok)
         IF ( allocok .GT. 0 ) THEN
            INFO(1)=-7
            INFO(2)=id%N
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'IWtemp(N)'
            END IF
            GOTO 92
         END IF
         CALL SMUMPS_83(
     &        id%N, id%MAPPING(1),
     &        id%NZ, id%IRN(1),id%JCN(1), id%PROCNODE_STEPS(1),
     &        id%STEP(1),
     &        id%NSLAVES, id%SYM_PERM(1),
     &        id%FILS(1), IWtemp, id%KEEP(1),id%KEEP8(1),
     &        id%root%MBLOCK, id%root%NBLOCK,
     &        id%root%NPROW, id%root%NPCOL )
         deallocate( IWtemp )
 92      CONTINUE
      END IF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      RETURN
 110  FORMAT(/' ****** ANALYSIS STEP ********'/)
 150  FORMAT(
     & /' ** FAILURE DURING SMUMPS_26, DYNAMIC ALLOCATION OF',
     &     A30)
      END SUBROUTINE SMUMPS_26
      SUBROUTINE SMUMPS_537(N,NSLAVES,
     &     ICNTL,INFOG, NE, NFSIZ,
     &     FRERE, FILS,
     &     KEEP,KEEP8,PROCNODE,
     &     SSARBR,NBSA,PEAK,IERR
     &     )
      USE MUMPS_STATIC_MAPPING
      IMPLICIT NONE
      INTEGER N, NSLAVES, NBSA, IERR
      INTEGER ICNTL(40),INFOG(40),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NE(N),NFSIZ(N),FRERE(N),FILS(N),PROCNODE(N)
      INTEGER SSARBR(N)
      REAL PEAK
      CALL MUMPS_369(N,NSLAVES,
     &     ICNTL,INFOG, NE, NFSIZ,
     &     FRERE, FILS,
     &     KEEP,KEEP8,PROCNODE,
     &     SSARBR,NBSA,dble(PEAK),IERR
     &     )
      RETURN
      END SUBROUTINE SMUMPS_537
      SUBROUTINE SMUMPS_564(INODE, PROCNODE, VALUE, FILS, N)
      INTEGER, intent(in) :: INODE, N, VALUE
      INTEGER, intent(in) :: FILS(N)
      INTEGER, intent(inout) :: PROCNODE(N)
      INTEGER IN
      IN=INODE
      DO WHILE ( IN > 0 )
         PROCNODE( IN ) = VALUE
         IN=FILS( IN )
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_564
      SUBROUTINE SMUMPS_647(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(SMUMPS_STRUC)  :: id
      INTEGER   :: LP, MP, MPG, I
      INTEGER   :: MASTER
      LOGICAL   :: PROK, PROKG
      PARAMETER( MASTER = 0 )
      LP  = id%ICNTL( 1 )
      MP  = id%ICNTL( 2 )
      MPG = id%ICNTL( 3 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IF (id%MYID.eq.MASTER) THEN
        id%KEEP(256) = id%ICNTL(7) 
        id%KEEP(252) = id%ICNTL(32)
        IF (id%KEEP(252) < 0 .OR. id%KEEP(252) > 1 ) THEN
          id%KEEP(252) = 0
        ENDIF
        id%KEEP(251) = id%ICNTL(31)
        IF (id%KEEP(251) < 0 .OR. id%KEEP(251) > 2 ) THEN
          id%KEEP(251)=0
        ENDIF
        IF (id%KEEP(50) .EQ. 0 .AND. id%KEEP(252).EQ.1) THEN
          IF (id%KEEP(251) .NE. 1) id%KEEP(251) = 2
        ENDIF
        IF (id%KEEP(50) .NE.0 .AND. id%KEEP(251) .EQ. 2) THEN
          id%KEEP(251) = 0
        ENDIF
        IF (id%KEEP(251) .EQ. 1) THEN
          id%KEEP(201) = -1
        ENDIF
        IF (id%KEEP(252).EQ.1) THEN
          id%KEEP(253) = id%NRHS
          IF (id%KEEP(253) .LE. 0) THEN
            id%INFO(1)=-42
            id%INFO(2)=id%NRHS
            RETURN
          ENDIF
        ELSE
          id%KEEP(253) = 0
        ENDIF
      ENDIF
      IF ( (id%KEEP(24).NE.0) .AND.
     &     id%NSLAVES.eq.1 ) THEN
         id%KEEP(24) = 0
         IF ( PROKG ) THEN
            WRITE(MPG, '(A)')
     &           ' Resetting candidate strategy to 0 because NSLAVES=1'
            WRITE(MPG, '(A)') ' '
         END IF
      END IF
      IF ( (id%KEEP(24).EQ.0) .AND.
     &     id%NSLAVES.GT.1 ) THEN
         id%KEEP(24) = 8
      ENDIF
      IF ( (id%KEEP(24).NE.0)  .AND. (id%KEEP(24).NE.1)  .AND.
     &     (id%KEEP(24).NE.8)  .AND. (id%KEEP(24).NE.10) .AND.
     &     (id%KEEP(24).NE.12) .AND. (id%KEEP(24).NE.14) .AND.
     &     (id%KEEP(24).NE.16) .AND. (id%KEEP(24).NE.18)) THEN
         id%KEEP(24) = 8
         IF ( PROKG ) THEN
            WRITE(MPG, '(A)')
     &           ' Resetting candidate strategy to 8 '
            WRITE(MPG, '(A)') ' '
         END IF
      END IF
      id%KEEP8(21) = int(id%KEEP(85),8)
      IF ( id%MYID .EQ. MASTER ) THEN
        IF (id%KEEP(201).NE.-1) THEN
          id%KEEP(201)=id%ICNTL(22)
          IF (id%KEEP(201) .GT. 0) THEN
#if defined(OLD_OOC_NOPANEL)
            id%KEEP(201)=2
#else
            id%KEEP(201)=1
#endif
          ENDIF
        ENDIF
         id%KEEP(54) = id%ICNTL(18)
         IF ( id%KEEP(54) .LT. 0 .or. id%KEEP(54).GT.3 ) THEN
            IF ( PROKG ) THEN
               WRITE(MPG, *) ' Out-of-range value for id%ICNTL(18).'
               WRITE(MPG, *) ' Used 0 ie matrix not distributed'
            END IF
            id%KEEP(54) = 0
         END IF
         id%KEEP(55) = id%ICNTL(5)
         IF ( id%KEEP(55) .LT. 0 .OR. id%KEEP(55) .GT. 1 ) THEN
            IF ( PROKG ) THEN
               WRITE(MPG, *) ' Out-of-range value for id%ICNTL(5).'
               WRITE(MPG, *) ' Used 0 ie matrix is assembled'
            END IF
            id%KEEP(55) = 0
         END IF
         id%KEEP(60) = id%ICNTL(19)
         IF ( id%KEEP( 60 ) .LE. 0 ) id%KEEP( 60 ) = 0
         IF ( id%KEEP( 60 ) .GT. 3 ) id%KEEP( 60 ) = 0
         IF (id%KEEP(60) .NE. 0 .AND. id%SIZE_SCHUR == 0 ) THEN
            WRITE(MPG,'(A)')
     &           ' ** Schur option ignored because SIZE_SCHUR=0'
            id%KEEP(60)=0
         END IF
         IF ( id%KEEP(60) .NE.0 ) THEN
            id%KEEP(116) = id%SIZE_SCHUR
            IF (id%SIZE_SCHUR .LT. 0 .OR. id%SIZE_SCHUR .GE. id%N) THEN
              id%INFO(1)=-49
              id%INFO(2)=id%SIZE_SCHUR
              RETURN
            ENDIF
            IF ( .NOT. associated( id%LISTVAR_SCHUR ) ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 8
               RETURN
            ELSE IF (size(id%LISTVAR_SCHUR)<id%SIZE_SCHUR) THEN
               id%INFO(1) = -22
               id%INFO(2) = 8
               RETURN
            END IF
         ENDIF
         IF (id%KEEP(60) .EQ. 3 .AND. id%KEEP(50).NE.0) THEN
            IF (id%MBLOCK > 0 .AND. id%NBLOCK > 0 .AND.
     &           id%NPROW > 0 .AND. id%NPCOL > 0 ) THEN
               IF (id%NPROW *id%NPCOL .LE. id%NSLAVES) THEN
                  IF (id%MBLOCK .NE. id%NBLOCK ) THEN
                     id%INFO(1)=-31
                     id%INFO(2)=id%MBLOCK - id%NBLOCK
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         id%KEEP(244) = id%ICNTL(28)
         id%KEEP(245) = id%ICNTL(29)
#if ! defined(parmetis)         
         IF ((id%KEEP(244) .EQ. 2) .AND. (id%KEEP(245) .EQ. 2)) THEN
            id%INFO(1)  = -38
            IF(id%MYID .EQ.0 ) THEN
               WRITE(LP,'("ParMETIS not available.")')
               WRITE(LP,'("Aborting.")')
               RETURN
            END IF
         END IF
#endif
#if ! defined(ptscotch)         
         IF ((id%KEEP(244) .EQ. 2) .AND. (id%KEEP(245) .EQ. 1)) THEN
            id%INFO(1)  = -38
            IF(id%MYID .EQ.0 ) THEN
               WRITE(LP,'("PT-SCOTCH not available.")')
               WRITE(LP,'("Aborting.")')
               RETURN
            END IF
         END IF
#endif
         IF((id%KEEP(244) .GT. 2) .OR.
     &        (id%KEEP(244) .LT. 0)) id%KEEP(244)=0
         IF(id%KEEP(244) .EQ. 0) THEN
            id%KEEP(244) = 1
         ELSE IF (id%KEEP(244) .EQ. 2) THEN
            IF(id%KEEP(55) .NE. 0) THEN
               id%INFO(1)  = -39
               WRITE(LP,
     &              '("Incompatible values for ICNTL(5), ICNTL(28)")')
               WRITE(LP,
     &              '("Parallel analysis is not possible if the")')
               WRITE(LP,
     &              '("matrix is not assembled")')
               RETURN
            ELSE IF(id%KEEP(60) .NE. 0) THEN
               id%INFO(1)  = -39
               WRITE(LP,
     &              '("Incompatible values for ICNTL(19), ICNTL(28)")')
               WRITE(LP,
     &              '("Parallel analysis is not possible if SCHUR")')
               WRITE(LP,
     &              '("complement must be returned")')
               RETURN
            END IF
            IF(id%NSLAVES .LT. 2) THEN
               id%KEEP(244) = 1
               IF(PROKG) WRITE(MPG,
     &              '("Too few processes.
     & Reverting to sequential analysis")',advance='no')
               IF(id%KEEP(245) .EQ. 1) THEN
                  IF(PROKG) WRITE(MPG, '(" with SCOTCH")')
                  id%KEEP(256) = 3
               ELSE IF(id%KEEP(245) .EQ. 2) THEN
                  IF(PROKG) WRITE(MPG, '(" with Metis")')
                  id%KEEP(256) = 5
               ELSE
                  IF(PROKG) WRITE(MPG, '(".")')
                  id%KEEP(256) = 0
               END IF
            END IF
         END IF
         id%INFOG(32) = id%KEEP(244)
         IF ( (id%KEEP(244) .EQ. 1) .AND.
     &        (id%KEEP(256) .EQ. 1) ) THEN
            IF ( .NOT. associated( id%PERM_IN ) ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 3
               RETURN
            ELSE IF ( size( id%PERM_IN ) < id%N ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 3
               RETURN
            END IF
         ENDIF
         IF (id%KEEP(9) .LE. 1 ) id%KEEP(9) = 500
         IF ( id%KEEP8(21) .GT. 0_8 ) THEN 
            IF ((id%KEEP8(21).LE.1_8) .OR.
     &          (id%KEEP8(21).GT.int(id%KEEP(9),8)))
     &         id%KEEP8(21) = int(min(id%KEEP(9),100),8)
         ENDIF
         IF (id%KEEP(48). EQ. 1 ) id%KEEP(48) = -12345
         IF ( (id%KEEP(48).LT.0) .OR. (id%KEEP(48).GT.5) ) THEN
            id%KEEP(48)=5
         ENDIF
         IF ( (id%KEEP(60) .NE. 0) .AND. (id%KEEP(256) .EQ. 1) ) THEN
            DO I = 1, id%SIZE_SCHUR
               IF (id%PERM_IN(id%LISTVAR_SCHUR(I))
     &              .EQ. id%N-id%SIZE_SCHUR+I)
     &              CYCLE
               id%INFO(1) = -22
               id%INFO(2) = 8
               RETURN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Ignoring user-ordering, because incompatible with Schur.'
                  WRITE(MPG,'(A)') ' ** id%ICNTL(7) treated as 0.'
               END IF
               EXIT
            ENDDO
         END IF
         id%KEEP(95) = id%ICNTL(12)
         IF (id%KEEP(50).NE.2) id%KEEP(95) = 1
         IF ((id%KEEP(95).GT.3).OR.(id%KEEP(95).LT.0)) id%KEEP(95) = 0
         id%KEEP(23) = id%ICNTL(6)
         IF (id%KEEP(23).LT.0.OR.id%KEEP(23).GT.7) id%KEEP(23) = 7
         IF ( id%KEEP(50) .EQ. 1 ) THEN
            IF (id%KEEP(23) .NE. 0) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Max-trans not compatible with LLT factorization'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(95) .GT. 1) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** ICNTL(12) ignored: not compatible with LLT factorization'
               END IF
            ENDIF
            id%KEEP(95) = 1
         END IF
         IF  (id%KEEP(60) .GT. 0) THEN
            IF (id%KEEP(23) .NE. 0) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     &                 ' ** Max-trans not allowed because of Schur'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(52).NE.0) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Scaling during analysis not allowed because of Schur'
               ENDIF
               id%KEEP(52) = 0
            ENDIF
            IF (id%KEEP(95) .GT. 1) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** ICNTL(12) option not allowed because of Schur'
               END IF
            ENDIF
            id%KEEP(95) = 1
         END IF
         IF ( (id%KEEP(23) .NE. 0) .AND. (id%KEEP(256).EQ.1)) THEN
            id%KEEP(23) = 0
            id%KEEP(95) = 1
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** Max-trans not allowed because ordering is given'
            END IF
         END IF
         IF ( id%KEEP(256) .EQ. 1 ) THEN
            IF (id%KEEP(95) > 1 .AND. MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** ICNTL(12) option incompatible with given ordering'
            END IF
            id%KEEP(95) = 1
         END IF
         IF (id%KEEP(54) .NE. 0) THEN
            IF( id%KEEP(23) .NE. 0 ) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Max-trans not allowed because matrix is distributed'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(52).EQ.-2) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Scaling during analysis not allowed (matrix is distributed)'
               ENDIF
            ENDIF
            id%KEEP(52) = 0
            IF (id%KEEP(95) .GT. 1 .AND. MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** ICNTL(12) option not allowed because matrix is
     &distributed'
            ENDIF
            id%KEEP(95) = 1
         END IF
         IF ( id%KEEP(55) .NE. 0 ) THEN
            IF( id%KEEP(23) .NE. 0 ) THEN
               IF (MPG.GT.0) THEN
                  WRITE(MPG,'(A)')
     & ' ** Max-trans not allowed for element matrix'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (MPG.GT.0 .AND. id%KEEP(52).EQ.-2) THEN
               WRITE(MPG,'(A)')
     & ' ** Scaling not allowed at analysis for element matrix'
            ENDIF
            id%KEEP(52) = 0
            id%KEEP(95) = 1
         ENDIF
         IF(id%KEEP(244) .EQ. 2) THEN
            IF(id%KEEP(23) .EQ. 7) THEN
               id%KEEP(23) = 0
            ELSE IF (id%KEEP(23) .GT. 0) THEN
               id%INFO(1)  = -39
               id%KEEP(23) = 0
               WRITE(LP,
     &              '("Incompatible values for ICNTL(6), ICNTL(28)")')
               WRITE(LP,
     &              '("Maximum transversal not allowed
     &                 in parallel analysis")')
               RETURN
            END IF
         END IF
         IF ( id%KEEP(54) .NE. 0 .AND. id%KEEP(55) .NE. 0 ) THEN
            id%KEEP(54) = 0
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     & ' ** Distributed entry not available for element matrix'
            END IF
         ENDIF
         IF (id%ICNTL(39).NE.1 .and. id%ICNTL(39).NE.2) THEN
            id%KEEP(106)=1
         ELSE
            id%KEEP(106)=id%ICNTL(39)
         ENDIF
         IF(id%KEEP(50) .EQ. 2) THEN
            IF( .NOT. associated(id%A) ) THEN
               IF(id%KEEP(95) .EQ. 3) THEN
                  id%KEEP(95) = 2
               ENDIF
            ENDIF
            IF(id%KEEP(95) .EQ. 3 .AND. id%KEEP(256) .NE. 2) THEN
               IF (PROK) WRITE(MP,*)
     &              'WARNING: SMUMPS_203 constrained ordering not ', 
     &              'available with selected ordering'
               id%KEEP(95) = 2
            ENDIF 
            IF(id%KEEP(95) .EQ. 3) THEN
               id%KEEP(23) = 5
               id%KEEP(52) = -2
            ELSE IF(id%KEEP(95) .EQ. 2 .AND. 
     &              (id%KEEP(23) .EQ. 0 .OR. id%KEEP(23) .EQ. 7) ) THEN
               IF( associated(id%A) ) THEN
                  id%KEEP(23) = 5
               ELSE
                  id%KEEP(23) = 1
               ENDIF
            ELSE IF(id%KEEP(95) .EQ. 1) THEN
               id%KEEP(23) = 0
            ELSE IF(id%KEEP(95) .EQ. 0 .AND. id%KEEP(23) .EQ. 0) THEN
               id%KEEP(95) = 1
            ENDIF
         ELSE
            id%KEEP(95) = 1
         ENDIF
         id%KEEP(53)=0
         IF(id%KEEP(86).EQ.1)THEN
            IF(id%KEEP(47).LT.2) id%KEEP(47)=2
         ENDIF
         IF(id%KEEP(48).EQ.5)THEN
            IF(id%KEEP(50).EQ.0)THEN
               id%KEEP(87)=50
               id%KEEP(88)=50
            ELSE
               id%KEEP(87)=70
               id%KEEP(88)=70
            ENDIF
         ENDIF
         IF((id%NSLAVES.EQ.1).AND.(id%KEEP(76).GT.3))THEN
            id%KEEP(76)=2
         ENDIF
         IF(id%KEEP(81).GT.0)THEN
            IF(id%KEEP(47).LT.2) id%KEEP(47)=2
         ENDIF
      END IF
      RETURN
      END SUBROUTINE SMUMPS_647
      SUBROUTINE SMUMPS_664(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE(SMUMPS_STRUC)  :: id
      INTEGER, ALLOCATABLE :: REQPTR(:,:)
      INTEGER              :: MASTER, IERR, INDX, NRECV
      INTEGER              :: STATUS( MPI_STATUS_SIZE )
      INTEGER              :: LP, MP, MPG, I
      LOGICAL              :: PROK, PROKG
      PARAMETER( MASTER = 0 )
      LP  = id%ICNTL( 1 )
      MP  = id%ICNTL( 2 )
      MPG = id%ICNTL( 3 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IF ( id%KEEP(46) .EQ. 0 .AND. id%MYID .EQ. MASTER ) THEN
         id%NZ_loc = 0
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
         allocate( REQPTR( id%NPROCS, 3 ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) = 3 * id%NPROCS
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'REQPTR'
            END IF
            GOTO 13
         END IF
         allocate( id%IRN( id%NZ ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%NZ
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'IRN'
            END IF
            GOTO 13
         END IF
         allocate( id%JCN( id%NZ ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%NZ
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'JCN'
            END IF
            GOTO 13
         END IF
      END IF
 13   CONTINUE
      CALL MUMPS_276( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) < 0 ) RETURN
      IF ( id%MYID .EQ. MASTER ) THEN
         DO I = 1, id%NPROCS - 1
            CALL MPI_RECV( REQPTR( I+1, 1 ), 1, 
     &           MPI_INTEGER, I,
     &           COLLECT_NZ, id%COMM, STATUS, IERR )
         END DO
         IF ( id%KEEP(46) .eq. 0 ) THEN
            REQPTR( 1, 1 ) = 1
         ELSE
            REQPTR( 1, 1 ) = id%NZ_loc + 1
         END IF
         DO I = 2, id%NPROCS
            REQPTR( I, 1 ) = REQPTR( I, 1 ) + REQPTR( I-1, 1 )
         END DO
      ELSE
         CALL MPI_SEND( id%NZ_loc, 1, MPI_INTEGER, MASTER,
     &        COLLECT_NZ, id%COMM, IERR )
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
         NRECV = 0
         DO I = 1, id%NPROCS - 1
            IF ( REQPTR( I + 1, 1 ) - REQPTR( I, 1 ) .NE. 0 ) THEN
               NRECV = NRECV + 2
               CALL MPI_IRECV( id%IRN( REQPTR( I, 1 ) ),
     &              REQPTR( I + 1, 1 ) - REQPTR( I, 1 ), 
     &              MPI_INTEGER,
     &              I, COLLECT_IRN, id%COMM, REQPTR(I, 2), IERR )
               CALL MPI_IRECV( id%JCN( REQPTR( I, 1 ) ),
     &              REQPTR( I + 1, 1 ) - REQPTR( I, 1 ),   
     &              MPI_INTEGER,
     &              I, COLLECT_JCN, id%COMM, REQPTR(I, 3), IERR )
            ELSE
               REQPTR(I, 2) = MPI_REQUEST_NULL
               REQPTR(I, 3) = MPI_REQUEST_NULL
            END IF
         END DO
      ELSE
         IF ( id%NZ_loc .NE. 0 ) THEN
            CALL MPI_SEND( id%IRN_loc(1), id%NZ_loc, 
     &           MPI_INTEGER, MASTER,
     &           COLLECT_IRN, id%COMM, IERR )
            CALL MPI_SEND( id%JCN_loc(1), id%NZ_loc, 
     &           MPI_INTEGER, MASTER,
     &           COLLECT_JCN, id%COMM, IERR )
         END IF
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
         IF ( id%NZ_loc .NE. 0 ) THEN
            DO I=1,id%NZ_loc
               id%IRN(I) = id%IRN_loc(I)
               id%JCN(I) = id%JCN_loc(I)
            ENDDO
         END IF
         REQPTR( id%NPROCS, 2 ) = MPI_REQUEST_NULL
         REQPTR( id%NPROCS, 3 ) = MPI_REQUEST_NULL
         DO I = 1, NRECV
            CALL MPI_WAITANY
     &           ( 2 * id%NPROCS, REQPTR( 1, 2 ), INDX, STATUS, IERR )
         END DO
         deallocate( REQPTR )
      END IF
      RETURN
 150  FORMAT(
     &/' ** FAILURE DURING SMUMPS_664, DYNAMIC ALLOCATION OF',
     &     A30)
      END SUBROUTINE SMUMPS_664
      SUBROUTINE SMUMPS_658(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE(SMUMPS_STRUC)  :: id
      INTEGER              :: MASTER, IERR
      INTEGER              :: IUNIT
      LOGICAL              :: IS_ELEMENTAL
      LOGICAL              :: IS_DISTRIBUTED
      INTEGER              :: MM_WRITE
      INTEGER              :: MM_WRITE_CHECK
      CHARACTER(LEN=20)    :: MM_IDSTR
      LOGICAL              :: I_AM_SLAVE, I_AM_MASTER
      PARAMETER( MASTER = 0 )
      IUNIT = 69
      I_AM_SLAVE = ( id%MYID .NE. MASTER  .OR.
     &     ( id%MYID .EQ. MASTER .AND.
     &     id%KEEP(46) .EQ. 1 ) )
      I_AM_MASTER = (id%MYID.EQ.MASTER)
      IS_DISTRIBUTED = (id%KEEP(54) .EQ. 3)
      IS_ELEMENTAL   = (id%KEEP(55) .NE. 0)
      IF (id%MYID.EQ.MASTER .AND. .NOT. IS_DISTRIBUTED) THEN
        IF (id%WRITE_PROBLEM(1:20) .NE. "NAME_NOT_INITIALIZED")THEN
          OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM))
          CALL SMUMPS_166( id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &           IS_DISTRIBUTED,        
     &           IS_ELEMENTAL )         
          CLOSE(IUNIT)
        ENDIF
      ELSE IF (id%KEEP(54).EQ.3) THEN
         IF (id%WRITE_PROBLEM(1:20) .EQ. "NAME_NOT_INITIALIZED"
     &        .OR. .NOT. I_AM_SLAVE )THEN
            MM_WRITE = 0
         ELSE
            MM_WRITE = 1
         ENDIF
         CALL MPI_ALLREDUCE(MM_WRITE, MM_WRITE_CHECK, 1,
     &        MPI_INTEGER, MPI_SUM, id%COMM, IERR)
         IF (MM_WRITE_CHECK.EQ.id%NSLAVES .AND. I_AM_SLAVE) THEN
            WRITE(MM_IDSTR,'(I7)') id%MYID_NODES
            OPEN(IUNIT,
     &           FILE=trim(id%WRITE_PROBLEM)//trim(adjustl(MM_IDSTR)))
            CALL SMUMPS_166(id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &           IS_DISTRIBUTED,           
     &           IS_ELEMENTAL )            
            CLOSE(IUNIT)
         ENDIF
      ENDIF
      IF ( id%MYID.EQ.MASTER .AND.
     &     associated(id%RHS) .AND.
     &     id%WRITE_PROBLEM(1:20)
     &     .NE. "NAME_NOT_INITIALIZED")THEN
        OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM) //".rhs")
        CALL SMUMPS_179(IUNIT, id)
        CLOSE(IUNIT)
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_658
      SUBROUTINE SMUMPS_166
     & (id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &  IS_DISTRIBUTED, IS_ELEMENTAL )
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      LOGICAL, intent(in) :: I_AM_SLAVE,
     &                       I_AM_MASTER,
     &                       IS_DISTRIBUTED,
     &                       IS_ELEMENTAL
      INTEGER, intent(in) :: IUNIT
      TYPE(SMUMPS_STRUC), intent(in)  :: id
      CHARACTER (LEN=10)   :: SYMM
      CHARACTER (LEN=8)    :: ARITH
      INTEGER              :: I
      IF (IS_ELEMENTAL) THEN
        RETURN
      ENDIF
      IF (I_AM_MASTER .AND. .NOT. IS_DISTRIBUTED) THEN
         IF (associated(id%A)) THEN
               ARITH='real'
         ELSE
            ARITH='pattern '
         ENDIF
         IF (id%KEEP(50) .eq. 0) THEN
            SYMM="general"
         ELSE
            SYMM="symmetric"
         END IF
         WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix coordinate ',
     &           trim(ARITH)," ",trim(SYMM)
         WRITE(IUNIT,*) id%N, id%N, id%NZ
         IF (associated(id%A)) THEN
            DO I=1,id%NZ
               IF (id%KEEP(50).NE.0 .AND. id%IRN(I).LT.id%JCN(I)) THEN
                     WRITE(IUNIT,*) id%JCN(I), id%IRN(I), id%A(I)  
               ELSE
                     WRITE(IUNIT,*) id%IRN(I), id%JCN(I), id%A(I)  
               ENDIF
            ENDDO
         ELSE
            DO I=1,id%NZ
               IF (id%KEEP(50).NE.0 .AND. id%IRN(I).LT.id%JCN(I)) THEN
                  WRITE(IUNIT,*) id%JCN(I), id%IRN(I)
               ELSE
                     WRITE(IUNIT,*) id%IRN(I), id%JCN(I)
               ENDIF
            ENDDO
         ENDIF
      ELSE IF ( IS_DISTRIBUTED .AND. I_AM_SLAVE ) THEN
         IF (associated(id%A_loc)) THEN
               ARITH='real'
         ELSE
               ARITH='pattern '
         ENDIF
         IF (id%KEEP(50) .eq. 0) THEN
            SYMM="general"
         ELSE
            SYMM="symmetric"
         END IF
         WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix coordinate ',
     &           trim(ARITH)," ",trim(SYMM)
         WRITE(IUNIT,*) id%N, id%N, id%NZ_loc
         IF (associated(id%A_loc)) THEN
            DO I=1,id%NZ_loc
               IF (id%KEEP(50).NE.0 .AND.
     &             id%IRN_loc(I).LT.id%JCN_loc(I)) THEN
                     WRITE(IUNIT,*) id%JCN_loc(I), id%IRN_loc(I),
     &                    id%A_loc(I)
               ELSE
                     WRITE(IUNIT,*) id%IRN_loc(I), id%JCN_loc(I),
     &                    id%A_loc(I)
               ENDIF
            ENDDO
         ELSE
            DO I=1,id%NZ_loc
               IF (id%KEEP(50).NE.0 .AND. 
     &            id%IRN_loc(I).LT.id%JCN_loc(I)) THEN
                  WRITE(IUNIT,*) id%JCN_loc(I), id%IRN_loc(I)
               ELSE
                  WRITE(IUNIT,*) id%IRN_loc(I), id%JCN_loc(I)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_166
      SUBROUTINE SMUMPS_179(IUNIT, id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(SMUMPS_STRUC), intent(in)  :: id
      INTEGER, intent(in)             :: IUNIT
      CHARACTER (LEN=8)    :: ARITH
      INTEGER              :: I, J, K, LD_RHS
      IF (associated(id%RHS)) THEN
               ARITH='real'
        WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix array ',
     &           trim(ARITH),
     &           ' general'
        WRITE(IUNIT,*) id%N, id%NRHS
        IF ( id%NRHS .EQ. 1 ) THEN
           LD_RHS = id%N
        ELSE
           LD_RHS = id%LRHS
        ENDIF
        DO J = 1, id%NRHS
           DO I = 1, id%N
              K=(J-1)*LD_RHS+I
                 WRITE(IUNIT,*) id%RHS(K)
        ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_179
      SUBROUTINE SMUMPS_649( NSLAVES, NB_NIV2, MYID_NODES,
     &     CANDIDATES, I_AM_CAND )
      IMPLICIT NONE
      INTEGER, intent(in) :: NSLAVES, NB_NIV2, MYID_NODES
      INTEGER, intent(in) :: CANDIDATES( NSLAVES+1, NB_NIV2 )
      LOGICAL, intent(out):: I_AM_CAND( NB_NIV2 )
      INTEGER I, INIV2, NCAND
      DO INIV2=1, NB_NIV2
         I_AM_CAND(INIV2)=.FALSE.
         NCAND = CANDIDATES(NSLAVES+1,INIV2)
         DO I=1, NCAND
            IF (CANDIDATES(I,INIV2).EQ.MYID_NODES) THEN
               I_AM_CAND(INIV2)=.TRUE.
               EXIT
            ENDIF
         ENDDO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_649
      SUBROUTINE SMUMPS_251(N,IW,LIW,A,LA,
     &             NSTK_STEPS, NBPROCFILS,IFLAG,ND,FILS,STEP,
     &             FRERE, DAD, CAND,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &             MAXFRT, NTOTPV, NMAXNPIV, PTRIST, PTRAST,
     &             PIMASTER, PAMASTER, PTRARW,PTRAIW,
     &             ITLOC, RHS_MUMPS, IERROR,IPOOL, LPOOL,
     &             RINFO, POSFAC ,IWPOS, LRLU, IPTRLU,
     &             LRLUS, LEAF, NBROOT, NBRTOT,
     &             UU, ICNTL, PTLUST_S, PTRFAC, NSTEPS, INFO,
     &             KEEP,KEEP8,
     &             PROCNODE_STEPS,SLAVEF,MYID, COMM_NODES,
     &             MYID_NODES,
     &             BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,root,
     &             PERM, NELT, FRTPTR, FRTELT, LPTRAR,
     &             COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     &             MEM_DISTRIB, NE,
     &     DKEEP,PIVNUL_LIST,LPN_LIST)
      USE SMUMPS_LOAD
      USE SMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER N,IFLAG,NTOTPV,MAXFRT,LIW, LPTRAR, NMAXNPIV,
     &        IERROR, NSTEPS, INFO(40)
      INTEGER(8) :: LA
      REAL, TARGET :: A(LA)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER LPOOL
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER ITLOC(N+KEEP(253))
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER IW(LIW), NSTK_STEPS(KEEP(28)), NBPROCFILS(KEEP(28))
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND(KEEP(28))
      INTEGER FILS(N),PTRIST(KEEP(28))
      INTEGER STEP(N), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER PIMASTER(KEEP(28))
      INTEGER PTLUST_S(KEEP(28)), PERM(N)
      INTEGER CAND(SLAVEF+1,max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)),
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER IPOOL(LPOOL)
      INTEGER NE(KEEP(28))
      REAL RINFO(40)
      INTEGER(8) :: PAMASTER(KEEP(28)), PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: POSFAC, LRLU, LRLUS, IPTRLU
      INTEGER IWPOS, LEAF, NBROOT
      INTEGER COMM_LOAD, ASS_IRECV
      REAL UU, SEUIL, SEUIL_LDLT_NIV2
      INTEGER NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER        INTARR( max(1,KEEP(14)) )
      REAL DBLARR( max(1,KEEP(13)) )
      LOGICAL IS_ISOLATED_NODE
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      DOUBLE PRECISION, PARAMETER :: DZERO = 0.0D0, DONE = 1.0D0
      INTEGER INODE
      INTEGER IWPOSCB
      INTEGER FPERE, TYPEF
      INTEGER MP, LP, DUMMY(1)
      INTEGER NBFIN, NBRTOT, NBROOT_TRAITEES
      INTEGER NFRONT, IOLDPS
      INTEGER(8) NFRONT8
      INTEGER(8) :: POSELT
      INTEGER IPOSROOT, IPOSROOTROWINDICES
      INTEGER GLOBK109
      INTEGER(8) :: LBUFRX
      REAL, POINTER, DIMENSION(:) :: BUFRX
      LOGICAL :: IS_BUFRX_ALLOCATED
      DOUBLE PRECISION FLOP1
      INTEGER TYPE
      LOGICAL SON_LEVEL2, SET_IRECV, BLOCKING,
     &        MESSAGE_RECEIVED
      LOGICAL AVOID_DELAYED
      LOGICAL LAST_CALL
      INTEGER MASTER_ROOT
      INTEGER LOCAL_M, LOCAL_N
      INTEGER LRHS_CNTR_MASTER_ROOT, FWD_LOCAL_N_RHS
      LOGICAL ROOT_OWNER
      EXTERNAL MUMPS_330, MUMPS_275
      INTEGER MUMPS_330, MUMPS_275
      LOGICAL MUMPS_167,MUMPS_283
      EXTERNAL MUMPS_167,MUMPS_283
      LOGICAL SMUMPS_508
      EXTERNAL SMUMPS_508, SMUMPS_509
      LOGICAL STACK_RIGHT_AUTHORIZED
      INTEGER numroc
      EXTERNAL numroc
      INTEGER MAXFRW, NPVW, NOFFW, NELVAW, COMP,
     &        JOBASS, ETATASS
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, NextPiv2beWritten,
     &        IDUMMY
      INTEGER(8) :: ITMP8
      TYPE(IO_BLOCK) :: MonBloc
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION    OPASSW, OPELIW
       ASS_IRECV = MPI_REQUEST_NULL
       ITLOC(1:N+KEEP(253)) =0
       PTRIST  (1:KEEP(28))=0
       PTLUST_S(1:KEEP(28))=0
       PTRAST(1:KEEP(28))=0_8
       PTRFAC(1:KEEP(28))=-99999_8
      MP = ICNTL(2)
      LP = ICNTL(1)
      MAXFRW = 0
      NPVW   = 0
      NOFFW  = 0
      NELVAW = 0
      COMP  = 0
      OPASSW = DZERO
      OPELIW = DZERO
      IWPOSCB = LIW
      STACK_RIGHT_AUTHORIZED = .TRUE.
      CALL SMUMPS_22( .FALSE., 0_8,
     &     .FALSE., .FALSE., MYID_NODES, N, KEEP, KEEP8,
     &     IW, LIW, A, LA, LRLU, IPTRLU, IWPOS, IWPOSCB,
     &     PTRIST, PTRAST, STEP, PIMASTER,
     &     PAMASTER, KEEP(IXSZ), 0_8, -444, -444, .true.,
     &     COMP, LRLUS,
     &     IFLAG, IERROR
     &     )
      JOBASS  = 0
      ETATASS = 0
      NBFIN = NBRTOT
      NBROOT_TRAITEES = 0
      NBPROCFILS(1:KEEP(28)) = 0
      IF ( KEEP(38).NE.0 ) THEN
        IF (root%yes) THEN
            CALL SMUMPS_284(
     &        root, KEEP(38), N, IW, LIW,
     &        A, LA,
     &        FILS, MYID_NODES, PTRAIW, PTRARW,
     &        INTARR, DBLARR,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &        COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
        ENDIF
        IF ( IFLAG .LT. 0 ) GOTO 635
      END IF
 20   CONTINUE
      NIV1_FLAG=0
      SET_IRECV = .TRUE.
      BLOCKING = .FALSE.
      MESSAGE_RECEIVED = .FALSE.
      CALL SMUMPS_329(
     &      COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV,
     &      MESSAGE_RECEIVED,
     &      MPI_ANY_SOURCE, MPI_ANY_TAG,
     &      STATUS, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS,
     &      COMP, IFLAG,
     &      IERROR, COMM_NODES,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID_NODES, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS,
     &      PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &      STACK_RIGHT_AUTHORIZED )
      CALL SMUMPS_467(COMM_LOAD, KEEP)
      IF (MESSAGE_RECEIVED) THEN
          IF ( IFLAG .LT. 0 ) GO TO 640
          IF ( NBFIN .eq. 0 ) GOTO 640
      ELSE
          IF ( .NOT. SMUMPS_508( IPOOL, LPOOL) )THEN
            CALL SMUMPS_509( N, IPOOL, LPOOL,
     &      PROCNODE_STEPS,
     &      SLAVEF, STEP, INODE, KEEP,KEEP8, MYID_NODES, ND,
     &      (.NOT. STACK_RIGHT_AUTHORIZED) )
            STACK_RIGHT_AUTHORIZED = .TRUE.
            IF (KEEP(47) .GE. 3) THEN
              CALL SMUMPS_500(
     &              IPOOL, LPOOL,
     &              PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &              MYID_NODES, STEP, N, ND, FILS )
            ENDIF
            IF (KEEP(47).EQ.4) THEN
               IF(INODE.GT.0.AND.INODE.LE.N)THEN
                  IF((NE(STEP(INODE)).EQ.0).AND.
     &                 (FRERE(STEP(INODE)).EQ.0))THEN
                     IS_ISOLATED_NODE=.TRUE.
                  ELSE
                     IS_ISOLATED_NODE=.FALSE.
                  ENDIF
               ENDIF
               CALL SMUMPS_501(
     &              IS_ISOLATED_NODE,INODE,IPOOL,LPOOL,
     &              MYID_NODES,SLAVEF,COMM_LOAD,KEEP,KEEP8)
            ENDIF
            IF ((( KEEP(80) == 2 .OR. KEEP(80)==3 ) .AND.
     &           ( KEEP(47) == 4 )).OR.
     &           (KEEP(80) == 1 .AND. KEEP(47) .GE. 1)) THEN
               CALL SMUMPS_512(INODE,STEP,KEEP(28),
     &         PROCNODE_STEPS,FRERE,ND,COMM_LOAD,SLAVEF,
     &         MYID_NODES,KEEP,KEEP8,N)
            END IF
            GOTO 30
          ENDIF
      ENDIF
      GO TO 20
 30   CONTINUE
      IF ( INODE .LT. 0 ) THEN
        INODE = -INODE
        FPERE = DAD(STEP(INODE))
        GOTO 130
      ELSE IF (INODE.GT.N) THEN
       INODE = INODE - N
       IF (INODE.EQ.KEEP(38)) THEN
         NBROOT_TRAITEES = NBROOT_TRAITEES + 1
         IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
            NBFIN = NBFIN - NBROOT
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID_NODES,
     &          COMM_NODES, RACINE, SLAVEF)
            END IF
         ENDIF
         IF (NBFIN.EQ.0) GOTO 640
         GOTO 20
       ENDIF
       TYPE = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
       IF (TYPE.EQ.1) GOTO 100
       FPERE = DAD(STEP(INODE))
       AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     &                   .AND. KEEP(60).ne.0 )
       IF ( KEEP(50) .eq. 0 ) THEN
         CALL  SMUMPS_144( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, A, LA,
     &             UU, NOFFW,
     &             NPVW,
     &             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     &             NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC,
     &             STEP, PIMASTER, PAMASTER,
     &             NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     &             root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &              DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
          IF ( IFLAG .LT. 0 ) GOTO 640
       ELSE
         CALL  SMUMPS_141( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, A, LA,
     &             UU, NOFFW,
     &             NPVW,
     &             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     &             NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC,
     &             STEP, PIMASTER, PAMASTER,
     &             NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     &             root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL_LDLT_NIV2,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &              DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
        IF ( IFLAG .LT. 0 ) GOTO 640
        IF ( IW( PTLUST_S(STEP(INODE)) + KEEP(IXSZ) + 5 ) .GT. 1 ) THEN
             GOTO 20
        END IF
       END IF
       GOTO 130
      ENDIF
      IF (INODE.EQ.KEEP(38)) THEN
         CALL  SMUMPS_176( COMM_LOAD, ASS_IRECV,
     &    root, FRERE,
     &    INODE,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS, COMP,
     &    IFLAG, IERROR, COMM_NODES,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID_NODES, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
         IF ( IFLAG .LT. 0 ) GOTO 640
        GOTO 20
      ENDIF
      TYPE = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF (TYPE.EQ.1) THEN
        IF (KEEP(55).NE.0) THEN
         CALL SMUMPS_36( COMM_LOAD, ASS_IRECV,
     &        NELT, FRTPTR, FRTELT,
     &        N,INODE,IW,LIW,A,LA,
     &        IFLAG,IERROR,ND,
     &        FILS,FRERE,DAD,MAXFRW,root,OPASSW, OPELIW,
     &     PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS,
     &        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     &    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM, ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        ELSE
         JOBASS = 0
         CALL SMUMPS_252(COMM_LOAD, ASS_IRECV,
     &        N,INODE,IW,LIW,A,LA,
     &        IFLAG,IERROR,ND,
     &        FILS,FRERE,DAD,MAXFRW,root,OPASSW, OPELIW,
     &      PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS,
     &        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     &    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    JOBASS,ETATASS )
        ENDIF
       IF ( IFLAG .LT. 0 ) GOTO 640
        IF ((NBPROCFILS(STEP(INODE)).GT.0).OR.(SON_LEVEL2)) GOTO 20
      ELSE
        IF ( KEEP(55) .eq. 0 ) THEN
          CALL SMUMPS_253(COMM_LOAD, ASS_IRECV,
     &    N, INODE, IW, LIW, A, LA,
     &    IFLAG, IERROR,
     &    ND, FILS, FRERE, DAD, CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW,
     &    root, OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP,KEEP8,INTARR,DBLARR,
     &    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     &    MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     &    MEM_DISTRIB(0)
     &    )
        ELSE
          CALL SMUMPS_37( COMM_LOAD, ASS_IRECV,
     &    NELT, FRTPTR, FRTELT,
     &    N, INODE, IW, LIW, A, LA, IFLAG, IERROR,
     &    ND, FILS, FRERE, DAD, CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW,
     &    root, OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP,KEEP8,INTARR,DBLARR,
     &    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     &    MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     &    MEM_DISTRIB(0))
        END IF
        IF (IFLAG.LT.0) GOTO 640
        GOTO 20
      ENDIF
 100  CONTINUE
       FPERE = DAD(STEP(INODE))
      IF ( INODE .eq. KEEP(20) ) THEN
        POSELT = PTRAST(STEP(INODE))
        IF (PTRFAC(STEP(INODE)).NE.POSELT) THEN
          WRITE(*,*) "ERROR 2 in SMUMPS_251", POSELT
          CALL MUMPS_ABORT()
        ENDIF
        CALL SMUMPS_87
     &       ( IW(PTLUST_S(STEP(INODE))+KEEP(IXSZ)), KEEP(253) )
        GOTO 200
      END IF
      POSELT = PTRAST(STEP(INODE))
      IOLDPS = PTLUST_S(STEP(INODE))
      AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     &                   .AND. KEEP(60).ne.0 )
      IF (KEEP(50).EQ.0) THEN
         CALL SMUMPS_143( N, INODE, IW, LIW, A, LA,
     &               IOLDPS, POSELT,
     &               IFLAG, UU, NOFFW, NPVW,
     &               KEEP,KEEP8,
     &               STEP, PROCNODE_STEPS, MYID_NODES, SLAVEF,
     &               SEUIL, AVOID_DELAYED, ETATASS,
     &              DKEEP(1),PIVNUL_LIST(1),LPN_LIST, IWPOS)
         JOBASS = ETATASS
         IF (JOBASS.EQ.1) THEN
           CALL SMUMPS_252(COMM_LOAD, ASS_IRECV,
     &        N,INODE,IW,LIW,A,LA,
     &        IFLAG,IERROR,ND,
     &        FILS,FRERE,DAD,MAXFRW,root,OPASSW, OPELIW,
     &      PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP,PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS,
     &        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     &    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    JOBASS,ETATASS )
         ENDIF
      ELSE
         IW( IOLDPS+4+KEEP(IXSZ) ) = 1
         CALL SMUMPS_140( N, INODE,
     &           IW, LIW, A, LA,
     &           IOLDPS, POSELT,
     &           IFLAG, UU, NOFFW, NPVW,
     &           KEEP,KEEP8, MYID_NODES, SEUIL, AVOID_DELAYED,
     &           ETATASS,
     &           DKEEP(1),PIVNUL_LIST(1),LPN_LIST, IWPOS)
         IW( IOLDPS+4+KEEP(IXSZ) ) = STEP(INODE)
         JOBASS = ETATASS  
         IF (JOBASS.EQ.1) THEN
           CALL SMUMPS_252(COMM_LOAD, ASS_IRECV,
     &        N,INODE,IW,LIW,A,LA,
     &        IFLAG,IERROR,ND,
     &        FILS,FRERE,DAD,MAXFRW,root,OPASSW, OPELIW,
     &      PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP,PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS,
     &        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     &    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    JOBASS,ETATASS )
         ENDIF
      ENDIF
      IF (IFLAG.LT.0) GOTO 635
 130  CONTINUE
      TYPE  = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF ( FPERE .NE. 0 ) THEN
        TYPEF = MUMPS_330(PROCNODE_STEPS(STEP(FPERE)),SLAVEF)
      ELSE
        TYPEF = -9999
      END IF
      CALL SMUMPS_254( COMM_LOAD, ASS_IRECV,
     &       N,INODE,TYPE,TYPEF,LA,IW,LIW,A,
     &       IFLAG,IERROR,OPELIW,NELVAW,NMAXNPIV,
     &       PTRIST,PTLUST_S,PTRFAC,
     &       PTRAST, STEP, PIMASTER, PAMASTER,
     &       NE, POSFAC,LRLU,
     &       LRLUS,IPTRLU,ICNTL,KEEP,KEEP8,COMP,IWPOS,IWPOSCB,
     &       PROCNODE_STEPS,SLAVEF,FPERE,COMM_NODES,MYID_NODES,
     &       IPOOL, LPOOL, LEAF,
     &       NSTK_STEPS, NBPROCFILS, BUFR, LBUFR, LBUFR_BYTES, NBFIN,
     &       root, OPASSW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &       INTARR, DBLARR,
     &       ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IF (IFLAG.LT.0) GOTO 640
 200  CONTINUE
      IF ( INODE .eq. KEEP(38) ) THEN
        WRITE(*,*) 'Error .. in SMUMPS_251: ',
     &             ' INODE == KEEP(38)'
        Stop
      END IF
      IF ( FPERE.EQ.0 ) THEN
        NBROOT_TRAITEES = NBROOT_TRAITEES + 1
        IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
           IF (KEEP(201).EQ.1) THEN 
              CALL SMUMPS_681(IERR)
           ELSE IF ( KEEP(201).EQ.2) THEN 
              CALL SMUMPS_580(IERR)
           ENDIF
            NBFIN = NBFIN - NBROOT
            IF ( NBFIN .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in SMUMPS_251: ',
     &                   ' NBFIN=', NBFIN
              CALL MUMPS_ABORT()
            END IF
            IF ( NBROOT .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in SMUMPS_251: ',
     &                   ' NBROOT=', NBROOT
              CALL MUMPS_ABORT()
            END IF
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL SMUMPS_242( DUMMY(1), 1, MPI_INTEGER,
     &          MYID_NODES, COMM_NODES, RACINE, SLAVEF)
            END IF
        ENDIF
        IF (NBFIN.EQ.0)THEN
           GOTO 640
        ENDIF
      ELSEIF ( FPERE.NE.KEEP(38) .AND.
     &         MUMPS_275(PROCNODE_STEPS(STEP(FPERE)),SLAVEF).EQ.
     &         MYID_NODES ) THEN
        NSTK_STEPS(STEP(FPERE)) = NSTK_STEPS(STEP(FPERE))-1
        IF ( NSTK_STEPS( STEP( FPERE )).EQ.0) THEN
          IF (KEEP(234).NE.0 .AND.
     &      MUMPS_167(PROCNODE_STEPS(STEP(INODE)),SLAVEF))
     &      THEN
            STACK_RIGHT_AUTHORIZED = .FALSE.
          ENDIF
          CALL SMUMPS_507(N, IPOOL, LPOOL,
     &         PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     &         KEEP(80), KEEP(47), STEP, FPERE )
          IF (KEEP(47) .GE. 3) THEN
             CALL SMUMPS_500(
     &            IPOOL, LPOOL,
     &            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &            MYID_NODES, STEP, N, ND, FILS )
          ENDIF
          CALL MUMPS_137( FPERE, N, PROCNODE_STEPS,SLAVEF,
     &           ND, FILS, FRERE, STEP, PIMASTER, KEEP(28),
     &           KEEP(50), KEEP(253), FLOP1,
     &           IW, LIW, KEEP(IXSZ) )
          IF (FPERE.NE.KEEP(20))
     &    CALL SMUMPS_190(1,.FALSE.,FLOP1,KEEP,KEEP8)
        ENDIF
      ENDIF
      GO TO 20
 635  CONTINUE
      CALL SMUMPS_44( MYID_NODES, SLAVEF, COMM_NODES )
 640  CONTINUE
        CALL SMUMPS_255( INFO(1),
     &       ASS_IRECV, BUFR, LBUFR,
     &       LBUFR_BYTES,
     &       COMM_NODES,
     &       MYID_NODES, SLAVEF)
       CALL SMUMPS_180( INFO(1),
     &      BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      COMM_NODES, COMM_LOAD, SLAVEF, MP)
      CALL MPI_BARRIER( COMM_NODES, IERR )
       IF ( INFO(1) .GE. 0 ) THEN
          IF( KEEP(38) .NE. 0 .OR. KEEP(20).NE.0) THEN
            MASTER_ROOT = MUMPS_275(
     &                  PROCNODE_STEPS(STEP(max(KEEP(38),KEEP(20)))),
     &                  SLAVEF)
            ROOT_OWNER  = (MASTER_ROOT .EQ. MYID_NODES)
            IF ( KEEP(38) .NE. 0 )THEN
               IF (KEEP(60).EQ.0) THEN
                 IOLDPS  = PTLUST_S(STEP(KEEP(38)))
                 LOCAL_M = IW(IOLDPS+2+KEEP(IXSZ))
                 LOCAL_N = IW(IOLDPS+1+KEEP(IXSZ))
               ELSE
                 IOLDPS  = -999
                 LOCAL_M = root%SCHUR_MLOC
                 LOCAL_N = root%SCHUR_NLOC
               ENDIF
               ITMP8   = int(LOCAL_M,8)*int(LOCAL_N,8)
               LBUFRX = min(int(root%MBLOCK,8)*int(root%NBLOCK,8),
     &            int(root%TOT_ROOT_SIZE,8)*int(root%TOT_ROOT_SIZE,8) )
               IF ( LRLU .GT. LBUFRX ) THEN
                   BUFRX => A(POSFAC:POSFAC+LRLU-1_8)
                   LBUFRX=LRLU
                   IS_BUFRX_ALLOCATED = .FALSE.
               ELSE
                   ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                   IF (IERR.gt.0) THEN
                         INFO(1) = -9
                         CALL MUMPS_731(LBUFRX, INFO(2) )
                         IF (LP > 0 )
     &                   write(LP,*) ' Error allocating, real array ',
     &                   'of size before SMUMPS_146',  LBUFRX
                         CALL MUMPS_ABORT()
                   ENDIF
                   IS_BUFRX_ALLOCATED = .FALSE.
               ENDIF
               CALL SMUMPS_146( MYID_NODES,
     &               root, N, KEEP(38),
     &               COMM_NODES, IW, LIW, IWPOS + 1,
     &               A, LA, PTRAST, PTLUST_S, PTRFAC, STEP,
     &               INFO(1), KEEP(50), KEEP(19),
     &               BUFRX(1), LBUFRX, KEEP,KEEP8, DKEEP )
               IF (IS_BUFRX_ALLOCATED) DEALLOCATE ( BUFRX )
               NULLIFY(BUFRX)
                IF ( MYID_NODES .eq. 
     &               MUMPS_275(PROCNODE_STEPS(STEP(KEEP(38))),
     &                              SLAVEF)
     &             ) THEN
                   IF ( INFO(1) .EQ. -10 .OR. INFO(1) .EQ. -40 ) THEN
                      NPVW = NPVW + INFO(2)
                   ELSE
                      NPVW = NPVW + root%TOT_ROOT_SIZE
                      NMAXNPIV = max(NMAXNPIV,root%TOT_ROOT_SIZE)
                   END IF
                END IF
                IF (root%yes.AND.KEEP(60).EQ.0) THEN
                  IF (KEEP(252).EQ.0) THEN
                  IF (KEEP(201).EQ.1) THEN 
                    CALL MUMPS_729(LAFAC, IW(IOLDPS+XXR))
                    LIWFAC    = IW(IOLDPS+XXI)
                    TYPEFile  = TYPEF_L
                    NextPiv2beWritten = 1 
                    MonBloc%INODE    = KEEP(38)   
                    MonBloc%MASTER   = .TRUE.
                    MonBloc%Typenode = 3
                    MonBloc%NROW     = LOCAL_M
                    MonBloc%NCOL     = LOCAL_N
                    MonBloc%NFS      = MonBloc%NCOL
                    MonBloc%Last     = .TRUE.   
                    MonBloc%LastPiv  =  MonBloc%NCOL
                    NULLIFY(MonBloc%INDICES)
                    STRAT        = STRAT_WRITE_MAX
                    MonBloc%Last = .TRUE.
                    LAST_CALL = .TRUE.
                    CALL SMUMPS_688
     &                                 ( STRAT, TYPEFile,
     &                                  A(PTRFAC(STEP(KEEP(38)))),
     &                                  LAFAC, MonBloc,
     &                                  NextPiv2beWritten, IDUMMY,
     &                                  IW(IOLDPS), LIWFAC,
     &                                  MYID, KEEP8(31), IERR,LAST_CALL)
                  ELSE IF (KEEP(201).EQ.2) THEN
                    KEEP8(31)=KEEP8(31)+ ITMP8
                    CALL SMUMPS_576(KEEP(38),PTRFAC,
     &              KEEP,KEEP8,A,LA, ITMP8, IERR)
                    IF(IERR.LT.0)THEN
                      WRITE(*,*)MYID,
     &                ': Internal error in SMUMPS_576'
                      CALL MUMPS_ABORT()
                    ENDIF
                  ENDIF 
                  ENDIF 
                  IF (KEEP(201).NE.0 .OR. KEEP(252).NE.0) THEN
                     LRLUS = LRLUS + ITMP8 
                     IF (KEEP(252).NE.0) THEN
                       CALL SMUMPS_471(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,0_8,-ITMP8,
     &                 KEEP,KEEP8,LRLU)
                     ELSE         
                       CALL SMUMPS_471(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,ITMP8,    
     &                 0_8,
     &                 KEEP,KEEP8,LRLU)
                     ENDIF
                     IF (PTRFAC(STEP(KEEP(38))).EQ.POSFAC-ITMP8) THEN
                       POSFAC = POSFAC  - ITMP8
                       LRLU   = LRLU    + ITMP8
                     ENDIF
                  ELSE
                       CALL SMUMPS_471(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,ITMP8,    
     &                 0_8,
     &                 KEEP,KEEP8,LRLU)
                  ENDIF
                ENDIF  
                IF (root%yes. AND. KEEP(252) .NE. 0 .AND.
     &              (KEEP(60).EQ.0 .OR. KEEP(221).EQ.1)) THEN
                  IF (MYID_NODES .EQ. MASTER_ROOT) THEN
                    LRHS_CNTR_MASTER_ROOT = root%TOT_ROOT_SIZE*KEEP(253)
                  ELSE
                    LRHS_CNTR_MASTER_ROOT = 1
                  ENDIF
                  ALLOCATE(root%RHS_CNTR_MASTER_ROOT( 
     &                     LRHS_CNTR_MASTER_ROOT), stat=IERR )
                  IF (IERR.gt.0) THEN
                    INFO(1) = -13
                    CALL MUMPS_731(LRHS_CNTR_MASTER_ROOT,INFO(2))
                    IF (LP > 0 )
     &              write(LP,*) ' Error allocating, real array ',
     &              'of size before SMUMPS_146',
     &              LRHS_CNTR_MASTER_ROOT
                    CALL MUMPS_ABORT()
                  ENDIF
                  FWD_LOCAL_N_RHS = numroc(KEEP(253), root%NBLOCK,
     &            root%MYCOL, 0, root%NPCOL)
                  FWD_LOCAL_N_RHS = max(1,FWD_LOCAL_N_RHS)
                  CALL SMUMPS_156( MYID_NODES,
     &            root%TOT_ROOT_SIZE, KEEP(253),
     &            root%RHS_CNTR_MASTER_ROOT(1), LOCAL_M,
     &            FWD_LOCAL_N_RHS, root%MBLOCK, root%NBLOCK,
     &            root%RHS_ROOT(1,1), MASTER_ROOT,
     &            root%NPROW, root%NPCOL, COMM_NODES )
     &             
                ENDIF
            ELSE
                IF (KEEP(19).NE.0) THEN
                  CALL MPI_REDUCE(KEEP(109), GLOBK109, 1,
     &                 MPI_INTEGER, MPI_SUM,
     &                 MASTER_ROOT,
     &                 COMM_NODES, IERR)
                ENDIF
                IF (ROOT_OWNER) THEN
                   IPOSROOT = PTLUST_S(STEP(KEEP(20)))
                   NFRONT   = IW(IPOSROOT+KEEP(IXSZ)+3)   
                   NFRONT8  = int(NFRONT,8)
                   IPOSROOTROWINDICES=IPOSROOT+6+KEEP(IXSZ)+ 
     &                             IW(IPOSROOT+5+KEEP(IXSZ)) 
                   NPVW = NPVW + NFRONT 
                   NMAXNPIV = max(NMAXNPIV,NFRONT)
                END IF
               IF (ROOT_OWNER.AND.KEEP(60).NE.0) THEN 
                  IF ( PTRFAC(STEP(KEEP(20))) .EQ. POSFAC -
     &                 NFRONT8*NFRONT8 ) THEN
                    POSFAC = POSFAC - NFRONT8*NFRONT8
                    LRLUS  = LRLUS  + NFRONT8*NFRONT8
                    LRLU   = LRLUS  + NFRONT8*NFRONT8
                    CALL SMUMPS_471(.FALSE.,.FALSE.,
     &              LA-LRLUS,0_8,-NFRONT8*NFRONT8,KEEP,KEEP8,LRLU)
                  ENDIF
               ENDIF
            END IF
          END IF
       END IF
       IF ( KEEP(38) .NE. 0 ) THEN
         IF (MYID_NODES.EQ. 
     &        MUMPS_275(PROCNODE_STEPS(STEP(KEEP(38))),SLAVEF)
     &      ) THEN
           MAXFRW = max ( MAXFRW, root%TOT_ROOT_SIZE)
         END IF
       END IF
       MAXFRT       = MAXFRW
       NTOTPV       = NPVW
       INFO(12)     = NOFFW
       RINFO(2)     = real(OPASSW)
       RINFO(3)     = real(OPELIW)
       INFO(13)     = NELVAW
       INFO(14)     = COMP
      RETURN
      END SUBROUTINE SMUMPS_251
      SUBROUTINE SMUMPS_87( HEADER, KEEP253 )
        INTEGER HEADER( 6 ), KEEP253
        INTEGER NFRONT, NASS
        NFRONT = HEADER(1)
        IF ( HEADER(2) .ne. 0 ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 1 :',HEADER(2)
          CALL MUMPS_ABORT()
        END IF
        NASS   = abs( HEADER( 3 ) )
        IF ( NASS .NE. abs( HEADER( 4 ) ) ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 2 :',HEADER(3:4)
          CALL MUMPS_ABORT()
        END IF
        IF ( NASS+KEEP253 .NE. NFRONT ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 3 : not root'
          CALL MUMPS_ABORT()
        END IF
        HEADER( 1 ) = KEEP253 
        HEADER( 2 ) = 0
        HEADER( 3 ) = NFRONT 
        HEADER( 4 ) = NFRONT-KEEP253    
        RETURN
      END SUBROUTINE SMUMPS_87
      SUBROUTINE SMUMPS_136( id )
      USE SMUMPS_OOC
      USE SMUMPS_STRUC_DEF
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      include 'mpif.h'
      TYPE( SMUMPS_STRUC ) :: id
      LOGICAL I_AM_SLAVE
      INTEGER IERR, MASTER
      PARAMETER ( MASTER = 0 )
      I_AM_SLAVE = ( id%MYID .ne. MASTER .OR. id%KEEP(46) .NE. 0 )
      IF (id%KEEP(201).GT.0 .AND. I_AM_SLAVE) THEN
        CALL SMUMPS_587(id,IERR)
        IF (IERR < 0) THEN
          id%INFO(1) = -90
          id%INFO(2) = 0
        ENDIF
      END IF
      CALL MUMPS_276(id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID)
      IF (id%root%gridinit_done) THEN
        IF ( id%KEEP(38).NE.0 .and. id%root%yes ) THEN
          CALL blacs_gridexit( id%root%CNTXT_BLACS )
          id%root%gridinit_done = .FALSE.
        END IF
      END IF
      IF ( id%MYID .NE. MASTER .OR. id%KEEP(46) .ne. 0 ) THEN
        CALL MPI_COMM_FREE( id%COMM_NODES, IERR )
        CALL MPI_COMM_FREE( id%COMM_LOAD, IERR )
      END IF
      IF (associated(id%MEM_DIST))  THEN
         DEALLOCATE(id%MEM_DIST)
         NULLIFY(id%MEM_DIST)
      ENDIF
      IF (associated(id%MAPPING)) THEN
        DEALLOCATE(id%MAPPING)
        NULLIFY(id%MAPPING)
      END IF
       NULLIFY(id%SCHUR_CINTERFACE)
      IF ( id%KEEP(52) .NE. -1 .or. id%MYID .ne. MASTER ) THEN
        IF (associated(id%COLSCA)) THEN
          DEALLOCATE(id%COLSCA)
          NULLIFY(id%COLSCA)
        ENDIF
        IF (associated(id%ROWSCA)) THEN
          DEALLOCATE(id%ROWSCA)
          NULLIFY(id%ROWSCA)
        ENDIF
      END IF
      IF (associated(id%PTLUST_S)) THEN
        DEALLOCATE(id%PTLUST_S)
        NULLIFY(id%PTLUST_S)
      END IF
      IF (associated(id%PTRFAC)) THEN
        DEALLOCATE(id%PTRFAC)
        NULLIFY(id%PTRFAC)
      END IF
      IF (associated(id%POIDS)) THEN
        DEALLOCATE(id%POIDS)
        NULLIFY(id%POIDS)
      ENDIF
      IF (associated(id%IS)) THEN
        DEALLOCATE(id%IS)
        NULLIFY(id%IS)
      ENDIF
      IF (associated(id%IS1)) THEN
        DEALLOCATE(id%IS1)
        NULLIFY(id%IS1)
      ENDIF
      IF (associated(id%STEP))      THEN
        DEALLOCATE(id%STEP)
        NULLIFY(id%STEP)
      ENDIF
      IF (associated(id%Step2node))      THEN
        DEALLOCATE(id%Step2node)
        NULLIFY(id%Step2node)
      ENDIF
      IF (associated(id%NE_STEPS))  THEN
        DEALLOCATE(id%NE_STEPS)
        NULLIFY(id%NE_STEPS)
      ENDIF
      IF (associated(id%ND_STEPS))  THEN
        DEALLOCATE(id%ND_STEPS)
        NULLIFY(id%ND_STEPS)
      ENDIF
      IF (associated(id%FRERE_STEPS))  THEN
        DEALLOCATE(id%FRERE_STEPS)
        NULLIFY(id%FRERE_STEPS)
      ENDIF
      IF (associated(id%DAD_STEPS))  THEN
        DEALLOCATE(id%DAD_STEPS)
        NULLIFY(id%DAD_STEPS)
      ENDIF
      IF (associated(id%SYM_PERM))  THEN
        DEALLOCATE(id%SYM_PERM)
        NULLIFY(id%SYM_PERM)
      ENDIF
      IF (associated(id%UNS_PERM))  THEN
        DEALLOCATE(id%UNS_PERM)
        NULLIFY(id%UNS_PERM)
      ENDIF
      IF (associated(id%PIVNUL_LIST))  THEN
        DEALLOCATE(id%PIVNUL_LIST)
        NULLIFY(id%PIVNUL_LIST)
      ENDIF
      IF (associated(id%FILS))      THEN
        DEALLOCATE(id%FILS)
        NULLIFY(id%FILS)
      ENDIF
      IF (associated(id%PTRAR))     THEN
        DEALLOCATE(id%PTRAR)
        NULLIFY(id%PTRAR)
      ENDIF
      IF (associated(id%FRTPTR))    THEN
        DEALLOCATE(id%FRTPTR)
        NULLIFY(id%FRTPTR)
      ENDIF
      IF (associated(id%FRTELT))    THEN
        DEALLOCATE(id%FRTELT)
        NULLIFY(id%FRTELT)
      ENDIF
      IF (associated(id%NA))        THEN
        DEALLOCATE(id%NA)
        NULLIFY(id%NA)
      ENDIF
      IF (associated(id%PROCNODE_STEPS)) THEN
        DEALLOCATE(id%PROCNODE_STEPS)
        NULLIFY(id%PROCNODE_STEPS)
      ENDIF
      IF (associated(id%PROCNODE)) THEN
        DEALLOCATE(id%PROCNODE)
        NULLIFY(id%PROCNODE)
      ENDIF
      IF (associated(id%RHSCOMP)) THEN
        DEALLOCATE(id%RHSCOMP)
        NULLIFY(id%RHSCOMP)
      ENDIF
      IF (associated(id%POSINRHSCOMP)) THEN
        DEALLOCATE(id%POSINRHSCOMP)
        NULLIFY(id%POSINRHSCOMP)
      ENDIF
      IF (id%KEEP(46).eq.1 .and.
     &    id%KEEP(55).ne.0 .and.
     &    id%MYID .eq. MASTER .and.
     &    id%KEEP(52) .eq. 0 ) THEN
        NULLIFY(id%DBLARR)
      ELSE
        IF (associated(id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
      END IF
      IF (associated(id%INTARR))       THEN
        DEALLOCATE(id%INTARR)
        NULLIFY(id%INTARR)
      ENDIF
      IF (associated(id%root%RG2L_ROW))THEN
        DEALLOCATE(id%root%RG2L_ROW)
        NULLIFY(id%root%RG2L_ROW)
      ENDIF
      IF (associated(id%root%RG2L_COL))THEN
        DEALLOCATE(id%root%RG2L_COL)
        NULLIFY(id%root%RG2L_COL)
      ENDIF
      IF (associated(id%root%IPIV))    THEN
        DEALLOCATE(id%root%IPIV)
        NULLIFY(id%root%IPIV)
      ENDIF
      IF (associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE(id%root%RHS_CNTR_MASTER_ROOT)
        NULLIFY(id%root%RHS_CNTR_MASTER_ROOT)
      ENDIF
      IF (associated(id%root%RHS_ROOT))THEN
        DEALLOCATE(id%root%RHS_ROOT)
        NULLIFY(id%root%RHS_ROOT)
      ENDIF
      CALL SMUMPS_636(id)
      IF (associated(id%ELTPROC))     THEN
        DEALLOCATE(id%ELTPROC)
        NULLIFY(id%ELTPROC)
      ENDIF
      IF (associated(id%CANDIDATES)) THEN
        DEALLOCATE(id%CANDIDATES)
        NULLIFY(id%CANDIDATES)
      ENDIF
      IF (associated(id%I_AM_CAND)) THEN
        DEALLOCATE(id%I_AM_CAND)
        NULLIFY(id%I_AM_CAND)
      ENDIF
      IF (associated(id%ISTEP_TO_INIV2)) THEN
        DEALLOCATE(id%ISTEP_TO_INIV2)
        NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
      IF (I_AM_SLAVE) THEN
       IF (associated(id%TAB_POS_IN_PERE)) THEN
        DEALLOCATE(id%TAB_POS_IN_PERE)
        NULLIFY(id%TAB_POS_IN_PERE)
       ENDIF
       IF (associated(id%FUTURE_NIV2)) THEN
        DEALLOCATE(id%FUTURE_NIV2)
        NULLIFY(id%FUTURE_NIV2)
       ENDIF
      ENDIF
      IF(associated(id%DEPTH_FIRST))THEN
        DEALLOCATE(id%DEPTH_FIRST)
        NULLIFY(id%DEPTH_FIRST)
      ENDIF
      IF(associated(id%DEPTH_FIRST_SEQ))THEN
        DEALLOCATE(id%DEPTH_FIRST_SEQ)
        NULLIFY(id%DEPTH_FIRST_SEQ)
      ENDIF
      IF(associated(id%SBTR_ID))THEN
        DEALLOCATE(id%SBTR_ID)
        NULLIFY(id%SBTR_ID)
      ENDIF
      IF (associated(id%MEM_SUBTREE)) THEN
        DEALLOCATE(id%MEM_SUBTREE)
        NULLIFY(id%MEM_SUBTREE)
      ENDIF
      IF (associated(id%MY_ROOT_SBTR)) THEN
        DEALLOCATE(id%MY_ROOT_SBTR)
        NULLIFY(id%MY_ROOT_SBTR)
      ENDIF
      IF (associated(id%MY_FIRST_LEAF)) THEN
        DEALLOCATE(id%MY_FIRST_LEAF)
        NULLIFY(id%MY_FIRST_LEAF)
      ENDIF
      IF (associated(id%MY_NB_LEAF)) THEN
        DEALLOCATE(id%MY_NB_LEAF)
        NULLIFY(id%MY_NB_LEAF)
      ENDIF
      IF (associated(id%COST_TRAV)) THEN
        DEALLOCATE(id%COST_TRAV)
        NULLIFY(id%COST_TRAV)
      ENDIF
        IF(associated (id%OOC_INODE_SEQUENCE))THEN
           DEALLOCATE(id%OOC_INODE_SEQUENCE)
           NULLIFY(id%OOC_INODE_SEQUENCE)
        ENDIF
        IF(associated (id%OOC_TOTAL_NB_NODES))THEN
           DEALLOCATE(id%OOC_TOTAL_NB_NODES)
           NULLIFY(id%OOC_TOTAL_NB_NODES)
        ENDIF
        IF(associated (id%OOC_SIZE_OF_BLOCK))THEN
           DEALLOCATE(id%OOC_SIZE_OF_BLOCK)
           NULLIFY(id%OOC_SIZE_OF_BLOCK)
        ENDIF
        IF(associated (id%OOC_VADDR))THEN
           DEALLOCATE(id%OOC_VADDR)
           NULLIFY(id%OOC_VADDR)
        ENDIF
        IF(associated (id%OOC_NB_FILES))THEN
           DEALLOCATE(id%OOC_NB_FILES)
           NULLIFY(id%OOC_NB_FILES)
        ENDIF
      IF (id%KEEP8(24).EQ.0_8) THEN
        IF (associated(id%S))        DEALLOCATE(id%S)
      ELSE
      ENDIF
      NULLIFY(id%S)
      IF (I_AM_SLAVE) THEN
        CALL SMUMPS_57( IERR )
        CALL SMUMPS_59( IERR )
      END IF
      IF ( associated( id%BUFR ) ) DEALLOCATE( id%BUFR )
      NULLIFY( id%BUFR )
      RETURN
      END SUBROUTINE SMUMPS_136
      SUBROUTINE SMUMPS_150(MYID,COMM,S,MAXS,MAXS_BYTES)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      INTEGER COMM, MYID, MAXS, MAXS_BYTES
      INTEGER S( MAXS )
      INTEGER MSGTAG, MSGSOU, MSGLEN
      LOGICAL FLAG
      FLAG = .TRUE.
      DO WHILE ( FLAG )
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &                   FLAG, STATUS, IERR )
        IF (FLAG) THEN
          MSGTAG=STATUS(MPI_TAG)
          MSGSOU=STATUS(MPI_SOURCE)
          CALL MPI_GET_COUNT(STATUS,MPI_PACKED,MSGLEN,IERR)
          IF (MSGLEN <= MAXS_BYTES) THEN
            CALL MPI_RECV(S(1),MAXS_BYTES,MPI_PACKED,
     &      MSGSOU, MSGTAG, COMM, STATUS, IERR)
          ELSE
            EXIT
          ENDIF
        END IF
      END DO
      CALL MPI_BARRIER( COMM, IERR )
      RETURN
      END SUBROUTINE SMUMPS_150
      SUBROUTINE SMUMPS_254(COMM_LOAD, ASS_IRECV, 
     &    N, INODE, TYPE, TYPEF, 
     &    LA, IW, LIW, A,
     &    IFLAG, IERROR, OPELIW, NELVAW, NMAXNPIV,
     &    PTRIST, PTLUST_S,
     &    PTRFAC, PTRAST,
     &    STEP, PIMASTER, PAMASTER, NE,
     &    POSFAC, LRLU, LRLUS, IPTRLU, ICNTL, KEEP,KEEP8,
     &    COMP, IWPOS, IWPOSCB, PROCNODE_STEPS, SLAVEF,
     &    FPERE, COMM, MYID,
     &    IPOOL, LPOOL, LEAF, NSTK_S,
     &    NBPROCFILS,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, root,
     &    OPASSW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE SMUMPS_COMM_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER COMM, MYID, TYPE, TYPEF
      INTEGER N, LIW, INODE,IFLAG,IERROR
      INTEGER ICNTL(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA, POSFAC, LRLU, LRLUS, IPTRLU
      INTEGER IWPOSCB, IWPOS,
     &        FPERE, SLAVEF, NELVAW, NMAXNPIV
      INTEGER IW(LIW),PROCNODE_STEPS(KEEP(28))
      INTEGER(8) :: PTRAST  (KEEP(28))
      INTEGER(8) :: PTRFAC  (KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28)), NE(KEEP(28))
      REAL    A(LA)
      DOUBLE PRECISION OPASSW, OPELIW
      REAL DBLARR(max(1,KEEP(13)))
      INTEGER INTARR(max(1,KEEP(14))) 
      INTEGER ITLOC( N + KEEP(253) ), FILS( N ),
     &        ND( KEEP(28) ), FRERE( KEEP(28) )
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER LPOOL, LEAF, COMP
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK_S( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NBFIN
      INTEGER NFRONT_ESTIM,NELIM_ESTIM
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LP
      INTEGER NBROWS_ALREADY_SENT 
      INTEGER(8) :: POSELT, OPSFAC
      INTEGER(8) :: IOLD, INEW, FACTOR_POS
      INTEGER NSLAVES, NCB,
     &        H_INODE, IERR, NBCOL, NBROW, NBROW_SEND,
     &        NBROW_STACK, NBCOL_STACK, NELIM
      INTEGER NCBROW_ALREADY_MOVED, NCBROW_PREVIOUSLY_MOVED,
     &NCBROW_NEWLY_MOVED
      INTEGER(8) :: LAST_ALLOWED_POS
      INTEGER(8) :: LREQCB, MIN_SPACE_IN_PLACE
      INTEGER(8) :: SHIFT_VAL_SON
      INTEGER SHIFT_LIST_ROW_SON,
     &        SHIFT_LIST_COL_SON,
     &        LIST_ROW_SON, LIST_COL_SON, LIST_SLAVES
      INTEGER IOLDPS,NFRONT,NPIV,NASS,IOLDP1,PTROWEND,
     &        LREQI, LCONT
      INTEGER I,LDA, INIV2
      INTEGER MSGDEST, MSGTAG, CHK_LOAD
      INCLUDE 'mumps_headers.h'
      LOGICAL  COMPRESSCB, MUST_COMPACT_FACTORS
      LOGICAL  INPLACE
      INTEGER(8) :: SIZE_INPLACE
      INTEGER INTSIZ
      DOUBLE PRECISION FLOP1, FLOP1_EFFECTIVE
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      LOGICAL SSARBR, SSARBR_ROOT, MUMPS_167,
     &MUMPS_170
      EXTERNAL MUMPS_167, MUMPS_170
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      INPLACE = .FALSE.
      MIN_SPACE_IN_PLACE = 0_8
      IOLDPS = PTLUST_S(STEP(INODE))
      INTSIZ = IW(IOLDPS+XXI)
      NFRONT = IW(IOLDPS+KEEP(IXSZ))
      NPIV   = IW(IOLDPS + 1+KEEP(IXSZ))         
      NMAXNPIV = max(NPIV, NMAXNPIV)
      NASS   = iabs(IW(IOLDPS + 2+KEEP(IXSZ)))   
      NSLAVES=  IW(IOLDPS+5+KEEP(IXSZ))
      H_INODE= 6 + NSLAVES + KEEP(IXSZ)
      LCONT = NFRONT - NPIV
      NBCOL = LCONT
      SSARBR = MUMPS_167(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      SSARBR_ROOT = MUMPS_170
     &              (PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      LREQCB = 0_8
      INPLACE = .FALSE.
      COMPRESSCB= ((KEEP(215).EQ.0)
     &             .AND.(KEEP(50).NE.0)
     &             .AND.(TYPEF.EQ.1
     &             .OR.TYPEF.EQ.2
     &              )
     &             .AND.(TYPE.EQ.1))
      MUST_COMPACT_FACTORS = .TRUE.
      IF (KEEP(201).EQ.1 .OR. KEEP(201).EQ.-1) THEN
            MUST_COMPACT_FACTORS = .FALSE.
      ENDIF
      IF ((FPERE.EQ.0).AND.(NASS.NE.NPIV)) THEN 
        IFLAG = -10
        GOTO 600
      ENDIF
      NBROW      = LCONT
      IF (TYPE.EQ.2) NBROW = NASS - NPIV
      IF ((KEEP(50).NE.0).AND.(TYPE.EQ.2)) THEN
        LDA = NASS
      ELSE
        LDA = NFRONT
      ENDIF
      NBROW_SEND = NBROW
      NELIM = NASS-NPIV
      IF (TYPEF.EQ.2) NBROW_SEND = NELIM 
      POSELT = PTRAST(STEP(INODE))
      IF (POSELT .ne. PTRFAC(STEP(INODE))) THEN
        WRITE(*,*) "Error 1 in G"
        CALL MUMPS_ABORT()
      END IF
      NELVAW = NELVAW + NASS - NPIV
      IF (KEEP(50) .eq. 0) THEN
        KEEP8(10) = KEEP8(10) + int(NPIV,8) * int(NFRONT,8)
      ELSE
        KEEP8(10) = KEEP8(10) + ( int(NPIV,8)*int(NPIV+1,8) )/ 2_8
      ENDIF
      KEEP8(10) = KEEP8(10) + int(NBROW,8) * int(NPIV,8)
      CALL MUMPS_511( NFRONT, NPIV, NASS,
     &     KEEP(50), TYPE,FLOP1 )
      IF ( (.NOT. SSARBR_ROOT) .and. TYPE == 1) THEN
        IF (NE(STEP(INODE))==0) THEN
          CHK_LOAD=0
        ELSE
          CHK_LOAD=1
        ENDIF
        CALL SMUMPS_190(CHK_LOAD, .FALSE., -FLOP1,
     &                      KEEP,KEEP8)
      ENDIF
      FLOP1_EFFECTIVE = FLOP1
      OPELIW = OPELIW + FLOP1
      IF ( NPIV .NE. NASS ) THEN
        CALL MUMPS_511( NFRONT, NASS, NASS,
     &       KEEP(50), TYPE,FLOP1 )
         IF (.NOT. SSARBR_ROOT ) THEN
            IF (NE(STEP(INODE))==0) THEN
              CHK_LOAD=0
            ELSE
              CHK_LOAD=1
            ENDIF
            CALL SMUMPS_190(CHK_LOAD, .FALSE.,
     &                         FLOP1_EFFECTIVE-FLOP1,
     &                         KEEP,KEEP8)
         ENDIF
      END IF
      IF ( SSARBR_ROOT ) THEN
        NFRONT_ESTIM=ND(STEP(INODE)) + KEEP(253)
        NELIM_ESTIM=NASS-(NFRONT-NFRONT_ESTIM)
        CALL MUMPS_511(NFRONT_ESTIM,NELIM_ESTIM,NELIM_ESTIM,
     &       KEEP(50),1,FLOP1)
      END IF
      FLOP1=-FLOP1
      IF (SSARBR_ROOT) THEN
        CALL SMUMPS_190(0,.FALSE.,FLOP1,KEEP,KEEP8)
      ELSE
        CALL SMUMPS_190(2,.FALSE.,FLOP1,KEEP,KEEP8)
      ENDIF
      IF ( FPERE .EQ. 0 ) THEN
        IF ( KEEP(253) .NE. 0 .AND. KEEP(201).NE.-1
     &       .AND. KEEP(201).NE.1 ) THEN
          MUST_COMPACT_FACTORS = .TRUE.
          GOTO 190 
        ELSE
          MUST_COMPACT_FACTORS = .FALSE.
          GOTO 190 
        ENDIF
      ENDIF
      IF ( FPERE.EQ.KEEP(38) ) THEN
       NCB   = NFRONT - NASS
       SHIFT_LIST_ROW_SON = H_INODE + NASS
       SHIFT_LIST_COL_SON = H_INODE + NFRONT + NASS
       SHIFT_VAL_SON      = int(NASS,8)*int(NFRONT+1,8)
       IF (TYPE.EQ.1) THEN
         CALL SMUMPS_80(
     &    COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE,
     &    PTLUST_S, PTRAST,
     &    root, NCB, NCB, SHIFT_LIST_ROW_SON,
     &    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, NFRONT,
     &    ROOT_CONT_STATIC, MYID, COMM,
     &
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, .FALSE., ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF (IFLAG < 0 ) GOTO 500
       ENDIF
       MSGDEST=  MUMPS_275(PROCNODE_STEPS(STEP(FPERE)),SLAVEF)
       IOLDPS = PTLUST_S(STEP(INODE))
       LIST_ROW_SON = IOLDPS + H_INODE + NPIV
       LIST_COL_SON = IOLDPS + H_INODE + NFRONT + NPIV
       LIST_SLAVES  = IOLDPS + 6 + KEEP(IXSZ)
       IF (MSGDEST.EQ.MYID) THEN
         CALL SMUMPS_273( root, 
     &      INODE, NELIM, NSLAVES, IW(LIST_ROW_SON), 
     &      IW(LIST_COL_SON), IW(LIST_SLAVES),
     &
     &      PROCNODE_STEPS, IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S,
     &      ITLOC, RHS_MUMPS, COMP,
     &      IFLAG, IERROR, 
     &      IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,KEEP8,
     &      COMM, COMM_LOAD, FILS, ND)
         IF (IFLAG.LT.0) GOTO 600
       ELSE
        IERR = -1
        DO WHILE (IERR.EQ.-1)
         CALL SMUMPS_76( INODE, NELIM, 
     &     IW(LIST_ROW_SON), IW(LIST_COL_SON), NSLAVES, 
     &     IW(LIST_SLAVES), MSGDEST, COMM, IERR)
         IF ( IERR .EQ. -1 ) THEN
            BLOCKING   =.FALSE.
            SET_IRECV  =.TRUE.
            MESSAGE_RECEIVED = .FALSE.
            CALL SMUMPS_329( COMM_LOAD, ASS_IRECV, 
     &       BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &       MPI_ANY_SOURCE, MPI_ANY_TAG,
     &       STATUS,
     &       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &       IWPOS, IWPOSCB, IPTRLU,
     &       LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &       PTLUST_S, PTRFAC,
     &       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &       IFLAG, IERROR, COMM,
     &       NBPROCFILS,
     &       IPOOL, LPOOL, LEAF,
     &       NBFIN, MYID, SLAVEF,
     &
     &       root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &       FILS, PTRARW, PTRAIW,
     &       INTARR, DBLARR, ICNTL, KEEP,KEEP8,
     &       ND, FRERE, LPTRAR, NELT,
     &       FRTPTR, FRTELT, 
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &       .TRUE.)
            IF ( IFLAG .LT. 0 ) GOTO 500
            IOLDPS = PTLUST_S(STEP(INODE))
            LIST_ROW_SON = IOLDPS + H_INODE + NPIV
            LIST_COL_SON = IOLDPS + H_INODE + NFRONT + NPIV
            LIST_SLAVES  = IOLDPS + 6 + KEEP(IXSZ)
         ENDIF
        ENDDO
        IF ( IERR .EQ. -2 ) THEN 
            IERROR = ( 3 + NSLAVES + 2 * NELIM ) * KEEP( 34 ) 
            IFLAG  = - 17
            GOTO 600
        ELSE IF ( IERR .EQ. -3 ) THEN
            IERROR = ( 3 + NSLAVES + 2 * NELIM ) * KEEP( 34 ) 
            IFLAG  = -20
            GOTO 600
        ENDIF
       ENDIF
       IF (NELIM.EQ.0) THEN 
          POSELT = PTRAST(STEP(INODE))
          OPSFAC = POSELT + int(NPIV,8) * int(NFRONT,8) + int(NPIV,8)
          GOTO 190
       ELSE
          GOTO 500
       ENDIF
      ENDIF
      OPSFAC = POSELT + int(NPIV,8) * int(LDA,8) + int(NPIV,8)
      IF ( MUMPS_275(PROCNODE_STEPS(STEP(FPERE)),
     &     SLAVEF) .NE. MYID ) THEN
        MSGTAG =NOEUD
        MSGDEST=MUMPS_275( PROCNODE_STEPS(STEP(FPERE)), SLAVEF )
        IERR = -1 
        NBROWS_ALREADY_SENT = 0
        DO WHILE (IERR.EQ.-1)
          IF ( (TYPE.EQ.1) .AND. (TYPEF.EQ.1) ) THEN
             CALL SMUMPS_66( NBROWS_ALREADY_SENT,
     &         INODE, FPERE, NFRONT, 
     &         LCONT, NASS, NPIV, IW( IOLDPS +  H_INODE + NPIV ),
     &         IW( IOLDPS +  H_INODE + NPIV + NFRONT ),
     &         A( OPSFAC ), COMPRESSCB,
     &         MSGDEST, MSGTAG, COMM, IERR )
          ELSE
             IF ( (TYPE.EQ.2).AND.(KEEP(48).NE.0)) THEN
              INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
             ELSE
              INIV2 = -9999
             ENDIF
             CALL SMUMPS_70( NBROWS_ALREADY_SENT,
     &             FPERE, INODE, 
     &             NBROW_SEND, IW(IOLDPS +  H_INODE + NPIV ), 
     &             NBCOL, IW(IOLDPS +  H_INODE + NPIV + NFRONT ),
     &             A(OPSFAC), LDA, NELIM, TYPE,
     &             NSLAVES, IW(IOLDPS+6+KEEP(IXSZ)), MSGDEST, 
     &             COMM, IERR, 
     & 
     &             SLAVEF, KEEP,KEEP8, INIV2, TAB_POS_IN_PERE )
          END IF
          IF ( IERR .EQ. -1 ) THEN
            BLOCKING  = .FALSE.
            SET_IRECV = .TRUE.
            MESSAGE_RECEIVED = .FALSE.
            CALL SMUMPS_329( COMM_LOAD, ASS_IRECV, 
     &       BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &       MPI_ANY_SOURCE, MPI_ANY_TAG,
     &       STATUS,
     &       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &       IWPOS, IWPOSCB, IPTRLU,
     &       LRLU, LRLUS, N, IW, LIW, A, LA,
     &       PTRIST, PTLUST_S, PTRFAC,
     &       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &       IFLAG, IERROR, COMM,
     &       NBPROCFILS,
     &       IPOOL, LPOOL, LEAF,
     &       NBFIN, MYID, SLAVEF,
     &
     &       root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &       FILS, PTRARW, PTRAIW,
     &       INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &       LPTRAR, NELT, FRTPTR, FRTELT,
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.  )
             IF ( IFLAG .LT. 0 ) GOTO 500
          ENDIF
          IOLDPS = PTLUST_S(STEP( INODE ))
          OPSFAC = POSELT + int(NPIV,8) * int(LDA,8) + int(NPIV,8)
        END DO
        IF ( IERR .EQ. -2 .OR. IERR .EQ. -3 ) THEN 
          IF ( (TYPE.EQ.1) .AND. (TYPEF.EQ.1) ) THEN
            IERROR = ( 2*LCONT + 9 ) * KEEP( 34 ) + 
     &          LCONT*LCONT * KEEP( 35 )
          ELSE IF (KEEP(50).ne.0 .AND. TYPE .eq. 2 ) THEN
            IERROR = ( NBROW_SEND + NBCOL+ 5 + NSLAVES) 
     &                 * KEEP( 34 ) + 
     &          NBROW_SEND*NBROW_SEND*KEEP( 35 )
          ELSE
            IERROR = ( NBROW_SEND + NBCOL+ 5 + NSLAVES) * KEEP( 34 ) + 
     &          NBROW_SEND*NBCOL*KEEP( 35 )
          ENDIF
          IF (IERR .EQ. -2) THEN
            IFLAG = -17
            IF ( LP  >  0 ) THEN
              WRITE(LP, *) MYID,
     & ": FAILURE, SEND BUFFER TOO SMALL DURING
     & SMUMPS_254", TYPE, TYPEF
            ENDIF
          ENDIF
          IF (IERR .EQ. -3) THEN
            IFLAG = -20
            IF ( LP  >  0 ) THEN
              WRITE(LP, *) MYID,
     & ": FAILURE, RECV BUFFER TOO SMALL DURING
     & SMUMPS_254", TYPE, TYPEF
            ENDIF
          ENDIF
          GOTO 600
        ENDIF
      ENDIF
      IF ( MUMPS_275(PROCNODE_STEPS(STEP(FPERE)),
     &     SLAVEF) .EQ. MYID ) THEN
        LREQI = 2 + KEEP(IXSZ)
        NBROW_STACK = NBROW
        NBROW_SEND = 0
        IF ((KEEP(50).NE.0).AND.(TYPE.EQ.2)) THEN
          NBCOL_STACK = NBROW
        ELSE
          NBCOL_STACK = NBCOL
        ENDIF
      ELSE
        NBROW_STACK = NBROW-NBROW_SEND
        NBCOL_STACK = NBCOL
        LREQI       = 6 + NBROW_STACK + NBCOL + KEEP(IXSZ)
        IF (.NOT. (TYPE.EQ.1 .AND. TYPEF.EQ.2 ) ) GOTO 190
        IF (FPERE.EQ.0) GOTO 190
      ENDIF
      IF (COMPRESSCB) THEN
        LREQCB = ( int(NBCOL_STACK,8) * int( NBCOL_STACK + 1, 8) ) / 2_8
     &         - ( int(NBROW_SEND ,8) * int( NBROW_SEND  + 1, 8) ) / 2_8
      ELSE
        LREQCB = int(NBROW_STACK,8) * int(NBCOL_STACK,8)
      ENDIF
      INPLACE = ( KEEP(234).NE.0 )
      IF (KEEP(50).NE.0 .AND. TYPE .EQ. 2) INPLACE = .FALSE.
      INPLACE = INPLACE .OR. .NOT. MUST_COMPACT_FACTORS
      INPLACE = INPLACE .AND.
     &            ( PTLUST_S(STEP(INODE)) + INTSIZ .EQ. IWPOS )
      MIN_SPACE_IN_PLACE = 0_8
      IF ( INPLACE .AND. KEEP(50).eq. 0 .AND.
     &     MUST_COMPACT_FACTORS) THEN
        MIN_SPACE_IN_PLACE = int(NBCOL_STACK,8)
      ENDIF
      IF ( MIN_SPACE_IN_PLACE .GT. LREQCB ) THEN
        INPLACE = .FALSE.
      ENDIF
      CALL SMUMPS_22( INPLACE, MIN_SPACE_IN_PLACE,
     &   SSARBR, .FALSE.,
     &   MYID,N,KEEP,KEEP8,IW, LIW, A, LA,
     &   LRLU, IPTRLU,IWPOS,IWPOSCB,
     &   PTRIST,PTRAST,STEP, PIMASTER,PAMASTER,
     &   LREQI, LREQCB, INODE, S_NOTFREE, .TRUE.,
     &   COMP, LRLUS, IFLAG, IERROR )
      IF (IFLAG.LT.0) GOTO 600
      PTRIST(STEP(INODE)) = IWPOSCB+1
      IF ( MUMPS_275(PROCNODE_STEPS(STEP(FPERE)),
     &     SLAVEF) .EQ. MYID ) THEN
        PIMASTER (STEP(INODE)) = PTLUST_S(STEP(INODE))
        PAMASTER(STEP(INODE)) = IPTRLU + 1_8
        PTRAST(STEP(INODE)) = -99999999_8
          IW(IWPOSCB+1+KEEP(IXSZ)) = min(-NBCOL_STACK,-1)
          IW(IWPOSCB+2+KEEP(IXSZ)) = NBROW_STACK
          IF (COMPRESSCB) IW(IWPOSCB+1+XXS) = S_CB1COMP
      ELSE
        PTRAST(STEP(INODE)) = IPTRLU+1_8
        IF (COMPRESSCB) IW(IWPOSCB+1+XXS)=S_CB1COMP
        IW(IWPOSCB+1+KEEP(IXSZ)) = NBCOL
        IW(IWPOSCB+2+KEEP(IXSZ)) = 0
        IW(IWPOSCB+3+KEEP(IXSZ)) = NBROW_STACK
        IW(IWPOSCB+4+KEEP(IXSZ)) = 0
        IW(IWPOSCB+5+KEEP(IXSZ)) = 1
        IW(IWPOSCB+6+KEEP(IXSZ)) = 0
        IOLDP1   = PTLUST_S(STEP(INODE))+H_INODE
        PTROWEND = IWPOSCB+6+NBROW_STACK+KEEP(IXSZ)
        DO I = 1, NBROW_STACK
          IW(IWPOSCB+7+KEEP(IXSZ)+I-1) =
     &    IW(IOLDP1+NFRONT-NBROW_STACK+I-1)
        ENDDO
        DO I = 1, NBCOL
          IW(PTROWEND+I)=IW(IOLDP1+NFRONT+NPIV+I-1)
        ENDDO
      END IF 
      IF ( KEEP(50).NE.0 .AND. TYPE .EQ. 1 
     &     .AND. MUST_COMPACT_FACTORS ) THEN
        POSELT = PTRFAC(STEP(INODE))
        CALL SMUMPS_324(A(POSELT), LDA, 
     &                         NPIV, NBROW, KEEP(50))
        MUST_COMPACT_FACTORS = .FALSE.
      ENDIF
      IF (  KEEP(50).EQ.0 .AND. MUST_COMPACT_FACTORS )
     &     THEN
        LAST_ALLOWED_POS = POSELT + int(LDA,8)*int(NPIV+NBROW-1,8)
     &                     + int(NPIV,8)
      ELSE
        LAST_ALLOWED_POS = -1_8
      ENDIF
      NCBROW_ALREADY_MOVED = 0
 10   CONTINUE
      NCBROW_PREVIOUSLY_MOVED = NCBROW_ALREADY_MOVED
      IF (IPTRLU .LT. POSFAC ) THEN
        CALL SMUMPS_652( A, LA, LDA,
     &  POSELT, IPTRLU, NPIV, NBCOL_STACK, NBROW_STACK,
     &  NBROW_SEND, LREQCB, KEEP, COMPRESSCB,
     &  LAST_ALLOWED_POS, NCBROW_ALREADY_MOVED )
      ELSE
        CALL SMUMPS_705( A, LA, LDA,
     &  POSELT, IPTRLU, NPIV, NBCOL_STACK, NBROW_STACK,
     &  NBROW_SEND, LREQCB, KEEP, COMPRESSCB )
        NCBROW_ALREADY_MOVED = NBROW_STACK
      ENDIF
      IF (LAST_ALLOWED_POS .NE. -1_8) THEN
      MUST_COMPACT_FACTORS =.FALSE.
        IF ( NCBROW_ALREADY_MOVED .EQ. NBROW_STACK ) THEN
          NCBROW_ALREADY_MOVED = NCBROW_ALREADY_MOVED + NBROW_SEND
        ENDIF
        NCBROW_NEWLY_MOVED = NCBROW_ALREADY_MOVED
     &                    - NCBROW_PREVIOUSLY_MOVED
        FACTOR_POS = POSELT +
     &         int(LDA,8)*int(NPIV+NBROW-NCBROW_ALREADY_MOVED,8)
        CALL SMUMPS_651( A(FACTOR_POS), LDA, NPIV,
     &       NCBROW_NEWLY_MOVED )
        INEW = FACTOR_POS + int(NPIV,8) * int(NCBROW_NEWLY_MOVED,8)
        IOLD = INEW + int(NCBROW_NEWLY_MOVED,8) * int(NBCOL_STACK,8)
        DO I = 1, NCBROW_PREVIOUSLY_MOVED*NPIV
            A(INEW) = A(IOLD)
            IOLD = IOLD + 1_8
            INEW = INEW + 1_8
        ENDDO
        KEEP8(8)=KEEP8(8) + int(NCBROW_PREVIOUSLY_MOVED,8)
     &                    * int(NPIV,8)
        LAST_ALLOWED_POS = INEW
        IF (NCBROW_ALREADY_MOVED.LT.NBROW_STACK) THEN
          GOTO 10
        ENDIF
      ENDIF
 190  CONTINUE
      IF (MUST_COMPACT_FACTORS) THEN
       POSELT = PTRFAC(STEP(INODE))
       CALL SMUMPS_324(A(POSELT), LDA, 
     &                         NPIV, NBROW, KEEP(50))
       MUST_COMPACT_FACTORS = .FALSE.
      ENDIF
      IOLDPS = PTLUST_S(STEP(INODE))
      IW(IOLDPS+KEEP(IXSZ))     = NBCOL
      IW(IOLDPS + 1+KEEP(IXSZ)) = NASS - NPIV
      IF (TYPE.EQ.2) THEN
        IW(IOLDPS + 2+KEEP(IXSZ)) = NASS
      ELSE
        IW(IOLDPS + 2+KEEP(IXSZ)) = NFRONT
      ENDIF
      IW(IOLDPS + 3+KEEP(IXSZ)) = NPIV
      IF (INPLACE) THEN
        SIZE_INPLACE = LREQCB - MIN_SPACE_IN_PLACE
      ELSE
        SIZE_INPLACE = 0_8
      ENDIF
      CALL SMUMPS_93(SIZE_INPLACE,MYID,N,IOLDPS,TYPE, IW, LIW, 
     &    A, LA, POSFAC, LRLU, LRLUS, 
     &    IWPOS, PTRAST, PTRFAC, STEP, KEEP,KEEP8, SSARBR,INODE,IERR)
      IF(IERR.LT.0)THEN
         IFLAG=IERR
         IERROR=0
         GOTO 600
      ENDIF
 500  CONTINUE
      RETURN
  600 CONTINUE
      IF (IFLAG .NE. -1) CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_254
      SUBROUTINE SMUMPS_142( id)
      USE SMUMPS_COMM_BUFFER
      USE SMUMPS_LOAD
      USE SMUMPS_OOC
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
#ifndef SUN_
      INTERFACE
      SUBROUTINE SMUMPS_27(id, ANORMINF, LSCAL)
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET :: id
      REAL, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      END SUBROUTINE SMUMPS_27
      END INTERFACE
#endif
      TYPE(SMUMPS_STRUC), TARGET :: id
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INCLUDE 'mumps_headers.h'
      INTEGER NSEND, NSEND_TOT, LDPTRAR, NELT
      INTEGER NLOCAL, NLOCAL_TOT, KEEP13_SAVE, ITMP
      INTEGER(8) K67
      INTEGER(8) ITMP8
      INTEGER  MUMPS_275
      EXTERNAL MUMPS_275
      INTEGER MP, LP, MPG, allocok
      LOGICAL PROK, PROKG, LSCAL
      INTEGER SMUMPS_LBUF, SMUMPS_LBUFR_BYTES, SMUMPS_LBUF_INT
      INTEGER PTRIST, PTRWB, MAXELT_SIZE,
     &     ITLOC, IPOOL, NSTEPS, K28, LPOOL, LIW
      INTEGER IRANK, ID_ROOT
      INTEGER KKKK, NZ_locMAX
      INTEGER(8) MEMORY_MD_ARG
      INTEGER(8) MAXS_BASE8, MAXS_BASE_RELAXED8
      REAL CNTL4
      INTEGER MIN_PERLU, MAXIS_ESTIM
      INTEGER   MAXIS
      INTEGER(8) :: MAXS
      DOUBLE PRECISION TIME
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER PERLU, TOTAL_MBYTES, K231, K232, K233
      INTEGER COLOUR, COMM_FOR_SCALING 
      INTEGER LIWK, LWK, LWK_REAL
      LOGICAL I_AM_SLAVE, PERLU_ON, WK_USER_PROVIDED
      REAL :: ANORMINF, SEUIL, SEUIL_LDLT_NIV2
      REAL :: CNTL1, CNTL3, CNTL5, CNTL6, EPS
      INTEGER N, LPN_LIST,POSBUF
      INTEGER, DIMENSION (:), ALLOCATABLE :: ITMP2
      INTEGER I,K
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWK
      REAL, DIMENSION(:), ALLOCATABLE :: WK
      REAL, DIMENSION(:), ALLOCATABLE :: WK_REAL
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IWK8
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURS
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCS
      INTEGER BUREGISTRE(12)
      INTEGER BUINTSZ, BURESZ, BUJOB
      INTEGER BUMAXMN, M, SCMYID, SCNPROCS
      REAL    SCONEERR, SCINFERR
      INTEGER, POINTER ::  JOB, NZ
      REAL,DIMENSION(:),POINTER::RINFO, RINFOG
      REAL,DIMENSION(:),POINTER::    CNTL
      INTEGER,DIMENSION(:),POINTER::INFO, INFOG, KEEP
      INTEGER, DIMENSION(:), POINTER :: MYIRN_loc, MYJCN_loc
      REAL, DIMENSION(:), POINTER :: MYA_loc
      INTEGER, TARGET :: DUMMYIRN_loc(1), DUMMYJCN_loc(1)
      REAL, TARGET :: DUMMYA_loc(1)
      INTEGER(8),DIMENSION(:),POINTER::KEEP8
      INTEGER,DIMENSION(:),POINTER::ICNTL
      EXTERNAL SMUMPS_505
      INTEGER SMUMPS_505
      INTEGER(8) TOTAL_BYTES
      INTEGER(8) :: I8TMP
      INTEGER numroc
      EXTERNAL numroc
      REAL, DIMENSION(:), POINTER :: RHS_MUMPS
      LOGICAL :: RHS_MUMPS_ALLOCATED
      JOB=>id%JOB
      NZ=>id%NZ
      RINFO=>id%RINFO
      RINFOG=>id%RINFOG
      CNTL=>id%CNTL
      INFO=>id%INFO
      INFOG=>id%INFOG
      KEEP=>id%KEEP
      KEEP8=>id%KEEP8
      ICNTL=>id%ICNTL
      IF (id%NZ_loc .NE. 0) THEN
        MYIRN_loc=>id%IRN_loc
        MYJCN_loc=>id%JCN_loc
        MYA_loc=>id%A_loc
      ELSE
        MYIRN_loc=>DUMMYIRN_loc
        MYJCN_loc=>DUMMYJCN_loc
        MYA_loc=>DUMMYA_loc
      ENDIF
      N = id%N
      EPS = epsilon ( ZERO )
      NULLIFY(RHS_MUMPS)
      RHS_MUMPS_ALLOCATED = .FALSE.
      IF (KEEP8(24).GT.0_8) THEN 
           NULLIFY(id%S)
      ENDIF
      WK_USER_PROVIDED = (id%LWK_USER.NE.0)
      IF (WK_USER_PROVIDED) THEN
          IF (id%LWK_USER.GT.0) THEN
            KEEP8(24) = int(id%LWK_USER,8)
          ELSE
            KEEP8(24) = -int(id%LWK_USER,8)* 1000000_8 
          ENDIF
      ELSE
          KEEP8(24) = 0_8
      ENDIF
      KEEP13_SAVE = KEEP(13)
      id%DKEEP(4)=-1.0E0
      id%DKEEP(5)=-1.0E0
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LP  = ICNTL( 1 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IF ( PROK ) WRITE( MP, 130 )
      IF ( PROKG ) WRITE( MPG, 130 )
      IF ( PROKG .and. KEEP(53).GT.0 ) THEN
        WRITE(MPG,'(/A,I3)') ' Null space option :', KEEP(19)
        IF ( KEEP(21) .ne. 0 ) THEN 
          WRITE( MPG, '(A,I10)') ' Max deficiency    : ', KEEP(21)
        END IF
        IF ( KEEP(22) .ne. 0 ) THEN 
          WRITE( MPG, '(A,I10)') ' Min deficiency    : ', KEEP(22)
        END IF
      END IF
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
        IF (id%MYID .EQ. MASTER .AND. KEEP(201) .NE. -1) THEN
          KEEP(201)=id%ICNTL(22)
          IF (KEEP(201) .NE. 0) THEN 
#           if defined(OLD_OOC_NOPANEL)
              KEEP(201)=2
#           else
              KEEP(201)=1
#           endif
          ENDIF
        ENDIF
        CALL MPI_BCAST( KEEP(12), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(19), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(21), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (id%MYID.EQ.MASTER) THEN
          IF (KEEP(217).GT.2.OR.KEEP(217).LT.0) THEN
            KEEP(217)=0
          ENDIF
          KEEP(214)=KEEP(217)
          IF (KEEP(214).EQ.0) THEN
            IF (KEEP(201).NE.0) THEN 
              KEEP(214)=1
            ELSE
              KEEP(214)=2
            ENDIF
          ENDIF
        ENDIF
        CALL MPI_BCAST( KEEP(214), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (KEEP(201).NE.0) THEN
          CALL MPI_BCAST( KEEP(99), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(205), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(211), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        ENDIF
        IF ( KEEP(50) .eq. 1 ) THEN
          IF (id%CNTL(1) .ne. ZERO ) THEN
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)')
     &' ** Warning : SPD solver called, resetting CNTL(1) to 0.0E0'
            END IF
          END IF
          id%CNTL(1) = ZERO
        END IF
      IF (KEEP(219).NE.0) THEN
       CALL SMUMPS_617(max(KEEP(108),1),IERR)
       IF (IERR .NE. 0) THEN
          INFO(1) = -13
          INFO(2) = max(KEEP(108),1)
       END IF
      ENDIF
      IF (id%KEEP(252).EQ.1 .AND. id%MYID.EQ.MASTER) THEN
        IF (id%ICNTL(20).EQ.1) THEN 
          id%INFO(1)=-43
          id%INFO(2)=20
          IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Sparse RHS is incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(30).NE.0) THEN 
          id%INFO(1)=-43
          id%INFO(2)=30
          IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(9) .NE. 1) THEN
          id%INFO(1)=-43
          id%INFO(2)=9
          IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: sparse RHS incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ENDIF
      ENDIF
      CALL MUMPS_276( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF (INFO(1).LT.0) GOTO 530
      IF ( PROKG ) THEN
          WRITE( MPG, 172 ) id%NSLAVES, id%ICNTL(22),
     &    KEEP8(111), KEEP(126), KEEP(127), KEEP(28)
          IF (KEEP(252).GT.0) 
     &    WRITE(MPG,173) KEEP(253)
      ENDIF
      IF (KEEP(201).LE.0) THEN
        KEEP(IXSZ)=XSIZE_IC
      ELSE IF (KEEP(201).EQ.2) THEN
        KEEP(IXSZ)=XSIZE_OOC_NOPANEL
      ELSE IF (KEEP(201).EQ.1) THEN
        IF (KEEP(50).EQ.0) THEN
          KEEP(IXSZ)=XSIZE_OOC_UNSYM
        ELSE
          KEEP(IXSZ)=XSIZE_OOC_SYM
        ENDIF
      ENDIF
      IF (id%MYID.EQ.MASTER) KEEP(258)=ICNTL(33)
      CALL MPI_BCAST(KEEP(258), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      IF (KEEP(258) .NE. 0) THEN
        KEEP(259) = 0      
        KEEP(260) = 1      
        id%DKEEP(6)  = 1.0E0  
      ENDIF
      CALL MPI_BCAST(KEEP(52), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      LSCAL = ((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 8))
      IF (LSCAL) THEN
        IF ( id%MYID.EQ.MASTER ) THEN
        ENDIF
        IF (KEEP(52) .EQ. 7) THEN 
           K231= KEEP(231)
           K232= KEEP(232)
           K233= KEEP(233)
        ELSEIF (KEEP(52) .EQ. 8) THEN
           K231= KEEP(239)
           K232= KEEP(240)
           K233= KEEP(241)
        ENDIF
        CALL MPI_BCAST(id%DKEEP(3),1,MPI_REAL,MASTER,
     &       id%COMM,IERR)
        IF ( ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) .AND. 
     &       KEEP(54).NE.0 ) THEN
           IF ( id%MYID .NE. MASTER ) THEN
              IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
              IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
         ENDIF
         M = N
         BUMAXMN=M
         IF(N > BUMAXMN) BUMAXMN = N
         LIWK = 4*BUMAXMN
         ALLOCATE (IWK(LIWK),BURP(M),BUCP(N),
     &            BURS(2* (id%NPROCS)),BUCS(2* (id%NPROCS)),
     &            stat=allocok)
         IF (allocok > 0) THEN
            INFO(1)=-13
            INFO(2)=LIWK+M+N+4* (id%NPROCS)
         ENDIF
         CALL MUMPS_276( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF (INFO(1).LT.0) GOTO 530
         BUJOB = 1
         LWK_REAL   = 1
         ALLOCATE(WK_REAL(LWK_REAL))
         CALL SMUMPS_693(
     &        MYIRN_loc(1), MYJCN_loc(1), MYA_loc(1),
     &        id%NZ_loc,
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         IF(LIWK < BUINTSZ) THEN
            DEALLOCATE(IWK)
            LIWK = BUINTSZ
            ALLOCATE(IWK(LIWK), stat=allocok)
            IF (allocok > 0) THEN
               INFO(1)=-13
               INFO(2)=LIWK
            ENDIF
         ENDIF
         LWK_REAL = BURESZ
         DEALLOCATE(WK_REAL)
         ALLOCATE (WK_REAL(LWK_REAL), stat=allocok)
         IF (allocok > 0) THEN
            INFO(1)=-13
            INFO(2)=LWK_REAL
         ENDIF
         CALL MUMPS_276( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF (INFO(1).LT.0) GOTO 530
         BUJOB = 2   
         CALL SMUMPS_693(
     &        MYIRN_loc(1), MYJCN_loc(1), MYA_loc(1),
     &        id%NZ_loc,
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         id%DKEEP(4) = SCONEERR
         id%DKEEP(5) = SCINFERR
         DEALLOCATE(IWK, WK_REAL,BURP,BUCP,BURS, BUCS)
        ELSE IF ( KEEP(54) .EQ. 0 ) THEN
          IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8))  THEN
              IF (id%MYID.EQ.MASTER) THEN
                COLOUR = 0
              ELSE
                COLOUR = MPI_UNDEFINED
              ENDIF
              CALL MPI_COMM_SPLIT( id%COMM, COLOUR, 0,
     &             COMM_FOR_SCALING, IERR )
              IF (id%MYID.EQ.MASTER) THEN
                 M = N
                 BUMAXMN=N
                 IF(N > BUMAXMN) BUMAXMN = N
                 LIWK = 1
                 ALLOCATE (IWK(LIWK),BURP(1),BUCP(1),
     &                BURS(1),BUCS(1),
     &                stat=allocok)
                 LWK_REAL = M + N  
                 ALLOCATE (WK_REAL(LWK_REAL), stat=allocok)
                 IF (allocok > 0) THEN
                    INFO(1)=-13
                    INFO(2)=1
                 ENDIF
                 IF (INFO(1) .LT. 0) GOTO 400
                 CALL MPI_COMM_RANK(COMM_FOR_SCALING, SCMYID, IERR)
                 CALL MPI_COMM_SIZE(COMM_FOR_SCALING, SCNPROCS, IERR)
                 BUJOB = 1
                 CALL SMUMPS_693(
     &                id%IRN(1), id%JCN(1), id%A(1),
     &                id%NZ,
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 IF(LWK_REAL < BURESZ) THEN
                    INFO(1) = -136
                    GOTO 400
                 ENDIF
                 BUJOB = 2
                 CALL SMUMPS_693(id%IRN(1),
     &                id%JCN(1), id%A(1),
     &                id%NZ,
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 id%DKEEP(4) = SCONEERR
                 id%DKEEP(5) = SCINFERR
                 DEALLOCATE(WK_REAL)                 
                 DEALLOCATE (IWK,BURP,BUCP,
     &                BURS,BUCS)
              ENDIF
              CALL MPI_BCAST( id%DKEEP(4),2,MPI_REAL,
     &                        MASTER, id%COMM, IERR )
  400         CONTINUE
              IF (id%MYID.EQ.MASTER) THEN
                CALL MPI_COMM_FREE(COMM_FOR_SCALING, IERR)
              ENDIF
              CALL MUMPS_276(ICNTL(1), INFO(1), id%COMM, id%MYID)
              IF (INFO(1).LT.0) GOTO 530
          ELSE IF (id%MYID.EQ.MASTER) THEN
            IF (KEEP(52).GT.0 .AND. KEEP(52).LE.6) THEN
              IF ( KEEP(52) .eq. 5 .or. 
     &          KEEP(52) .eq. 6 ) THEN
                LWK = NZ
              ELSE
                LWK = 1
              END IF
              LWK_REAL = 5 * N
              ALLOCATE( WK_REAL( LWK_REAL ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                INFO(1) = -13
                INFO(2) = LWK_REAL
                GOTO 137
              END IF
              ALLOCATE( WK( LWK ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                INFO(1) = -13
                INFO(2) = LWK
                GOTO 137
              END IF
              CALL SMUMPS_217(N, NZ, KEEP(52), id%A(1),
     &             id%IRN(1), id%JCN(1),
     &             id%COLSCA(1), id%ROWSCA(1),
     &             WK, LWK, WK_REAL, LWK_REAL, ICNTL(1), INFO(1) )
              DEALLOCATE( WK_REAL )
              DEALLOCATE( WK )
            ENDIF
          ENDIF
        ENDIF 
        IF (id%MYID.EQ.MASTER) THEN
          IF (PROKG.AND.(KEEP(52).EQ.7.OR.KEEP(52).EQ.8) 
     &             .AND. (K233+K231+K232).GT.0) THEN
           IF (K232.GT.0) WRITE(MPG, 166) id%DKEEP(4)
          ENDIF
        ENDIF
      ENDIF 
        LSCAL = (LSCAL .OR. (KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        IF (LSCAL .AND. KEEP(258).NE.0 .AND. id%MYID .EQ. MASTER) THEN
          DO I = 1, id%N
            CALL SMUMPS_761(id%ROWSCA(I),
     &           id%DKEEP(6),    
     &           KEEP(259))   
          ENDDO
          IF (KEEP(50) .EQ. 0) THEN 
            DO I = 1, id%N
              CALL SMUMPS_761(id%COLSCA(I),
     &           id%DKEEP(6),    
     &           KEEP(259))   
            ENDDO
          ELSE
            CALL SMUMPS_765(id%DKEEP(6), KEEP(259))
          ENDIF
          CALL SMUMPS_766(id%DKEEP(6), KEEP(259))
        ENDIF
 137  CONTINUE
      IF (associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE (id%root%RHS_CNTR_MASTER_ROOT)
        NULLIFY (id%root%RHS_CNTR_MASTER_ROOT)
      ENDIF
      IF ( id%MYID.EQ.MASTER.AND. KEEP(252).EQ.1 .AND.
     &      id%NRHS .NE. id%KEEP(253) ) THEN
          id%INFO(1)=-42
          id%INFO(2)=id%KEEP(253)
      ENDIF
      IF (id%KEEP(252) .EQ. 1) THEN
          IF ( id%MYID.NE.MASTER ) THEN
            id%KEEP(254) = N              
            id%KEEP(255) = N*id%KEEP(253) 
            ALLOCATE(RHS_MUMPS(id%KEEP(255)),stat=IERR)
            IF (IERR > 0) THEN
               INFO(1)=-13
               INFO(2)=id%KEEP(255)
               IF (LP > 0)
     &         WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
               NULLIFY(RHS_MUMPS)
            ENDIF
            RHS_MUMPS_ALLOCATED = .TRUE.
          ELSE 
            id%KEEP(254)=id%LRHS              
            id%KEEP(255)=id%LRHS*(id%KEEP(253)-1)+id%N 
            RHS_MUMPS=>id%RHS
            RHS_MUMPS_ALLOCATED = .FALSE.
            IF (LSCAL) THEN
              DO K=1, id%KEEP(253)
                DO I=1, N
                  RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          = RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          * id%ROWSCA(I)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          DO I= 1, id%KEEP(253)
            CALL MPI_BCAST(RHS_MUMPS((I-1)*id%KEEP(254)+1), N,
     &           MPI_REAL, MASTER,id%COMM,IERR)
          END DO
      ELSE
          id%KEEP(255)=1
          ALLOCATE(RHS_MUMPS(1))
          RHS_MUMPS_ALLOCATED = .TRUE.
      ENDIF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 530
      KEEP(110)=ICNTL(24)
      CALL MPI_BCAST(KEEP(110), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      IF (KEEP(110).NE.1) KEEP(110)=0
      IF (id%MYID .EQ. MASTER) CNTL3 = id%CNTL(3)
      CALL MPI_BCAST(CNTL3, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL5 = id%CNTL(5)
      CALL MPI_BCAST(CNTL5, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL6 = id%CNTL(6)
      CALL MPI_BCAST(CNTL6, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL1 = id%CNTL(1)
      CALL MPI_BCAST(CNTL1, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      ANORMINF = ZERO
      IF (KEEP(19).EQ.0) THEN 
         SEUIL = ZERO
      ELSE
         CALL SMUMPS_27(  id , ANORMINF, LSCAL )
         IF (CNTL6 .LT. ZERO) THEN
           SEUIL = EPS*ANORMINF
         ELSE
           SEUIL = CNTL6*ANORMINF
         ENDIF
         IF (PROKG) WRITE(MPG,*)
     &   ' ABSOLUTE PIVOT THRESHOLD for rank revealing =',SEUIL
      ENDIF
      SEUIL_LDLT_NIV2 = SEUIL
      IF (KEEP(110).EQ.0) THEN 
         id%DKEEP(1) = -1.0E0
         id%DKEEP(2) = ZERO
      ELSE
         IF (ANORMINF.EQ.ZERO) 
     &       CALL SMUMPS_27(  id , ANORMINF, LSCAL )
         IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(1)  = abs(CNTL(3))
         ELSE IF  (CNTL3 .GT. ZERO) THEN
           id%DKEEP(1)  = CNTL3*ANORMINF
         ELSE 
           id%DKEEP(1)  = 1.0E-5*EPS*ANORMINF
         ENDIF
         IF (PROKG) WRITE(MPG,*)
     &    ' ZERO PIVOT DETECTION ON, THRESHOLD          =',id%DKEEP(1)
         IF (CNTL5.GT.ZERO) THEN
            id%DKEEP(2) = CNTL5 * ANORMINF
            IF (PROKG) WRITE(MPG,*) 
     &    ' FIXATION FOR NULL PIVOTS                    =',id%DKEEP(2)
         ELSE
            IF (PROKG) WRITE(MPG,*) 'INFINITE FIXATION '
            IF (id%KEEP(50).EQ.0) THEN
             id%DKEEP(2) = -max(1.0E10*ANORMINF, 
     &                sqrt(huge(ANORMINF))/1.0E8)
            ELSE
            id%DKEEP(2) = ZERO
            ENDIF
         ENDIF
      ENDIF
      IF (KEEP(53).NE.0) THEN
        ID_ROOT =MUMPS_275(id%PROCNODE_STEPS(id%STEP(KEEP(20))),
     &                          id%NSLAVES)
        IF ( KEEP( 46 )  .NE. 1 ) THEN
          ID_ROOT = ID_ROOT + 1
        END IF
      ENDIF
      IF ( associated( id%PIVNUL_LIST) ) DEALLOCATE(id%PIVNUL_LIST)
      IF(KEEP(110) .EQ. 1) THEN
         LPN_LIST = N
      ELSE
         LPN_LIST = 1
      ENDIF
      IF (KEEP(19).NE.0 .AND.
     &   (ID_ROOT.EQ.id%MYID .OR. id%MYID.EQ.MASTER)) THEN
         LPN_LIST = N
      ENDIF
      ALLOCATE( id%PIVNUL_LIST(LPN_LIST),stat = IERR )
      IF ( IERR .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=LPN_LIST
      END IF
      id%PIVNUL_LIST(1:LPN_LIST) = 0
      KEEP(109) = 0
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 530
      IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
        IF (id%MYID .EQ. MASTER) CNTL4 = id%CNTL(4)
        CALL MPI_BCAST( CNTL4, 1, MPI_REAL,
     &                MASTER, id%COMM, IERR )
        IF ( CNTL4 .GE. ZERO ) THEN
         KEEP(97) = 1
         IF ( CNTL4 .EQ. ZERO ) THEN
            IF(ANORMINF .EQ. ZERO) THEN
               CALL SMUMPS_27(  id , ANORMINF, LSCAL )
            ENDIF
            SEUIL = sqrt(EPS) * ANORMINF
         ELSE
            SEUIL = CNTL4
         ENDIF
         SEUIL_LDLT_NIV2 = SEUIL
        ELSE 
         SEUIL = ZERO
        ENDIF
      ENDIF
      KEEP(98)  = 0
      KEEP(103) = 0
      KEEP(105) = 0
      MAXS      = 1_8
      IF ( id%MYID.EQ.MASTER ) THEN
        ITMP = ICNTL(23)
      END IF
      CALL MPI_BCAST( ITMP, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      IF (WK_USER_PROVIDED) ITMP = 0
      ITMP8 = int(ITMP, 8)
      KEEP8(4) = ITMP8 * 1000000_8   
      PERLU = KEEP(12)
      IF (KEEP(201) .EQ. 0) THEN
        MAXS_BASE8=KEEP8(12)
       ELSE
        MAXS_BASE8=KEEP8(14)
      ENDIF
      IF (WK_USER_PROVIDED) THEN
        MAXS = KEEP8(24)
      ELSE
       IF ( MAXS_BASE8 .GT. 0_8 ) THEN
          MAXS_BASE_RELAXED8 =
     &         MAXS_BASE8 + int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8)
          IF (MAXS_BASE_RELAXED8 > huge(MAXS)) THEN
            WRITE(*,*) "Internal error: I8 overflow"
            CALL MUMPS_ABORT()
          ENDIF
          MAXS_BASE_RELAXED8 = max(MAXS_BASE_RELAXED8, 1_8)
          MAXS = MAXS_BASE_RELAXED8
       ELSE
        MAXS = 1_8
        MAXS_BASE_RELAXED8 = 1_8
       END IF
      ENDIF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (INFO(1) .LT. 0) THEN
        GOTO 530
      ENDIF
      IF ((.NOT.WK_USER_PROVIDED).AND.(I_AM_SLAVE)) THEN
          IF (KEEP(96).GT.0) THEN
            MAXS=int(KEEP(96),8)
          ELSE
            IF (KEEP8(4) .NE. 0_8) THEN
              PERLU_ON = .TRUE.
              CALL SMUMPS_214( id%KEEP(1), id%KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%LNA, id%NZ, id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE., KEEP(201), 
     &        PERLU_ON, TOTAL_BYTES)
              MAXS_BASE_RELAXED8=MAXS_BASE_RELAXED8 +
     &        (KEEP8(4)-TOTAL_BYTES)/int(KEEP(35),8) 
              IF (MAXS_BASE_RELAXED8 > int(huge(MAXS),8)) THEN
                WRITE(*,*) "Internal error: I8 overflow"
                CALL MUMPS_ABORT()
              ELSE IF (MAXS_BASE_RELAXED8 .LE. 0_8) THEN
                id%INFO(1)=-9
                IF ( -MAXS_BASE_RELAXED8 .GT.
     &               int(huge(id%INFO(1)),8) ) THEN
                  WRITE(*,*) "I8: OVERFLOW"
                  CALL MUMPS_ABORT()
                ENDIF
                id%INFO(2)=-int(MAXS_BASE_RELAXED8,4)
              ELSE
                MAXS=MAXS_BASE_RELAXED8
              ENDIF
            ENDIF
          ENDIF
      ENDIF 
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (INFO(1) .LT. 0) THEN
        GOTO 530
      ENDIF
      CALL SMUMPS_713(PROKG, MPG, MAXS, id%NSLAVES,
     & id%COMM, "effective relaxed size of S              =")
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 530
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        CALL SMUMPS_188( dble(id%COST_SUBTREES),
     &        KEEP(64), KEEP(66),MAXS )
        K28=KEEP(28)
        MEMORY_MD_ARG = min(int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8 ),
     &                      max(0_8, MAXS-MAXS_BASE8))
        CALL SMUMPS_185( id, MEMORY_MD_ARG, MAXS )
        CALL SMUMPS_587(id, IERR)
        IF (IERR < 0) THEN
          INFO(1) = -90
          INFO(2) = 0
          GOTO 112
        ENDIF
        IF (KEEP(201) .GT. 0) THEN
           IF (KEEP(201).EQ.1 
     &         .AND.KEEP(50).EQ.0 
     &         .AND.KEEP(251).NE.2 
     &         ) THEN
             OOC_NB_FILE_TYPE=2 
           ELSE
             OOC_NB_FILE_TYPE=1 
           ENDIF
           IF (KEEP(205) .GT. 0) THEN
             KEEP(100) = KEEP(205)
           ELSE
             IF (KEEP(201).EQ.1) THEN 
               I8TMP = int(OOC_NB_FILE_TYPE,8) * 2_8 * int(KEEP(226),8)
             ELSE
               I8TMP = 2_8 * KEEP8(119)
             ENDIF
             I8TMP = I8TMP +  int(max(KEEP(12),0),8) *
     &               (I8TMP/100_8+1_8)
             I8TMP = min(I8TMP, 12000000_8)
             KEEP(100)=int(I8TMP)
           ENDIF
           IF (KEEP(201).EQ.1) THEN
             IF ( KEEP(99) < 3 ) THEN
               KEEP(99) = KEEP(99) + 3
             ENDIF
             IF (id%MYID_NODES .eq. MASTER) THEN
               write(6,*) ' PANEL: INIT and force STRAT_IO= ',
     &         id%KEEP(99)
             ENDIF
           ENDIF
           IF (KEEP(99) .LT.3) KEEP(100)=0
           IF((dble(KEEP(100))*dble(KEEP(35))/dble(2)).GT.
     &       (dble(1999999999)))THEN
             IF (PROKG) THEN
               WRITE(MPG,*)id%MYID,': Warning: DIM_BUF_IO might be
     &  too big for Filesystem'
             ENDIF
           ENDIF
           ALLOCATE (id%OOC_INODE_SEQUENCE(KEEP(28),
     &          OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_INODE_SEQUENCE)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_TOTAL_NB_NODES(OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = OOC_NB_FILE_TYPE
              NULLIFY(id%OOC_TOTAL_NB_NODES)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_SIZE_OF_BLOCK(KEEP(28),
     &          OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_SIZE_OF_BLOCK)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_VADDR(KEEP(28),OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_VADDR)
              GOTO 112
           ENDIF
        ENDIF
      ENDIF
 112  CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (INFO(1) < 0) THEN
        GOTO 513
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (KEEP(201) .GT. 0) THEN
           IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
             CALL SMUMPS_575(id,MAXS)
           ELSE
             WRITE(*,*) "Internal error in SMUMPS_142"
             CALL MUMPS_ABORT()
           ENDIF
           IF(INFO(1).LT.0)THEN
              GOTO 111
           ENDIF
        ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
        CALL SMUMPS_190(0,.FALSE.,dble(id%COST_SUBTREES),
     &          id%KEEP(1),id%KEEP8(1))
#endif
        IF (INFO(1).LT.0) GOTO 111
#if defined(stephinfo)
        write(*,*) 'proc ',id%MYID,' array of dist : ',
     &       id%MEM_DIST(0:id%NSLAVES - 1)
#endif
      END IF
      IF ( associated (id%S) ) THEN
        DEALLOCATE(id%S)
        NULLIFY(id%S)
        KEEP8(23)=0_8  
      ENDIF
#if defined (LARGEMATRICES)
      IF ( id%MYID .ne. MASTER ) THEN
#endif
      IF (.NOT.WK_USER_PROVIDED) THEN
        ALLOCATE (id%S(MAXS),stat=IERR)
        KEEP8(23) = MAXS
        IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          CALL MUMPS_735(MAXS, INFO(2))
          NULLIFY(id%S)
          KEEP8(23)=0_8
        ENDIF
      ELSE
       id%S => id%WK_USER(1:KEEP8(24))
      ENDIF
#if defined (LARGEMATRICES)
      END IF
#endif
 111  CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP(55) .eq. 0 ) THEN
        IF (associated( id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
        IF ( I_AM_SLAVE .and. KEEP(13) .ne. 0 ) THEN
          ALLOCATE( id%DBLARR( KEEP(13) ), stat = IERR )
        ELSE
          ALLOCATE( id%DBLARR( 1 ), stat =IERR )
        END IF
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     &       ':Error allocating DBLARR : IERR = ', IERR
          INFO(1)=-13
          INFO(2)=KEEP(13)
          NULLIFY(id%DBLARR)
          GOTO 100
        END IF
      ELSE
         IF ( associated( id%INTARR ) ) THEN
           DEALLOCATE( id%INTARR )
           NULLIFY( id%INTARR )
         END IF
         IF ( I_AM_SLAVE .and. KEEP(14) .ne. 0 ) THEN
           ALLOCATE( id%INTARR( KEEP(14) ), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = KEEP(14)
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         ELSE
           ALLOCATE( id%INTARR(1),stat=allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         END IF
         IF (associated( id%DBLARR)) THEN
           DEALLOCATE(id%DBLARR)
           NULLIFY(id%DBLARR)
         ENDIF
         IF ( I_AM_SLAVE ) THEN
           IF (      id%MYID_NODES .eq. MASTER
     &       .AND.   KEEP(46)   .eq. 1
     &       .AND.   KEEP(52)   .eq. 0 ) THEN
             id%DBLARR => id%A_ELT
           ELSE
             IF ( KEEP(13) .ne. 0 ) THEN
               ALLOCATE( id%DBLARR( KEEP(13) ), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = KEEP(13)
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF
             ELSE
               ALLOCATE( id%DBLARR(1), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = 1
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF 
             END IF
           END IF
         ELSE
           ALLOCATE( id%DBLARR(1), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%DBLARR)
             GOTO 100
           END IF
         END IF
      END IF
      IF ( KEEP(38).NE.0 .AND.  I_AM_SLAVE ) THEN
         CALL SMUMPS_165( id%N,
     &   id%root, id%FILS(1), KEEP(38), id%KEEP(1), id%INFO(1) )
      END IF
 100  CONTINUE
      CALL MUMPS_276( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP( 55 ) .eq. 0 ) THEN
      IF (KEEP(38).NE.0 .AND. I_AM_SLAVE) THEN
        LWK = numroc( id%root%ROOT_SIZE, id%root%MBLOCK,
     &             id%root%MYROW, 0, id%root%NPROW ) 
        LWK = max( 1, LWK )
        LWK = LWK*
     &        numroc( id%root%ROOT_SIZE, id%root%NBLOCK,
     &        id%root%MYCOL, 0, id%root%NPCOL )
        LWK = max( 1, LWK )
      ELSE
        LWK = 1
      ENDIF
      IF (MAXS .LT. int(LWK,8)) THEN
           INFO(1) = -9
           INFO(2) = LWK
      ENDIF
      CALL MUMPS_276( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP(54) .eq. 0 ) THEN
        IF ( id%MYID .eq. MASTER ) THEN
          ALLOCATE(IWK(id%N), stat=allocok)
          IF ( allocok .NE. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
          END IF
#if defined(LARGEMATRICES)
          IF ( associated (id%S) ) THEN
            DEALLOCATE(id%S)
            NULLIFY(id%S)
            KEEP8(23)=0_8
          ENDIF
          ALLOCATE (WK(LWK),stat=IERR)
          IF ( IERR .GT. 0 ) THEN
            INFO(1) = -13
            INFO(2) = LWK
            write(6,*) ' PB1 ALLOC LARGEMAT'
          ENDIF
#endif
        ENDIF
        CALL MUMPS_276( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
        IF ( INFO(1).LT.0 ) GOTO 500
        IF ( id%MYID .eq. MASTER ) THEN
          IF (PROKG ) THEN
            CALL MUMPS_291(TIME)
          END IF
          IF ( .not. associated( id%INTARR ) ) THEN
            ALLOCATE( id%INTARR( 1 ) )
          ENDIF
#if defined(LARGEMATRICES)
          CALL SMUMPS_148(id%N, NZ, id%A(1),
     &      id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &      LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &      id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &      min(KEEP(39),id%NZ),
     &      LP, id%COMM, id%root, KEEP,KEEP8,
     &      id%FILS(1), IWK(1), 
     &
     &      id%INTARR(1), id%DBLARR(1),
     &      id%PTRAR(1), id%PTRAR(id%N+1),
     &      id%FRERE_STEPS(1), id%STEP(1), WK(1), int(LWK,8),
     &      id%ISTEP_TO_INIV2, id%I_AM_CAND,
     &      id%CANDIDATES) 
          write(6,*) '!!! A,IRN,JCN are freed during facto '
          DEALLOCATE (id%A)
          NULLIFY(id%A)
          DEALLOCATE (id%IRN)
          NULLIFY (id%IRN)
          DEALLOCATE (id%JCN)
          NULLIFY (id%JCN)
          IF (.NOT.WK_USER_PROVIDED) THEN
            ALLOCATE (id%S(MAXS),stat=IERR)
            KEEP8(23) = MAXS
            IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = MAXS
              NULLIFY(id%S)
              KEEP8(23)=0_8
              write(6,*) ' PB2 ALLOC LARGEMAT',MAXS
              CALL MUMPS_ABORT()
            ENDIF
          ELSE
            id%S => id%WK_USER(1:KEEP8(24))
          ENDIF
          id%S(MAXS-LWK+1:MAXS) = WK(1:LWK)
          DEALLOCATE (WK)
#else
          CALL SMUMPS_148(id%N, NZ, id%A(1),
     &    id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &    LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &    id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &    min(KEEP(39),id%NZ),
     &    LP, id%COMM, id%root, KEEP(1),KEEP8(1),
     &    id%FILS(1), IWK(1),
     &
     &    id%INTARR(1), id%DBLARR(1),
     &    id%PTRAR(1), id%PTRAR(id%N+1),
     &    id%FRERE_STEPS(1), id%STEP(1), id%S(1), MAXS,
     &    id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &    id%CANDIDATES(1,1) ) 
#endif
          DEALLOCATE(IWK)
          IF ( PROKG ) THEN
            CALL MUMPS_292(TIME)
            WRITE(MPG,160) TIME
            CALL MUMPS_291(TIME)
          END IF
        ELSE
          CALL SMUMPS_145( id%N,
     &       id%DBLARR( 1 ), max(1,KEEP( 13 )),
     &       id%INTARR( 1 ), max(1,KEEP( 14 )),
     &       id%PTRAR( 1 ),
     &       id%PTRAR(id%N+1),
     &       KEEP( 1 ), KEEP8(1), id%MYID, id%COMM,
     &       min(id%KEEP(39),id%NZ),
     &
     &       id%S(1), MAXS,
     &       id%root,
     &       id%PROCNODE_STEPS(1), id%NSLAVES,
     &       id%SYM_PERM(1), id%FRERE_STEPS(1), id%STEP(1),
     &       id%INFO(1), id%INFO(2) )
        ENDIF
      ELSE
      IF (PROKG ) THEN
        CALL MUMPS_291(TIME)
      END IF
      IF ( I_AM_SLAVE ) THEN
       NZ_locMAX = 0
       CALL MPI_ALLREDUCE(id%NZ_loc, NZ_locMAX, 1, MPI_INTEGER, 
     &                   MPI_MAX, id%COMM_NODES, IERR)
        CALL SMUMPS_282( id%N,
     &  id%NZ_loc,
     &  id,
     &  id%DBLARR(1), KEEP(13), id%INTARR(1),
     &  KEEP(14), id%PTRAR(1), id%PTRAR(id%N+1),
     &  KEEP(1), KEEP8(1), id%MYID_NODES,
     &  id%COMM_NODES, min(id%KEEP(39),NZ_locMAX),
     &  id%S(1), MAXS, id%root, id%PROCNODE_STEPS(1),
     &  id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     &  id%ICNTL(1), id%INFO(1), NSEND, NLOCAL,
     &  id%ISTEP_TO_INIV2(1),
     &  id%CANDIDATES(1,1) )
        IF ( ( KEEP(52).EQ.7 ).OR. (KEEP(52).EQ.8) ) THEN
          IF ( id%MYID > 0 ) THEN
            IF (associated(id%ROWSCA)) THEN
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
            IF (associated(id%COLSCA)) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ENDIF
          ENDIF
        ENDIF
#if defined(LARGEMATRICES)
         IF (associated(id%IRN_loc)) THEN
            DEALLOCATE(id%IRN_loc)
            NULLIFY(id%IRN_loc)
         ENDIF
         IF (associated(id%JCN_loc)) THEN
            DEALLOCATE(id%JCN_loc)
            NULLIFY(id%JCN_loc)
         ENDIF
         IF (associated(id%A_loc)) THEN
            DEALLOCATE(id%A_loc)
            NULLIFY(id%A_loc)
         ENDIF
       write(6,*) ' Warning :', 
     &        ' id%A_loc, IRN_loc, JCN_loc deallocated !!! '
#endif
      IF (PROK) THEN
        WRITE(MP,120) NLOCAL, NSEND
      END IF
      END IF
      IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
        NSEND  = 0
        NLOCAL = 0
      END IF
      CALL MPI_REDUCE( NSEND, NSEND_TOT, 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( NLOCAL, NLOCAL_TOT, 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
        WRITE(MPG,125) NLOCAL_TOT, NSEND_TOT
      END IF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( PROKG ) THEN
          CALL MUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL MUMPS_291(TIME)
      END IF
      END IF
      ELSE
        IF (PROKG ) THEN
          CALL MUMPS_291(TIME)
        END IF
      IF ( id%MYID.eq.MASTER)
     &CALL SMUMPS_213( id%ELTPTR(1),
     &                        id%NELT,
     &                        MAXELT_SIZE )
      CALL SMUMPS_126( id%N, id%NELT, id%NA_ELT,
     &     id%COMM, id%MYID,
     &     id%NSLAVES, id%PTRAR(1),
     &     id%PTRAR(id%NELT+2),
     &     id%INTARR(1), id%DBLARR(1),
     &     id%KEEP(1), id%KEEP8(1), MAXELT_SIZE,
     &     id%FRTPTR(1), id%FRTELT(1),
     &     id%S(1), MAXS, id%FILS(1),
     &     id, id%root )
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( PROKG ) THEN
          CALL MUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL MUMPS_291(TIME)
      END IF
      END IF
      IF ( I_AM_SLAVE )  THEN
        CALL SMUMPS_528(id%MYID_NODES)
        SMUMPS_LBUFR_BYTES = KEEP( 44 ) * KEEP( 35 )
        SMUMPS_LBUFR_BYTES = max( SMUMPS_LBUFR_BYTES,
     &                      100000 )
        PERLU = KEEP( 12 )
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU=2
        ELSE
          MIN_PERLU=0
        ENDIF
        SMUMPS_LBUFR_BYTES = SMUMPS_LBUFR_BYTES
     &        + int( 2.0E0 * real(max(PERLU,MIN_PERLU))*
     &        real(SMUMPS_LBUFR_BYTES)/100E0)
        IF (KEEP(48)==5) THEN
           KEEP8(21) = KEEP8(22) + int( real(max(PERLU,MIN_PERLU))*
     &        real(KEEP8(22))/100E0,8)
        ENDIF
        SMUMPS_LBUF = int( real(KEEP(213)) / 100.0E0 *
     &                      real(KEEP(43)) * real(KEEP(35))  )
        SMUMPS_LBUF = max( SMUMPS_LBUF, 100000 )
        SMUMPS_LBUF = SMUMPS_LBUF
     &                 + int( 2.0E0 * real(max(PERLU,MIN_PERLU))*
     &                   real(SMUMPS_LBUF)/100E0)
        SMUMPS_LBUF = max(SMUMPS_LBUF, SMUMPS_LBUFR_BYTES+3*KEEP(34))
        IF(id%KEEP(48).EQ.4)THEN
           SMUMPS_LBUFR_BYTES=SMUMPS_LBUFR_BYTES*5
           SMUMPS_LBUF=SMUMPS_LBUF*5
        ENDIF
        SMUMPS_LBUF_INT = ( KEEP(56) + id%NSLAVES * id%NSLAVES ) * 5
     &               * KEEP(34)
        IF ( KEEP( 38 ) .NE. 0 ) THEN
          KKKK = MUMPS_275( id%PROCNODE_STEPS(id%STEP(KEEP(38))),
     &                           id%NSLAVES )
          IF ( KKKK .EQ. id%MYID_NODES ) THEN
             SMUMPS_LBUF_INT = SMUMPS_LBUF_INT + 
     &     10 *  
     &      2 * ( id%NE_STEPS(id%STEP(KEEP(38))) + 1 ) * id%NSLAVES
     &                      * KEEP(34)
          END IF
        END IF
        IF ( MP .GT. 0 ) THEN
          WRITE( MP, 9999 ) SMUMPS_LBUFR_BYTES,
     &                      SMUMPS_LBUF, SMUMPS_LBUF_INT
        END IF
 9999   FORMAT( /,' Allocated buffers',/,' ------------------',/,
     &  ' Size of reception buffer in bytes ...... = ', I10,
     &  /,
     &  ' Size of async. emission buffer (bytes).. = ', I10,/,
     &  ' Small emission buffer (bytes) .......... = ', I10)
        CALL SMUMPS_55( SMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     &   ':Error allocating small Send buffer:IERR='
     &   ,IERR
          INFO(1)= -13
          INFO(2)= (SMUMPS_LBUF_INT+KEEP(34)-1)/KEEP(34)
          GO TO 110
        END IF
        CALL SMUMPS_53( SMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating Send buffer:IERR='
     &   ,IERR
          INFO(1)= -13
          INFO(2)= (SMUMPS_LBUF+KEEP(34)-1)/KEEP(34)
          GO TO 110
        END IF
        id%LBUFR_BYTES = SMUMPS_LBUFR_BYTES
        id%LBUFR = (SMUMPS_LBUFR_BYTES+KEEP(34)-1)/KEEP(34)
        IF (associated(id%BUFR)) DEALLOCATE(id%BUFR)
        ALLOCATE( id%BUFR( id%LBUFR ),stat=IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating BUFR:IERR='
     &   ,IERR
          INFO(1)=-13
          INFO(2)=id%LBUFR
          NULLIFY(id%BUFR)
          GO TO 110
        END IF
        PERLU          = KEEP( 12 )
        IF (KEEP(201).GT.0) THEN
          MAXIS_ESTIM   = KEEP(225)
        ELSE
          MAXIS_ESTIM   = KEEP(15)
        ENDIF
        MAXIS = max( 1,
     &       MAXIS_ESTIM + 2 * max(PERLU,10) * 
     &          ( MAXIS_ESTIM / 100 + 1 )
     &  )
        IF (associated(id%IS)) DEALLOCATE( id%IS )
        ALLOCATE( id%IS( MAXIS  ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) id%MYID,':Error allocating IS:IERR=',IERR
         INFO(1)=-13
         INFO(2)=MAXIS
         NULLIFY(id%IS)
         GO TO 110
        END IF
        LIW = MAXIS
        IF (associated( id%PTLUST_S )) DEALLOCATE(id%PTLUST_S)
        ALLOCATE( id%PTLUST_S( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocatingPTLUST:IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=id%KEEP(28)
          NULLIFY(id%PTLUST_S)
          GOTO 100
        END IF
        IF (associated( id%PTRFAC )) DEALLOCATE(id%PTRFAC)
        ALLOCATE( id%PTRFAC( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocatingPTRFAC:IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=id%KEEP(28)
          NULLIFY(id%PTRFAC)
          GOTO 100
        END IF
        PTRIST = 1
        PTRWB  = PTRIST + id%KEEP(28)
        ITLOC  = PTRWB  + 3 * id%KEEP(28)
        IPOOL  = ITLOC  + id%N + id%KEEP(253)
        LPOOL  = SMUMPS_505(id%KEEP(1),id%KEEP8(1))
        ALLOCATE( IWK(  IPOOL + LPOOL - 1 ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating IWK : IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=IPOOL + LPOOL - 1
          GOTO 110
        END IF
        ALLOCATE(IWK8( 2 * id%KEEP(28)), stat = IERR)
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating IWK : IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=2 * id%KEEP(28)
          GOTO 110
        END IF
      ENDIF
 110  CONTINUE
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( I_AM_SLAVE )  THEN
        CALL SMUMPS_60( id%LBUFR_BYTES )
        IF (MP .GT. 0) THEN
          WRITE( MP, 170 ) MAXS, MAXIS, KEEP8(12), KEEP(15), KEEP(13),
     &    KEEP(14), KEEP8(11), KEEP(26), KEEP(27)
        ENDIF
      END IF
      PERLU_ON = .TRUE.
      CALL SMUMPS_214( id%KEEP(1), id%KEEP8(1),
     &     id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     &     id%NA_ELT,
     &     id%NSLAVES, TOTAL_MBYTES, .FALSE., id%KEEP(201),
     &     PERLU_ON, TOTAL_BYTES)
      id%INFO(16) = TOTAL_MBYTES
      IF ( MP .gt. 0 ) THEN
          WRITE(MP,'(A,I10) ')
     &    ' ** Space in MBYTES used during factorization  :',
     &                id%INFO(16)
      END IF
      CALL MUMPS_243( id%MYID, id%COMM,
     &                           id%INFO(16), id%INFOG(18), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Memory relaxation parameter ( ICNTL(14)  )            :',
     &  KEEP(12)
        WRITE( MPG,'(A,I10) ')
     &  ' ** Rank of processor needing largest memory in facto     :',
     &  IRANK
        WRITE( MPG,'(A,I10) ')
     &  ' ** Space in MBYTES used by this processor for facto      :',
     &  id%INFOG(18)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during facto    :',
     &  ( id%INFOG(19)-id%INFO(16) ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during facto    :',
     &  id%INFOG(19) / id%NSLAVES
        END IF
      END IF
      KEEP8(31)= 0_8
      KEEP8(10) = 0_8
      KEEP8(8)=0_8
      INFO(9:14)=0
      RINFO(2:3)=ZERO
      IF ( I_AM_SLAVE ) THEN
        IF ( KEEP(55) .eq. 0 ) THEN
          LDPTRAR = id%N
        ELSE
          LDPTRAR = id%NELT + 1
        END IF
        IF ( id%KEEP(55) .NE. 0 ) THEN
          NELT = id%NELT
        ELSE
          NELT = 1
        END IF
        CALL SMUMPS_244( id%N, NSTEPS, id%S(1),
     &      MAXS, id%IS( 1 ), LIW,
     &      id%SYM_PERM(1), id%NA(1), id%LNA, id%NE_STEPS(1),
     &      id%ND_STEPS(1), id%FILS(1), id%STEP(1),
     &      id%FRERE_STEPS(1), id%DAD_STEPS(1), id%CANDIDATES(1,1), 
     &      id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     &      id%PTRAR(1), LDPTRAR, IWK( PTRIST ),
     &      id%PTLUST_S(1), id%PTRFAC(1), IWK( PTRWB ),
     &      IWK8,
     &      IWK( ITLOC ), RHS_MUMPS(1), IWK( IPOOL ), LPOOL,
     &      CNTL1, ICNTL(1), INFO(1), RINFO(1), KEEP(1),KEEP8(1),
     &      id%PROCNODE_STEPS(1),
     &      id%NSLAVES, id%COMM_NODES,
     &      id%MYID, id%MYID_NODES,
     &      id%BUFR(1),id%LBUFR,id%LBUFR_BYTES,
     &      id%INTARR(1), id%DBLARR(1), id%root,
     &      NELT, id%FRTPTR(1), 
     &      id%FRTELT(1), id%COMM_LOAD, id%ASS_IRECV, SEUIL,
     &      SEUIL_LDLT_NIV2, id%MEM_DIST(0),
     &       id%DKEEP(1),id%PIVNUL_LIST(1),LPN_LIST)
        IF ( MP . GT. 0 .and. KEEP(38) .ne. 0 ) THEN
          WRITE( MP, 175 ) KEEP(49)
        END IF
        DEALLOCATE( IWK  )
        DEALLOCATE( IWK8 )
      ENDIF
        IF ( KEEP(55) .eq. 0 ) THEN
          IF (associated( id%DBLARR)) THEN
            DEALLOCATE(id%DBLARR)
            NULLIFY(id%DBLARR)
          ENDIF
        ELSE
          DEALLOCATE( id%INTARR)
          NULLIFY( id%INTARR )
          IF (      id%MYID_NODES .eq. MASTER
     &      .AND.   KEEP(46)   .eq. 1
     &      .AND.   KEEP(52)   .eq. 0 ) THEN
            NULLIFY( id%DBLARR )
          ELSE
            IF (associated( id%DBLARR)) THEN
              DEALLOCATE(id%DBLARR)
              NULLIFY(id%DBLARR)
            ENDIF
          END IF
        END IF
      IF ( KEEP(19) .NE. 0 ) THEN
        IF ( KEEP(46) .NE. 1 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( KEEP(17), 1, MPI_INTEGER, 1, DEFIC_TAG,
     &                   id%COMM, STATUS, IERR )
          ELSE IF ( id%MYID .EQ. 1 ) THEN
            CALL MPI_SEND( KEEP(17), 1, MPI_INTEGER, 0, DEFIC_TAG,
     &                   id%COMM, IERR )
          END IF
        END IF
      END IF
      IF (associated(id%BUFR)) THEN
        DEALLOCATE(id%BUFR)
        NULLIFY(id%BUFR)
      END IF
      CALL SMUMPS_57( IERR )
      CALL SMUMPS_59( IERR )
      IF (KEEP(219).NE.0) THEN
      CALL SMUMPS_620()
      ENDIF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      CALL SMUMPS_770(id)
      IF (KEEP(201) .GT. 0) THEN
         IF ((KEEP(201).EQ.1) .OR. (KEEP(201).EQ.2)) THEN
            IF ( I_AM_SLAVE ) THEN
               CALL SMUMPS_591(IERR)
               IF(IERR.LT.0)THEN
                  INFO(1)=IERR
                  INFO(2)=0
               ENDIF
            ENDIF
            CALL MUMPS_276( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
         END IF
      END IF
      IF ( PROKG ) THEN
         CALL MUMPS_292(TIME)
         WRITE(MPG,180) TIME
      END IF
      PERLU_ON = .TRUE.
      CALL SMUMPS_214( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%LNA, id%NZ,
     &     id%NA_ELT,
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     PERLU_ON, TOTAL_BYTES)
      KEEP8(7) = TOTAL_BYTES
      id%INFO(22) = TOTAL_MBYTES
      IF ( MP .gt. 0 ) THEN
          WRITE(MP,'(A,I10) ')
     &    ' ** Effective minimum Space in MBYTES for facto  :',
     &                TOTAL_MBYTES
      ENDIF
      IF (I_AM_SLAVE) THEN
       K67 = KEEP8(67)
      ELSE
       K67 = 0_8
      ENDIF
      CALL MUMPS_735(K67,id%INFO(21))
      CALL SMUMPS_713(PROKG, MPG, K67, id%NSLAVES,
     & id%COMM, "effective space used in S   (KEEP8(67)   =")
      CALL MUMPS_243( id%MYID, id%COMM,
     &                    TOTAL_MBYTES, id%INFOG(21), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Rank of processor needing largest memory :',
     &  IRANK
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Space in MBYTES used by this processor   :',
     &  id%INFOG(21)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Avg. Space in MBYTES per working proc    :',
     &  ( id%INFOG(22)-TOTAL_MBYTES ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Avg. Space in MBYTES per working proc    :',
     &  id%INFOG(22) / id%NSLAVES
        END IF
      END IF
      KEEP(33) = INFO(11) 
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_REAL,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(247) = 0
      CALL MPI_REDUCE( KEEP(246), KEEP(247), 1, MPI_INTEGER, 
     &                 MPI_MAX, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_REAL,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MUMPS_646( KEEP8(31),KEEP8(6), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_735(KEEP8(6), INFOG(9))
      CALL MPI_REDUCE( INFO(10), INFOG(10), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_ALLREDUCE( INFO(11), INFOG(11), 1, MPI_INTEGER,
     &                 MPI_MAX, id%COMM, IERR)
      KEEP(133) = INFOG(11)
      CALL MPI_REDUCE( INFO(12), INFOG(12), 3, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( KEEP(103), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(229) = INFOG(25)
      CALL MPI_REDUCE( KEEP(105), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(230) = INFOG(25)
      INFO(25) = KEEP(98)
      CALL MPI_ALLREDUCE( INFO(25), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, id%COMM, IERR)
      CALL MUMPS_646( KEEP8(8), KEEP8(108), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_735(KEEP8(10), INFO(27))
      CALL MUMPS_646( KEEP8(10),KEEP8(110), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_735(KEEP8(110), INFOG(29))
      IF (KEEP(258).NE.0) THEN
        RINFOG(13)=0.0E0
        IF (KEEP(260).EQ.-1) THEN
          id%DKEEP(6)=-id%DKEEP(6)
        ENDIF
        CALL SMUMPS_764(
     &           id%COMM, id%DKEEP(6), KEEP(259),
     &           RINFOG(12), INFOG(34), id%NPROCS)
        IF (id%KEEP(50).EQ.0 .AND. id%MYID.EQ. MASTER) THEN
          IF (id%KEEP(23).NE.0) THEN
            CALL SMUMPS_767(
     &           RINFOG(12), id%N,
     &           id%STEP(1),
     &           id%UNS_PERM(1) )
          ENDIF
        ENDIF
      ENDIF
      IF(KEEP(110) .EQ. 1) THEN
         INFO(18) = KEEP(109)
         CALL MPI_ALLREDUCE( KEEP(109), KEEP(112), 1, MPI_INTEGER,
     &        MPI_SUM, id%COMM, IERR)
      ELSE
         INFO(18)  = 0
         KEEP(109) = 0
         KEEP(112) = 0
      ENDIF
      INFOG(28)=KEEP(112)+KEEP(17)
      IF (KEEP(17) .NE. 0) THEN
        IF (id%MYID .EQ. ID_ROOT) THEN
          INFO(18)=INFO(18)+KEEP(17)
        ENDIF
        IF (ID_ROOT .EQ. MASTER) THEN
          IF (id%MYID.EQ.MASTER) THEN
            DO I=1, KEEP(17)
              id%PIVNUL_LIST(KEEP(112)+I)=id%PIVNUL_LIST(KEEP(109)+I)
            ENDDO
          ENDIF
        ELSE
          IF (id%MYID .EQ. ID_ROOT) THEN
            CALL MPI_SEND(id%PIVNUL_LIST(KEEP(109)+1), KEEP(17),
     &                    MPI_INTEGER, MASTER, ZERO_PIV,
     &                    id%COMM, IERR)
          ELSE IF (id%MYID .EQ. MASTER) THEN
            CALL MPI_RECV(id%PIVNUL_LIST(KEEP(112)+1), KEEP(17),
     &                    MPI_INTEGER, ID_ROOT, ZERO_PIV,
     &                    id%COMM, STATUS, IERR )
          ENDIF
        ENDIF
      ENDIF
      IF(KEEP(110) .EQ. 1) THEN
         ALLOCATE(ITMP2(id%NPROCS),stat = IERR )  
         IF ( IERR .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%NPROCS
         END IF
         CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
         IF (INFO(1).LT.0) GOTO 490
         CALL MPI_GATHER ( KEEP(109),1, MPI_INTEGER, 
     &        ITMP2(1), 1, MPI_INTEGER, 
     &        MASTER, id%COMM, IERR)
         IF(id%MYID .EQ. MASTER) THEN
            POSBUF = ITMP2(1)+1
            KEEP(220)=1
            DO I = 1,id%NPROCS-1
               CALL MPI_RECV(id%PIVNUL_LIST(POSBUF), ITMP2(I+1), 
     &              MPI_INTEGER,I, 
     &              ZERO_PIV, id%COMM, STATUS, IERR)
               CALL MPI_SEND(POSBUF, 1, MPI_INTEGER, I, ZERO_PIV,
     &              id%COMM, IERR)
               POSBUF = POSBUF + ITMP2(I+1)
            ENDDO
         ELSE
            CALL MPI_SEND( id%PIVNUL_LIST(1), KEEP(109), MPI_INTEGER,
     &           MASTER,ZERO_PIV, id%COMM, IERR)
            CALL MPI_RECV( KEEP(220), 1, MPI_INTEGER, MASTER, ZERO_PIV,
     &           id%COMM, STATUS, IERR )
         ENDIF
      ENDIF
 490  IF (allocated(ITMP2)) DEALLOCATE(ITMP2)
      IF ( PROKG ) THEN
          WRITE(MPG,99984) RINFOG(2),RINFOG(3),KEEP8(6),INFOG(10),
     &                    INFOG(11), KEEP8(110)
          IF (id%KEEP(50) == 1 .OR. id%KEEP(50) == 2) THEN
            WRITE(MPG, 99987) INFOG(12)
          END IF
          IF (id%KEEP(50) == 0) THEN
            WRITE(MPG, 99985) INFOG(12)
          END IF
          IF (id%KEEP(50) .NE. 1) THEN
            WRITE(MPG, 99982) INFOG(13)
          END IF
          IF (KEEP(97) .NE. 0) THEN
            WRITE(MPG, 99986) KEEP(98)
          ENDIF
          IF (id%KEEP(50) == 2) THEN
             WRITE(MPG, 99988) KEEP(229)
             WRITE(MPG, 99989) KEEP(230)
          ENDIF
          IF (KEEP(110) .NE.0) THEN
              WRITE(MPG, 99991) KEEP(112)
          ENDIF
          IF ( KEEP(17) .ne. 0 )
     &    WRITE(MPG, 99983) KEEP(17)
          IF (KEEP(110).NE.0.OR.KEEP(17).NE.0)
     &    WRITE(MPG, 99992) KEEP(17)+KEEP(112)
          WRITE(MPG, 99981) INFOG(14)
          IF ((KEEP(201).EQ.0.OR.KEEP(201).EQ.2).AND.
     &        KEEP(50).EQ.0) THEN
          WRITE(MPG, 99980) KEEP8(108)
          ENDIF
          IF  ((KEEP(60).NE.0) .AND. INFOG(25).GT.0) THEN
           WRITE(MPG, '(A)') 
     & " ** Warning Static pivoting was necessary"
           WRITE(MPG, '(A)') 
     & " ** to factor interior variables with Schur ON"
          ENDIF
          IF (KEEP(258).NE.0) THEN
            WRITE(MPG,99978) RINFOG(12)
            WRITE(MPG,99977) INFOG(34)
          ENDIF
      END IF
 500  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
            CALL SMUMPS_592(id,IERR)
            IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
         ENDIF
         IF (WK_USER_PROVIDED) THEN
            NULLIFY(id%S)
         ELSE IF (KEEP(201).NE.0) THEN
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   
            KEEP8(23)=0_8
         ENDIF
      ELSE  
         IF (WK_USER_PROVIDED) THEN
            NULLIFY(id%S)
         ELSE
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   
            KEEP8(23)=0_8
         END IF
      END IF
 513  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         CALL SMUMPS_183( INFO(1), IERR )
         IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
      ENDIF
      CALL MUMPS_276( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
 530  CONTINUE
      IF (RHS_MUMPS_ALLOCATED) DEALLOCATE(RHS_MUMPS)
      NULLIFY(RHS_MUMPS)
      id%KEEP(13) = KEEP13_SAVE
      RETURN
 120  FORMAT(/' LOCAL REDISTRIB: DATA LOCAL/SENT         =',I12,I12)
 125  FORMAT(/' REDISTRIB: TOTAL DATA LOCAL/SENT         =',I12,I12)
 130  FORMAT(/' ****** FACTORIZATION STEP ********'/)
 160  FORMAT(' GLOBAL TIME FOR MATRIX DISTRIBUTION       =',F12.4)
 165  FORMAT(' Convergence error after scaling for INF-NORM',
     &       ' (option 7/8)   =',D9.2)
 166  FORMAT(' Convergence error after scaling for ONE-NORM',
     &       ' (option 7/8)   =',D9.2)
 170  FORMAT(/' STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Size of internal working array S         =',I12/
     &        ' Size of internal working array IS        =',I12/
     &        ' MINIMUM (ICNTL(14)=0) size of S          =',I12/
     &        ' MINIMUM (ICNTL(14)=0) size of IS         =',I12/
     &        ' REAL SPACE FOR ORIGINAL MATRIX           =',I12/
     &        ' INTEGER SPACE FOR ORIGINAL MATRIX        =',I12/
     &        ' REAL SPACE FOR FACTORS                   =',I12/
     &        ' INTEGER SPACE FOR FACTORS                =',I12/
     &        ' MAXIMUM FRONTAL SIZE (ESTIMATED)         =',I12)
 172  FORMAT(/' GLOBAL STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' NUMBER OF WORKING PROCESSES              =',I12/
     &        ' OUT-OF-CORE OPTION (ICNTL(22))           =',I12/
     &        ' REAL SPACE FOR FACTORS                   =',I12/
     &        ' INTEGER SPACE FOR FACTORS                =',I12/
     &        ' MAXIMUM FRONTAL SIZE (ESTIMATED)         =',I12/
     &        ' NUMBER OF NODES IN THE TREE              =',I12)
 173  FORMAT( ' PERFORM FORWARD DURING FACTO, NRHS       =',I12)
 175  FORMAT(/' NUMBER OF ENTRIES FOR // ROOT            =',I12)
 180  FORMAT(/' ELAPSED TIME FOR FACTORIZATION           =',F12.4)
99977 FORMAT( ' INFOG(34)  DETERMINANT (base 2 exponent) =',I12)
99978 FORMAT( ' RINFOG(12) DETERMINANT (real part)       =',F12.4)
99980 FORMAT( ' KEEP8(108) Extra copies IP stacking      =',I12)
99981 FORMAT( ' INFOG(14)  NUMBER OF MEMORY COMPRESS     =',I12)
99982 FORMAT( ' INFOG(13)  NUMBER OF DELAYED PIVOTS      =',I12)
99983 FORMAT( ' NB OF NULL PIVOTS DETECTED BY ICNTL(16)  =',I12)
99991 FORMAT( ' NB OF NULL PIVOTS DETECTED BY ICNTL(24)  =',I12)
99992 FORMAT( ' INFOG(28)  ESTIMATED DEFICIENCY          =',I12)
99984 FORMAT(/' GLOBAL STATISTICS '/
     &        ' RINFOG(2)  OPERATIONS IN NODE ASSEMBLY   =',1PD10.3/
     &        ' ------(3)  OPERATIONS IN NODE ELIMINATION=',1PD10.3/
     &        ' INFOG (9)  REAL SPACE FOR FACTORS        =',I12/
     &        ' INFOG(10)  INTEGER SPACE FOR FACTORS     =',I12/
     &        ' INFOG(11)  MAXIMUM FRONT SIZE            =',I12/
     &        ' INFOG(29)  NUMBER OF ENTRIES IN FACTORS  =',I12)
99985 FORMAT( ' INFOG(12) NB OF OFF DIAGONAL PIVOTS      =',I12)
99986 FORMAT( ' INFOG(25) NB TINY PIVOTS/STATIC PIVOTING =',I12)
99987 FORMAT( ' INFOG(12) NB OF NEGATIVE PIVOTS          =',I12)
99988 FORMAT( ' NUMBER OF 2x2 PIVOTS in type 1 nodes     =',I12)
99989 FORMAT( ' NUMBER OF 2x2 PIVOTS in type 2 nodes     =',I12)
      END SUBROUTINE SMUMPS_142
      SUBROUTINE SMUMPS_713(PROKG, MPG, VAL, NSLAVES,
     &     COMM, MSG)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL PROKG
      INTEGER MPG
      INTEGER(8) VAL
      INTEGER NSLAVES
      INTEGER COMM
      CHARACTER*42 MSG 
      INTEGER(8) MAX_VAL
      INTEGER IERR, MASTER
      REAL LOC_VAL, AVG_VAL
      PARAMETER(MASTER=0)
      CALL MUMPS_646( VAL, MAX_VAL, MPI_MAX, MASTER, COMM)
      LOC_VAL = real(VAL)/real(NSLAVES)
      CALL MPI_REDUCE( LOC_VAL, AVG_VAL, 1, MPI_REAL,
     &                 MPI_SUM, MASTER, COMM, IERR )
      IF (PROKG) THEN
        WRITE(MPG,100) " Maximum ", MSG, MAX_VAL
        WRITE(MPG,100) " Average ", MSG, int(AVG_VAL,8)
      ENDIF
      RETURN
 100  FORMAT(A9,A42,I12)
      END SUBROUTINE SMUMPS_713
      SUBROUTINE SMUMPS_770(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(SMUMPS_STRUC) :: id
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER :: ID_SCHUR, SIZE_SCHUR, LD_SCHUR, IB, BL4
      INTEGER :: ROW_LENGTH, I
      INTEGER(8) :: SURFSCHUR8, BL8, SHIFT8
      INTEGER(8) :: ISCHUR_SRC, ISCHUR_DEST, ISCHUR_SYM, ISCHUR_UNS
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      IF (id%INFO(1) .LT. 0) RETURN
      IF (id%KEEP(60) .EQ. 0) RETURN
      ID_SCHUR =MUMPS_275(
     &    id%PROCNODE_STEPS(id%STEP(max(id%KEEP(20),id%KEEP(38)))),
     &    id%NSLAVES)
      IF ( id%KEEP( 46 )  .NE. 1 ) THEN
        ID_SCHUR = ID_SCHUR + 1
      END IF
      IF (id%MYID.EQ.ID_SCHUR) THEN
        IF (id%KEEP(60).EQ.1) THEN
          LD_SCHUR =
     &    id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+2+id%KEEP(IXSZ))
          SIZE_SCHUR = LD_SCHUR - id%KEEP(253)
        ELSE
          LD_SCHUR   = -999999 
          SIZE_SCHUR = id%root%TOT_ROOT_SIZE
        ENDIF
      ELSE IF (id%MYID .EQ. MASTER) THEN
        SIZE_SCHUR = id%KEEP(116)
        LD_SCHUR = -44444 
      ELSE
        RETURN
      ENDIF
      SURFSCHUR8 = int(SIZE_SCHUR,8)*int(SIZE_SCHUR,8)
      IF (id%KEEP(60) .GT. 1) THEN
        IF (id%KEEP(221).EQ.1) THEN 
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR.EQ.MASTER) THEN 
              CALL scopy(SIZE_SCHUR,
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1), 1,
     &             id%REDRHS((I-1)*id%LREDRHS+1), 1)
            ELSE
              IF (id%MYID.EQ.ID_SCHUR) THEN
                CALL MPI_SEND(
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1),
     &             SIZE_SCHUR,
     &             MPI_REAL,
     &             MASTER, TAG_SCHUR,
     &             id%COMM, IERR )
              ELSE 
                CALL MPI_RECV( id%REDRHS((I-1)*id%LREDRHS+1),
     &             SIZE_SCHUR,
     &             MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &             id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
          ENDDO
          IF (id%MYID.EQ.ID_SCHUR) THEN
            DEALLOCATE(id%root%RHS_CNTR_MASTER_ROOT)
            NULLIFY   (id%root%RHS_CNTR_MASTER_ROOT)
          ENDIF
        ENDIF
        RETURN
      ENDIF
      IF (id%KEEP(252).EQ.0) THEN
        IF ( ID_SCHUR .EQ. MASTER ) THEN 
          CALL SMUMPS_756( SURFSCHUR8,
     &      id%S(id%PTRFAC(id%STEP(id%KEEP(20)))),
     &      id%SCHUR(1) )
        ELSE
          BL8=int(huge(BL4)/id%KEEP(35)/10,8)
          DO IB=1, int((SURFSCHUR8+BL8-1_8) / BL8)
            SHIFT8 = int(IB-1,8) * BL8                
            BL4    = int(min(BL8,SURFSCHUR8-SHIFT8)) 
            IF ( id%MYID .eq. ID_SCHUR ) THEN
              CALL MPI_SEND( id%S( SHIFT8 +
     &          id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ)))),
     &          BL4,
     &          MPI_REAL,
     &          MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
            ELSE IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_RECV( id%SCHUR(1_8 + SHIFT8),
     &                     BL4,
     &                     MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &                     id%COMM, STATUS, IERR )
            END IF
          ENDDO
        END IF
      ELSE
        ISCHUR_SRC = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &               +4+id%KEEP(IXSZ)))
        ISCHUR_DEST= 1_8
        DO I=1, SIZE_SCHUR
          ROW_LENGTH = SIZE_SCHUR
          IF (ID_SCHUR.EQ.MASTER) THEN 
            CALL scopy(ROW_LENGTH, id%S(ISCHUR_SRC), 1,
     &                 id%SCHUR(ISCHUR_DEST),1)
          ELSE
            IF (id%MYID.EQ.ID_SCHUR) THEN
              CALL MPI_SEND( id%S(ISCHUR_SRC), ROW_LENGTH,
     &        MPI_REAL,
     &        MASTER, TAG_SCHUR,
     &        id%COMM, IERR )
            ELSE
              CALL MPI_RECV( id%SCHUR(ISCHUR_DEST),
     &                   ROW_LENGTH,
     &                   MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &                   id%COMM, STATUS, IERR )
            ENDIF
          ENDIF
          ISCHUR_SRC = ISCHUR_SRC+int(LD_SCHUR,8)
          ISCHUR_DEST= ISCHUR_DEST+int(SIZE_SCHUR,8)
        ENDDO
        IF (id%KEEP(221).EQ.1) THEN 
          ISCHUR_SYM = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8) *
     &                    int(LD_SCHUR,8)
          ISCHUR_UNS =
     &                 id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8)
          ISCHUR_DEST = 1_8
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR .EQ. MASTER) THEN 
              IF (id%KEEP(50) .EQ. 0) THEN
                CALL scopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ELSE
                CALL scopy(SIZE_SCHUR, id%S(ISCHUR_SYM), 1,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ENDIF
            ELSE
              IF (id%MYID .NE. MASTER) THEN
                IF (id%KEEP(50) .EQ. 0) THEN
                  CALL scopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &            id%S(ISCHUR_SYM), 1)
                ENDIF
                CALL MPI_SEND(id%S(ISCHUR_SYM), SIZE_SCHUR,
     &          MPI_REAL, MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
              ELSE
                CALL MPI_RECV(id%REDRHS(ISCHUR_DEST),
     &          SIZE_SCHUR, MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &          id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
            IF (id%KEEP(50).EQ.0) THEN
              ISCHUR_UNS = ISCHUR_UNS + int(LD_SCHUR,8)
            ELSE
              ISCHUR_SYM = ISCHUR_SYM + int(LD_SCHUR,8)
            ENDIF
            ISCHUR_DEST = ISCHUR_DEST + int(id%LREDRHS,8)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_770
      SUBROUTINE SMUMPS_83
     & ( N, MAPPING, NZ, IRN, JCN, PROCNODE, STEP,
     &   SLAVEF, PERM, FILS,
     &   RG2L, KEEP,KEEP8, MBLOCK, NBLOCK, NPROW, NPCOL )
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, NZ, SLAVEF, MBLOCK, NBLOCK, NPROW, NPCOL
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IRN( NZ ), JCN( NZ ) 
      INTEGER MAPPING( NZ ), STEP( N )
      INTEGER PROCNODE( KEEP(28) ), PERM( N ), FILS( N ), RG2L( N )
      INTEGER MUMPS_275, MUMPS_330
      EXTERNAL MUMPS_275, MUMPS_330
      INTEGER K, IOLD, JOLD, INEW, JNEW, ISEND, JSEND, IARR, INODE
      INTEGER TYPE_NODE, DEST
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INODE = KEEP(38)
      K = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = K
        INODE = FILS( INODE )
        K = K + 1
      END DO
      DO K = 1, NZ
        IOLD = IRN( K )
        JOLD = JCN( K )
        IF ( IOLD .GT. N .OR. IOLD .LT. 1 .OR.
     &       JOLD .GT. N .OR. JOLD .LT. 1 ) THEN
           MAPPING( K ) = -1
           CYCLE
        END IF
        IF ( IOLD .eq. JOLD ) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM( IOLD )
          JNEW = PERM( JOLD )
          IF ( INEW .LT. JNEW ) THEN
            ISEND = IOLD
            IF ( KEEP(50) .ne. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          END IF
        END IF
        IARR = abs( ISEND )
        TYPE_NODE = MUMPS_330( PROCNODE(abs(STEP(IARR))),
     &                              SLAVEF )
        IF ( TYPE_NODE .eq. 1 .or. TYPE_NODE .eq. 2 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = MUMPS_275( PROCNODE(abs(STEP(IARR))),
     &                             SLAVEF ) + 1
          ELSE
            DEST = MUMPS_275( PROCNODE(abs(STEP(IARR))),
     &                             SLAVEF )
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L( JSEND )
            JPOSROOT = RG2L( IARR  )
          ELSE
            IPOSROOT = RG2L( IARR  )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/MBLOCK, NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/NBLOCK, NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * NPCOL + JCOL_GRID
          END IF
        END IF
        MAPPING( K ) = DEST
      END DO
      RETURN
      END SUBROUTINE SMUMPS_83
      SUBROUTINE SMUMPS_282(
     & N, NZ_loc, id,
     & DBLARR, LDBLARR, INTARR, LINTARR,
     & PTRAIW, PTRARW, KEEP,KEEP8, MYID, COMM, NBRECORDS,
     &
     & A, LA, root, PROCNODE_STEPS, SLAVEF, PERM, STEP,
     & ICNTL, INFO, NSEND, NLOCAL,
     & ISTEP_TO_INIV2, CANDIDATES
     & )
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, NZ_loc
      TYPE (SMUMPS_STRUC) :: id
      INTEGER LDBLARR, LINTARR
      REAL DBLARR( LDBLARR )
      INTEGER INTARR( LINTARR )
      INTEGER PTRAIW( N ), PTRARW( N )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER MYID, COMM, NBRECORDS
      INTEGER(8) :: LA
      INTEGER SLAVEF
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      INTEGER CANDIDATES(SLAVEF+1, max(1,KEEP(56)))
      REAL A( LA )
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER PROCNODE_STEPS(KEEP(28)), PERM( N ), STEP( N )
      INTEGER INFO( 40 ), ICNTL(40)
      INTEGER MUMPS_275, MUMPS_330, numroc, 
     &        MUMPS_810
      EXTERNAL MUMPS_275, MUMPS_330, numroc, 
     &        MUMPS_810
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE ), MSGSOU
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IW4
      INTEGER END_MSG_2_RECV
      INTEGER I, K, I1, IA
      INTEGER TYPE_NODE, DEST
      INTEGER IOLD, JOLD, IARR, ISEND, JSEND, INEW, JNEW
      INTEGER allocok,  TYPESPLIT, T4MASTER, INIV2
      LOGICAL T4_MASTER_CONCERNED
      REAL VAL
      INTEGER(8) :: PTR_ROOT
      INTEGER LOCAL_M, LOCAL_N, ARROW_ROOT
      INTEGER IROW_GRID, JCOL_GRID, IPOSROOT, JPOSROOT
      INTEGER MP,LP
      INTEGER KPROBE, FREQPROBE
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFI
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BUFR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI
      REAL, ALLOCATABLE, DIMENSION(:) :: BUFRECR
      INTEGER IACT( SLAVEF ), IREQI( SLAVEF ), IREQR( SLAVEF )
      LOGICAL SEND_ACTIVE( SLAVEF )
      LOGICAL FLAG
      INTEGER NSEND, NLOCAL
      INTEGER MASTER_NODE, ISTEP
      NSEND = 0
      NLOCAL = 0
      LP = ICNTL(1)
      MP = ICNTL(2)
      END_MSG_2_RECV = SLAVEF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, 2, SLAVEF ), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating int buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = ( NBRECORDS * 2 + 1 ) * SLAVEF * 2
      END IF
      ALLOCATE( BUFR( NBRECORDS, 2, SLAVEF), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating real buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFRECI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * 2 + 1
        GOTO 20
      END IF
      ALLOCATE( BUFRECR( NBRECORDS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS
        GOTO 20
      END IF
      ALLOCATE( IW4( N, 2 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(LP,*) '** Error allocating IW4 for matrix distribution'
        INFO(1) = -13
        INFO(2) = N * 2
      END IF
 20   CONTINUE
      CALL MUMPS_276( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      ARROW_ROOT = 0
      DO I = 1, N
          I1 = PTRAIW( I )
          IA = PTRARW( I )
          IF ( IA .GT. 0 ) THEN
            DBLARR( IA ) = ZERO
            IW4( I, 1 ) = INTARR( I1 )
            IW4( I, 2 ) = -INTARR( I1 + 1 )
            INTARR( I1 + 2 ) = I
          END IF
      END DO
      IF ( KEEP(38) .NE. 0 ) THEN
          IF (KEEP(60)==0) THEN
          LOCAL_M = numroc( root%ROOT_SIZE, root%MBLOCK,
     &               root%MYROW, 0, root%NPROW )
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
     &        (I-1)*root%SCHUR_LLD+root%SCHUR_MLOC) = ZERO
            ENDDO
          ENDIF
      END IF
      DO I = 1, SLAVEF
        BUFI( 1, 1, I ) = 0
      END DO
      DO I = 1, SLAVEF
        BUFI( 1, 2, I ) = 0
      END DO
      DO I = 1, SLAVEF
        SEND_ACTIVE( I ) = .FALSE.
        IACT( I ) = 1
      END DO
      KPROBE = 0
      FREQPROBE = max(1,NBRECORDS/10)
      DO K = 1, NZ_loc
        KPROBE = KPROBE + 1
        IF ( KPROBE .eq. FREQPROBE ) THEN
          KPROBE = 0
          CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     &                     FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            MSGSOU = STATUS( MPI_SOURCE )
            CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, 
     &                 MPI_INTEGER,
     &                 MSGSOU, ARR_INT, COMM, STATUS, IERR )
            CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_REAL,
     &                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
            CALL SMUMPS_102(
     &             BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &             KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &             A, LA,
     &             END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &             ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR
     &             )
          END IF
        END IF
        IOLD = id%IRN_loc(K)
        JOLD = id%JCN_loc(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     &                 .OR.(JOLD.LT.1) ) CYCLE
        VAL = id%A_loc(K)
        IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) THEN
          VAL = VAL * id%ROWSCA(IOLD)*id%COLSCA(JOLD)
        ENDIF
        IF (IOLD.EQ.JOLD) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM(IOLD)
          JNEW = PERM(JOLD)
          IF (INEW.LT.JNEW) THEN
            ISEND = IOLD
            IF ( KEEP(50) .NE. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          ENDIF
        ENDIF
        IARR = abs( ISEND )
        ISTEP = abs(STEP(IARR))
        TYPE_NODE = MUMPS_330(   PROCNODE_STEPS(ISTEP),
     &                                SLAVEF )
        MASTER_NODE= MUMPS_275(  PROCNODE_STEPS(ISTEP),
     &                                SLAVEF )
        TYPESPLIT  = MUMPS_810( PROCNODE_STEPS(ISTEP),
     &                               SLAVEF )
        T4_MASTER_CONCERNED = .FALSE.
        T4MASTER               = -9999
        IF (TYPE_NODE.EQ.2) THEN
         INIV2         = ISTEP_TO_INIV2(ISTEP)
         IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
          T4_MASTER_CONCERNED = .TRUE.
          T4MASTER=CANDIDATES(CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
         ENDIF
        ENDIF
        IF ( TYPE_NODE .eq. 1 ) THEN
          DEST = MASTER_NODE
        ELSE IF ( TYPE_NODE .eq. 2 ) THEN
          IF ( ISEND .LT. 0 ) THEN
            DEST = -1
          ELSE
            DEST = MASTER_NODE
          END IF
        ELSE
          IF ( ISEND < 0 ) THEN
            IPOSROOT = root%RG2L_ROW(JSEND)
            JPOSROOT = root%RG2L_ROW(IARR )
          ELSE
            IPOSROOT = root%RG2L_ROW(IARR )
            JPOSROOT = root%RG2L_ROW(JSEND)
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          DEST = IROW_GRID * root%NPCOL + JCOL_GRID
        END IF
        if (DEST .eq. -1) then
          NLOCAL = NLOCAL + 1
          NSEND = NSEND + SLAVEF -1
        else
          if (DEST .eq.MYID ) then
            NLOCAL = NLOCAL + 1
          else
            NSEND = NSEND + 1
          endif
        end if
        IF ( DEST.EQ.-1) THEN
         DO I=1, CANDIDATES(SLAVEF+1,ISTEP_TO_INIV2(ISTEP))
            DEST=CANDIDATES(I,ISTEP_TO_INIV2(ISTEP))
            CALL SMUMPS_101( DEST, ISEND, JSEND, VAL,
     &   BUFI, BUFR, BUFRECI, BUFRECR,
     &   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1), 
     &   root, KEEP,KEEP8 )
         ENDDO
         DEST=MASTER_NODE
         CALL SMUMPS_101( DEST, ISEND, JSEND, VAL,
     &   BUFI, BUFR, BUFRECI, BUFRECR,
     &   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1),
     &   root, KEEP,KEEP8 )
         IF (T4_MASTER_CONCERNED) THEN
          DEST = T4MASTER
          CALL SMUMPS_101( DEST, ISEND, JSEND, VAL,
     &    BUFI, BUFR, BUFRECI, BUFRECR,
     &    NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &    SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &    N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &    PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1),
     &    root, KEEP,KEEP8 )
         ENDIF
        ELSE
         CALL SMUMPS_101( DEST, ISEND, JSEND, VAL,
     &   BUFI, BUFR, BUFRECI, BUFRECR,
     &   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1), 
     &   root, KEEP,KEEP8 )
         IF (T4_MASTER_CONCERNED) THEN
          DEST = T4MASTER
          CALL SMUMPS_101( DEST, ISEND, JSEND, VAL,
     &    BUFI, BUFR, BUFRECI, BUFRECR,
     &    NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &    SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &    N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &    PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1),
     &    root, KEEP,KEEP8 )
         ENDIF
        ENDIF
      END DO
      DEST = -2
        CALL SMUMPS_101( DEST, ISEND, JSEND, VAL,
     &  BUFI, BUFR, BUFRECI, BUFRECR,
     &  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &  N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, 
     &  IW4(1,1), root, KEEP,KEEP8 )
      DO WHILE ( END_MSG_2_RECV .NE. 0 )
        CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, MPI_INTEGER,
     &                 MPI_ANY_SOURCE, ARR_INT, COMM, STATUS, IERR )
        MSGSOU = STATUS( MPI_SOURCE )
        CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_REAL,
     &                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
        CALL SMUMPS_102(
     &           BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &           KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &           A, LA,
     &           END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &           ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &           INTARR, LINTARR, DBLARR, LDBLARR
     &           )
      END DO
      DO I = 1, SLAVEF
        IF ( SEND_ACTIVE( I ) ) THEN
          CALL MPI_WAIT( IREQI( I ), STATUS, IERR )
          CALL MPI_WAIT( IREQR( I ), STATUS, IERR )
        END IF
      END DO
      KEEP(49) = ARROW_ROOT
      DEALLOCATE( IW4 )
      DEALLOCATE( BUFI )
      DEALLOCATE( BUFR )
      DEALLOCATE( BUFRECI )
      DEALLOCATE( BUFRECR )
      RETURN
      END SUBROUTINE SMUMPS_282
      SUBROUTINE SMUMPS_101( DEST, ISEND, JSEND, VAL,
     &  BUFI, BUFR, BUFRECI, BUFRECR,
     &  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR, N,
     &  PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4, root,
     &  KEEP,KEEP8 )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER ISEND, JSEND, DEST, NBRECORDS, SLAVEF, COMM, MYID, N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER ARROW_ROOT, END_MSG_2_RECV, LOCAL_M, LOCAL_N
      INTEGER LINTARR, LDBLARR
      INTEGER(8) :: LA, PTR_ROOT
      INTEGER BUFI( NBRECORDS * 2 + 1, 2, SLAVEF )
      INTEGER BUFRECI( NBRECORDS * 2 + 1 )
      INTEGER IREQI(SLAVEF), IREQR(SLAVEF), IACT(SLAVEF)
      INTEGER IW4( N, 2 )
      INTEGER PTRAIW( N ), PTRARW( N ), PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR( LINTARR )
      REAL DBLARR( LDBLARR ), A( LA )
      LOGICAL SEND_ACTIVE(SLAVEF)
      REAL BUFR( NBRECORDS, 2, SLAVEF )
      REAL BUFRECR( NBRECORDS )
      REAL VAL
      INTEGER ISLAVE, IBEG, IEND, NBREC, IREQ
      INTEGER TAILLE_SEND_I, TAILLE_SEND_R, MSGSOU
      LOGICAL FLAG, SEND_LOCAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, STATUS(MPI_STATUS_SIZE)
      IF ( DEST .eq. -2 ) THEN
        IBEG = 1
        IEND = SLAVEF
      ELSE
        IBEG = DEST + 1
        IEND = DEST + 1
      END IF
      SEND_LOCAL = .FALSE.
      DO ISLAVE = IBEG, IEND
        NBREC = BUFI(1,IACT(ISLAVE),ISLAVE)
        IF ( DEST .eq. -2 ) THEN
          BUFI(1,IACT(ISLAVE),ISLAVE) = - NBREC
        END IF
        IF ( DEST .eq. -2 .or. NBREC + 1 > NBRECORDS ) THEN
          DO WHILE ( SEND_ACTIVE( ISLAVE ) )
            CALL MPI_TEST( IREQR( ISLAVE ), FLAG, STATUS, IERR )
            IF ( .NOT. FLAG ) THEN
                CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     &                           FLAG, STATUS, IERR )
                IF ( FLAG ) THEN
                  MSGSOU = STATUS(MPI_SOURCE)
                  CALL MPI_RECV( BUFRECI(1), 2*NBRECORDS+1,
     &                  MPI_INTEGER, MSGSOU, ARR_INT, COMM,
     &                  STATUS, IERR )
                  CALL MPI_RECV( BUFRECR(1), NBRECORDS,
     &                  MPI_REAL, MSGSOU,
     &                  ARR_REAL, COMM, STATUS, IERR )
                  CALL SMUMPS_102(
     &              BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &              A, LA,
     &              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &              INTARR, LINTARR, DBLARR, LDBLARR
     &              )
                END IF
            ELSE
                CALL MPI_WAIT( IREQI( ISLAVE ), STATUS, IERR )
                SEND_ACTIVE( ISLAVE ) = .FALSE.
            END IF
          END DO
          IF ( ISLAVE - 1 .ne. MYID ) THEN
            TAILLE_SEND_I = NBREC * 2 + 1
            TAILLE_SEND_R = NBREC
            CALL MPI_ISEND( BUFI(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_I,
     &          MPI_INTEGER, ISLAVE - 1, ARR_INT, COMM,
     &          IREQI( ISLAVE ), IERR )
            CALL MPI_ISEND( BUFR(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_R,
     &          MPI_REAL, ISLAVE - 1, ARR_REAL, COMM,
     &          IREQR( ISLAVE ), IERR )
            SEND_ACTIVE( ISLAVE ) = .TRUE.
          ELSE
            SEND_LOCAL = .TRUE.
          END IF
          IACT( ISLAVE ) = 3 - IACT( ISLAVE )
          BUFI( 1, IACT( ISLAVE ), ISLAVE ) = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          IREQ = BUFI(1,IACT(ISLAVE),ISLAVE) + 1
          BUFI(1,IACT(ISLAVE),ISLAVE) = IREQ
          BUFI(IREQ*2,IACT(ISLAVE),ISLAVE)  = ISEND
          BUFI(IREQ*2+1,IACT(ISLAVE),ISLAVE) = JSEND
          BUFR(IREQ,IACT(ISLAVE),ISLAVE )    = VAL
        END IF
      END DO
      IF ( SEND_LOCAL ) THEN
            ISLAVE = MYID + 1
            CALL SMUMPS_102(
     &              BUFI(1,3-IACT(ISLAVE),ISLAVE),
     &              BUFR(1,3-IACT(ISLAVE),ISLAVE),
     &              NBRECORDS, N, IW4(1,1),
     &              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &              A, LA,
     &              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &              INTARR, LINTARR, DBLARR, LDBLARR
     &              )
      END IF
      RETURN
      END SUBROUTINE SMUMPS_101
      SUBROUTINE SMUMPS_102
     &           ( BUFI, BUFR, NBRECORDS, N, IW4,
     &             KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT, A, LA,
     &             END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     &             SLAVEF, ARROW_ROOT,
     &             PTRAIW, PTRARW, PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER NBRECORDS, N, ARROW_ROOT, MYID, SLAVEF
      INTEGER BUFI( NBRECORDS * 2 + 1 )
      REAL BUFR( NBRECORDS )
      INTEGER IW4( N, 2 )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER END_MSG_2_RECV
      INTEGER PTRAIW( N ), PTRARW( N ), PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER LINTARR, LDBLARR
      INTEGER INTARR( LINTARR )
      INTEGER LOCAL_M, LOCAL_N
      INTEGER(8) :: PTR_ROOT, LA
      REAL A( LA ), DBLARR( LDBLARR )
      INTEGER MUMPS_330, MUMPS_275
      EXTERNAL MUMPS_330, MUMPS_275
      INTEGER IREC, NB_REC, NODE_TYPE, IPROC
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID,
     &        ILOCROOT, JLOCROOT
      INTEGER IA, IS1, ISHIFT, IIW, IS, IAS, IARR, JARR
      INTEGER TAILLE
      REAL VAL
      NB_REC = BUFI( 1 )
      IF ( NB_REC .LE. 0 ) THEN
        END_MSG_2_RECV = END_MSG_2_RECV - 1
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) GOTO 100
      DO IREC = 1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        NODE_TYPE = MUMPS_330( 
     &              PROCNODE_STEPS(abs(STEP(abs( IARR )))),
     &              SLAVEF )
        IF ( NODE_TYPE .eq. 3 ) THEN
          ARROW_ROOT = ARROW_ROOT + 1
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          IF ( IROW_GRID .NE. root%MYROW .OR.
     &       JCOL_GRID .NE. root%MYCOL ) THEN
            WRITE(*,*) MYID,':INTERNAL Error: recvd root arrowhead '
            WRITE(*,*) MYID,':not belonging to me. IARR,JARR=',IARR,JARR
            WRITE(*,*) MYID,':IROW_GRID,JCOL_GRID=',IROW_GRID,JCOL_GRID
            WRITE(*,*) MYID,':MYROW, MYCOL=', root%MYROW, root%MYCOL
            WRITE(*,*) MYID,':IPOSROOT,JPOSROOT=', IPOSROOT, JPOSROOT
            CALL MUMPS_ABORT()
          END IF
          ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
          JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
          IF (KEEP(60)==0) THEN
            A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT-1,8)) =  A( PTR_ROOT
     &        + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT - 1,8) )
     &      + VAL
          ELSE
            root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8) )
     &      = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8))
     &      + VAL
          ENDIF
        ELSE IF (IARR.GE.0) THEN
         IF (IARR.EQ.JARR) THEN
          IA = PTRARW(IARR)
          DBLARR(IA) = DBLARR(IA) + VAL
         ELSE
          IS1 =  PTRAIW(IARR)
          ISHIFT      = INTARR(IS1) + IW4(IARR,2)
          IW4(IARR,2) = IW4(IARR,2) - 1
          IIW         = IS1 + ISHIFT + 2
          INTARR(IIW)     = JARR
          IS          = PTRARW(IARR)
          IAS         = IS + ISHIFT
          DBLARR(IAS) = VAL
         ENDIF
        ELSE
           IARR = -IARR
           ISHIFT      = PTRAIW(IARR)+IW4(IARR,1)+2
           INTARR(ISHIFT)  = JARR
           IAS         = PTRARW(IARR)+IW4(IARR,1)
           IW4(IARR,1) = IW4(IARR,1) - 1
           DBLARR(IAS)      = VAL
           IPROC = MUMPS_275( PROCNODE_STEPS(abs(STEP(IARR))),
     &                             SLAVEF )
           IF ( (KEEP(50) .NE. 0 .OR. KEEP(234).NE.0)
     &          .AND.
     &          IW4(IARR,1) .EQ. 0 .AND.
     &          IPROC .EQ. MYID
     &          .AND. STEP(IARR) > 0 ) THEN
             TAILLE = INTARR( PTRAIW(IARR) )
             CALL SMUMPS_310( N, PERM,
     &            INTARR( PTRAIW(IARR) + 3 ),
     &            DBLARR( PTRARW(IARR) + 1 ),
     &            TAILLE, 1, TAILLE )
           END IF
        ENDIF
      ENDDO
 100  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_102
      SUBROUTINE SMUMPS_151( NRHS, N, KEEP28, IWCB, LIWW,
     &       W, LWC,
     &       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER NRHS, N,LIWW,LWC,POSWCB,IWPOSCB, KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28),PTRACB(KEEP28)
      REAL W(LWC)
      INTEGER SIZFI, SIZFR
      IF ( IWPOSCB .eq. LIWW ) RETURN
      DO WHILE ( IWCB( IWPOSCB + 2 ) .eq. 0 )
        SIZFR = IWCB( IWPOSCB + 1 )
        SIZFI =  2  
        SIZFR = SIZFR * NRHS
        IWPOSCB = IWPOSCB + SIZFI
        POSWCB  = POSWCB  + SIZFR
        IF ( IWPOSCB .eq. LIWW ) RETURN
      END DO
      RETURN
      END SUBROUTINE SMUMPS_151
      SUBROUTINE SMUMPS_95(NRHS,N,KEEP28,IWCB,LIWW,W,LWC,
     &       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER NRHS,N,LIWW,LWC,POSWCB,IWPOSCB,KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28),PTRACB(KEEP28)
      REAL W(LWC)
      INTEGER IPTIW,IPTA,SIZFI,SIZFR,LONGI,LONGR
      INTEGER I
      IPTIW = IWPOSCB
      IPTA  = POSWCB
      LONGI = 0
      LONGR = 0
      IF ( IPTIW .EQ. LIWW ) RETURN
10    CONTINUE
      IF (IWCB(IPTIW+2).EQ.0) THEN
        SIZFR  = IWCB(IPTIW+1)
        SIZFI =  2  
        SIZFR  = SIZFR * NRHS
        IF (LONGI.NE.0) THEN
          DO 20 I=0,LONGI-1
            IWCB(IPTIW + SIZFI - I) = IWCB (IPTIW - I )
 20       CONTINUE 
          DO 30 I=0,LONGR-1
            W(IPTA + SIZFR - I)   = W(IPTA - I )
 30       CONTINUE
        ENDIF
        DO 40 I=1,KEEP28
          IF ((PTRICB(I).LE.(IPTIW+1)).AND.
     &        (PTRICB(I).GT.IWPOSCB) ) THEN
            PTRICB(I) = PTRICB(I) + SIZFI
            PTRACB(I) = PTRACB(I) + SIZFR
          ENDIF 
40      CONTINUE 
        IWPOSCB = IWPOSCB + SIZFI
        IPTIW   = IPTIW + SIZFI
        POSWCB = POSWCB + SIZFR
        IPTA   = IPTA + SIZFR     
       ELSE
        SIZFR  = IWCB(IPTIW+1)
        SIZFI  = 2
        SIZFR  = SIZFR * NRHS
        IPTIW = IPTIW + SIZFI
        LONGI = LONGI + SIZFI
        IPTA  = IPTA + SIZFR
        LONGR = LONGR + SIZFR
       ENDIF
       IF (IPTIW.NE.LIWW) GOTO 10
       RETURN
       END SUBROUTINE SMUMPS_95
      SUBROUTINE SMUMPS_205(MTYPE, IFLAG, N, NZ,
     &    LHS, WRHS, W, RHS, GIVSOL, SOL, ANORM, XNORM, SCLNRM,
     &    MPRINT, ICNTL, KEEP,KEEP8)
      INTEGER MTYPE,N,NZ,IFLAG,ICNTL(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      REAL RHS(N),LHS(N)
      REAL WRHS(N),SOL(*)
      REAL W(N)
      REAL RESMAX,RESL2,XNORM, ERMAX,MAXSOL,
     &     COMAX, SCLNRM, ERL2, ERREL
      REAL ANORM,DZERO,EPSI
      LOGICAL GIVSOL,PROK
      INTEGER MPRINT, MP
      INTEGER K
      INTRINSIC abs, max, sqrt
      MP = ICNTL(2)
      PROK = (MPRINT .GT. 0)
      DZERO = 0.0E0
      EPSI = 0.1E-9
      ANORM = DZERO
      RESMAX = DZERO
      RESL2 = DZERO
      DO 40 K = 1, N
        RESMAX = max(RESMAX, abs(RHS(K)))
        RESL2 = RESL2 + abs(RHS(K)) * abs(RHS(K))
        ANORM = max(ANORM, W(K))
   40 CONTINUE
      XNORM = DZERO
      DO 50 K = 1, N
        XNORM = max(XNORM, abs(LHS(K)))
   50 CONTINUE
      IF (XNORM .GT. EPSI) THEN
        SCLNRM = RESMAX / (ANORM * XNORM)
      ELSE
        IFLAG = IFLAG + 2
        IF ((MP .GT. 0) .AND. (ICNTL(4) .GE. 2)) WRITE( MP, * ) 
     &' max-NORM of computed solut. is zero'
        SCLNRM = RESMAX / ANORM
      ENDIF
      RESL2 = sqrt(RESL2)
      ERMAX = DZERO
      COMAX = DZERO
      ERL2 = DZERO
      IF (.NOT.GIVSOL) THEN
        IF (PROK) WRITE( MPRINT, 90 ) RESMAX, RESL2, ANORM, XNORM, 
     &      SCLNRM
      ELSE
        MAXSOL = DZERO
        DO 60 K = 1, N
          MAXSOL = max(MAXSOL, abs(SOL(K)))
   60   CONTINUE
        DO 70 K = 1, N
          ERL2 = abs(LHS(K) - SOL(K)) ** 2 + ERL2
          ERMAX = max(ERMAX, abs(LHS(K) - SOL(K)))
   70   CONTINUE
        DO 80 K = 1, N
          IF (abs(SOL(K)) .GT. EPSI) THEN
            COMAX = max(COMAX, (abs(LHS(K) - SOL(K)) / abs(SOL(K))))
          ENDIF
   80   CONTINUE
        ERL2 = sqrt(ERL2)
        IF (MAXSOL .GT. EPSI) THEN
          ERREL = ERMAX / MAXSOL
        ELSE
          IFLAG = IFLAG + 2
          IF ((MP .GT. 0) .AND. (ICNTL(4) .GE. 2)) WRITE( MP, * ) 
     &' MAX-NORM of exact solution is zero'
          ERREL = ERMAX
        ENDIF
        IF (PROK) WRITE( MPRINT, 100 ) ERMAX, ERL2, ERREL, COMAX, RESMAX
     &      , RESL2, ANORM, XNORM, SCLNRM
      ENDIF
   90  FORMAT (/' RESIDUAL IS ............ (MAX-NORM)        =',1PD9.2/
     &       '                       .. (2-NORM)          =',1PD9.2/
     &       ' RINFOG(4):NORM OF input  Matrix  (MAX-NORM)=',1PD9.2/
     &       ' RINFOG(5):NORM OF Computed SOLUT (MAX-NORM)=',1PD9.2/
     &       ' RINFOG(6):SCALED RESIDUAL ...... (MAX-NORM)=',1PD9.2)
      RETURN
  100  FORMAT (/' ERROR IS     ............ (MAX-NORM)       =',1PD9.2/
     &       '              ............ (2-NORM)         =',1PD9.2/
     &       ' RELATIVE ERROR........... (MAX-NORM)       =',1PD9.2/
     &       ' Comp. Wise ERROR......... (MAX-NORM)       =',1PD9.2/
     &       ' AND RESIDUAL IS ......... (MAX-NORM)       =',1PD9.2/
     &       '                        .. (2-NORM)         =',1PD9.2/
     &       ' NORM OF input  MATRIX ... (MAX-NORM)       =',1PD9.2/
     &       ' NORM of computed SOLUT... (MAX-NORM)       =',1PD9.2/
     &       ' SCALED RESIDUAL ......... (MAX-NORM)       =',1PD9.2)
      END SUBROUTINE SMUMPS_205
      SUBROUTINE SMUMPS_206(NZ, N, RHS,
     &    X, Y, D, R_W, C_W, IW, KASE,
     &    OMEGA, ERX, JOB, COND, MAXIT, NOITER, LP, KEEP,KEEP8,
     &    ARRET )
      IMPLICIT NONE
      INTEGER NZ, N, KASE, KEEP(500), JOB
      INTEGER(8) KEEP8(150)
      INTEGER IW(N,2)
      REAL RHS(N)
      REAL X(N), Y(N)
      REAL D(N)
      REAL R_W(N,2)
      REAL C_W(N)
      INTEGER LP, MAXIT, NOITER
      REAL COND(2),OMEGA(2)
      REAL ARRET
      REAL CGCE, CTAU
      DATA  CTAU /1.0E3/, CGCE /0.2E0/
      LOGICAL LCOND1, LCOND2
      INTEGER IFLAG, JUMP, I, IMAX
      REAL ERX, DXMAX
      REAL CONVER, OM1, OM2, DXIMAX
      REAL ZERO, ONE,TAU, DD
      REAL OLDOMG(2)
      INTEGER SMUMPS_IXAMAX
      INTRINSIC     abs, max
      SAVE LCOND1, LCOND2, JUMP,  DXIMAX, DXMAX, CONVER,
     &     OM1, OLDOMG, IFLAG
      DATA ZERO /0.0E0/, ONE /1.0E0/
      IF (KASE .EQ. 0) THEN
        LCOND1 = .FALSE.
        LCOND2 = .FALSE.
        COND(1) = ONE
        COND(2) = ONE
        ERX = ZERO
        OM1 = ZERO
        IFLAG = 0
        NOITER = 0
        JUMP = 1
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 30
      CASE(2)
        GOTO 10
      CASE(3)
        GOTO 110
      CASE(4)
        GOTO 150
      CASE(5)
        GOTO 35
      CASE DEFAULT
      END SELECT
   10 CONTINUE
      DO 20 I = 1, N
        X(I) = X(I) + Y(I)
   20 CONTINUE
      IF (NOITER .GT. MAXIT) THEN
        IFLAG = IFLAG + 8
        GOTO 70
      ENDIF
   30 CONTINUE
      KASE = 14
      JUMP = 5
      RETURN
   35 CONTINUE
      IMAX = SMUMPS_IXAMAX(N, X, 1)
      DXMAX = abs(X(IMAX))
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO 40 I = 1, N
        TAU = (R_W(I, 2) * DXMAX + abs(RHS(I))) * real(N) * CTAU
        DD = R_W(I, 1) + abs(RHS(I))
        IF ((DD + TAU) .GT. TAU) THEN
          OMEGA(1) = max(OMEGA(1), abs(Y(I)) / DD)
          IW(I, 1) = 1
        ELSE
          IF (TAU .GT. ZERO) THEN
            OMEGA(2) = max(OMEGA(2),
     &                     abs(Y(I)) / (DD + R_W(I, 2) * DXMAX))
          ENDIF
          IW(I, 1) = 2
        ENDIF
   40 CONTINUE
      OM2 = OMEGA(1) + OMEGA(2)
      IF (OM2 .LT. ARRET ) GOTO 70
      IF (MAXIT .EQ. 0) GOTO 70
      IF (NOITER .GT. 1 .AND. OM2 .GT. OM1 * CGCE) THEN
        CONVER = OM2 / OM1
        IF (OM2 .GT. OM1) THEN
          OMEGA(1) = OLDOMG(1)
          OMEGA(2) = OLDOMG(2)
          DO 50 I = 1, N
            X(I) = C_W(I)
   50     CONTINUE
        ENDIF
        GOTO 70
      ENDIF
      DO 60 I = 1, N
        C_W(I) = X(I)
   60 CONTINUE
      OLDOMG(1) = OMEGA(1)
      OLDOMG(2) = OMEGA(2)
      OM1 = OM2
      NOITER = NOITER + 1
      KASE = 2
      JUMP = 2
      RETURN
   70 KASE = 0
      IF (JOB .LE. 0) GOTO 170
      DO 80 I = 1, N
        IF (IW(I, 1) .EQ. 1) THEN
          R_W(I, 1) = R_W(I, 1) + abs(RHS(I))
          R_W(I, 2) = ZERO
          LCOND1 = .TRUE.
        ELSE
          R_W(I, 2) = R_W(I, 2) * DXMAX + R_W(I, 1)
          R_W(I, 1) = ZERO
          LCOND2 = .TRUE.
        ENDIF
   80 CONTINUE
      DO 90 I = 1, N
        C_W(I) = X(I) * D(I)
   90 CONTINUE
      IMAX = SMUMPS_IXAMAX(N, C_W(1), 1)
      DXIMAX = abs(C_W(IMAX))
      IF (.NOT.LCOND1) GOTO 130
  100 CALL SMUMPS_218(N, KASE, Y, COND(1), C_W, IW(1, 2))
      IF (KASE .EQ. 0) GOTO 120
      IF (KASE .EQ. 1) CALL SMUMPS_204(N, Y, D)
      IF (KASE .EQ. 2) CALL SMUMPS_204(N, Y, R_W)
      JUMP = 3
      RETURN
  110 CONTINUE
      IF (KASE .EQ. 1) CALL SMUMPS_204(N, Y, R_W)
      IF (KASE .EQ. 2) CALL SMUMPS_204(N, Y, D)
      GOTO 100
  120 IF (DXIMAX .GT. ZERO) COND(1) = COND(1) / DXIMAX
      ERX = OMEGA(1) * COND(1)
  130 IF (.NOT.LCOND2) GOTO 170
      KASE = 0
  140 CALL SMUMPS_218(N, KASE, Y, COND(2), C_W, IW(1, 2))
      IF (KASE .EQ. 0) GOTO 160
      IF (KASE .EQ. 1) CALL SMUMPS_204(N, Y, D)
      IF (KASE .EQ. 2) CALL SMUMPS_204(N, Y, R_W(1, 2))
      JUMP = 4
      RETURN
  150 CONTINUE
      IF (KASE .EQ. 1) CALL SMUMPS_204(N, Y, R_W(1, 2))
      IF (KASE .EQ. 2) CALL SMUMPS_204(N, Y, D)
      GOTO 140
  160 IF (DXIMAX .GT. ZERO) THEN
        COND(2) = COND(2) / DXIMAX
      ENDIF
      ERX = ERX + OMEGA(2) * COND(2)
  170 KASE = -IFLAG
      RETURN
      END SUBROUTINE SMUMPS_206
      SUBROUTINE SMUMPS_207(A, NZ, N, IRN, ICN, Z, KEEP,KEEP8)
      INTEGER NZ, N, I, J, K, KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IRN(NZ), ICN(NZ)
      REAL A(NZ)
      REAL Z(N)
      REAL ZERO
      INTRINSIC     abs
      DATA ZERO /0.0E0/
      DO 10 I = 1, N
        Z(I) = ZERO
   10 CONTINUE
      IF (KEEP(50) .EQ.0) THEN
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + abs(A(K))
       ENDDO
      ELSE
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + abs(A(K))
        IF (J.NE.I) THEN 
          Z(J) = Z(J) + abs(A(K))
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_207
      SUBROUTINE SMUMPS_289(A, NZ, N, IRN, ICN, Z,
     &            KEEP, KEEP8, COLSCA)
      INTEGER,   intent(in)  :: NZ, N, KEEP(500)
      INTEGER(8), intent(in)  :: KEEP8(150)
      INTEGER,   intent(in)  :: IRN(NZ), ICN(NZ)
      REAL,   intent(in)  :: A(NZ)
      REAL,      intent(in)  :: COLSCA(N)
      REAL,      intent(out) :: Z(N)
      REAL  ZERO
      DATA ZERO /0.0E0/
      INTEGER I, J, K
      DO 10 I = 1, N
        Z(I) = ZERO
   10 CONTINUE
      IF (KEEP(50) .EQ.0) THEN
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + abs(A(K)*COLSCA(J))
       ENDDO
      ELSE
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + abs(A(K)*COLSCA(J))
        IF (J.NE.I) THEN
          Z(J) = Z(J) + abs(A(K)*COLSCA(I))
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_289
      SUBROUTINE SMUMPS_208(A, NZ, N, IRN, ICN, RHS, X, R, W,
     &           KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER,   intent(in)   :: NZ, N, KEEP(500)
      INTEGER(8), intent(in)   ::  KEEP8(150)
      INTEGER,   intent(in)   :: IRN(NZ), ICN(NZ)
      REAL,   intent(in)   :: A(NZ), RHS(N), X(N)
      REAL,      intent(out)  :: W(N)
      REAL,   intent(out)  :: R(N)
      INTEGER I, K, J
      REAL ZERO
      DATA ZERO /0.0E0/
      REAL D
      DO I = 1, N
        R(I) = RHS(I)
        W(I) = ZERO
      ENDDO
      DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .GT. N) .OR. (J .GT. N) .OR. (I .LT. 1) .OR. (J .LT. 1))
     &      CYCLE
        D = A(K) * X(J)
        R(I) = R(I) - D
        W(I) = W(I) + abs(D)
        IF ((I.NE.J) .AND. (KEEP(50).NE.0) ) THEN
          D = A(K) * X(I)
          R(J) = R(J) - D
          W(J) = W(J) + abs(D)
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_208
      SUBROUTINE SMUMPS_204(N, R, W)
      INTEGER, intent(in)  :: N
      REAL,    intent(in)  :: W(N)
      REAL, intent(inout) :: R(N)
      INTEGER I
      DO 10 I = 1, N
        R(I) = R(I) * W(I)
   10 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_204
      SUBROUTINE SMUMPS_218(N, KASE, X, EST, W, IW)
      INTEGER, intent(in)    :: N
      INTEGER, intent(inout) :: KASE
      INTEGER IW(N)
      REAL W(N), X(N)
      REAL EST
      INTRINSIC abs, nint, real, sign
      INTEGER SMUMPS_IXAMAX
      EXTERNAL SMUMPS_IXAMAX
      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      INTEGER I, ITER, J, JLAST, JUMP
      REAL ALTSGN
      REAL TEMP
      SAVE ITER, J, JLAST, JUMP
      REAL ZERO, ONE
      PARAMETER( ZERO = 0.0E0 )
      PARAMETER( ONE = 1.0E0 )
      REAL, PARAMETER :: RZERO = 0.0E0
      REAL, PARAMETER :: RONE = 1.0E0
      IF (KASE .EQ. 0) THEN
        DO 10 I = 1, N
          X(I) = ONE / real(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        RETURN
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 20
      CASE(2)
        GOTO 40
      CASE(3)
        GOTO 70
      CASE(4)
        GOTO 120
      CASE(5)
        GOTO 160
      CASE DEFAULT
      END SELECT
   20 CONTINUE
      IF (N .EQ. 1) THEN
        W(1) = X(1)
        EST = abs(W(1))
        GOTO 190
      ENDIF
      DO 30 I = 1, N
        X(I)  = sign( RONE,real(X(I)) )
        IW(I) = nint(real(X(I)))
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
   40 CONTINUE
      J = SMUMPS_IXAMAX(N, X, 1)
      ITER = 2
   50 CONTINUE
      DO 60 I = 1, N
        X(I) = ZERO
   60 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
   70 CONTINUE
      DO 80 I = 1, N
        W(I) = X(I)
   80 CONTINUE
      DO 90 I = 1, N
        IF (nint(sign(RONE, real(X(I)))) .NE. IW(I)) GOTO 100
   90 CONTINUE
      GOTO 130
  100 CONTINUE
      DO 110 I = 1, N
        X(I) = sign(RONE, real(X(I)))
        IW(I) = nint(real(X(I)))
  110 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
  120 CONTINUE
      JLAST = J
      J = SMUMPS_IXAMAX(N, X, 1)
      IF ((abs(X(JLAST)) .NE. abs(X(J))) .AND. (ITER .LT. ITMAX)) THEN
        ITER = ITER + 1
        GOTO 50
      ENDIF
  130 CONTINUE
      EST = RZERO
      DO 140 I = 1, N
        EST = EST + abs(W(I))
  140 CONTINUE
      ALTSGN = RONE
      DO 150 I = 1, N
        X(I) = ALTSGN * (RONE + real(I - 1) / real(N - 1))
        ALTSGN = -ALTSGN
  150 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
  160 CONTINUE
      TEMP = RZERO
      DO 170 I = 1, N
        TEMP = TEMP + abs(X(I))
  170 CONTINUE
      TEMP = 2.0E0 * TEMP / real(3 * N)
      IF (TEMP .GT. EST) THEN
        DO 180 I = 1, N
          W(I) = X(I)
  180   CONTINUE
        EST = TEMP
      ENDIF
  190 KASE = 0
      RETURN
      END SUBROUTINE SMUMPS_218
      SUBROUTINE SMUMPS_278( MTYPE, N, NZ, ASPK, IRN, ICN,
     &    LHS, WRHS, W, RHS, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER MTYPE, N, NZ
      INTEGER IRN( NZ ), ICN( NZ )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      REAL, intent(in) :: ASPK( NZ )
      REAL, intent(in) :: LHS( N ), WRHS( N )
      REAL, intent(out):: RHS( N )
      REAL,    intent(out):: W( N )
      INTEGER K, I, J
      REAL DZERO
      PARAMETER(DZERO = 0.0E0)
      DO 10 K = 1, N
        W(K) = DZERO
        RHS(K) = WRHS(K)
   10 CONTINUE
      IF ( KEEP(50) .EQ. 0 ) THEN
       IF (MTYPE .EQ. 1) THEN
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          RHS(I) = RHS(I) - ASPK(K) * LHS(J)
          W(I) = W(I) + abs(ASPK(K))
        ENDDO
       ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          RHS(J) = RHS(J) - ASPK(K) * LHS(I)
          W(J) = W(J) + abs(ASPK(K))
        ENDDO
       ENDIF
      ELSE
       DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          RHS(I) = RHS(I) - ASPK(K) * LHS(J)
          W(I) = W(I) + abs(ASPK(K))
          IF (J.NE.I) THEN
            RHS(J) = RHS(J) - ASPK(K) * LHS(I)
            W(J) = W(J) + abs(ASPK(K))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_278
      SUBROUTINE SMUMPS_121( MTYPE, N, 
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     &    LHS, WRHS, W, RHS, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      REAL A_ELT(NA_ELT)
      REAL LHS( N ), WRHS( N ), RHS( N )
      REAL W(N)
      CALL SMUMPS_257(N, NELT, ELTPTR, ELTVAR, A_ELT,
     &                         LHS, RHS, KEEP(50), MTYPE )
      RHS = WRHS - RHS
      CALL SMUMPS_119( MTYPE, N, 
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     &    W, KEEP,KEEP8 )
      RETURN
      END SUBROUTINE SMUMPS_121
      SUBROUTINE SMUMPS_119( MTYPE, N, 
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     &    W, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      REAL A_ELT(NA_ELT)
      REAL TEMP
      REAL W(N)
      INTEGER K, I, J, IEL, SIZEI, IELPTR
      REAL DZERO
      PARAMETER(DZERO = 0.0E0)
      W = DZERO
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) = 
     &           W( ELTVAR( IELPTR + I) )
     &           + abs(A_ELT( K ))
               K = K + 1
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
               TEMP = TEMP + abs( A_ELT(K))
               K = K + 1
              END DO
              W(ELTVAR( IELPTR + J )) = 
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) = 
     &        W(ELTVAR( IELPTR + J )) + abs(A_ELT( K ))
          K = K + 1
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) = 
     &           W(ELTVAR( IELPTR + J )) + abs(A_ELT( K ))
              W(ELTVAR( IELPTR + I ) ) = 
     &           W(ELTVAR( IELPTR + I )) + abs(A_ELT( K ))
              K = K + 1
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_119
      SUBROUTINE SMUMPS_135(MTYPE, N, 
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     &    W, KEEP,KEEP8, COLSCA )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      REAL COLSCA(N)
      REAL A_ELT(NA_ELT)
      REAL W(N)
      REAL TEMP, TEMP2
      INTEGER K, I, J, IEL, SIZEI, IELPTR
      REAL DZERO
      PARAMETER(DZERO = 0.0E0)
      W = DZERO
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              TEMP2 = abs(COLSCA(ELTVAR( IELPTR + J) ))
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) =
     &           W( ELTVAR( IELPTR + I) )
     &           + abs(A_ELT( K )) * TEMP2
               K = K + 1
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              TEMP2= abs(COLSCA(ELTVAR( IELPTR + J) ))
              DO I = 1, SIZEI
               TEMP = TEMP + abs(A_ELT( K )) * TEMP2
               K = K + 1
              END DO
              W(ELTVAR( IELPTR + J )) =
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) =
     &        W(ELTVAR( IELPTR + J )) + 
     &        abs( A_ELT( K )*COLSCA(ELTVAR( IELPTR + J)) )
          K = K + 1
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) =
     &           W(ELTVAR( IELPTR + J )) + 
     &           abs(A_ELT( K )*COLSCA(ELTVAR( IELPTR + J)))
              W(ELTVAR( IELPTR + I ) ) =
     &           W(ELTVAR( IELPTR + I )) + 
     &           abs(A_ELT( K )*COLSCA(ELTVAR( IELPTR + I)))
              K = K + 1
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_135
      SUBROUTINE SMUMPS_122( MTYPE, N, NELT, ELTPTR, 
     &                     LELTVAR, ELTVAR, NA_ELT, A_ELT,
     &                     SAVERHS, X, Y, W, K50 )
      IMPLICIT NONE
      INTEGER N, NELT, K50, MTYPE, LELTVAR, NA_ELT
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( LELTVAR )
      REAL A_ELT( NA_ELT ), X( N ), Y( N ), 
     &                 SAVERHS(N)
      REAL W(N)
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      REAL ZERO
      REAL TEMP
      REAL TEMP2
      PARAMETER( ZERO = 0.0E0 )
      Y = SAVERHS
      W = ZERO
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( K50 .eq. 0 ) THEN
          IF ( MTYPE .eq. 1 ) THEN
            DO J = 1, SIZEI
              TEMP = X( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                Y( ELTVAR( IELPTR + I ) ) =
     &          Y( ELTVAR( IELPTR + I ) ) -
     &             A_ELT( K ) * TEMP
                W( ELTVAR( IELPTR + I ) ) =
     &          W( ELTVAR( IELPTR + I ) ) +
     &             abs( A_ELT( K ) * TEMP )
                K = K + 1
              END DO
            END DO
          ELSE
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              TEMP2 = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP - 
     &          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                TEMP2 = TEMP2 +  abs(
     &          A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
              W( ELTVAR( IELPTR + J ) ) = TEMP2
            END DO
          END IF
        ELSE
          DO J = 1, SIZEI
            Y( ELTVAR( IELPTR + J ) ) =
     &      Y( ELTVAR( IELPTR + J ) ) -
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            W( ELTVAR( IELPTR + J ) ) =
     &      W( ELTVAR( IELPTR + J ) ) + abs(
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
            K = K + 1
            DO I = J+1, SIZEI
              Y( ELTVAR( IELPTR + I ) ) =
     &        Y( ELTVAR( IELPTR + I ) ) -
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     &        Y( ELTVAR( IELPTR + J ) ) -
     &           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              W( ELTVAR( IELPTR + I ) ) =
     &        W( ELTVAR( IELPTR + I ) ) + abs(
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
              W( ELTVAR( IELPTR + J ) ) =
     &        W( ELTVAR( IELPTR + J ) ) + abs(
     &           A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_122
      SUBROUTINE SMUMPS_643(
     &     INODE,PTRFAC,KEEP,A,LA,STEP,
     &     KEEP8,N,MUST_BE_PERMUTED,IERR)
      USE SMUMPS_OOC
      IMPLICIT NONE
      INTEGER INODE,KEEP(500),N
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER IERR
      REAL A(LA)      
      INTEGER RETURN_VALUE
      LOGICAL MUST_BE_PERMUTED
      RETURN_VALUE=SMUMPS_726(INODE,PTRFAC,
     &     KEEP(28),A,LA,IERR)
      IF(RETURN_VALUE.EQ.OOC_NODE_NOT_IN_MEM)THEN
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
         CALL SMUMPS_578(INODE,PTRFAC,
     &        KEEP,KEEP8,A,IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
         CALL SMUMPS_577(
     &        A(PTRFAC(STEP(INODE))),
     &        INODE,IERR
     &        )
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ELSE
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ENDIF
      IF(RETURN_VALUE.NE.OOC_NODE_PERMUTED)THEN
         MUST_BE_PERMUTED=.TRUE.
         CALL SMUMPS_682(INODE)
      ELSE
         MUST_BE_PERMUTED=.FALSE.
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_643
      SUBROUTINE SMUMPS_257( N, NELT, ELTPTR, ELTVAR, A_ELT,
     &                         X, Y, K50, MTYPE )
      IMPLICIT NONE
      INTEGER N, NELT, K50, MTYPE
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( * )
      REAL A_ELT( * ), X( N ), Y( N )
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      REAL TEMP
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      Y = ZERO
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( K50 .eq. 0 ) THEN
          IF ( MTYPE .eq. 1 ) THEN
            DO J = 1, SIZEI
              TEMP = X( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                Y( ELTVAR( IELPTR + I ) ) =
     &          Y( ELTVAR( IELPTR + I ) ) +
     &             A_ELT( K ) * TEMP
                K = K + 1
              END DO
            END DO
          ELSE
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP + 
     &          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
            END DO
          END IF
        ELSE
          DO J = 1, SIZEI
            Y( ELTVAR( IELPTR + J ) ) =
     &      Y( ELTVAR( IELPTR + J ) ) +
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            K = K + 1
            DO I = J+1, SIZEI
              Y( ELTVAR( IELPTR + I ) ) =
     &        Y( ELTVAR( IELPTR + I ) ) +
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     &        Y( ELTVAR( IELPTR + J ) ) +
     &           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_257
      SUBROUTINE SMUMPS_192
     &( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     &  LDLT, MTYPE)
      IMPLICIT NONE
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      REAL A_loc( NZ_loc ), X( N ), Y_loc( N )
      INTEGER LDLT, MTYPE
      INTEGER I, J, K
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      Y_loc = ZERO
      IF ( LDLT .eq. 0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &          (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
        ENDDO
        ELSE
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     &        .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
        ENDDO
        END IF
      ELSE
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &        (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SMUMPS_192
      SUBROUTINE SMUMPS_256( N, NZ, IRN, ICN, ASPK, X, Y,
     &                         LDLT, MTYPE, MAXTRANS, PERM )
      INTEGER N, NZ, LDLT, MTYPE, MAXTRANS
      INTEGER IRN( NZ ), ICN( NZ ) 
      INTEGER PERM( N )
      REAL ASPK( NZ ), X( N ), Y( N )
      INTEGER K, I, J
      REAL PX( N )
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      Y = ZERO
      IF ( MAXTRANS .eq. 1 .and. MTYPE .eq. 1) THEN
        DO I = 1, N
          PX(I) = X( PERM( I ) )
        END DO
      ELSE
        PX = X
      END IF
      IF ( LDLT .eq. 0 ) THEN
       IF (MTYPE .EQ. 1) THEN
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
        ENDDO
       ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(J) = Y(J) + ASPK(K) * PX(I)
        ENDDO
       ENDIF
      ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
          IF (J.NE.I) THEN
            Y(J) = Y(J) + ASPK(K) * PX(I)
          ENDIF
        ENDDO
      END IF
      IF ( MAXTRANS .EQ. 1 .AND. MTYPE .eq. 0 ) THEN
      PX = Y
      DO I = 1, N
        Y( PERM( I ) ) = PX( I )
      END DO
      END IF
      RETURN
      END SUBROUTINE SMUMPS_256
      SUBROUTINE SMUMPS_193
     &( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     &  LDLT, MTYPE)
      IMPLICIT NONE
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      REAL A_loc( NZ_loc ), X( N )
      REAL Y_loc( N )
      INTEGER LDLT, MTYPE
      INTEGER I, J, K
      REAL RZERO
      PARAMETER( RZERO = 0.0E0 )
      Y_loc = RZERO
      IF ( LDLT .eq. 0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &          (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + abs( A_loc(K) * X(J) )
        ENDDO
        ELSE
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     &        .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(J) = Y_loc(J) + abs( A_loc(K) * X(I) )
        ENDDO
        END IF
      ELSE
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &        (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + abs( A_loc(K) * X(J) )
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + abs( A_loc(K) * X(I) )
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SMUMPS_193
