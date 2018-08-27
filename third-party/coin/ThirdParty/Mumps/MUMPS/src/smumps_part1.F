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
      SUBROUTINE SMUMPS( id )
      USE SMUMPS_OOC
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
C       matrix in assembled format (ICNTL(5)=0, and ICNTL(18) $\neq$ 3),
      INTERFACE
      SUBROUTINE SMUMPS_758
     &(idRHS, idINFO, idN, idNRHS, idLRHS)
      REAL, DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      END SUBROUTINE SMUMPS_758
      SUBROUTINE SMUMPS_26( id )
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET :: id
      END SUBROUTINE SMUMPS_26
      SUBROUTINE SMUMPS_142( id )
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET :: id
      END SUBROUTINE SMUMPS_142
      SUBROUTINE SMUMPS_301( id )
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET :: id
      END SUBROUTINE SMUMPS_301
      SUBROUTINE SMUMPS_349(id, LP)
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER  :: LP
      END SUBROUTINE SMUMPS_349
      END INTERFACE
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      TYPE (SMUMPS_STRUC) :: id
      INTEGER JOBMIN, JOBMAX, OLDJOB
      INTEGER I, J, MP, LP, MPG, KEEP235SAVE, KEEP242SAVE,
     &        KEEP243SAVE
      LOGICAL LANAL, LFACTO, LSOLVE, PROK, FLAG, PROKG
      LOGICAL NOERRORBEFOREPERM
      LOGICAL UNS_PERM_DONE
      INTEGER COMM_SAVE
      INTEGER JOB, N, NZ, NELT
      INTEGER, PARAMETER :: ICNTL18DIST_MIN = 1
      INTEGER, PARAMETER :: ICNTL18DIST_MAX = 3
      INTEGER, DIMENSION(:), ALLOCATABLE :: UNS_PERM_INV
      NOERRORBEFOREPERM = .FALSE.
      UNS_PERM_DONE = .FALSE.
      JOB  = id%JOB
      N    = id%N
      NZ   = id%NZ
      NELT = id%NELT
      id%INFO(1) = 0
      id%INFO(2) = 0
      IF ( JOB .NE. -1 ) THEN
        LP      = id%ICNTL(1)
        MP      = id%ICNTL(2)
        MPG     = id%ICNTL(3)
        PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.3))
        PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
        IF (PROKG) THEN
           IF (id%ICNTL(5) .NE. 1) THEN
              WRITE(MPG,'(A,I4,I12,I15)') 
     &             'Entering SMUMPS driver with JOB, N, NZ =', JOB,N,NZ
           ELSE
              WRITE(MPG,'(A,I4,I12,I15)') 
     &             'Entering SMUMPS driver with JOB, N, NELT =', JOB,N
     &             ,NELT
           ENDIF
        ENDIF
      ELSE
        MPG = 0
        PROK = .FALSE.
        PROKG = .FALSE.
        LP = 6
        MP = 6
      END IF
      CALL MPI_INITIALIZED( FLAG, IERR )
      IF ( .NOT. FLAG ) THEN
        WRITE(LP,990)
 990  FORMAT(' Error in SMUMPS initialization: MPI is not running.')
        id%INFO(1) = -23
        id%INFO(2) =   0
        GOTO 500
      END IF
       COMM_SAVE = id%COMM
       CALL MPI_COMM_DUP( COMM_SAVE, id%COMM, IERR )
      CALL MPI_ALLREDUCE(JOB,JOBMIN,1,MPI_INTEGER,MPI_MAX,
     &                   id%COMM,IERR)
      CALL MPI_ALLREDUCE(JOB,JOBMAX,1,MPI_INTEGER,MPI_MIN,
     &                   id%COMM,IERR)
      IF ( JOBMIN .NE. JOBMAX ) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
      IF ( JOB .EQ. -1 ) THEN
        id%INFO(1)=0
        id%INFO(2)=0
        IF ( id%KEEP(40) .EQ. 1 - 456789 .OR.
     &      id%KEEP(40) .EQ. 2 - 456789 .OR.
     &      id%KEEP(40) .EQ. 3 -456789 ) THEN
        IF ( id%N > 0 ) THEN
          id%INFO(1)=-3
          id%INFO(2)=JOB
        ENDIF
        ENDIF
        CALL MPI_COMM_RANK(id%COMM, id%MYID, IERR)
        CALL MUMPS_276( id%ICNTL(1),
     &                       id%INFO(1),
     &                       id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) THEN
           IF (id%KEEP(201).GT.0) THEN
             CALL SMUMPS_587(id, IERR)
           ENDIF
           GOTO 499
        ENDIF
        CALL SMUMPS_163( id )
        GOTO 500
      END IF
      IF ( JOB .EQ. -2 ) THEN
        id%KEEP(40)= -2 - 456789
        CALL SMUMPS_136( id )
        GOTO 500
      END IF
      IF ((JOB.LT.1).OR.(JOB.GT.6)) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        IF ( id%ICNTL(18) .LT. ICNTL18DIST_MIN
     &       .OR. id%ICNTL(18) .GT. ICNTL18DIST_MAX ) THEN
         IF ((N.LE.0).OR.((N+N+N)/3.NE.N)) THEN
          id%INFO(1) = -16
          id%INFO(2) = N
         END IF
         IF (id%ICNTL(5).NE.1) THEN
           IF (NZ.LE.0) THEN
             id%INFO(1) = -2
             id%INFO(2) = NZ
           END IF
         ELSE
           IF (NELT.LE.0) THEN
             id%INFO(1) = -24
             id%INFO(2) = NELT
           END IF
         ENDIF
        END IF
        IF ( (id%KEEP(46).EQ.0).AND.(id%NPROCS.LE.1) ) 
     &     THEN
          id%INFO(1) = -21
          id%INFO(2) = id%NPROCS
        ENDIF
      END IF
      CALL MUMPS_276( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 499
      LANAL  = .FALSE.
      LFACTO = .FALSE.
      LSOLVE = .FALSE.
      IF ((JOB.EQ.1).OR.(JOB.EQ.4).OR.
     &    (JOB.EQ.6))               LANAL  = .TRUE.
      IF ((JOB.EQ.2).OR.(JOB.EQ.4).OR.
     &    (JOB.EQ.5).OR.(JOB.EQ.6)) LFACTO = .TRUE.
      IF ((JOB.EQ.3).OR.(JOB.EQ.5).OR.
     &    (JOB.EQ.6))               LSOLVE = .TRUE.
      IF (MP.GT.0) CALL SMUMPS_349(id, MP)
      OLDJOB = id%KEEP( 40 ) + 456789
      IF ( LANAL ) THEN
        IF ( OLDJOB .EQ. 0 .OR. OLDJOB .GT. 3 .OR. OLDJOB .LT. -1 ) THEN
          id%INFO(1) = -3
          id%INFO(2) = JOB
          GOTO 499
        END IF
        IF ( OLDJOB .GE. 2 ) THEN
          IF (associated(id%IS)) THEN
            DEALLOCATE  (id%IS)
            NULLIFY     (id%IS)
          END IF
          IF (associated(id%S)) THEN
            DEALLOCATE  (id%S)
            NULLIFY     (id%S)
          END IF
        END IF   
      END IF
      IF ( LFACTO ) THEN
         IF ( OLDJOB .LT. 1 .and. .NOT. LANAL ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
      IF ( LSOLVE ) THEN
         IF ( OLDJOB .LT. 2 .AND. .NOT. LFACTO ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
#if ! defined (LARGEMATRICES)
      NOERRORBEFOREPERM =.TRUE.
      UNS_PERM_DONE=.FALSE.
      IF (id%MYID .eq. MASTER .AND. id%KEEP(23) .NE. 0) THEN
        IF ( id%JOB .EQ. 2 .OR. id%JOB .EQ. 5 .OR.
     &       (id%JOB .EQ. 3 .AND. (id%ICNTL(10) .NE.0 .OR.
     &        id%ICNTL(11).NE. 0))) THEN
          UNS_PERM_DONE = .TRUE.
          ALLOCATE(UNS_PERM_INV(id%N),stat=IERR)
          IF (IERR .GT. 0) THEN
              id%INFO(1)=-13
              id%INFO(2)=id%N
              IF (id%ICNTL(1) .GT.  0 .AND. id%ICNTL(4) .GE.1) THEN
                WRITE(id%ICNTL(2),99993)
              END IF
            GOTO 510
          ENDIF
          DO I = 1, id%N
            UNS_PERM_INV(id%UNS_PERM(I))=I
          END DO
          DO I = 1, id%NZ
            J = id%JCN(I)
            IF (J.LE.0.OR.J.GT.id%N) CYCLE
            id%JCN(I)=UNS_PERM_INV(J)
          END DO
          DEALLOCATE(UNS_PERM_INV)
        END IF
      END IF
#endif
        CALL MUMPS_276( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
      IF (LANAL) THEN
        id%KEEP(40)=-1 -456789
        IF (id%MYID.EQ.MASTER) THEN
          id%INFOG(7) = -9999
          id%INFOG(23) = 0
          id%INFOG(24) = 1
          IF (associated(id%IS1)) DEALLOCATE(id%IS1)
          IF ( id%ICNTL(5) .NE. 1 ) THEN  
             IF ( id%KEEP(50) .NE. 1 
     &            .AND. (
     &            (id%ICNTL(6) .NE. 0 .AND. id%ICNTL(7) .NE.1)
     &            .OR.
     &            id%ICNTL(12) .NE. 1) ) THEN
                id%MAXIS1 = 11 * N
             ELSE
              id%MAXIS1 = 10 * N
             END IF
          ELSE
            id%MAXIS1 = 6 * N + 2 * NELT + 2
          ENDIF
          ALLOCATE( id%IS1(id%MAXIS1), stat=IERR )
          IF (IERR.gt.0) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%MAXIS1
            IF ( LP .GT.0 ) 
     &      WRITE(LP,*) 'Problem in allocating work array for analysis.'
            GO TO 100
          END IF
          IF ( associated( id%PROCNODE ) )
     &          DEALLOCATE( id%PROCNODE )
          ALLOCATE( id%PROCNODE(id%N), stat=IERR )
          IF (IERR.gt.0) THEN
            id%INFO(1) = -7
            id%INFO(2) = id%N
            IF ( LP .GT. 0 ) THEN
              WRITE(LP,*) 'Problem in allocating work array PROCNODE'
            END IF
            GOTO 100
          END IF
          id%PROCNODE(1:id%N) = 0
          IF ( id%ICNTL(5) .EQ. 1 ) THEN 
            IF ( associated( id%ELTPROC ) )
     &            DEALLOCATE( id%ELTPROC )
            ALLOCATE( id%ELTPROC(id%NELT), stat=IERR )
            IF (IERR.gt.0) THEN
              id%INFO(1) = -7
              id%INFO(2) = id%NELT
              IF ( LP .GT. 0 ) THEN
                WRITE(LP,*) 'Problem in allocating work array ELTPROC'
              END IF
              GOTO 100
            END IF
          END IF
          IF ( id%ICNTL(5) .NE. 1 ) THEN 
            id%NA_ELT=0
            IF ( id%ICNTL(18) .LT. ICNTL18DIST_MIN
     &        .OR. id%ICNTL(18) .GT. ICNTL18DIST_MAX ) THEN
             IF ( .not. associated( id%IRN ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
             ELSE IF ( size( id%IRN ) < id%NZ ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
             ELSE IF ( .not. associated( id%JCN ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
             ELSE IF ( size( id%JCN ) < id%NZ ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
             END IF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF (LP.GT.0) WRITE(LP,*)
     &           'Error in analysis: IRN/JCN badly allocated.'
            END IF
          ELSE
            IF ( .not. associated( id%ELTPTR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( size( id%ELTPTR ) < id%NELT+1 ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( .not. associated( id%ELTVAR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
            ELSE 
              id%LELTVAR = id%ELTPTR( id%NELT+1 ) - 1
              IF ( size( id%ELTVAR ) < id%LELTVAR ) THEN 
                id%INFO(1) = -22
                id%INFO(2) = 2
              ELSE
                id%NA_ELT = 0
                IF ( id%KEEP(50) .EQ. 0 ) THEN
                  DO I = 1,NELT
                    J = id%ELTPTR(I+1) - id%ELTPTR(I)
                    J = (J * J)
                    id%NA_ELT = id%NA_ELT + J
                  ENDDO
                ELSE
                  DO I = 1,NELT
                    J = id%ELTPTR(I+1) - id%ELTPTR(I)
                    J = (J * (J+1))/2
                    id%NA_ELT = id%NA_ELT + J
                  ENDDO
                ENDIF
              ENDIF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF (LP.GT.0) WRITE(LP,*) 
     &           'Error in analysis: ELTPTR/ELTVAR badly allocated.'
            END IF
          ENDIF
 100    CONTINUE
        END IF
        CALL MUMPS_276( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
         id%KEEP(52) = id%ICNTL(8)
         IF ( id%KEEP(52) .GT. 8 .OR. id%KEEP(52).LT.-2)
     &        id%KEEP(52) = 77
         IF ((id%KEEP(52).EQ.77).AND.(id%KEEP(50).EQ.1)) THEN
           id%KEEP(52) = 0
         ENDIF
         IF ( id%KEEP(52).EQ.77 .OR. id%KEEP(52).LE.-2) THEN 
           IF (.not.associated(id%A)) id%KEEP(52) = 0
         ENDIF
         IF(id%KEEP(52) .EQ. -1) id%KEEP(52) = 0
         CALL SMUMPS_26( id )
        IF (id%MYID .eq. MASTER) THEN
           IF (id%KEEP(52) .NE. 0) THEN
             id%INFOG(33)=id%KEEP(52)
           ELSE
             id%INFOG(33)=id%ICNTL(8)
           ENDIF
        ENDIF
        IF (id%MYID .eq. MASTER) id%INFOG(24)=id%KEEP(95)
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
        id%KEEP(40) = 1 -456789
      END IF
      IF (LFACTO) THEN
         id%KEEP(40) = 1 - 456789
        IF ( id%MYID .EQ. MASTER ) THEN
          IF (id%KEEP(60).EQ.1) THEN
             IF ( associated( id%SCHUR_CINTERFACE)) THEN
               id%SCHUR=>id%SCHUR_CINTERFACE
     &          (1:id%SIZE_SCHUR*id%SIZE_SCHUR)
             ENDIF
             IF ( .NOT. associated (id%SCHUR)) THEN
              IF (LP.GT.0) 
     &        write(LP,'(A)') 
     &                      ' SCHUR not associated'
              id%INFO(1)=-22
              id%INFO(2)=9
             ELSE IF ( size(id%SCHUR) .LT.
     &                id%SIZE_SCHUR * id%SIZE_SCHUR ) THEN
                IF (LP.GT.0) 
     &          write(LP,'(A)') 
     &                ' SCHUR allocated but too small' 
                id%INFO(1)=-22
                id%INFO(2)=9
             END IF
          END IF
          IF ( id%KEEP(55) .EQ. 0 ) THEN
           IF ( id%KEEP(54).eq.0 ) THEN
            IF ( .not. associated( id%A ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            ELSE IF ( size( id%A ) < id%NZ ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            END IF
           END IF
          ELSE
            IF ( .not. associated( id%A_ELT ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            ELSE 
              IF ( size( id%A_ELT ) < id%NA_ELT ) THEN
                id%INFO( 1 ) = -22
                id%INFO( 2 ) = 4
              ENDIF
            END IF
          ENDIF
          CALL MUMPS_633(id%KEEP(12),id%ICNTL(14),
     &         id%KEEP(50),id%KEEP(54),id%ICNTL(6),id%ICNTL(8))
          CALL SMUMPS_635(N,id%KEEP(1),id%ICNTL(1),MPG)
          IF( id%KEEP(52) .EQ. -2 .AND. id%ICNTL(8) .NE. -2 .AND.
     &        id%ICNTL(8).NE. 77 ) THEN
             IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)') ' ** WARNING : SCALING'
                WRITE(MPG,'(A)') 
     &               ' ** scaling already computed during analysis'
                WRITE(MPG,'(A)') 
     &               ' ** keeping the scaling from the analysis'
             ENDIF
          ENDIF
          IF (id%KEEP(52) .NE. -2) THEN
            id%KEEP(52)=id%ICNTL(8)
          ENDIF
          IF ( id%KEEP(52) .GT. 8 .OR. id%KEEP(52).LT.-2)
     &    id%KEEP(52) = 77
          IF (id%KEEP(52).EQ.77) THEN
            IF (id%KEEP(50).EQ.1) THEN
              id%KEEP(52) = 0
            ELSE
              id%KEEP(52) = 7 
            ENDIF
          ENDIF
          IF (id%KEEP(23) .NE. 0 .AND. id%ICNTL(8) .EQ. -1) THEN
             IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)') ' ** WARNING : SCALING'
                WRITE(MPG,'(A)') 
     &               ' ** column permutation applied:'
                WRITE(MPG,'(A)') 
     &               ' ** column scaling has to be permuted'
             ENDIF 
          ENDIF
          IF ( id%KEEP( 19 ) .ne. 0 .and. id%KEEP( 52 ).ne. 0 ) THEN
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with null space)'
            END IF
            id%KEEP(52) = 0
          END IF
          IF ( id%KEEP(60) .ne. 0 .and. id%KEEP(52) .ne. 0 ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 .AND. id%ICNTL(8) .NE. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with Schur)'
            END IF
          END IF
          IF (id%KEEP(54) .NE. 0 .AND. 
     &        id%KEEP(52).NE.7 .AND. id%KEEP(52).NE.8 .AND.
     &        id%KEEP(52) .NE. 0 ) THEN
             id%KEEP(52) = 0
             IF ( MPG .GT. 0 .and. id%ICNTL(8) .ne. 0 ) THEN
               WRITE(MPG,'(A)')
     &         ' ** Warning: This scaling option not available'
               WRITE(MPG,'(A)') ' ** for distributed matrix entry'
             END IF
          END IF
          IF ( id%KEEP(50) .NE. 0 ) THEN
             IF ( id%KEEP(52).ne.  1 .and.
     &            id%KEEP(52).ne. -1 .and.
     &            id%KEEP(52).ne.  0 .and.
     &            id%KEEP(52).ne.  7 .and.
     &            id%KEEP(52).ne.  8 .and.
     &            id%KEEP(52).ne. -2 .and.
     &            id%KEEP(52).ne. 77) THEN
              IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)')
     &  ' ** Warning: Scaling option n.a. for symmetric matrix'
              END IF
              id%KEEP(52) = 0
            END IF
          END IF
          IF (id%KEEP(55) .NE. 0 .AND. 
     &        ( id%KEEP(52) .gt. 0 ) ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)')
     &        ' ** (only user scaling av. for elt. entry)'
            END IF
          END IF
          IF ( id%KEEP(52) .eq. -1 ) THEN
            IF ( .not. associated( id%ROWSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( size( id%ROWSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( .not. associated( id%COLSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            ELSE IF ( size( id%COLSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            END IF
          END IF
          IF (id%KEEP(52).GT.0 .AND.
     &        id%KEEP(52) .LE.8) THEN
            IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=IERR)
            IF (IERR .GT.0) id%INFO(1)=-13
            ALLOCATE( id%ROWSCA(N), stat=IERR)
            IF (IERR .GT.0) id%INFO(1)=-13
          END IF
          IF (.NOT. associated(id%COLSCA)) THEN
            ALLOCATE( id%COLSCA(1), stat=IERR)
          END IF
          IF (IERR .GT.0) id%INFO(1)=-13
          IF (.NOT. associated(id%ROWSCA))
     &    ALLOCATE( id%ROWSCA(1), stat=IERR)
          IF (IERR .GT.0) id%INFO(1)=-13
          IF ( id%INFO(1) .eq. -13 ) THEN
            IF ( LP .GT. 0 )
     &      WRITE(LP,*) 'Problems in allocations before facto'
            GOTO 200
          END IF
          IF (id%KEEP(252) .EQ. 1) THEN
             CALL SMUMPS_758
     &       (id%RHS,id%INFO,id%N,id%NRHS,id%LRHS)
             CALL SMUMPS_807(id)
             CALL SMUMPS_769(id)
          ENDIF
 200      CONTINUE
        END IF        
        CALL MPI_BCAST( id%KEEP(221), 1, MPI_INTEGER, MASTER, id%COMM,
     &                  IERR )
        IF (id%KEEP(60).EQ.2.OR.id%KEEP(60).EQ.3) THEN
          IF ( id%root%yes ) THEN
            IF ( associated( id%SCHUR_CINTERFACE )) THEN
              id%SCHUR=>id%SCHUR_CINTERFACE
     &          (1:id%SCHUR_LLD*(id%root%SCHUR_NLOC-1)+
     &          id%root%SCHUR_MLOC)
            ENDIF
            IF (id%SCHUR_LLD < id%root%SCHUR_MLOC) THEN
              IF (LP.GT.0) write(LP,*) 
     &          ' SCHUR leading dimension SCHUR_LLD ', 
     &          id%SCHUR_LLD, 'too small with respect to', 
     &          id%root%SCHUR_MLOC
              id%INFO(1)=-30
              id%INFO(2)=id%SCHUR_LLD
            ELSE IF ( .NOT. associated (id%SCHUR)) THEN
              IF (LP.GT.0) write(LP,'(A)') 
     &                      ' SCHUR not associated'
              id%INFO(1)=-22
              id%INFO(2)=9
            ELSE IF (size(id%SCHUR) <
     &          id%SCHUR_LLD*(id%root%SCHUR_NLOC-1)+
     &          id%root%SCHUR_MLOC) THEN
              IF (LP.GT.0) THEN 
                write(LP,'(A)') 
     &                      ' SCHUR allocated but too small'
                write(LP,*) id%MYID, ' : Size Schur=', 
     &          size(id%SCHUR), 
     &          ' SCHUR_LLD= ', id%SCHUR_LLD, 
     &          ' SCHUR_MLOC=', id%root%SCHUR_NLOC, 
     &          ' SCHUR_NLOC=', id%root%SCHUR_NLOC
              ENDIF
              id%INFO(1)=-22
              id%INFO(2)= 9
            ELSE
               id%root%SCHUR_LLD=id%SCHUR_LLD
               IF (id%root%SCHUR_NLOC==0) THEN
                 ALLOCATE(id%root%SCHUR_POINTER(1))
               ELSE
                id%root%SCHUR_POINTER=>id%SCHUR
               ENDIF
            ENDIF
          ENDIF
        ENDIF
        CALL MUMPS_276( id%ICNTL(1),
     &                      id%INFO(1),
     &                      id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
        CALL SMUMPS_142(id)
        IF (id%MYID .eq. MASTER) id%INFOG(33)=id%KEEP(52)
        IF (id%KEEP(60).EQ.2.OR.id%KEEP(60).EQ.3) THEN
          IF (id%root%yes) THEN
            IF (id%root%SCHUR_NLOC==0) THEN
               DEALLOCATE(id%root%SCHUR_POINTER)
               NULLIFY(id%root%SCHUR_POINTER)
            ELSE
               NULLIFY(id%root%SCHUR_POINTER)
            ENDIF
          ENDIF
        ENDIF
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
        id%KEEP(40) = 2 - 456789
      END IF
      IF (LSOLVE) THEN
        id%KEEP(40) = 2 -456789
        IF (id%MYID .eq. MASTER) THEN
           KEEP235SAVE = id%KEEP(235)
           KEEP242SAVE = id%KEEP(242)
           KEEP243SAVE = id%KEEP(243)
           IF (id%KEEP(242).EQ.0) id%KEEP(243)=0
        ENDIF
        CALL SMUMPS_301(id)
        IF (id%MYID .eq. MASTER) THEN
           id%KEEP(235) = KEEP235SAVE
           id%KEEP(242) = KEEP242SAVE
           id%KEEP(243) = KEEP243SAVE
        ENDIF
        IF (id%INFO(1).LT.0) GOTO 499
        id%KEEP(40) = 3 -456789
      ENDIF
      IF (MP.GT.0) CALL SMUMPS_349(id, MP)
      GOTO 500
  499 PROK  = ((id%ICNTL(1).GT.0).AND.
     &         (id%ICNTL(4).GE.1))
      IF (PROK) WRITE (id%ICNTL(1),99995) id%INFO(1)
      IF (PROK) WRITE (id%ICNTL(1),99994) id%INFO(2)
500   CONTINUE
#if ! defined(LARGEMATRICES)
      IF (id%MYID .eq. MASTER .AND. id%KEEP(23) .NE. 0
     &    .AND. NOERRORBEFOREPERM) THEN
        IF (id%JOB .NE. 3 .OR. UNS_PERM_DONE) THEN
          DO I = 1, id%NZ
            J=id%JCN(I)
            IF (J.LE.0.OR.J.GT.id%N) CYCLE
            id%JCN(I)=id%UNS_PERM(J)
          END DO
        END IF
      END IF
#endif
 510  CONTINUE
      CALL SMUMPS_300(id%INFO(1), id%INFOG(1), id%COMM, id%MYID)
      CALL MPI_BCAST( id%RINFOG(1), 40, MPI_REAL, MASTER,
     &                    id%COMM, IERR )
      IF (id%MYID.EQ.MASTER.and.MPG.GT.0.and.
     & id%INFOG(1).lt.0) THEN
        WRITE(MPG,'(A,I12)') ' On return from SMUMPS, INFOG(1)=',
     &      id%INFOG(1)
        WRITE(MPG,'(A,I12)') ' On return from SMUMPS, INFOG(2)=',
     &      id%INFOG(2)
      END IF
       CALL MPI_COMM_FREE( id%COMM, IERR )
       id%COMM = COMM_SAVE
      RETURN
99995 FORMAT (' ** ERROR RETURN ** FROM SMUMPS INFO(1)=', I3)
99994 FORMAT (' ** INFO(2)=', I10)
99993 FORMAT (' ** Allocation error: could not permute JCN.')
      END SUBROUTINE SMUMPS
      SUBROUTINE SMUMPS_300( INFO, INFOG, COMM, MYID )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INFO(40), INFOG(40), COMM, MYID
      INTEGER TMP1(2),TMP(2)
      INTEGER ROOT, IERR
      INTEGER MASTER
      PARAMETER (MASTER=0)
      IF ( INFO(1) .ge. 0  .and. INFO(2) .ge. 0 ) THEN
        INFOG(1) = INFO(1)
        INFOG(2) = INFO(2)
      ELSE
        INFOG(1) = INFO(1)
        TMP1(1) = INFO(1)
        TMP1(2) = MYID
        CALL MPI_ALLREDUCE(TMP1,TMP,1,MPI_2INTEGER,
     &                     MPI_MINLOC,COMM,IERR )
        INFOG(2) = INFO(2)
        ROOT = TMP(2)
        CALL MPI_BCAST( INFOG(1), 1, MPI_INTEGER, ROOT, COMM, IERR )
        CALL MPI_BCAST( INFOG(2), 1, MPI_INTEGER, ROOT, COMM, IERR )
      END IF
      CALL MPI_BCAST(INFOG(3), 38, MPI_INTEGER, MASTER, COMM, IERR )
      RETURN
      END SUBROUTINE SMUMPS_300
      SUBROUTINE SMUMPS_349(id, LP)
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER  :: LP
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LT.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      IF (id%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19),ICNTL(22)
           IF ((ICNTL(6).EQ.5).OR.(ICNTL(6).EQ.6).OR.
     &          (ICNTL(12).NE.1) )  THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,993) ICNTL(14)
         CASE(2);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,993) ICNTL(14)
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,992) ICNTL(8)
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,993) ICNTL(14)
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19),ICNTL(22)
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,993) ICNTL(14)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19),ICNTL(22)
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
           WRITE (LP,993) ICNTL(14)
        END SELECT
      ENDIF
 980  FORMAT (/'***********CONTROL PARAMETERS (ICNTL)**************'/)
 990  FORMAT (
     &     'ICNTL(1)   Output stream for error messages        =',I10/
     &     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     &     'ICNTL(3)   Output stream for global information    =',I10/
     &     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     &     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     &     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     &     'ICNTL(7)   Ordering                                =',I10/
     &     'ICNTL(12)  LDLT ordering strat ( keep(95) )        =',I10/
     &     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     &     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     &     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10/
     &     'ICNTL(22)  Out-off-core option (0=Off, >0=ON)      =',I10)
 992  FORMAT (
     &     'ICNTL(8)   Scaling strategy                        =',I10)
 993  FORMAT (
     &     'ICNTL(14)  Percent of memory increase              =',I10)
 995  FORMAT (
     &     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)       =',I10/
     &     'ICNTL(10)  Max steps iterative refinement          =',I10/
     &     'ICNTL(11)  Error analysis ( 0= off, else=on)       =',I10)
 998  FORMAT (
     &     '      Size of SCHUR matrix (SIZE_SHUR)             =',I10)
      END SUBROUTINE SMUMPS_349
      SUBROUTINE SMUMPS_350(id, LP)
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER ::LP
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL, KEEP
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LT.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      KEEP=>id%KEEP
      IF (id%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6))THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(2);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4) 
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).NE.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF  
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
           WRITE (LP,993) KEEP(12)
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF              
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),KEEP(248),ICNTL(21)
           WRITE (LP,993) KEEP(12)
        END SELECT
      ENDIF
 980  FORMAT (/'******INTERNAL VALUE OF PARAMETERS (ICNTL/KEEP)****'/)
 990  FORMAT (
     &     'ICNTL(1)   Output stream for error messages        =',I10/
     &     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     &     'ICNTL(3)   Output stream for global information    =',I10/
     &     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     &     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     &     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     &     'ICNTL(7)   Ordering                                =',I10/
     &     'ICNTL(12)  LDLT ordering strat ( keep(95) )        =',I10/
     &     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     &     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     &     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10/
     &     'ICNTL(22)  Out-off-core option (0=Off, >0=ON)      =',I10)
 992  FORMAT (
     &     'ICNTL(8)   Scaling strategy ( keep(52) )           =',I10)
 993  FORMAT (
     &     'ICNTL(14)  Percent of memory increase ( keep(12) ) =',I10)
 995  FORMAT (
     &     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)      =',I10/
     &     'ICNTL(10)  Max steps iterative refinement          =',I10/
     &     'ICNTL(11)  Error analysis ( 0= off, else=on)       =',I10/
     &     'ICNTL(20)  Dense (0) or sparse (1) RHS             =',I10/
     &     'ICNTL(21)  Gathered (0) or distributed(1) solution =',I10)
      END SUBROUTINE SMUMPS_350
      SUBROUTINE SMUMPS_758
     &       (idRHS, idINFO, idN, idNRHS, idLRHS)
      IMPLICIT NONE
      REAL, DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      IF ( .not. associated( idRHS ) ) THEN
              idINFO( 1 ) = -22
              idINFO( 2 ) = 7
      ELSE IF (idNRHS.EQ.1) THEN
               IF ( size( idRHS ) < idN ) THEN
                  idINFO( 1 ) = -22
                  idINFO( 2 ) = 7
               ENDIF
      ELSE IF (idLRHS < idN) 
     &            THEN
                  idINFO( 1 ) = -26
                  idINFO( 2 ) = idLRHS
      ELSE IF 
     &      (size(idRHS)<(idNRHS*idLRHS-idLRHS+idN)) 
     &            THEN
                  idINFO( 1 ) = -22
                  idINFO( 2 ) = 7
      END IF
      RETURN
      END SUBROUTINE SMUMPS_758
      SUBROUTINE SMUMPS_807(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (SMUMPS_STRUC) :: id
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (id%MYID.EQ.MASTER) THEN
        id%KEEP(221)=id%ICNTL(26)
        IF (id%KEEP(221).ne.0 .and. id%KEEP(221) .NE.1
     &      .AND.id%KEEP(221).ne.2) id%KEEP(221)=0
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_807
      SUBROUTINE SMUMPS_769(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (SMUMPS_STRUC) :: id
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (id%MYID .EQ. MASTER) THEN
          IF ( id%KEEP(221) == 1 .or. id%KEEP(221) == 2 ) THEN
            IF (id%KEEP(221) == 2 .and. id%JOB == 2) THEN
              id%INFO(1)=-35
              id%INFO(2)=id%KEEP(221)
              GOTO 333
            ENDIF
            IF (id%KEEP(221) == 1 .and. id%KEEP(252) == 1
     &          .and. id%JOB == 3) THEN
              id%INFO(1)=-35
              id%INFO(2)=id%KEEP(221)
            ENDIF
            IF ( id%KEEP(60).eq. 0 .or. id%SIZE_SCHUR.EQ.0 ) THEN
              id%INFO(1)=-33
              id%INFO(2)=id%KEEP(221)
              GOTO 333
            ENDIF
            IF ( .NOT. associated( id%REDRHS)) THEN
              id%INFO(1)=-22
              id%INFO(2)=15
              GOTO 333
            ELSE IF (id%NRHS.EQ.1) THEN
              IF (size(id%REDRHS) < id%SIZE_SCHUR ) THEN
                id%INFO(1)=-22
                id%INFO(2)=15
                GOTO 333
              ENDIF
            ELSE IF (id%LREDRHS < id%SIZE_SCHUR) THEN
              id%INFO(1)=-34
              id%INFO(2)=id%LREDRHS
              GOTO 333
            ELSE IF
     &      (size(id%REDRHS)<
     &         id%NRHS*id%LREDRHS-id%LREDRHS+id%SIZE_SCHUR)
     &      THEN
              id%INFO(1)=-22
              id%INFO(2)=15
              GOTO 333
            ENDIF
          ENDIF
      ENDIF
 333  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_769
      SUBROUTINE SMUMPS_24( MYID, SLAVEF, N,
     &           PROCNODE, STEP, PTRAIW, PTRARW, ISTEP_TO_INIV2,
     &           I_AM_CAND,
     &           KEEP, KEEP8, ICNTL, id )
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (SMUMPS_STRUC) :: id
      INTEGER MYID, N, SLAVEF
      INTEGER KEEP( 500 ), ICNTL( 40 )
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE( KEEP(28) ), STEP( N ),
     &        PTRAIW( N ), PTRARW( N )
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      LOGICAL I_AM_CAND(max(1,KEEP(56)))
      LOGICAL I_AM_SLAVE
      LOGICAL I_AM_CAND_LOC
      INTEGER MUMPS_330, MUMPS_275, MUMPS_810
      EXTERNAL MUMPS_330, MUMPS_275, MUMPS_810
      INTEGER ISTEP, I, IPTRI, IPTRR, NCOL, NROW, allocok
      INTEGER TYPE_PARALL, ITYPE, IRANK, INIV2, TYPESPLIT 
      LOGICAL T4_MASTER_CONCERNED
      TYPE_PARALL = KEEP(46)
      I_AM_SLAVE = (KEEP(46).EQ.1 .OR. MYID.NE.0)
      KEEP(14) = 0
      KEEP(13) = 0
      DO I = 1, N
        ISTEP=abs(STEP(I))
        ITYPE = MUMPS_330( PROCNODE(ISTEP), SLAVEF )
        IRANK = MUMPS_275( PROCNODE(ISTEP), SLAVEF )
        TYPESPLIT = MUMPS_810 ( PROCNODE(ISTEP), SLAVEF )
        I_AM_CAND_LOC = .FALSE.
        T4_MASTER_CONCERNED = .FALSE.
        IF (ITYPE.EQ.2) THEN
         INIV2         = ISTEP_TO_INIV2(ISTEP)
         IF (I_AM_SLAVE)  THEN 
           I_AM_CAND_LOC = I_AM_CAND(INIV2)
          IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
           IF ( TYPE_PARALL .eq. 0 ) THEN
            T4_MASTER_CONCERNED = 
     &     ( id%CANDIDATES (id%CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
     &       .EQ.MYID-1 )
           ELSE
            T4_MASTER_CONCERNED = 
     &     ( id%CANDIDATES (id%CANDIDATES(SLAVEF+1, INIV2)+1,INIV2 ) 
     &       .EQ.MYID )
           ENDIF
          ENDIF
         ENDIF
        ENDIF
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK = IRANK + 1
        END IF
        IF (
     &       ( (ITYPE .EQ. 1.OR.ITYPE.EQ.2) .AND.
     &            IRANK .EQ. MYID ) 
     &       .OR.
     &       ( T4_MASTER_CONCERNED ) 
     &     ) THEN
          KEEP( 14 ) = KEEP( 14 ) + 3 + PTRAIW( I ) + PTRARW( I )
          KEEP( 13 ) = KEEP( 13 ) + 1 + PTRAIW( I ) + PTRARW( I )
        ELSE IF ( ITYPE .EQ. 3 ) THEN
        ELSE IF ( ITYPE .EQ. 2 .AND. I_AM_CAND_LOC ) THEN
           PTRARW( I ) = 0
           KEEP(14) = KEEP(14) + 3 + PTRAIW( I ) + PTRARW( I )
           KEEP(13) = KEEP(13) + 1 + PTRAIW( I ) + PTRARW( I )
        END IF
      END DO
      IF ( associated( id%INTARR ) ) THEN
        DEALLOCATE( id%INTARR )
        NULLIFY( id%INTARR )
      END IF
      IF ( KEEP(14) > 0 ) THEN
      ALLOCATE( id%INTARR( KEEP(14) ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        id%INFO(1) = -7
        id%INFO(2) = KEEP(14)
        RETURN
      END IF
      ELSE
      ALLOCATE( id%INTARR( 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        id%INFO(1) = -7
        id%INFO(2) = 1
        RETURN
      END IF
      END IF
      IPTRI = 1
      IPTRR = 1
      DO I = 1, N
        ISTEP = abs(STEP(I))
        ITYPE = MUMPS_330( PROCNODE(ISTEP), SLAVEF )
        IRANK = MUMPS_275( PROCNODE(ISTEP), SLAVEF )
        TYPESPLIT = MUMPS_810 ( PROCNODE(ISTEP), SLAVEF )
        I_AM_CAND_LOC = .FALSE.
        T4_MASTER_CONCERNED = .FALSE.
        IF (ITYPE.EQ.2) THEN
          INIV2         = ISTEP_TO_INIV2(ISTEP)
          IF (I_AM_SLAVE)  THEN
           I_AM_CAND_LOC = I_AM_CAND(INIV2)
           IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
            IF ( TYPE_PARALL .eq. 0 ) THEN
             T4_MASTER_CONCERNED = 
     &       (id%CANDIDATES (id%CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
     &         .EQ.MYID-1 )
            ELSE
              T4_MASTER_CONCERNED = 
     &        (id%CANDIDATES (id%CANDIDATES(SLAVEF+1,INIV2)+1,INIV2) 
     &         .EQ.MYID )
            ENDIF
           ENDIF
          ENDIF
        ENDIF
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK =IRANK + 1
        END IF
        IF (
     &      ( ITYPE .eq. 2 .and.
     &        IRANK .eq. MYID )
     & .or.
     &      ( ITYPE .eq. 1 .and.
     &        IRANK .eq. MYID )
     & .or.
     &      ( T4_MASTER_CONCERNED )
     &     )  THEN
          NCOL = PTRAIW( I )
          NROW = PTRARW( I )
          id%INTARR( IPTRI     ) = NCOL
          id%INTARR( IPTRI + 1 ) = -NROW
          id%INTARR( IPTRI + 2 ) = I
          PTRAIW( I ) = IPTRI
          PTRARW( I ) = IPTRR
          IPTRI = IPTRI + NCOL + NROW + 3
          IPTRR = IPTRR + NCOL + NROW + 1
        ELSE IF ( ITYPE .eq. 2  .AND. I_AM_CAND_LOC ) THEN
           NCOL = PTRAIW( I )
           NROW = 0
           id%INTARR( IPTRI     ) = NCOL
           id%INTARR( IPTRI + 1 ) = -NROW
           id%INTARR( IPTRI + 2 ) = I
           PTRAIW( I ) = IPTRI
           PTRARW( I ) = IPTRR
           IPTRI = IPTRI + NCOL + NROW + 3
           IPTRR = IPTRR + NCOL + NROW + 1
        ELSE
          PTRAIW(I) = 0
          PTRARW(I) = 0
        END IF
      END DO
      IF ( IPTRI - 1 .NE. KEEP(14) ) THEN
        WRITE(*,*) 'Error 1 in anal_arrowheads',  
     &      ' IPTRI - 1, KEEP(14)=', IPTRI - 1, KEEP(14)
        CALL MUMPS_ABORT()
      END IF
      IF ( IPTRR - 1 .NE. KEEP(13) ) THEN
        WRITE(*,*) 'Error 2 in anal_arrowheads'
        CALL MUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE SMUMPS_24
      SUBROUTINE SMUMPS_148(N, NZ, ASPK, 
     &   IRN, ICN, PERM,
     &   LSCAL,COLSCA,ROWSCA,
     &   MYID, SLAVEF, PROCNODE_STEPS, NBRECORDS,
     &   LP, COMM, root, KEEP, KEEP8, FILS, RG2L,
     &   INTARR, DBLARR, PTRAIW, PTRARW, FRERE_STEPS,
     &   STEP, A, LA, ISTEP_TO_INIV2, I_AM_CAND, CANDIDATES )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INTEGER N,NZ, COMM, NBRECORDS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      REAL ASPK(NZ)
      REAL COLSCA(*), ROWSCA(*)
      INTEGER IRN(NZ), ICN(NZ) 
      INTEGER PERM(N), PROCNODE_STEPS(KEEP(28))
      INTEGER RG2L( N ), FILS( N )
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      LOGICAL I_AM_CAND(max(1,KEEP(56)))
      INTEGER LP, SLAVEF, MYID
      INTEGER CANDIDATES(SLAVEF+1, max(1,KEEP(56)))
      LOGICAL LSCAL
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER(8) :: LA
      INTEGER PTRAIW( N ), PTRARW( N ), FRERE_STEPS( KEEP(28) )
      INTEGER STEP(N)
      INTEGER INTARR( max(1,KEEP(14)) )
      REAL A( LA ), DBLARR(max(1,KEEP(13)))
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: BUFI
      REAL, DIMENSION(:,:), ALLOCATABLE :: BUFR
      INTEGER MUMPS_275, MUMPS_330, numroc, 
     &        MUMPS_810
      EXTERNAL MUMPS_275, MUMPS_330, numroc, 
     &        MUMPS_810
      REAL VAL
      INTEGER IOLD,JOLD,INEW,JNEW,ISEND,JSEND,DEST,I,K,IARR
      INTEGER IPOSROOT, JPOSROOT
      INTEGER IROW_GRID, JCOL_GRID
      INTEGER INODE, ISTEP
      INTEGER NBUFS
      INTEGER ARROW_ROOT, TAILLE
      INTEGER LOCAL_M, LOCAL_N
      INTEGER(8) :: PTR_ROOT
      INTEGER TYPENODE_TMP, MASTER_NODE
      LOGICAL I_AM_CAND_LOC, I_AM_SLAVE
      INTEGER I1, IA, JARR, ILOCROOT, JLOCROOT
      INTEGER IS1, ISHIFT, IIW, IS, IAS
      INTEGER allocok, TYPESPLIT, T4MASTER, INIV2
      LOGICAL T4_MASTER_CONCERNED
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      ARROW_ROOT = 0
      I_AM_SLAVE=(MYID.NE.0.OR.KEEP(46).EQ.1)
      IF ( KEEP(46) .eq. 0 ) THEN
        NBUFS = SLAVEF
      ELSE
        NBUFS = SLAVEF - 1
        ALLOCATE( IW4( N, 2 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          WRITE(*,*) 'Error allocating IW4'
          CALL MUMPS_ABORT()
        END IF
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
              root%SCHUR_POINTER(int(I-1,8)*int(root%SCHUR_LLD,8)+1_8:
     &        int(I-1,8)*int(root%SCHUR_LLD,8)+int(root%SCHUR_MLOC,8))=
     &        ZERO
            ENDDO
          ENDIF
        END IF
      END IF
      IF (NBUFS.GT.0) THEN
       ALLOCATE( BUFI(NBRECORDS*2+1,NBUFS),stat=allocok )
       IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) 'Error allocating BUFI'
        CALL MUMPS_ABORT()
       END IF
       ALLOCATE( BUFR( NBRECORDS, NBUFS ), stat=allocok )
       IF ( allocok .GT. 0 ) THEN
         WRITE(*,*) 'Error allocating BUFR'
         CALL MUMPS_ABORT()
       END IF
       DO I = 1, NBUFS
        BUFI( 1, I ) = 0
       ENDDO
      ENDIF
      INODE = KEEP(38)
      I     = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = I
        INODE = FILS( INODE )
        I = I + 1
      END DO
      DO 120 K=1,NZ
        IOLD = IRN(K)
        JOLD = ICN(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     &                 .OR.(JOLD.LT.1) ) THEN
           GOTO 120
        END IF
        IF (LSCAL) THEN
          VAL = ASPK(K)*ROWSCA(IOLD)*COLSCA(JOLD)
        ELSE
          VAL = ASPK(K)
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
        IARR  = abs( ISEND )
        ISTEP = abs( STEP(IARR) )
        TYPENODE_TMP = MUMPS_330( PROCNODE_STEPS(ISTEP),
     &                                 SLAVEF ) 
        MASTER_NODE  = MUMPS_275( PROCNODE_STEPS(ISTEP),
     &                                 SLAVEF )
        TYPESPLIT  = MUMPS_810(  PROCNODE_STEPS(ISTEP),
     &                                 SLAVEF )
        I_AM_CAND_LOC          = .FALSE.
        T4_MASTER_CONCERNED = .FALSE.
        T4MASTER               = -9999
        IF (TYPENODE_TMP.EQ.2) THEN
         INIV2         = ISTEP_TO_INIV2(ISTEP)
         IF (I_AM_SLAVE)  I_AM_CAND_LOC = I_AM_CAND(INIV2)
         IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
          T4_MASTER_CONCERNED = .TRUE.
          T4MASTER=CANDIDATES(CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
          IF ( KEEP(46) .eq. 0 ) THEN
           T4MASTER=T4MASTER+1
          ENDIF
         ENDIF
        ENDIF
        IF ( TYPENODE_TMP .EQ. 1 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = MASTER_NODE + 1
          ELSE
            DEST = MASTER_NODE
          END IF
        ELSE IF ( TYPENODE_TMP .EQ. 2 ) THEN
          IF ( ISEND .LT. 0  ) THEN
            DEST = -1
          ELSE
            IF ( KEEP( 46 ) .eq. 0 ) THEN
              DEST = MASTER_NODE + 1
            ELSE 
              DEST = MASTER_NODE
            END IF
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L(JSEND)
            JPOSROOT = RG2L(IARR)
          ELSE
            IPOSROOT = RG2L( IARR )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * root%NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * root%NPCOL + JCOL_GRID
          END IF
        END IF
        IF ( DEST .eq. 0 .or.
     &      ( DEST .eq. -1 .and. KEEP( 46 ) .eq. 1 .AND.
     &       ( I_AM_CAND_LOC .OR. MASTER_NODE .EQ. 0 ) )
     &     .or. 
     &      ( T4MASTER.EQ.0 )
     &     ) THEN
          IARR = ISEND  
          JARR = JSEND
          IF ( TYPENODE_TMP .eq. 3 ) THEN
            ARROW_ROOT = ARROW_ROOT + 1
            IF ( IROW_GRID .EQ. root%MYROW .AND.
     &         JCOL_GRID .EQ. root%MYCOL ) THEN
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
             IF (KEEP(60)==0) THEN
               A( PTR_ROOT
     &           + int(JLOCROOT - 1,8) * int(LOCAL_M,8) 
     &           + int(ILOCROOT - 1,8) )
     &         =  A( PTR_ROOT
     &           + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &           + int(ILOCROOT - 1,8) )
     &         + VAL
             ELSE
               root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                           * int(root%SCHUR_LLD,8)
     &                           + int(ILOCROOT,8) )
     &          = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                           *    int(root%SCHUR_LLD,8)
     &                           +    int(ILOCROOT,8))
     &          + VAL
             ENDIF
            ELSE
              WRITE(*,*) MYID,':INTERNAL Error: root arrowhead '
              WRITE(*,*) MYID,':is not belonging to me. IARR,JARR='
     &        ,IARR,JARR
              CALL MUMPS_ABORT()
            END IF
          ELSE IF ( IARR .GE. 0 ) THEN
            IF ( IARR .eq. JARR ) THEN
              IA = PTRARW( IARR )
              DBLARR( IA ) = DBLARR( IA ) + VAL
            ELSE
              IS1 =  PTRAIW(IARR)
              ISHIFT      = INTARR(IS1) + IW4(IARR,2)
              IW4(IARR,2) = IW4(IARR,2) - 1
              IIW         = IS1 + ISHIFT + 2
              INTARR(IIW)     = JARR
              IS          = PTRARW(IARR)
              IAS         = IS + ISHIFT
              DBLARR(IAS) = VAL
            END IF
          ELSE
            IARR = -IARR
            ISHIFT      = PTRAIW(IARR)+IW4(IARR,1)+2
            INTARR(ISHIFT)  = JARR
            IAS         = PTRARW(IARR)+IW4(IARR,1)
            IW4(IARR,1) = IW4(IARR,1) - 1
            DBLARR(IAS)      = VAL
            IF ( (KEEP(50) .NE. 0 .OR. KEEP(234).NE.0 )
     &           .AND.  IW4(IARR,1) .EQ. 0 .AND.
     &           STEP( IARR) > 0 ) THEN
              IF (MUMPS_275( PROCNODE_STEPS(abs(STEP(IARR))),
     &                            SLAVEF ) == MYID) THEN
                TAILLE = INTARR( PTRAIW(IARR) )
                CALL SMUMPS_310( N, PERM,
     &             INTARR( PTRAIW(IARR) + 3 ),
     &             DBLARR( PTRARW(IARR) + 1 ),
     &             TAILLE, 1, TAILLE )
              END IF
            END IF
          ENDIF
        END IF
        IF ( DEST.EQ. -1 ) THEN
         DO I=1, CANDIDATES(SLAVEF+1,ISTEP_TO_INIV2(ISTEP))
           DEST=CANDIDATES(I,ISTEP_TO_INIV2(ISTEP))
           IF (KEEP(46).EQ.0) DEST=DEST+1
           IF (DEST.NE.0)
     &     CALL SMUMPS_34( ISEND, JSEND, VAL,
     &     DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     &     LP, COMM, KEEP(46))
         ENDDO
         DEST = MASTER_NODE
         IF (KEEP(46).EQ.0) DEST=DEST+1
         IF ( DEST .NE. 0 ) THEN
           CALL SMUMPS_34( ISEND, JSEND, VAL,
     &     DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     &     LP, COMM, KEEP(46))
         ENDIF
         IF ((T4_MASTER_CONCERNED).AND.(T4MASTER.GT.0)) THEN 
          CALL SMUMPS_34( ISEND, JSEND, VAL,
     &     T4MASTER, BUFI, BUFR, NBRECORDS, NBUFS, 
     &     LP, COMM, KEEP(46))
         ENDIF 
        ELSE IF ( DEST .GT. 0 ) THEN
         CALL SMUMPS_34( ISEND, JSEND, VAL,
     &    DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     &    LP, COMM, KEEP(46))
         IF ( T4MASTER.GT.0 ) THEN
          CALL SMUMPS_34( ISEND, JSEND, VAL,
     &    T4MASTER, BUFI, BUFR, NBRECORDS, NBUFS, 
     &    LP, COMM, KEEP(46))
         ENDIF
        ELSE IF ( T4MASTER.GT.0 ) THEN
         CALL SMUMPS_34( ISEND, JSEND, VAL,
     &    T4MASTER, BUFI, BUFR, NBRECORDS, NBUFS, 
     &    LP, COMM, KEEP(46))
        END IF
  120 CONTINUE
      KEEP(49) = ARROW_ROOT
      IF (NBUFS.GT.0) THEN
       CALL SMUMPS_18(
     &   BUFI, BUFR, NBRECORDS, NBUFS,
     &   LP, COMM, KEEP( 46 ) )
      ENDIF
      IF ( KEEP( 46 ) .NE. 0 ) DEALLOCATE( IW4 )
      IF (NBUFS.GT.0) THEN
        DEALLOCATE( BUFI )
        DEALLOCATE( BUFR )
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_148
      SUBROUTINE SMUMPS_34(ISEND, JSEND, VAL,
     &   DEST, BUFI, BUFR, NBRECORDS, NBUFS, LP, COMM,
     &   TYPE_PARALL )
      IMPLICIT NONE
      INTEGER ISEND, JSEND, DEST, NBUFS, NBRECORDS, TYPE_PARALL
      INTEGER BUFI( NBRECORDS * 2 + 1, NBUFS )
      REAL BUFR( NBRECORDS, NBUFS )
      INTEGER COMM
      INTEGER LP
      REAL VAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER TAILLE_SENDI, TAILLE_SENDR, IREQ
         IF (BUFI(1,DEST)+1.GT.NBRECORDS) THEN
          TAILLE_SENDI = BUFI(1,DEST) * 2 + 1
          TAILLE_SENDR = BUFI(1,DEST)
          CALL MPI_SEND(BUFI(1,DEST),TAILLE_SENDI,
     &                   MPI_INTEGER,
     &                   DEST, ARROWHEAD, COMM, IERR )
          CALL MPI_SEND( BUFR(1,DEST), TAILLE_SENDR,
     &                   MPI_REAL, DEST,
     &                   ARROWHEAD, COMM, IERR )
          BUFI(1,DEST) = 0
         ENDIF
         IREQ = BUFI(1,DEST) + 1
         BUFI(1,DEST) = IREQ
         BUFI( IREQ * 2, DEST )     = ISEND
         BUFI( IREQ * 2 + 1, DEST ) = JSEND
         BUFR( IREQ, DEST )         = VAL
      RETURN
      END SUBROUTINE SMUMPS_34
      SUBROUTINE SMUMPS_18(
     &   BUFI, BUFR, NBRECORDS, NBUFS, LP, COMM,
     &   TYPE_PARALL )
      IMPLICIT NONE
      INTEGER NBUFS, NBRECORDS, TYPE_PARALL
      INTEGER BUFI( NBRECORDS * 2 + 1, NBUFS )
      REAL BUFR( NBRECORDS, NBUFS )
      INTEGER COMM
      INTEGER LP
      INTEGER ISLAVE, TAILLE_SENDI, TAILLE_SENDR, IERR
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        DO ISLAVE = 1,NBUFS 
          TAILLE_SENDI = BUFI(1,ISLAVE) * 2 + 1
          TAILLE_SENDR = BUFI(1,ISLAVE)
          BUFI(1,ISLAVE) = - BUFI(1,ISLAVE)
          CALL MPI_SEND(BUFI(1,ISLAVE),TAILLE_SENDI,
     &                   MPI_INTEGER,
     &                   ISLAVE, ARROWHEAD, COMM, IERR )
          IF ( TAILLE_SENDR .NE. 0 ) THEN
            CALL MPI_SEND( BUFR(1,ISLAVE), TAILLE_SENDR,
     &                     MPI_REAL, ISLAVE,
     &                     ARROWHEAD, COMM, IERR )
          END IF
        ENDDO
      RETURN
      END SUBROUTINE SMUMPS_18
      RECURSIVE SUBROUTINE SMUMPS_310( N, PERM, 
     &            INTLIST, DBLLIST, TAILLE, LO, HI )
      IMPLICIT NONE
      INTEGER N, TAILLE
      INTEGER PERM( N ) 
      INTEGER INTLIST( TAILLE )
      REAL DBLLIST( TAILLE )
      INTEGER LO, HI
      INTEGER I,J
      INTEGER ISWAP, PIVOT
      REAL sswap
      I = LO
      J = HI
      PIVOT = PERM(INTLIST((I+J)/2))
 10   IF (PERM(INTLIST(I)) < PIVOT) THEN
        I=I+1
        GOTO 10
      ENDIF
 20   IF (PERM(INTLIST(J)) > PIVOT) THEN
        J=J-1
        GOTO 20
      ENDIF
      IF (I < J) THEN
        ISWAP = INTLIST(I)
        INTLIST(I) = INTLIST(J)
        INTLIST(J)=ISWAP
        sswap = DBLLIST(I)
        DBLLIST(I) = DBLLIST(J)
        DBLLIST(J) = sswap
      ENDIF
      IF ( I <= J) THEN
        I = I+1
        J = J-1
      ENDIF
      IF ( I <= J ) GOTO 10
      IF ( LO < J ) CALL SMUMPS_310(N, PERM,
     &              INTLIST, DBLLIST, TAILLE, LO, J)
      IF ( I < HI ) CALL SMUMPS_310(N, PERM,
     &              INTLIST, DBLLIST, TAILLE, I, HI)
      RETURN
      END SUBROUTINE SMUMPS_310
      SUBROUTINE SMUMPS_145(  N,
     &    DBLARR, LDBLARR, INTARR, LINTARR, PTRAIW, PTRARW, 
     &    KEEP, KEEP8, MYID,  COMM, NBRECORDS,
     &    A, LA, root,
     &    PROCNODE_STEPS,
     &    SLAVEF, PERM, FRERE_STEPS, STEP, INFO1, INFO2
     &   )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INTEGER N, MYID, LDBLARR, LINTARR,
     &        COMM
      INTEGER INTARR(LINTARR) 
      INTEGER PTRAIW(N), PTRARW(N) 
      INTEGER   KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8), intent(IN) :: LA
      INTEGER PROCNODE_STEPS( KEEP(28) ), PERM( N )
      INTEGER SLAVEF, NBRECORDS
      REAL A( LA )
      INTEGER INFO1, INFO2
      REAL DBLARR(LDBLARR)
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER, POINTER, DIMENSION(:) :: BUFI
      REAL, POINTER, DIMENSION(:) :: BUFR
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      LOGICAL FINI 
      INTEGER IREC, NB_REC, IARR, JARR, IA, I1, I, allocok
      INTEGER IS, IS1, ISHIFT, IIW, IAS
      INTEGER LOCAL_M, LOCAL_N, ILOCROOT, JLOCROOT, 
     &        IPOSROOT, JPOSROOT, TAILLE,
     &        IPROC
      INTEGER FRERE_STEPS( KEEP(28) ), STEP(N)
      INTEGER(8) :: PTR_ROOT
      INTEGER ARROW_ROOT, TYPE_PARALL
      INTEGER MUMPS_330, MUMPS_275
      EXTERNAL MUMPS_330, MUMPS_275
      REAL VAL
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MASTER
      PARAMETER(MASTER=0)
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR
      INTEGER numroc
      EXTERNAL numroc
      TYPE_PARALL = KEEP(46)
      ARROW_ROOT=0
      ALLOCATE( BUFI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = NBRECORDS * 2 + 1
        WRITE(*,*) MYID,': Could not allocate BUFI: goto 500'
        GOTO 500
      END IF
      ALLOCATE( BUFR( NBRECORDS )        , stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = NBRECORDS
        WRITE(*,*) MYID,': Could not allocate BUFR: goto 500'
        GOTO 500
      END IF
      ALLOCATE( IW4(N,2), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = 2 * N
        WRITE(*,*) MYID,': Could not allocate IW4: goto 500'
        GOTO 500
      END IF
      IF ( KEEP(38).NE.0) THEN
        IF (KEEP(60)==0) THEN
         LOCAL_M = numroc( root%ROOT_SIZE, root%MBLOCK,
     &             root%MYROW, 0, root%NPROW )
         LOCAL_M = max( 1, LOCAL_M )
         LOCAL_N = numroc( root%ROOT_SIZE, root%NBLOCK,
     &             root%MYCOL, 0, root%NPCOL )
         PTR_ROOT = LA - int(LOCAL_M,8) * int(LOCAL_N,8) + 1_8
         IF ( PTR_ROOT .LE. LA ) THEN
           A( PTR_ROOT:LA ) = ZERO
         END IF
        ELSE
         DO I=1, root%SCHUR_NLOC
           root%SCHUR_POINTER((I-1)*root%SCHUR_LLD+1:
     &     (I-1)*root%SCHUR_LLD+root%SCHUR_MLOC)=ZERO
         ENDDO
        ENDIF
      END IF
      FINI = .FALSE.
      DO I=1,N
       I1 = PTRAIW(I)
       IA = PTRARW(I)
       IF (IA.GT.0) THEN
        DBLARR(IA) = ZERO
        IW4(I,1) = INTARR(I1)
        IW4(I,2) = -INTARR(I1+1)
        INTARR(I1+2)=I
       ENDIF
      ENDDO
      DO WHILE (.NOT.FINI) 
       CALL MPI_RECV( BUFI(1), 2*NBRECORDS+1, 
     &                MPI_INTEGER, MASTER, 
     &                ARROWHEAD,
     &                COMM, STATUS, IERR )
       NB_REC = BUFI(1)
       IF (NB_REC.LE.0) THEN
         FINI = .TRUE.
         NB_REC = -NB_REC 
       ENDIF
       IF (NB_REC.EQ.0) EXIT
       CALL MPI_RECV( BUFR(1), NBRECORDS, MPI_REAL,
     &                  MASTER, ARROWHEAD,
     &                COMM, STATUS, IERR )
       DO IREC=1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        IF ( MUMPS_330( PROCNODE_STEPS(abs(STEP(abs(IARR)))),
     &                       SLAVEF ) .eq. 3 ) THEN
          ARROW_ROOT = ARROW_ROOT + 1
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
          END IF
            ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
            JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
           IF (KEEP(60)==0) THEN
             A( PTR_ROOT + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &                   + int(ILOCROOT - 1,8) )
     &       =  A( PTR_ROOT + int(JLOCROOT - 1,8)
     &                      * int(LOCAL_M,8)
     &                      + int(ILOCROOT - 1,8))
     &        + VAL
           ELSE
             root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8) )
     &       = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8))
     &       + VAL
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
           IF ( (KEEP(50) .NE. 0 .OR. KEEP(234).NE.0)
     &          .AND.  IW4(IARR,1) .EQ. 0 
     &          .AND. STEP(IARR) > 0 ) THEN
              IPROC = MUMPS_275( PROCNODE_STEPS(abs(STEP(IARR))),
     &                                SLAVEF )
              IF ( TYPE_PARALL .eq. 0 ) THEN
                IPROC = IPROC + 1
              END IF 
              IF (IPROC .EQ. MYID) THEN
                TAILLE = INTARR( PTRAIW(IARR) )
                CALL SMUMPS_310( N, PERM,
     &            INTARR( PTRAIW(IARR) + 3 ),
     &            DBLARR( PTRARW(IARR) + 1 ),
     &            TAILLE, 1, TAILLE )
              END IF
           END IF
        ENDIF
       ENDDO
      END DO
      DEALLOCATE( BUFI )
      DEALLOCATE( BUFR )
      DEALLOCATE( IW4 )
 500  CONTINUE
      KEEP(49) = ARROW_ROOT
      RETURN 
      END SUBROUTINE SMUMPS_145
      SUBROUTINE SMUMPS_266( MYID, BUFR, LBUFR, 
     &     LBUFR_BYTES,
     &     IWPOS, IWPOSCB,
     &     IPTRLU, LRLU, LRLUS,
     &     TNBPROCFILS, N, IW, LIW, A, LA,
     &     PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, COMP,
     &     KEEP,KEEP8, ITLOC, RHS_MUMPS,
     &     IFLAG, IERROR )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INTEGER MYID
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA
      INTEGER IWPOS, IWPOSCB, N, LIW
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER PTRIST(KEEP(28)), STEP(N), 
     & PIMASTER(KEEP(28)), 
     & TNBPROCFILS( KEEP(28) ), ITLOC( N + KEEP(253) )
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER COMP, IFLAG, IERROR
      INTEGER INODE, NBPROCFILS, NCOL, NROW, NASS, NSLAVES
      INTEGER NSLAVES_RECU, NFRONT
      INTEGER LREQ
      INTEGER(8) :: LREQCB
      DOUBLE PRECISION FLOP1
      INCLUDE 'mumps_headers.h'
      INODE = BUFR( 1 )
      NBPROCFILS = BUFR( 2 )
      NROW = BUFR( 3 )
      NCOL = BUFR( 4 )
      NASS = BUFR( 5 )
      NFRONT = BUFR( 6 )
      NSLAVES_RECU = BUFR( 7 )
      IF ( KEEP(50) .eq. 0 ) THEN
         FLOP1 = dble( NASS * NROW ) +
     &     dble(NROW*NASS)*dble(2*NCOL-NASS-1)
      ELSE
         FLOP1 = dble( NASS ) * dble( NROW )
     &            * dble( 2 * NCOL - NROW - NASS + 1)
      END IF
      CALL SMUMPS_190(1,.TRUE.,FLOP1, KEEP,KEEP8)
      IF ( KEEP(50) .eq. 0 ) THEN
        NSLAVES = NSLAVES_RECU + XTRA_SLAVES_UNSYM
      ELSE
        NSLAVES = NSLAVES_RECU + XTRA_SLAVES_SYM   
      END IF
      LREQ   = NROW + NCOL + 6 + NSLAVES + KEEP(IXSZ)
      LREQCB = int(NCOL,8) * int(NROW,8)
      CALL SMUMPS_22(.FALSE., 0_8, .FALSE.,.TRUE.,
     &   MYID,N, KEEP,KEEP8, IW, LIW, A, LA,
     &   LRLU, IPTRLU,IWPOS,IWPOSCB,
     &   PTRIST,PTRAST, STEP, PIMASTER,PAMASTER,
     &   LREQ, LREQCB, INODE, S_ACTIVE, .TRUE.,
     &   COMP, LRLUS, IFLAG, IERROR
     &     )
      IF ( IFLAG .LT. 0 ) RETURN
      PTRIST(STEP(INODE)) = IWPOSCB + 1
      PTRAST(STEP(INODE)) = IPTRLU  + 1_8
      IW( IWPOSCB + 1+KEEP(IXSZ) ) = NCOL
      IW( IWPOSCB + 2+KEEP(IXSZ) ) = - NASS
      IW( IWPOSCB + 3+KEEP(IXSZ) ) = NROW
      IW( IWPOSCB + 4+KEEP(IXSZ) ) = 0
      IW( IWPOSCB + 5+KEEP(IXSZ) ) = NASS
      IW( IWPOSCB + 6+KEEP(IXSZ) ) = NSLAVES
      IW( IWPOSCB + 7+KEEP(IXSZ)+NSLAVES : 
     &           IWPOSCB + 6+KEEP(IXSZ)+NSLAVES + NROW + NCOL )
     &= BUFR( 8 + NSLAVES_RECU : 7 + NSLAVES_RECU + NROW + NCOL )
      IF ( KEEP(50) .eq. 0 ) THEN
        IW( IWPOSCB + 7+KEEP(IXSZ) ) = S_ROOTBAND_INIT
        IF (NSLAVES_RECU.GT.0) 
     &  IW( IWPOSCB + 7+XTRA_SLAVES_UNSYM+KEEP(IXSZ): 
     &       IWPOSCB+6+KEEP(IXSZ)+NSLAVES_RECU ) =
     &       BUFR( 8: 7 + NSLAVES_RECU )
      ELSE
        IW( IWPOSCB + 7+KEEP(IXSZ) ) = 0
        IW( IWPOSCB + 8+KEEP(IXSZ) ) = NFRONT
        IW( IWPOSCB + 9+KEEP(IXSZ) ) = S_ROOTBAND_INIT
        IW( IWPOSCB + 7+XTRA_SLAVES_SYM+KEEP(IXSZ):
     &      IWPOSCB + 6+XTRA_SLAVES_SYM+KEEP(IXSZ)+NSLAVES_RECU ) =
     &       BUFR( 8: 7 + NSLAVES_RECU )
      END IF
      TNBPROCFILS(STEP( INODE )) = NBPROCFILS
      RETURN
      END SUBROUTINE SMUMPS_266
      SUBROUTINE SMUMPS_163( id )
      USE SMUMPS_STRUC_DEF
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (SMUMPS_STRUC) id
      INTEGER MASTER, IERR,PAR_loc,SYM_loc
      PARAMETER( MASTER = 0 )
      INTEGER color
      CALL MPI_COMM_SIZE(id%COMM, id%NPROCS, IERR )
      PAR_loc=id%PAR
      SYM_loc=id%SYM
      CALL MPI_BCAST(PAR_loc,1,MPI_INTEGER,MASTER,id%COMM,IERR)
      CALL MPI_BCAST(SYM_loc,1,MPI_INTEGER,MASTER,id%COMM,IERR)
      IF ( PAR_loc .eq. 0 ) THEN
        IF ( id%MYID .eq. MASTER ) THEN
          color = MPI_UNDEFINED
        ELSE
          color = 0
        END IF
        CALL MPI_COMM_SPLIT( id%COMM, color, 0,
     &                       id%COMM_NODES, IERR )
        id%NSLAVES = id%NPROCS - 1
      ELSE
        CALL MPI_COMM_DUP( id%COMM, id%COMM_NODES, IERR )
        id%NSLAVES = id%NPROCS
      END IF
      IF (PAR_loc .ne. 0 .or. id%MYID .NE. MASTER) THEN
        CALL MPI_COMM_DUP( id%COMM_NODES, id%COMM_LOAD, IERR )
      ENDIF
      CALL SMUMPS_20( id%NSLAVES, id%LWK_USER,
     &    id%CNTL(1), id%ICNTL(1),
     &    id%KEEP(1), id%KEEP8(1), id%INFO(1), id%INFOG(1),
     &    id%RINFO(1), id%RINFOG(1),
     &    SYM_loc, PAR_loc, id%DKEEP(1) )
      id%WRITE_PROBLEM="NAME_NOT_INITIALIZED"
      CALL MUMPS_SET_VERSION( id%VERSION_NUMBER )
      id%OOC_TMPDIR="NAME_NOT_INITIALIZED"
      id%OOC_PREFIX="NAME_NOT_INITIALIZED"
      id%NRHS = 1
      id%LRHS = 0 
      id%LREDRHS = 0
      CALL SMUMPS_61( id%KEEP( 34 ), id%KEEP(35) )
      NULLIFY(id%BUFR)
      id%MAXIS1 = 0
      id%INST_Number = -1
      id%N = 0; id%NZ = 0
      NULLIFY(id%IRN)
      NULLIFY(id%JCN)
      NULLIFY(id%A)
      id%NZ_loc = 0 
      NULLIFY(id%IRN_loc)
      NULLIFY(id%JCN_loc)
      NULLIFY(id%A_loc)
      NULLIFY(id%MAPPING)
      NULLIFY(id%RHS)
      NULLIFY(id%REDRHS)
      id%NZ_RHS=0
      NULLIFY(id%RHS_SPARSE)
      NULLIFY(id%IRHS_SPARSE)
      NULLIFY(id%IRHS_PTR)
      NULLIFY(id%ISOL_loc)
      id%LSOL_loc=0
      NULLIFY(id%SOL_loc)
      NULLIFY(id%COLSCA)
      NULLIFY(id%ROWSCA)
      NULLIFY(id%PERM_IN)
      NULLIFY(id%IS)
      NULLIFY(id%IS1)
      NULLIFY(id%STEP)
      NULLIFY(id%Step2node)
      NULLIFY(id%DAD_STEPS)
      NULLIFY(id%NE_STEPS)
      NULLIFY(id%ND_STEPS)
      NULLIFY(id%FRERE_STEPS)
      NULLIFY(id%SYM_PERM)
      NULLIFY(id%UNS_PERM)
      NULLIFY(id%PIVNUL_LIST)
      NULLIFY(id%FILS)
      NULLIFY(id%PTRAR)
      NULLIFY(id%FRTPTR)
      NULLIFY(id%FRTELT)
      NULLIFY(id%NA)
      id%LNA=0
      NULLIFY(id%PROCNODE_STEPS)
      NULLIFY(id%S)
      NULLIFY(id%PROCNODE)
      NULLIFY(id%POIDS)
      NULLIFY(id%PTLUST_S)
      NULLIFY(id%PTRFAC)
      NULLIFY(id%INTARR) 
      NULLIFY(id%DBLARR)
      NULLIFY(id%DEPTH_FIRST)
      NULLIFY(id%DEPTH_FIRST_SEQ)
      NULLIFY(id%SBTR_ID)
      NULLIFY(id%MEM_SUBTREE)
      NULLIFY(id%MEM_SUBTREE)
      NULLIFY(id%MY_ROOT_SBTR)
      NULLIFY(id%MY_FIRST_LEAF)
      NULLIFY(id%MY_NB_LEAF)
      NULLIFY(id%COST_TRAV)
      NULLIFY(id%RHSCOMP)
      NULLIFY(id%POSINRHSCOMP)
      NULLIFY(id%OOC_INODE_SEQUENCE)
      NULLIFY(id%OOC_TOTAL_NB_NODES)
      NULLIFY(id%OOC_SIZE_OF_BLOCK)
      NULLIFY(id%OOC_FILE_NAME_LENGTH)
      NULLIFY(id%OOC_FILE_NAMES)
      NULLIFY(id%OOC_VADDR)
      NULLIFY(id%OOC_NB_FILES)
      NULLIFY(id%CB_SON_SIZE)
      NULLIFY(id%root%RHS_CNTR_MASTER_ROOT)
      NULLIFY(id%root%RHS_ROOT)
      NULLIFY(id%root%RG2L_ROW)
      NULLIFY(id%root%RG2L_COL)
      NULLIFY(id%root%IPIV)
      NULLIFY(id%root%SCHUR_POINTER)
      NULLIFY(id%SCHUR_CINTERFACE)
      id%NELT=0
      NULLIFY(id%ELTPTR)
      NULLIFY(id%ELTVAR)
      NULLIFY(id%A_ELT)
      NULLIFY(id%ELTPROC)
      id%SIZE_SCHUR = 0
      NULLIFY( id%LISTVAR_SCHUR )
      NULLIFY( id%SCHUR )
      id%NPROW      = 0
      id%NPCOL      = 0
      id%MBLOCK     = 0
      id%NBLOCK     = 0
      id%SCHUR_MLOC = 0 
      id%SCHUR_NLOC = 0 
      id%SCHUR_LLD  = 0
      NULLIFY(id%ISTEP_TO_INIV2)
      NULLIFY(id%I_AM_CAND)
      NULLIFY(id%FUTURE_NIV2)
      NULLIFY(id%TAB_POS_IN_PERE)
      NULLIFY(id%CANDIDATES)
      CALL SMUMPS_637(id)
      NULLIFY(id%MEM_DIST)
      NULLIFY(id%SUP_PROC)
      id%Deficiency = 0
      id%root%LPIV = -1
      id%root%yes  = .FALSE.
      id%root%gridinit_done  = .FALSE.
        IF ( id%KEEP( 46 ) .ne. 0  .OR.
     &     id%MYID .ne. MASTER ) THEN
          CALL MPI_COMM_RANK
     &         (id%COMM_NODES, id%MYID_NODES, IERR )
        ELSE
          id%MYID_NODES = -464646
        ENDIF
      RETURN
      END SUBROUTINE SMUMPS_163
      SUBROUTINE SMUMPS_252( COMM_LOAD, ASS_IRECV,
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
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, JOBASS, ETATASS 
     &    )
      USE SMUMPS_COMM_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER IZERO 
      PARAMETER (IZERO=0)
      INTEGER N,LIW,NSTEPS
      INTEGER(8) LA, LRLU, LRLUS, IPTRLU, POSFAC
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     &        IWPOS, IWPOSCB, COMP
      INTEGER JOBASS,ETATASS 
      LOGICAL SON_LEVEL2
      REAL A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER NBPROCFILS(KEEP(28)) 
      INTEGER NSTK_S(KEEP(28)),PROCNODE_STEPS(KEEP(28))
      INTEGER IPOOL( LPOOL )
      INTEGER BUFR( LBUFR )
      INTEGER IDUMMY(1)
      INTEGER IW(LIW), ITLOC(N+KEEP(253)),
     &        PTRARW(N), PTRAIW(N), ND(KEEP(28)), PERM(N), 
     &        FILS(N), FRERE(KEEP(28)), DAD(KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28)),
     &              PAMASTER(KEEP(28))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER      INTARR(max(1,KEEP(14)))
      REAL DBLARR(max(1,KEEP(13)))
      INTEGER  MUMPS_330
      EXTERNAL MUMPS_330
      INTEGER LP, HS, HF
      INTEGER NBPANELS_L, NBPANELS_U
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFS4FATHER
      INTEGER(8) NFRONT8, LAELL8, LAELL_REQ8
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER LREQ_OOC
      INTEGER(8) :: SIZFR
      INTEGER SIZFI, NCB
      INTEGER J1,J2
      INTEGER NCOLS, NROWS, LDA_SON
      INTEGER(8) :: JJ2, ICT13
#if defined(ALLOW_NON_INIT)
      INTEGER(8) :: NUMROWS, JJ3, JJ8, APOS_ini
#endif
      INTEGER NELIM,JJ,JJ1,J3,
     &        IBROT,IORG
      INTEGER JPOS,ICT11
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4
      INTEGER(8) IACHK, POSELT, LAPOS2, IACHK_ini
      INTEGER(8) APOS, APOS2, APOS3, POSEL1, ICT12
      INTEGER AINPUT
      INTEGER NSLAVES, NSLSON, NPIVS,NPIV_ANA,NPIV
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      INTEGER ISON_IN_PLACE 
      INTEGER ISON_TOP 
      INTEGER(8) SIZE_ISON_TOP8
      LOGICAL RESET_TO_ZERO, RISK_OF_SAME_POS,
     &        RISK_OF_SAME_POS_THIS_LINE, OMP_PARALLEL_FLAG
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE
      INTEGER INDX, FIRST_INDEX, SHIFT_INDEX
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INCLUDE 'mumps_headers.h'
      INTEGER NCBSON
      LOGICAL SAME_PROC
      INTRINSIC real
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER NELT, LPTRAR
      EXTERNAL MUMPS_167
      LOGICAL MUMPS_167
      LOGICAL SSARBR
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      DOUBLE PRECISION FLOP1,FLOP1_EFF
      EXTERNAL MUMPS_170
      LOGICAL MUMPS_170
      COMPRESSCB =.FALSE.
      NELT       = 1
      LPTRAR     = N
      NFS4FATHER = -1
      IN         = INODE
      NBPROCFILS(STEP(IN)) = 0
      LEVEL = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF (LEVEL.NE.1) THEN 
       write(6,*) 'Error1 in mpi51f_niv1 '
       CALL MUMPS_ABORT()
      ENDIF
      NSLAVES = 0
      HF =  6 + NSLAVES + KEEP(IXSZ)
      IF (JOBASS.EQ.0) THEN
        ETATASS= 0 
      ELSE
        ETATASS= 2 
        IOLDPS = PTLUST_S(STEP(INODE)) 
        NFRONT = IW(IOLDPS + KEEP(IXSZ)) 
        NASS1  = iabs(IW(IOLDPS + 2 + KEEP(IXSZ)))
        ICT11 = IOLDPS + HF - 1 + NFRONT
        SSARBR=MUMPS_167(PROCNODE_STEPS(STEP(INODE)),
     &                        SLAVEF)
        NUMORG = 0
        DO WHILE (IN.GT.0)
          NUMORG = NUMORG + 1
          IN = FILS(IN)
        ENDDO
        NUMSTK = 0
        IFSON = -IN
        ISON = IFSON
        IF (ISON .NE. 0) THEN
         DO WHILE (ISON .GT. 0)
           NUMSTK = NUMSTK + 1
           ISON = FRERE(STEP(ISON))
         ENDDO
        ENDIF
        GOTO 123
      ENDIF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      ENDDO
      NPIV_ANA=NUMORG
      NSTEPS = NSTEPS + 1
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .NE. 0) THEN
        DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         NASS = NASS + IW(PIMASTER(STEP(ISON))+1+KEEP(IXSZ))
         ISON = FRERE(STEP(ISON))
         ENDDO
      ENDIF
      NFRONT = ND(STEP(INODE)) + NASS + KEEP(253) 
      NASS1 = NASS + NUMORG
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1) THEN 
        CALL SMUMPS_684( KEEP(50), NFRONT, NFRONT, NASS1,
     &                                NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL SMUMPS_94(N, KEEP(28),
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
          ENDIF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      ENDIF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      ISON_TOP      = -9999
      ISON_IN_PLACE = -9999
      SIZE_ISON_TOP8 = 0_8
      IF (KEEP(234).NE.0) THEN
        IF ( IWPOSCB .NE. LIW ) THEN 
        IF ( IWPOSCB+IW(IWPOSCB+1+XXI).NE.LIW) THEN
          ISON = IW( IWPOSCB + 1 + XXN )
          IF ( DAD( STEP( ISON ) ) .EQ. INODE .AND.
     &    MUMPS_330(PROCNODE_STEPS(STEP(ISON)),SLAVEF)
     &    .EQ. 1 )
     &    THEN
            ISON_TOP = ISON
            CALL MUMPS_729(SIZE_ISON_TOP8,IW(IWPOSCB + 1 + XXR))
            IF (LRLU .LT. int(NFRONT,8) * int(NFRONT,8)) THEN
              ISON_IN_PLACE = ISON
            ENDIF
          END IF
        END IF
        END IF
      END IF
      NIV1 = .TRUE.
      IF (KEEP(50).EQ.0 .AND. KEEP(234) .EQ. 0) THEN
        CALL  MUMPS_81(MYID, INODE, N, IOLDPS, HF, NFRONT, 
     &        NFRONT_EFF, DAD,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     &        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     &        INTARR, ITLOC, RHS_MUMPS, FILS, FRERE,
     &        SON_LEVEL2, NIV1, NBPROCFILS, KEEP, KEEP8, IFLAG,
     &        PROCNODE_STEPS, SLAVEF )
      ELSE
        CALL MUMPS_86( MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM, DAD,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     &        INTARR, ITLOC, RHS_MUMPS, FILS, FRERE,
     &        SON_LEVEL2, NIV1, NBPROCFILS, KEEP, KEEP8, IFLAG,
     &        ISON_IN_PLACE, 
     &        PROCNODE_STEPS, SLAVEF)
        IF (IFLAG.LT.0) GOTO 300
      ENDIF
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
              CALL SMUMPS_190(0,.FALSE.,FLOP1-FLOP1_EFF,
     &             KEEP,KEEP8)
           ENDIF
           IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF))
           NFRONT = NFRONT_EFF
           LREQ = HF + 2 * NFRONT + LREQ_OOC
        ELSE
           Write(*,*) ' ERROR 1 during ass_niv1', NFRONT, NFRONT_EFF
           GOTO 270
        ENDIF
      ENDIF
      NFRONT8=int(NFRONT,8)
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
        CALL SMUMPS_691(KEEP(50),
     &       NBPANELS_L, NBPANELS_U, NASS1, 
     &       IOLDPS + HF + 2 * NFRONT, IW, LIW)
      ENDIF
      NCB   = NFRONT - NASS1
      MAXFRW = max0(MAXFRW, NFRONT)
      ICT11 = IOLDPS + HF - 1 + NFRONT 
      LAELL8 = NFRONT8 * NFRONT8
      LAELL_REQ8 = LAELL8
      IF ( ISON_IN_PLACE > 0 ) THEN
        LAELL_REQ8 = LAELL8 - SIZE_ISON_TOP8
      ENDIF
      IF (LRLU .LT. LAELL_REQ8) THEN
        IF (LRLUS .LT. LAELL_REQ8) THEN
          GOTO 280
        ELSE
          CALL SMUMPS_94
     &        (N, KEEP(28), IW, LIW, A, LA, LRLU, IPTRLU,
     &         IWPOS, IWPOSCB, PTRIST, PTRAST, STEP, PIMASTER,
     &         PAMASTER,KEEP(216),LRLUS,KEEP(IXSZ))
          COMP = COMP + 1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          ENDIF
        ENDIF
      ENDIF
      LRLU = LRLU - LAELL8 
      LRLUS = LRLUS - LAELL8 + SIZE_ISON_TOP8
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL8
      SSARBR=MUMPS_167(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      CALL SMUMPS_471(SSARBR,.FALSE.,
     &     LA-LRLUS, 
     &     0_8,
     &     LAELL8-SIZE_ISON_TOP8, 
     &     KEEP,KEEP8,
     &     LRLU)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = min(POSELT + LAELL8 - 1_8, IPTRLU)
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = min(POSELT + LAELL8 - 1_8, IPTRLU)
        DO JJ8 = POSELT, LAPOS2
           A( JJ8 ) = ZERO
        ENDDO
      ELSE
        IF (ETATASS.EQ.1) THEN
         APOS_ini = POSELT
         DO JJ8 = 0_8, NFRONT8 - 1_8
          JJ3 = min(JJ8,int(NASS1-1,8)) 
          APOS = APOS_ini + JJ8 * NFRONT8
          A(APOS:APOS+JJ3) = ZERO
         END DO
        ELSE
          APOS_ini = POSELT
          NUMROWS = min(NFRONT8, (IPTRLU-APOS_ini) / NFRONT8 )
          DO JJ8 = 0_8, NUMROWS - 1_8
             APOS = APOS_ini + JJ8 * NFRONT8
             A(APOS:APOS + JJ8) = ZERO
          ENDDO
          IF( NUMROWS .LT. NFRONT8 ) THEN
            APOS = APOS_ini + NFRONT8*NUMROWS
            A(APOS : min(IPTRLU,APOS+NUMROWS)) = ZERO
          ENDIF
        ENDIF
      END IF
#endif
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      PTLUST_S(STEP(INODE)) = IOLDPS
      IW(IOLDPS+XXI)   = LREQ  
      CALL MUMPS_730(LAELL8,IW(IOLDPS+XXR)) 
      IW(IOLDPS+XXS) =-9999
      IW(IOLDPS+XXS+1:IOLDPS+KEEP(IXSZ)-1)=-99999
      IW(IOLDPS + KEEP(IXSZ))   = NFRONT
      IW(IOLDPS + KEEP(IXSZ) + 1) = 0
      IW(IOLDPS + KEEP(IXSZ) + 2) = -NASS1
      IW(IOLDPS + KEEP(IXSZ) + 3) = -NASS1
      IW(IOLDPS + KEEP(IXSZ) + 4) = STEP(INODE)
      IW(IOLDPS + KEEP(IXSZ) + 5)   = NSLAVES
 123  CONTINUE  
      IF (NUMSTK.NE.0) THEN
        IF (ISON_TOP > 0) THEN
          ISON = ISON_TOP
        ELSE
          ISON = IFSON
        ENDIF
        DO 220 IELL = 1, NUMSTK
          ISTCHK    = PIMASTER(STEP(ISON))
          LSTK      = IW(ISTCHK + KEEP(IXSZ))
          NELIM     = IW(ISTCHK + KEEP(IXSZ) + 1)
          NPIVS     = IW(ISTCHK + KEEP(IXSZ) + 3)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON    = IW(ISTCHK + KEEP(IXSZ) + 5)
          HS        = 6 + KEEP(IXSZ) + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          IF ( SAME_PROC ) THEN
            COMPRESSCB=( IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
            COMPRESSCB=( IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          LEVEL1    = NSLSON.EQ.0
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW( ISTCHK + KEEP(IXSZ) + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          SIZFI   = HS + NROWS + NCOLS 
          J1 = ISTCHK + HS + NROWS + NPIVS
          IF ( .NOT. LEVEL1 .AND. NELIM.EQ.0 ) GOTO 205
          IF (LEVEL1) THEN
           J2 = J1 + LSTK - 1
           SIZFR  = int(LSTK,8)*int(LSTK,8)
           IF (COMPRESSCB) SIZFR = (int(LSTK,8)*int(LSTK+1,8))/2_8
          ELSE
           IF ( KEEP(50).eq.0 ) THEN
             SIZFR = int(NELIM,8) * int(LSTK,8)
           ELSE
             SIZFR = int(NELIM,8) * int(NELIM,8)
           END IF
           J2 = J1 + NELIM - 1
          ENDIF
          IF (JOBASS.EQ.0) OPASSW = OPASSW + dble(SIZFR)
          IACHK = PAMASTER(STEP(ISON))
          IF ( KEEP(50) .eq. 0 ) THEN
            POSEL1 = PTRAST(STEP(INODE)) - NFRONT8
            IF (NFRONT .EQ. LSTK.AND. ISON.EQ.ISON_IN_PLACE
     &          .AND.IACHK + SIZFR - 1_8 .EQ. POSFAC - 1_8 ) THEN
               GOTO 205
            ENDIF
            IF (J2.GE.J1) THEN
              RESET_TO_ZERO = (IACHK .LT. POSFAC) 
              RISK_OF_SAME_POS = IACHK + SIZFR - 1_8 .EQ. POSFAC - 1_8
              RISK_OF_SAME_POS_THIS_LINE = .FALSE.
              IACHK_ini = IACHK
              OMP_PARALLEL_FLAG = (RESET_TO_ZERO.EQV..FALSE.).AND.
     &            ((J2-J1).GT.300)
              DO 170 JJ = J1, J2
                APOS = POSEL1 + int(IW(JJ),8) * int(NFRONT,8)
                IACHK = IACHK_ini + int(JJ-J1,8)*int(LSTK,8)
                IF (RISK_OF_SAME_POS) THEN
                  IF (JJ.EQ.J2) THEN
                    RISK_OF_SAME_POS_THIS_LINE =
     &                  (ISON .EQ. ISON_IN_PLACE)
     &                  .AND. ( APOS + int(IW(J1+LSTK-1)-1,8).EQ.
     &                          IACHK+int(LSTK-1,8) )
                  ENDIF
                ENDIF
                IF ((IACHK .GE. POSFAC).AND.(JJ>J1))THEN
                   RESET_TO_ZERO =.FALSE.
                ENDIF
                IF (RESET_TO_ZERO) THEN
                  IF (RISK_OF_SAME_POS_THIS_LINE) THEN
                    DO JJ1 = 1, LSTK
                      JJ2 = APOS + int(IW(J1 + JJ1 - 1) - 1,8)
                      IF ( IACHK+int(JJ1-1,8) .NE. JJ2 ) THEN
                        A(JJ2) = A(IACHK + int(JJ1 - 1,8))
                        A(IACHK + int(JJ1 -1,8)) = ZERO
                      ENDIF
                    ENDDO
                  ELSE
                    DO JJ1 = 1, LSTK
                      JJ2 = APOS + int(IW(J1+JJ1-1),8) - 1_8
                      A(JJ2) = A(IACHK + int(JJ1 - 1,8))
                      A(IACHK + int(JJ1 -1,8)) = ZERO
                    ENDDO
                  ENDIF
                ELSE 
                  DO JJ1 = 1, LSTK
                    JJ2 = APOS + int(IW(J1+JJ1-1),8) - 1_8
                    A(JJ2) = A(JJ2) + A(IACHK + int(JJ1 - 1,8))
                  ENDDO
                ENDIF
  170         CONTINUE
            END IF
          ELSE
            IF (LEVEL1) THEN
             LDA_SON = LSTK  
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = SIZFR
            ELSE
              LCB = int(LDA_SON,8)* int(J2-J1+1,8)
            ENDIF
            CALL SMUMPS_178(A, LA,
     &           PTRAST(STEP( INODE )), NFRONT, NASS1,
     &           IACHK, LDA_SON, LCB,
     &           IW( J1 ), J2 - J1 + 1, NELIM, ETATASS, 
     &           COMPRESSCB, (ISON.EQ.ISON_IN_PLACE)
     &          )
          ENDIF
  205     IF (LEVEL1) THEN 
           IF (SAME_PROC) ISTCHK = PTRIST(STEP(ISON))
           IF ((SAME_PROC).AND.ETATASS.NE.1) THEN
             IF (KEEP(50).NE.0) THEN
              J2 = J1 + LSTK - 1
              DO JJ = J1, J2
               IW(JJ) = IW(JJ - NROWS)
              ENDDO
             ELSE
              J2 = J1 + LSTK - 1
              J3 = J1 + NELIM
              DO JJ = J3, J2
               IW(JJ) = IW(JJ - NROWS)
              ENDDO
              IF (NELIM .NE. 0) THEN
                J3 = J3 - 1
                DO JJ = J1, J3
                 JPOS = IW(JJ) + ICT11
                 IW(JJ) = IW(JPOS)
                ENDDO
              ENDIF
             ENDIF
           ENDIF
           IF (ETATASS.NE.1) THEN
             IF ( SAME_PROC ) THEN 
               PTRIST(STEP(ISON))   = -99999999
             ELSE
               PIMASTER(STEP( ISON )) = -99999999
             ENDIF
             CALL SMUMPS_152(SSARBR, MYID, N, ISTCHK,
     &          PAMASTER(STEP(ISON)),
     &          IW, LIW, LRLU, LRLUS, IPTRLU,
     &          IWPOSCB, LA, KEEP,KEEP8,
     &          (ISON .EQ. ISON_TOP)
     &          )
           ENDIF
          ELSE
           PDEST = ISTCHK + 6 + KEEP(IXSZ)
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL MUMPS_49( 
     &                KEEP, KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
              CALL SMUMPS_210( COMM_LOAD, ASS_IRECV, 
     &             BUFR, LBUFR, LBUFR_BYTES,
     &             INODE, ISON, NSLAVES, IDUMMY,
     &             NFRONT, NASS1,NFS4FATHER,
     &             TROW_SIZE, IW( INDX ),
     &         PROCNODE_STEPS,
     &         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &         LRLUS, N, IW,
     &         LIW, A, LA,
     &         PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &         PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL,
     &         LEAF, NBFIN, ICNTL, KEEP, KEEP8,  root,
     &         OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ND, FRERE,
     &         LPTRAR, NELT, IW, IW, 
     &
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 500
              EXIT
             ENDIF
           ENDDO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  SMUMPS_71( 
     &           INODE, NFRONT, NASS1, NFS4FATHER, 
     &           ISON, MYID,
     &       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     &       COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     &       KEEP, KEEP8, STEP, N, 
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &        )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL SMUMPS_329( 
     &         COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS,
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA,
     &         PTRIST, PTLUST_S, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, COMM,
     &         NBPROCFILS,
     &         IPOOL, LPOOL, LEAF,
     &         NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ICNTL, KEEP, KEEP8, ND, FRERE,
     &         LPTRAR, NELT, IW, IW,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
               IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           ENDDO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
           ENDIF
          ENDIF
        ISON = FRERE(STEP(ISON))
        IF (ISON .LE. 0) THEN
          ISON = IFSON
        ENDIF
  220 CONTINUE
      END IF
      IF (ETATASS.EQ.2) GOTO 500
      POSELT = PTRAST(STEP(INODE))
      IBROT = INODE
      DO 260 IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        IJROW = INTARR(J1)
        ICT12 = POSELT + int(IJROW - NFRONT - 1,8)
Cduplicates --> CVD$ DEPCHK
        DO 240 JJ = J1, J2
           APOS2 = ICT12 + int(INTARR(JJ),8) * NFRONT8
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
          AINPUT = AINPUT + 1
  240   CONTINUE
        IF (J3 .LE. J4) THEN
          ICT13 = POSELT + int(IJROW - 1,8) * NFRONT8
          NBCOL = J4 - J3 + 1
Cduplicates--> CVD$ DEPCHK
CduplicatesCVD$ NODEPCHK
          DO 250 JJ = 1, NBCOL
            APOS3 = ICT13 + int(INTARR(J3 + JJ - 1) - 1,8)
            A(APOS3) = A(APOS3) + DBLARR(AINPUT + JJ - 1)
  250     CONTINUE
        ENDIF
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
  260 CONTINUE
      GOTO 500
  270 CONTINUE
      IFLAG = -8
      IERROR = LREQ
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE IN INTEGER ALLOCATION DURING SMUMPS_252'
      ENDIF
      GOTO 490
  280 CONTINUE
      IFLAG = -9
      CALL MUMPS_731(LAELL_REQ8 - LRLUS, IERROR)
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, WORKSPACE TOO SMALL DURING SMUMPS_252'
      ENDIF
      GOTO 490
  290 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &  ' FAILURE, SEND BUFFER TOO SMALL DURING SMUMPS_252'
      ENDIF
      IFLAG = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &  ' FAILURE, SEND BUFFER TOO SMALL DURING SMUMPS_252'
      ENDIF
      IFLAG = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     & ' FAILURE IN INTEGER DYNAMIC ALLOCATION DURING SMUMPS_252'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK + 1
  490 CALL  SMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_252
      SUBROUTINE SMUMPS_253(COMM_LOAD, ASS_IRECV,
     &    N, INODE, IW, LIW, A, LA, IFLAG,
     &    IERROR, ND, FILS, FRERE, DAD,
     &    CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW, root,
     &    OPASSW, OPELIW, PTRIST, PTLUST_S,  PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_S,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS, 
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP, KEEP8,INTARR,DBLARR,
     &    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, LEAF, IPOOL, LPOOL,
     &    PERM , MEM_DISTRIB)
      USE SMUMPS_COMM_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N,LIW,NSTEPS, NBFIN
      INTEGER(8) :: LA
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     &        LPOOL, LEAF, IWPOS, IWPOSCB, COMP
      INTEGER(8) :: LRLUS, LRLU, IPTRLU, POSFAC
      REAL A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER COMM, SLAVEF, MYID,  LBUFR, LBUFR_BYTES
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER IW(LIW), ITLOC(N+KEEP(253)),
     &        PTRARW(N), PTRAIW(N), ND(KEEP(28)),
     &        FILS(N), FRERE(KEEP(28)), DAD (KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), 
     & PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), PERM(N)
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER   CAND(SLAVEF+1, max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER NBPROCFILS(KEEP(28)),
     &        PROCNODE_STEPS(KEEP(28)), BUFR(LBUFR)
      INTEGER      INTARR(max(1,KEEP(14)))
      REAL DBLARR(max(1,KEEP(13)))
      INTEGER LP, HS, HF, HF_OLD,NCBSON, NSLAVES_OLD,
     &        NBSPLIT
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFS4FATHER,I
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER(8) :: LAELL8
      INTEGER LREQ_OOC
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      INTEGER NCB
      INTEGER J1,J2,J3,MP
      INTEGER(8) :: JJ8, LAPOS2, JJ2, JJ3
      INTEGER NELIM,JJ,JJ1,NPIVS,NCOLS,NROWS,
     &        IBROT,IORG
      INTEGER LDAFS, LDA_SON
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4, NUMORG_SPLIT
      INTEGER(8) :: ICT13
      INTEGER(8) :: IACHK, APOS, APOS2, POSELT, ICT12, POSEL1
      INTEGER AINPUT
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, ISLAVE, PDEST
      INTEGER PDEST1(1)
      INTEGER TYPESPLIT
      INTEGER ISON_IN_PLACE 
      LOGICAL IS_ofType5or6
      LOGICAL SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER IZERO
      INTEGER IDUMMY(1)
      PARAMETER( IZERO = 0 )
      INTEGER MUMPS_275, MUMPS_330, MUMPS_810
      EXTERNAL MUMPS_275, MUMPS_330, MUMPS_810
      REAL ZERO
      REAL RZERO
      PARAMETER(RZERO = 0.0E0 )
      PARAMETER( ZERO = 0.0E0 )
      INTEGER NELT, LPTRAR, NCBSON_MAX
      logical :: force_cand
      INTEGER ETATASS
      INCLUDE 'mumps_headers.h'
      INTEGER (8) :: APOSMAX
      REAL  MAXARR
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok, 
     &        NCB_SPLIT, SIZE_LIST_SPLIT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST, COPY_CAND
      INTEGER NBPANELS_L, NBPANELS_U
      MP         = ICNTL(2)
      IS_ofType5or6    = .FALSE.
      COMPRESSCB = .FALSE.
      ETATASS    = 0  
      IN         = INODE
      NBPROCFILS(STEP(IN)) = 0
      NSTEPS = NSTEPS + 1
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      ENDDO
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      NCBSON_MAX = 0
      NELT = 1
      LPTRAR = 1
      DO WHILE (ISON .GT. 0)
        NUMSTK = NUMSTK + 1
        IF ( KEEP(48)==5 .AND. 
     &       MUMPS_330(PROCNODE_STEPS(STEP(ISON)),
     &       SLAVEF) .EQ. 1) THEN
          NCBSON_MAX = max
     &      (
     &       IW(PIMASTER(STEP(ISON))+KEEP(IXSZ)), NCBSON_MAX
     &       )
        ENDIF
        NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 + KEEP(IXSZ))
        ISON = FRERE(STEP(ISON))
      ENDDO
      NFRONT = ND(STEP(INODE)) + NASS + KEEP(253)
      NASS1 = NASS + NUMORG
      NCB   = NFRONT - NASS1
      if((KEEP(24).eq.0).or.(KEEP(24).eq.1)) then
         force_cand=.FALSE.
      else
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
       IS_ofType5or6 =    (TYPESPLIT.EQ.5 .OR. TYPESPLIT.EQ.6)
       IF  ( (TYPESPLIT.EQ.4) 
     &               .OR.(TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) 
     &     )  THEN
        IF (TYPESPLIT.EQ.4) THEN
         ALLOCATE(COPY_CAND(SLAVEF+1),stat=allocok)
         IF (allocok > 0 ) THEN
           GOTO 265
         ENDIF
         CALL SMUMPS_791 (
     &      INODE, STEP, N, SLAVEF, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      CAND(1,INIV2), ICNTL, COPY_CAND,
     &      NBSPLIT, NUMORG_SPLIT, TMP_SLAVES_LIST(1),
     &      SIZE_TMP_SLAVES_LIST 
     &                                    )
         NCB_SPLIT = NCB-NUMORG_SPLIT
         SIZE_LIST_SPLIT = SIZE_TMP_SLAVES_LIST - NBSPLIT
         CALL SMUMPS_472( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, COPY_CAND,
     &     MEM_DISTRIB(0), NCB_SPLIT, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST(NBSPLIT+1),
     &     SIZE_LIST_SPLIT,INODE )
         DEALLOCATE (COPY_CAND)
         CALL SMUMPS_790 (
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
         CALL SMUMPS_792 (
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
        CALL SMUMPS_472( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, CAND(1,INIV2),
     &     MEM_DISTRIB(0), NCB, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST,
     &     SIZE_TMP_SLAVES_LIST,INODE )
       ENDIF
      HF   = NSLAVES + 6 + KEEP(IXSZ)
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1) THEN 
        CALL SMUMPS_684(KEEP(50), NASS1, NFRONT, NASS1,
     &                               NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL SMUMPS_94(N, KEEP(28),
     &        IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER,
     &        KEEP(216),LRLUS,KEEP(IXSZ))
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
      IF (KEEP(50).EQ.0 .AND. KEEP(234).EQ.0) THEN
        CALL  MUMPS_81(MYID, INODE, N, IOLDPS, HF, NFRONT,
     &        NFRONT_EFF, DAD,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     &        INTARR, ITLOC, RHS_MUMPS, FILS, FRERE,
     &        SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG,
     &        PROCNODE_STEPS, SLAVEF )
      ELSE
        ISON_IN_PLACE = -9999
        CALL MUMPS_86( MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM, DAD,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     &        INTARR, ITLOC, RHS_MUMPS, FILS, FRERE,
     &        SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG,
     &        ISON_IN_PLACE,
     &        PROCNODE_STEPS, SLAVEF)
        IF (IFLAG.LT.0) GOTO 250
      ENDIF
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (
     &        (TYPESPLIT.EQ.5) .OR. (TYPESPLIT.EQ.6)  ) THEN
          WRITE(6,*) ' Internal error 1 in fac_ass due to plitting ',
     &     ' INODE, NFRONT, NFRONT_EFF =', INODE, NFRONT, NFRONT_EFF 
          CALL MUMPS_ABORT()
        ENDIF
        IF (NFRONT.GT.NFRONT_EFF) THEN
            NCB    = NFRONT_EFF - NASS1
            NSLAVES_OLD = NSLAVES
            HF_OLD      = HF
            IF (TYPESPLIT.EQ.4) THEN
             WRITE(6,*) ' Internal error 2 in fac_ass due',
     &     ' to splitting ',
     &     ' INODE, NFRONT, NFRONT_EFF =', INODE, NFRONT, NFRONT_EFF 
             CALL MUMPS_ABORT()
            ELSE
             CALL SMUMPS_472( NCBSON_MAX,
     &       SLAVEF, KEEP,KEEP8, ICNTL,
     &       CAND(1,INIV2),
     &       MEM_DISTRIB(0), NCB, NFRONT_EFF, NSLAVES,
     &       TAB_POS_IN_PERE(1,INIV2),
     &       TMP_SLAVES_LIST, SIZE_TMP_SLAVES_LIST,INODE )
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
          Write(*,*) MYID,': ERROR 2 during ass_niv2, INODE=', INODE, 
     &     ' NFRONT, NFRONT_EFF=', NFRONT, NFRONT_EFF
          GOTO 270
        ENDIF
      ENDIF
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
        CALL SMUMPS_691(KEEP(50),
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
      IW(IOLDPS+6+KEEP(IXSZ):IOLDPS+5+KEEP(IXSZ)+NSLAVES)=
     &             TMP_SLAVES_LIST(1:NSLAVES)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE) 
      IF ( KEEP(73) .EQ. 0 ) THEN
#endif
#endif
      CALL SMUMPS_461(MYID, SLAVEF, COMM_LOAD,
     &     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &     NASS1, KEEP,KEEP8, IW(IOLDPS+6+KEEP(IXSZ)), NSLAVES,INODE)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE) 
      ENDIF
#endif
#endif
      IF(KEEP(86).EQ.1)THEN
         IF(mod(KEEP(24),2).eq.0)THEN
            CALL SMUMPS_533(SLAVEF,CAND(SLAVEF+1,INIV2),
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, NSLAVES,INODE)
         ELSEIF((KEEP(24).EQ.0).OR.(KEEP(24).EQ.1))THEN
            CALL SMUMPS_533(SLAVEF,SLAVEF-1,
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, NSLAVES,INODE)
         ENDIF
      ENDIF
      DEALLOCATE(TMP_SLAVES_LIST)
      IF (KEEP(50).EQ.0) THEN
        LAELL8 = int(NASS1,8) * int(NFRONT,8)
        LDAFS = NFRONT
      ELSE
        LAELL8 = int(NASS1,8)*int(NASS1,8)
        IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2)
     &     LAELL8 = LAELL8+int(NASS1,8)
        LDAFS = NASS1
      ENDIF
      IF (LRLU .LT. LAELL8) THEN
        IF (LRLUS .LT. LAELL8) THEN
          GOTO 280
        ELSE
         CALL SMUMPS_94(N, KEEP(28),
     &      IW, LIW, A, LA,
     &      LRLU, IPTRLU,
     &      IWPOS, IWPOSCB, PTRIST, PTRAST,
     &      STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &      KEEP(IXSZ))
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
      CALL SMUMPS_471(.FALSE.,.FALSE.,LA-LRLUS,0_8,LAELL8,
     &     KEEP,KEEP8,LRLU)
      POSEL1 = POSELT - int(LDAFS,8)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
        A(POSELT:LAPOS2) = ZERO
      ELSE
        APOS = POSELT
        DO JJ8 = 0_8, int(LDAFS-1,8)
          A(APOS:APOS+JJ8) = ZERO
          APOS = APOS + int(LDAFS,8)
        END DO
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
          A(APOS:APOS+int(LDAFS,8)-1_8)=ZERO
        ENDIF
      END IF
#endif
      IF ((NUMSTK.NE.0).AND.(NASS.NE.0)) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK = PIMASTER(STEP(ISON))
          NELIM = IW(ISTCHK + 1+KEEP(IXSZ))
          IF (NELIM.EQ.0) GOTO 210
          LSTK    = IW(ISTCHK+KEEP(IXSZ))
          NPIVS   = IW(ISTCHK + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS=0
          NSLSON  = IW(ISTCHK + 5+KEEP(IXSZ))
          HS      = 6 + NSLSON  + KEEP(IXSZ)
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          IF ( SAME_PROC ) THEN
           COMPRESSCB=( IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
           COMPRESSCB=( IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW(ISTCHK + 2+KEEP(IXSZ))
          ELSE
           NROWS = NCOLS
          ENDIF
          OPASSW = OPASSW + dble(NELIM*LSTK)
          J1 = ISTCHK + HS + NROWS + NPIVS
          J2 = J1 + NELIM - 1
          IACHK = PAMASTER(STEP(ISON))
          IF (KEEP(50).eq.0) THEN
           IF (IS_ofType5or6) THEN
            APOS = POSELT  
            DO JJ8 = 1_8, int(NELIM,8)*int(LSTK,8)
             A(APOS+JJ8-1_8) = A(APOS+JJ8-1_8) + A(IACHK+JJ8-1_8)
            ENDDO
           ELSE
            DO 170 JJ = J1, J2
             APOS = POSEL1 + int(IW(JJ),8) * int(LDAFS,8)
             DO 160 JJ1 = 1, LSTK
              JJ2 = APOS + int(IW(J1 + JJ1 - 1),8) - 1_8
              A(JJ2) = A(JJ2) + A(IACHK + int(JJ1 - 1,8))
  160        CONTINUE
             IACHK = IACHK + int(LSTK,8)
  170       CONTINUE
           ENDIF
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
            CALL SMUMPS_178( A, LA,
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
      IBROT = INODE
      APOSMAX = POSELT + int(NASS1,8)*int(NASS1,8)
      DO 260 IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        IJROW = INTARR(J1) 
        ICT12 = POSELT + int(IJROW - 1 - LDAFS, 8)
        MAXARR = RZERO
CduplicatesCVD$ NODEPCHK
        DO 240 JJ = J1, J2
          IF (KEEP(219).NE.0) THEN
            IF (INTARR(JJ).LE.NASS1) THEN
              APOS2 = ICT12 + int(INTARR(JJ),8) * int(LDAFS,8)
              A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            ELSEIF (KEEP(50).EQ.2) THEN
              MAXARR = max(MAXARR,abs(DBLARR(AINPUT)))
            ENDIF
          ELSE
            IF (INTARR(JJ).LE.NASS1) THEN
              APOS2 = ICT12 + int(INTARR(JJ),8) * int(LDAFS,8)
              A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            ENDIF
          ENDIF
          AINPUT = AINPUT + 1
  240   CONTINUE
        IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
           A(APOSMAX+int(IJROW-1,8)) = MAXARR
        ENDIF
        IF (J3 .GT. J4) GOTO 255
        ICT13 = POSELT + int(IJROW - 1,8) * int(LDAFS,8)
        NBCOL = J4 - J3 + 1
CduplicatesCVD$ NODEPCHK
CduplicatesCVD$ NODEPCHK
        DO JJ = 1, NBCOL
          JJ3 = ICT13 + int(INTARR(J3 + JJ - 1),8) - 1_8
          A(JJ3) = A(JJ3) + DBLARR(AINPUT + JJ - 1)
        ENDDO
  255   CONTINUE
        IF (KEEP(50).EQ.0) THEN
          DO JJ = 1, KEEP(253)
            APOS = POSELT +
     &             int(IJROW-1,8) * int(LDAFS,8) +
     &             int(LDAFS-KEEP(253)+JJ-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ENDIF
        IBROT = FILS(IBROT)
  260 CONTINUE
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
           CALL SMUMPS_68( INODE,
     &      NBPROCFILS(STEP(INODE)),
     &      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     &      IZERO, IDUMMY,
     &      IW(PDEST), NFRONT, COMM, IERR)
         ELSE
           NBCOL = NASS1+SHIFT_INDEX+NBLIG
           CALL SMUMPS_68( INODE,
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
          CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     &     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &     MPI_ANY_SOURCE, MPI_ANY_TAG,
     &     STATUS, BUFR, LBUFR,
     &     LBUFR_BYTES,
     &     PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &     LRLU, LRLUS, N, IW, LIW, A, LA,
     &     PTRIST, PTLUST_S, PTRFAC,
     &     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &     IERROR, COMM,
     &     NBPROCFILS,
     &     IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &     root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &     FILS, PTRARW, PTRAIW,
     &     INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &     LPTRAR, NELT, IW, IW,
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
          IF ( IFLAG .LT. 0 ) GOTO 500
          IF (MESSAGE_RECEIVED) THEN
           IOLDPS = PTLUST_S(STEP(INODE))
           PTRCOL = IOLDPS + HF + NFRONT
           PTRROW = IOLDPS + HF + NASS1 + SHIFT_INDEX
          ENDIF
         ENDIF
        ENDDO
        IF (IERR .EQ. -2) GOTO 300
        IF (IERR .EQ. -3) GOTO 305
        PTRROW = PTRROW + NBLIG
        PDEST  = PDEST + 1
      ENDDO
      IF (NUMSTK.EQ.0) GOTO 500
      ISON = IFSON
      DO IELL = 1, NUMSTK
        ISTCHK = PIMASTER(STEP(ISON))
        NELIM = IW(ISTCHK + 1+KEEP(IXSZ))
        LSTK    = IW(ISTCHK+KEEP(IXSZ))
        NPIVS   = IW(ISTCHK + 3+KEEP(IXSZ))
        IF ( NPIVS .LT. 0 ) NPIVS = 0
        NSLSON  = IW(ISTCHK + 5+KEEP(IXSZ))
        HS      = 6 + NSLSON + KEEP(IXSZ)
        NCOLS     = NPIVS + LSTK
        SAME_PROC     = (ISTCHK.LE.IWPOS)
        IF (.NOT.SAME_PROC) THEN
         NROWS = IW(ISTCHK + 2+KEEP(IXSZ))
        ELSE
         NROWS = NCOLS
        ENDIF
        PDEST   = ISTCHK + 6 + KEEP(IXSZ)
        NCBSON  = LSTK - NELIM
        PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
        IF (KEEP(219).NE.0) THEN
          IF(KEEP(50) .EQ. 2) THEN
           NFS4FATHER = NCBSON
           DO I=0,NCBSON-1
              IF(IW(PTRCOL+I) .GT. NASS1) THEN
                 NFS4FATHER = I
                 EXIT
              ENDIF
           ENDDO
           NFS4FATHER = NFS4FATHER+NELIM
          ELSE
           NFS4FATHER = 0
          ENDIF
        ELSE
          NFS4FATHER = 0
        ENDIF
        IF (NSLSON.EQ.0) THEN
          NSLSON = 1
          PDEST1(1)  = MUMPS_275(PROCNODE_STEPS(STEP(ISON)),
     &                                SLAVEF)
          IF (PDEST1(1).EQ.MYID) THEN
            CALL SMUMPS_211( COMM_LOAD, ASS_IRECV, 
     &      BUFR, LBUFR, LBUFR_BYTES,
     &      INODE, ISON, NSLAVES, 
     &      IW( PTLUST_S(STEP(INODE)) + 6 +KEEP(IXSZ)),
     &      NFRONT, NASS1, NFS4FATHER, NCBSON,
     &           IW( PTRCOL ),
     &      PROCNODE_STEPS,
     &      SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &      LRLUS, N, IW,
     &      LIW, A, LA,
     &      PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &      PIMASTER, PAMASTER, NSTK_S, COMP,
     &      IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     &      NBFIN, ICNTL, KEEP,KEEP8, root,
     &      OPASSW, OPELIW,
     &      ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &      ND, FRERE, LPTRAR, NELT, IW, IW,
     &
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
           IF ( IFLAG .LT. 0 ) GOTO 500
          ELSE
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON))+HS+NROWS+NPIVS+NELIM
            CALL  SMUMPS_71( 
     &           INODE, NFRONT,NASS1,NFS4FATHER, 
     &           ISON, MYID,
     &      NSLAVES, IW( PTLUST_S(STEP(INODE))+6+KEEP(IXSZ) ),
     &      IW(PTRCOL), NCBSON,
     &      COMM, IERR, PDEST1, NSLSON, SLAVEF, 
     &      KEEP,KEEP8, STEP, N, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
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
     &        root,OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &        FILS, PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     &        NELT, IW, IW, 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           ENDDO
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
              CALL SMUMPS_210( COMM_LOAD, ASS_IRECV, 
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
     &        INTARR, DBLARR, ND, FRERE, LPTRAR, NELT, IW,
     &        IW, 
     &        
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &        )
              IF ( IFLAG .LT. 0 ) GOTO 500
              EXIT
            ENDIF
          ENDDO
          IF (PIMASTER(STEP(ISON)).GT.0) THEN
          IERR = -1
          DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  =  PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  SMUMPS_71( 
     &           INODE, NFRONT,NASS1, NFS4FATHER,
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
             CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
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
     &        LPTRAR, NELT, IW, IW, 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
          ENDDO
          IF (IERR .EQ. -2) GOTO 290
          IF (IERR .EQ. -3) GOTO 295
          ENDIF
        ENDIF
       ISON = FRERE(STEP(ISON))
      ENDDO
      GOTO 500
  250 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE IN INTEGER DYNAMIC ALLOCATION DURING
     & SMUMPS_253'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK + 1
      GOTO 490
  265 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE ALLOCATING TMP_SLAVES_LIST',
     &                 ' DURING SMUMPS_253'
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
     &' FAILURE IN INTEGER ALLOCATION DURING SMUMPS_253'
      ENDIF
      GOTO 490
  280 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, WORKSPACE TOO SMALL DURING SMUMPS_253'
      ENDIF
      IFLAG = -9
      CALL MUMPS_731(LAELL8-LRLUS, IERROR)
      GOTO 490
  290 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (1) DURING SMUMPS_253'
      ENDIF
      IFLAG = -17
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (1) DURING SMUMPS_253'
      ENDIF
      IFLAG = -20
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (2) DURING SMUMPS_253'
      ENDIF
      IFLAG = -17
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  305 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (2) DURING SMUMPS_253'
      ENDIF
      IFLAG = -17
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      IERROR =  LREQ  * KEEP( 34 ) 
  490 CALL SMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_253
      SUBROUTINE SMUMPS_39(N, INODE, IW, LIW, A, LA, 
     &    ISON, NBROWS, NBCOLS, ROWLIST,
     &    VALSON, PTLUST_S, PTRAST, STEP, PIMASTER,
     &    OPASSW, IWPOSCB, MYID, KEEP,KEEP8, IS_ofType5or6,
     &    LDA_VALSON )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA
      INTEGER N,LIW,MYID
      INTEGER INODE,ISON, IWPOSCB
      INTEGER NBROWS, NBCOLS, LDA_VALSON
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER IW(LIW), STEP(N), PIMASTER(KEEP(28)),
     &        PTLUST_S(KEEP(28)), ROWLIST(NBROWS)
      REAL A(LA), VALSON(LDA_VALSON,NBROWS)
      DOUBLE PRECISION OPASSW
      LOGICAL, INTENT(IN) :: IS_ofType5or6
      INTEGER(8) :: POSELT, POSEL1, APOS, JJ2
      INTEGER HF,HS, NSLAVES, NFRONT, NASS1,
     &        IOLDPS, ISTCHK, LSTK, NSLSON,NELIM,
     &        NPIVS,NCOLS,J1,JJ,JJ1,NROWS,
     &        LDAFS_PERE, IBEG, DIAG
      INCLUDE 'mumps_headers.h'
      LOGICAL SAME_PROC
      INTRINSIC real
      IOLDPS = PTLUST_S(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      NFRONT = IW(IOLDPS+KEEP(IXSZ))
      NASS1  = iabs(IW(IOLDPS + 2+KEEP(IXSZ)))
      NSLAVES= IW(IOLDPS+5+KEEP(IXSZ))
      IF (KEEP(50).EQ.0) THEN
        LDAFS_PERE = NFRONT
      ELSE
        IF ( NSLAVES .eq. 0 ) THEN
          LDAFS_PERE = NFRONT
        ELSE
          LDAFS_PERE = NASS1
        ENDIF
      ENDIF
      HF      = 6 + NSLAVES + KEEP(IXSZ)
      POSEL1 = POSELT - int(LDAFS_PERE,8)
      ISTCHK = PIMASTER(STEP(ISON))
      LSTK = IW(ISTCHK+KEEP(IXSZ))
      NSLSON  = IW(ISTCHK + 5+KEEP(IXSZ))
      HS      = 6 + NSLSON + KEEP(IXSZ)
      OPASSW = OPASSW + dble(NBROWS*NBCOLS)
      NELIM = IW(ISTCHK + 1+KEEP(IXSZ))
      NPIVS = IW(ISTCHK + 3+KEEP(IXSZ))
      IF (NPIVS.LT.0) NPIVS = 0
      NCOLS = NPIVS + LSTK
      SAME_PROC = (ISTCHK.LT.IWPOSCB)
      IF (SAME_PROC) THEN
       NROWS = NCOLS
      ELSE
       NROWS = IW(ISTCHK+2+KEEP(IXSZ))
      ENDIF
      J1 = ISTCHK + NROWS + HS + NPIVS
      IF (KEEP(50).EQ.0) THEN
       IF (IS_ofType5or6) THEN
         APOS = POSEL1 + int(ROWLIST(1),8) * int(LDAFS_PERE,8)
         DO JJ = 1, NBROWS
           DO JJ1 = 1, NBCOLS
             JJ2 = APOS + int(JJ1-1,8)
             A(JJ2)=A(JJ2)+VALSON(JJ1,JJ)
           ENDDO
           APOS = APOS + int(LDAFS_PERE,8)
         ENDDO
       ELSE
        DO 170 JJ = 1, NBROWS
         APOS = POSEL1 + int(ROWLIST(JJ),8) * int(LDAFS_PERE,8)
         DO 160 JJ1 = 1, NBCOLS
           JJ2 = APOS + int(IW(J1 + JJ1 - 1) - 1,8)
           A(JJ2) = A(JJ2) + VALSON(JJ1,JJ) 
  160    CONTINUE
  170   CONTINUE
       ENDIF
      ELSE
       IF (IS_ofType5or6) THEN
        APOS = POSEL1 + int(ROWLIST(1),8) * int(LDAFS_PERE,8)
        DIAG = ROWLIST(1)
        DO JJ = 1, NBROWS
          DO JJ1 = 1, DIAG
            JJ2 = APOS+int(JJ1-1,8)
            A(JJ2) = A(JJ2) + VALSON(JJ1,JJ)
          ENDDO
          DIAG = DIAG+1
          APOS = APOS + int(LDAFS_PERE,8)
        ENDDO
       ELSE
        DO JJ = 1, NBROWS
         IF (ROWLIST(JJ).LE.NASS1.and..NOT.IS_ofType5or6) THEN
          APOS = POSEL1 + int(ROWLIST(JJ) - 1,8)
          DO JJ1 = 1, NELIM
           JJ2 = APOS + int(IW(J1+JJ1-1),8)*int(LDAFS_PERE,8)
            A(JJ2) = A(JJ2) + VALSON(JJ1,JJ)
          ENDDO
          IBEG = NELIM+1
         ELSE
          IBEG = 1
         ENDIF
         APOS = POSEL1 + int(ROWLIST(JJ),8) * int(LDAFS_PERE,8)
         DO JJ1 = IBEG, NBCOLS
          IF (ROWLIST(JJ).LT.IW(J1 + JJ1 - 1)) EXIT
          JJ2 = APOS + int(IW(J1 + JJ1 - 1) - 1,8)
          A(JJ2) = A(JJ2) + VALSON(JJ1,JJ)
         ENDDO
        ENDDO
       ENDIF  
      ENDIF   
      RETURN
      END SUBROUTINE SMUMPS_39
      SUBROUTINE SMUMPS_539
     &    (N, INODE, IW, LIW, A, LA, 
     &    NBROWS, NBCOLS,
     &    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     &    RHS_MUMPS, FILS, PTRARW, PTRAIW, INTARR, DBLARR, 
     &    ICNTL, KEEP,KEEP8, MYID)
      IMPLICIT NONE
      INTEGER N,LIW
      INTEGER(8) :: LA
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER INODE, MYID
      INTEGER NBROWS, NBCOLS 
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER IW(LIW), ITLOC(N+KEEP(253)), STEP(N),
     &        PTRIST(KEEP(28)), FILS(N), PTRARW(N), PTRAIW(N)
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      REAL A(LA),
     &        DBLARR(max(1,KEEP(13)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER(8) :: POSELT, ICT12, APOS
      INTEGER IOLDPS, NBCOLF, NBROWF, NSLAVES, HF,
     &        K1,K2,K,J,JPOS,NASS,JJ,
     &        IN,AINPUT,JK,J1,J2,IJROW, ILOC
      INTEGER    :: K1RHS, K2RHS, JFirstRHS
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INCLUDE 'mumps_headers.h'
      IOLDPS  = PTRIST(STEP(INODE))
      POSELT  = PTRAST(STEP(INODE))
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NASS    = IW(IOLDPS+1+KEEP(IXSZ))
      NSLAVES = IW(IOLDPS+5+KEEP(IXSZ))
      HF      = 6 + NSLAVES + KEEP(IXSZ)
      IF (NASS.LT.0) THEN
          NASS         = -NASS
          IW(IOLDPS+1+KEEP(IXSZ)) = NASS
          A(POSELT:POSELT+int(NBROWF,8)*int(NBCOLF,8)-1_8) =
     &    ZERO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NASS - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -JPOS
           JPOS     = JPOS + 1
          ENDDO
          K1 = IOLDPS + HF 
          K2 = K1 + NBROWF - 1
          JPOS = 1
          IF ((KEEP(253).GT.0).AND.(KEEP(50).NE.0)) THEN
           K1RHS = 0
           K2RHS = -1
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = JPOS
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
                ILOC = ITLOC(J)    
                APOS = POSELT+int(ILOC-1,8)*int(NBCOLF,8) + 
     &                 int(IJROW-1,8) 
                A(APOS) = A(APOS) + RHS_MUMPS(
     &                    (JFirstRHS+(K-K1RHS)-1)*KEEP(254)+IN)
              ENDDO
              IN = FILS(IN)
             ENDDO
            ENDIF
          ELSE  
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = JPOS
            JPOS     = JPOS + 1
           ENDDO
          ENDIF
          IN = INODE
          DO WHILE (IN.GT.0) 
           AINPUT = PTRARW(IN)
           JK     = PTRAIW(IN)
           JJ     = JK + 1
           J1     = JJ + 1
           J2 = J1 + INTARR(JK)
           IJROW = -ITLOC(INTARR(J1))
           ICT12 = POSELT +int(- NBCOLF + IJROW - 1,8)
           DO JJ= J1,J2
            ILOC = ITLOC(INTARR(JJ))
            IF (ILOC.GT.0) THEN
              APOS = ICT12 + int(ILOC,8)*int(NBCOLF,8)
              A(APOS) = A(APOS) + DBLARR(AINPUT)
            ENDIF
            AINPUT  = AINPUT + 1
           ENDDO
           IN = FILS(IN)
          ENDDO
          K1 = IOLDPS + HF
          K2 = K1 + NBROWF + NASS - 1
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          ENDDO
      ENDIF
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_539
      SUBROUTINE SMUMPS_531
     & (N, INODE, IW, LIW, NBROWS, STEP, PTRIST,
     & ITLOC, RHS_MUMPS, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER N, LIW
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER INODE
      INTEGER NBROWS
      INTEGER IW(LIW), ITLOC(N+KEEP(253)), STEP(N),
     &        PTRIST(KEEP(28))
      REAL :: RHS_MUMPS(KEEP(255))
      INCLUDE 'mumps_headers.h'
      INTEGER IOLDPS, NBCOLF, NBROWF, NSLAVES, HF,
     &        K1,K2,K,J
      IOLDPS  = PTRIST(STEP(INODE))
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NSLAVES = IW(IOLDPS+5+KEEP(IXSZ))
      HF      = 6 + NSLAVES+KEEP(IXSZ)
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = 0
          ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_531
      SUBROUTINE SMUMPS_40(N, INODE, IW, LIW, A, LA, 
     &    NBROWS, NBCOLS, ROWLIST, COLLIST, VALSON, 
     &    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     &    RHS_MUMPS, FILS,
     &    ICNTL, KEEP,KEEP8, MYID, IS_ofType5or6, LDA_VALSON)
      IMPLICIT NONE
      INTEGER N,LIW
      INTEGER(8) :: LA
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER INODE, MYID
      LOGICAL, intent(in) :: IS_ofType5or6
      INTEGER NBROWS, NBCOLS, LDA_VALSON
      INTEGER ROWLIST(NBROWS), COLLIST(NBCOLS)
      INTEGER IW(LIW), ITLOC(N+KEEP(253)), STEP(N),
     &        PTRIST(KEEP(28)), FILS(N)
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRAST(KEEP(28))
      REAL A(LA), VALSON(LDA_VALSON,NBROWS)
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER(8) :: POSEL1, POSELT, APOS, K8
      INTEGER IOLDPS, NBCOLF, NBROWF, NSLAVES, HF,
     &        I,J,NASS,IDIAG
      INCLUDE 'mumps_headers.h'
      INTRINSIC real
      IOLDPS  = PTRIST(STEP(INODE))
      POSELT  = PTRAST(STEP(INODE))
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NASS    = IW(IOLDPS+1+KEEP(IXSZ))
       IF ( NBROWS .GT. NBROWF ) THEN
          WRITE(*,*) ' ERR: ERROR : NBROWS > NBROWF'
          WRITE(*,*) ' ERR: INODE =', INODE
          WRITE(*,*) ' ERR: NBROW=',NBROWS,'NBROWF=',NBROWF
          WRITE(*,*) ' ERR: ROW_LIST=', ROWLIST
          CALL MUMPS_ABORT()
       END IF
      NSLAVES = IW(IOLDPS+5+KEEP(IXSZ))
      HF      = 6 + NSLAVES+KEEP(IXSZ)
      IF (NBROWS.GT.0) THEN
          POSEL1 = POSELT - int(NBCOLF,8)
          IF (KEEP(50).EQ.0) THEN
           IF (IS_ofType5or6) THEN
            APOS = POSEL1 + int(ROWLIST(1),8) * int(NBCOLF,8)
            DO I=1, NBROWS
                DO J = 1, NBCOLS
                  A(APOS+int(J-1,8)) = A( APOS+int(J-1,8)) + VALSON(J,I)
                ENDDO
                APOS = APOS + int(NBCOLF,8)
            END DO
           ELSE
            DO I=1,NBROWS
             APOS = POSEL1 + int(ROWLIST(I),8) * int(NBCOLF,8)
             DO J=1,NBCOLS
              K8 = APOS + int(ITLOC(COLLIST(J)),8) - 1_8
              A(K8) = A(K8) + VALSON(J,I)
             ENDDO
            ENDDO
           ENDIF
          ELSE
           IF (IS_ofType5or6) THEN
            APOS = POSEL1 + int(ROWLIST(1),8) * int(NBCOLF,8)
     &              + int((NBROWS-1),8)*int(NBCOLF,8)
            IDIAG = 0
            DO I=NBROWS,1,-1
             A(APOS:APOS+int(NBCOLS-IDIAG-1,8))= 
     &         A(APOS:APOS+int(NBCOLS-IDIAG-1,8)) +
     &         VALSON(1:NBCOLS-IDIAG,I)
             APOS = APOS - int(NBCOLF,8)
             IDIAG = IDIAG + 1
            ENDDO
           ELSE
            DO I=1,NBROWS
             APOS = POSEL1 + int(ROWLIST(I),8) * int(NBCOLF,8)
             DO J=1,NBCOLS
              IF (ITLOC(COLLIST(J)) .EQ. 0) THEN 
                  write(6,*) ' .. exit for col =', J
                  EXIT
              ENDIF
              K8 = APOS + int(ITLOC(COLLIST(J)),8) - 1_8
              A(K8) = A(K8) + VALSON(J,I)
             ENDDO
            ENDDO
           ENDIF
          ENDIF
          OPASSW = OPASSW + dble(NBROWS*NBCOLS)
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_40
      SUBROUTINE SMUMPS_178( A, LA,
     &             IAFATH, NFRONT, NASS1,
     &             IACB, NCOLS, LCB,
     &             IW, NROWS, NELIM, ETATASS,
     &             CB_IS_COMPRESSED, IS_INPLACE
     &             )
      IMPLICIT NONE
      INTEGER NFRONT, NASS1
      INTEGER(8) :: LA
      INTEGER NCOLS, NROWS, NELIM
      INTEGER(8) :: LCB
      REAL A( LA )
      INTEGER(8) :: IAFATH, IACB
      INTEGER IW( NCOLS )
      INTEGER ETATASS
      LOGICAL CB_IS_COMPRESSED, IS_INPLACE
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      LOGICAL RESET_TO_ZERO, RISK_OF_SAME_POS,
     &        RISK_OF_SAME_POS_THIS_LINE, OMP_FLAG
      INTEGER I, J
      INTEGER(8) :: APOS, POSELT
      INTEGER(8) :: IPOSCB, IBEGCBROW, IENDFRONT
      IENDFRONT =  IAFATH+int(NFRONT,8)*int(NFRONT,8)-1_8
      IF ( IS_INPLACE ) THEN
        IPOSCB=1_8
        RESET_TO_ZERO    = IACB .LT. IENDFRONT + 1_8
        RISK_OF_SAME_POS = IACB + LCB .EQ. IENDFRONT + 1_8
        RISK_OF_SAME_POS_THIS_LINE = .FALSE.
        DO I=1, NROWS
          POSELT = int(IW(I)-1,8) * int(NFRONT,8)
          IF (.NOT. CB_IS_COMPRESSED ) THEN
            IPOSCB = 1_8 + int(I - 1,8) * int(NCOLS,8)
            IF (IACB+IPOSCB-1_8 .GE. IENDFRONT + 1_8) THEN
              RESET_TO_ZERO = .FALSE.
            ENDIF
          ENDIF
          IF ( RISK_OF_SAME_POS ) THEN
            IF (I.EQ.NROWS .OR. .NOT. CB_IS_COMPRESSED) THEN
              IF ( IAFATH + POSELT + int(IW(I)-1,8) .EQ.
     &             IACB+IPOSCB+int(I-1-1,8)) THEN
                 RISK_OF_SAME_POS_THIS_LINE = .TRUE.
              ENDIF
            ENDIF
          ENDIF
          IF (RESET_TO_ZERO) THEN
            IF ( RISK_OF_SAME_POS_THIS_LINE ) THEN
              DO J=1, I
                APOS = POSELT + int(IW( J ),8)
                IF (IAFATH + APOS - 1_8.NE. IACB+IPOSCB-1_8) THEN
                  A(IAFATH+ APOS -1_8) = A(IACB+IPOSCB-1_8)
                  A(IACB+IPOSCB-1_8) = ZERO
                ENDIF
                IPOSCB = IPOSCB + 1_8
              ENDDO
            ELSE
              DO J=1, I
                APOS = POSELT + int(IW( J ),8)
                A(IAFATH+ APOS -1_8) = A(IACB+IPOSCB-1_8)
                A(IACB+IPOSCB-1_8) = ZERO
                IPOSCB = IPOSCB + 1_8
              ENDDO
            ENDIF
          ELSE
            DO J=1, I
              APOS = POSELT + int(IW( J ),8)
              A(IAFATH+ APOS -1_8) = A(IACB+IPOSCB-1_8)
              IPOSCB = IPOSCB + 1_8
            ENDDO
          ENDIF
          IF (.NOT. CB_IS_COMPRESSED ) THEN
            IBEGCBROW = IACB+IPOSCB-1_8
            IF ( IBEGCBROW .LE. IENDFRONT ) THEN
              A(IBEGCBROW:IBEGCBROW+int(NCOLS-I,8)-1_8)=ZERO
            ENDIF
          ENDIF
          IF (IACB+IPOSCB-1_8 .GE. IENDFRONT + 1_8) THEN
            RESET_TO_ZERO = .FALSE.
          ENDIF
        ENDDO
        RETURN
      ENDIF
      IF ((ETATASS.EQ.0) .OR. (ETATASS.EQ.1)) THEN
        IPOSCB = 1_8
        DO I = 1, NELIM
          POSELT = int( IW( I ) - 1, 8 ) * int(NFRONT, 8)
          IF (.NOT. CB_IS_COMPRESSED) THEN
            IPOSCB = 1_8 + int( I - 1, 8 ) * int(NCOLS,8)
          ENDIF
          DO J = 1, I
            APOS = POSELT + int(IW( J ),8)
            A(IAFATH+ APOS -1_8) = A(IAFATH+ APOS -1_8)
     &                           + A(IACB+IPOSCB-1_8)
            IPOSCB = IPOSCB + 1_8
          END DO
        END DO
      ENDIF
      IF ((ETATASS.EQ.0).OR.(ETATASS.EQ.1)) THEN
        OMP_FLAG = (NROWS-NELIM).GE.300
       DO I = NELIM + 1, NROWS
          IF (CB_IS_COMPRESSED) THEN
            IPOSCB = (int(I,8) * int(I-1,8)) / 2_8 + 1_8
          ELSE
            IPOSCB = int(I-1,8) * int(NCOLS,8) + 1_8
          ENDIF
          POSELT = int(IW( I ),8)
          IF (POSELT.LE. int(NASS1,8) .AND. .NOT. IS_INPLACE) THEN 
            DO J = 1, NELIM
              APOS = POSELT + int( IW( J ) - 1, 8 ) * int(NFRONT,8)
              A(IAFATH+APOS-1_8) = A(IAFATH+APOS-1_8) +
     &                             A(IACB+IPOSCB-1_8)
              IPOSCB = IPOSCB + 1_8
            END DO
          ELSE
            POSELT = int( IW( I ) - 1, 8 ) * int(NFRONT, 8)
            DO J = 1, NELIM
             APOS = POSELT + int(IW( J ), 8)
             A(IAFATH+APOS-1_8) = A(IAFATH+APOS-1_8)
     &                          + A(IACB+IPOSCB-1_8)
             IPOSCB = IPOSCB + 1_8
            END DO
          ENDIF
          IF (ETATASS.EQ.1) THEN
            POSELT = int( IW( I ) - 1, 8 ) * int(NFRONT, 8)
            DO J = NELIM + 1, I
                 IF (IW(J).GT.NASS1) EXIT
                 APOS = POSELT + int(IW( J ), 8)
                 A(IAFATH+APOS-1_8) = A(IAFATH+APOS-1_8)
     &                              + A(IACB+IPOSCB-1_8)
                 IPOSCB = IPOSCB +1_8
            END DO
          ELSE
            POSELT = int( IW( I ) - 1, 8 ) * int(NFRONT, 8)
            DO J = NELIM + 1, I
             APOS = POSELT + int(IW( J ), 8)
             A(IAFATH+APOS-1_8) = A(IAFATH+APOS-1_8)
     &                          + A(IACB+IPOSCB-1_8)
             IPOSCB = IPOSCB + 1_8
            END DO
          ENDIF
        END DO
      ELSE  
        DO I= NROWS, NELIM+1, -1
          IF (CB_IS_COMPRESSED) THEN
            IPOSCB = (int(I,8)*int(I+1,8))/2_8 
          ELSE
            IPOSCB = int(I-1,8) * int(NCOLS,8) + int(I,8)
          ENDIF
          POSELT = int(IW( I ),8)
          IF (POSELT.LE.int(NASS1,8)) EXIT
          POSELT = int( IW( I ) - 1, 8 ) * int(NFRONT, 8)
          DO J=I,NELIM+1, -1
            IF (IW(J).LE.NASS1) EXIT
            APOS = POSELT + int(IW( J ), 8)
            A(IAFATH+APOS-1_8) = A(IAFATH+APOS-1_8)
     &                         + A(IACB+IPOSCB-1_8)
            IPOSCB = IPOSCB - 1_8
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_178
      SUBROUTINE SMUMPS_530(N, ISON, INODE, IWPOSCB,
     &           PIMASTER, PTLUST_S, IW, LIW, STEP, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER N, ISON, INODE, IWPOSCB
      INTEGER KEEP(500), STEP(N)
      INTEGER(8) KEEP8(150)
      INTEGER PIMASTER(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER LIW
      INTEGER IW(LIW)
      INTEGER ISTCHK, LSTK, NSLSON, HS, NROWS, NCOLS, NPIVS, NELIM
      INTEGER IOLDPS, NFRONT, NSLAVES, ICT11, HF
      INTEGER J1, J2, J3, JJ, JPOS
      LOGICAL SAME_PROC
      INCLUDE 'mumps_headers.h'
      ISTCHK = PIMASTER(STEP(ISON))
      LSTK   = IW(ISTCHK+KEEP(IXSZ))
      NSLSON = IW(ISTCHK+5+KEEP(IXSZ))
      HS     = 6 + NSLSON + KEEP(IXSZ)
      NELIM  = IW(ISTCHK + 1+KEEP(IXSZ))
      NPIVS  = IW(ISTCHK + 3+KEEP(IXSZ))
      NCOLS  = NPIVS + LSTK
      IF ( NPIVS < 0 ) NPIVS = 0
      SAME_PROC = ISTCHK < IWPOSCB
      IF (SAME_PROC) THEN
       NROWS = NCOLS
      ELSE
       NROWS = IW(ISTCHK+2+KEEP(IXSZ))
      ENDIF
      J1 = ISTCHK + NROWS + HS + NPIVS
      IF (KEEP(50).NE.0) THEN
          J2 = J1 +  LSTK - 1
          DO JJ = J1, J2
            IW(JJ) = IW(JJ - NROWS)
          ENDDO
      ELSE
            J2 = J1 + LSTK - 1
            J3 = J1 + NELIM
            DO JJ = J3, J2
             IW(JJ) = IW(JJ - NROWS)
            ENDDO
            IF (NELIM .NE. 0) THEN
              IOLDPS = PTLUST_S(STEP(INODE))
              NFRONT = IW(IOLDPS+KEEP(IXSZ))
              NSLAVES= IW(IOLDPS+5+KEEP(IXSZ))
              HF     = 6 + NSLAVES+KEEP(IXSZ)
              ICT11 = IOLDPS + HF - 1 + NFRONT
              J3 = J3 - 1
              DO 190 JJ = J1, J3
               JPOS = IW(JJ) + ICT11
               IW(JJ) = IW(JPOS)
  190         CONTINUE
            ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_530
      SUBROUTINE SMUMPS_619(
     &     N, INODE, IW, LIW, A, LA, 
     &     ISON, NBCOLS,
     &     VALSON, PTLUST_S, PTRAST, STEP, PIMASTER,
     &     OPASSW, IWPOSCB,MYID, KEEP,KEEP8 )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA
      INTEGER N,LIW,MYID
      INTEGER INODE,ISON,IWPOSCB
      INTEGER NBCOLS
      INTEGER IW(LIW), STEP(N), 
     &     PIMASTER(KEEP(28)),
     &     PTLUST_S(KEEP(28))
      INTEGER(8) PTRAST(KEEP(28))
      REAL A(LA)
      REAL VALSON(NBCOLS)
      DOUBLE PRECISION OPASSW
      INTEGER HF,HS, NSLAVES, NASS1,
     &     IOLDPS, ISTCHK,
     &     LSTK, NSLSON,NELIM,NPIVS,NCOLS, J1,
     &     JJ1,NROWS
      INTEGER(8) POSELT, APOS, JJ2
      INCLUDE 'mumps_headers.h'
      LOGICAL SAME_PROC
      INTRINSIC real
      IOLDPS = PTLUST_S(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      NASS1  = iabs(IW(IOLDPS + 2 + KEEP(IXSZ)))
      NSLAVES= IW(IOLDPS+5 + KEEP(IXSZ))
      HF      = 6 + NSLAVES + KEEP(IXSZ)
      ISTCHK = PIMASTER(STEP(ISON))
      LSTK = IW(ISTCHK + KEEP(IXSZ))
      NSLSON  = IW(ISTCHK + 5 + KEEP(IXSZ))
      HS      = 6 + NSLSON + KEEP(IXSZ)
      NELIM = IW(ISTCHK + 1 + KEEP(IXSZ))
      NPIVS = IW(ISTCHK + 3 + KEEP(IXSZ))
      IF (NPIVS.LT.0) NPIVS = 0
      NCOLS = NPIVS + LSTK
      SAME_PROC = (ISTCHK.LT.IWPOSCB)
      IF (SAME_PROC) THEN
       NROWS = NCOLS
      ELSE
       NROWS = IW(ISTCHK+2 + KEEP(IXSZ))
      ENDIF
      J1 = ISTCHK + NROWS + HS + NPIVS
      APOS = POSELT + int(NASS1,8)*int(NASS1,8) - 1_8
      DO JJ1 = 1, NBCOLS
         JJ2 = APOS+int(IW(J1 + JJ1 - 1),8)
         IF(abs(A(JJ2)) .LT. VALSON(JJ1))
     &         A(JJ2) = VALSON(JJ1)
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_619
      RECURSIVE SUBROUTINE SMUMPS_264(
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
     &    )
      USE SMUMPS_OOC
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA
      INTEGER(8) :: POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER COMM, MYID
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTLUST_S(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), ND(KEEP(28))
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      REAL  DBLARR(max(1,KEEP(13)))
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER INODE, POSITION, NPIV, IERR, LP
      INTEGER NCOL
      INTEGER(8) :: POSBLOCFACTO
      INTEGER(8) :: LAELL
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW
      INTEGER (8) :: LPOS, LPOS1, LPOS2, IPOS, KPOS
      INTEGER ICT11
      INTEGER I, IPIV, FPERE
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      REAL ONE,ALPHA
      PARAMETER (ONE = 1.0E0, ALPHA=-1.0E0)
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPivDummy
      TYPE(IO_BLOCK) :: MonBloc
      LOGICAL LAST_CALL
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      FPERE    = -1
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LAELL = int(NPIV,8) * int(NCOL,8)
      IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          CALL MUMPS_731(LAELL - LRLUS, IERROR)
          IF (ICNTL(1).GT.0 .AND. ICNTL(4).GE.1) THEN
            LP=ICNTL(1)
            WRITE(LP,*)
     &" FAILURE, WORKSPACE TOO SMALL DURING SMUMPS_264"
          ENDIF
          GOTO 700
        END IF
        CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     &      LRLU, IPTRLU,
     &      IWPOS, IWPOSCB, PTRIST, PTRAST,
     &      STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &      KEEP(IXSZ))
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     &       ,LRLU,LRLUS
             IFLAG = -9
             CALL MUMPS_731( LAELL-LRLUS, IERROR )
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IF (ICNTL(1).GT.0 .AND. ICNTL(4).GE.1) THEN
            LP=ICNTL(1)
            WRITE(LP,*)
     &" FAILURE IN INTEGER ALLOCATION DURING SMUMPS_264"
          ENDIF
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
      END IF
      LRLU  = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL SMUMPS_471(.FALSE., .FALSE.,
     &               LA-LRLUS,0_8,LAELL,KEEP,KEEP8,LRLU)
      IPIV = IWPOS
      IWPOS = IWPOS + NPIV
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*NCOL, 
     &                 MPI_REAL,
     &                 COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
         DO WHILE ( PTRIST(STEP(INODE)) .EQ. 0 )
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL SMUMPS_329( COMM_LOAD,
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
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &    )
          IF ( IFLAG .LT. 0 ) GOTO 600
        END DO
      ENDIF
      DO WHILE ( NBPROCFILS( STEP(INODE)) .NE. 0 ) 
        BLOCKING = .TRUE.
        SET_IRECV = .FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL SMUMPS_329( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
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
        CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
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
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.  )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS +KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 +KEEP(IXSZ))
      NROW1  = IW( IOLDPS + 2 +KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ))
      NSLAV1 = IW( IOLDPS + 5 +KEEP(IXSZ))
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      IF (NPIV.GT.0) THEN
        ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
        DO I = 1, NPIV
          IF (IW(IPIV+I-1).EQ.I) CYCLE
          ISW = IW(ICT11+I)
          IW(ICT11+I) = IW(ICT11+IW(IPIV+I-1))
          IW(ICT11+IW(IPIV+I-1)) = ISW
          IPOS = POSELT + int(NPIV1 + I - 1,8)
          KPOS = POSELT + int(NPIV1 + IW(IPIV+I-1) - 1,8)
          CALL sswap(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
        ENDDO
        LPOS2 = POSELT + int(NPIV1,8)
        CALL strsm('L','L','N','N',NPIV, NROW1, ONE, 
     &           A(POSBLOCFACTO), NCOL, A(LPOS2), NCOL1)
        LPOS1 = POSBLOCFACTO+int(NPIV,8)
        LPOS  = LPOS2 + int(NPIV,8)
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
        LAST_CALL = .FALSE.
        CALL SMUMPS_688( STRAT, TYPEF_L, A(POSELT),
     &       LAFAC, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
      ENDIF
      IF ( NPIV .GT. 0 ) THEN
        CALL sgemm('N','N', NCOL-NPIV,NROW1,NPIV,
     &             ALPHA,A(LPOS1),NCOL,
     &             A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
      ENDIF
      IW(IOLDPS+KEEP(IXSZ) ) = IW(IOLDPS+KEEP(IXSZ) ) - NPIV
      IW(IOLDPS + 3+KEEP(IXSZ) ) = IW(IOLDPS+3+KEEP(IXSZ) ) + NPIV
      IF (LASTBL) IW(IOLDPS+1+KEEP(IXSZ) ) = IW(IOLDPS + 3+KEEP(IXSZ) )
      IF ( .not. LASTBL .AND. 
     &  (IW(IOLDPS+1+KEEP(IXSZ)) .EQ. IW(IOLDPS + 3+KEEP(IXSZ))) ) THEN
        write(*,*) ' ERROR 1 **** IN BLACFACTO '
        CALL MUMPS_ABORT()
      ENDIF
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      POSFAC = POSFAC - LAELL
      CALL SMUMPS_471(.FALSE.,.FALSE.,
     &             LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLU)
      IWPOS = IWPOS - NPIV
      FLOP1 = dble( NPIV1*NROW1 ) +
     &        dble(NROW1*NPIV1)*dble(2*NCOL1-NPIV1-1)
     &   -
     &        dble((NPIV1+NPIV)*NROW1 ) -
     &        dble(NROW1*(NPIV1+NPIV))*dble(2*NCOL1-NPIV1-NPIV-1)
      CALL SMUMPS_190( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      IF (LASTBL) THEN
         CALL SMUMPS_759(
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
      ENDIF 
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_264
      SUBROUTINE SMUMPS_699( COMM_LOAD, ASS_IRECV, 
     &   MSGLEN, BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, POSFAC,
     &   N, IW, LIW, A, LA, PTRIST, PTLUST_S, PTRFAC, PTRAST,
     &   STEP, PIMASTER, PAMASTER, NBPROCFILS,
     &   COMP, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, NSTK_S,
     &   FILS, PTRARW, PTRAIW, INTARR, DBLARR, NBFIN,
     &   MYID, COMM, ICNTL, KEEP,KEEP8, IFLAG, IERROR,
     &   IPOOL, LPOOL, LEAF, ND, FRERE_STEPS, LPTRAR, NELT,
     &   FRTPTR, FRTELT, 
     &   ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE SMUMPS_LOAD
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV, MSGLEN
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: LRLU, IPTRLU, LRLUS, LA, POSFAC
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER NBFIN
      INTEGER COMP
      INTEGER NELT, LPTRAR
      INTEGER PROCNODE_STEPS( KEEP(28) ), PTRIST(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER PTLUST_S( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER ITLOC( N + KEEP(253)), NSTK_S( KEEP(28) ), FILS( N )
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER ND(KEEP(28)), FRERE_STEPS( KEEP(28) )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER INTARR( max(1,KEEP(14)) )
      REAL DBLARR( max( 1,KEEP(13)) )
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER COMM, MYID, IFLAG, IERROR
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER FRTPTR(N+1), FRTELT( NELT )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER NFS4FATHER
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER MUMPS_275, MUMPS_810
      EXTERNAL MUMPS_275, MUMPS_810
      INTEGER IERR
      INTEGER NBROWS_ALREADY_SENT, NBROWS_PACKET
      INTEGER I, INODE, ISON, POSITION, NBROW, LROW, IROW, INDCOL
      INTEGER LREQI
      INTEGER(8) :: LREQA, POSCONTRIB
      INTEGER ROW_LENGTH
      INTEGER MASTER
      INTEGER ISTCHK
      LOGICAL SAME_PROC
      LOGICAL SLAVE_NODE
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED, IS_ofType5or6
      INTEGER ISHIFT_BUFR, LBUFR_LOC, LBUFR_BYTES_LOC
      INTEGER TYPESPLIT
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, ISON, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NBROW, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, LROW, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NBROWS_ALREADY_SENT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NBROWS_PACKET, 1,
     &                 MPI_INTEGER, COMM, IERR )
      MASTER     = MUMPS_275(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      SLAVE_NODE = MASTER .NE. MYID
      TYPESPLIT = MUMPS_810(PROCNODE_STEPS(STEP(INODE)),
     &                  SLAVEF)
      IS_ofType5or6 = ((TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6))
      IF (SLAVE_NODE .AND. PTRIST(STEP(INODE)) ==0) THEN
        ISHIFT_BUFR     = ( MSGLEN + KEEP(34) ) / KEEP(34)
        LBUFR_LOC       = LBUFR - ISHIFT_BUFR + 1
        LBUFR_BYTES_LOC = LBUFR_LOC * KEEP(34)
        DO WHILE ( PTRIST( STEP(INODE) ) .EQ. 0 )
          MASTER = MUMPS_275(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
          BLOCKING = .TRUE.
          SET_IRECV = .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     &     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &     MASTER, MAITRE_DESC_BANDE,
     &     STATUS, 
     &     BUFR(ISHIFT_BUFR), LBUFR_LOC, LBUFR_BYTES_LOC,
     &     PROCNODE_STEPS, POSFAC,
     &     IWPOS, IWPOSCB, IPTRLU,
     &     LRLU, LRLUS, N, IW, LIW, A, LA,
     &     PTRIST, PTLUST_S, PTRFAC,
     &     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &     IFLAG, IERROR, COMM,
     &     NBPROCFILS, IPOOL, LPOOL, LEAF,
     &     NBFIN, MYID, SLAVEF,
     &
     &     root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, 
     &     PTRARW, PTRAIW,
     &     INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     &     LPTRAR, NELT, FRTPTR, FRTELT, 
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
          IF (IFLAG.LT.0) RETURN
        END DO
      ENDIF
      IF ( SLAVE_NODE ) THEN
         LREQI = LROW + NBROWS_PACKET
      ELSE
         LREQI = NBROWS_PACKET
      END IF
         LREQA = int(LROW,8)
         IF ( LRLU .LT. LREQA .OR. IWPOS + LREQI
     &        - 1 .GT. IWPOSCB ) THEN
            IF ( LRLUS .LT. LREQA ) THEN
               IFLAG = -9
               CALL MUMPS_731( LREQA - LRLUS, IERROR )
               CALL SMUMPS_44( MYID, SLAVEF, COMM )
               RETURN
            END IF
            CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &           KEEP(IXSZ))
            COMP = COMP+1
            IF ( LRLU .NE. LRLUS ) THEN
               WRITE(*,*) 'PB compress ass..process_contrib'
               WRITE(*,*) 'LRLU,LRLUS=',LRLU,LRLUS
               IFLAG = -9
               CALL MUMPS_731( LREQA - LRLUS, IERROR )
               CALL SMUMPS_44( MYID, SLAVEF, COMM )
               RETURN
            END IF
            IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
               IFLAG  = -8
               IERROR = IWPOS + LREQI - 1 - IWPOSCB
               CALL SMUMPS_44( MYID, SLAVEF, COMM )
               RETURN
            END IF
         END IF
         LRLU  = LRLU - LREQA
         LRLUS = LRLUS - LREQA
         POSCONTRIB = POSFAC
         POSFAC = POSFAC + LREQA
         KEEP8(67) = min(LRLUS, KEEP8(67))
         CALL SMUMPS_471(.FALSE.,.FALSE.,
     &        LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLU)
         IF  ( SLAVE_NODE ) THEN
            IROW   = IWPOS
            INDCOL = IWPOS + NBROWS_PACKET
         ELSE
            IROW   = IWPOS
            INDCOL = -1
         END IF
         IWPOS = IWPOS + LREQI
         IF ( SLAVE_NODE ) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           IW( INDCOL ), LROW, MPI_INTEGER,
     &           COMM, IERR )
         END IF
         DO I = 1, NBROWS_PACKET
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           IW( IROW + I - 1 ), 1, MPI_INTEGER,
     &           COMM, IERR )
         END DO
         IF ( SLAVE_NODE ) THEN
            IF ( NBROWS_ALREADY_SENT + NBROWS_PACKET == NBROW ) THEN
              NBPROCFILS(STEP(INODE))=NBPROCFILS(STEP(INODE))-1
            ENDIF
            IF ( KEEP(55) .eq. 0 ) THEN               
               CALL SMUMPS_539
     &              (N, INODE, IW, LIW, A, LA,
     &              NBROW, LROW,
     &              OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &              ITLOC, RHS_MUMPS,
     &              FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &              KEEP,KEEP8, MYID )
            ELSE
               CALL SMUMPS_123(
     &              NELT, FRTPTR, FRTELT,
     &              N, INODE, IW, LIW, A, LA,
     &              NBROW, LROW,
     &              OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &              ITLOC, RHS_MUMPS,
     &              FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &              KEEP,KEEP8, MYID )
            ENDIF
            DO I=1,NBROWS_PACKET
               IF(KEEP(50).NE.0)THEN
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ROW_LENGTH,
     &                 1,
     &                 MPI_INTEGER,
     &                 COMM, IERR )
               ELSE
                 ROW_LENGTH=LROW
               ENDIF
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              A(POSCONTRIB),
     &              ROW_LENGTH,
     &              MPI_REAL,
     &              COMM, IERR )
               CALL SMUMPS_40(N, INODE, IW, LIW, A, LA,
     &              1, ROW_LENGTH, IW( IROW+I-1 ),IW(INDCOL),
     &              A(POSCONTRIB),
     &              OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &              ITLOC, RHS_MUMPS,
     &              FILS, ICNTL, KEEP,KEEP8, MYID, IS_ofType5or6, 
     &              ROW_LENGTH )
            ENDDO
            CALL SMUMPS_531
     &           (N, INODE, IW, LIW,
     &           NBROWS_PACKET, STEP, PTRIST,
     &           ITLOC, RHS_MUMPS,KEEP,KEEP8)
         ELSE
            DO I=1,NBROWS_PACKET
               IF(KEEP(50).NE.0)THEN
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ROW_LENGTH,
     &                 1,
     &                 MPI_INTEGER,
     &                 COMM, IERR )
               ELSE
                 ROW_LENGTH=LROW
               ENDIF
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              A(POSCONTRIB),
     &              ROW_LENGTH,
     &              MPI_REAL,
     &              COMM, IERR )
               CALL SMUMPS_39(N, INODE, IW, LIW, A, LA,
     &              ISON, 1, ROW_LENGTH, IW( IROW +I-1 ),
     &              A(POSCONTRIB), PTLUST_S, PTRAST,
     &              STEP, PIMASTER, OPASSW,
     &              IWPOSCB, MYID, KEEP,KEEP8,
     &              IS_ofType5or6, ROW_LENGTH
     &)
            ENDDO
          IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          IF (KEEP(219).NE.0) THEN
            IF(KEEP(50) .EQ. 2) THEN
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              NFS4FATHER,
     &              1,
     &              MPI_INTEGER,
     &              COMM, IERR )
               IF(NFS4FATHER .GT. 0) THEN
                  CALL SMUMPS_617(NFS4FATHER,IERR)
                  IF (IERR .NE. 0) THEN
                        IERROR         = BUF_LMAX_ARRAY
                        IFLAG          = -13
                        CALL SMUMPS_44( MYID, SLAVEF, COMM )
                        RETURN
                  ENDIF
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 BUF_MAX_ARRAY,
     &                 NFS4FATHER,
     &                 MPI_REAL,
     &                 COMM, IERR )
                  CALL SMUMPS_619(N, INODE, IW, LIW, A, LA,
     &                 ISON, NFS4FATHER,
     &                 BUF_MAX_ARRAY, PTLUST_S, PTRAST,
     &                 STEP, PIMASTER, OPASSW,
     &                 IWPOSCB, MYID, KEEP,KEEP8)
               ENDIF
            ENDIF
          ENDIF
          ENDIF
          IF (NBROWS_ALREADY_SENT + NBROWS_PACKET == NBROW ) THEN
            NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) - 1
            NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
            IF ( NBPROCFILS(STEP(ISON)) .EQ. 0) THEN
               ISTCHK = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK .LT. IWPOSCB
               IF (SAME_PROC) THEN
                  CALL SMUMPS_530(N, ISON, INODE, IWPOSCB,
     &                 PIMASTER, PTLUST_S, IW, LIW, STEP, KEEP,KEEP8)
               ENDIF
               IF (SAME_PROC) THEN
                  ISTCHK = PTRIST(STEP(ISON))
                  PTRIST(STEP( ISON) ) = -99999999
               ELSE
                  PIMASTER(STEP( ISON )) = -99999999
               ENDIF
               CALL SMUMPS_152(.FALSE., MYID, N, ISTCHK,
     &              PAMASTER(STEP(ISON)),
     &              IW, LIW, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &              LA, KEEP,KEEP8, .FALSE.
     &              )
            ENDIF
            IF ( NBPROCFILS(STEP(INODE)) .EQ. 0 ) THEN
               CALL SMUMPS_507( N, IPOOL, LPOOL,
     &              PROCNODE_STEPS,
     &              SLAVEF, KEEP(28), KEEP(76), KEEP(80),
     &              KEEP(47), STEP, INODE+N )
               IF (KEEP(47) .GE. 3) THEN
                  CALL SMUMPS_500(
     &          IPOOL, LPOOL, 
     &                 PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &                 MYID, STEP, N, ND, FILS )
               ENDIF
            ENDIF
          ENDIF 
      END IF 
         IWPOS = IWPOS - LREQI
         LRLU = LRLU + LREQA
         LRLUS = LRLUS + LREQA
         POSFAC = POSFAC - LREQA
         CALL SMUMPS_471(.FALSE.,.FALSE.,
     &        LA-LRLUS,0_8,-LREQA,KEEP,KEEP8,LRLU)
      RETURN
      END SUBROUTINE SMUMPS_699
      SUBROUTINE SMUMPS_143( N, INODE, IW, LIW, A, LA,
     &                           IOLDPS, POSELT, IFLAG, UU, NOFFW,
     &                           NPVW,
     &                           KEEP,KEEP8, STEP,
     &                           PROCNODE_STEPS, MYID, SLAVEF, SEUIL,
     &                           AVOID_DELAYED, ETATASS,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, 
     &     IWPOS )
      USE SMUMPS_OOC      
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, NOFFW, NPVW
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER MYID, SLAVEF, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) ), STEP(N)
      REAL UU, SEUIL
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK
      INTEGER NASS, NEL1, NPIVB, NPIVE, NBOLKJ, NBTLKJ
      REAL UUTEMP
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, LNextPiv2beWritten, 
     &        UNextPiv2beWritten, IFLAG_OOC,
     &        PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &        PP_LastPIVRPTRFilled_L, 
     &        PP_LastPIVRPTRFilled_U
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INCLUDE 'mumps_headers.h'
      EXTERNAL MUMPS_330, SMUMPS_221, SMUMPS_233, 
     &         SMUMPS_229,
     &         SMUMPS_225, SMUMPS_232, SMUMPS_231,
     &         SMUMPS_220,
     &         SMUMPS_228, SMUMPS_236
      INTEGER  MUMPS_330
      LOGICAL STATICMODE
      REAL SEUIL_LOC
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
      IBEG_BLOCK=1
      NFRONT = IW(IOLDPS+KEEP(IXSZ))
      NASS   = iabs(IW(IOLDPS+2+KEEP(IXSZ)))
      IF (NASS .GT. KEEP(3)) THEN
        NBOLKJ = min( KEEP(6), NASS )
      ELSE
        NBOLKJ = min( KEEP(5), NASS )
      ENDIF
      NBTLKJ = NBOLKJ
        IF (KEEP(201).EQ.1) THEN 
          CALL MUMPS_729(LAFAC,IW(IOLDPS+XXR))
          LIWFAC    = IW(IOLDPS+XXI)
          TYPEFile     = TYPEF_BOTH_LU  
          LNextPiv2beWritten = 1 
          UNextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = LNextPiv2beWritten 
          PP_FIRST2SWAP_U = UNextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          MonBloc%LastPanelWritten_U = 0 
          PP_LastPIVRPTRFilled_L = 0 
          PP_LastPIVRPTRFilled_U = 0 
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 1
          MonBloc%NROW     = NFRONT
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -88877    
          NULLIFY(MonBloc%INDICES)   
        ENDIF
 50   CONTINUE
      CALL SMUMPS_221(NFRONT,NASS,N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     &     IFLAG,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     &     DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &     PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &     PP_LastPIVRPTRFilled_U)
      IF (IFLAG.LT.0) GOTO 500  
      IF (INOPV.EQ.1) THEN
         IF(STATICMODE) THEN
            INOPV = -1
            GOTO 50  
         ENDIF
         GOTO 80
      ENDIF
      IF (INOPV.EQ.2) THEN
         CALL SMUMPS_233(IBEG_BLOCK,
     &            NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &            IOLDPS,POSELT,NBOLKJ, NBTLKJ,KEEP(4),KEEP(IXSZ))
         GOTO 50
      ENDIF
      NPVW = NPVW + 1
      IF (NASS.LE.1) THEN
       CALL SMUMPS_229(NFRONT,N,INODE,IW,LIW,A,LA,
     &                 IOLDPS,POSELT,KEEP(IXSZ))
       IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+1+KEEP(IXSZ)) + 1
       GO TO 500
      ENDIF
       CALL SMUMPS_225(IBEG_BLOCK,NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     &             IOLDPS,POSELT,IFINB,
     &             NBTLKJ,KEEP(4),KEEP(IXSZ))
       IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+1+KEEP(IXSZ)) + 1
       IF (IFINB.EQ.0) GOTO 50
       IF (KEEP(201).EQ.1) THEN  
           MonBloc%LastPiv= IW(IOLDPS+1+KEEP(IXSZ))
           STRAT          = STRAT_TRY_WRITE
           TYPEFile       = TYPEF_U  
           LAST_CALL      = .FALSE.
           CALL SMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
        ENDIF
       IF (IFINB.EQ.(-1)) GOTO 80
       NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
       NEL1   = NASS - NPIV
      CALL SMUMPS_232(A,LA,
     &           NFRONT,NPIV,NASS,POSELT,NBTLKJ)
      GO TO 50
 80   CONTINUE
      NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
      IF (NPIV.LE.0) GO TO 110
      NEL1   = NFRONT - NASS
      IF (NEL1.LE.0) GO TO 110
        IF (KEEP(201).EQ.1) THEN 
         STRAT          = STRAT_TRY_WRITE
         TYPEFile     = TYPEF_BOTH_LU  
         MonBloc%LastPiv= NPIV
         CALL SMUMPS_642(A(POSELT), LAFAC, NFRONT, 
     &      NPIV, NASS, IW(IOLDPS), LIWFAC, 
     &      MonBloc, TYPEFile, MYID, KEEP8,
     &      STRAT, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten)
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
        ELSE
          CALL SMUMPS_231(A,LA,NFRONT, NPIV,NASS,POSELT)
        ENDIF
 110  CONTINUE
      IF (MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
     &                   .EQ.1) THEN
        NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
        IBEG_BLOCK = NPIV
        IF (NASS.EQ.NPIV) GOTO 500
 120    CALL SMUMPS_220(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &     INOPV,NOFFW,IOLDPS,POSELT,UU,SEUIL,
     &     KEEP, DKEEP,
     &     PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &     PP_LastPIVRPTRFilled_U)
        IF (INOPV.NE.1) THEN
         NPVW = NPVW + 1
         CALL SMUMPS_228(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &                 IOLDPS,POSELT,IFINB,KEEP(IXSZ))
         IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS+1+KEEP(IXSZ)) + 1
       IF (IFINB.EQ.0) GOTO 120
        ENDIF
        NPIV   = IW(IOLDPS+1+KEEP(IXSZ))
        NPIVB  = IBEG_BLOCK
        NPIVE  = NPIV - NPIVB
        NEL1   = NFRONT - NASS
        IF ((NPIVE.LE.0).OR.(NEL1.EQ.0)) GO TO 500
        CALL SMUMPS_236(A,LA,NPIVB,
     &                NFRONT,NPIV,NASS,POSELT)
      ENDIF
 500  CONTINUE
       IF (KEEP(201).EQ.1) THEN 
          STRAT            = STRAT_WRITE_MAX   
          MonBloc%Last     = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+KEEP(IXSZ))
          TYPEFile     = TYPEF_BOTH_LU  
          LAST_CALL    = .TRUE.
          CALL SMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
          CALL SMUMPS_644 (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
       ENDIF
      RETURN
      END SUBROUTINE SMUMPS_143
      RECURSIVE SUBROUTINE SMUMPS_322(
     &    COMM_LOAD, ASS_IRECV,
     &    MSGSOU, MSGTAG, MSGLEN,
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
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER MSGSOU, MSGTAG, MSGLEN
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER KEEP(500), ICNTL( 40 )
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS, LA
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER INTARR( max(1,KEEP(14)) )
      REAL DBLARR( max(1,KEEP(13)) )
      INTEGER INIV2, ISHIFT, IBEG
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      LOGICAL FLAG
      INTEGER MP, LP
      INTEGER TMP( 2 )
      INTEGER NBRECU, POSITION, INODE, ISON, IROOT
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE,
     &     LMAP, FPERE, NELIM,
     &     HDMAPLIG,NFS4FATHER,
     &     TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      DOUBLE PRECISION FLOP1
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      CHARACTER(LEN=35)::SUBNAME
      MP = ICNTL(2)
      LP = ICNTL(1)
      SUBNAME="??????"
      CALL SMUMPS_467(COMM_LOAD, KEEP)
      IF ( MSGTAG .EQ. RACINE ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NBRECU,
     &     1, MPI_INTEGER, COMM, IERR)
          NBRECU = BUFR( 1 )
          NBFIN =  NBFIN - NBRECU
      ELSEIF ( MSGTAG .EQ. NOEUD ) THEN
          CALL SMUMPS_269( MYID,KEEP,KEEP8,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST, PTRAST,
     &    STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP, FPERE, FLAG, IFLAG, IERROR, COMM,
     &    ITLOC, RHS_MUMPS )
          SUBNAME="SMUMPS_269"
          IF ( IFLAG .LT. 0 ) GO TO 500
          IF ( FLAG ) THEN
            CALL SMUMPS_507(N, IPOOL, LPOOL,
     &           PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     &           KEEP(80), KEEP(47), STEP, FPERE )
            IF (KEEP(47) .GE. 3) THEN
               CALL SMUMPS_500(
     &              IPOOL, LPOOL,
     &              PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &              MYID, STEP, N, ND, FILS )
            ENDIF
            CALL MUMPS_137( FPERE, N,
     &           PROCNODE_STEPS,SLAVEF,
     &           ND, FILS, FRERE, STEP, PIMASTER,
     &           KEEP(28), KEEP(50), KEEP(253), FLOP1,
     &           IW, LIW, KEEP(IXSZ) )
            IF (FPERE.NE.KEEP(20))
     &        CALL SMUMPS_190(1,.FALSE.,FLOP1,KEEP,KEEP8)
          ENDIF
      ELSEIF ( MSGTAG .EQ. END_NIV2_LDLT ) THEN
          INODE = BUFR( 1 )
          CALL SMUMPS_507(N, IPOOL, LPOOL,
     &         PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     &         KEEP(80), KEEP(47),
     &         STEP, -INODE )
          IF (KEEP(47) .GE. 3) THEN
             CALL SMUMPS_500(
     &            IPOOL, LPOOL,
     &            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &            MYID, STEP, N, ND, FILS )
          ENDIF
      ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          IFLAG  = -001
          IERROR = MSGSOU
          GOTO 100
      ELSEIF ( MSGTAG .EQ. MAITRE_DESC_BANDE ) THEN
        CALL SMUMPS_266( MYID,BUFR, LBUFR,
     &    LBUFR_BYTES, IWPOS,
     &    IWPOSCB,
     &    IPTRLU, LRLU, LRLUS, NBPROCFILS,
     &    N, IW, LIW, A, LA,
     &    PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, COMP,
     &    KEEP,KEEP8, ITLOC, RHS_MUMPS,
     &    IFLAG, IERROR )
          SUBNAME="SMUMPS_266"
        IF ( IFLAG .LT. 0 ) GO to 500
      ELSEIF ( MSGTAG .EQ. MAITRE2           ) THEN
        CALL SMUMPS_268(MYID,BUFR, LBUFR, LBUFR_BYTES,
     &    PROCNODE_STEPS, SLAVEF, IWPOS, IWPOSCB,
     &    IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA,
     &    PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    KEEP,KEEP8, ND, FILS, FRERE, ITLOC, RHS_MUMPS,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          SUBNAME="SMUMPS_268"
        IF ( IFLAG .LT. 0 ) GO to 500
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO        ) THEN
        CALL SMUMPS_264( COMM_LOAD, ASS_IRECV,
     &   BUFR,  LBUFR, LBUFR_BYTES,
     &   PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM , IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO_SYM_SLAVE    ) THEN
        CALL SMUMPS_263( COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO_SYM    ) THEN
        CALL SMUMPS_274( COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      ELSEIF ( MSGTAG .EQ. CONTRIB_TYPE2    ) THEN
        CALL SMUMPS_699( COMM_LOAD, ASS_IRECV,
     &       MSGLEN, BUFR, LBUFR,
     &       LBUFR_BYTES, PROCNODE_STEPS,
     &       SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, POSFAC,
     &       N, IW, LIW, A, LA, PTRIST,
     &       PTLUST_S, PTRFAC, PTRAST,
     &       STEP, PIMASTER, PAMASTER, NBPROCFILS, COMP, root,
     &       OPASSW, OPELIW, ITLOC, RHS_MUMPS, NSTK_S,
     &       FILS, PTRARW, PTRAIW, INTARR, DBLARR, NBFIN, MYID, COMM,
     &       ICNTL, KEEP,KEEP8, IFLAG, IERROR, IPOOL, LPOOL, LEAF,
     &       ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF ( IFLAG .LT. 0 ) GO TO 100
      ELSEIF ( MSGTAG .EQ. MAPLIG            ) THEN
         HDMAPLIG = 7
         INODE        = BUFR( 1 )
         ISON         = BUFR( 2 )
         NSLAVES_PERE = BUFR( 3 )
         NFRONT_PERE  = BUFR( 4 )
         NASS_PERE    = BUFR( 5 )
         LMAP         = BUFR( 6 )
         NFS4FATHER = BUFR(7)
         IF ( (NSLAVES_PERE.NE.0).AND.(KEEP(48).NE.0) ) THEN
            INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
            ISHIFT = NSLAVES_PERE+1
            TAB_POS_IN_PERE(1:NSLAVES_PERE+1, INIV2) =
     &           BUFR(HDMAPLIG+1:HDMAPLIG+1+NSLAVES_PERE)
            TAB_POS_IN_PERE(SLAVEF+2, INIV2) = NSLAVES_PERE
         ELSE
            ISHIFT = 0
         ENDIF
         IBEG = HDMAPLIG+1+ISHIFT
         CALL SMUMPS_210( COMM_LOAD, ASS_IRECV,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    INODE, ISON, NSLAVES_PERE,
     &    BUFR(IBEG),
     &    NFRONT_PERE, NASS_PERE, NFS4FATHER,LMAP,
     &    BUFR(IBEG+NSLAVES_PERE),
     &    PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB,
     &    IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP,
     &    IFLAG, IERROR, MYID, COMM, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, ICNTL, KEEP,KEEP8, root,
     &    OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &    ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
         IF ( IFLAG .LT. 0 ) GO TO 100
      ELSE IF ( MSGTAG .EQ. ROOT_CONT_STATIC ) THEN
        CALL SMUMPS_700(
     &        BUFR, LBUFR, LBUFR_BYTES,
     &        root, N, IW, LIW, A, LA, NBPROCFILS,
     &        LRLU, IPTRLU, IWPOS, IWPOSCB,
     &        PTRIST, PTLUST_S, PTRFAC, PTRAST,
     &        STEP, PIMASTER, PAMASTER,
     &        COMP, LRLUS, IPOOL, LPOOL, LEAF,
     &        FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &        KEEP,KEEP8, IFLAG, IERROR, COMM, COMM_LOAD,
     &        ITLOC, RHS_MUMPS,
     &        ND, PROCNODE_STEPS, SLAVEF)
        SUBNAME="SMUMPS_700"
        IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. ROOT_NON_ELIM_CB ) THEN
        IROOT  = KEEP( 38 )
        MSGSOU = MUMPS_275( PROCNODE_STEPS(STEP(IROOT)), 
     &           SLAVEF )
        IF ( PTLUST_S( STEP(IROOT)) .EQ. 0 ) THEN
          CALL MPI_RECV( TMP, 2 * KEEP(34), MPI_PACKED,
     &                   MSGSOU, ROOT_2SLAVE,
     &                   COMM, STATUS, IERR )
          CALL SMUMPS_270( TMP( 1 ), TMP( 2 ),
     &    root,
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
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND )
          SUBNAME="SMUMPS_270"
          IF ( IFLAG .LT. 0 ) GOTO 500
        END IF
        CALL SMUMPS_700(
     &       BUFR, LBUFR, LBUFR_BYTES,
     &       root, N, IW, LIW, A, LA, NBPROCFILS,
     &       LRLU, IPTRLU, IWPOS, IWPOSCB,
     &       PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     &       COMP, LRLUS, IPOOL, LPOOL, LEAF,
     &       FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &       KEEP,KEEP8, IFLAG, IERROR, COMM, COMM_LOAD,
     &       ITLOC, RHS_MUMPS,
     &       ND, PROCNODE_STEPS, SLAVEF )
          SUBNAME="SMUMPS_700"
        IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. ROOT_2SON ) THEN
         ISON = BUFR( 1 )
         NELIM = BUFR( 2 )
         CALL SMUMPS_271( COMM_LOAD, ASS_IRECV,
     &    ISON, NELIM, root,
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
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GO TO 100
         IF (MYID.NE.MUMPS_275(PROCNODE_STEPS(STEP(ISON)), 
     &          SLAVEF)) THEN
          IF (KEEP(50).EQ.0) THEN
            IF (IW(PTRIST(STEP(ISON))+6+KEEP(IXSZ)).EQ.
     &                                 S_REC_CONTSTATIC) THEN
             IW(PTRIST(STEP(ISON))+6+KEEP(IXSZ)) = S_ROOT2SON_CALLED
            ELSE
             CALL SMUMPS_626( N, ISON, PTRIST, PTRAST,
     &       IW, LIW, A, LA, LRLU, LRLUS, IWPOSCB,
     &       IPTRLU, STEP, MYID, KEEP
     &    )
            ENDIF
          ELSE
           IF (IW(PTRIST(STEP(ISON))+8+KEEP(IXSZ)).EQ.
     &                                 S_REC_CONTSTATIC) THEN
             IW(PTRIST(STEP(ISON))+8+KEEP(IXSZ)) = S_ROOT2SON_CALLED
           ELSE
             CALL SMUMPS_626( N, ISON, PTRIST, PTRAST,
     &       IW, LIW, A, LA, LRLU, LRLUS, IWPOSCB,
     &       IPTRLU, STEP, MYID, KEEP
     &    )
           ENDIF
          ENDIF
         ENDIF
      ELSE IF ( MSGTAG .EQ. ROOT_2SLAVE ) THEN
          TOT_ROOT_SIZE    = BUFR( 1 )
          TOT_CONT_TO_RECV = BUFR( 2 )
          CALL SMUMPS_270( TOT_ROOT_SIZE,
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
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND )
          IF ( IFLAG .LT. 0 ) GO TO 100
      ELSE IF ( MSGTAG .EQ. ROOT_NELIM_INDICES ) THEN
         ISON         = BUFR( 1 )
         NELIM        = BUFR( 2 )
         NSLAVES_PERE = BUFR( 3 )
         CALL SMUMPS_273( root,
     &    ISON, NELIM, NSLAVES_PERE, BUFR(4), BUFR(4+BUFR(2)),
     &    BUFR(4+2*BUFR(2)),
     &
     &    PROCNODE_STEPS,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S,
     &    ITLOC, RHS_MUMPS, COMP,
     &    IFLAG, IERROR,
     &    IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,KEEP8,
     &    COMM, COMM_LOAD, FILS, ND)
          SUBNAME="SMUMPS_273"
         IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. UPDATE_LOAD ) THEN
         WRITE(*,*) "Internal error 3 in SMUMPS_322"
         CALL MUMPS_ABORT()
      ELSE IF ( MSGTAG .EQ. TAG_DUMMY   ) THEN
      ELSE
         IF ( LP > 0 )
     &     WRITE(LP,*) MYID,
     &': Internal error, routine SMUMPS_322.',MSGTAG
         IFLAG = -100
         IERROR= MSGTAG
         GOTO 500
      ENDIF
 100  CONTINUE
      RETURN
 500  CONTINUE
      IF ( ICNTL(1) .GT. 0 .AND. ICNTL(4).GE.1 ) THEN
        LP=ICNTL(1)
        IF (IFLAG.EQ.-9) THEN
         WRITE(LP,*) 'FAILURE, WORKSPACE TOO SMALL DURING ',SUBNAME
        ENDIF
        IF (IFLAG.EQ.-8) THEN
         WRITE(LP,*) 'FAILURE IN INTEGER ALLOCATION DURING ',SUBNAME
        ENDIF
        IF (IFLAG.EQ.-13) THEN
         WRITE(LP,*) 'FAILURE IN DYNAMIC ALLOCATION DURING ',SUBNAME
        ENDIF
      ENDIF
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_322
      RECURSIVE SUBROUTINE SMUMPS_280(
     &    COMM_LOAD, ASS_IRECV,
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
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT ,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: POSFAC, LA, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) ),
     &        PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER INTARR( max(1,KEEP(14)) )
      REAL DBLARR( max(1,KEEP(13)) )
      INTEGER MSGSOU, MSGTAG, MSGLEN, IERR
      MSGSOU = STATUS( MPI_SOURCE )
      MSGTAG = STATUS( MPI_TAG )
      CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
      IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
        IFLAG  = -20
        IERROR = MSGLEN
         WRITE(*,*) ' RECEPTION BUF TOO SMALL, Msgtag/len=',
     &                MSGTAG,MSGLEN
        CALL SMUMPS_44( MYID, SLAVEF, COMM )
        RETURN
       ENDIF
       CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     &                 MSGTAG,
     &                 COMM, STATUS, IERR )
       CALL SMUMPS_322(
     &      COMM_LOAD, ASS_IRECV,
     &      MSGSOU, MSGTAG, MSGLEN, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &      IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
      RETURN
      END SUBROUTINE SMUMPS_280
      RECURSIVE SUBROUTINE SMUMPS_329(
     &    COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV,
     &    MESSAGE_RECEIVED, MSGSOU, MSGTAG,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    STACK_RIGHT_AUTHORIZED )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL, INTENT (IN)  :: BLOCKING
      LOGICAL, INTENT (IN)  :: SET_IRECV
      LOGICAL, INTENT (INOUT) :: MESSAGE_RECEIVED
      INTEGER, INTENT (IN) :: MSGSOU, MSGTAG
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: LA, POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) ),
     &        PTLUST_S(KEEP(28))
      INTEGER STEP(N),
     & PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N + KEEP(253) ), FILS( N )
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER INTARR( max(1,KEEP(14)) )
      REAL DBLARR( max(1,KEEP(13)) )
      LOGICAL, intent(in) :: STACK_RIGHT_AUTHORIZED
       LOGICAL FLAG, RIGHT_MESS, FLAGbis
       INTEGER LP, MSGSOU_LOC, MSGTAG_LOC, MSGLEN_LOC
       INTEGER IERR
       INTEGER STATUS_BIS( MPI_STATUS_SIZE )
       INTEGER, SAVE :: RECURS = 0
      CALL SMUMPS_467(COMM_LOAD, KEEP)
      IF ( .NOT. STACK_RIGHT_AUTHORIZED ) THEN
          RETURN
      ENDIF
      RECURS = RECURS + 1
      LP = ICNTL(1)
      IF (ICNTL(4).LT.1) LP=-1
      IF ( MESSAGE_RECEIVED ) THEN
        MSGSOU_LOC = MPI_ANY_SOURCE
        MSGTAG_LOC = MPI_ANY_TAG
        GOTO 250
      ENDIF
      IF ( ASS_IRECV .NE. MPI_REQUEST_NULL) THEN
      RIGHT_MESS = .TRUE.
       IF (BLOCKING) THEN
         CALL MPI_WAIT(ASS_IRECV,
     &                STATUS, IERR)
         FLAG = .TRUE.
         IF ( ( (MSGSOU.NE.MPI_ANY_SOURCE) .OR.
     &      (MSGTAG.NE.MPI_ANY_TAG) )  ) THEN
           IF ( MSGSOU.NE.MPI_ANY_SOURCE) THEN
             RIGHT_MESS = MSGSOU.EQ.STATUS(MPI_SOURCE)
           ENDIF
           IF ( MSGTAG.NE.MPI_ANY_TAG) THEN
             RIGHT_MESS =
     &       ( (MSGTAG.EQ.STATUS(MPI_TAG)).AND.RIGHT_MESS )
           ENDIF
           IF (.NOT.RIGHT_MESS) THEN
             CALL MPI_PROBE(MSGSOU,MSGTAG,
     &           COMM, STATUS_BIS, IERR)
           ENDIF
         ENDIF
       ELSE
        CALL MPI_TEST(ASS_IRECV,
     &             FLAG, STATUS, IERR)
       ENDIF
       IF (IERR.LT.0) THEN
        IFLAG = -20
        IF (LP.GT.0)
     &  write(LP,*) ' Error return from MPI_TEST ',
     &     IFLAG, ' in SMUMPS_329'
        CALL SMUMPS_44( MYID, SLAVEF, COMM )
        RETURN
       ENDIF
       IF ( FLAG ) THEN
         MESSAGE_RECEIVED = .TRUE.
         MSGSOU_LOC = STATUS( MPI_SOURCE )
         MSGTAG_LOC = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN_LOC, IERR )
           IF (.NOT.RIGHT_MESS) RECURS = RECURS + 10
         CALL SMUMPS_322( COMM_LOAD, ASS_IRECV,
     &      MSGSOU_LOC, MSGTAG_LOC, MSGLEN_LOC, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &      IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS,
     &      PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
           IF (.NOT.RIGHT_MESS) RECURS = RECURS - 10
          IF ( IFLAG .LT. 0 ) RETURN
           IF (.NOT.RIGHT_MESS) THEN
              IF (ASS_IRECV .NE. MPI_REQUEST_NULL) THEN
                CALL MUMPS_ABORT()
              ENDIF
             CALL MPI_IPROBE(MSGSOU,MSGTAG,
     &           COMM, FLAGbis, STATUS, IERR)
             IF (FLAGbis) THEN
               MSGSOU_LOC = STATUS( MPI_SOURCE )
               MSGTAG_LOC = STATUS( MPI_TAG )
               CALL SMUMPS_280( COMM_LOAD, ASS_IRECV,
     &            STATUS, BUFR, LBUFR,
     &            LBUFR_BYTES,
     &            PROCNODE_STEPS, POSFAC,
     &            IWPOS, IWPOSCB, IPTRLU,
     &            LRLU, LRLUS, N, IW, LIW, A, LA,
     &            PTRIST, PTLUST_S, PTRFAC,
     &            PTRAST, STEP, PIMASTER, PAMASTER,
     &            NSTK_S, COMP, IFLAG,
     &            IERROR, COMM,
     &            NBPROCFILS,
     &            IPOOL, LPOOL, LEAF,
     &            NBFIN, MYID, SLAVEF,
     &
     &            root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &            FILS, PTRARW, PTRAIW,
     &            INTARR, DBLARR, ICNTL,
     &            KEEP,KEEP8, ND, FRERE,
     &            LPTRAR, NELT, FRTPTR, FRTELT,
     &            ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
                  IF ( IFLAG .LT. 0 ) RETURN
             ENDIF
           ENDIF
       ENDIF
      ELSE
         IF (BLOCKING) THEN
           CALL MPI_PROBE(MSGSOU,MSGTAG,
     &           COMM, STATUS, IERR)
           FLAG = .TRUE.
         ELSE
           CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,
     &           COMM, FLAG, STATUS, IERR)
         ENDIF
         IF (FLAG) THEN
          MSGSOU_LOC = STATUS( MPI_SOURCE )
          MSGTAG_LOC = STATUS( MPI_TAG )
          MESSAGE_RECEIVED = .TRUE.
          CALL SMUMPS_280( COMM_LOAD, ASS_IRECV,
     &      STATUS, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &      IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
          IF ( IFLAG .LT. 0 ) RETURN
         ENDIF
      ENDIF
 250  CONTINUE
      RECURS  = RECURS - 1
      IF ( NBFIN .EQ. 0 ) RETURN
      IF ( RECURS .GT. 3 ) RETURN
      IF ( KEEP(36).EQ.1 .AND. SET_IRECV  .AND.
     &      (ASS_IRECV.EQ.MPI_REQUEST_NULL) .AND.
     &    MESSAGE_RECEIVED ) THEN
       CALL MPI_IRECV ( BUFR(1),
     &      LBUFR_BYTES, MPI_PACKED, MPI_ANY_SOURCE,
     &      MPI_ANY_TAG, COMM,
     &      ASS_IRECV, IERR )
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_329
      SUBROUTINE SMUMPS_255( INFO1,
     &    ASS_IRECV,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    COMM,
     &    MYID, SLAVEF)
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER COMM
      INTEGER MYID, SLAVEF, INFO1, DEST
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL NO_ACTIVE_IRECV
      INTEGER MSGSOU_LOC, MSGTAG_LOC
      INTEGER IERR, DUMMY
      INTRINSIC mod
      IF (SLAVEF .EQ. 1) RETURN
      IF (ASS_IRECV.EQ.MPI_REQUEST_NULL) THEN
        NO_ACTIVE_IRECV=.TRUE.
      ELSE
        CALL MPI_TEST(ASS_IRECV, NO_ACTIVE_IRECV,
     &                STATUS, IERR)
      ENDIF
      CALL MPI_BARRIER(COMM,IERR)
      DUMMY = 1
      DEST = mod(MYID+1, SLAVEF)
      CALL SMUMPS_62
     &    (DUMMY, DEST, TAG_DUMMY, COMM, IERR)
      IF (NO_ACTIVE_IRECV) THEN
        CALL MPI_RECV( BUFR, LBUFR,
     &             MPI_INTEGER, MPI_ANY_SOURCE,
     &             TAG_DUMMY, COMM, STATUS, IERR )
      ELSE
        CALL MPI_WAIT(ASS_IRECV,
     &                STATUS, IERR)
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_255
      SUBROUTINE SMUMPS_180(
     &    INFO1, BUFR, LBUFR, LBUFR_BYTES,
     &    COMM_NODES, COMM_LOAD, SLAVEF, MP )
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER COMM_NODES, COMM_LOAD, SLAVEF, INFO1, MP
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL FLAG, BUFFERS_EMPTY, BUFFERS_EMPTY_ON_ALL_PROCS
      INTEGER MSGSOU_LOC, MSGTAG_LOC, COMM_EFF
      INTEGER IERR
      INTEGER IBUF_EMPTY, IBUF_EMPTY_ON_ALL_PROCS
      IF (SLAVEF.EQ.1) RETURN
      BUFFERS_EMPTY_ON_ALL_PROCS = .FALSE.
 10   CONTINUE
      FLAG = .TRUE.
      DO WHILE ( FLAG )
        COMM_EFF = COMM_NODES
        CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,
     &       COMM_NODES, FLAG, STATUS, IERR)
        IF ( .NOT. FLAG ) THEN
          COMM_EFF = COMM_LOAD
          CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         COMM_LOAD, FLAG, STATUS, IERR)
        END IF
        IF (FLAG) THEN
            MSGSOU_LOC = STATUS( MPI_SOURCE )
            MSGTAG_LOC = STATUS( MPI_TAG )
               CALL MPI_RECV( BUFR, LBUFR_BYTES,
     &             MPI_PACKED, MSGSOU_LOC,
     &             MSGTAG_LOC, COMM_EFF, STATUS, IERR )
           ENDIF
         END DO
        IF (BUFFERS_EMPTY_ON_ALL_PROCS) THEN
        RETURN
        ENDIF
        CALL SMUMPS_469(BUFFERS_EMPTY)
        IF ( BUFFERS_EMPTY ) THEN
          IBUF_EMPTY = 0
        ELSE
          IBUF_EMPTY = 1
        ENDIF
        CALL MPI_ALLREDUCE(IBUF_EMPTY,
     &                     IBUF_EMPTY_ON_ALL_PROCS,
     &                     1, MPI_INTEGER, MPI_MAX,
     &                     COMM_NODES, IERR)
        IF ( IBUF_EMPTY_ON_ALL_PROCS == 0) THEN
          BUFFERS_EMPTY_ON_ALL_PROCS = .TRUE.
        ELSE
          BUFFERS_EMPTY_ON_ALL_PROCS = .FALSE.
        ENDIF
        GOTO 10
      END SUBROUTINE SMUMPS_180
      INTEGER FUNCTION SMUMPS_748
     &     ( HBUF_SIZE, NNMAX, K227, K50 )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NNMAX, K227, K50
      INTEGER(8), INTENT(IN) :: HBUF_SIZE
      INTEGER K227_LOC
      INTEGER NBCOL_MAX
      INTEGER EFFECTIVE_SIZE
      NBCOL_MAX=int(HBUF_SIZE / int(NNMAX,8))
      K227_LOC = abs(K227)
      IF (K50.EQ.2) THEN
         K227_LOC=max(K227_LOC,2)
         EFFECTIVE_SIZE =  min(NBCOL_MAX-1, K227_LOC-1)
      ELSE
         EFFECTIVE_SIZE =  min(NBCOL_MAX, K227_LOC)
      ENDIF
      IF (EFFECTIVE_SIZE.LE.0) THEN   
         write(6,*) 'Internal buffers too small to store ', 
     &        ' ONE col/row of size', NNMAX
         CALL MUMPS_ABORT()
      ENDIF
      SMUMPS_748 = EFFECTIVE_SIZE
      RETURN
      END FUNCTION SMUMPS_748
      SUBROUTINE SMUMPS_698( IPIV, LPIV, ISHIFT,
     &     THE_PANEL, NBROW, NBCOL, KbeforePanel )
      IMPLICIT NONE
      INTEGER LPIV, ISHIFT, NBROW, NBCOL, KbeforePanel
      INTEGER IPIV(LPIV)
      REAL THE_PANEL(NBROW, NBCOL)
      INTEGER I, IPERM
      DO I = 1, LPIV
         IPERM=IPIV(I)
         IF ( I+ISHIFT.NE.IPERM) THEN
            CALL sswap(NBCOL,
     &           THE_PANEL(I+ISHIFT-KbeforePanel,1), NBROW,
     &           THE_PANEL(IPERM-KbeforePanel,1), NBROW)
         ENDIF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_698
      SUBROUTINE SMUMPS_667(TYPEF,
     &     NBPANELS,
     &     I_PIVPTR, I_PIV, IPOS, IW, LIW)
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER, intent(out) :: NBPANELS, I_PIVPTR, I_PIV
      INTEGER, intent(in) :: TYPEF 
      INTEGER, intent(in) :: LIW, IPOS
      INTEGER IW(LIW)
      INTEGER I_NBPANELS, I_NASS
      I_NASS       = IPOS
      I_NBPANELS   = I_NASS + 1 
      NBPANELS     = IW(I_NBPANELS) 
      I_PIVPTR     = I_NBPANELS + 1 
      I_PIV        = I_PIVPTR + NBPANELS 
      IF (TYPEF==TYPEF_U) THEN
         I_NBPANELS   = I_PIV+IW(I_NASS) 
         NBPANELS     = IW(I_NBPANELS) 
         I_PIVPTR     = I_NBPANELS + 1 
         I_PIV        = I_PIVPTR + NBPANELS 
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_667
      SUBROUTINE SMUMPS_691(K50,NBPANELS_L,NBPANELS_U,
     &     NASS, IPOS, IW, LIW )
      IMPLICIT NONE
      INTEGER K50
      INTEGER IPOS, NASS, NBPANELS_L, NBPANELS_U, LIW
      INTEGER IW(LIW)
      INTEGER IPOS_U
      IF (K50.EQ.1) THEN
         WRITE(*,*) "Internal error: SMUMPS_691 called"
      ENDIF
      IW(IPOS)=NASS
      IW(IPOS+1)=NBPANELS_L
      IW(IPOS+2:IPOS+1+NBPANELS_L)=NASS+1
      IF (K50 == 0) THEN
         IPOS_U=IPOS+2+NASS+NBPANELS_L
         IW(IPOS_U)=NBPANELS_U
         IW(IPOS_U+1:IPOS_U+NBPANELS_U)=NASS+1
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_691
      SUBROUTINE SMUMPS_644 ( 
     &     IWPOS, IOLDPS, IW, LIW, MonBloc, NFRONT, KEEP
     &     )
      USE SMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER, INTENT(IN)    :: IOLDPS, LIW, NFRONT,
     &     KEEP(500)
      INTEGER, INTENT(INOUT) :: IWPOS, IW(LIW)
      TYPE(IO_BLOCK), INTENT(IN):: MonBloc
      INTEGER :: NBPANELS_L,I_PIVRPTR_L, I_PIVR_L, NBPANELS_U, 
     &     I_PIVRPTR_U, I_PIVR_U, XSIZE, IBEGOOC
      LOGICAL FREESPACE         
      IF (KEEP(50).EQ.1) RETURN 
      IF ((IOLDPS+IW(IOLDPS+XXI)).NE.IWPOS) RETURN
      XSIZE   = KEEP(IXSZ)
      IBEGOOC = IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE
      CALL SMUMPS_667(TYPEF_L, NBPANELS_L, 
     &     I_PIVRPTR_L, I_PIVR_L, 
     &     IBEGOOC, IW, LIW)
      FREESPACE = 
     &     (MonBloc%LastPiv.EQ.(IW(I_PIVRPTR_L)-1))
      IF (KEEP(50).EQ.0) THEN
         CALL SMUMPS_667(TYPEF_U, NBPANELS_U, 
     &        I_PIVRPTR_U, I_PIVR_U, 
     &        IBEGOOC, IW, LIW)
         FREESPACE =  FREESPACE .AND.
     &        (MonBloc%LastPiv.EQ.(IW(I_PIVRPTR_U)-1)) 
      ENDIF
      IF (FREESPACE) THEN
         IW(IBEGOOC) = -7777    
         IW(IOLDPS+XXI) = IBEGOOC - IOLDPS + 1 
         IWPOS = IBEGOOC+1      
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_644
      SUBROUTINE SMUMPS_684(K50, NBROW_L, NBCOL_U, NASS,
     &     NBPANELS_L, NBPANELS_U, LREQ)
      USE SMUMPS_OOC       
      IMPLICIT NONE
      INTEGER, intent(IN)  :: K50, NBROW_L, NBCOL_U, NASS
      INTEGER, intent(OUT) :: NBPANELS_L, NBPANELS_U, LREQ
      NBPANELS_L=-99999
      NBPANELS_U=-99999
      IF (K50.EQ.1) THEN
         LREQ = 0
         RETURN
      ENDIF
      NBPANELS_L = (NASS / SMUMPS_690(NBROW_L))+1
      LREQ =    1               
     &     + 1                  
     &     + NASS               
     &     + NBPANELS_L         
      IF (K50.eq.0) THEN
         NBPANELS_U = (NASS / SMUMPS_690(NBCOL_U) ) +1
         LREQ = LREQ + 1        
     &        + NASS            
     &        + NBPANELS_U      
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_684
      SUBROUTINE SMUMPS_755
     &           (IW_LOCATION, MUST_BE_PERMUTED)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IW_LOCATION
      LOGICAL, INTENT(INOUT) :: MUST_BE_PERMUTED
      IF (IW_LOCATION .EQ. -7777) THEN
        MUST_BE_PERMUTED = .FALSE.
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_755
