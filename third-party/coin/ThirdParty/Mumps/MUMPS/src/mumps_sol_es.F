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
       MODULE MUMPS_SOL_ES
      PRIVATE
      PUBLIC:: PRUNED_SIZE_LOADED 
      PUBLIC:: MUMPS_797
      PUBLIC:: MUMPS_802
      PUBLIC:: MUMPS_798
      PUBLIC:: MUMPS_803
      PUBLIC:: MUMPS_804
      INTEGER(8), POINTER, DIMENSION(:,:) :: SIZE_OF_BLOCK
      INTEGER(8) :: PRUNED_SIZE_LOADED
      CONTAINS
      SUBROUTINE MUMPS_804(SIZE_OF_BLOCK_ARG, KEEP201)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KEEP201
      INTEGER(8), POINTER, DIMENSION(:,:) :: SIZE_OF_BLOCK_ARG
      IF (KEEP201 > 0) THEN
        SIZE_OF_BLOCK => SIZE_OF_BLOCK_ARG
      ELSE
        NULLIFY(SIZE_OF_BLOCK)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_804
      SUBROUTINE MUMPS_798( 
     &     fill,
     &     DAD, NE_STEPS, FRERE, KEEP28,
     &     FILS, STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: fill
      INTEGER, INTENT(IN) :: N, KEEP28
      INTEGER, INTENT(IN) :: DAD(KEEP28),NE_STEPS(KEEP28),FRERE(KEEP28)
      INTEGER, INTENT(IN) :: FILS(N), STEP(N)
      INTEGER, INTENT(IN) :: nodes_RHS(KEEP28),  nb_nodes_RHS
      INTEGER :: nb_prun_nodes
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_List(nb_prun_nodes) 
      INTEGER :: nb_prun_roots
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Roots(nb_prun_roots) 
      INTEGER :: nb_prun_leaves
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Leaves(nb_prun_leaves) 
      LOGICAL :: TO_PROCESS(KEEP28)  
      INTEGER :: IN, I, ISTEP, TMP, TMPsave
      nb_prun_nodes = 0
      nb_prun_leaves = 0
      TO_PROCESS(:) = .FALSE.
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         TMPsave = TMP
         ISTEP = STEP(TMP)
         DO WHILE(.NOT.TO_PROCESS(ISTEP))
            TO_PROCESS(ISTEP) = .TRUE.
            nb_prun_nodes = nb_prun_nodes + 1
            IF(fill) THEN
               Pruned_List(nb_prun_nodes) = TMP
            END IF
            IN = FILS(TMP) 
            DO WHILE(IN.GT.0) 
               IN = FILS(IN)
            END DO
            IF (IN.LT.0) THEN 
               TMP = -IN
               ISTEP = STEP(TMP)
            ELSE 
               nb_prun_leaves = nb_prun_leaves + 1
               IF(fill) THEN
                  Pruned_Leaves(nb_prun_leaves) = TMP
               END IF
               IF(TMP.NE.TMPsave) THEN
                  TMP = abs(FRERE(ISTEP))
                  IF(TMP.NE.0) THEN 
                     ISTEP = STEP(TMP) 
                  END IF
               END IF
            END IF
         END DO
      END DO
      nb_prun_roots = 0
      DO I=1,nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         IF(DAD(ISTEP).NE.0) THEN 
            IF(.NOT.TO_PROCESS(STEP(DAD(ISTEP)))) THEN
               nb_prun_roots = nb_prun_roots + 1
               IF(fill) THEN
                  Pruned_Roots(nb_prun_roots) = TMP
               END IF
            END IF
         ELSE 
            nb_prun_roots = nb_prun_roots + 1
            IF(fill) THEN
               Pruned_Roots(nb_prun_roots) = TMP
            END IF          
         END IF
      END DO
      RETURN
      END SUBROUTINE MUMPS_798
      SUBROUTINE MUMPS_797(
     &     fill,
     &     DAD, KEEP28,
     &     STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     Pruned_SONS, TO_PROCESS,
     &     nb_prun_nodes,nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: fill
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: STEP(N)
      INTEGER, INTENT(IN) :: KEEP28
      INTEGER, INTENT(IN) :: DAD(KEEP28)
      INTEGER, INTENT(IN) :: nb_nodes_RHS
      INTEGER, INTENT(IN) :: nodes_RHS(nb_nodes_RHS)
      INTEGER :: nb_prun_nodes
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_List(nb_prun_nodes) 
      INTEGER :: nb_prun_roots
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Roots(nb_prun_roots) 
      INTEGER :: nb_prun_leaves
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Leaves(nb_prun_leaves) 
      INTEGER :: Pruned_SONS(KEEP28) 
      LOGICAL :: TO_PROCESS(KEEP28)  
      INTEGER :: IN, I, ISTEP, TMP
      nb_prun_nodes = 0
      nb_prun_roots = 0
      TO_PROCESS(:) = .FALSE.
      Pruned_SONS(:) = -1
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         TO_PROCESS(ISTEP) = .TRUE.
         IF (Pruned_SONS(ISTEP) .eq. -1) THEN
            Pruned_SONS(ISTEP) = 0
            nb_prun_nodes = nb_prun_nodes + 1
            IF(fill) THEN
               Pruned_List(nb_prun_nodes) = nodes_RHS(I)
            END IF
            IN = nodes_RHS(I)
            IN = DAD(STEP(IN))
            DO WHILE (IN.NE.0)
               TO_PROCESS(STEP(IN)) = .TRUE.
               IF (Pruned_SONS(STEP(IN)).eq.-1) THEN 
                  nb_prun_nodes = nb_prun_nodes + 1
                  IF(fill) THEN
                     Pruned_List(nb_prun_nodes) = IN
                  END IF
                  Pruned_SONS(STEP(IN)) = 1
                  TMP = IN
                  IN = DAD(STEP(IN))
               ELSE 
                  Pruned_SONS(STEP(IN)) = Pruned_SONS(STEP(IN)) + 1
                  GOTO 201
               ENDIF
            ENDDO
            nb_prun_roots = nb_prun_roots +1
            IF(fill) THEN
               Pruned_Roots(nb_prun_roots) = TMP
            END IF
         ENDIF
  201    CONTINUE
      ENDDO
      nb_prun_leaves = 0
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         IF (Pruned_SONS(ISTEP).EQ.0) THEN
            nb_prun_leaves = nb_prun_leaves +1
            IF(fill) THEN
              Pruned_Leaves(nb_prun_leaves) = TMP
            END IF
         END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_797
      SUBROUTINE MUMPS_803(MYID, N, KEEP28, KEEP201,
     &           KEEP8_31,
     &           STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_LOC)
      INTEGER, intent(in) :: KEEP28, KEEP201, OOC_FCT_TYPE_LOC, MYID, N
      INTEGER(8), intent(in) :: KEEP8_31
      INTEGER, intent(in) :: nb_prun_nodes
      INTEGER, intent(in) :: Pruned_List(nb_prun_nodes)
      INTEGER, intent(in) :: STEP(N)
      INTEGER I, ISTEP
      INTEGER(8) :: Pruned_Size
#if defined(Mila_Print)
      write(*,*) ' in Pruned List nodes:',nb_prun_nodes
      write(*,*) Pruned_List(1:nb_prun_nodes)
#endif
      IF (KEEP201 .GT. 0) THEN
        Pruned_Size = 0_8
        DO I = 1, nb_prun_nodes
          ISTEP = STEP(Pruned_List(I))
          Pruned_Size = Pruned_Size + SIZE_OF_BLOCK
     &                  (ISTEP, OOC_FCT_TYPE_LOC)
        ENDDO
        PRUNED_SIZE_LOADED = PRUNED_SIZE_LOADED +Pruned_Size
#if defined(Mila_Print)
        write(*,*) 'Pruned_Size Total_Size:',
     &        Pruned_Size, KEEP8_31
        write(*,*) MYID,'Gain (%) = ', dble(100)
     &        - (dble(Pruned_Size)*dble(100)) /dble(KEEP8_31)
        IF (Pruned_Size .EQ. 0) THEN
            WRITE(*,*) "NOT NORMAL BEHAVIOUR !!"
            DO I = 1, nb_nodes_RHS
               WRITE(*,*) "starting_node node_size",
     &              nodes_RHS(I), 
     &              SIZE_OF_BLOCK(STEP(nodes_RHS(I)),OOC_FCT_TYPE_LOC)
            ENDDO
        ENDIF
        write(*,*) '============================='
#endif
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_803
      SUBROUTINE MUMPS_802
     &                (MYID, N, KEEP28, KEEP201, KEEP8_31,
     &                STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_LOC
     & )
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP28, KEEP201, OOC_FCT_TYPE_LOC, N
      INTEGER(8), intent(in) :: KEEP8_31
      INTEGER, intent(in) :: nb_prun_nodes, MYID
      INTEGER, intent(in) :: Pruned_List(nb_prun_nodes)
      INTEGER, intent(in) :: STEP(N)
      INCLUDE 'mpif.h'
      INTEGER I, ISTEP
      INTEGER(8) :: Pruned_Size
      Pruned_Size = 0_8
      DO I = 1, nb_prun_nodes
        ISTEP = STEP(Pruned_List(I))
        IF (KEEP201 .GT. 0) THEN
            Pruned_Size = Pruned_Size + SIZE_OF_BLOCK
     &                    (ISTEP, OOC_FCT_TYPE_LOC)
        ENDIF
      ENDDO
      IF (KEEP201.GT.0) THEN
#       if defined(Mila_Print)
          write(*,*) MYID,'PR leaves NODES',nb_prun_leaves,
     &     Pruned_Leaves(1:nb_prun_leaves)
          write(*,*) MYID,'PR  NODES',Pos_List,
     &     Pruned_List(1:Pos_List)
          write(*,*) 'PR root NODES',
     &     Pruned_Roots(nb_prun_roots)
#       endif
        IF (KEEP8_31 .NE. 0_8) THEN
#         if defined(Mila_Print)
            write(*,*) MYID,'PRUNED and TOTAL Size:',
     &        Pruned_Size, KEEP8_31
            write(*,*) MYID,'Gain (%) = ', dble(100)
     &        - ((dble(Pruned_Size)*dble(100))/dble(KEEP8_31))
            IF (MYID.EQ.0) 
     &        write(*,*) '============================='
#         endif
          PRUNED_SIZE_LOADED = PRUNED_SIZE_LOADED +Pruned_Size
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_802
      END MODULE MUMPS_SOL_ES
      SUBROUTINE MUMPS_780
     &          (PERM_STRAT, SYM_PERM,
     &           IRHS_PTR, NHRS,
     &           PERM_RHS, SIZEPERM, IERR
     &         )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: PERM_STRAT, NHRS, SIZEPERM
      INTEGER, INTENT(IN) :: SYM_PERM(SIZEPERM)
      INTEGER, INTENT(IN) :: IRHS_PTR(NHRS)
      INTEGER, INTENT(OUT):: IERR   
      INTEGER, INTENT(OUT):: PERM_RHS(SIZEPERM)
      DOUBLE PRECISION :: RAND_NUM
      INTEGER  I, J, STRAT
      IERR = 0
      STRAT = PERM_STRAT 
      IF( (STRAT.NE.-3).AND.
     &    (STRAT.NE.-2).AND.
     &    (STRAT.NE.-1).AND.
     &    (STRAT.NE. 1).AND.
     &    (STRAT.NE. 2).AND.
     &    (STRAT.NE. 6) ) THEN
        WRITE(*,*)"Warning: incorrect value for the RHS permutation; ",
     &            "defaulting to post-order"
        STRAT = 1
      END IF
      IF (STRAT .EQ. -3) THEN
         WRITE(*,*) "Processing the RHS in random order"
         PERM_RHS(1:SIZEPERM)=0
         DO I=1, SIZEPERM 
           CALL random_number(RAND_NUM) 
           RAND_NUM = RAND_NUM*dble(SIZEPERM) 
           J = ceiling(RAND_NUM) 
           DO WHILE (PERM_RHS(J).NE.0) 
             CALL random_number(RAND_NUM)
             RAND_NUM = RAND_NUM*dble(SIZEPERM)
             J = ceiling(RAND_NUM)
           ENDDO
           PERM_RHS(J)=I
         ENDDO
      ELSEIF (STRAT .EQ. -2) THEN
         WRITE(*,*) "Processing the RHS in inverse order"
         DO I=1, SIZEPERM
            PERM_RHS(SIZEPERM -I +1) = I
         ENDDO
      ELSEIF (STRAT .EQ. -1) THEN
         WRITE(*,*) "Processing the RHS in natural order"
         DO I=1, SIZEPERM
            PERM_RHS(I) = I
         ENDDO
      ELSEIF (STRAT .EQ.  1) THEN
         WRITE(*,*) "Processing the RHS in post-order"
         DO I=1, SIZEPERM
            PERM_RHS(SYM_PERM(I)) = I
         ENDDO
      ELSEIF (STRAT .EQ.  2) THEN
         WRITE(*,*) "Processing the RHS in pre-order"
         DO I=1, SIZEPERM
            PERM_RHS(SIZEPERM-SYM_PERM(I)+1) = I
         ENDDO
      ENDIF
      END SUBROUTINE MUMPS_780
      SUBROUTINE MUMPS_772
     &    (PERM_RHS, SIZEPERM, N, KEEP_28,
     &     PROCNODE, STEP_S, Nprocs, Step2node,
     &     IERR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: SIZEPERM
      INTEGER, intent(in) :: N, KEEP_28, Nprocs
      INTEGER, intent(in) :: PROCNODE(KEEP_28), STEP_S(N)
      INTEGER, intent(in) :: Step2node(KEEP_28)
      INTEGER, INTENT(OUT):: IERR
      INTEGER, INTENT(INOUT):: PERM_RHS(SIZEPERM)
      INTEGER I, TMP_RHS, TMP2, proc_num
      INTEGER , ALLOCATABLE :: TEMP_LOC_ARRAY(:)
      INTEGER PTR(0:Nprocs-1)
      INTEGER  MUMPS_275, MUMPS_330
      EXTERNAL MUMPS_275, MUMPS_330
      IERR = 0
      ALLOCATE(TEMP_LOC_ARRAY(SIZEPERM), stat=IERR)
      IF (IERR.GT.0) THEN
       WRITE(6,*) " Not enough memory to allocate working ",
     &            " arrays in MUMPS_772 "
       CALL MUMPS_ABORT()
      ENDIF
      proc_num = 0
      PTR(:) = 1
      DO I = 1, SIZEPERM
 555     CONTINUE
         IF ( PTR(proc_num).LE.SIZEPERM ) THEN
            TMP_RHS = PERM_RHS(PTR(proc_num)) 
            TMP2 = Step2node(abs (STEP_S(TMP_RHS) )) 
              IF (proc_num .EQ. MUMPS_275
     &           (PROCNODE(STEP_S(TMP2)),Nprocs)) THEN
               TEMP_LOC_ARRAY(I) = TMP_RHS
               PTR(proc_num) = PTR(proc_num)+1
               IF ( (MUMPS_330(PROCNODE(STEP_S(TMP2)),
     &              Nprocs).EQ.1)
     &             ) THEN
                  proc_num = mod(proc_num+1,Nprocs)
                  proc_num = mod(proc_num+1,Nprocs)
               ENDIF
            ELSE
               PTR(proc_num) = PTR(proc_num)+1
               GOTO 555
            ENDIF
         ELSE
            proc_num = mod(proc_num+1,Nprocs)
            GOTO 555
         ENDIF
      ENDDO
      WRITE(*,*) "Used interleaving of the RHS"
      DO I = 1, SIZEPERM
         PERM_RHS(I) = TEMP_LOC_ARRAY(I)
      ENDDO
      IF (allocated(TEMP_LOC_ARRAY)) DEALLOCATE (TEMP_LOC_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_772
