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
      SUBROUTINE DMUMPS_246(MYID, N, STEP, FRERE, FILS,
     &     NA, LNA, NE, DAD, ND, PROCNODE, SLAVEF,
     &     NRLADU, NIRADU, NIRNEC, NRLNEC,
     &     NRLNEC_ACTIVE, 
     &     NIRADU_OOC, NIRNEC_OOC,
     &     MAXFR, OPSA,
     &     KEEP,KEEP8, LOCAL_M, LOCAL_N, SBUF_RECOLD,
     &     SBUF_SEND, SBUF_REC, OPS_SUBTREE, NSTEPS,
     &     I_AM_CAND,NMB_PAR2, ISTEP_TO_INIV2, CANDIDATES, 
     &     IFLAG, IERROR
     &     ,MAX_FRONT_SURFACE_LOCAL
     &     ,MAX_SIZE_FACTOR, ENTRIES_IN_FACTORS_LOC
     &     ,ENTRIES_IN_FACTORS_LOC_MASTERS
     &     )
      IMPLICIT NONE
      INTEGER  MYID, N, LNA, IFLAG, IERROR
      INTEGER  NIRADU, NIRNEC
      INTEGER(8) NRLADU, NRLNEC, NRLNEC_ACTIVE
      INTEGER(8) NRLADU_CURRENT, NRLADU_ROOT_3
      INTEGER NIRADU_OOC, NIRNEC_OOC
      INTEGER MAXFR, NSTEPS
      INTEGER(8) MAX_FRONT_SURFACE_LOCAL
      INTEGER STEP(N)
      INTEGER FRERE(NSTEPS), FILS(N), NA(LNA), NE(NSTEPS),
     &        ND(NSTEPS), PROCNODE(NSTEPS), DAD(NSTEPS)
      INTEGER  SLAVEF, KEEP(500), LOCAL_M, LOCAL_N
      INTEGER(8) KEEP8(150)
      INTEGER(8) ENTRIES_IN_FACTORS_LOC,
     &          ENTRIES_IN_FACTORS_LOC_MASTERS
      INTEGER  SBUF_SEND, SBUF_REC
      INTEGER(8) SBUF_RECOLD
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      LOGICAL  I_AM_CAND(NMB_PAR2)
      INTEGER  CANDIDATES( SLAVEF+1, NMB_PAR2 )
      DOUBLE PRECISION OPSA
      DOUBLE PRECISION OPSA_LOC 
      INTEGER(8) MAX_SIZE_FACTOR
      DOUBLE PRECISION OPS_SUBTREE
      DOUBLE PRECISION OPS_SBTR_LOC 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TNSTK, IPOOL, LSTKI 
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR
      INTEGER(8) SBUFS_CB, SBUFR_CB
      INTEGER SBUFR, SBUFS
      INTEGER BLOCKING_RHS
      INTEGER ITOP,NELIM,NFR
      INTEGER(8) ISTKR, LSTK
      INTEGER ISTKI,  STKI, ISTKI_OOC
      INTEGER K,NSTK, IFATH
      INTEGER INODE, LEAF, NBROOT, IN
      INTEGER LEVEL, MAXITEMPCB
      INTEGER(8) CURRENT_ACTIVE_MEM, MAXTEMPCB
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF, INSSARBR
      INTEGER LEVELF, NCB, SIZECBI
      INTEGER(8) NCB8
      INTEGER(8) NFR8, NELIM8
      INTEGER(8) SIZECB, SIZECBINFR, SIZECB_SLAVE
      INTEGER SIZEHEADER, SIZEHEADER_OOC, XSIZE_OOC
      INTEGER EXTRA_PERM_INFO_OOC
      INTEGER NBROWMAX, NSLAVES_LOC, NSLAVES_PASSED,
     &         NELIMF, NFRF, NCBF,
     &         NBROWMAXF, LKJIB,
     &         LKJIBT, NBR, NBCOLFAC
      INTEGER(8) LEV3MAXREC, CBMAXR, CBMAXS
      INTEGER ALLOCOK
      INTEGER PANEL_SIZE
      LOGICAL COMPRESSCB
      DOUBLE PRECISION OPS_NODE, OPS_NODE_MASTER, OPS_NODE_SLAVE
      INTEGER(8) ENTRIES_NODE_UPPER_PART, ENTRIES_NODE_LOWER_PART
      INCLUDE 'mumps_headers.h'
      INTEGER WHAT
      INTEGER(8) IDUMMY8
      INTRINSIC min, int, real
      INTEGER DMUMPS_748
      EXTERNAL DMUMPS_748
      INTEGER MUMPS_275, MUMPS_330
      LOGICAL MUMPS_170
      INTEGER MUMPS_52
      EXTERNAL MUMPS_503, MUMPS_52
      EXTERNAL MUMPS_275, MUMPS_330, 
     &         MUMPS_170
      logical :: FORCE_CAND, CONCERNED, UPDATES, STACKCB, MASTERSON
      integer :: IFSON, LEVELSON
      IF (KEEP(50).eq.2) THEN
        EXTRA_PERM_INFO_OOC = 1
      ELSE IF (KEEP(50).eq.0) THEN
        EXTRA_PERM_INFO_OOC = 2
      ELSE
        EXTRA_PERM_INFO_OOC = 0
      ENDIF
      COMPRESSCB=( KEEP(215).EQ.0 .AND. KEEP(50).NE.0 )
      MAX_FRONT_SURFACE_LOCAL=0_8
      MAX_SIZE_FACTOR=0_8
      ALLOCATE( LSTKR(NSTEPS), TNSTK(NSTEPS), IPOOL(NSTEPS),
     &          LSTKI(NSTEPS) , stat=ALLOCOK)
      if (ALLOCOK .GT. 0) THEN
        IFLAG  =-7
        IERROR = 4*NSTEPS
        RETURN
      endif
      LKJIB = max(KEEP(5),KEEP(6))
      TNSTK = NE
      LEAF = NA(1)+1
      IPOOL(1:LEAF-1) = NA(3:3+LEAF-2)
      NBROOT = NA(2)
#if defined(OLD_OOC_NOPANEL)
      XSIZE_OOC=XSIZE_OOC_NOPANEL
#else
      IF (KEEP(50).EQ.0) THEN
              XSIZE_OOC=XSIZE_OOC_UNSYM
      ELSE
              XSIZE_OOC=XSIZE_OOC_SYM
      ENDIF
#endif
      SIZEHEADER_OOC = XSIZE_OOC+6 
      SIZEHEADER = XSIZE_IC + 6  
      ISTKR      = 0_8
      ISTKI      = 0
      ISTKI_OOC  = 0
      OPSA_LOC   = dble(0.0D0)
      ENTRIES_IN_FACTORS_LOC = 0_8
      ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
      OPS_SBTR_LOC = dble(0.0D0)
      NRLADU     = 0_8
      NIRADU     = 0
      NIRADU_OOC = 0
      NRLADU_CURRENT = 0_8
      NRLADU_ROOT_3 = 0_8
      NRLNEC_ACTIVE = 0_8
      NRLNEC     = 0_8
      NIRNEC     = 0
      NIRNEC_OOC = 0
      MAXFR      = 0
      ITOP       = 0
      MAXTEMPCB  = 0_8
      MAXITEMPCB = 0
      SBUFS_CB   = 1_8
      SBUFS      = 1
      SBUFR_CB   = 1_8
      SBUFR      = 1
      IF (KEEP(38) .NE. 0 .AND. KEEP(60).EQ.0) THEN
        INODE  = KEEP(38)
        NRLADU_ROOT_3 = int(LOCAL_M,8)*int(LOCAL_N,8)
        NRLADU = NRLADU_ROOT_3
        NRLNEC_ACTIVE = NRLADU_CURRENT
        MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_ROOT_3)
        NRLNEC = NRLADU
        IF (MUMPS_275(PROCNODE(STEP(INODE)),SLAVEF)
     &                                       .EQ. MYID) THEN
          NIRADU     = SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = SIZEHEADER_OOC+2*(ND(STEP(INODE))+KEEP(253))
        ELSE
          NIRADU     = SIZEHEADER
          NIRADU_OOC = SIZEHEADER_OOC
        ENDIF
        NIRNEC     = NIRADU
        NIRNEC_OOC = NIRADU_OOC
      ENDIF
      IF((KEEP(24).eq.0).OR.(KEEP(24).eq.1)) THEN
         FORCE_CAND=.FALSE.           
      ELSE
         FORCE_CAND=(mod(KEEP(24),2).eq.0)
      END IF
 90   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF - 1
         INODE = IPOOL(LEAF)
      ELSE 
         WRITE(MYID+6,*) ' ERROR 1 in DMUMPS_246 '
         CALL MUMPS_ABORT()
      ENDIF
 95   CONTINUE 
      NFR    = ND(STEP(INODE))+KEEP(253)
      NFR8   = int(NFR,8)
      NSTK   = NE(STEP(INODE))
      NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
      NELIM8=int(NELIM,8)
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
      IFSON = -IN
      IFATH = DAD(STEP(INODE))
      MASTER = MUMPS_275(PROCNODE(STEP(INODE)),SLAVEF)
     &           .EQ. MYID
      LEVEL  = MUMPS_330(PROCNODE(STEP(INODE)),SLAVEF)
      INSSARBR = MUMPS_170(PROCNODE(STEP(INODE)),
     &        SLAVEF)
      UPDATE=.FALSE.
       if(.NOT.FORCE_CAND) then
         UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
       else
         if(MASTER.and.(LEVEL.ne.3)) then
            UPDATE = .TRUE.
         else if(LEVEL.eq.2) then
            if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(INODE)))) THEN
              UPDATE = .TRUE.
            end if
         end if
       end if
      NCB      = NFR-NELIM
      NCB8     = int(NCB,8)
      SIZECBINFR = NCB8*NCB8
      IF (KEEP(50).EQ.0) THEN
        SIZECB = SIZECBINFR
      ELSE
        IFATH = DAD(STEP(INODE))
        IF ( IFATH.NE.KEEP(38) .AND. COMPRESSCB ) THEN
          SIZECB    = (NCB8*(NCB8+1_8))/2_8
        ELSE
          SIZECB    = SIZECBINFR
        ENDIF
      ENDIF
      SIZECBI      = 2* NCB  + SIZEHEADER
      IF (LEVEL.NE.2) THEN
        NSLAVES_LOC     = -99999999
        SIZECB_SLAVE = -99999997_8
        NBROWMAX        = NCB
      ELSE
        IF (KEEP(48) .EQ. 5) THEN
          WHAT = 5 
          IF (FORCE_CAND) THEN
            NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &                    ISTEP_TO_INIV2(STEP(INODE)))
          ELSE
            NSLAVES_LOC=SLAVEF-1
          ENDIF
          NSLAVES_PASSED=NSLAVES_LOC
        ELSE
          WHAT = 2 
          NSLAVES_PASSED=SLAVEF
          NSLAVES_LOC   =SLAVEF-1
        ENDIF
         CALL MUMPS_503(WHAT, KEEP,KEEP8,
     &     NCB, NFR, NSLAVES_PASSED, NBROWMAX, SIZECB_SLAVE
     &    )
      ENDIF
      IF (KEEP(60).GT.1) THEN
         IF (MASTER .AND. INODE.EQ.KEEP(38)) THEN
          NIRADU     = NIRADU+SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC+
     &                 2*(ND(STEP(INODE))+KEEP(253))
         ENDIF
      ENDIF
      IF (LEVEL.EQ.3) THEN
         IF ( 
     &     KEEP(60).LE.1
     &      ) THEN
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+
     &                 int(LOCAL_M,8)*int(LOCAL_N,8))
           NRLADU_CURRENT = int(LOCAL_M,8)*int(LOCAL_N,8)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_ROOT_3 + 
     &                        NRLADU_CURRENT+ISTKR)
         ENDIF
         IF (MASTER) THEN 
            IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
      ENDIF
      IF(KEEP(86).EQ.1)THEN
         IF(MASTER.AND.(.NOT.MUMPS_170(
     &        PROCNODE(STEP(INODE)), SLAVEF))
     &     )THEN
            IF(LEVEL.EQ.1)THEN
               MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &              NFR8*NFR8)
            ELSEIF(LEVEL.EQ.2)THEN
               IF(KEEP(50).EQ.0)THEN
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NFR8*NELIM8)
               ELSE
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NELIM8*NELIM8)
                 IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                  MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                  NELIM8*(NELIM8+1_8))
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF (LEVEL.EQ.2) THEN
        IF (MASTER) THEN
          IF (KEEP(50).EQ.0) THEN
             SBUFS = max(SBUFS, NFR*LKJIB+LKJIB+4)
          ELSE
             SBUFS = max(SBUFS, NELIM*LKJIB+NELIM+6)
          ENDIF
        ELSEIF (UPDATE) THEN
            if (KEEP(50).EQ.0) THEN
              SBUFR   = max(SBUFR, NFR*LKJIB+LKJIB+4)
            else
              SBUFR = max( SBUFR, NELIM*LKJIB+NELIM+6 )
              IF (KEEP(50).EQ.1) THEN
                LKJIBT  = LKJIB
              ELSE
                LKJIBT  = min( NELIM, LKJIB * 2 )
              ENDIF
              SBUFS = max(SBUFS,
     &                        LKJIBT*NBROWMAX+6)
              SBUFR = max( SBUFR, NBROWMAX*LKJIBT+6 )
            endif
        ENDIF
      ENDIF
      IF ( UPDATE ) THEN
          IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
            NIRADU     = NIRADU + 2*NFR + SIZEHEADER
            NIRADU_OOC = NIRADU_OOC + 2*NFR + SIZEHEADER_OOC
            PANEL_SIZE = DMUMPS_748(
     &      2_8*int(KEEP(226),8), NFR, KEEP(227), KEEP(50))
            NIRADU_OOC = NIRADU_OOC +
     &      EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
            IF (KEEP(50).EQ.0) THEN
             NRLADU_CURRENT = int(NELIM,8)*int(2*NFR-NELIM,8)
             NRLADU = NRLADU + NRLADU_CURRENT
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
            ELSE
             NRLADU_CURRENT = int(NELIM,8)*int(NFR,8)
             NRLADU = NRLADU + NRLADU_CURRENT
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
            ENDIF
            SIZECBI        = 2* NCB  + 6 + 3
          ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
              NIRADU     = NIRADU+SIZEHEADER +SLAVEF-1+2*NFR 
              NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC +SLAVEF-1+2*NFR 
              IF (KEEP(50).EQ.0) THEN
                NBCOLFAC=NFR
              ELSE
                NBCOLFAC=NELIM
              ENDIF
              PANEL_SIZE = DMUMPS_748(
     &        2_8*int(KEEP(226),8), NBCOLFAC, KEEP(227), KEEP(50))
              NIRADU_OOC = NIRADU_OOC +
     &        EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
              NRLADU_CURRENT = int(NBCOLFAC,8)*int(NELIM,8)
              NRLADU = NRLADU + NRLADU_CURRENT
              MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
               SIZECB     = 0_8
               SIZECBINFR = 0_8
               SIZECBI    = NCB + 5 +  SLAVEF - 1
            ELSE
             SIZECB=SIZECB_SLAVE
             SIZECBINFR = SIZECB
             NIRADU       = NIRADU+4+NELIM+NBROWMAX
             NIRADU_OOC   = NIRADU_OOC+4+NELIM+NBROWMAX
             IF (KEEP(50).EQ.0) THEN
               NRLADU   = NRLADU + int(NELIM,8)*int(NBROWMAX,8)
             ELSE 
               NRLADU   = NRLADU + int(NELIM,8)*int(NCB/NSLAVES_LOC,8)
             ENDIF
             NRLADU_CURRENT = int(NELIM,8)*int(NBROWMAX,8)
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
             SIZECBI       = 4 + NBROWMAX + NCB
             IF (KEEP(50).NE.0) THEN 
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_SYM
             ELSE
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_UNSYM 
             ENDIF
            ENDIF
         ENDIF
         NIRNEC = max0(NIRNEC,
     &             NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
         NIRNEC_OOC = max0(NIRNEC_OOC,
     &             NIRADU_OOC+ISTKI_OOC+SIZECBI+MAXITEMPCB +
     &             (XSIZE_OOC-XSIZE_IC) )  
         CURRENT_ACTIVE_MEM = ISTKR+SIZECBINFR
         IF (NSTK .NE. 0 .AND. INSSARBR .AND.
     &     KEEP(234).NE.0 .AND. KEEP(55).EQ.0) THEN
           CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTKR(ITOP)
         ENDIF
         IF (KEEP(50).NE.0.AND.UPDATE.AND.LEVEL.EQ.1) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM +
     &              int(NELIM,8)*int(NCB,8)
         ENDIF
         IF (MASTER .AND.  KEEP(219).NE.0.AND.
     &       KEEP(50).EQ.2.AND.LEVEL.EQ.2) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + int(NELIM,8)
         ENDIF
         IF (SLAVEF.EQ.1) THEN
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
         ELSE
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+MAXTEMPCB)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+MAXTEMPCB)
         ENDIF
         IF (NFR.GT.MAXFR) MAXFR = NFR
         IF (NSTK.GT.0) THEN
            DO 70 K=1,NSTK
               LSTK = LSTKR(ITOP)
               ISTKR = ISTKR - LSTK
               IF (K==1 .AND. INSSARBR.AND.KEEP(234).NE.0
     &            .AND.KEEP(55).EQ.0) THEN
               ELSE
                 CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTK
               ENDIF
               STKI = LSTKI( ITOP )
               ISTKI = ISTKI - STKI
               ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
               ITOP = ITOP - 1
               IF (ITOP.LT.0) THEN
                  write(*,*) MYID,
     &            ': ERROR 2 in DMUMPS_246. ITOP = ',ITOP
                  CALL MUMPS_ABORT()
               ENDIF
 70         CONTINUE
         ENDIF 
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = MUMPS_275(PROCNODE(STEP(IFSON)),SLAVEF)
     &                  .EQ.MYID
            LEVELSON  = MUMPS_330(PROCNODE(STEP(IFSON)),SLAVEF)
            if(.NOT.FORCE_CAND) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     &                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFSON)))) then
                    UPDATES = .TRUE.
                  end if
               end if
            end if
            IF (UPDATES) THEN
              LSTK = LSTKR(ITOP)
              ISTKR = ISTKR - LSTK
              STKI = LSTKI( ITOP )
              ISTKI = ISTKI - STKI
              ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
              ITOP = ITOP - 1
              IF (ITOP.LT.0) THEN
                write(*,*) MYID,
     &          ': ERROR 2 in DMUMPS_246. ITOP = ',ITOP
                CALL MUMPS_ABORT()
              ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) 
         END DO
      ENDIF
      IF (
     &        ( (INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0) ) 
     &       .AND.
     &        ( (INODE.NE.KEEP(38)).OR.(KEEP(60).LE.1) ) 
     &      )
     &THEN
            ENTRIES_NODE_LOWER_PART = int(NFR-NELIM,8) * int(NELIM,8)
            IF ( KEEP(50).EQ.0 ) THEN
              ENTRIES_NODE_UPPER_PART = int(NFR,8) * int(NELIM,8)
            ELSE
              ENTRIES_NODE_UPPER_PART =
     &        (int(NELIM,8)*int(NELIM+1,8))/2_8
            ENDIF
            IF (KEEP(50).EQ.2 .AND. LEVEL.EQ.3) THEN
              CALL MUMPS_511(NFR, 
     &           NELIM, NELIM,0,
     &           1,OPS_NODE)
            ELSE
              CALL MUMPS_511(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
            ENDIF
            IF (LEVEL.EQ.2) THEN
              CALL MUMPS_511(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           2,OPS_NODE_MASTER)
              OPS_NODE_SLAVE=OPS_NODE-OPS_NODE_MASTER
            ENDIF
      ELSE
           OPS_NODE = 0.0D0
           ENTRIES_NODE_UPPER_PART = 0_8
           ENTRIES_NODE_LOWER_PART = 0_8
      ENDIF
      IF ( MASTER ) 
     &  ENTRIES_IN_FACTORS_LOC_MASTERS = 
     &                     ENTRIES_IN_FACTORS_LOC_MASTERS +
     &                            ENTRIES_NODE_UPPER_PART +
     &                            ENTRIES_NODE_LOWER_PART
      IF (UPDATE.OR.LEVEL.EQ.3) THEN
         IF ( LEVEL .EQ. 3 ) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE / dble( SLAVEF )
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                            ENTRIES_NODE_UPPER_PART /
     &                            int(SLAVEF,8)
            IF (MASTER)
     &      ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                               mod(ENTRIES_NODE_UPPER_PART,
     &                                   int(SLAVEF,8))
         ELSE IF (MASTER .AND. LEVEL.EQ.2) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE_MASTER
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                      ENTRIES_NODE_UPPER_PART +
     &                      mod(ENTRIES_NODE_LOWER_PART,
     &                          int(NSLAVES_LOC,8))
         ELSE IF (MASTER .AND. LEVEL.EQ.1) THEN
            OPSA_LOC = OPSA_LOC + dble(OPS_NODE)
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                               ENTRIES_NODE_UPPER_PART +
     &                               ENTRIES_NODE_LOWER_PART
         ELSE IF (UPDATE) THEN 
            OPSA_LOC = OPSA_LOC + 
     &            dble(OPS_NODE_SLAVE)/dble(NSLAVES_LOC)
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC
     &                 + ENTRIES_NODE_LOWER_PART / 
     &                 int(NSLAVES_LOC,8)
         ENDIF
         IF (MUMPS_170(PROCNODE(STEP(INODE)),
     &        SLAVEF) .OR. NE(STEP(INODE))==0) THEN
           IF (LEVEL == 1) THEN
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ELSE
             CALL MUMPS_511(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ENDIF
         ENDIF
        ENDIF
      IF (IFATH .EQ. 0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
      ELSE
         NFRF = ND(STEP(IFATH))+KEEP(253)
         IF (DAD(STEP(IFATH)).EQ.0) THEN
           NELIMF = NFRF
         ELSE
           NELIMF = 0
           IN = IFATH
           DO WHILE (IN.GT.0)
              IN = FILS(IN)
              NELIMF = NELIMF+1
           ENDDO
         ENDIF
         NCBF = NFRF - NELIMF
         LEVELF = MUMPS_330(PROCNODE(STEP(IFATH)),SLAVEF)
         MASTERF= MUMPS_275(PROCNODE(STEP(IFATH)),SLAVEF).EQ.MYID
         UPDATEF= .FALSE.
         if(.NOT.FORCE_CAND) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(MASTERF.and.(LEVELF.ne.3)) then
               UPDATEF = .TRUE.
            else if (LEVELF.eq.2) then
               if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFATH)))) THEN
                 UPDATEF = .TRUE.
               end if
            end if
         end if
         CONCERNED  = UPDATEF .OR. UPDATE
         IF (LEVELF .NE. 2) THEN
           NBROWMAXF = -999999
         ELSE
           IF (KEEP(48) .EQ. 5) THEN
               WHAT = 4
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &               ISTEP_TO_INIV2(STEP(IFATH)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
           ELSE
               WHAT = 1 
               NSLAVES_LOC=SLAVEF
           ENDIF
           CALL MUMPS_503( WHAT, KEEP, KEEP8,
     &     NCBF, NFRF, NSLAVES_LOC, NBROWMAXF, IDUMMY8 
     &          )
         ENDIF
         IF(LEVEL.EQ.1.AND.UPDATE.AND.
     &      (UPDATEF.OR.LEVELF.EQ.2)
     &      .AND.LEVELF.NE.3) THEN
             IF ( INSSARBR .AND. KEEP(234).NE.0) THEN
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
             ELSE
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+SIZECB)
             ENDIF
         ENDIF
         IF (UPDATE .AND. LEVEL.EQ.2 .AND. .NOT. MASTER) THEN
             NRLNEC =
     &         max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+NRLADU_CURRENT)
             NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,2_8*NRLADU_CURRENT+
     &         NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
         ENDIF
        IF (LEVELF.EQ.3) THEN
          IF (LEVEL.EQ.1) THEN
            LEV3MAXREC = int(min(NCB,LOCAL_M),8) *
     &                   int(min(NCB,LOCAL_N),8)
          ELSE
            LEV3MAXREC = min(SIZECB,
     &                 int(min(NBROWMAX,LOCAL_M),8)
     &                *int(min(NCB,LOCAL_N),8))
          ENDIF
          MAXTEMPCB  = max(MAXTEMPCB, LEV3MAXREC)
          MAXITEMPCB = max(MAXITEMPCB,SIZECBI+SIZEHEADER)
          SBUFR_CB   = max(SBUFR_CB, LEV3MAXREC+int(SIZECBI,8))
          NIRNEC = max(NIRNEC,NIRADU+ISTKI+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
          NIRNEC_OOC = max(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
        ENDIF
        IF (CONCERNED) THEN
         IF (LEVELF.EQ.2) THEN
           IF (UPDATE.AND.(LEVEL.NE.2.OR..NOT.MASTER)) THEN
             IF(MASTERF)THEN
                 NBR = min(NBROWMAXF,NBROWMAX)
             ELSE
                 NBR = min(max(NELIMF,NBROWMAXF),NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
               CBMAXS = int(NBR,8)*int(NCB,8)
             ELSE
               CBMAXS = int(NBR,8)*int(NCB,8) -
     &                  (int(NBR,8)*int(NBR-1,8))/2_8
             ENDIF
           ELSE
              CBMAXS = 0_8
           END IF
           IF (MASTERF) THEN
             IF (LEVEL.EQ.1) THEN
                IF (.NOT.UPDATE) THEN
                  NBR = min(NELIMF, NCB)
                ELSE
                  NBR = 0
                ENDIF
             ELSE
                NBR = min(NELIMF, NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
                CBMAXR = int(NBR,8)*NCB8
             ELSE
                CBMAXR = int(NBR,8)*int(min(NCB,NELIMF),8)-
     &                   (int(NBR,8)*int(NBR-1,8))/2_8
                CBMAXR = min(CBMAXR, int(NELIMF,8)*int(NELIMF+1,8)/2_8)
                CBMAXR = min(CBMAXR, SIZECB)
                IF ((LEVEL.EQ.1).AND.(.NOT. COMPRESSCB)) THEN
                  CBMAXR = min(CBMAXR,(NCB8*(NCB8+1_8))/2_8)
                ENDIF
             ENDIF
           ELSE IF (UPDATEF) THEN
              NBR = min(NBROWMAXF,NBROWMAX)
              CBMAXR = int(NBR,8) * NCB8
              IF (KEEP(50).NE.0) THEN
                CBMAXR = CBMAXR - (int(NBR,8)*(int(NBR-1,8)))/2_8
              ENDIF
           ELSE
              CBMAXR = 0_8
           ENDIF
         ELSEIF (LEVELF.EQ.3) THEN
           CBMAXR = LEV3MAXREC
           IF (UPDATE.AND. .NOT. (MASTER.AND.LEVEL.EQ.2)) THEN
             CBMAXS = LEV3MAXREC
           ELSE
             CBMAXS = 0_8
           ENDIF
         ELSE
           IF (MASTERF) THEN
             CBMAXS = 0_8
             NBR = min(NFRF,NBROWMAX)
             IF ((LEVEL.EQ.1).AND.UPDATE) THEN
                NBR = 0
             ENDIF
             CBMAXR = int(NBR,8)*int(min(NFRF,NCB),8)
             IF (LEVEL.EQ.2)
     &       CBMAXR = min(CBMAXR, SIZECB_SLAVE)
             IF ( KEEP(50).NE.0 )  THEN
              CBMAXR = min(CBMAXR,(int(NFRF,8)*int(NFRF+1,8))/2_8)
             ELSE
              CBMAXR = min(CBMAXR,int(NFRF,8)*int(NFRF,8))
             ENDIF
           ELSE
             CBMAXR = 0_8
             CBMAXS = SIZECB
           ENDIF
         ENDIF
         IF (UPDATE) THEN
           CBMAXS = min(CBMAXS, SIZECB)
           IF ( .not. ( LEVELF .eq. 1 .AND. UPDATEF ) )THEN
              SBUFS_CB = max(SBUFS_CB, CBMAXS+int(SIZECBI,8))
           ENDIF
         ENDIF
         STACKCB = .FALSE.
         IF (UPDATEF) THEN
          STACKCB = .TRUE.
          SIZECBI     = 2 * NFR + SIZEHEADER
          IF (LEVEL.EQ.1) THEN
             IF (KEEP(50).NE.0.AND.LEVELF.NE.3
     &           .AND.COMPRESSCB) THEN
                 SIZECB = (NCB8*(NCB8+1_8))/2_8
             ELSE
                 SIZECB = NCB8*NCB8
             ENDIF
             IF (MASTER) THEN
               SIZECBI     = 2+ XSIZE_IC
             ELSE IF (LEVELF.EQ.1) THEN
               SIZECB  = min(CBMAXR,SIZECB)
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, int(SIZECBI,8)+SIZECB)
               SIZECBI    =  2 * NCB + SIZEHEADER     
             ELSE 
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, 
     &                      min(SIZECB,CBMAXR) + int(SIZECBI,8))
               MAXTEMPCB  = max(MAXTEMPCB, min(SIZECB,CBMAXR)) 
               SIZECBI    =  2 * NCB + SIZEHEADER 
               MAXITEMPCB = max(MAXITEMPCB, SIZECBI)
               SIZECBI     = 0
               SIZECB      = 0_8
             ENDIF
          ELSE 
             SIZECB = SIZECB_SLAVE
             MAXTEMPCB  = max(MAXTEMPCB, min(CBMAXR,SIZECB) )
             MAXITEMPCB = max(MAXITEMPCB,NBROWMAX+NCB+SIZEHEADER)
             IF (.NOT. 
     &        (UPDATE.AND.(.NOT.MASTER).AND.(NSLAVES_LOC.EQ.1))
     &          ) 
     &       SBUFR_CB = max(SBUFR_CB, 
     &            min(CBMAXR,SIZECB) + int(NBROWMAX + NCB + 6,8))
             IF (MASTER) THEN
              SIZECBI     =  NCB + 5 +  SLAVEF - 1 + XSIZE_IC
              SIZECB  = 0_8
             ELSE IF (UPDATE) THEN
              SIZECBI      =  NFR + 6 + SLAVEF - 1 + XSIZE_IC
              IF (KEEP(50).EQ.0) THEN
                SIZECBI = SIZECBI + NBROWMAX + NFR + 
     &                    SIZEHEADER
              ELSE
                SIZECBI = SIZECBI + NBROWMAX + NFR +
     &                    SIZEHEADER+ NSLAVES_LOC
              ENDIF
             ELSE
              SIZECB      = 0_8
              SIZECBI     = 0
             ENDIF
          ENDIF
         ELSE
           IF (LEVELF.NE.3) THEN
               STACKCB     = .TRUE.
               SIZECB      = 0_8
               SIZECBI     = 0
               IF ( (LEVEL.EQ.1) .AND. (LEVELF.NE.1) ) THEN
                  IF (COMPRESSCB) THEN 
                      SIZECB  = (NCB8*(NCB8+1_8))/2_8
                  ELSE
                      SIZECB  = NCB8*NCB8
                  ENDIF
                  SIZECBI     = 2 * NCB + SIZEHEADER
               ELSE IF (LEVEL.EQ.2) THEN
                 IF (MASTER) THEN
                   SIZECBI     =  NCB + 5 +  SLAVEF - 1 + XSIZE_IC
                 ELSE 
                   SIZECB  = SIZECB_SLAVE
                   SIZECBI = SIZECBI + NBROWMAX + NFR + SIZEHEADER
                 ENDIF 
               ENDIF
           ENDIF
         ENDIF
         IF (STACKCB) THEN
           IF (FRERE(STEP(INODE)).EQ.0) THEN
                  write(*,*) ' ERROR 3 in DMUMPS_246'
                  CALL MUMPS_ABORT()
           ENDIF
           ITOP = ITOP + 1
           IF ( ITOP .GT. NSTEPS ) THEN
             WRITE(*,*) 'ERROR 4 in DMUMPS_246 '
           ENDIF
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI + SIZECBI
           ISTKI_OOC = ISTKI_OOC + SIZECBI + (XSIZE_OOC-XSIZE_IC)
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
           NIRNEC = max0(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
           NIRNEC_OOC = max0(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &                  MAXITEMPCB + 
     &                    (XSIZE_OOC-XSIZE_IC) )  
         ENDIF 
        ENDIF 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
      ENDIF 
 115  CONTINUE
      BLOCKING_RHS = KEEP(84)
      IF (KEEP(84).EQ.0) BLOCKING_RHS=1
      NRLNEC = max(NRLNEC, 
     &         NRLADU+int(4*KEEP(127)*abs(BLOCKING_RHS),8))
      IF (BLOCKING_RHS .LT. 0) THEN
        BLOCKING_RHS = - 2 * BLOCKING_RHS
      ENDIF
      NRLNEC_ACTIVE = max(NRLNEC_ACTIVE, MAX_SIZE_FACTOR+
     &                    int(4*KEEP(127)*BLOCKING_RHS,8))
      SBUF_RECOLD = max(int(SBUFR,8),SBUFR_CB)
      SBUF_RECOLD = max(SBUF_RECOLD,
     &        MAXTEMPCB+int(MAXITEMPCB,8)) + 10_8
      SBUF_REC = max(SBUFR, int(min(100000_8,SBUFR_CB)))
      SBUF_REC = SBUF_REC   + 17
      SBUF_REC = SBUF_REC + 2 * KEEP(127) + SLAVEF - 1 + 7
      SBUF_SEND = max(SBUFS, int(min(100000_8,SBUFR_CB)))
      SBUF_SEND = SBUF_SEND + 17 
      IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
         SBUF_RECOLD = SBUF_RECOLD+int(KEEP(108)+1,8)
         SBUF_REC = SBUF_REC+KEEP(108)+1
         SBUF_SEND = SBUF_SEND+KEEP(108)+1
      ENDIF
      IF (SLAVEF.EQ.1) THEN 
         SBUF_RECOLD = 1_8
         SBUF_REC = 1
         SBUF_SEND= 1
      ENDIF
      DEALLOCATE( LSTKR, TNSTK, IPOOL,
     &          LSTKI )
      OPS_SUBTREE = dble(OPS_SBTR_LOC)
      OPSA        = dble(OPSA_LOC)
      KEEP(66)    = int(OPSA_LOC/1000000.d0)
      RETURN
      END SUBROUTINE DMUMPS_246
      RECURSIVE SUBROUTINE 
     &    DMUMPS_271( COMM_LOAD, ASS_IRECV, 
     &    INODE, NELIM_ROOT, root, 
     &
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
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER KEEP(500), ICNTL( 40 )
      INTEGER(8) KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER INODE, NELIM_ROOT
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: LA, POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS(KEEP(28))
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC( N + KEEP(253) ), FILS( N )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) ) 
      INTEGER INTARR(max(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mumps_tags.h'
      INTEGER I, LCONT, NCOL_TO_SEND, LDA
      INTEGER(8) :: SHIFT_VAL_SON, POSELT
      INTEGER FPERE, IOLDPS, NFRONT, NPIV, NASS, NSLAVES,
     &        H_INODE, NELIM, NBCOL, LIST_NELIM_ROW, 
     &        LIST_NELIM_COL, NELIM_LOCAL, TYPE_SON, 
     &        NROW, NCOL, NBROW, SHIFT_LIST_ROW_SON,
     &        SHIFT_LIST_COL_SON, LDAFS, IERR,
     &        STATUS( MPI_STATUS_SIZE ), ISON, PDEST_MASTER_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER MSGSOU, MSGTAG
      LOGICAL INVERT
      INCLUDE 'mumps_headers.h'
      INTEGER  MUMPS_275, MUMPS_330
      EXTERNAL MUMPS_275, MUMPS_330
      FPERE = KEEP(38)
      TYPE_SON = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF ( MUMPS_275( PROCNODE_STEPS(STEP(INODE)), 
     &     SLAVEF ).EQ.MYID) THEN
       IOLDPS   = PTLUST_S(STEP(INODE))
       NFRONT   = IW(IOLDPS+KEEP(IXSZ))
       NPIV     = IW(IOLDPS+1+KEEP(IXSZ))
       NASS     = iabs(IW(IOLDPS + 2+KEEP(IXSZ)))
       NSLAVES  =  IW(IOLDPS+5+KEEP(IXSZ))
       H_INODE  = 6 + NSLAVES + KEEP(IXSZ)
       NELIM    = NASS - NPIV
       NBCOL = NFRONT - NPIV
       LIST_NELIM_ROW = IOLDPS + H_INODE + NPIV
       LIST_NELIM_COL = LIST_NELIM_ROW + NFRONT
           IF (NELIM.LE.0) THEN
            write(6,*) ' ERROR 1 in DMUMPS_271 ', NELIM
            write(6,*) MYID,':Process root2son: INODE=',INODE,
     & 'Header=',IW(PTLUST_S(STEP(INODE)):PTLUST_S(STEP(INODE))
     &  +5+KEEP(IXSZ))
            CALL MUMPS_ABORT()
           ENDIF
       NELIM_LOCAL = NELIM_ROOT
       DO I=1, NELIM
        root%RG2L_ROW(IW(LIST_NELIM_ROW)) = NELIM_LOCAL
        root%RG2L_COL(IW(LIST_NELIM_COL)) = NELIM_LOCAL
        NELIM_LOCAL = NELIM_LOCAL + 1
        LIST_NELIM_ROW = LIST_NELIM_ROW + 1
        LIST_NELIM_COL = LIST_NELIM_COL + 1
       ENDDO
       NBROW = NFRONT - NPIV
       NROW = NELIM
       IF ( KEEP( 50 ) .eq. 0 ) THEN
         NCOL = NFRONT - NPIV
       ELSE
         NCOL = NELIM
       END IF
       SHIFT_LIST_ROW_SON = H_INODE + NPIV
       SHIFT_LIST_COL_SON = H_INODE + NFRONT + NPIV
       IF ( KEEP(50).eq.0 .OR. TYPE_SON .eq. 1 ) THEN
         LDAFS = NFRONT
       ELSE
         LDAFS = NASS
       END IF
       SHIFT_VAL_SON = int(NPIV,8) * int(LDAFS,8) + int(NPIV,8)
       CALL DMUMPS_80( COMM_LOAD,
     &   ASS_IRECV, 
     &   N, INODE, FPERE,
     &   PTLUST_S(1), PTRAST(1),
     &   root, NROW, NCOL, SHIFT_LIST_ROW_SON,
     &   SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDAFS,
     &   ROOT_NON_ELIM_CB, MYID, COMM,
     &   BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &   IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &   PTRIST, PTLUST_S(1), PTRFAC(1), PTRAST(1),
     &   STEP, PIMASTER, PAMASTER,
     &   NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &   IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &   OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &   FILS, PTRARW, PTRAIW,
     &   INTARR, DBLARR, ICNTL, KEEP,KEEP8, .FALSE., ND, FRERE,
     &   LPTRAR, NELT, FRTPTR, FRTELT, 
     &   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
       IF (IFLAG.LT.0 ) RETURN
       IF (TYPE_SON.EQ.1) THEN
        NROW = NFRONT - NASS
        NCOL = NELIM
        SHIFT_LIST_ROW_SON = H_INODE + NASS
        SHIFT_LIST_COL_SON = H_INODE + NFRONT + NPIV
        SHIFT_VAL_SON      = int(NASS,8) * int(NFRONT,8) + int(NPIV,8)
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          INVERT = .FALSE.
        ELSE
          INVERT = .TRUE.
        END IF
        CALL DMUMPS_80( COMM_LOAD, ASS_IRECV,
     &    N, INODE, FPERE,
     &    PTLUST_S, PTRAST,
     &    root, NROW, NCOL, SHIFT_LIST_ROW_SON,
     &    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, NFRONT,
     &    ROOT_NON_ELIM_CB, MYID, COMM,
     &
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, INVERT, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF (IFLAG.LT.0 ) RETURN
       ENDIF
       IOLDPS = PTLUST_S(STEP(INODE))
       POSELT = PTRAST(STEP(INODE))
       IW(IOLDPS + 4+KEEP(IXSZ)) = STEP(INODE)
       PTRFAC(STEP(INODE))=POSELT
       IF ( TYPE_SON .eq. 1 ) THEN
         NBROW = NFRONT - NPIV
       ELSE
         NBROW = NELIM
       END IF
       IF ( TYPE_SON .eq. 1 .OR. KEEP(50).EQ.0) THEN
         LDA = NFRONT
       ELSE
         LDA = NPIV+NBROW
       ENDIF
       CALL DMUMPS_324(A(POSELT), LDA,
     &          NPIV, NBROW, KEEP(50))
       IW(IOLDPS + KEEP(IXSZ))     = NBCOL
       IW(IOLDPS + 1 +KEEP(IXSZ)) = NASS - NPIV
       IF (TYPE_SON.EQ.2) THEN
        IW(IOLDPS + 2 +KEEP(IXSZ)) = NASS
       ELSE
        IW(IOLDPS + 2 +KEEP(IXSZ)) = NFRONT
       ENDIF
       IW(IOLDPS + 3 +KEEP(IXSZ)) = NPIV
      CALL DMUMPS_93(0_8,MYID,N,IOLDPS,TYPE_SON,IW,LIW,
     &    A, LA, POSFAC, LRLU, LRLUS,
     &    IWPOS, PTRAST,PTRFAC,STEP, KEEP,KEEP8, .FALSE.,INODE,IERR)
      IF(IERR.LT.0)THEN
         IFLAG=IERR
         IERROR=0
         RETURN
      ENDIF
      ELSE 
        ISON = INODE
        PDEST_MASTER_ISON = 
     &      MUMPS_275(PROCNODE_STEPS(STEP(ISON)), SLAVEF)
        DO WHILE ( PTRIST(STEP(ISON)) .EQ. 0)
          BLOCKING = .TRUE.
          SET_IRECV = .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    PDEST_MASTER_ISON, MAITRE_DESC_BANDE,
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
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
          IF ( IFLAG .LT. 0 ) RETURN
        ENDDO
        DO WHILE (
     &     ( IW( PTRIST(STEP(ISON)) + 1  +KEEP(IXSZ)) .NE.
     &       IW( PTRIST(STEP(ISON)) + 3  +KEEP(IXSZ)) ) .OR.
     &     ( KEEP(50) .NE. 0 .AND.
     &       IW( PTRIST(STEP(ISON)) + 6  +KEEP(IXSZ)) .NE. 0 ) )
          IF ( KEEP(50).eq.0) THEN
            MSGSOU = PDEST_MASTER_ISON
            MSGTAG = BLOC_FACTO
          ELSE
            IF ( IW( PTRIST(STEP(ISON)) + 1  +KEEP(IXSZ)) .NE.
     &           IW( PTRIST(STEP(ISON)) + 3  +KEEP(IXSZ)) ) THEN
              MSGSOU = PDEST_MASTER_ISON
              MSGTAG = BLOC_FACTO_SYM
            ELSE
              MSGSOU = MPI_ANY_SOURCE
              MSGTAG = BLOC_FACTO_SYM_SLAVE
            END IF
          END IF
          BLOCKING  = .TRUE.
          SET_IRECV = .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MSGSOU, MSGTAG,
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
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
          IF ( IFLAG .LT. 0 ) RETURN
        END DO
       IOLDPS = PTRIST(STEP(INODE))
       LCONT  = IW(IOLDPS+KEEP(IXSZ))
       NROW   = IW(IOLDPS+2+KEEP(IXSZ))
       NPIV   = IW(IOLDPS+3+KEEP(IXSZ))
       NASS   = IW(IOLDPS+4+KEEP(IXSZ))
       NELIM  = NASS-NPIV
       IF (NELIM.LE.0) THEN
         write(6,*) MYID,': INODE,LCONT, NROW, NPIV, NASS, NELIM=',
     &   INODE,LCONT, NROW, NPIV, NASS, NELIM
         write(6,*) MYID,': IOLDPS=',IOLDPS
         write(6,*) MYID,': ERROR 2 in DMUMPS_271 '
         CALL MUMPS_ABORT()
       ENDIF
       NSLAVES= IW(IOLDPS+5+KEEP(IXSZ))
       H_INODE = 6 + NSLAVES + KEEP(IXSZ)
       LIST_NELIM_COL = IOLDPS + H_INODE + NROW + NPIV
       NELIM_LOCAL = NELIM_ROOT
       DO I = 1, NELIM
        root%RG2L_COL(IW(LIST_NELIM_COL)) = NELIM_LOCAL
        root%RG2L_ROW(IW(LIST_NELIM_COL)) = NELIM_LOCAL
        NELIM_LOCAL = NELIM_LOCAL + 1
        LIST_NELIM_COL = LIST_NELIM_COL + 1
       ENDDO
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5+KEEP(IXSZ)) + KEEP(IXSZ)
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NPIV
       NCOL_TO_SEND       = NELIM
       IF (IW(IOLDPS+XXS).EQ.S_NOLCBNOCONTIG38.OR.
     &     IW(IOLDPS+XXS).EQ.S_ALL) THEN
         SHIFT_VAL_SON      = int(NPIV,8)
         LDA                = LCONT + NPIV
       ELSE IF (IW(IOLDPS+XXS).EQ.S_NOLCBCONTIG38) THEN
         SHIFT_VAL_SON = int(NROW,8)*int(LCONT+NPIV-NELIM,8)
         LDA           = NELIM
       ELSE IF (IW(IOLDPS+XXS).EQ.S_NOLCLEANED38) THEN
         SHIFT_VAL_SON=0_8
         LDA = NELIM
       ELSE
         write(*,*) MYID,": internal error in DMUMPS_271",
     &   IW(IOLDPS+XXS), "INODE=",INODE
         CALL MUMPS_ABORT()
       ENDIF
       IF ( KEEP( 50 ) .eq. 0 ) THEN
         INVERT = .FALSE.
       ELSE
         INVERT = .TRUE.
       END IF
       CALL DMUMPS_80( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE,
     &    PTRIST, PTRAST,
     &    root, NROW, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     &    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA,
     &    ROOT_NON_ELIM_CB, MYID, COMM,
     &
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8, INVERT, ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF (IFLAG.LT.0 ) RETURN
       IF (KEEP(214).EQ.2) THEN
        CALL DMUMPS_314( N, INODE,
     &      PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     &      LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     &      IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &      IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, TYPE_SON
     &      )
       ENDIF
        IF (IFLAG.LT.0) THEN
           CALL DMUMPS_44( MYID, SLAVEF, COMM )
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_271
      SUBROUTINE DMUMPS_221(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,
     &     DKEEP,PIVNUL_LIST,LPN_LIST,
     &
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U)
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,IFLAG,INOPV,NOFFW
      INTEGER(8) :: LA
      DOUBLE PRECISION A(LA) 
      DOUBLE PRECISION UU, SEUIL
      INTEGER IW(LIW)
      INTEGER(8) :: POSELT
      INTEGER  IOLDPS
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION    DKEEP(30)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION SWOP
      INTEGER XSIZE
      INTEGER(8) :: APOS, IDIAG
      INTEGER(8) :: J1, J2, J3, JJ
      INTEGER(8) :: NFRONT8
      DOUBLE PRECISION AMROW
      DOUBLE PRECISION RMAX
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER DMUMPS_IXAMAX
      INTRINSIC max
      DOUBLE PRECISION, PARAMETER :: RZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
        PIVNUL  = DKEEP(1)
        FIXA    = DKEEP(2)
        CSEUIL  = SEUIL
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        NFRONT8 = int(NFRONT,8)
        IF (KEEP(201).EQ.1) THEN
          CALL DMUMPS_667(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
          CALL DMUMPS_667(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
        NASSW   = iabs(IW(IOLDPS+3+XSIZE))
        IF(INOPV .EQ. -1) THEN
           APOS = POSELT + NFRONT8*int(NPIVP1-1,8) + int(NPIV,8)
           IDIAG = APOS
           IF(abs(A(APOS)).LT.SEUIL) THEN
              IF(dble(A(APOS)) .GE. RZERO) THEN
                 A(APOS) = CSEUIL
              ELSE
                 A(APOS) = -CSEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ELSE IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_762(A(APOS), DKEEP(6), KEEP(259))
           ENDIF
           IF (KEEP(201).EQ.1) THEN
             IF (KEEP(251).EQ.0) THEN 
               CALL DMUMPS_680( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
             ENDIF
             CALL DMUMPS_680( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
          DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + NFRONT8*int(IPIV-1,8) + int(NPIV,8)
            JMAX = 1
            IF (UU.GT.RZERO) GO TO 340
            IF (abs(A(APOS)).EQ.RZERO) GO TO 630
            GO TO 380
  340       AMROW = RZERO
            J1 = APOS
            J2 = APOS + int(- NPIV + NASS - 1,8)
             J     = NASS -NPIV
             JMAX  = DMUMPS_IXAMAX(J,A(J1),1)
             JJ    = J1 + int(JMAX - 1,8)
             AMROW = abs(A(JJ))
            RMAX = AMROW
            J1 = J2 + 1_8
            J2 = APOS +int(- NPIV + NFRONT - 1 - KEEP(253),8)
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = max(abs(A(JJ)),RMAX)
  360       CONTINUE
  370       IDIAG = APOS + int(IPIV - NPIVP1,8)
            IF ( RMAX .LE. PIVNUL ) THEN
               KEEP(109) = KEEP(109)+1
               ISW = IOLDPS+IW(IOLDPS+1+XSIZE)+6+XSIZE+
     &                      IW(IOLDPS+5+XSIZE)+IPIV-NPIVP1
               PIVNUL_LIST(KEEP(109)) = IW(ISW)
               IF(dble(FIXA).GT.RZERO) THEN
                  IF(dble(A(IDIAG)) .GE. RZERO) THEN
                     A(IDIAG) = FIXA
                  ELSE
                     A(IDIAG) = -FIXA
                  ENDIF
               ELSE
                 J1 = APOS
                 J2 = APOS + int(- NPIV + NFRONT - 1 - KEEP(253),8)
                 DO JJ=J1,J2
                   A(JJ) = ZERO
                 ENDDO
                 A(IDIAG) = -FIXA
               ENDIF
               JMAX = IPIV - NPIV
               GOTO 385   
            ENDIF
            IF (abs(A(IDIAG)).GT. max(UU*RMAX,SEUIL)) THEN
              JMAX = IPIV - NPIV
              GO TO 380
            ENDIF
            IF (AMROW.LE. max(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_762(
     &        A( APOS+int(JMAX-1,8) ),
     &        DKEEP(6),
     &        KEEP(259) )
            ENDIF
  385       CONTINUE
            IF (IPIV.EQ.NPIVP1) GO TO 400
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8)*NFRONT8
            J2 = J1 + NFRONT8 - 1_8
            J3 = POSELT + int(IPIV-1,8)*NFRONT8
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1_8
  390       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIVP1 + XSIZE
            ISWPS2 = IOLDPS + 5 + IPIV + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8)
            J2 = POSELT + int(NPIV + JMAX - 1,8)
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT8
              J2 = J2 + NFRONT8
  410       CONTINUE
            ISWPS1 = IOLDPS + 5 + NFRONT + NPIV + 1 +XSIZE
            ISWPS2 = IOLDPS + 5 + NFRONT + NPIV + JMAX +XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 430
  630 CONTINUE
      IFLAG = -10
      WRITE(*,*) 'Detected a null pivot, INODE/NPIV=',INODE,NPIV
      GOTO 430
  420 CONTINUE
              IF (KEEP(201).EQ.1) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_680( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, IPIV, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_680( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIV+JMAX, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
 430  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_221
      SUBROUTINE DMUMPS_220(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &   INOPV,NOFFW,IOLDPS,POSELT,UU,SEUIL,KEEP, DKEEP,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U)
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,INOPV
      INTEGER(8) :: LA
      INTEGER KEEP(500)
      DOUBLE PRECISION    DKEEP(30)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION AMROW
      DOUBLE PRECISION RMAX
      DOUBLE PRECISION  SWOP
      INTEGER(8) :: APOS, POSELT
      INTEGER(8) :: J1, J2, J3_8, JJ, IDIAG
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER NOFFW,NPIV,IPIV
      INTEGER J, J3
      INTEGER NPIVP1,JMAX,ISW,ISWPS1
      INTEGER ISWPS2,KSW,XSIZE
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INTEGER DMUMPS_IXAMAX
      INCLUDE 'mumps_headers.h'
      INTRINSIC max
      DOUBLE PRECISION, PARAMETER :: RZERO = 0.0D0
        NFRONT8 = int(NFRONT,8)
        INOPV   = 0
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
          CALL DMUMPS_667(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)
     &              +KEEP(IXSZ),
     &       IW, LIW)
          CALL DMUMPS_667(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
          DO 460 IPIV=NPIVP1,NASS
            APOS = POSELT + NFRONT8*int(NPIV,8) + int(IPIV-1,8)
            JMAX = 1
            AMROW = RZERO
            J1    = APOS
            J3    = NASS -NPIV
            JMAX  = DMUMPS_IXAMAX(J3,A(J1),NFRONT)
            JJ    = J1 + int(JMAX-1,8)*NFRONT8
            AMROW = abs(A(JJ))
            RMAX  = AMROW
            J1    = APOS +  int(NASS-NPIV,8) * NFRONT8
            J3 = NFRONT - NASS - KEEP(253)
            IF (J3.EQ.0) GOTO 370
            DO 360 J=1,J3
              RMAX = max(abs(A(J1)),RMAX)
              J1 = J1 + NFRONT8
  360       CONTINUE
  370       IF (RMAX.EQ.RZERO) GO TO 460
            IDIAG = APOS + int(IPIV - NPIVP1,8)*NFRONT8
            IF (abs(A(IDIAG)).GE.max(UU*RMAX,SEUIL)) THEN
              JMAX = IPIV - NPIV
              GO TO 380
            ENDIF
            IF (AMROW.LT.max(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_762(
     &             A(APOS + int(JMAX - 1,8) * NFRONT8 ),
     &             DKEEP(6),
     &             KEEP(259) )
            ENDIF
            IF (IPIV.EQ.NPIVP1) GO TO 400
            KEEP(260)=-KEEP(260)
            J1   = POSELT + int(NPIV,8)
            J3_8 = POSELT + int(IPIV-1,8)
            DO 390 J= 1,NFRONT
              SWOP  = A(J1)
              A(J1) = A(J3_8)
              A(J3_8) = SWOP
              J1 = J1 + NFRONT8
              J3_8 = J3_8 + NFRONT8
  390       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIVP1 + NFRONT + XSIZE
            ISWPS2 = IOLDPS + 5 + IPIV + NFRONT + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8) * NFRONT8
            J2 = POSELT + int(NPIV + JMAX - 1,8) * NFRONT8
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + 1_8
              J2 = J2 + 1_8
  410       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIV + 1 + XSIZE
            ISWPS2 = IOLDPS + 5 + NPIV + JMAX + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
       INOPV = 1
       GOTO 430
  420 CONTINUE
              IF (KEEP(201).EQ.1) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_680( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIV+JMAX,
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_680( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, IPIV,
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
 430  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_220
      SUBROUTINE DMUMPS_225(IBEG_BLOCK,
     &     NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &     IOLDPS,POSELT,IFINB,LKJIB,LKJIT,XSIZE)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,IFINB,LKJIB,IBEG_BLOCK
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER LKJIT, XSIZE
      DOUBLE PRECISION ONE, ALPHA
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,KROW,NEL
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8= int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
        IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
          IF (NASS.LT.LKJIT) THEN
           IW(IOLDPS+3+XSIZE) = NASS
          ELSE
           IW(IOLDPS+3+XSIZE) = min0(NASS,LKJIB)
          ENDIF
        ENDIF
        JROW2 = IW(IOLDPS+3+XSIZE)
        NEL2   = JROW2 - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
          IW(IOLDPS+3+XSIZE) = min0(JROW2+LKJIB,NASS)
          IBEG_BLOCK = NPIVP1+1
         ENDIF
        ELSE
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT8
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541     CONTINUE
         LPOS   = APOS + NFRONT8
         UUPOS  = APOS + 1_8
         CALL dger(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     &              A(LPOS+1_8),NFRONT)
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_225
      SUBROUTINE DMUMPS_229(NFRONT,N,INODE,IW,LIW,A,LA,IOLDPS,
     &          POSELT,XSIZE)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LIW,XSIZE
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA,VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS
      INTEGER(8) :: NFRONT8, LPOS, IRWPOS
      INTEGER IOLDPS,NPIV,NEL
      INTEGER JROW
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
        NFRONT8= int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + int(NPIV,8) * NFRONT8 + int(NPIV,8)
        IF (NEL.EQ.0) GO TO 650
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT8
        DO 340 JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT8
  340   CONTINUE
        LPOS   = APOS + NFRONT8
        UUPOS  = APOS+1_8
        DO 440 JROW = 1,NEL
             IRWPOS  = LPOS + 1_8
             ALPHA   = -A(LPOS)
             CALL daxpy(NEL,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT8
  440   CONTINUE
  650   RETURN
        END SUBROUTINE DMUMPS_229
      SUBROUTINE DMUMPS_228(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &       IOLDPS,POSELT,IFINB,XSIZE)
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER NFRONT,NASS,N,LIW,INODE,IFINB
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA,VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS, IRWPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS,NPIV,KROW, XSIZE
      INTEGER NEL,ICOL,NEL2
      INTEGER NPIVP1
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
        NFRONT8=int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        NEL2   = NASS - NPIVP1
        IFINB  = 0
        IF (NPIVP1.EQ.NASS) IFINB = 1
        APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT8
        DO 541 KROW = 1,NEL
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541    CONTINUE
        LPOS   = APOS + NFRONT8
        UUPOS  = APOS + 1_8
        DO 440 ICOL = 1,NEL
             IRWPOS  = LPOS + 1_8
             ALPHA   = -A(LPOS)
             CALL daxpy(NEL2,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT8
  440   CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_228
      SUBROUTINE DMUMPS_231(A,LA,NFRONT,
     &       NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER(8) :: LA,POSELT
      DOUBLE PRECISION    A(LA)
      INTEGER NFRONT, NPIV, NASS
      INTEGER(8) :: LPOS, LPOS1, LPOS2
      INTEGER NEL1,NEL11
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT + int(NASS,8)*int(NFRONT,8)
        CALL dtrsm('L','L','N','N',NPIV,NEL1,ONE,A(POSELT),NFRONT,
     &              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + int(NPIV,8)
        LPOS1  = POSELT + int(NPIV,8)
        CALL dgemm('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     &          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_231
      SUBROUTINE DMUMPS_642(A,LAFAC,NFRONT,
     &      NPIV,NASS, IW, LIWFAC,
     &      MonBloc, TYPEFile, MYID, KEEP8,
     &      STRAT, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten)
      USE DMUMPS_OOC   
      IMPLICIT NONE
      INTEGER NFRONT, NPIV, NASS
      INTEGER(8) :: LAFAC
      INTEGER  LIWFAC, TYPEFile, MYID, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten, STRAT
      DOUBLE PRECISION  A(LAFAC)
      INTEGER  IW(LIWFAC)
      INTEGER(8) KEEP8(150)
      TYPE(IO_BLOCK) :: MonBloc 
      INTEGER(8) :: LPOS2,LPOS1,LPOS
      INTEGER NEL1,NEL11
      DOUBLE PRECISION ALPHA, ONE
      LOGICAL LAST_CALL
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = 1_8 + int(NASS,8) * int(NFRONT,8)
        CALL dtrsm('L','L','N','N',NPIV,NEL1,ONE,A(1),NFRONT,
     &              A(LPOS2),NFRONT)
        LAST_CALL=.FALSE.
           CALL DMUMPS_688
     &          ( STRAT, TYPEFile, 
     &           A, LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW, LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
        LPOS   = LPOS2 + int(NPIV,8)
        LPOS1  = int(1 + NPIV,8)
        CALL dgemm('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     &          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_642
      SUBROUTINE DMUMPS_232(A,LA,NFRONT,NPIV,NASS,POSELT,LKJIB)
      INTEGER NFRONT, NPIV, NASS, LKJIB
      INTEGER (8) :: POSELT, LA
      DOUBLE PRECISION    A(LA)
      INTEGER(8) :: POSELT_LOCAL, LPOS, LPOS1, LPOS2
      INTEGER NEL1, NEL11, NPBEG
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        POSELT_LOCAL = POSELT
        NEL1   = NASS - NPIV
        NPBEG  = NPIV - LKJIB + 1
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT_LOCAL + int(NPIV,8)*int(NFRONT,8)
     &                        + int(NPBEG - 1,8)
        POSELT_LOCAL = POSELT_LOCAL + int(NPBEG-1,8)*int(NFRONT,8)
     &                              + int(NPBEG-1,8)
        CALL dtrsm('L','L','N','N',LKJIB,NEL1,ONE,A(POSELT_LOCAL),
     &               NFRONT,A(LPOS2),NFRONT)
        LPOS   = LPOS2 + int(LKJIB,8)
        LPOS1  = POSELT_LOCAL + int(LKJIB,8)
        CALL dgemm('N','N',NEL11,NEL1,LKJIB,ALPHA,A(LPOS1),
     &       NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_232
      SUBROUTINE DMUMPS_233(IBEG_BLOCK,
     &    NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &    IOLDPS,POSELT,LKJIB_ORIG,LKJIB,LKJIT,XSIZE )
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,LIW
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB_ORIG,LKJIB, INODE, IBEG_BLOCK
      INTEGER(8) :: POSELT, LPOS, LPOS1, LPOS2, POSLOCAL
      INTEGER(8) :: IPOS, KPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1, NEL11
      INTEGER LBP, HF
      INTEGER LBPT,I1,K1,II,ISWOP,LBP1
      INTEGER LKJIT, XSIZE
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8=int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        JROW2  = iabs(IW(IOLDPS+3+XSIZE))
        NPBEG  = IBEG_BLOCK
        HF     = 6 + IW(IOLDPS+5+XSIZE) +XSIZE
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
         LKJIB       = LKJIB_ORIG + NONEL
         IW(IOLDPS+3+XSIZE)= min0(NPIV+LKJIB,NASS)
        ELSE
          IW(IOLDPS+3+XSIZE) = NASS
        ENDIF
        IBEG_BLOCK = NPIV + 1
        NEL1   = NASS - JROW2
        LKJIW  = NPIV - NPBEG + 1
        NEL11  = NFRONT - NPIV
        IF ((NEL1.NE.0).AND.(LKJIW.NE.0)) THEN
          LPOS2  = POSELT + int(JROW2,8)*NFRONT8 +
     &             int(NPBEG - 1,8)
          POSLOCAL = POSELT + int(NPBEG-1,8)*NFRONT8 + int(NPBEG - 1,8)
          CALL dtrsm('L','L','N','N',LKJIW,NEL1,ONE,
     &               A(POSLOCAL),NFRONT,
     &               A(LPOS2),NFRONT)
          LPOS   = LPOS2 + int(LKJIW,8)
          LPOS1  = POSLOCAL + int(LKJIW,8)
          CALL dgemm('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     &          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_233
      SUBROUTINE DMUMPS_236(A,LA,NPIVB,NFRONT,
     &                             NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER NPIVB,NASS
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER(8) :: APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER(8) :: LPOS, LPOS1, LPOS2
      INTEGER NEL1, NEL11, NPIVE
      DOUBLE PRECISION    ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + int(NPIVB,8)*int(NFRONT,8)
     &                  + int(NPIVB,8)
        LPOS2  = APOS + int(NASSL,8)
        CALL dtrsm('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     &              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + int(NFRONT,8)*int(NPIVE,8)
        LPOS1  = APOS  + int(NFRONT,8)*int(NPIVE,8)
        CALL dgemm('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     &          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_236
       SUBROUTINE DMUMPS_217(N, NZ, NSCA, 
     &      ASPK, IRN, ICN, COLSCA, ROWSCA, WK, LWK, WK_REAL,
     &      LWK_REAL, ICNTL, INFO)
       IMPLICIT NONE
      INTEGER N, NZ, NSCA
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER ICNTL(40), INFO(40)
      DOUBLE PRECISION    ASPK(NZ)
      DOUBLE PRECISION COLSCA(*), ROWSCA(*)
      INTEGER LWK, LWK_REAL
      DOUBLE PRECISION    WK(LWK)
      DOUBLE PRECISION WK_REAL(LWK_REAL)
      INTEGER MPG,LP
      INTEGER IWNOR
      INTEGER I, K
      LOGICAL PROK
      DOUBLE PRECISION ONE
      PARAMETER( ONE = 1.0D0 )
      LP      = ICNTL(1)
      MPG     = ICNTL(2)
      MPG    = ICNTL(3)
      PROK   = (MPG.GT.0)
      IF (PROK) WRITE(MPG,101)
 101    FORMAT(/' ****** SCALING OF ORIGINAL MATRIX '/)
        IF (NSCA.EQ.1) THEN
         IF (PROK)
     &    WRITE (MPG,*) ' DIAGONAL SCALING '
        ELSEIF (NSCA.EQ.2) THEN
         IF (PROK)
     &   WRITE (MPG,*) ' SCALING BASED ON (MC29)'
        ELSEIF (NSCA.EQ.3) THEN
         IF (PROK)
     &   WRITE (MPG,*) ' COLUMN SCALING'
        ELSEIF (NSCA.EQ.4) THEN
         IF (PROK)
     &   WRITE (MPG,*) ' ROW AND COLUMN SCALING (1 Pass)'
        ELSEIF (NSCA.EQ.5) THEN
         IF (PROK)
     &   WRITE (MPG,*) ' MC29 FOLLOWED BY ROW &COL SCALING'
        ELSEIF (NSCA.EQ.6) THEN
         IF (PROK)
     &   WRITE (MPG,*) ' MC29 FOLLOWED BY COLUMN SCALING'
        ENDIF
        DO 10 I=1,N
            COLSCA(I) = ONE
            ROWSCA(I) = ONE
 10     CONTINUE
        IF ((NSCA.EQ.5).OR.
     &      (NSCA.EQ.6))                   THEN
          IF (NZ.GT.LWK) GOTO 400
          DO 15 K=1,NZ
           WK(K) = ASPK(K)
  15      CONTINUE
        ENDIF
        IF (5*N.GT.LWK_REAL) GOTO 410
        IWNOR = 1
          IF (NSCA.EQ.1) THEN
            CALL DMUMPS_238(N,NZ,ASPK,IRN,ICN,
     &        COLSCA,ROWSCA,MPG)
          ELSEIF (NSCA.EQ.2) THEN
            CALL DMUMPS_239(N,NZ,ASPK,IRN,ICN,
     &      ROWSCA,COLSCA,WK_REAL(IWNOR),MPG,MPG,NSCA)
          ELSEIF (NSCA.EQ.3) THEN
            CALL DMUMPS_241(N,NZ,ASPK,IRN,ICN,WK_REAL(IWNOR),
     &      COLSCA, MPG)
          ELSEIF (NSCA.EQ.4) THEN
            CALL DMUMPS_287(N,NZ,IRN,ICN,ASPK,
     &      WK_REAL(IWNOR),WK_REAL(IWNOR+N),COLSCA,ROWSCA,MPG)
          ELSEIF (NSCA.EQ.5) THEN
            CALL DMUMPS_239(N,NZ,WK,IRN,ICN,
     &           ROWSCA,COLSCA,WK_REAL(IWNOR),MPG,MPG,NSCA)
            CALL DMUMPS_241(N,NZ,WK,IRN,ICN,WK_REAL(IWNOR),
     &           COLSCA, MPG)
          ELSEIF (NSCA.EQ.6) THEN
            CALL DMUMPS_239(N,NZ,WK,IRN,ICN,
     &           ROWSCA,COLSCA,WK_REAL(IWNOR),MPG,MPG,NSCA)
            CALL DMUMPS_240(NSCA,N,NZ,IRN,ICN,WK,
     &           WK_REAL(IWNOR+N),ROWSCA,MPG)
            CALL DMUMPS_241(N,NZ,WK,IRN,ICN,
     &           WK_REAL(IWNOR), COLSCA, MPG)
          ENDIF
      GOTO 500
 400  INFO(1) = -5
      INFO(2) = NZ-LWK
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1))
     & WRITE(LP,*) '*** ERROR: Not enough space to scale matrix'
      GOTO 500
 410  INFO(1) = -5
      INFO(2) = 5*N-LWK_REAL
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1))
     & WRITE(LP,*) '*** ERROR: Not enough space to scale matrix'
      GOTO 500
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_217
      SUBROUTINE DMUMPS_287(N,NZ,IRN,ICN,VAL,
     &    RNOR,CNOR,COLSCA,ROWSCA,MPRINT)
      INTEGER N, NZ
      DOUBLE PRECISION    VAL(NZ)
      DOUBLE PRECISION    RNOR(N),CNOR(N)
      DOUBLE PRECISION    COLSCA(N),ROWSCA(N)
      DOUBLE PRECISION    CMIN,CMAX,RMIN,ARNOR,ACNOR
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION    VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO, ONE
      PARAMETER(ZERO=0.0D0, ONE=1.0D0)
      DO 50 J=1,N
       CNOR(J)   = ZERO
       RNOR(J)   = ZERO
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     &        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = abs(VAL(K))
            IF (VDIAG.GT.CNOR(J)) THEN
              CNOR(J) =     VDIAG
            ENDIF
            IF (VDIAG.GT.RNOR(I)) THEN
              RNOR(I) =     VDIAG
            ENDIF
 100   CONTINUE
      IF (MPRINT.GT.0) THEN
       CMIN = CNOR(1)
       CMAX = CNOR(1)
       RMIN = RNOR(1)
       DO 111 I=1,N
        ARNOR = RNOR(I)
        ACNOR = CNOR(I)
        IF (ACNOR.GT.CMAX) CMAX=ACNOR
        IF (ACNOR.LT.CMIN) CMIN=ACNOR
        IF (ARNOR.LT.RMIN) RMIN=ARNOR
 111   CONTINUE
       WRITE(MPRINT,*) '**** STAT. OF MATRIX PRIOR ROW&COL SCALING'
       WRITE(MPRINT,*) ' MAXIMUM NORM-MAX OF COLUMNS:',CMAX
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF COLUMNS:',CMIN
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF ROWS   :',RMIN
      ENDIF
      DO 120 J=1,N
       IF (CNOR(J).LE.ZERO) THEN
         CNOR(J)   = ONE
       ELSE
         CNOR(J)   = ONE / CNOR(J)
       ENDIF
 120  CONTINUE
      DO 130 J=1,N
       IF (RNOR(J).LE.ZERO) THEN
         RNOR(J)   = ONE
       ELSE
         RNOR(J)   = ONE / RNOR(J)
       ENDIF
 130  CONTINUE
       DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I) * RNOR(I)
        COLSCA(I) = COLSCA(I) * CNOR(I)
 110   CONTINUE
      IF (MPRINT.GT.0)
     &  WRITE(MPRINT,*) ' END OF SCALING BY MAX IN ROW AND COL'
      RETURN
      END SUBROUTINE DMUMPS_287
      SUBROUTINE DMUMPS_239(N,NZ,VAL,ROWIND,COLIND,
     &       RNOR,CNOR,WNOR,MPRINT,MP,
     &       NSCA)
      INTEGER N, NZ
      DOUBLE PRECISION    VAL(NZ)
      DOUBLE PRECISION WNOR(5*N)
      DOUBLE PRECISION RNOR(N), CNOR(N)
      INTEGER COLIND(NZ),ROWIND(NZ)
      INTEGER J,I,K
      INTEGER MPRINT,MP,NSCA
      INTEGER IFAIL9
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0)
      DO 15 I=1,N
       RNOR(I) = ZERO
       CNOR(I) = ZERO
  15  CONTINUE
      CALL DMUMPS_216(N,N,NZ,VAL,ROWIND,COLIND,
     &     RNOR,CNOR,WNOR, MP,IFAIL9)
*CVD$ NODEPCHK
*CVD$ VECTOR
*CVD$ CONCUR
      DO 30 I=1,N
       CNOR(I) = exp(CNOR(I))
       RNOR(I) = exp(RNOR(I))
  30  CONTINUE
      IF ((NSCA.EQ.5).OR.(NSCA.EQ.6)) THEN
        DO 100 K=1,NZ
          I   = ROWIND(K)
          J   = COLIND(K)
          IF (min(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 100
          VAL(K) = VAL(K) * CNOR(J) * RNOR(I)
 100    CONTINUE
      ENDIF
      IF (MPRINT.GT.0) 
     &   WRITE(MPRINT,*) ' END OF SCALING USING MC29'
      RETURN
      END SUBROUTINE DMUMPS_239
      SUBROUTINE DMUMPS_241(N,NZ,VAL,IRN,ICN,
     &       CNOR,COLSCA,MPRINT)
      INTEGER N,NZ
      DOUBLE PRECISION VAL(NZ)
      DOUBLE PRECISION CNOR(N)
      DOUBLE PRECISION COLSCA(N)
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DO 10 J=1,N
       CNOR(J)   = ZERO
  10  CONTINUE
      DO 100 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.LE.0).OR.(I.GT.N).OR.
     &      (J.LE.0).OR.(J.GT.N)) GOTO 100
        VDIAG = abs(VAL(K))
        IF (VDIAG.GT.CNOR(J)) THEN
           CNOR(J) =     VDIAG
        ENDIF
 100  CONTINUE
      DO 110 J=1,N
       IF (CNOR(J).LE.ZERO) THEN
         CNOR(J)   = ONE
       ELSE
         CNOR(J)   = ONE/CNOR(J)
       ENDIF
 110  CONTINUE
       DO 215 I=1,N
        COLSCA(I) = COLSCA(I) * CNOR(I)
 215   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF COLUMN SCALING'
      RETURN
      END SUBROUTINE DMUMPS_241
      SUBROUTINE DMUMPS_238(N,NZ,VAL,IRN,ICN,
     &      COLSCA,ROWSCA,MPRINT)
      INTEGER   N, NZ
      DOUBLE PRECISION  VAL(NZ)
      DOUBLE PRECISION ROWSCA(N),COLSCA(N)
      INTEGER   IRN(NZ),ICN(NZ)
      DOUBLE PRECISION      VDIAG
      INTEGER   MPRINT,I,J,K
      INTRINSIC sqrt
      DOUBLE PRECISION ZERO, ONE
      PARAMETER(ZERO=0.0D0, ONE=1.0D0)
      DO 10 I=1,N
       ROWSCA(I)   = ONE
  10  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          IF ((I.GT.N).OR.(I.LE.0)) GOTO 100
          J = ICN(K)
          IF (I.EQ.J) THEN
            VDIAG = abs(VAL(K))
            IF (VDIAG.GT.ZERO) THEN
              ROWSCA(J) = ONE/(sqrt(VDIAG))
            ENDIF
          ENDIF
 100   CONTINUE
       DO 110 I=1,N
        COLSCA(I) = ROWSCA(I)
 110   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF DIAGONAL SCALING'
      RETURN
      END SUBROUTINE DMUMPS_238
      SUBROUTINE DMUMPS_240(NSCA,N,NZ,IRN,ICN,VAL,
     &    RNOR,ROWSCA,MPRINT)
      INTEGER N, NZ, NSCA
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION VAL(NZ)
      DOUBLE PRECISION RNOR(N)
      DOUBLE PRECISION ROWSCA(N)
      DOUBLE PRECISION VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      DO 50 J=1,N
       RNOR(J)   = ZERO
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     &        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = abs(VAL(K))
            IF (VDIAG.GT.RNOR(I)) THEN
              RNOR(I) =  VDIAG
            ENDIF
 100   CONTINUE
      DO 130 J=1,N
       IF (RNOR(J).LE.ZERO) THEN
         RNOR(J)   = ONE
       ELSE
         RNOR(J)   = ONE/RNOR(J)
       ENDIF
 130  CONTINUE
      DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I)* RNOR(I)
 110  CONTINUE
      IF ( (NSCA.EQ.4) .OR. (NSCA.EQ.6) ) THEN
        DO 150 K=1,NZ
          I   = IRN(K)
          J   = ICN(K)
          IF (min(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 150
          VAL(K) = VAL(K) * RNOR(I)
 150    CONTINUE
      ENDIF
      IF (MPRINT.GT.0)
     &  WRITE(MPRINT,'(A)') '  END OF ROW SCALING'
      RETURN
      END SUBROUTINE DMUMPS_240
      SUBROUTINE DMUMPS_216(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION    R(M),C(N)
      DOUBLE PRECISION W(M*2+N*3)
      INTEGER LP,IFAIL
      INTRINSIC log,abs,min
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      DOUBLE PRECISION ONE
      DOUBLE PRECISION SMIN,ZERO
      PARAMETER (ONE=1.0D0,SMIN=0.1D0,ZERO=0.0D0)
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      DOUBLE PRECISION E,E1,EM,Q,Q1,QM,S,S1,SM,U,V
      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF
      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3
      DO 10 I = 1,M
         R(I) = ZERO
         W(I1+I) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = ZERO
         W(I2+J) = ZERO
         W(I3+J) = ZERO
         W(I4+J) = ZERO
   20 CONTINUE
      DO 30 K = 1,NE
         U = abs(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (min(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = log(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*dble(NE)
      DO 60 K = 1,NE
         IF (abs(A(K)).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (min(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
      E = ZERO
      Q = ONE
      S = ZERO
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (abs(S).LE.abs(SM)) GO TO 160
      DO 150 ITER = 1,MAXIT
         DO 80 K = 1,NE
            IF (abs(A(K)).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (min(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = ZERO
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
         IF (abs(S).LE.abs(SM)) E = ZERO
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (abs(S).LE.abs(SM)) GO TO 180
         EM = E*E1
         DO 110 K = 1,NE
            IF (abs(A(K)).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (min(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = ZERO
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
         IF (abs(S).LE.abs(SM)) Q = ONE
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
         IF (abs(S).LE.abs(SM)) GO TO 160
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
  180 DO 190 K = 1,NE
         IF (abs(A(K)).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (min(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN
  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     &    ' **** Error return from DMUMPS_216 ****',' IFAIL =',IFAIL
      END SUBROUTINE DMUMPS_216
      SUBROUTINE DMUMPS_27( id,  ANORMINF, LSCAL )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      TYPE(DMUMPS_STRUC), TARGET :: id
      DOUBLE PRECISION, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      INTEGER, DIMENSION (:), POINTER :: KEEP,INFO
      INTEGER(8), DIMENSION (:), POINTER :: KEEP8
      LOGICAL :: I_AM_SLAVE
      DOUBLE PRECISION DUMMY(1)
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0)
      DOUBLE PRECISION, ALLOCATABLE :: SUMR(:), SUMR_LOC(:)
      INTEGER :: allocok, MTYPE, I
      INFO =>id%INFO
      KEEP =>id%KEEP
      KEEP8 =>id%KEEP8
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
      IF (id%MYID .EQ. MASTER) THEN
       ALLOCATE( SUMR( id%N ), stat =allocok )
       IF (allocok .GT.0 ) THEN
        id%INFO(1)=-13
        id%INFO(2)=id%N
        RETURN
       ENDIF
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
             IF (.NOT.LSCAL) THEN
              CALL DMUMPS_207(id%A(1),
     &          id%NZ, id%N,
     &          id%IRN(1), id%JCN(1),
     &          SUMR, KEEP(1),KEEP8(1) )
             ELSE
              CALL DMUMPS_289(id%A(1),
     &          id%NZ, id%N,
     &          id%IRN(1), id%JCN(1), 
     &          SUMR, KEEP(1), KEEP8(1),
     &          id%COLSCA(1))
             ENDIF
            ELSE
             MTYPE = 1
             IF (.NOT.LSCAL) THEN
              CALL DMUMPS_119(MTYPE, id%N,
     &           id%NELT, id%ELTPTR(1),
     &           id%LELTVAR, id%ELTVAR(1),
     &           id%NA_ELT, id%A_ELT(1),
     &           SUMR, KEEP(1),KEEP8(1) )
             ELSE
              CALL DMUMPS_135(MTYPE, id%N,
     &           id%NELT, id%ELTPTR(1),
     &           id%LELTVAR, id%ELTVAR(1),
     &           id%NA_ELT, id%A_ELT(1),
     &           SUMR, KEEP(1),KEEP8(1), id%COLSCA(1))
             ENDIF
            ENDIF
          ENDIF
      ELSE
          ALLOCATE( SUMR_LOC( id%N ), stat =allocok )
          IF (allocok .GT.0 ) THEN
             id%INFO(1)=-13
             id%INFO(2)=id%N
             RETURN
          ENDIF
          IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
           IF (.NOT.LSCAL) THEN
              CALL DMUMPS_207(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%IRN_loc(1), id%JCN_loc(1), 
     &          SUMR_LOC, id%KEEP(1),id%KEEP8(1) )
           ELSE
              CALL DMUMPS_289(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%IRN_loc(1), id%JCN_loc(1), 
     &          SUMR_LOC, id%KEEP(1),id%KEEP8(1),
     &          id%COLSCA(1))
           ENDIF
          ELSE
           SUMR_LOC = ZERO
          ENDIF
          IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SUMR_LOC, SUMR,
     &        id%N, MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MASTER,id%COMM, IERR)
          ELSE
              CALL MPI_REDUCE( SUMR_LOC, DUMMY,
     &        id%N, MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MASTER,id%COMM, IERR)
          END IF
        DEALLOCATE (SUMR_LOC)
      ENDIF
      IF ( id%MYID .eq. MASTER ) THEN
       ANORMINF = dble(ZERO)
        IF (LSCAL) THEN
         DO I = 1, id%N
          ANORMINF = max(abs(id%ROWSCA(I) * SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ELSE
         DO I = 1, id%N
          ANORMINF = max(abs(SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ENDIF
      ENDIF
      CALL MPI_BCAST(ANORMINF, 1,
     &              MPI_DOUBLE_PRECISION, MASTER,
     &              id%COMM, IERR )
      IF (id%MYID .eq. MASTER) DEALLOCATE (SUMR)
      RETURN
      END SUBROUTINE DMUMPS_27
      SUBROUTINE DMUMPS_693(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     M, N, NUMPROCS, MYID, COMM,
     &     RPARTVEC, CPARTVEC,
     &     RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     ROWSCA, COLSCA, WRKRC, ISZWRKRC,
     &     SYM, NB1, NB2, NB3, EPS,
     &     ONENORMERR,INFNORMERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER NZ_loc, M, N, IWRKSZ, OP
      INTEGER NUMPROCS, MYID, COMM
      INTEGER INTSZ, RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      DOUBLE PRECISION A_loc(NZ_loc)
      INTEGER RPARTVEC(M), RSNDRCVSZ(2*NUMPROCS)
      INTEGER CPARTVEC(N), CSNDRCVSZ(2*NUMPROCS)
      INTEGER IWRK(IWRKSZ)
      INTEGER REGISTRE(12)
      DOUBLE PRECISION ROWSCA(M)
      DOUBLE PRECISION COLSCA(N)
      INTEGER ISZWRKRC
      DOUBLE PRECISION WRKRC(ISZWRKRC)
      DOUBLE PRECISION ONENORMERR,INFNORMERR
      INTEGER SYM, NB1, NB2, NB3
      DOUBLE PRECISION EPS
      EXTERNAL DMUMPS_694,DMUMPS_687, 
     &     DMUMPS_670
      INTEGER I
      IF(SYM.EQ.0) THEN
         CALL DMUMPS_694(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &        M, N, NUMPROCS, MYID, COMM,
     &        RPARTVEC, CPARTVEC,
     &        RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &        IWRK, IWRKSZ,
     &        INTSZ, RESZ, OP,
     &        ROWSCA, COLSCA, WRKRC, ISZWRKRC,
     &        NB1, NB2, NB3, EPS,
     &        ONENORMERR, INFNORMERR)  
      ELSE
         CALL DMUMPS_687(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &        N, NUMPROCS, MYID, COMM,
     &        RPARTVEC, 
     &        RSNDRCVSZ, REGISTRE,
     &        IWRK, IWRKSZ,
     &        INTSZ, RESZ, OP,
     &        ROWSCA, WRKRC, ISZWRKRC,
     &        NB1, NB2, NB3, EPS,
     &        ONENORMERR, INFNORMERR)  
         DO I=1,N
            COLSCA(I) = ROWSCA(I)
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_693
      SUBROUTINE DMUMPS_694(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     M, N, NUMPROCS, MYID, COMM,
     &     RPARTVEC, CPARTVEC,
     &     RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     ROWSCA, COLSCA, WRKRC, ISZWRKRC,
     &     NB1, NB2, NB3, EPS,
     &     ONENORMERR, INFNORMERR)    
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER NZ_loc, M, N, IWRKSZ, OP
      INTEGER NUMPROCS, MYID, COMM
      INTEGER INTSZ, RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      DOUBLE PRECISION A_loc(NZ_loc)
      INTEGER RPARTVEC(M), RSNDRCVSZ(2*NUMPROCS)
      INTEGER CPARTVEC(N), CSNDRCVSZ(2*NUMPROCS)
      INTEGER IWRK(IWRKSZ)
      INTEGER REGISTRE(12)
      DOUBLE PRECISION ROWSCA(M)
      DOUBLE PRECISION COLSCA(N)
      INTEGER ISZWRKRC
      DOUBLE PRECISION WRKRC(ISZWRKRC)  
      DOUBLE PRECISION ONENORMERR,INFNORMERR
      INTEGER IRSNDRCVNUM, ORSNDRCVNUM
      INTEGER IRSNDRCVVOL, ORSNDRCVVOL
      INTEGER ICSNDRCVNUM, OCSNDRCVNUM
      INTEGER ICSNDRCVVOL, OCSNDRCVVOL
      INTEGER  INUMMYR, INUMMYC
      INTEGER IMYRPTR,IMYCPTR 
      INTEGER IRNGHBPRCS, IRSNDRCVIA,IRSNDRCVJA
      INTEGER ORNGHBPRCS, ORSNDRCVIA,ORSNDRCVJA
      INTEGER ICNGHBPRCS, ICSNDRCVIA,ICSNDRCVJA
      INTEGER OCNGHBPRCS, OCSNDRCVIA,OCSNDRCVJA
      INTEGER ISTATUS, REQUESTS, TMPWORK
      INTEGER ITDRPTR, ITDCPTR, ISRRPTR
      INTEGER OSRRPTR, ISRCPTR, OSRCPTR
      INTEGER NB1, NB2, NB3
      DOUBLE PRECISION EPS
      INTEGER ITER, NZIND, IR, IC
      DOUBLE PRECISION ELM
      INTEGER TAG_COMM_COL
      PARAMETER(TAG_COMM_COL=100)
      INTEGER TAG_COMM_ROW
      PARAMETER(TAG_COMM_ROW=101)
      INTEGER TAG_ITERS
      PARAMETER(TAG_ITERS=102)
      EXTERNAL DMUMPS_654,
     &     DMUMPS_672, 
     &     DMUMPS_674,
     &     DMUMPS_662, 
     &     DMUMPS_743,
     &     DMUMPS_745,
     &     DMUMPS_660,
     &     DMUMPS_670,
     &     DMUMPS_671,
     &     DMUMPS_657,
     &     DMUMPS_656
      INTEGER DMUMPS_743 
      INTEGER DMUMPS_745
      DOUBLE PRECISION DMUMPS_737
      DOUBLE PRECISION DMUMPS_738
      INTRINSIC abs
      DOUBLE PRECISION RONE, RZERO
      PARAMETER(RONE=1.0D0,RZERO=0.0D0)
      INTEGER RESZR, RESZC
      INTEGER INTSZR, INTSZC
      INTEGER MAXMN
      INTEGER I, IERROR
      DOUBLE PRECISION ONEERRROW, ONEERRCOL, ONEERRL, ONEERRG
      DOUBLE PRECISION INFERRROW, INFERRCOL, INFERRL, INFERRG
      INTEGER OORANGEIND
      INFERRG = -RONE
      ONEERRG = -RONE
      OORANGEIND = 0
      MAXMN = M
      IF(MAXMN < N) MAXMN = N
      IF(OP == 1) THEN
         IF(NUMPROCS > 1) THEN
            CALL DMUMPS_654(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           RPARTVEC, M, N,
     &           IWRK, IWRKSZ)
            CALL DMUMPS_654(MYID, NUMPROCS, COMM,
     &           JCN_loc, IRN_loc,  NZ_loc,
     &           CPARTVEC, N, M,
     &           IWRK, IWRKSZ)         
            CALL DMUMPS_672(MYID, NUMPROCS, M, RPARTVEC,
     &           NZ_loc, IRN_loc, N, JCN_loc,
     &           IRSNDRCVNUM,IRSNDRCVVOL,
     &           ORSNDRCVNUM,ORSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS), COMM)
            CALL DMUMPS_672(MYID, NUMPROCS, N, CPARTVEC,
     &           NZ_loc, JCN_loc, M, IRN_loc, 
     &           ICSNDRCVNUM,ICSNDRCVVOL,
     &           OCSNDRCVNUM,OCSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           CSNDRCVSZ(1), CSNDRCVSZ(1+NUMPROCS), COMM)         
            CALL DMUMPS_662(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           RPARTVEC, CPARTVEC, M, N,
     &           INUMMYR,
     &           INUMMYC,     
     &           IWRK, IWRKSZ)
            INTSZR =  IRSNDRCVNUM + ORSNDRCVNUM + 
     &           IRSNDRCVVOL + ORSNDRCVVOL +
     &           2*(NUMPROCS+1) + INUMMYR
            INTSZC = ICSNDRCVNUM + OCSNDRCVNUM + 
     &           ICSNDRCVVOL + OCSNDRCVVOL +
     &           2*(NUMPROCS+1) + INUMMYC
            INTSZ = INTSZR + INTSZC + MAXMN + 
     &           (MPI_STATUS_SIZE +1) *  NUMPROCS
         ELSE
             IRSNDRCVNUM = 0
             ORSNDRCVNUM = 0
             IRSNDRCVVOL = 0 
             ORSNDRCVVOL = 0
             INUMMYR = 0
             ICSNDRCVNUM  = 0 
             OCSNDRCVNUM  = 0
             ICSNDRCVVOL = 0  
             OCSNDRCVVOL  = 0
             INUMMYC = 0
             INTSZ = 0
          ENDIF
          RESZR = M + IRSNDRCVVOL + ORSNDRCVVOL
          RESZC = N + ICSNDRCVVOL + OCSNDRCVVOL
          RESZ = RESZR  + RESZC  
          REGISTRE(1) = IRSNDRCVNUM 
          REGISTRE(2) = ORSNDRCVNUM
          REGISTRE(3) = IRSNDRCVVOL 
          REGISTRE(4) = ORSNDRCVVOL
          REGISTRE(5) = ICSNDRCVNUM 
          REGISTRE(6) = OCSNDRCVNUM
          REGISTRE(7) = ICSNDRCVVOL
          REGISTRE(8) = OCSNDRCVVOL
          REGISTRE(9) = INUMMYR
          REGISTRE(10) = INUMMYC
          REGISTRE(11) = INTSZ
          REGISTRE(12) = RESZ
       ELSE
          IRSNDRCVNUM = REGISTRE(1) 
          ORSNDRCVNUM = REGISTRE(2) 
          IRSNDRCVVOL = REGISTRE(3)
          ORSNDRCVVOL = REGISTRE(4) 
          ICSNDRCVNUM = REGISTRE(5) 
          OCSNDRCVNUM = REGISTRE(6) 
          ICSNDRCVVOL = REGISTRE(7) 
          OCSNDRCVVOL = REGISTRE(8) 
          INUMMYR = REGISTRE(9) 
          INUMMYC = REGISTRE(10)
          IF(NUMPROCS > 1) THEN
             CALL DMUMPS_660(MYID, NUMPROCS,COMM,    
     &            IRN_loc, JCN_loc, NZ_loc,
     &            RPARTVEC, CPARTVEC, M, N,
     &            IWRK(1), INUMMYR,
     &            IWRK(1+INUMMYR), INUMMYC,     
     &            IWRK(1+INUMMYR+INUMMYC), IWRKSZ-INUMMYR-INUMMYC )
             IMYRPTR = 1
             IMYCPTR = IMYRPTR + INUMMYR
             IRNGHBPRCS = IMYCPTR+ INUMMYC
             IRSNDRCVIA = IRNGHBPRCS+IRSNDRCVNUM
             IRSNDRCVJA = IRSNDRCVIA + NUMPROCS+1
             ORNGHBPRCS = IRSNDRCVJA + IRSNDRCVVOL
             ORSNDRCVIA = ORNGHBPRCS + ORSNDRCVNUM
             ORSNDRCVJA = ORSNDRCVIA + NUMPROCS + 1
             ICNGHBPRCS = ORSNDRCVJA + ORSNDRCVVOL
             ICSNDRCVIA = ICNGHBPRCS + ICSNDRCVNUM
             ICSNDRCVJA = ICSNDRCVIA + NUMPROCS+1
             OCNGHBPRCS = ICSNDRCVJA + ICSNDRCVVOL
             OCSNDRCVIA = OCNGHBPRCS + OCSNDRCVNUM
             OCSNDRCVJA = OCSNDRCVIA +  NUMPROCS + 1
             REQUESTS = OCSNDRCVJA + OCSNDRCVVOL
             ISTATUS = REQUESTS + NUMPROCS
             TMPWORK = ISTATUS + MPI_STATUS_SIZE *  NUMPROCS
             CALL DMUMPS_674(MYID, NUMPROCS, M, RPARTVEC,
     &            NZ_loc, IRN_loc,N, JCN_loc,
     &            IRSNDRCVNUM, IRSNDRCVVOL, 
     &            IWRK(IRNGHBPRCS),IWRK(IRSNDRCVIA),IWRK(IRSNDRCVJA),
     &            ORSNDRCVNUM, ORSNDRCVVOL, 
     &            IWRK(ORNGHBPRCS),IWRK(ORSNDRCVIA),IWRK(ORSNDRCVJA),
     &            RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS), IWRK(REQUESTS),
     &            TAG_COMM_ROW, COMM)
             CALL DMUMPS_674(MYID, NUMPROCS, N, CPARTVEC,
     &            NZ_loc, JCN_loc, M, IRN_loc,
     &            ICSNDRCVNUM, ICSNDRCVVOL, 
     &            IWRK(ICNGHBPRCS),
     &            IWRK(ICSNDRCVIA),
     &            IWRK(ICSNDRCVJA),
     &            OCSNDRCVNUM, OCSNDRCVVOL, 
     &            IWRK(OCNGHBPRCS),IWRK(OCSNDRCVIA),IWRK(OCSNDRCVJA),
     &            CSNDRCVSZ(1), CSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS),  IWRK(REQUESTS),
     &            TAG_COMM_COL, COMM)
             CALL DMUMPS_670(ROWSCA, M, RZERO)
             CALL DMUMPS_670(COLSCA, N, RZERO)
             CALL DMUMPS_671(ROWSCA, M, 
     &            IWRK(IMYRPTR),INUMMYR, RONE)
             CALL DMUMPS_671(COLSCA, N, 
     &            IWRK(IMYCPTR),INUMMYC, RONE)   
          ELSE
             CALL DMUMPS_670(ROWSCA, M, RONE)
             CALL DMUMPS_670(COLSCA, N, RONE)            
          ENDIF
          ITDRPTR = 1
          ITDCPTR = ITDRPTR + M
          ISRRPTR = ITDCPTR + N
          OSRRPTR = ISRRPTR + IRSNDRCVVOL
          ISRCPTR = OSRRPTR + ORSNDRCVVOL
          OSRCPTR = ISRCPTR + ICSNDRCVVOL
          IF(NUMPROCS == 1)THEN
             OSRCPTR = OSRCPTR - 1
             ISRCPTR = ISRCPTR - 1
             OSRRPTR = OSRRPTR - 1
             ISRRPTR = ISRRPTR - 1
          ELSE
             IF(IRSNDRCVVOL == 0) ISRRPTR = ISRRPTR - 1
             IF(ORSNDRCVVOL == 0) OSRRPTR = OSRRPTR - 1
             IF(ICSNDRCVVOL == 0) ISRCPTR = ISRCPTR - 1
             IF(OCSNDRCVVOL == 0) OSRCPTR = OSRCPTR - 1
          ENDIF
          ITER = 1
          DO WHILE (ITER.LE.NB1+NB2+NB3)
             IF(NUMPROCS > 1) THEN
                CALL DMUMPS_650(WRKRC(ITDRPTR),M,
     &               IWRK(IMYRPTR),INUMMYR)
                CALL DMUMPS_650(WRKRC(ITDCPTR),N,
     &               IWRK(IMYCPTR),INUMMYC)
             ELSE
                CALL DMUMPS_670(WRKRC(ITDRPTR),M, RZERO)
                CALL DMUMPS_670(WRKRC(ITDCPTR),N, RZERO)
             ENDIF
             IF((ITER.LE.NB1).OR.(ITER > NB1+NB2)) THEN
                IF((ITER.EQ.1).OR.(OORANGEIND.EQ.1)) THEN
                   DO NZIND=1,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                         IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                            WRKRC(ITDRPTR-1+IR)= ELM
                         ENDIF
                         IF(WRKRC(ITDCPTR-1+IC)<ELM) THEN
                            WRKRC(ITDCPTR-1+IC)= ELM
                         ENDIF
                      ELSE
                         OORANGEIND = 1
                      ENDIF
                   ENDDO
                ELSEIF(OORANGEIND.EQ.0) THEN
                   DO NZIND=1,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                      IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                         WRKRC(ITDRPTR-1+IR)= ELM
                      ENDIF
                      IF(WRKRC(ITDCPTR-1+IC)<ELM) THEN
                         WRKRC(ITDCPTR-1+IC)= ELM
                      ENDIF
                   ENDDO
                ENDIF
                IF(NUMPROCS > 1) THEN
                   CALL DMUMPS_657(MYID, NUMPROCS,
     &                  WRKRC(ITDCPTR), N, TAG_ITERS+ITER, 
     &                  ICSNDRCVNUM,IWRK(ICNGHBPRCS),
     &                  ICSNDRCVVOL,IWRK(ICSNDRCVIA), IWRK(ICSNDRCVJA), 
     &                  WRKRC(ISRCPTR),
     &                  OCSNDRCVNUM,IWRK(OCNGHBPRCS),
     &                  OCSNDRCVVOL,IWRK(OCSNDRCVIA), IWRK(OCSNDRCVJA),
     &                  WRKRC( OSRCPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
                  CALL DMUMPS_657(MYID, NUMPROCS,
     &                  WRKRC(ITDRPTR), M, TAG_ITERS+2+ITER, 
     &                  IRSNDRCVNUM,IWRK(IRNGHBPRCS),
     &                  IRSNDRCVVOL,IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                  WRKRC(ISRRPTR),
     &                  ORSNDRCVNUM,IWRK(ORNGHBPRCS),
     &                  ORSNDRCVVOL,IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                  WRKRC( OSRRPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRROW = DMUMPS_737(ROWSCA, 
     &                    WRKRC(ITDRPTR), M,
     &                    IWRK(IMYRPTR),INUMMYR)
                     INFERRCOL = DMUMPS_737(COLSCA,  
     &                    WRKRC(ITDCPTR), N,
     &                    IWRK(IMYCPTR),INUMMYC)
                     INFERRL = INFERRCOL
                     IF(INFERRROW > INFERRL ) THEN
                        INFERRL = INFERRROW                   
                     ENDIF
                     CALL MPI_ALLREDUCE(INFERRL, INFERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_665(COLSCA,  WRKRC(ITDCPTR), 
     &                       N,
     &                       IWRK(IMYCPTR),INUMMYC)
                        CALL DMUMPS_665(ROWSCA,  WRKRC(ITDRPTR), 
     &                       M,
     &                       IWRK(IMYRPTR),INUMMYR)         
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF                  
               ELSE
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRROW = DMUMPS_738(ROWSCA, 
     &                    WRKRC(ITDRPTR), M)
                     INFERRCOL = DMUMPS_738(COLSCA,  
     &                    WRKRC(ITDCPTR), N)
                     INFERRL = INFERRCOL
                     IF(INFERRROW > INFERRL) THEN
                        INFERRL = INFERRROW                    
                     ENDIF                     
                     INFERRG = INFERRL
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_666(COLSCA,  WRKRC(ITDCPTR), N)
                        CALL DMUMPS_666(ROWSCA,  WRKRC(ITDRPTR), M)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF 
                  ENDIF
               ENDIF
            ELSE
               IF((ITER .EQ.1).OR.(OORANGEIND.EQ.1))THEN
                  DO NZIND=1,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                        WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                        WRKRC(ITDCPTR-1+IC) = WRKRC(ITDCPTR-1+IC) + ELM
                     ELSE
                        OORANGEIND = 1
                     ENDIF
                  ENDDO
               ELSEIF(OORANGEIND.EQ.0) THEN
                  DO NZIND=1,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                     WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                     WRKRC(ITDCPTR-1+IC) = WRKRC(ITDCPTR-1+IC) + ELM
                  ENDDO
               ENDIF
               IF(NUMPROCS > 1) THEN                 
                  CALL DMUMPS_656(MYID, NUMPROCS,
     &                 WRKRC(ITDCPTR), N, TAG_ITERS+ITER, 
     &                 ICSNDRCVNUM, IWRK(ICNGHBPRCS),
     &                 ICSNDRCVVOL, IWRK(ICSNDRCVIA), IWRK(ICSNDRCVJA), 
     &                 WRKRC(ISRCPTR),
     &                 OCSNDRCVNUM, IWRK(OCNGHBPRCS),
     &                 OCSNDRCVVOL, IWRK(OCSNDRCVIA), IWRK(OCSNDRCVJA),
     &                 WRKRC( OSRCPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
                  CALL DMUMPS_656(MYID, NUMPROCS,
     &                 WRKRC(ITDRPTR), M, TAG_ITERS+2+ITER, 
     &                 IRSNDRCVNUM, IWRK(IRNGHBPRCS),
     &                 IRSNDRCVVOL, IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                 WRKRC(ISRRPTR),
     &                 ORSNDRCVNUM, IWRK(ORNGHBPRCS),
     &                 ORSNDRCVVOL, IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                 WRKRC( OSRRPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRROW = DMUMPS_737(ROWSCA, 
     &                    WRKRC(ITDRPTR), M,
     &                    IWRK(IMYRPTR),INUMMYR)
                     ONEERRCOL = DMUMPS_737(COLSCA,  
     &                    WRKRC(ITDCPTR), N,
     &                    IWRK(IMYCPTR),INUMMYC)
                     ONEERRL = ONEERRCOL
                     IF(ONEERRROW > ONEERRL ) THEN
                        ONEERRL = ONEERRROW                   
                     ENDIF
                     CALL MPI_ALLREDUCE(ONEERRL, ONEERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_665(COLSCA,  WRKRC(ITDCPTR), 
     &                       N,
     &                       IWRK(IMYCPTR),INUMMYC)
                        CALL DMUMPS_665(ROWSCA,  WRKRC(ITDRPTR), 
     &                       M,
     &                       IWRK(IMYRPTR),INUMMYR)          
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF                            
               ELSE
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRROW = DMUMPS_738(ROWSCA, 
     &                    WRKRC(ITDRPTR), M)
                     ONEERRCOL = DMUMPS_738(COLSCA,  
     &                    WRKRC(ITDCPTR), N)
                     ONEERRL = ONEERRCOL
                     IF(ONEERRROW > ONEERRL) THEN
                        ONEERRL = ONEERRROW                    
                     ENDIF                     
                     ONEERRG = ONEERRL
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_666(COLSCA,  WRKRC(ITDCPTR), N)
                        CALL DMUMPS_666(ROWSCA,  WRKRC(ITDRPTR), M)
                        ITER = NB1+NB2+1                        
                        CYCLE
                     ENDIF
                  ENDIF                  
               ENDIF 
            ENDIF
            IF(NUMPROCS > 1) THEN               
               CALL DMUMPS_665(COLSCA,  WRKRC(ITDCPTR), N,
     &              IWRK(IMYCPTR),INUMMYC)
               CALL DMUMPS_665(ROWSCA,  WRKRC(ITDRPTR), M,
     &              IWRK(IMYRPTR),INUMMYR)               
            ELSE
               CALL DMUMPS_666(COLSCA,  WRKRC(ITDCPTR), N)
               CALL DMUMPS_666(ROWSCA,  WRKRC(ITDRPTR), M)
            ENDIF
            ITER = ITER + 1
         ENDDO
         ONENORMERR = ONEERRG 
         INFNORMERR = INFERRG 
         IF(NUMPROCS > 1) THEN
            CALL MPI_REDUCE(ROWSCA, WRKRC(1), M, MPI_DOUBLE_PRECISION,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            IF(MYID.EQ.0) THEN
               DO I=1, M
                  ROWSCA(I) = WRKRC(I)
               ENDDO
            ENDIF
            CALL MPI_REDUCE(COLSCA, WRKRC(1+M), N, MPI_DOUBLE_PRECISION,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            If(MYID.EQ.0) THEN
               DO I=1, N
                  COLSCA(I) = WRKRC(I+M)
               ENDDO
            ENDIF         
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_694
      SUBROUTINE DMUMPS_687(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     N, NUMPROCS, MYID, COMM,
     &     PARTVEC, 
     &     RSNDRCVSZ, 
     &     REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     SCA, WRKRC, ISZWRKRC,
     &     NB1, NB2, NB3, EPS,
     &     ONENORMERR, INFNORMERR)    
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER NZ_loc, N, IWRKSZ, OP
      INTEGER NUMPROCS, MYID, COMM
      INTEGER INTSZ, RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      DOUBLE PRECISION A_loc(NZ_loc)
      INTEGER PARTVEC(N), RSNDRCVSZ(2*NUMPROCS)
      INTEGER IWRK(IWRKSZ)
      INTEGER REGISTRE(12)
      DOUBLE PRECISION SCA(N)
      INTEGER ISZWRKRC
      DOUBLE PRECISION WRKRC(ISZWRKRC)
      INTEGER IRSNDRCVNUM, ORSNDRCVNUM
      INTEGER IRSNDRCVVOL, ORSNDRCVVOL
      INTEGER  INUMMYR
      INTEGER IMYRPTR,IMYCPTR 
      INTEGER IRNGHBPRCS, IRSNDRCVIA,IRSNDRCVJA
      INTEGER ORNGHBPRCS, ORSNDRCVIA,ORSNDRCVJA
      INTEGER ISTATUS, REQUESTS, TMPWORK
      INTEGER ITDRPTR, ISRRPTR, OSRRPTR
      DOUBLE PRECISION ONENORMERR,INFNORMERR
      INTEGER NB1, NB2, NB3
      DOUBLE PRECISION EPS
      INTEGER ITER, NZIND, IR, IC
      DOUBLE PRECISION ELM
      INTEGER TAG_COMM_ROW
      PARAMETER(TAG_COMM_ROW=101)
      INTEGER TAG_ITERS
      PARAMETER(TAG_ITERS=102)
      EXTERNAL DMUMPS_655,
     &     DMUMPS_673, 
     &     DMUMPS_692,
     &     DMUMPS_663, 
     &     DMUMPS_742,
     &     DMUMPS_745,
     &     DMUMPS_661,
     &     DMUMPS_657,
     &     DMUMPS_656,
     &     DMUMPS_670,
     &     DMUMPS_671
      INTEGER DMUMPS_742 
      INTEGER DMUMPS_745
      DOUBLE PRECISION DMUMPS_737
      DOUBLE PRECISION DMUMPS_738
      INTRINSIC abs
      DOUBLE PRECISION RONE, RZERO
      PARAMETER(RONE=1.0D0,RZERO=0.0D0)
      INTEGER INTSZR
      INTEGER MAXMN
      INTEGER I, IERROR
      DOUBLE PRECISION ONEERRL, ONEERRG
      DOUBLE PRECISION INFERRL, INFERRG
      INTEGER OORANGEIND
      OORANGEIND = 0
      INFERRG = -RONE
      ONEERRG = -RONE
      MAXMN = N
      IF(OP == 1) THEN
         IF(NUMPROCS > 1) THEN
            CALL DMUMPS_655(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           PARTVEC, N,
     &           IWRK, IWRKSZ)
            CALL DMUMPS_673(MYID, NUMPROCS, N, PARTVEC,
     &           NZ_loc, IRN_loc,JCN_loc, IRSNDRCVNUM,IRSNDRCVVOL,
     &           ORSNDRCVNUM, ORSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS), COMM)
            CALL DMUMPS_663(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           PARTVEC, N,
     &           INUMMYR,
     &           IWRK, IWRKSZ)
            INTSZR =  IRSNDRCVNUM + ORSNDRCVNUM + 
     &           IRSNDRCVVOL + ORSNDRCVVOL +
     &           2*(NUMPROCS+1) + INUMMYR
            INTSZ = INTSZR + N + 
     &           (MPI_STATUS_SIZE +1) *  NUMPROCS
         ELSE
            IRSNDRCVNUM = 0
            ORSNDRCVNUM = 0
            IRSNDRCVVOL = 0 
            ORSNDRCVVOL = 0
            INUMMYR = 0
            INTSZ = 0
         ENDIF
         RESZ = N + IRSNDRCVVOL + ORSNDRCVVOL
         REGISTRE(1) = IRSNDRCVNUM 
         REGISTRE(2) = ORSNDRCVNUM
         REGISTRE(3) = IRSNDRCVVOL 
         REGISTRE(4) = ORSNDRCVVOL
         REGISTRE(9) = INUMMYR
         REGISTRE(11) = INTSZ
         REGISTRE(12) = RESZ
      ELSE
         IRSNDRCVNUM = REGISTRE(1) 
         ORSNDRCVNUM = REGISTRE(2) 
         IRSNDRCVVOL = REGISTRE(3)
         ORSNDRCVVOL = REGISTRE(4) 
         INUMMYR = REGISTRE(9) 
          IF(NUMPROCS > 1) THEN
             CALL DMUMPS_661(MYID, NUMPROCS,COMM,    
     &            IRN_loc, JCN_loc, NZ_loc,
     &            PARTVEC, N,
     &            IWRK(1), INUMMYR,
     &            IWRK(1+INUMMYR), IWRKSZ-INUMMYR)
             IMYRPTR = 1
             IMYCPTR = IMYRPTR + INUMMYR
             IRNGHBPRCS = IMYCPTR 
             IRSNDRCVIA = IRNGHBPRCS+IRSNDRCVNUM
             IRSNDRCVJA = IRSNDRCVIA + NUMPROCS+1
             ORNGHBPRCS = IRSNDRCVJA + IRSNDRCVVOL
             ORSNDRCVIA = ORNGHBPRCS + ORSNDRCVNUM
             ORSNDRCVJA = ORSNDRCVIA + NUMPROCS + 1
             REQUESTS = ORSNDRCVJA + ORSNDRCVVOL 
             ISTATUS = REQUESTS + NUMPROCS
             TMPWORK = ISTATUS + MPI_STATUS_SIZE *  NUMPROCS
             CALL DMUMPS_692(MYID, NUMPROCS, N, PARTVEC,
     &            NZ_loc, IRN_loc, JCN_loc,
     &            IRSNDRCVNUM, IRSNDRCVVOL, 
     &            IWRK(IRNGHBPRCS),IWRK(IRSNDRCVIA),IWRK(IRSNDRCVJA),
     &            ORSNDRCVNUM, ORSNDRCVVOL, 
     &            IWRK(ORNGHBPRCS),IWRK(ORSNDRCVIA),IWRK(ORSNDRCVJA),
     &            RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS), IWRK(REQUESTS),
     &            TAG_COMM_ROW, COMM)
             CALL DMUMPS_670(SCA, N, RZERO)
             CALL DMUMPS_671(SCA, N, 
     &            IWRK(IMYRPTR),INUMMYR, RONE)
          ELSE
             CALL DMUMPS_670(SCA, N, RONE)
          ENDIF
          ITDRPTR = 1
          ISRRPTR = ITDRPTR + N
          OSRRPTR = ISRRPTR + IRSNDRCVVOL
          IF(NUMPROCS == 1)THEN
             OSRRPTR = OSRRPTR - 1
             ISRRPTR = ISRRPTR - 1
          ELSE
             IF(IRSNDRCVVOL == 0) ISRRPTR = ISRRPTR - 1
             IF(ORSNDRCVVOL == 0) OSRRPTR = OSRRPTR - 1
          ENDIF
          ITER = 1
          DO WHILE(ITER.LE.NB1+NB2+NB3)
             IF(NUMPROCS > 1) THEN
                CALL DMUMPS_650(WRKRC(ITDRPTR),N,
     &               IWRK(IMYRPTR),INUMMYR)
             ELSE
                CALL DMUMPS_670(WRKRC(ITDRPTR),N, RZERO)
             ENDIF
             IF((ITER.LE.NB1).OR.(ITER > NB1+NB2)) THEN
                IF((ITER.EQ.1).OR.(OORANGEIND.EQ.1)) THEN
                   DO NZIND=1,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                         IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                            WRKRC(ITDRPTR-1+IR)= ELM
                         ENDIF
                         IF(WRKRC(ITDRPTR-1+IC)<ELM) THEN
                            WRKRC(ITDRPTR-1+IC)= ELM
                         ENDIF
                      ELSE
                         OORANGEIND = 1
                      ENDIF
                   ENDDO
                ELSEIF(OORANGEIND.EQ.0) THEN
                   DO NZIND=1,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                      IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                         WRKRC(ITDRPTR-1+IR)= ELM
                      ENDIF
                      IF(WRKRC(ITDRPTR-1+IC)<ELM) THEN
                         WRKRC(ITDRPTR-1+IC)= ELM
                      ENDIF
                   ENDDO
                ENDIF                      
                IF(NUMPROCS > 1) THEN
                  CALL DMUMPS_657(MYID, NUMPROCS,
     &                  WRKRC(ITDRPTR), N, TAG_ITERS+2+ITER, 
     &                  IRSNDRCVNUM,IWRK(IRNGHBPRCS),
     &                  IRSNDRCVVOL,IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                  WRKRC(ISRRPTR),
     &                  ORSNDRCVNUM,IWRK(ORNGHBPRCS),
     &                  ORSNDRCVVOL,IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                  WRKRC( OSRRPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRL = DMUMPS_737(SCA,  
     &                    WRKRC(ITDRPTR), N,
     &                    IWRK(IMYRPTR),INUMMYR)                  
                     CALL MPI_ALLREDUCE(INFERRL, INFERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_665(SCA,  WRKRC(ITDRPTR), N,
     &                       IWRK(IMYRPTR),INUMMYR)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRL = DMUMPS_738(SCA, 
     &                    WRKRC(ITDRPTR), N)
                     INFERRG = INFERRL
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_666(SCA,  WRKRC(ITDRPTR), N)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF 
                  ENDIF
               ENDIF
            ELSE
               IF((ITER.EQ.1).OR.(OORANGEIND.EQ.1))THEN
                  DO NZIND=1,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                        WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                        IF(IR.NE.IC) THEN
                           WRKRC(ITDRPTR-1+IC) = 
     &                          WRKRC(ITDRPTR-1+IC) + ELM
                        ENDIF
                     ELSE
                        OORANGEIND = 1
                     ENDIF
                  ENDDO
               ELSEIF(OORANGEIND.EQ.0)THEN
                  DO NZIND=1,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                     WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                     IF(IR.NE.IC) THEN
                        WRKRC(ITDRPTR-1+IC) = WRKRC(ITDRPTR-1+IC) + ELM
                     ENDIF
                  ENDDO
               ENDIF
               IF(NUMPROCS > 1) THEN
                  CALL DMUMPS_656(MYID, NUMPROCS,
     &                 WRKRC(ITDRPTR), N, TAG_ITERS+2+ITER, 
     &                 IRSNDRCVNUM, IWRK(IRNGHBPRCS),
     &                 IRSNDRCVVOL, IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                 WRKRC(ISRRPTR),
     &                 ORSNDRCVNUM, IWRK(ORNGHBPRCS),
     &                 ORSNDRCVVOL, IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                 WRKRC( OSRRPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRL = DMUMPS_737(SCA,  
     &                    WRKRC(ITDRPTR), N,
     &                    IWRK(IMYRPTR),INUMMYR) 
                     CALL MPI_ALLREDUCE(ONEERRL, ONEERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_665(SCA,  WRKRC(ITDRPTR), N,
     &                       IWRK(IMYRPTR),INUMMYR)
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF
               ELSE
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRL = DMUMPS_738(SCA, 
     &                    WRKRC(ITDRPTR), N)
                     ONEERRG = ONEERRL
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_666(SCA,  WRKRC(ITDRPTR), N)
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF(NUMPROCS > 1) THEN
               CALL DMUMPS_665(SCA,  WRKRC(ITDRPTR), N,
     &              IWRK(IMYRPTR),INUMMYR)
            ELSE
               CALL DMUMPS_666(SCA,  WRKRC(ITDRPTR), N)
            ENDIF     
            ITER = ITER + 1
         ENDDO
         ONENORMERR = ONEERRG 
         INFNORMERR = INFERRG 
         IF(NUMPROCS > 1) THEN
            CALL MPI_REDUCE(SCA, WRKRC(1), N, MPI_DOUBLE_PRECISION,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            IF(MYID.EQ.0) THEN
               DO I=1, N
                  SCA(I) = WRKRC(I)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_687
      SUBROUTINE DMUMPS_654(MYID, NUMPROCS, COMM,
     & IRN_loc, JCN_loc, NZ_loc,
     & IPARTVEC, ISZ, OSZ,
     & IWRK, IWSZ)
      IMPLICIT NONE
      EXTERNAL DMUMPS_703
      INTEGER MYID, NUMPROCS, COMM
      INTEGER NZ_loc, ISZ, IWSZ, OSZ
      INTEGER IRN_loc(NZ_loc), JCN_loc(NZ_loc)
      INTEGER IPARTVEC(ISZ)
      INTEGER IWRK(IWSZ)
      INCLUDE 'mpif.h'
      INTEGER I
      INTEGER OP, IERROR
      INTEGER IR, IC
      IF(NUMPROCS.NE.1) THEN
         CALL MPI_OP_CREATE(DMUMPS_703, .TRUE., OP, IERROR)
         CALL DMUMPS_668(IWRK, 4*ISZ, ISZ)
         DO I=1,ISZ
            IWRK(2*I-1) = 0
            IWRK(2*I) = MYID
         ENDDO
         DO I=1,NZ_loc
            IR = IRN_loc(I)
            IC = JCN_loc(I)
            IF((IR.GE.1).AND.(IR.LE.ISZ).AND.
     &           (IC.GE.1).AND.(IC.LE.OSZ)) THEN
               IWRK(2*IR-1) = IWRK(2*IR-1) + 1
            ENDIF
         ENDDO
         CALL MPI_ALLREDUCE(IWRK(1), IWRK(1+2*ISZ), ISZ,
     &        MPI_2INTEGER, OP, COMM, IERROR)      
         DO I=1,ISZ
            IPARTVEC(I) = IWRK(2*I+2*ISZ)
         ENDDO
         CALL MPI_OP_FREE(OP, IERROR)
      ELSE
         DO I=1,ISZ
            IPARTVEC(I) = 0
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_654
      SUBROUTINE DMUMPS_662(MYID, NUMPROCS, COMM,
     &     IRN_loc, JCN_loc, NZ_loc,
     &     ROWPARTVEC, COLPARTVEC, M, N,
     &     INUMMYR,
     &     INUMMYC,     
     &     IWRK, IWSZ)
      IMPLICIT NONE
      INTEGER MYID, NUMPROCS, NZ_loc, M, N
      INTEGER IRN_loc(NZ_loc), JCN_loc(NZ_loc)
      INTEGER ROWPARTVEC(M)
      INTEGER COLPARTVEC(N)
      INTEGER INUMMYR, INUMMYC
      INTEGER IWSZ
      INTEGER IWRK(IWSZ)
      INTEGER COMM
      INTEGER I, IR, IC
      INUMMYR = 0
      INUMMYC = 0
      DO I=1,M
         IWRK(I) = 0
         IF(ROWPARTVEC(I).EQ.MYID) THEN
            IWRK(I)=1
            INUMMYR = INUMMYR + 1
         ENDIF
      ENDDO
      DO I=1,NZ_loc
         IR = IRN_loc(I)
         IC = JCN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.M).AND.
     &        ((IC.GE.1).AND.(IC.LE.N)) ) THEN
            IF(IWRK(IR) .EQ. 0) THEN
               IWRK(IR)= 1
               INUMMYR = INUMMYR + 1
            ENDIF
         ENDIF
      ENDDO
      DO I=1,N
         IWRK(I) = 0
         IF(COLPARTVEC(I).EQ.MYID) THEN
            IWRK(I)= 1
            INUMMYC = INUMMYC + 1
         ENDIF
      ENDDO
      DO I=1,NZ_loc
         IC = JCN_loc(I)
         IR = IRN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.M).AND.
     &        ((IC.GE.1).AND.(IC.LE.N)) ) THEN
            IF(IWRK(IC) .EQ. 0) THEN
               IWRK(IC)= 1
               INUMMYC = INUMMYC + 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_662
      SUBROUTINE DMUMPS_660(MYID, NUMPROCS,COMM,    
     &     IRN_loc, JCN_loc, NZ_loc,
     &     ROWPARTVEC, COLPARTVEC, M, N,
     &     MYROWINDICES, INUMMYR,
     &     MYCOLINDICES, INUMMYC,     
     &     IWRK, IWSZ  )
      IMPLICIT NONE
      INTEGER MYID, NUMPROCS, NZ_loc, M, N
      INTEGER INUMMYR, INUMMYC, IWSZ
      INTEGER IRN_loc(NZ_loc), JCN_loc(NZ_loc)
      INTEGER ROWPARTVEC(M)
      INTEGER COLPARTVEC(N)
      INTEGER MYROWINDICES(INUMMYR)
      INTEGER MYCOLINDICES(INUMMYC)
      INTEGER IWRK(IWSZ)
      INTEGER COMM
      INTEGER I, IR, IC, ITMP, MAXMN
      MAXMN = M
      IF(N > MAXMN) MAXMN = N
      DO I=1,M
         IWRK(I) = 0
         IF(ROWPARTVEC(I).EQ.MYID) IWRK(I)=1
      ENDDO
      DO I=1,NZ_loc
         IR = IRN_loc(I)         
         IC = JCN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.M).AND.
     &      ((IC.GE.1).AND.(IC.LE.N))  ) THEN
            IF(IWRK(IR) .EQ. 0) IWRK(IR)= 1
         ENDIF
      ENDDO
      ITMP = 1
      DO I=1,M
         IF(IWRK(I).EQ.1) THEN
            MYROWINDICES(ITMP) = I
            ITMP  = ITMP + 1
         ENDIF
      ENDDO
      DO I=1,N
         IWRK(I) = 0
         IF(COLPARTVEC(I).EQ.MYID) IWRK(I)= 1
      ENDDO
      DO I=1,NZ_loc
         IR = IRN_loc(I)         
         IC = JCN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.M).AND.
     &      ((IC.GE.1).AND.(IC.LE.N))  ) THEN
            IF(IWRK(IC) .EQ. 0) IWRK(IC)= 1
         ENDIF
      ENDDO
      ITMP = 1
      DO I=1,N
         IF(IWRK(I).EQ.1) THEN
            MYCOLINDICES(ITMP) = I
            ITMP  = ITMP + 1
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_660
      INTEGER FUNCTION DMUMPS_744(D, DSZ, INDX, INDXSZ, EPS)
      IMPLICIT NONE
      INTEGER DSZ, INDXSZ
      DOUBLE PRECISION D(DSZ)
      INTEGER INDX(INDXSZ)
      DOUBLE PRECISION EPS
      INTEGER I, IID
      DOUBLE PRECISION RONE
      PARAMETER(RONE=1.0D0)
      DMUMPS_744 = 1
      DO I=1, INDXSZ
         IID = INDX(I)
         IF (.NOT.( (D(IID).LE.(RONE+EPS)).AND.
     &        ((RONE-EPS).LE.D(IID)) )) THEN
            DMUMPS_744 = 0         
         ENDIF
      ENDDO
      RETURN
      END FUNCTION DMUMPS_744
      INTEGER FUNCTION DMUMPS_745(D, DSZ, EPS)
      IMPLICIT NONE
      INTEGER DSZ
      DOUBLE PRECISION D(DSZ)
      DOUBLE PRECISION EPS
      INTEGER I
      DOUBLE PRECISION RONE
      PARAMETER(RONE=1.0D0)
      DMUMPS_745 = 1
      DO I=1, DSZ
         IF (.NOT.( (D(I).LE.(RONE+EPS)).AND.
     &        ((RONE-EPS).LE.D(I)) )) THEN
            DMUMPS_745 = 0         
         ENDIF
      ENDDO
      RETURN
      END FUNCTION DMUMPS_745
      INTEGER FUNCTION DMUMPS_743(DR, M, INDXR, INDXRSZ,
     &     DC, N, INDXC, INDXCSZ, EPS, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER M, N, INDXRSZ, INDXCSZ
      DOUBLE PRECISION DR(M), DC(N)
      INTEGER INDXR(INDXRSZ), INDXC(INDXCSZ)
      DOUBLE PRECISION EPS
      INTEGER COMM
      EXTERNAL DMUMPS_744
      INTEGER  DMUMPS_744
      INTEGER GLORES, MYRESR, MYRESC, MYRES
      INTEGER IERR
      MYRESR =  DMUMPS_744(DR, M, INDXR, INDXRSZ, EPS)
      MYRESC =  DMUMPS_744(DC, N, INDXC, INDXCSZ, EPS)
      MYRES = MYRESR + MYRESC
      CALL MPI_ALLREDUCE(MYRES, GLORES, 1, MPI_INTEGER,
     &     MPI_SUM, COMM, IERR)
      DMUMPS_743 = GLORES
      RETURN
      END FUNCTION DMUMPS_743
      DOUBLE PRECISION FUNCTION DMUMPS_737(D, TMPD, DSZ,
     &     INDX, INDXSZ)
      IMPLICIT NONE 
      INTEGER DSZ, INDXSZ
      DOUBLE PRECISION D(DSZ)
      DOUBLE PRECISION TMPD(DSZ)
      INTEGER INDX(INDXSZ)
      DOUBLE PRECISION RONE
      PARAMETER(RONE=1.0D0)
      INTEGER I, IIND
      DOUBLE PRECISION ERRMAX
      INTRINSIC abs
      ERRMAX = -RONE
      DO I=1,INDXSZ
         IIND = INDX(I)
         IF(abs(RONE-TMPD(IIND)).GT.ERRMAX) THEN
            ERRMAX = abs(RONE-TMPD(IIND))
         ENDIF
      ENDDO           
      DMUMPS_737 = ERRMAX
      RETURN
      END FUNCTION DMUMPS_737
      DOUBLE PRECISION FUNCTION DMUMPS_738(D, TMPD, DSZ)
      IMPLICIT NONE 
      INTEGER DSZ
      DOUBLE PRECISION D(DSZ)
      DOUBLE PRECISION TMPD(DSZ)
      DOUBLE PRECISION RONE
      PARAMETER(RONE=1.0D0)
      INTEGER I
      DOUBLE PRECISION ERRMAX1
      INTRINSIC abs
      ERRMAX1 = -RONE
      DO I=1,DSZ
         IF(abs(RONE-TMPD(I)).GT.ERRMAX1) THEN
            ERRMAX1 = abs(RONE-TMPD(I))
         ENDIF
      ENDDO
      DMUMPS_738 = ERRMAX1
      RETURN
      END FUNCTION DMUMPS_738
      SUBROUTINE DMUMPS_665(D,  TMPD, DSZ,
     &        INDX, INDXSZ)
      IMPLICIT NONE
      INTEGER DSZ, INDXSZ
      DOUBLE PRECISION D(DSZ)
      DOUBLE PRECISION TMPD(DSZ)
      INTEGER INDX(INDXSZ)
      INTRINSIC sqrt
      INTEGER I, IIND
      DOUBLE PRECISION RZERO
      PARAMETER(RZERO=0.0D0)
      DO I=1,INDXSZ
         IIND = INDX(I)
         IF (TMPD(IIND).NE.RZERO) D(IIND) = D(IIND)/sqrt(TMPD(IIND))
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_665
      SUBROUTINE DMUMPS_666(D,  TMPD, DSZ)
      IMPLICIT NONE
      INTEGER DSZ
      DOUBLE PRECISION D(DSZ)
      DOUBLE PRECISION TMPD(DSZ)
      INTRINSIC sqrt
      INTEGER I
      DOUBLE PRECISION RZERO
      PARAMETER(RZERO=0.0D0)
      DO I=1,DSZ
         IF (TMPD(I) .NE. RZERO) D(I) = D(I)/sqrt(TMPD(I))
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_666
      SUBROUTINE DMUMPS_671(D, DSZ, INDX, INDXSZ, VAL)
      IMPLICIT NONE
      INTEGER DSZ, INDXSZ
      DOUBLE PRECISION D(DSZ)
      INTEGER INDX(INDXSZ)
      DOUBLE PRECISION VAL
      INTEGER I, IIND
      DO I=1,INDXSZ
         IIND = INDX(I)
         D(IIND) = VAL
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_671
      SUBROUTINE DMUMPS_702(D, DSZ, INDX, INDXSZ)
      IMPLICIT NONE
      INTEGER DSZ, INDXSZ
      DOUBLE PRECISION D(DSZ)
      INTEGER INDX(INDXSZ)
      INTEGER I, IIND
      DO I=1,INDXSZ
         IIND  = INDX(I)
         D(IIND) = 1.0D0/D(IIND)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_702
      SUBROUTINE DMUMPS_670(D, DSZ, VAL)
      IMPLICIT NONE
      INTEGER DSZ
      DOUBLE PRECISION D(DSZ)
      DOUBLE PRECISION VAL
      INTEGER I
      DO I=1,DSZ
         D(I) = VAL
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_670
      SUBROUTINE DMUMPS_650(TMPD, TMPSZ, INDX, INDXSZ)
      IMPLICIT NONE
      INTEGER TMPSZ,INDXSZ 
      DOUBLE PRECISION TMPD(TMPSZ)
      INTEGER INDX(INDXSZ)
      INTEGER I
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO=0.0D0)
      DO I=1,INDXSZ
         TMPD(INDX(I)) = DZERO
      ENDDO      
      RETURN
      END SUBROUTINE DMUMPS_650
      SUBROUTINE DMUMPS_703(INV, INOUTV, LEN, DTYPE)
      IMPLICIT NONE
      INTEGER LEN
      INTEGER INV(2*LEN) 
      INTEGER INOUTV(2*LEN)
      INTEGER DTYPE
      INTEGER I
      INTEGER DIN, DINOUT, PIN, PINOUT
      DO I=1,2*LEN-1,2
         DIN = INV(I)     
         PIN = INV(I+1)   
         DINOUT = INOUTV(I)
         PINOUT = INOUTV(I+1)
         IF (DINOUT < DIN) THEN
            INOUTV(I) = DIN
            INOUTV(I+1) = PIN
         ELSE IF (DINOUT == DIN) THEN
            IF ((mod(DINOUT,2).EQ.0).AND.(PIN<PINOUT)) THEN
              INOUTV(I+1) = PIN
            ELSE IF ((mod(DINOUT,2).EQ.1).AND.(PIN>PINOUT)) THEN
              INOUTV(I+1) = PIN
            ENDIF
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_703
      SUBROUTINE DMUMPS_668(IW, IWSZ, IVAL)
      IMPLICIT NONE
      INTEGER IWSZ
      INTEGER IW(IWSZ)
      INTEGER IVAL
      INTEGER I
      DO I=1,IWSZ
         IW(I)=IVAL
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_668
      SUBROUTINE DMUMPS_704(MYID, NUMPROCS,
     & IRN_loc, JCN_loc, NZ_loc,
     & ROWPARTVEC, COLPARTVEC, M, N,
     & MYROWINDICES, INUMMYR,
     & MYCOLINDICES, INUMMYC,     
     & IWRKROW, IWRKCOL, IWSZR, IWSZC, COMM    )
      IMPLICIT NONE
      INTEGER MYID, NUMPROCS, NZ_loc, M, N
      INTEGER IRN_loc(NZ_loc), JCN_loc(NZ_loc)
      INTEGER ROWPARTVEC(M)
      INTEGER COLPARTVEC(N)
      INTEGER MYROWINDICES(M)
      INTEGER MYCOLINDICES(N)
      INTEGER INUMMYR, INUMMYC
      INTEGER IWSZR, IWSZC
      INTEGER IWRKROW(IWSZR)
      INTEGER IWRKCOL(IWSZC)
      INTEGER COMM
      INTEGER I, IR, IC, ITMP
      INUMMYR = 0
      INUMMYC = 0
      DO I=1,M
         IWRKROW(I) = 0
         IF(ROWPARTVEC(I).EQ.MYID) THEN
            IWRKROW(I)=1
            INUMMYR = INUMMYR + 1
         ENDIF
      ENDDO
      DO I=1,NZ_loc
         IR = IRN_loc(I)
         IC = JCN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.M).AND.
     &        ((IC.GE.1).AND.(IC.LE.N))) THEN
            IF(IWRKROW(IR) .EQ. 0) THEN
               IWRKROW(IR)= 1
               INUMMYR = INUMMYR + 1
            ENDIF
         ENDIF
      ENDDO
      ITMP = 1
      DO I=1,M
         IF(IWRKROW(I).EQ.1) THEN
            MYROWINDICES(ITMP) = I
            ITMP  = ITMP + 1
         ENDIF
      ENDDO
      DO I=1,N
         IWRKCOL(I) = 0
         IF(COLPARTVEC(I).EQ.MYID) THEN
            IWRKCOL(I)= 1
            INUMMYC = INUMMYC + 1
         ENDIF
      ENDDO
      DO I=1,NZ_loc
         IR = IRN_loc(I)
         IC = JCN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.M).AND.
     &        ((IC.GE.1).AND.(IC.LE.N))) THEN
            IF(IWRKCOL(IC) .EQ. 0) THEN
               IWRKCOL(IC)= 1
               INUMMYC = INUMMYC + 1
            ENDIF
         ENDIF
      ENDDO
      ITMP = 1
      DO I=1,N
         IF(IWRKCOL(I).EQ.1) THEN
            MYCOLINDICES(ITMP) = I
            ITMP  = ITMP + 1
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_704
      SUBROUTINE DMUMPS_672(MYID, NUMPROCS, ISZ, IPARTVEC,
     &     NZ_loc, INDX, OSZ, OINDX,ISNDRCVNUM,ISNDRCVVOL,
     &     OSNDRCVNUM,OSNDRCVVOL,
     &     IWRK,IWRKSZ, SNDSZ, RCVSZ, COMM)
      IMPLICIT NONE
      INTEGER MYID, NUMPROCS, NZ_loc, ISZ, IWRKSZ, OSZ
      INTEGER ISNDRCVNUM, ISNDRCVVOL
      INTEGER OSNDRCVNUM, OSNDRCVVOL
      INTEGER COMM
      INTEGER INDX(NZ_loc)
      INTEGER OINDX(NZ_loc)
      INTEGER IPARTVEC(ISZ)
      INTEGER IWRK(IWRKSZ)
      INTEGER SNDSZ(NUMPROCS)
      INTEGER RCVSZ(NUMPROCS)
      INCLUDE 'mpif.h'
      INTEGER I
      INTEGER IIND, IIND2, PIND
      INTEGER IERROR
      DO I=1,NUMPROCS
         SNDSZ(I) = 0
         RCVSZ(I) = 0
      ENDDO
      DO I=1,IWRKSZ
         IWRK(I) = 0
      ENDDO
      DO I=1,NZ_loc
         IIND = INDX(I)
         IIND2 = OINDX(I)
         IF((IIND.GE.1).AND.(IIND.LE.ISZ).AND.
     &        (IIND2.GE.1).AND.(IIND2.LE.OSZ))THEN
            PIND = IPARTVEC(IIND)
            IF(PIND .NE. MYID) THEN
               IF(IWRK(IIND).EQ.0) THEN
                  IWRK(IIND) = 1
                  SNDSZ(PIND+1) = SNDSZ(PIND+1)+1
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      CALL MPI_ALLTOALL(SNDSZ, 1, MPI_INTEGER,
     & RCVSZ, 1, MPI_INTEGER, COMM, IERROR)
      ISNDRCVNUM = 0 
      ISNDRCVVOL = 0
      OSNDRCVNUM = 0
      OSNDRCVVOL = 0
      DO I=1, NUMPROCS
         IF(SNDSZ(I) > 0) OSNDRCVNUM = OSNDRCVNUM + 1
         OSNDRCVVOL = OSNDRCVVOL + SNDSZ(I)
         IF(RCVSZ(I) > 0) ISNDRCVNUM = ISNDRCVNUM + 1
         ISNDRCVVOL = ISNDRCVVOL + RCVSZ(I)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_672
      SUBROUTINE DMUMPS_674(MYID, NUMPROCS, ISZ, IPARTVEC,
     &     NZ_loc, INDX, OSZ, OINDX,
     &     ISNDRCVNUM, ISNDVOL, INGHBPRCS, ISNDRCVIA, ISNDRCVJA,
     &     OSNDRCVNUM, OSNDVOL, ONGHBPRCS, OSNDRCVIA, OSNDRCVJA,
     &     SNDSZ, RCVSZ, IWRK, 
     &     ISTATUS, REQUESTS,
     &     ITAGCOMM, COMM )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MYID, NUMPROCS, NZ_loc, ISZ, ISNDVOL, OSNDVOL, OSZ
      INTEGER INDX(NZ_loc)
      INTEGER OINDX(NZ_loc)
      INTEGER IPARTVEC(ISZ)
      INTEGER ISNDRCVNUM, INGHBPRCS(ISNDRCVNUM)
      INTEGER ISNDRCVIA(NUMPROCS+1)
      INTEGER ISNDRCVJA(ISNDVOL)
      INTEGER OSNDRCVNUM, ONGHBPRCS(OSNDRCVNUM)
      INTEGER OSNDRCVIA(NUMPROCS+1)
      INTEGER OSNDRCVJA(OSNDVOL)
      INTEGER SNDSZ(NUMPROCS)
      INTEGER RCVSZ(NUMPROCS)
      INTEGER IWRK(ISZ)
      INTEGER ISTATUS(MPI_STATUS_SIZE, ISNDRCVNUM)
      INTEGER REQUESTS(ISNDRCVNUM)
      INTEGER ITAGCOMM, COMM
      INTEGER I, IIND, IIND2, IPID, OFFS 
      INTEGER IWHERETO, POFFS, ITMP, IERROR
      DO I=1,ISZ
         IWRK(I) = 0
      ENDDO
      OFFS = 1
      POFFS = 1
      DO I=1,NUMPROCS
         OSNDRCVIA(I) = OFFS + SNDSZ(I)
         IF(SNDSZ(I) > 0) THEN
            ONGHBPRCS(POFFS)=I
            POFFS = POFFS + 1
         ENDIF         
         OFFS  = OFFS +  SNDSZ(I)
      ENDDO
      OSNDRCVIA(NUMPROCS+1) = OFFS
      DO I=1,NZ_loc
         IIND=INDX(I)
         IIND2 = OINDX(I)
         IF((IIND.GE.1).AND.(IIND.LE.ISZ).AND.
     &        (IIND2.GE.1).AND.(IIND2.LE.OSZ) ) THEN
            IPID=IPARTVEC(IIND)
            IF(IPID.NE.MYID) THEN
               IF(IWRK(IIND).EQ.0) THEN
                  IWHERETO = OSNDRCVIA(IPID+1)-1
                  OSNDRCVIA(IPID+1) = OSNDRCVIA(IPID+1)-1
                  OSNDRCVJA(IWHERETO) = IIND
                  IWRK(IIND) = 1
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      CALL MPI_BARRIER(COMM,IERROR)
      OFFS = 1
      POFFS = 1
      ISNDRCVIA(1) = 1
      DO I=2,NUMPROCS+1
         ISNDRCVIA(I) = OFFS + RCVSZ(I-1)
         IF(RCVSZ(I-1) > 0) THEN
            INGHBPRCS(POFFS)=I-1
            POFFS = POFFS + 1
         ENDIF         
         OFFS  = OFFS +  RCVSZ(I-1)
      ENDDO
      CALL MPI_BARRIER(COMM,IERROR)      
      DO I=1, ISNDRCVNUM
         IPID = INGHBPRCS(I)
         OFFS = ISNDRCVIA(IPID)
         ITMP = ISNDRCVIA(IPID+1) - ISNDRCVIA(IPID)
         CALL MPI_IRECV(ISNDRCVJA(OFFS), ITMP, MPI_INTEGER,IPID-1,
     &     ITAGCOMM, COMM, REQUESTS(I),IERROR)   
      ENDDO
      DO I=1,OSNDRCVNUM
         IPID = ONGHBPRCS(I)
         OFFS = OSNDRCVIA(IPID)
         ITMP = OSNDRCVIA(IPID+1)-OSNDRCVIA(IPID)
         CALL MPI_SEND(OSNDRCVJA(OFFS), ITMP, MPI_INTEGER, IPID-1,
     &        ITAGCOMM, COMM,IERROR)
      ENDDO
      IF(ISNDRCVNUM > 0) THEN
         CALL MPI_WAITALL(ISNDRCVNUM, REQUESTS(1),ISTATUS(1,1),IERROR)
      ENDIF
      CALL MPI_BARRIER(COMM,IERROR)
      RETURN
      END SUBROUTINE DMUMPS_674
      SUBROUTINE DMUMPS_657(MYID, NUMPROCS,TMPD, IDSZ, ITAGCOMM, 
     &     ISNDRCVNUM, INGHBPRCS,
     &     ISNDRCVVOL, ISNDRCVIA, ISNDRCVJA, ISNDRCVA,
     &     OSNDRCVNUM, ONGHBPRCS,
     &     OSNDRCVVOL, OSNDRCVIA, OSNDRCVJA, OSNDRCVA,
     &     ISTATUS, REQUESTS,
     &     COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MYID, NUMPROCS, IDSZ, ITAGCOMM
      INTEGER ISNDRCVNUM,OSNDRCVNUM, ISNDRCVVOL, OSNDRCVVOL
      DOUBLE PRECISION TMPD(IDSZ)
      INTEGER INGHBPRCS(ISNDRCVNUM), ONGHBPRCS(OSNDRCVNUM)
      INTEGER ISNDRCVIA(NUMPROCS+1), ISNDRCVJA(ISNDRCVVOL)
      DOUBLE PRECISION ISNDRCVA(ISNDRCVVOL)
      INTEGER OSNDRCVIA(NUMPROCS+1), OSNDRCVJA(OSNDRCVVOL)
      DOUBLE PRECISION OSNDRCVA(OSNDRCVVOL)
      INTEGER ISTATUS(MPI_STATUS_SIZE, max(ISNDRCVNUM,OSNDRCVNUM))
      INTEGER REQUESTS(max(ISNDRCVNUM,OSNDRCVNUM))
      INTEGER COMM, IERROR
      INTEGER I, PID, OFFS, SZ, J, JS, JE, IID
      DO I=1,ISNDRCVNUM
         PID = INGHBPRCS(I)
         OFFS = ISNDRCVIA(PID)
         SZ = ISNDRCVIA(PID+1) -  ISNDRCVIA(PID) 
         CALL MPI_IRECV(ISNDRCVA(OFFS), SZ, 
     &        MPI_DOUBLE_PRECISION, PID-1, 
     &        ITAGCOMM,COMM,REQUESTS(I), IERROR)
      ENDDO
      DO I=1,OSNDRCVNUM
         PID = ONGHBPRCS(I)
         OFFS = OSNDRCVIA(PID)
         SZ = OSNDRCVIA(PID+1) - OSNDRCVIA(PID) 
         JS = OSNDRCVIA(PID) 
         JE =  OSNDRCVIA(PID+1) - 1
         DO J=JS, JE
            IID = OSNDRCVJA(J)
            OSNDRCVA(J) = TMPD(IID)
         ENDDO
         CALL MPI_SEND(OSNDRCVA(OFFS), SZ, MPI_DOUBLE_PRECISION, PID-1,
     &        ITAGCOMM, COMM, IERROR)
      ENDDO
      IF(ISNDRCVNUM > 0) THEN
         CALL MPI_WAITALL(ISNDRCVNUM, REQUESTS(1),ISTATUS(1,1),IERROR)
      ENDIF
      DO I=1,ISNDRCVNUM
         PID = INGHBPRCS(I)
         JS = ISNDRCVIA(PID)
         JE = ISNDRCVIA(PID+1)-1
         DO J=JS,JE
            IID = ISNDRCVJA(J)
            IF(TMPD(IID) < ISNDRCVA(J)) TMPD(IID)= ISNDRCVA(J)
         ENDDO
      ENDDO
      DO I=1,OSNDRCVNUM
         PID = ONGHBPRCS(I)
         OFFS = OSNDRCVIA(PID)
         SZ = OSNDRCVIA(PID+1) -  OSNDRCVIA(PID) 
         CALL MPI_IRECV(OSNDRCVA(OFFS), SZ, 
     &        MPI_DOUBLE_PRECISION, PID-1, 
     &        ITAGCOMM+1,COMM,REQUESTS(I), IERROR)
      ENDDO
      DO I=1,ISNDRCVNUM
         PID = INGHBPRCS(I)
         OFFS = ISNDRCVIA(PID)
         SZ = ISNDRCVIA(PID+1)-ISNDRCVIA(PID)
         JS = ISNDRCVIA(PID)
         JE = ISNDRCVIA(PID+1) -1
         DO J=JS, JE
            IID = ISNDRCVJA(J)
            ISNDRCVA(J) = TMPD(IID)
         ENDDO
         CALL MPI_SEND(ISNDRCVA(OFFS), SZ, MPI_DOUBLE_PRECISION, PID-1,
     &        ITAGCOMM+1, COMM, IERROR)
      ENDDO
      IF(OSNDRCVNUM > 0) THEN
         CALL MPI_WAITALL(OSNDRCVNUM, REQUESTS(1),ISTATUS(1,1),IERROR)
      ENDIF
      DO I=1,OSNDRCVNUM
         PID = ONGHBPRCS(I)
         JS = OSNDRCVIA(PID) 
         JE = OSNDRCVIA(PID+1) - 1
         DO J=JS,JE
            IID = OSNDRCVJA(J)
            TMPD(IID)=OSNDRCVA(J)
         ENDDO
      ENDDO
      RETURN
      END  SUBROUTINE DMUMPS_657
      SUBROUTINE DMUMPS_656(MYID, NUMPROCS,TMPD, IDSZ, ITAGCOMM, 
     &     ISNDRCVNUM, INGHBPRCS,
     &     ISNDRCVVOL, ISNDRCVIA, ISNDRCVJA, ISNDRCVA,
     &     OSNDRCVNUM, ONGHBPRCS,
     &     OSNDRCVVOL, OSNDRCVIA, OSNDRCVJA, OSNDRCVA,
     &     ISTATUS, REQUESTS,
     &     COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MYID, NUMPROCS, IDSZ, ITAGCOMM
      INTEGER ISNDRCVNUM,OSNDRCVNUM, ISNDRCVVOL, OSNDRCVVOL
      DOUBLE PRECISION TMPD(IDSZ)
      INTEGER INGHBPRCS(ISNDRCVNUM), ONGHBPRCS(OSNDRCVNUM)
      INTEGER ISNDRCVIA(NUMPROCS+1), ISNDRCVJA(ISNDRCVVOL)
      DOUBLE PRECISION ISNDRCVA(ISNDRCVVOL)
      INTEGER OSNDRCVIA(NUMPROCS+1), OSNDRCVJA(OSNDRCVVOL)
      DOUBLE PRECISION OSNDRCVA(OSNDRCVVOL)
      INTEGER ISTATUS(MPI_STATUS_SIZE, max(ISNDRCVNUM,OSNDRCVNUM))
      INTEGER REQUESTS(max(ISNDRCVNUM,OSNDRCVNUM))
      INTEGER COMM, IERROR
      INTEGER I, PID, OFFS, SZ, J, JS, JE, IID
      DO I=1,ISNDRCVNUM
         PID = INGHBPRCS(I)
         OFFS = ISNDRCVIA(PID)
         SZ = ISNDRCVIA(PID+1) -  ISNDRCVIA(PID) 
         CALL MPI_IRECV(ISNDRCVA(OFFS), SZ, 
     &        MPI_DOUBLE_PRECISION, PID-1, 
     &        ITAGCOMM,COMM,REQUESTS(I), IERROR)
      ENDDO
      DO I=1,OSNDRCVNUM
         PID = ONGHBPRCS(I)
         OFFS = OSNDRCVIA(PID)
         SZ = OSNDRCVIA(PID+1) - OSNDRCVIA(PID) 
         JS = OSNDRCVIA(PID) 
         JE =  OSNDRCVIA(PID+1) - 1
         DO J=JS, JE
            IID = OSNDRCVJA(J)
            OSNDRCVA(J) = TMPD(IID)
         ENDDO
         CALL MPI_SEND(OSNDRCVA(OFFS), SZ, MPI_DOUBLE_PRECISION, PID-1,
     &        ITAGCOMM, COMM, IERROR)
      ENDDO
      IF(ISNDRCVNUM > 0) THEN
         CALL MPI_WAITALL(ISNDRCVNUM, REQUESTS(1),ISTATUS(1,1),IERROR)
      ENDIF
      DO I=1,ISNDRCVNUM
         PID = INGHBPRCS(I)
         JS = ISNDRCVIA(PID)
         JE = ISNDRCVIA(PID+1)-1
         DO J=JS,JE
            IID = ISNDRCVJA(J)
            TMPD(IID)  = TMPD(IID)+ ISNDRCVA(J)
         ENDDO
      ENDDO
      DO I=1,OSNDRCVNUM
         PID = ONGHBPRCS(I)
         OFFS = OSNDRCVIA(PID)
         SZ = OSNDRCVIA(PID+1) -  OSNDRCVIA(PID) 
         CALL MPI_IRECV(OSNDRCVA(OFFS), SZ, 
     &        MPI_DOUBLE_PRECISION, PID-1, 
     &        ITAGCOMM+1,COMM,REQUESTS(I), IERROR)
      ENDDO
      DO I=1,ISNDRCVNUM
         PID = INGHBPRCS(I)
         OFFS = ISNDRCVIA(PID)
         SZ = ISNDRCVIA(PID+1)-ISNDRCVIA(PID)
         JS = ISNDRCVIA(PID)
         JE = ISNDRCVIA(PID+1) -1
         DO J=JS, JE
            IID = ISNDRCVJA(J)
            ISNDRCVA(J) = TMPD(IID)
         ENDDO
         CALL MPI_SEND(ISNDRCVA(OFFS), SZ, MPI_DOUBLE_PRECISION, PID-1,
     &        ITAGCOMM+1, COMM, IERROR)
      ENDDO
      IF(OSNDRCVNUM > 0) THEN
         CALL MPI_WAITALL(OSNDRCVNUM, REQUESTS(1),ISTATUS(1,1),IERROR)
      ENDIF
      DO I=1,OSNDRCVNUM
         PID = ONGHBPRCS(I)
         JS = OSNDRCVIA(PID) 
         JE = OSNDRCVIA(PID+1) - 1
         DO J=JS,JE
            IID = OSNDRCVJA(J)
            TMPD(IID)=OSNDRCVA(J)
         ENDDO
      ENDDO
      RETURN
      END  SUBROUTINE DMUMPS_656
      SUBROUTINE DMUMPS_655(MYID, NUMPROCS, COMM,
     & IRN_loc, JCN_loc, NZ_loc,
     & IPARTVEC, ISZ,
     & IWRK, IWSZ)
      IMPLICIT NONE
      EXTERNAL DMUMPS_703
      INTEGER MYID, NUMPROCS, COMM
      INTEGER NZ_loc, ISZ, IWSZ
      INTEGER IRN_loc(NZ_loc), JCN_loc(NZ_loc)
      INTEGER IPARTVEC(ISZ)
      INTEGER IWRK(IWSZ)
      INCLUDE 'mpif.h'
      INTEGER I
      INTEGER OP, IERROR
      INTEGER IR, IC
      IF(NUMPROCS.NE.1) THEN
         CALL MPI_OP_CREATE(DMUMPS_703, .TRUE., OP, IERROR)
         CALL DMUMPS_668(IWRK, 4*ISZ, ISZ)
         DO I=1,ISZ
            IWRK(2*I-1) = 0
            IWRK(2*I) = MYID
         ENDDO
         DO I=1,NZ_loc
            IR = IRN_loc(I)
            IC = JCN_loc(I)
            IF((IR.GE.1).AND.(IR.LE.ISZ).AND.
     &           (IC.GE.1).AND.(IC.LE.ISZ)) THEN
               IWRK(2*IR-1) = IWRK(2*IR-1) + 1
               IWRK(2*IC-1) = IWRK(2*IC-1) + 1
            ENDIF
         ENDDO
         CALL MPI_ALLREDUCE(IWRK(1), IWRK(1+2*ISZ), ISZ,
     &        MPI_2INTEGER, OP, COMM, IERROR)      
         DO I=1,ISZ
            IPARTVEC(I) = IWRK(2*I+2*ISZ)
         ENDDO
         CALL MPI_OP_FREE(OP, IERROR)
      ELSE
         DO I=1,ISZ
            IPARTVEC(I) = 0
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_655
      SUBROUTINE DMUMPS_673(MYID, NUMPROCS, ISZ, IPARTVEC,
     & NZ_loc, INDX,OINDX,ISNDRCVNUM,ISNDRCVVOL,OSNDRCVNUM,OSNDRCVVOL,
     & IWRK,IWRKSZ, SNDSZ, RCVSZ, COMM)
      IMPLICIT NONE
      INTEGER MYID, NUMPROCS, NZ_loc, ISZ, IWRKSZ
      INTEGER ISNDRCVNUM, ISNDRCVVOL
      INTEGER OSNDRCVNUM, OSNDRCVVOL
      INTEGER COMM
      INTEGER INDX(NZ_loc), OINDX(NZ_loc)
      INTEGER IPARTVEC(ISZ)
      INTEGER IWRK(IWRKSZ)
      INTEGER SNDSZ(NUMPROCS)
      INTEGER RCVSZ(NUMPROCS)
      INCLUDE 'mpif.h'
      INTEGER I
      INTEGER IIND, IIND2, PIND
      INTEGER IERROR
      DO I=1,NUMPROCS
         SNDSZ(I) = 0
         RCVSZ(I) = 0
      ENDDO
      DO I=1,IWRKSZ
         IWRK(I) = 0
      ENDDO
      DO I=1,NZ_loc
         IIND = INDX(I)
         IIND2 = OINDX(I)
         IF((IIND.GE.1).AND.(IIND.LE.ISZ).AND.(IIND2.GE.1)
     &        .AND.(IIND2.LE.ISZ)) THEN
            PIND = IPARTVEC(IIND)
            IF(PIND .NE. MYID) THEN
               IF(IWRK(IIND).EQ.0) THEN
                  IWRK(IIND) = 1
                  SNDSZ(PIND+1) = SNDSZ(PIND+1)+1
               ENDIF
            ENDIF
            IIND = OINDX(I)
            PIND = IPARTVEC(IIND)
            IF(PIND .NE. MYID) THEN
               IF(IWRK(IIND).EQ.0) THEN
                  IWRK(IIND) = 1
                  SNDSZ(PIND+1) = SNDSZ(PIND+1)+1
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      CALL MPI_ALLTOALL(SNDSZ, 1, MPI_INTEGER,
     &     RCVSZ, 1, MPI_INTEGER, COMM, IERROR)
      ISNDRCVNUM = 0 
      ISNDRCVVOL = 0
      OSNDRCVNUM = 0
      OSNDRCVVOL = 0
      DO I=1, NUMPROCS
         IF(SNDSZ(I) > 0) OSNDRCVNUM = OSNDRCVNUM + 1
         OSNDRCVVOL = OSNDRCVVOL + SNDSZ(I)
         IF(RCVSZ(I) > 0) ISNDRCVNUM = ISNDRCVNUM + 1
         ISNDRCVVOL = ISNDRCVVOL + RCVSZ(I)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_673
      SUBROUTINE DMUMPS_663(MYID, NUMPROCS, COMM,
     &     IRN_loc, JCN_loc, NZ_loc,
     &     PARTVEC, N,
     &     INUMMYR,
     &     IWRK, IWSZ)
      IMPLICIT NONE
      INTEGER MYID, NUMPROCS, NZ_loc, N
      INTEGER IRN_loc(NZ_loc), JCN_loc(NZ_loc)
      INTEGER PARTVEC(N)
      INTEGER INUMMYR
      INTEGER IWSZ
      INTEGER IWRK(IWSZ)
      INTEGER COMM
      INTEGER I, IR, IC
      INUMMYR = 0
      DO I=1,N
         IWRK(I) = 0
         IF(PARTVEC(I).EQ.MYID) THEN
            IWRK(I)=1
            INUMMYR = INUMMYR + 1
         ENDIF
      ENDDO
      DO I=1,NZ_loc
         IR = IRN_loc(I)
         IC = JCN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.N).AND.
     &        ((IC.GE.1).AND.(IC.LE.N))) THEN
            IF(IWRK(IR) .EQ. 0) THEN
               IWRK(IR)= 1
               INUMMYR = INUMMYR + 1
            ENDIF
         ENDIF
         IF((IR.GE.1).AND.(IR.LE.N).AND.
     &        ((IC.GE.1).AND.(IC.LE.N))) THEN
            IF(IWRK(IC).EQ.0) THEN
               IWRK(IC)= 1
               INUMMYR = INUMMYR + 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_663
      INTEGER FUNCTION DMUMPS_742(D, N, INDXR, INDXRSZ,
     &     EPS, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER N, INDXRSZ
      DOUBLE PRECISION D(N)
      INTEGER INDXR(INDXRSZ)
      DOUBLE PRECISION EPS
      INTEGER COMM
      EXTERNAL DMUMPS_744
      INTEGER  DMUMPS_744
      INTEGER GLORES, MYRESR, MYRES
      INTEGER IERR
      MYRESR =  DMUMPS_744(D, N, INDXR, INDXRSZ, EPS)
      MYRES = 2*MYRESR 
      CALL MPI_ALLREDUCE(MYRES, GLORES, 1, MPI_INTEGER,
     &     MPI_SUM, COMM, IERR)
      DMUMPS_742 = GLORES
      RETURN
      END FUNCTION DMUMPS_742
      SUBROUTINE DMUMPS_661(MYID, NUMPROCS,COMM,    
     &     IRN_loc, JCN_loc, NZ_loc,
     &     PARTVEC, N,
     &     MYROWINDICES, INUMMYR,
     &     IWRK, IWSZ  )
      IMPLICIT NONE
      INTEGER MYID, NUMPROCS, NZ_loc, N
      INTEGER INUMMYR, IWSZ
      INTEGER IRN_loc(NZ_loc), JCN_loc(NZ_loc)
      INTEGER PARTVEC(N)
      INTEGER MYROWINDICES(INUMMYR)
      INTEGER IWRK(IWSZ)
      INTEGER COMM
      INTEGER I, IR, IC, ITMP, MAXMN
      MAXMN = N
      DO I=1,N
         IWRK(I) = 0
         IF(PARTVEC(I).EQ.MYID) IWRK(I)=1
      ENDDO
      DO I=1,NZ_loc
         IR = IRN_loc(I)   
         IC = JCN_loc(I)
         IF((IR.GE.1).AND.(IR.LE.N).AND.
     &        ((IC.GE.1).AND.(IC.LE.N))) THEN
            IF(IWRK(IR) .EQ. 0) IWRK(IR)= 1
         ENDIF
         IF((IR.GE.1).AND.(IR.LE.N).AND.
     &        ((IC.GE.1).AND.(IC.LE.N))) THEN
            IF(IWRK(IC) .EQ.0) IWRK(IC)=1
         ENDIF
      ENDDO
      ITMP = 1
      DO I=1,N
         IF(IWRK(I).EQ.1) THEN
            MYROWINDICES(ITMP) = I
            ITMP  = ITMP + 1
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_661
      SUBROUTINE DMUMPS_692(MYID, NUMPROCS, ISZ, IPARTVEC,
     & NZ_loc, INDX, OINDX,
     & ISNDRCVNUM, ISNDVOL, INGHBPRCS, ISNDRCVIA, ISNDRCVJA,
     & OSNDRCVNUM, OSNDVOL, ONGHBPRCS, OSNDRCVIA, OSNDRCVJA,
     & SNDSZ, RCVSZ, IWRK, 
     & ISTATUS, REQUESTS,
     &  ITAGCOMM, COMM )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MYID, NUMPROCS, NZ_loc, ISZ, ISNDVOL, OSNDVOL
      INTEGER INDX(NZ_loc), OINDX(NZ_loc)
      INTEGER IPARTVEC(ISZ)
      INTEGER ISNDRCVNUM, INGHBPRCS(ISNDRCVNUM)
      INTEGER ISNDRCVIA(NUMPROCS+1)
      INTEGER ISNDRCVJA(ISNDVOL)
      INTEGER OSNDRCVNUM, ONGHBPRCS(OSNDRCVNUM)
      INTEGER OSNDRCVIA(NUMPROCS+1)
      INTEGER OSNDRCVJA(OSNDVOL)
      INTEGER SNDSZ(NUMPROCS)
      INTEGER RCVSZ(NUMPROCS)
      INTEGER IWRK(ISZ)
      INTEGER ISTATUS(MPI_STATUS_SIZE, ISNDRCVNUM)
      INTEGER REQUESTS(ISNDRCVNUM)
      INTEGER ITAGCOMM, COMM
      INTEGER I, IIND,IIND2,IPID,OFFS,IWHERETO,POFFS, ITMP, IERROR
      DO I=1,ISZ
         IWRK(I) = 0
      ENDDO
      OFFS = 1
      POFFS = 1
      DO I=1,NUMPROCS
         OSNDRCVIA(I) = OFFS + SNDSZ(I)
         IF(SNDSZ(I) > 0) THEN
            ONGHBPRCS(POFFS)=I
            POFFS = POFFS + 1
         ENDIF         
         OFFS  = OFFS +  SNDSZ(I)
      ENDDO
      OSNDRCVIA(NUMPROCS+1) = OFFS
      DO I=1,NZ_loc
         IIND=INDX(I)
         IIND2 = OINDX(I)
         IF((IIND.GE.1).AND.(IIND.LE.ISZ).AND.(IIND2.GE.1)
     &        .AND.(IIND2.LE.ISZ)) THEN            
            IPID=IPARTVEC(IIND)
            IF(IPID.NE.MYID) THEN
               IF(IWRK(IIND).EQ.0) THEN
                  IWHERETO = OSNDRCVIA(IPID+1)-1
                  OSNDRCVIA(IPID+1) = OSNDRCVIA(IPID+1)-1
                  OSNDRCVJA(IWHERETO) = IIND
                  IWRK(IIND) = 1
               ENDIF
            ENDIF
            IIND = OINDX(I)
            IPID=IPARTVEC(IIND)
            IF(IPID.NE.MYID) THEN
               IF(IWRK(IIND).EQ.0) THEN
                  IWHERETO = OSNDRCVIA(IPID+1)-1
                  OSNDRCVIA(IPID+1) = OSNDRCVIA(IPID+1)-1
                  OSNDRCVJA(IWHERETO) = IIND
                  IWRK(IIND) = 1
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      CALL MPI_BARRIER(COMM,IERROR)
      OFFS = 1
      POFFS = 1
      ISNDRCVIA(1) = 1
      DO I=2,NUMPROCS+1
         ISNDRCVIA(I) = OFFS + RCVSZ(I-1)
         IF(RCVSZ(I-1) > 0) THEN
            INGHBPRCS(POFFS)=I-1
            POFFS = POFFS + 1
         ENDIF         
         OFFS  = OFFS +  RCVSZ(I-1)
      ENDDO
      CALL MPI_BARRIER(COMM,IERROR)      
      DO I=1, ISNDRCVNUM
         IPID = INGHBPRCS(I)
         OFFS = ISNDRCVIA(IPID)
         ITMP = ISNDRCVIA(IPID+1) - ISNDRCVIA(IPID)
         CALL MPI_IRECV(ISNDRCVJA(OFFS), ITMP, MPI_INTEGER,IPID-1,
     &     ITAGCOMM, COMM, REQUESTS(I),IERROR)   
      ENDDO
      DO I=1,OSNDRCVNUM
         IPID = ONGHBPRCS(I)
         OFFS = OSNDRCVIA(IPID)
         ITMP = OSNDRCVIA(IPID+1)-OSNDRCVIA(IPID)
         CALL MPI_SEND(OSNDRCVJA(OFFS), ITMP, MPI_INTEGER, IPID-1,
     &        ITAGCOMM, COMM,IERROR)
      ENDDO
      IF(ISNDRCVNUM > 0) THEN
         CALL MPI_WAITALL(ISNDRCVNUM, REQUESTS(1),ISTATUS(1,1),IERROR)
      ENDIF
      CALL MPI_BARRIER(COMM,IERROR)
      RETURN
      END SUBROUTINE DMUMPS_692
      SUBROUTINE DMUMPS_628(IW,LREC,SIZE_FREE,XSIZE)
      INTEGER, intent(in) :: LREC, XSIZE
      INTEGER, intent(in) :: IW(LREC)
      INTEGER(8), intent(out):: SIZE_FREE
      INCLUDE 'mumps_headers.h'
      IF (IW(1+XXS).EQ.S_NOLCBCONTIG .OR.
     &    IW(1+XXS).EQ.S_NOLCBNOCONTIG) THEN
        SIZE_FREE=int(IW(1+XSIZE+2),8)*int(IW(1+XSIZE+3),8)
      ELSE IF (IW(1+XXS).EQ.S_NOLCBCONTIG38 .OR.
     &         IW(1+XXS).EQ.S_NOLCBNOCONTIG38) THEN
        SIZE_FREE=int(IW(1+XSIZE+2),8)*int(IW(1+XSIZE)+
     &            IW(1+XSIZE + 3) -
     &          ( IW(1+XSIZE + 4)
     &          - IW(1+XSIZE + 3) ), 8)
      ELSE
        SIZE_FREE=0_8
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_628
      SUBROUTINE DMUMPS_629
     &(IW,LIW,IXXP,ICURRENT,NEXT, RCURRENT,ISIZE2SHIFT)
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: RCURRENT
      INTEGER LIW,IXXP,ICURRENT,NEXT,ISIZE2SHIFT
      INTEGER IW(LIW)
      INTEGER(8) :: RSIZE
      ICURRENT=NEXT
      CALL MUMPS_729( RSIZE, IW(ICURRENT + XXR) )
      RCURRENT = RCURRENT - RSIZE
      NEXT=IW(ICURRENT+XXP)
      IW(IXXP)=ICURRENT+ISIZE2SHIFT
      IXXP=ICURRENT+XXP
      RETURN
      END SUBROUTINE DMUMPS_629
      SUBROUTINE DMUMPS_630(IW,LIW,BEG2SHIFT,END2SHIFT,ISIZE2SHIFT)
      IMPLICIT NONE
      INTEGER LIW, BEG2SHIFT, END2SHIFT, ISIZE2SHIFT
      INTEGER IW(LIW)
      INTEGER I
      IF (ISIZE2SHIFT.GT.0) THEN
        DO I=END2SHIFT,BEG2SHIFT,-1
          IW(I+ISIZE2SHIFT)=IW(I)
        ENDDO
      ELSE IF (ISIZE2SHIFT.LT.0) THEN
        DO I=BEG2SHIFT,END2SHIFT
          IW(I+ISIZE2SHIFT)=IW(I)
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_630
      SUBROUTINE DMUMPS_631(A, LA, BEG2SHIFT, END2SHIFT, RSIZE2SHIFT)
      IMPLICIT NONE
      INTEGER(8) :: LA, BEG2SHIFT, END2SHIFT, RSIZE2SHIFT
      DOUBLE PRECISION A(LA)
      INTEGER(8) :: I
      IF (RSIZE2SHIFT.GT.0_8) THEN
        DO I=END2SHIFT,BEG2SHIFT,-1_8
          A(I+RSIZE2SHIFT)=A(I)
        ENDDO
      ELSE IF (RSIZE2SHIFT.LT.0_8) THEN
        DO I=BEG2SHIFT,END2SHIFT
          A(I+RSIZE2SHIFT)=A(I)
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_631
      SUBROUTINE DMUMPS_94(N,KEEP28,IW,LIW,A,LA,
     &       LRLU,IPTRLU,IWPOS,
     &       IWPOSCB,PTRIST,PTRAST,STEP,PIMASTER,PAMASTER,
     &       KEEP216,LRLUS,XSIZE)
      IMPLICIT NONE
      INTEGER(8) :: LA, LRLU, IPTRLU, LRLUS
      INTEGER N,LIW,KEEP28,
     &        IWPOS,IWPOSCB,KEEP216,XSIZE
      INTEGER(8) :: PTRAST(KEEP28), PAMASTER(KEEP28)
      INTEGER IW(LIW),PTRIST(KEEP28),
     &        STEP(N), PIMASTER(KEEP28)
      DOUBLE PRECISION A(LA)
      INCLUDE 'mumps_headers.h' 
      INTEGER ICURRENT, NEXT, STATE_NEXT
      INTEGER(8) :: RCURRENT
      INTEGER ISIZE2SHIFT
      INTEGER(8) :: RSIZE2SHIFT
      INTEGER IBEGCONTIG
      INTEGER(8) :: RBEGCONTIG
      INTEGER(8) :: RBEG2SHIFT, REND2SHIFT
      INTEGER INODE
      INTEGER(8) :: FREE_IN_REC
      INTEGER(8) :: RCURRENT_SIZE
      INTEGER IXXP
      ISIZE2SHIFT=0
      RSIZE2SHIFT=0_8
      ICURRENT  = LIW-XSIZE+1
      RCURRENT = LA+1_8
      IBEGCONTIG = -999999 
      RBEGCONTIG = -999999_8 
      NEXT = IW(ICURRENT+XXP)
      IF (NEXT.EQ.TOP_OF_STACK) RETURN
      STATE_NEXT = IW(NEXT+XXS)
      IXXP = ICURRENT+XXP
  10     CONTINUE
         IF ( STATE_NEXT .NE. S_FREE .AND.
     &        (KEEP216.EQ.3.OR.
     &         (STATE_NEXT .NE. S_NOLCBNOCONTIG .AND.
     &          STATE_NEXT .NE. S_NOLCBCONTIG .AND.
     &          STATE_NEXT .NE. S_NOLCBNOCONTIG38 .AND.
     &          STATE_NEXT .NE. S_NOLCBCONTIG38))) THEN
            CALL DMUMPS_629(IW,LIW,
     &           IXXP, ICURRENT, NEXT, RCURRENT, ISIZE2SHIFT)
            CALL MUMPS_729(RCURRENT_SIZE, IW(ICURRENT+XXR))
            IF (IBEGCONTIG < 0) THEN
              IBEGCONTIG=ICURRENT+IW(ICURRENT+XXI)-1
            ENDIF
            IF (RBEGCONTIG < 0_8) THEN
              RBEGCONTIG=RCURRENT+RCURRENT_SIZE-1_8
            ENDIF
            INODE=IW(ICURRENT+XXN)
            IF (RSIZE2SHIFT .NE. 0_8) THEN
                IF (PTRAST(STEP(INODE)).EQ.RCURRENT)
     &            PTRAST(STEP(INODE))=
     &            PTRAST(STEP(INODE))+RSIZE2SHIFT
                IF (PAMASTER(STEP(INODE)).EQ.RCURRENT)
     &            PAMASTER(STEP(INODE))=
     &            PAMASTER(STEP(INODE))+RSIZE2SHIFT
            ENDIF
            IF (ISIZE2SHIFT .NE. 0) THEN
                IF (PTRIST(STEP(INODE)).EQ.ICURRENT)
     &            PTRIST(STEP(INODE))=
     &            PTRIST(STEP(INODE))+ISIZE2SHIFT
                IF (PIMASTER(STEP(INODE)).EQ.ICURRENT)
     &            PIMASTER(STEP(INODE))=
     &            PIMASTER(STEP(INODE))+ISIZE2SHIFT
            ENDIF
            IF (NEXT .NE. TOP_OF_STACK) THEN
              STATE_NEXT=IW(NEXT+XXS)
              GOTO 10
            ENDIF
         ENDIF
  20     CONTINUE
         IF (IBEGCONTIG.NE.0 .AND. ISIZE2SHIFT .NE. 0) THEN
           CALL DMUMPS_630(IW,LIW,ICURRENT,IBEGCONTIG,ISIZE2SHIFT)
           IF (IXXP .LE.IBEGCONTIG) THEN
           IXXP=IXXP+ISIZE2SHIFT
           ENDIF
         ENDIF
         IBEGCONTIG=-9999
  25     CONTINUE
         IF (RBEGCONTIG .GT.0_8 .AND. RSIZE2SHIFT .NE. 0_8) THEN
           CALL DMUMPS_631(A,LA,RCURRENT,RBEGCONTIG,RSIZE2SHIFT)
         ENDIF
         RBEGCONTIG=-99999_8
  30     CONTINUE
         IF (NEXT.EQ. TOP_OF_STACK) GOTO 100
         IF (STATE_NEXT .EQ. S_NOLCBCONTIG .OR.
     &       STATE_NEXT .EQ. S_NOLCBNOCONTIG .OR.
     &       STATE_NEXT .EQ. S_NOLCBCONTIG38 .OR.
     &       STATE_NEXT .EQ. S_NOLCBNOCONTIG38) THEN
           IF ( KEEP216.eq.3) THEN
             WRITE(*,*) "Internal error 2 in DMUMPS_94"
           ENDIF
           IF (RBEGCONTIG > 0_8) GOTO 25
           CALL DMUMPS_629
     &       (IW,LIW,IXXP,ICURRENT,NEXT, RCURRENT,ISIZE2SHIFT)
           IF (IBEGCONTIG < 0 ) THEN
             IBEGCONTIG=ICURRENT+IW(ICURRENT+XXI)-1
           ENDIF
           CALL DMUMPS_628(IW(ICURRENT),
     &                              LIW-ICURRENT+1,
     &                              FREE_IN_REC,
     &                              XSIZE)
           IF (STATE_NEXT .EQ. S_NOLCBNOCONTIG) THEN
             CALL DMUMPS_627(A,LA,RCURRENT,
     &            IW(ICURRENT+XSIZE+2),
     &            IW(ICURRENT+XSIZE),
     &            IW(ICURRENT+XSIZE)+IW(ICURRENT+XSIZE+3), 0,
     &            IW(ICURRENT+XXS),RSIZE2SHIFT) 
           ELSE IF (STATE_NEXT .EQ. S_NOLCBNOCONTIG38) THEN
             CALL DMUMPS_627(A,LA,RCURRENT,
     &            IW(ICURRENT+XSIZE+2),
     &            IW(ICURRENT+XSIZE),
     &            IW(ICURRENT+XSIZE)+IW(ICURRENT+XSIZE+3),
     &            IW(ICURRENT+XSIZE+4)-IW(ICURRENT+XSIZE+3), 
     &            IW(ICURRENT+XXS),RSIZE2SHIFT) 
           ELSE IF (RSIZE2SHIFT .GT.0_8) THEN
             RBEG2SHIFT = RCURRENT + FREE_IN_REC
             CALL MUMPS_729(RCURRENT_SIZE, IW(ICURRENT+XXR))
             REND2SHIFT = RCURRENT + RCURRENT_SIZE - 1_8
             CALL DMUMPS_631(A, LA,
     &                          RBEG2SHIFT, REND2SHIFT,
     &                          RSIZE2SHIFT)
           ENDIF
           INODE=IW(ICURRENT+XXN)
           IF (ISIZE2SHIFT.NE.0) THEN
             PTRIST(STEP(INODE))=PTRIST(STEP(INODE))+ISIZE2SHIFT
           ENDIF
           PTRAST(STEP(INODE))=PTRAST(STEP(INODE))+RSIZE2SHIFT+
     &                         FREE_IN_REC
           CALL MUMPS_724(IW(ICURRENT+XXR),FREE_IN_REC)
           IF (STATE_NEXT.EQ.S_NOLCBCONTIG.OR.
     &         STATE_NEXT.EQ.S_NOLCBNOCONTIG) THEN
             IW(ICURRENT+XXS)=S_NOLCLEANED
           ELSE
             IW(ICURRENT+XXS)=S_NOLCLEANED38
           ENDIF
           RSIZE2SHIFT=RSIZE2SHIFT+FREE_IN_REC
           RBEGCONTIG=-9999_8
           IF (NEXT.EQ.TOP_OF_STACK) THEN
             GOTO 20
           ELSE
             STATE_NEXT=IW(NEXT+XXS)
           ENDIF
           GOTO 30
         ENDIF
         IF (IBEGCONTIG.GT.0) THEN
           GOTO 20
         ENDIF
  40     CONTINUE
         IF (STATE_NEXT == S_FREE) THEN
            ICURRENT = NEXT
            CALL MUMPS_729( RCURRENT_SIZE, IW(ICURRENT + XXR) )
            ISIZE2SHIFT = ISIZE2SHIFT + IW(ICURRENT+XXI)
            RSIZE2SHIFT = RSIZE2SHIFT + RCURRENT_SIZE
            RCURRENT    = RCURRENT    - RCURRENT_SIZE
            NEXT=IW(ICURRENT+XXP)
            IF (NEXT.EQ.TOP_OF_STACK) THEN
              WRITE(*,*) "Internal error 1 in DMUMPS_94"
              CALL MUMPS_ABORT()
            ENDIF
            STATE_NEXT  = IW(NEXT+XXS)
            GOTO 40
         ENDIF
      GOTO 10
 100  CONTINUE
      IWPOSCB = IWPOSCB + ISIZE2SHIFT
      LRLU    = LRLU    + RSIZE2SHIFT
      IPTRLU  = IPTRLU  + RSIZE2SHIFT
      RETURN
      END SUBROUTINE DMUMPS_94
      SUBROUTINE DMUMPS_632(IREC, IW, LIW,
     &            ISIZEHOLE, RSIZEHOLE)
      IMPLICIT NONE
      INTEGER, intent(in) :: IREC, LIW
      INTEGER, intent(in) :: IW(LIW)
      INTEGER, intent(out):: ISIZEHOLE
      INTEGER(8), intent(out) :: RSIZEHOLE
      INTEGER IRECLOC
      INTEGER(8) :: RECLOC_SIZE
      INCLUDE 'mumps_headers.h'
      ISIZEHOLE=0
      RSIZEHOLE=0_8
      IRECLOC = IREC + IW( IREC+XXI )
 10   CONTINUE
      CALL MUMPS_729(RECLOC_SIZE, IW(IRECLOC+XXR))
      IF (IW(IRECLOC+XXS).EQ.S_FREE) THEN
        ISIZEHOLE=ISIZEHOLE+IW(IRECLOC+XXI)
        RSIZEHOLE=RSIZEHOLE+RECLOC_SIZE
        IRECLOC=IRECLOC+IW(IRECLOC+XXI)
        GOTO 10
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_632
      SUBROUTINE DMUMPS_627(A, LA, RCURRENT,
     &           NROW, NCB, LD, NELIM, NODESTATE, ISHIFT)
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER LD, NROW, NCB, NELIM, NODESTATE
      INTEGER(8) :: ISHIFT
      INTEGER(8) :: LA, RCURRENT
      DOUBLE PRECISION A(LA)
      INTEGER I,J
      INTEGER(8) :: IOLD,INEW
      LOGICAL NELIM_ROOT
      NELIM_ROOT=.TRUE.
      IF (NODESTATE.EQ. S_NOLCBNOCONTIG) THEN
         NELIM_ROOT=.FALSE.
         IF (NELIM.NE.0)  THEN
           WRITE(*,*) "Internal error 1 IN DMUMPS_627"
           CALL MUMPS_ABORT()
         ENDIF
      ELSE IF (NODESTATE .NE. S_NOLCBNOCONTIG38) THEN
           WRITE(*,*) "Internal error 2 in DMUMPS_627"
     &                ,NODESTATE
           CALL MUMPS_ABORT()
      ENDIF
      IF (ISHIFT .LT.0_8) THEN
        WRITE(*,*) "Internal error 3 in DMUMPS_627",ISHIFT
        CALL MUMPS_ABORT()
      ENDIF
      IF (NELIM_ROOT) THEN
        IOLD=RCURRENT+int(LD,8)*int(NROW,8)+int(NELIM-1-NCB,8)
      ELSE
        IOLD = RCURRENT+int(LD,8)*int(NROW,8)-1_8
      ENDIF
      INEW = RCURRENT+int(LD,8)*int(NROW,8)+ISHIFT-1_8
      DO I = NROW, 1, -1
        IF (I.EQ.NROW .AND. ISHIFT.EQ.0_8.AND.
     &    .NOT. NELIM_ROOT) THEN
          IOLD=IOLD-int(LD,8)
          INEW=INEW-int(NCB,8)
          CYCLE
        ENDIF
        IF (NELIM_ROOT) THEN
          DO J=1,NELIM
            A( INEW ) = A( IOLD + int(- J + 1,8))
            INEW = INEW - 1_8
          ENDDO
        ELSE
          DO J=1, NCB
            A( INEW ) = A( IOLD + int(- J + 1, 8))
            INEW = INEW - 1_8
          ENDDO
        ENDIF
        IOLD = IOLD - int(LD,8)
      ENDDO
      IF (NELIM_ROOT) THEN
        NODESTATE=S_NOLCBCONTIG38
      ELSE
        NODESTATE=S_NOLCBCONTIG
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_627
      SUBROUTINE DMUMPS_700(BUFR,LBUFR,
     &     LBUFR_BYTES,
     &     root, N, IW, LIW, A, LA,
     &     NBPROCFILS, LRLU, IPTRLU, IWPOS, IWPOSCB,
     &     PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     &     COMP, LRLUS, IPOOL, LPOOL, LEAF,
     &     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &     KEEP,KEEP8, IFLAG, IERROR, COMM, COMM_LOAD,
     &     ITLOC, RHS_MUMPS,
     &     ND,PROCNODE_STEPS,SLAVEF )
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC ) :: root
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA, LRLU, IPTRLU, LRLUS
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER LBUFR, LBUFR_BYTES, N, LIW,
     &        IWPOS, IWPOSCB, COMP, COMM, COMM_LOAD, IFLAG,
     &        IERROR
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LEAF )
      INTEGER PTRIST(KEEP(28))
      INTEGER PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28)), ITLOC( N+KEEP(253) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER BUFR( LBUFR_BYTES ), NBPROCFILS( KEEP(28) )
      INTEGER IW( LIW )
      INTEGER ND(KEEP(28)), PROCNODE_STEPS(KEEP(28)),SLAVEF
      DOUBLE PRECISION A( LA )
      INTEGER   MYID
      INTEGER FILS( N ), PTRAIW(N), PTRARW( N )
      INTEGER INTARR(max(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
        INCLUDE 'mpif.h'
        INTEGER IERR
        INTEGER POSITION, LOCAL_M, LOCAL_N, LREQI
        INTEGER(8) :: LREQA, POS_ROOT
        INTEGER NSUBSET_ROW, NSUBSET_COL, IROOT, ISON, NSUBSET_COL_EFF
        INTEGER NSUPCOL_EFF
        INTEGER NBROWS_ALREADY_SENT, NBROWS_PACKET
        INTEGER NSUPROW, NSUPCOL, BBPCBP
        INCLUDE 'mumps_headers.h'
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   ISON, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NSUBSET_ROW, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NSUPROW, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NSUBSET_COL, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NSUPCOL, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                   COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NBROWS_PACKET, 1, MPI_INTEGER,
     &                   COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   BBPCBP, 1, MPI_INTEGER,
     &                   COMM, IERR )
        IF (BBPCBP .EQ. 1) THEN
          NSUBSET_COL_EFF = NSUBSET_COL - NSUPCOL
          NSUPCOL_EFF = 0
        ELSE
          NSUBSET_COL_EFF = NSUBSET_COL
          NSUPCOL_EFF = NSUPCOL
        ENDIF
        IROOT = KEEP( 38 )
        IF ( PTRIST( STEP(IROOT) ) .NE. 0 .OR.
     &       PTLUST_S( STEP(IROOT)) .NE. 0 ) THEN
          IF (NBROWS_ALREADY_SENT + NBROWS_PACKET .EQ. NSUBSET_ROW
     &       - NSUPROW .OR.  NSUBSET_ROW - NSUPROW.EQ.0 .OR.
     &       NSUBSET_COL_EFF .EQ. 0)THEN
            NBPROCFILS(STEP(IROOT)) = NBPROCFILS(STEP(IROOT))-1
            IF ( NBPROCFILS( STEP(IROOT) ) .eq. 0 ) THEN
              IF (KEEP(201).EQ.1) THEN 
                 CALL DMUMPS_681(IERR)
              ELSEIF (KEEP(201).EQ.2) THEN 
                 CALL DMUMPS_580(IERR)              
              ENDIF
              CALL DMUMPS_507( N, IPOOL, LPOOL,
     &             PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     &             KEEP(80), KEEP(47),
     &             STEP, IROOT + N)
              IF (KEEP(47) .GE. 3) THEN
                 CALL DMUMPS_500(
     &                IPOOL, LPOOL, 
     &                PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &                MYID, STEP, N, ND, FILS )
              ENDIF
            ENDIF
          ENDIF
        ELSE
           IF (NBROWS_ALREADY_SENT + NBROWS_PACKET .EQ.
     &       NSUBSET_ROW - NSUPROW .OR.
     &        NSUBSET_ROW - NSUPROW.EQ.0 .OR.
     &        NSUBSET_COL_EFF .EQ. 0)THEN
             NBPROCFILS(STEP( IROOT ) ) = -1
           ENDIF
           IF (KEEP(60) == 0) THEN
            CALL DMUMPS_284( root, IROOT, N,
     &                     IW, LIW, A, LA,
     &                     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     &                     LRLU, IPTRLU,
     &                     IWPOS, IWPOSCB, PTRIST, PTRAST,
     &                     STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &                     COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
            IF ( IFLAG .LT. 0 ) RETURN
           ELSE
             PTRIST(STEP(IROOT)) = -55555
           ENDIF
        END IF
      IF (KEEP(60) .EQ.0) THEN
        IF ( PTRIST(STEP(IROOT)) .GE. 0 ) THEN
          IF ( PTRIST(STEP(IROOT)) .NE. 0 ) THEN
               LOCAL_N  = -IW( PTRIST(STEP( IROOT )) + KEEP(IXSZ)    )
               LOCAL_M  =  IW( PTRIST(STEP( IROOT )) + 1 + KEEP(IXSZ))
               POS_ROOT = PAMASTER(STEP( IROOT ))
          ELSE
               LOCAL_N = IW( PTLUST_S(STEP( IROOT ) ) + 1 + KEEP(IXSZ))
               LOCAL_M = IW( PTLUST_S(STEP( IROOT ) ) + 2 + KEEP(IXSZ))
               POS_ROOT = PTRFAC(IW(PTLUST_S(STEP(IROOT))+4+
     &                    KEEP(IXSZ)))
          END IF
         ENDIF
      ELSE
          LOCAL_M = root%SCHUR_LLD
          LOCAL_N = root%SCHUR_NLOC
      ENDIF
        IF ( (BBPCBP.EQ.1).AND. (NBROWS_ALREADY_SENT.EQ.0).AND.
     &     (min(NSUPROW, NSUPCOL) .GT. 0)
     &     ) THEN
         LREQI = NSUPROW+NSUPCOL
         LREQA = int(NSUPROW,8) * int(NSUPCOL,8)
         IF ( (LREQA.NE.0_8) .AND.
     &       (PTRIST(STEP(IROOT)).LT.0).AND.
     &       KEEP(60)==0) THEN
          WRITE(*,*) ' Error in DMUMPS_700'
          CALL MUMPS_ABORT()
         ENDIF
         CALL DMUMPS_22(.FALSE.,0_8,.FALSE.,.FALSE.,
     &     MYID,N,KEEP,KEEP8,IW,LIW,A, LA,
     &     LRLU, IPTRLU, IWPOS, IWPOSCB, PTRIST,
     &     PTRAST, STEP, PIMASTER, PAMASTER,
     &     LREQI, LREQA, -1234, S_NOTFREE, .FALSE.,
     &     COMP, LRLUS, IFLAG, IERROR
     &          )
         IF ( IFLAG .LT. 0 ) RETURN
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   IW( IWPOSCB + 1 ), LREQI,
     &                   MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   A( IPTRLU + 1_8 ), int(LREQA),
     &                   MPI_DOUBLE_PRECISION, COMM, IERR )
         CALL DMUMPS_38( NSUPROW, NSUPCOL,
     &                     IW( IWPOSCB + 1 ), 
     &                     IW( IWPOSCB + NSUPROW + 1 ), NSUPCOL,
     &                     A( IPTRLU + 1_8 ),
     &                     A( 1 ), 
     &                     LOCAL_M, LOCAL_N,
     &                  root%RHS_ROOT(1,1), root%RHS_NLOC,
     &                  1)
         IWPOSCB = IWPOSCB + LREQI
         IPTRLU  = IPTRLU  + LREQA
         LRLU    = LRLU    + LREQA
         LRLUS   = LRLUS   + LREQA
         CALL DMUMPS_471(.FALSE.,.FALSE.,
     &                    LA-LRLUS,0_8,-LREQA,KEEP,KEEP8,LRLU)
        ENDIF  
        LREQI = NBROWS_PACKET + NSUBSET_COL_EFF
        LREQA = int(NBROWS_PACKET,8) * int(NSUBSET_COL_EFF,8)
        IF ( (LREQA.NE.0_8) .AND.
     &       (PTRIST(STEP(IROOT)).LT.0).AND.
     &       KEEP(60)==0) THEN
         WRITE(*,*) ' Error in DMUMPS_700'
         CALL MUMPS_ABORT()
        ENDIF
        IF (LREQA.NE.0_8) THEN
          CALL DMUMPS_22(.FALSE.,0_8,.FALSE.,.FALSE.,
     &     MYID,N,KEEP,KEEP8,IW,LIW,A, LA,
     &     LRLU, IPTRLU, IWPOS, IWPOSCB, PTRIST,
     &     PTRAST, STEP, PIMASTER, PAMASTER,
     &     LREQI, LREQA, -1234, S_NOTFREE, .FALSE.,
     &     COMP, LRLUS, IFLAG, IERROR
     &          )
          IF ( IFLAG .LT. 0 ) RETURN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   IW( IWPOSCB + 1 ), LREQI,
     &                   MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   A( IPTRLU + 1_8 ), int(LREQA),
     &                   MPI_DOUBLE_PRECISION, COMM, IERR )
          IF (KEEP(60).EQ.0) THEN
            CALL DMUMPS_38( NBROWS_PACKET, NSUBSET_COL_EFF,
     &                     IW( IWPOSCB + 1 ),
     &                     IW( IWPOSCB + NBROWS_PACKET + 1 ),
     &                     NSUPCOL_EFF,
     &                     A( IPTRLU + 1_8 ),
     &                     A( POS_ROOT ), LOCAL_M, LOCAL_N,
     &                  root%RHS_ROOT(1,1), root%RHS_NLOC,
     &                  0)   
          ELSE
            CALL DMUMPS_38( NBROWS_PACKET, NSUBSET_COL_EFF,
     &                     IW( IWPOSCB + 1 ),
     &                     IW( IWPOSCB + NBROWS_PACKET + 1 ),
     &                     NSUPCOL_EFF,
     &                     A( IPTRLU + 1_8 ),
     &                     root%SCHUR_POINTER(1),
     &                     root%SCHUR_LLD , root%SCHUR_NLOC,
     &                  root%RHS_ROOT(1,1), root%RHS_NLOC,
     &                  0)  
          ENDIF
          IWPOSCB = IWPOSCB + LREQI
          IPTRLU  = IPTRLU  + LREQA
          LRLU    = LRLU    + LREQA
          LRLUS   = LRLUS   + LREQA
          CALL DMUMPS_471(.FALSE.,.FALSE.,
     &                    LA-LRLUS,0_8,-LREQA,KEEP,KEEP8,LRLU)
        ENDIF
      RETURN
      END SUBROUTINE DMUMPS_700
      SUBROUTINE DMUMPS_762(PIV, DETER, NEXP)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: PIV
      DOUBLE PRECISION, intent(inout) :: DETER
      INTEGER, intent(inout) :: NEXP
      DETER=DETER*fraction(PIV)
      NEXP=NEXP+exponent(PIV)+exponent(DETER)
      DETER=fraction(DETER)
      RETURN
      END SUBROUTINE DMUMPS_762
      SUBROUTINE DMUMPS_761(PIV, DETER, NEXP)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: PIV
      DOUBLE PRECISION, intent(inout) :: DETER
      INTEGER, intent(inout) :: NEXP
      DETER=DETER*fraction(PIV)
      NEXP=NEXP+exponent(PIV)+exponent(DETER)
      DETER=fraction(DETER)
      RETURN
      END SUBROUTINE DMUMPS_761
      SUBROUTINE DMUMPS_763(BLOCK_SIZE,IPIV,
     &                      MYROW, MYCOL, NPROW, NPCOL,
     &                      A, LOCAL_M, LOCAL_N, N, MYID,
     &                      DETER,NEXP,SYM)
      IMPLICIT NONE
      INTEGER, intent (in)    :: SYM
      INTEGER, intent (inout) :: NEXP
      DOUBLE PRECISION, intent (inout) :: DETER
      INTEGER, intent (in)    :: BLOCK_SIZE, NPROW, NPCOL,
     &                           LOCAL_M, LOCAL_N, N
      INTEGER, intent (in)    :: MYROW, MYCOL, MYID, IPIV(LOCAL_M)
      DOUBLE PRECISION, intent(in) :: A(*)
      INTEGER  I,IMX,DI,NBLOCK,IBLOCK,ILOC,JLOC,
     &         ROW_PROC,COL_PROC, K
      DI = LOCAL_M + 1 
      NBLOCK = ( N - 1 ) / BLOCK_SIZE
      DO IBLOCK = 0, NBLOCK
        ROW_PROC = mod( IBLOCK, NPROW ) 
        IF ( MYROW.EQ.ROW_PROC ) THEN
          COL_PROC = mod( IBLOCK, NPCOL )
          IF ( MYCOL.EQ.COL_PROC ) THEN
            ILOC = ( IBLOCK / NPROW ) * BLOCK_SIZE
            JLOC = ( IBLOCK / NPCOL ) * BLOCK_SIZE
            I   =   ILOC + JLOC *  LOCAL_M + 1
            IMX = min(ILOC+BLOCK_SIZE,LOCAL_M)
     &            + (min(JLOC+BLOCK_SIZE,LOCAL_N)-1)*LOCAL_M
     &            + 1
            K=1
            DO WHILE ( I .LT. IMX )
              CALL DMUMPS_762(A(I),DETER,NEXP)
              IF (SYM.NE.1) THEN 
                IF (IPIV(ILOC+K) .NE. IBLOCK*BLOCK_SIZE+K) THEN
                  DETER = -DETER
                ENDIF
              ENDIF
              K = K + 1
              I = I + DI
            END DO
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_763
      SUBROUTINE DMUMPS_764(
     &           COMM, DETER_IN, NEXP_IN,
     &           DETER_OUT, NEXP_OUT, NPROCS)
      IMPLICIT NONE
      INTEGER, intent(in) :: COMM, NPROCS
      DOUBLE PRECISION, intent(in) :: DETER_IN
      INTEGER,intent(in) :: NEXP_IN
      DOUBLE PRECISION,intent(out):: DETER_OUT
      INTEGER,intent(out):: NEXP_OUT
      INTEGER            :: IERR_MPI
      EXTERNAL DMUMPS_771
      INTEGER TWO_SCALARS_TYPE, DETERREDUCE_OP
      DOUBLE PRECISION :: INV(2)
      DOUBLE PRECISION :: OUTV(2)
      INCLUDE 'mpif.h'
      IF (NPROCS .EQ. 1) THEN
        DETER_OUT = DETER_IN
        NEXP_OUT  = NEXP_IN
        RETURN
      ENDIF
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_DOUBLE_PRECISION,
     &                         TWO_SCALARS_TYPE,
     &                         IERR_MPI)
      CALL MPI_TYPE_COMMIT(TWO_SCALARS_TYPE, IERR_MPI)
      CALL MPI_OP_CREATE(DMUMPS_771,  
     &                   .TRUE.,             
     &                   DETERREDUCE_OP,     
     &                   IERR_MPI)
      INV(1)=DETER_IN
      INV(2)=dble(NEXP_IN)
      CALL MPI_ALLREDUCE( INV, OUTV, 1, TWO_SCALARS_TYPE,
     &                    DETERREDUCE_OP, COMM, IERR_MPI)
      CALL MPI_OP_FREE(DETERREDUCE_OP, IERR_MPI)
      CALL MPI_TYPE_FREE(TWO_SCALARS_TYPE, IERR_MPI)
      DETER_OUT = OUTV(1)
      NEXP_OUT  = int(OUTV(2))
      RETURN
      END SUBROUTINE DMUMPS_764
      SUBROUTINE DMUMPS_771(INV, INOUTV, NEL, DATATYPE)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NEL, DATATYPE
      DOUBLE PRECISION, INTENT(IN)    :: INV    ( 2 * NEL )
      DOUBLE PRECISION, INTENT(INOUT) :: INOUTV ( 2 * NEL )
      INTEGER I, TMPEXPIN, TMPEXPINOUT
      DO I = 1, NEL
        TMPEXPIN    = int(INV   (I*2))
        TMPEXPINOUT = int(INOUTV(I*2))
        CALL DMUMPS_762(INV(I*2-1), 
     &                          INOUTV(I*2-1), 
     &                          TMPEXPINOUT)   
        TMPEXPINOUT = TMPEXPINOUT + TMPEXPIN
        INOUTV(I*2) = dble(TMPEXPINOUT)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_771
      SUBROUTINE DMUMPS_765(DETER, NEXP)
      IMPLICIT NONE
      INTEGER, intent (inout) :: NEXP
      DOUBLE PRECISION, intent (inout) :: DETER
      DETER=DETER*DETER
      NEXP=NEXP+NEXP
      RETURN
      END SUBROUTINE DMUMPS_765
      SUBROUTINE DMUMPS_766(DETER, NEXP)
      IMPLICIT NONE
      INTEGER, intent (inout) :: NEXP
      DOUBLE PRECISION, intent (inout) :: DETER
      DETER=1.0D0/DETER
      NEXP=-NEXP
      RETURN
      END SUBROUTINE DMUMPS_766
      SUBROUTINE DMUMPS_767(DETER, N, VISITED, PERM)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(inout) :: DETER
      INTEGER, intent(in)    :: N
      INTEGER, intent(inout) :: VISITED(N)
      INTEGER, intent(in)    :: PERM(N)
      INTEGER I, J, K
      K = 0
      DO I = 1, N
        IF (VISITED(I) .GT. N) THEN
          VISITED(I)=VISITED(I)-N-N-1
          CYCLE
        ENDIF
        J = PERM(I)
        DO WHILE (J.NE.I)
          VISITED(J) = VISITED(J) + N + N + 1
          K = K + 1
          J = PERM(J)
        ENDDO
      ENDDO
      IF (mod(K,2).EQ.1) THEN
        DETER = -DETER
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_767
      SUBROUTINE DMUMPS_224(NFRONT,NASS,IBEGKJI, LPIV, TIPIV,
     &    N,INODE,IW,LIW,A,LA,
     &    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU,SEUIL,KEEP,KEEP8,
     &     DKEEP,PIVNUL_LIST,LPN_LIST,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U)
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER IBEGKJI, LPIV 
      INTEGER TIPIV(LPIV)
      INTEGER(8) :: LA
      DOUBLE PRECISION A(LA)
      INTEGER NFRONT,NASS,N,LIW,INODE,IFLAG,INOPV,NOFFW
      DOUBLE PRECISION UU, SEUIL
      INTEGER IW(LIW) 
      INTEGER IOLDPS
      INTEGER(8) :: POSELT
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      DOUBLE PRECISION SWOP
      INTEGER(8) :: APOS, IDIAG
      INTEGER(8) :: J1, J2, JJ, J3_8
      INTEGER(8) :: NFRONT8
      INTEGER ILOC
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      DOUBLE PRECISION RZERO, RMAX, AMROW, ONE
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J3,ISW,ISWPS1
      INTEGER ISWPS2,KSW, HF
      INCLUDE 'mumps_headers.h'
      INTEGER DMUMPS_IXAMAX
      INTRINSIC max
      DATA RZERO /0.0D0/
      DATA ONE /1.0D0/
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
      INTEGER XSIZE
        PIVNUL  = DKEEP(1)
        FIXA    = DKEEP(2)
        CSEUIL  = SEUIL
        NFRONT8=int(NFRONT,8)
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
        NPIVP1  = NPIV + 1
        IF (KEEP(201).EQ.1) THEN
          CALL DMUMPS_667(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
          CALL DMUMPS_667(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
        ILOC    = NPIVP1 - IBEGKJI + 1
        TIPIV(ILOC) = ILOC
        NASSW   = iabs(IW(IOLDPS+3+XSIZE))
        IF(INOPV .EQ. -1) THEN
           APOS = POSELT + NFRONT8*int(NPIVP1-1,8) + int(NPIV,8)
           IDIAG = APOS
           IF(abs(A(APOS)).LT.SEUIL) THEN
              IF (dble(A(APOS)) .GE. RZERO) THEN
                A(APOS) = CSEUIL
              ELSE
                A(APOS) = -CSEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ELSE IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_762(A(APOS), DKEEP(6), KEEP(259))
           ENDIF
           IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
             IF (KEEP(251).EQ.0) THEN
               CALL DMUMPS_680( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
             ENDIF
             CALL DMUMPS_680( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
          DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + NFRONT8*int(IPIV-1,8) + int(NPIV,8)
            JMAX = 1
            IF (UU.GT.RZERO) GO TO 340
            IF (A(APOS).EQ.ZERO) GO TO 630
            GO TO 380
  340       AMROW = RZERO
            J1 = APOS
            J2 = APOS +int(- NPIV + NASS - 1,8)
             J3    = NASS -NPIV
             JMAX  = DMUMPS_IXAMAX(J3,A(J1),1)
             JJ    = int(JMAX,8) + J1 - 1_8
             AMROW = abs(A(JJ))
            RMAX = AMROW
            J1 = J2 + 1_8
            J2 = APOS +int(- NPIV + NFRONT - KEEP(253)- 1,8)
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = max(abs(A(JJ)),RMAX)
  360       CONTINUE
  370       IDIAG = APOS + int(IPIV - NPIVP1,8)
            IF (RMAX.LE.PIVNUL) THEN
               KEEP(109) = KEEP(109)+1
               ISW = IOLDPS+IW(IOLDPS+1+KEEP(IXSZ))+6+KEEP(IXSZ)+
     &                      IW(IOLDPS+5+KEEP(IXSZ))+IPIV-NPIVP1
               PIVNUL_LIST(KEEP(109)) = IW(ISW)
               IF(dble(FIXA).GT.RZERO) THEN
                  IF(dble(A(IDIAG)) .GE. RZERO) THEN
                     A(IDIAG) = FIXA
                  ELSE
                     A(IDIAG) = -FIXA
                  ENDIF
               ELSE
                 J1 = APOS
                 J2 = APOS +int(- NPIV + NFRONT - KEEP(253) - 1,8)
                 DO JJ=J1,J2
                   A(JJ)= ZERO
                 ENDDO
                 A(IDIAG) = -FIXA
               ENDIF
               JMAX = IPIV - NPIV
               GOTO 385   
            ENDIF
            IF (abs(A(IDIAG)).GT.max(UU*RMAX,SEUIL)) THEN
               JMAX = IPIV - NPIV
               GO TO 380
            ENDIF
            IF (AMROW.LE.max(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
            IF (KEEP(258).NE.0) THEN
              CALL DMUMPS_762( A(APOS+int(JMAX-1,8)),
     &                                 DKEEP(6),
     &                                 KEEP(259))
            ENDIF
  385       CONTINUE
            IF (IPIV.EQ.NPIVP1) GO TO 400
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8)*NFRONT8
            J2 = J1 + NFRONT8 - 1_8
            J3_8 = POSELT + int(IPIV-1,8)*NFRONT8
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3_8)
              A(J3_8) = SWOP
              J3_8 = J3_8 + 1_8
  390       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NPIVP1
            ISWPS2 = IOLDPS + HF - 1 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            KEEP(260)=-KEEP(260)
            TIPIV(ILOC) = ILOC + JMAX - 1
            J1 = POSELT + int(NPIV,8)
            J2 = POSELT + int(NPIV + JMAX - 1,8)
            DO 410 KSW=1,NASS
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT8
              J2 = J2 + NFRONT8
  410       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + HF - 1 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420  
  460     CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 430
  630 CONTINUE
      IFLAG = -10
      WRITE(*,*) 'NIV2:Detected 0 pivot, INODE,NPIV=',INODE,NPIV
      GOTO 430
  420 CONTINUE
              IF (KEEP(201).EQ.1) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_680( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, IPIV, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_680( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIV+JMAX, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
  430 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_224
      SUBROUTINE  DMUMPS_294( COMM_LOAD, ASS_IRECV, 
     &             N, INODE, FPERE,
     &             IW, LIW, 
     &             IOLDPS, POSELT, A, LA, LDA_FS, 
     &             IBEGKJI, IEND, TIPIV, LPIV, LASTBL, NB_BLOC_FAC,
     &
     &             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL, 
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC,
     &             STEP, PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE DMUMPS_COMM_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, IBEGKJI, IEND, LPIV, 
     &        IOLDPS, LDA_FS, NB_BLOC_FAC
      INTEGER(8) :: POSELT, LA
      INTEGER IW(LIW), TIPIV(LPIV)
      LOGICAL LASTBL
      DOUBLE PRECISION A(LA)
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT ) 
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NBFIN, IFLAG, IERROR, LEAF, LPOOL,
     &        SLAVEF, ICNTL(40)
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB, COMP
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     &        ITLOC(N+KEEP(253)), FILS(N),
     &        PTRARW(LPTRAR), PTRAIW(LPTRAR), 
     &        ND( KEEP(28) ), FRERE( KEEP(28) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER(8) :: PTRAST  (KEEP(28))
      INTEGER(8) :: PTRFAC  (KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)),
     &        NBPROCFILS(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(max(1,KEEP(13)))
      EXTERNAL  DMUMPS_329
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: APOS, LREQA
      INTEGER NPIV, NCOL, PDEST, NSLAVES
      INTEGER IERR, LREQI
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION FLOP1,FLOP2
      NSLAVES= IW(IOLDPS+5+KEEP(IXSZ))
          IF (NSLAVES.EQ.0) THEN
           WRITE(6,*) ' ERROR 1 in DMUMPS_294 '
           CALL MUMPS_ABORT()
          ENDIF
      NPIV   = IEND - IBEGKJI + 1
      NCOL   = LDA_FS - IBEGKJI + 1
      APOS   = POSELT + int(LDA_FS,8)*int(IBEGKJI-1,8) +
     &                  int(IBEGKJI - 1,8)
      IF (IBEGKJI > 0) THEN
       CALL MUMPS_511( LDA_FS, IBEGKJI-1, LPIV,
     &                            KEEP(50),2,FLOP1)
      ELSE
        FLOP1=0.0D0
      ENDIF
      CALL MUMPS_511( LDA_FS, IEND, LPIV,
     &                           KEEP(50),2,FLOP2)
      FLOP2 = FLOP1 - FLOP2
      CALL DMUMPS_190(1, .FALSE., FLOP2, KEEP,KEEP8)
      IF ((NPIV.GT.0) .OR. 
     &    ((NPIV.EQ.0).AND.(LASTBL)) ) THEN
        PDEST  = IOLDPS + 6 + KEEP(IXSZ)
        IERR = -1
        IF ( NPIV .NE. 0 ) THEN
          NB_BLOC_FAC = NB_BLOC_FAC + 1
        END IF
        DO WHILE (IERR .EQ.-1)
          CALL DMUMPS_65( INODE, LDA_FS, NCOL, 
     &               NPIV, FPERE, LASTBL, TIPIV, A(APOS),
     &               IW(PDEST), NSLAVES, KEEP(50), NB_BLOC_FAC,
     &               COMM, IERR )
        IF (IERR.EQ.-1) THEN
           BLOCKING  = .FALSE.
           SET_IRECV = .TRUE.
           MESSAGE_RECEIVED = .FALSE.
           CALL DMUMPS_329( COMM_LOAD, ASS_IRECV, 
     &      BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &      MPI_ANY_SOURCE, MPI_ANY_TAG,
     &      STATUS, BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     &      IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. )
           IF (MESSAGE_RECEIVED) POSELT = PTRAST(STEP(INODE))
           IF ( IFLAG .LT. 0 ) GOTO 500
         ENDIF
        ENDDO
        IF (IERR .EQ. -2 .OR. IERR.EQ.-3 ) THEN
          IF (IERR.EQ.-2) IFLAG = -17
          IF (IERR.EQ.-3) IFLAG = -20
          LREQA = int(NCOL,8)*int(NPIV,8)
          LREQI = NPIV + 6 + 2*NSLAVES
          CALL MUMPS_731(
     &    int(LREQI,8) * int(KEEP(34),8) + LREQA * int(KEEP(35),8),
     &    IERROR)
          GOTO 300
        ENDIF
      ENDIF
      GOTO 500
  300 CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
 500  RETURN
      END SUBROUTINE  DMUMPS_294
      SUBROUTINE DMUMPS_273( ROOT, 
     &    INODE, NELIM, NSLAVES, ROW_LIST,
     &    COL_LIST, SLAVE_LIST, 
     &
     &    PROCNODE_STEPS, IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S,
     &    ITLOC, RHS_MUMPS, COMP,
     &    IFLAG, IERROR, 
     &    IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,KEEP8,
     &    COMM,COMM_LOAD,FILS,ND )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: ROOT
      INTEGER INODE, NELIM, NSLAVES 
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER ROW_LIST(*), COL_LIST(*), 
     &        SLAVE_LIST(*)
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER PTRIST( KEEP(28) ), PTLUST_S(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), ITLOC( N + KEEP(253) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER IFLAG, IERROR
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF
      INTEGER COMM,COMM_LOAD,ND(KEEP(28)),FILS(N)
      INTEGER IROOT, TYPE_INODE, DEB_ROW, DEB_COL,
     &        NOINT
      INTEGER(8) :: NOREAL
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MUMPS_330
      EXTERNAL MUMPS_330
      IROOT        = KEEP(38)
      NSTK_S(STEP(IROOT))= NSTK_S(STEP(IROOT)) - 1
      KEEP(42) = KEEP(42) + NELIM
      TYPE_INODE= MUMPS_330(PROCNODE_STEPS(STEP(INODE)), SLAVEF)
      IF (TYPE_INODE.EQ.1) THEN 
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + 1
        ELSE 
         KEEP(41) = KEEP(41) + 3
        ENDIF
      ELSE
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + NSLAVES
        ELSE 
         KEEP(41) = KEEP(41) + 2*NSLAVES + 1
        ENDIF
      ENDIF
      IF  (NELIM.EQ.0) THEN
        PIMASTER(STEP(INODE)) = 0 
      ELSE
       NOINT = 6 + NSLAVES + NELIM  + NELIM + KEEP(IXSZ)
       NOREAL= 0_8
       CALL DMUMPS_22(.FALSE.,0_8,.FALSE.,.FALSE.,
     &   MYID,N,KEEP,KEEP8,IW,LIW, A, LA,
     &   LRLU, IPTRLU,IWPOS,IWPOSCB,
     &   PTRIST,PTRAST,STEP, PIMASTER, PAMASTER,
     &   NOINT, NOREAL, INODE, S_NOTFREE, .TRUE.,
     &   COMP, LRLUS, IFLAG, IERROR
     &      )
       IF ( IFLAG .LT. 0 ) THEN
         WRITE(*,*) ' Failure in int space allocation in CB area ',
     &    ' during assembly of root : DMUMPS_273',
     &    ' size required was :', NOINT,
     &    'INODE=',INODE,' NELIM=',NELIM, ' NSLAVES=', NSLAVES
         RETURN
        ENDIF
        PIMASTER(STEP( INODE )) = IWPOSCB + 1
        PAMASTER(STEP( INODE )) = IPTRLU  + 1_8
        IW( IWPOSCB + 1+KEEP(IXSZ) ) = 2*NELIM
        IW( IWPOSCB + 2+KEEP(IXSZ) ) = NELIM
        IW( IWPOSCB + 3+KEEP(IXSZ) ) = 0
        IW( IWPOSCB + 4+KEEP(IXSZ) ) = 0
        IW( IWPOSCB + 5+KEEP(IXSZ) ) = 1
        IW( IWPOSCB + 6+KEEP(IXSZ) ) = NSLAVES
        IF (NSLAVES.GT.0) THEN
         IW( IWPOSCB+7+KEEP(IXSZ):IWPOSCB+7+KEEP(IXSZ)+NSLAVES-1) = 
     &                   SLAVE_LIST(1:NSLAVES)
        ENDIF
        DEB_ROW = IWPOSCB+7+NSLAVES+KEEP(IXSZ)
        IW(DEB_ROW : DEB_ROW+NELIM -1) = ROW_LIST(1:NELIM)
        DEB_COL = DEB_ROW + NELIM
        IW(DEB_COL : DEB_COL+NELIM -1) = COL_LIST(1:NELIM)
      ENDIF
      IF (NSTK_S(STEP(IROOT)) .EQ. 0 ) THEN
          CALL DMUMPS_507(N, IPOOL, LPOOL, PROCNODE_STEPS,
     &         SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     &         STEP, IROOT )
          IF (KEEP(47) .GE. 3) THEN
             CALL DMUMPS_500(
     &            IPOOL, LPOOL, 
     &            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &            MYID, STEP, N, ND, FILS )
          ENDIF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_273
      SUBROUTINE DMUMPS_363(N,FRERE, STEP, FILS,
     &     NA,LNA,NE,ND, DAD, LDAD, USE_DAD,
     &     NSTEPS,PERM,SYM,INFO,LP,K215,K234,K55,
     &     PROCNODE,SLAVEF, PEAK,SBTR_WHICH_M
     &     )
      IMPLICIT NONE
      INTEGER N,PERM,SYM, NSTEPS, LNA, LP,LDAD
      INTEGER FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER NA(LNA), NE(NSTEPS), ND(NSTEPS)
      INTEGER K215,K234,K55
      INTEGER DAD(LDAD)
      LOGICAL USE_DAD
      INTEGER INFO(40)
      INTEGER SLAVEF,PROCNODE(NSTEPS)
      INTEGER :: SBTR_WHICH_M
      EXTERNAL MUMPS_275
      INTEGER MUMPS_275
      DOUBLE PRECISION PEAK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COST_TRAV
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH
      INTEGER IFATH,IN,NSTK,INODE,I,allocok,LOCAL_PERM
      INTEGER(8) NCB
      INTEGER(8) NELIM,NFR
      INTEGER NFR4,NELIM4
      INTEGER LEAF,NBLEAF,NBROOT, SIZE_TAB
      INTEGER, DIMENSION (:), ALLOCATABLE :: IPOOL,TNSTK
      INTEGER, DIMENSION (:), ALLOCATABLE,TARGET :: SON,TEMP
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: M,M_TOTAL, fact
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: TAB1,TAB2
      INTEGER, DIMENSION (:), POINTER :: TAB
      INTEGER dernier,fin
      INTEGER cour,II
      INTEGER CB_current,CB_MAX,ROOT_OF_CUR_SBTR
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: T1,T2
      INTEGER, DIMENSION (:), ALLOCATABLE :: RESULT
      INTEGER(8) MEM_SIZE,FACT_SIZE,SUM,MEM_SEC_PERM,FACT_SIZE_T,
     &     MEM_SIZE_T,TOTAL_MEM_SIZE,TMP_TOTAL_MEM_SIZE,TMP_SUM,
     &     SIZECB, SIZECB_LASTSON
      INTEGER(8) TMP8
      LOGICAL   SBTR_M
      INTEGER FIRST_LEAF,SIZE_SBTR
      EXTERNAL MUMPS_170,MUMPS_167
      LOGICAL MUMPS_170,MUMPS_167
      DOUBLE PRECISION COST_NODE
      INCLUDE 'mumps_headers.h'
      TOTAL_MEM_SIZE=0_8
      ROOT_OF_CUR_SBTR=0
      IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.
     &     (PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4).OR.
     &     (PERM.EQ.5).OR.(PERM.EQ.6))THEN
         LOCAL_PERM=0
      ENDIF
      SBTR_M=.FALSE.
      MEM_SIZE=0_8
      FACT_SIZE=0_8
      IF ((PERM.LT.0 .OR. PERM.GT.7)) THEN
         WRITE(*,*) "Internal Error in DMUMPS_363",PERM
         CALL MUMPS_ABORT()
      END IF
      NBLEAF = NA(1)
      NBROOT = NA(2)
      IF((PERM.EQ.0).AND.(NBROOT.EQ.NBLEAF)) RETURN
      IF ((PERM.NE.7).AND.(SBTR_M.OR.(PERM.EQ.2))) THEN
         IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
            ALLOCATE(M_TOTAL(NSTEPS), stat=allocok )
            IF (allocok > 0) THEN
               IF ( LP .GT. 0 )
     &              WRITE(LP,*)'Memory allocation error in
     &              DMUMPS_363'
               INFO(1)=-7
               INFO(2)=NSTEPS
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IF(PERM.NE.7)THEN
         ALLOCATE(M(NSTEPS),stat=allocok )
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error
     &in DMUMPS_363'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
      ENDIF
      ALLOCATE( IPOOL(NBLEAF), fact(NSTEPS),TNSTK(NSTEPS),
     &     stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_363'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      II=0
      DO I=1,NSTEPS
         TNSTK(I) = NE(I)
         IF(NE(I).GE.II) II=NE(I)
      ENDDO
      SIZE_TAB=max(II,NBROOT)
      ALLOCATE(SON(II), TEMP(II),
     &         TAB1(SIZE_TAB), TAB2(SIZE_TAB), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_363'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      ALLOCATE(T1(SIZE_TAB),T2(SIZE_TAB),
     &         RESULT(SIZE_TAB),stat=allocok)
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_363'
        INFO(1)=-7
        INFO(2)=SIZE_TAB
        RETURN
      ENDIF
      IF(PERM.EQ.7) THEN
         GOTO 001
      ENDIF
      IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
         ALLOCATE(COST_TRAV(NSTEPS), stat=allocok )
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error
     &       in DMUMPS_363'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
         COST_TRAV=0.0D0
         COST_NODE=0.0d0
      ENDIF
      IF(NBROOT.EQ.NBLEAF)THEN
        IF((PERM.NE.1).OR.(PERM.EQ.4).OR.(PERM.EQ.6))THEN
          WRITE(*,*)'Internal Error in reordertree:'
          WRITE(*,*)'  problem with perm parameter in reordertree'
          CALL MUMPS_ABORT()
        ENDIF
        DO I=1,NBROOT
          TAB1(I)=int(ND(STEP(NA(I+2+NBLEAF))),8)
          IPOOL(I)=NA(I+2+NBLEAF)
          M(STEP(IPOOL(I)))=TAB1(I)*TAB1(I)
        ENDDO
        CALL DMUMPS_462(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &    RESULT,T1,T2)
        GOTO 789
      ENDIF
      IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
         ALLOCATE(DEPTH(NSTEPS),stat=allocok)
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error in
     &           DMUMPS_363'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
         DEPTH=0
         NBROOT = NA(2)
         IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
         fin=NBROOT
         LEAF=NA(1)
 499     CONTINUE
         INODE=IPOOL(fin)
         IF(INODE.LT.0)THEN
            WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
            CALL MUMPS_ABORT()
         ENDIF
         IN=INODE
 4602    IN = FILS(IN)
         IF (IN .GT. 0 ) THEN
            GOTO 4602
         ENDIF
         IN=-IN
         DO I=1,NE(STEP(INODE))
            SON(I)=IN
            IN=FRERE(STEP(IN))
         ENDDO
         DO I=1,NE(STEP(INODE))
            IPOOL(fin)=SON(I)
            DEPTH(STEP(SON(I)))=DEPTH(STEP(INODE))+1
            SON(I)=0
            fin=fin+1
         ENDDO
         IF(NE(STEP(INODE)).EQ.0)THEN
            LEAF=LEAF-1
         ELSE
            fin=fin-1
            GOTO 499
         ENDIF
         fin=fin-1
         IF(fin.EQ.0) GOTO 489
         GOTO 499
 489     CONTINUE
      ENDIF
      DO I=1,NSTEPS
         M(I)=0_8
         IF (SBTR_M.OR.(PERM.EQ.2))  THEN
            IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
               M_TOTAL(I)=0_8
            ENDIF
         ENDIF
      ENDDO
      DO I=1,NSTEPS
         fact(I)=0_8
      ENDDO
      IPOOL(1:NBLEAF)=NA(3:2+NBLEAF)
      LEAF = NBLEAF + 1
 91   CONTINUE
        IF (LEAF.NE.1) THEN
           LEAF = LEAF -1
           INODE = IPOOL(LEAF)
        ENDIF
 96     CONTINUE
        NFR    = int(ND(STEP(INODE)),8)
        NSTK   = NE(STEP(INODE))
        NELIM4 = 0
        IN = INODE
 101    NELIM4 = NELIM4 + 1
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 101
        NELIM=int(NELIM4,8)
        IF(NE(STEP(INODE)).EQ.0) THEN
           M(STEP(INODE))=NFR*NFR
           IF (SBTR_M.OR.(PERM.EQ.2))  THEN
                 M_TOTAL(STEP(INODE))=NFR*NFR
           ENDIF
        ENDIF
        IF((PERM.EQ.4).OR.(PERM.EQ.3))THEN
           IF(MUMPS_170(PROCNODE(STEP(INODE)),
     &         SLAVEF))THEN
              DEPTH(STEP(INODE))=0
           ENDIF
        ENDIF
        IF ( SYM .eq. 0 ) THEN
          fact(STEP(INODE))=fact(STEP(INODE))+
     &      (2_8*NFR*NELIM)-(NELIM*NELIM)
        ELSE
          fact(STEP(INODE))=fact(STEP(INODE))+NFR*NELIM
        ENDIF
        IF (USE_DAD) THEN
          IFATH = DAD( STEP(INODE) )
        ELSE
          IN = INODE
 113      IN = FRERE(IN)
          IF (IN.GT.0) GO TO 113
          IFATH = -IN
        ENDIF
        IF (IFATH.EQ.0) THEN
           NBROOT = NBROOT - 1
           IF (NBROOT.EQ.0) GOTO 116
           GOTO 91
        ELSE
           fact(STEP(IFATH))=fact(STEP(IFATH))+fact(STEP(INODE))
           IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
              DEPTH(STEP(IFATH))=max(DEPTH(STEP(INODE)),
     &             DEPTH(STEP(IFATH)))
           ENDIF
        ENDIF
        TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
        IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
           INODE = IFATH        
           IN=INODE
           dernier=IN
           I=1
 5700      IN = FILS(IN)
           IF (IN .GT. 0 ) THEN
             dernier=IN
             I=I+1
             GOTO 5700
           ENDIF
           NCB=int(ND(STEP(INODE))-I,8)
           IN=-IN
           IF(PERM.NE.7)THEN
              DO I=1,NE(STEP(INODE))
                 SON(I)=IN
                 TEMP(I)=IN
                 IF(IN.GT.0) IN=FRERE(STEP(IN))
              ENDDO
           ELSE
              DO I=NE(STEP(INODE)),1,-1
                 SON(I)=IN
                 TEMP(I)=IN
                 IF(IN.GT.0) IN=FRERE(STEP(IN))
              ENDDO
           ENDIF
           NFR = int(ND(STEP(INODE)),8)
           DO II=1,NE(STEP(INODE))
             TAB1(II)=0_8
             TAB2(II)=0_8
             cour=SON(II)
             NELIM4=1
 151         cour=FILS(cour)
             IF(cour.GT.0) THEN
                NELIM4=NELIM4+1
                GOTO 151
             ENDIF
             NELIM=int(NELIM4,8)
             IF((SYM.EQ.0).OR.(K215.NE.0)) THEN
                SIZECB=(int(ND(STEP(SON(II))),8)-NELIM)
     &                *(int(ND(STEP(SON(II))),8)-NELIM)
             ELSE
                SIZECB=(int(ND(STEP(SON(II))),8)-NELIM)
     &                *(int(ND(STEP(SON(II))),8)-
     &               NELIM+1_8)/2_8
             ENDIF
             IF((PERM.EQ.0).OR.(PERM.EQ.5))THEN
                IF (K234 .NE. 0 .AND. K55.EQ.0 ) THEN
                   TMP8=NFR
                   TMP8=TMP8*TMP8
                   TAB1(II)=max(TMP8, M(STEP(SON(II)))) - SIZECB
                   TAB2(II)=SIZECB
                ELSE
                   TAB1(II)=M(STEP(SON(II)))- SIZECB
                   TAB2(II)=SIZECB
                ENDIF
             ENDIF
             IF((PERM.EQ.1).OR.(PERM.EQ.6)) THEN
                TAB1(II)=M(STEP(SON(II)))-SIZECB
                TAB1(II)=TAB1(II)-fact(STEP(SON(II)))
                TAB2(II)=SIZECB+fact(STEP(SON(II)))
             ENDIF
             IF(PERM.EQ.2)THEN
                IF (MUMPS_170(PROCNODE(STEP(INODE)),
     &               SLAVEF))THEN
                   TAB1(II)=M_TOTAL(STEP(SON(II)))-SIZECB
     &                  -fact(STEP(SON(II)))
                   TAB2(II)=SIZECB
                ELSE
                   TAB1(II)=M(STEP(SON(II)))-SIZECB
                   TAB2(II)=SIZECB             
                ENDIF
             ENDIF
             IF(PERM.EQ.3)THEN
                IF (MUMPS_170(PROCNODE(STEP(INODE)),
     &               SLAVEF))THEN
                   TAB1(II)=M(STEP(SON(II)))-SIZECB
                   TAB2(II)=SIZECB               
                ELSE
                   TAB1(II)=int(DEPTH(STEP(SON(II))),8)
                   TAB2(II)=M(STEP(SON(II)))
                ENDIF
             ENDIF
             IF(PERM.EQ.4)THEN
                IF (MUMPS_170(PROCNODE(STEP(INODE)),
     &               SLAVEF))THEN
                   TAB1(II)=M(STEP(SON(II)))-
     &                  SIZECB-fact(STEP(SON(II)))
                   TAB2(II)=SIZECB             
                ELSE
                   TAB1(II)=int(DEPTH(STEP(SON(II))),8)
                   TAB2(II)=M(STEP(SON(II)))
                ENDIF
             ENDIF
          ENDDO
          CALL DMUMPS_462(SON,NE(STEP(INODE)),TAB1,TAB2,
     &         LOCAL_PERM
     &           ,RESULT,T1,T2)
          IF(PERM.EQ.0) THEN
             DO II=1,NE(STEP(INODE))
               cour=TEMP(II)
               NELIM4=1
 153           cour=FILS(cour)
               IF(cour.GT.0) THEN
                  NELIM4=NELIM4+1
                  GOTO 153
               ENDIF
               NELIM=int(NELIM4,8)
               IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                  SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                 (int(ND(STEP(TEMP(II))),8)-NELIM)
               ELSE
                  SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                    (int(ND(STEP(TEMP(II))),8)-NELIM+1_8)/2_8
               ENDIF
               TAB1(II)=SIZECB
             ENDDO
             CALL DMUMPS_462(TEMP,NE(STEP(INODE)),TAB1,TAB2,3,
     &         RESULT,T1,T2)             
           ENDIF
           IF(PERM.EQ.1) THEN
              DO II=1,NE(STEP(INODE))
                cour=TEMP(II)
                NELIM4=1
 187            cour=FILS(cour)
                IF(cour.GT.0) THEN
                   NELIM4=NELIM4+1
                   GOTO 187
                ENDIF    
                NELIM=int(NELIM4,8)
                IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                   SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                    (int(ND(STEP(TEMP(II))),8)-NELIM)
                ELSE
                   SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                    (int(ND(STEP(TEMP(II))),8)-NELIM+1_8)/2_8
                ENDIF
                TAB1(II)=SIZECB+fact(STEP(TEMP(II)))
             ENDDO
             CALL DMUMPS_462(TEMP,NE(STEP(INODE)),TAB1,TAB2,3,
     &         RESULT,T1,T2)
           ENDIF
           CONTINUE
           IFATH=INODE
           DO II=1,2
              SUM=0_8
              FACT_SIZE=0_8
              FACT_SIZE_T=0_8
              MEM_SIZE=0_8
              MEM_SIZE_T=0_8
              CB_MAX=0
              CB_current=0
              TMP_SUM=0_8
              IF(II.EQ.1) TAB=>SON 
              IF(II.EQ.2) TAB=>TEMP
              DO I=1,NE(STEP(INODE))
                 cour=TAB(I)
                 NELIM4=1
 149             cour=FILS(cour)
                 IF(cour.GT.0) THEN
                    NELIM4=NELIM4+1
                    GOTO 149
                 ENDIF    
                 NELIM=int(NELIM4, 8)
                 NFR=int(ND(STEP(TAB(I))),8)
                 IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                    SIZECB=(NFR-NELIM)*(NFR-NELIM)
                 ELSE
                    SIZECB=(NFR-NELIM)*(NFR-NELIM+1_8)/2_8
                 ENDIF
                 MEM_SIZE=max(MEM_SIZE,(M(STEP(TAB(I)))+SUM+FACT_SIZE))
                 IF (SBTR_M.OR.(PERM.EQ.2)) THEN
                       MEM_SIZE_T=max(MEM_SIZE_T,(M_TOTAL(STEP(TAB(I)))+
     &                      SUM+
     &                      FACT_SIZE_T))
                       FACT_SIZE_T=FACT_SIZE_T+fact(STEP(TAB(I)))
                 ENDIF
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                (M(STEP(TAB(I)))+SUM+FACT_SIZE))
                 TMP_SUM=TMP_SUM+fact(STEP(TAB(I)))
                 SUM=SUM+SIZECB
                 SIZECB_LASTSON = SIZECB
                 IF((PERM.EQ.1).OR.(PERM.EQ.4))THEN
                    FACT_SIZE=FACT_SIZE+fact(STEP(TAB(I)))
                 ENDIF
              ENDDO
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=NCB*NCB
              ELSE
                 SIZECB=(NCB*(NCB+1_8))/2_8
              ENDIF
              IF (K234.NE.0 .AND. K55.EQ.0) THEN
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                ( (   int(ND(STEP(IFATH)),8)
     &                    * int(ND(STEP(IFATH)),8) )
     &                  + SUM-SIZECB_LASTSON+TMP_SUM )
     &           )
              ELSE IF (K234.NE.0 .AND. K55.NE.0) THEN
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                ( ( int(ND(STEP(IFATH)),8)
     &                  * int(ND(STEP(IFATH)),8) )
     &                  + SUM + TMP_SUM )
     &           )
              ELSE
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                ( ( int(ND(STEP(IFATH)),8)
     &                  * int(ND(STEP(IFATH)),8))
     &                  + max(SUM,SIZECB) + TMP_SUM )
     &                )
              ENDIF
              IF(II.EQ.1)THEN
                 TMP_TOTAL_MEM_SIZE=TOTAL_MEM_SIZE
              ENDIF
              IF(II.EQ.1)THEN
                 IF (K234.NE.0 .AND. K55.EQ.0) THEN
                   M(STEP(IFATH))=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+SUM-SIZECB_LASTSON+
     &             FACT_SIZE))
                 ELSE IF (K234.NE.0 .AND. K55.NE.0) THEN
                   M(STEP(IFATH))=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+SUM+FACT_SIZE))
                 ELSE
                   M(STEP(IFATH))=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+max(SUM,SIZECB)+FACT_SIZE))
                 ENDIF
                 IF (SBTR_M.OR.(PERM.EQ.2))  THEN
                       M_TOTAL(STEP(IFATH))=max(MEM_SIZE_T,
     &                      ((int(ND(STEP(IFATH)),8)
     &                      *int(ND(STEP(IFATH)),8))+max(SUM,SIZECB)+
     &                      FACT_SIZE_T))
                 ENDIF
              ENDIF
              IF((II.EQ.2).AND.(PERM.EQ.1).OR.(PERM.EQ.0).OR.
     &             (PERM.EQ.5).OR.(PERM.EQ.6).OR.
     &             (.NOT.SBTR_M.OR.(SBTR_WHICH_M.NE.1)))THEN
                 MEM_SEC_PERM=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+max(SUM,SIZECB)+FACT_SIZE))
              ENDIF
              IF((PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4))THEN
                 MEM_SEC_PERM=huge(MEM_SEC_PERM)
              ENDIF
           ENDDO
           IF(MEM_SEC_PERM.EQ.M(STEP(IFATH))) THEN
              TAB=>TEMP
           ELSE IF (MEM_SEC_PERM.LT.M(STEP(IFATH))) THEN
              WRITE(*,*)'Probleme dans reorder!!!!'
              CALL MUMPS_ABORT()
           ELSE 
              TOTAL_MEM_SIZE=TMP_TOTAL_MEM_SIZE
              TAB=>SON
           ENDIF
           DO I=NE(STEP(INODE)),1,-1
              IF(I.EQ.NE(STEP(INODE))) THEN
                 FILS(dernier)=-TAB(I)
                 dernier=TAB(I)
                 GOTO 222
              ENDIF
              IF(I.EQ.1) THEN
                 FRERE(STEP(dernier))=TAB(I)
                 FRERE(STEP(TAB(I)))=-INODE
                 GOTO 222
              ENDIF
              IF(I.GT.1) THEN
                 FRERE(STEP(dernier))=TAB(I)
                 dernier=TAB(I)
                 GOTO 222
              ENDIF
 222          CONTINUE
           ENDDO
           GOTO 96
        ELSE
           GOTO 91
        ENDIF
 116    CONTINUE
        NBROOT = NA(2)
        IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        IF (PERM.eq.1) THEN
          DO I=1,NBROOT
            TAB1(I)=M(STEP(NA(I+2+NBLEAF)))-fact(STEP(NA(I+2+NBLEAF)))
            TAB1(I)=-TAB1(I)
          ENDDO
          CALL DMUMPS_462(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &      RESULT,T1,T2)
          IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        ENDIF
 001    CONTINUE
        fin=NBROOT
        LEAF=NA(1)
        FIRST_LEAF=-9999
        SIZE_SBTR=0
 999    CONTINUE
        INODE=IPOOL(fin)
        IF(INODE.LT.0)THEN
           WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
           CALL MUMPS_ABORT()
        ENDIF
        IN=INODE
 5602   IN = FILS(IN)
        IF (IN .GT. 0 ) THEN
           dernier=IN
           GOTO 5602
        ENDIF
        IN=-IN
        IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
           IF(SLAVEF.NE.1)THEN
              IF (USE_DAD) THEN
                 IFATH=DAD(INODE)
              ELSE
                IN = INODE
 395            IN = FRERE(IN)
                IF (IN.GT.0) GO TO 395
                IFATH = -IN
              ENDIF
              NFR4   = ND(STEP(INODE))
              NFR    = int(NFR4,8)
              NELIM4 = 0
              IN = INODE
 396          NELIM4 = NELIM4 + 1
              IN = FILS(IN)
              IF (IN .GT. 0 ) GOTO 396
              NELIM=int(NELIM4,8)
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=(NFR-NELIM)*(NFR-NELIM)
              ELSE
                 SIZECB=(NFR-NELIM)*(NFR-NELIM+1_8)/2_8
              ENDIF
              CALL MUMPS_511(NFR4,NELIM4,NELIM4,
     &             SYM,1,COST_NODE)
              IF(IFATH.NE.0)THEN
                 IF(MUMPS_167(PROCNODE(STEP(INODE)),SLAVEF))THEN
                    COST_TRAV(STEP(INODE))=COST_TRAV(STEP(
     &                   ROOT_OF_CUR_SBTR))
                 ELSE
                    COST_TRAV(STEP(INODE))=dble(COST_NODE)+
     &                   COST_TRAV(STEP(IFATH))+
     &                   dble(SIZECB*18_8)  
                 ENDIF
              ELSE
                 COST_TRAV(STEP(INODE))=dble(COST_NODE)
              ENDIF
           ENDIF
        ENDIF
        DO I=1,NE(STEP(INODE))
           TEMP(I)=IN
           IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
              IF((SLAVEF.NE.1).AND.(.NOT.MUMPS_170(
     &             PROCNODE(STEP(INODE)),SLAVEF)))THEN
                 NFR4   = ND(STEP(INODE))
                 NFR    = int(NFR4,8)
                 NELIM4 = 0
                 II = TEMP(I)
 845             NELIM4 = NELIM4 + 1
                 II = FILS(II)
                 IF (II .GT. 0 ) GOTO 845
                 NELIM=int(NELIM4,8)
                 CALL MUMPS_511(NFR4,NELIM4,NELIM4,
     &                SYM,1,COST_NODE)
                 TAB1(I)=int(dble(COST_NODE)+
     &                COST_TRAV(STEP(INODE)),8)
                 TAB2(I)=0_8
              ELSE
                 SON(I)=IN
              ENDIF
           ELSE
              SON(I)=IN
           ENDIF
           IN=FRERE(STEP(IN))
        ENDDO
        IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
           IF((SLAVEF.NE.1).AND.(.NOT.MUMPS_170(
     &          PROCNODE(STEP(INODE)),SLAVEF)))THEN
              CALL DMUMPS_462(TEMP,NE(STEP(INODE)),TAB1,TAB2,
     &             LOCAL_PERM
     &             ,RESULT,T1,T2)
              TAB=>TEMP
              DO I=NE(STEP(INODE)),1,-1
                 IF(I.EQ.NE(STEP(INODE))) THEN
                    FILS(dernier)=-TAB(I)
                    dernier=TAB(I)
                    GOTO 221
                 ENDIF
                 IF(I.EQ.1) THEN
                    FRERE(STEP(dernier))=TAB(I)
                    FRERE(STEP(TAB(I)))=-INODE
                    GOTO 221
                 ENDIF
                 IF(I.GT.1) THEN
                    FRERE(STEP(dernier))=TAB(I)
                    dernier=TAB(I)
                    GOTO 221
                 ENDIF
 221             CONTINUE
                 SON(NE(STEP(INODE))-I+1)=TAB(I)
              ENDDO
           ENDIF
        ENDIF
        DO I=1,NE(STEP(INODE))
           IPOOL(fin)=SON(I)
           SON(I)=0
           fin=fin+1
        ENDDO
        IF(NE(STEP(INODE)).EQ.0)THEN
           IF(PERM.NE.7)THEN
              NA(LEAF+2)=INODE
           ENDIF
           LEAF=LEAF-1
        ELSE
           fin=fin-1
           GOTO 999
        ENDIF
        fin=fin-1
        IF(fin.EQ.0) THEN
           GOTO 789
        ENDIF
        GOTO 999
 789    CONTINUE
        IF(PERM.EQ.7) GOTO 5483
        NBROOT=NA(2)
        NBLEAF=NA(1)
        PEAK=0.0D0
        FACT_SIZE=0_8
        DO I=1,NBROOT
           PEAK=max(PEAK,dble(M(STEP(NA(2+NBLEAF+I)))))
           FACT_SIZE=FACT_SIZE+fact(STEP(NA(2+NBLEAF+I)))
        ENDDO
 5483   CONTINUE
        DEALLOCATE(IPOOL)
        DEALLOCATE(fact)
        DEALLOCATE(TNSTK)
        DEALLOCATE(SON)
        DEALLOCATE(TAB2)
        DEALLOCATE(TAB1)
        DEALLOCATE(T1)
        DEALLOCATE(T2)
        DEALLOCATE(RESULT)
        DEALLOCATE(TEMP)
        IF(PERM.NE.7)THEN
           DEALLOCATE(M)
        ENDIF
        IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
           DEALLOCATE(DEPTH)
        ENDIF
        IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
           DEALLOCATE(COST_TRAV)
        ENDIF
        IF ((PERM.NE.7).AND.(SBTR_M.OR.(PERM.EQ.2))) THEN
           IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1).OR.(PERM.EQ.2))THEN
              DEALLOCATE(M_TOTAL)
           ENDIF
        ENDIF
      RETURN
      END SUBROUTINE DMUMPS_363
      SUBROUTINE DMUMPS_364(N,FRERE, STEP, FILS,
     &     NA,LNA,NE,ND, DAD, LDAD, USE_DAD,
     &     NSTEPS,PERM,SYM,INFO,LP,K47,K81,K76,K215,K234,K55,
     &     PROCNODE,MEM_SUBTREE,SLAVEF, SIZE_MEM_SBTR, PEAK
     &     ,SBTR_WHICH_M,SIZE_DEPTH_FIRST,SIZE_COST_TRAV,
     &     DEPTH_FIRST_TRAV,DEPTH_FIRST_SEQ,COST_TRAV,MY_FIRST_LEAF,
     &     MY_NB_LEAF,MY_ROOT_SBTR,SBTR_ID
     &     )
      IMPLICIT NONE
      INTEGER N,PERM,SYM, NSTEPS, LNA, LP, SIZE_MEM_SBTR,LDAD
      INTEGER FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER NA(LNA), NE(NSTEPS), ND(NSTEPS)
      INTEGER K47,K81,K76,K215,K234,K55
      INTEGER DAD(LDAD)
      LOGICAL USE_DAD
      INTEGER INFO(40)
      INTEGER SLAVEF,PROCNODE(NSTEPS)
      DOUBLE PRECISION, intent(out) :: MEM_SUBTREE(SIZE_MEM_SBTR,SLAVEF)
      INTEGER :: SBTR_WHICH_M
      INTEGER MY_FIRST_LEAF(SIZE_MEM_SBTR,SLAVEF),
     &     MY_ROOT_SBTR(SIZE_MEM_SBTR,SLAVEF),
     &     MY_NB_LEAF(SIZE_MEM_SBTR,SLAVEF)
      EXTERNAL MUMPS_283,MUMPS_275
      LOGICAL MUMPS_283
      INTEGER MUMPS_275
      DOUBLE PRECISION PEAK
      INTEGER SIZE_DEPTH_FIRST,DEPTH_FIRST_TRAV(SIZE_DEPTH_FIRST),
     &     DEPTH_FIRST_SEQ(SIZE_DEPTH_FIRST)
      INTEGER SIZE_COST_TRAV
      INTEGER SBTR_ID(SIZE_DEPTH_FIRST),OOC_CUR_SBTR
      DOUBLE PRECISION COST_TRAV(SIZE_COST_TRAV)
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH
      INTEGER IFATH,IN,INODE,I,allocok,LOCAL_PERM
      INTEGER(8) NELIM,NFR
      INTEGER NFR4,NELIM4
      INTEGER LEAF,NBLEAF,NBROOT, SIZE_TAB
      INTEGER, DIMENSION (:), ALLOCATABLE :: IPOOL,TNSTK
      INTEGER, DIMENSION (:), ALLOCATABLE,TARGET :: SON,TEMP
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: M,M_TOTAL, fact
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: TAB1,TAB2
      INTEGER x,dernier,fin,RANK_TRAV
      INTEGER II
      INTEGER ROOT_OF_CUR_SBTR
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: T1,T2
      INTEGER, DIMENSION (:), ALLOCATABLE :: RESULT
      INTEGER(8) MEM_SIZE,FACT_SIZE,
     &     TOTAL_MEM_SIZE,
     &     SIZECB
      LOGICAL   SBTR_M
      INTEGER INDICE(SLAVEF),ID,FIRST_LEAF,SIZE_SBTR
      EXTERNAL MUMPS_170,MUMPS_167
      LOGICAL MUMPS_170,MUMPS_167
      DOUBLE PRECISION COST_NODE
      INTEGER CUR_DEPTH_FIRST_RANK
      INCLUDE 'mumps_headers.h'
      TOTAL_MEM_SIZE=0_8
      ROOT_OF_CUR_SBTR=0
      IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.
     &     (PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4).OR.
     &     (PERM.EQ.5).OR.(PERM.EQ.6))THEN
         LOCAL_PERM=0
      ENDIF
      IF (K47 == 4 .OR. ((K47.GE.2).AND.(K81.GE. 1))) THEN
        DO I=1,SLAVEF
          INDICE(I)=1
        ENDDO
        DO I=1,SLAVEF
          DO x=1,SIZE_MEM_SBTR
            MEM_SUBTREE(x,I)=-1.0D0
          ENDDO
        ENDDO
      ENDIF
      SBTR_M=((K47 == 4 .OR. ((K47.GE.2).AND.(K81 .GE. 1))))
      MEM_SIZE=0_8
      FACT_SIZE=0_8
      IF ((PERM.GT.7).AND.
     & (.NOT.(K47 == 4 .OR. ((K47.GE.2).AND.(K81 .GE. 1))))) THEN
         WRITE(*,*) "Internal Error in DMUMPS_363",PERM
         CALL MUMPS_ABORT()
      END IF
      NBLEAF = NA(1)
      NBROOT = NA(2)
      CUR_DEPTH_FIRST_RANK=1
      IF((PERM.EQ.0).AND.(NBROOT.EQ.NBLEAF)) RETURN
      IF (SBTR_M.OR.(PERM.EQ.2))  THEN
         IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
            ALLOCATE(M_TOTAL(NSTEPS), stat=allocok )
            IF (allocok > 0) THEN
               IF ( LP .GT. 0 )
     &              WRITE(LP,*)'Memory allocation error in
     &              DMUMPS_363'
               INFO(1)=-7
               INFO(2)=NSTEPS
               RETURN
            ENDIF
         ENDIF
      ENDIF
      ALLOCATE( IPOOL(NBLEAF), M(NSTEPS), fact(NSTEPS),
     &          TNSTK(NSTEPS), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_363'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      II=0
      DO I=1,NSTEPS
         TNSTK(I) = NE(I)
         IF(NE(I).GE.II) II=NE(I)
      ENDDO
      SIZE_TAB=max(II,NBROOT)
      ALLOCATE(SON(II), TEMP(II),
     &         TAB1(SIZE_TAB), TAB2(SIZE_TAB), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_363'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      ALLOCATE(T1(SIZE_TAB),T2(SIZE_TAB),
     &         RESULT(SIZE_TAB),stat=allocok)
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_363'
        INFO(1)=-7
        INFO(2)=SIZE_TAB
        RETURN
      ENDIF
      IF(NBROOT.EQ.NBLEAF)THEN
        IF((PERM.NE.1).OR.(PERM.EQ.4).OR.(PERM.EQ.6))THEN
          WRITE(*,*)'Internal Error in reordertree:'
          WRITE(*,*)'  problem with perm parameter in reordertree'
          CALL MUMPS_ABORT()
        ENDIF
        DO I=1,NBROOT
          TAB1(I)=int(ND(STEP(NA(I+2+NBLEAF))),8)
          IPOOL(I)=NA(I+2+NBLEAF)
          M(STEP(IPOOL(I)))=TAB1(I)*TAB1(I)
        ENDDO
        CALL DMUMPS_462(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &    RESULT,T1,T2)
        GOTO 789
      ENDIF
      IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
         ALLOCATE(DEPTH(NSTEPS),stat=allocok)
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error in
     &           DMUMPS_363'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
         DEPTH=0
         NBROOT = NA(2)
         IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
         fin=NBROOT
         LEAF=NA(1)
 499     CONTINUE
         INODE=IPOOL(fin)
         IF(INODE.LT.0)THEN
            WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
            CALL MUMPS_ABORT()
         ENDIF
         IN=INODE
 4602    IN = FILS(IN)
         IF (IN .GT. 0 ) THEN
            GOTO 4602
         ENDIF
         IN=-IN
         DO I=1,NE(STEP(INODE))
            SON(I)=IN
            IN=FRERE(STEP(IN))
         ENDDO
         DO I=1,NE(STEP(INODE))
            IPOOL(fin)=SON(I)
            DEPTH(STEP(SON(I)))=DEPTH(STEP(INODE))+1
            SON(I)=0
            fin=fin+1
         ENDDO
         IF(NE(STEP(INODE)).EQ.0)THEN
            LEAF=LEAF-1
         ELSE
            fin=fin-1
            GOTO 499
         ENDIF
         fin=fin-1
         IF(fin.EQ.0) GOTO 489
         GOTO 499
 489     CONTINUE
      ENDIF
      IF(K76.EQ.4.OR.(K76.EQ.6))THEN
         RANK_TRAV=NSTEPS
         DEPTH_FIRST_TRAV=0
         DEPTH_FIRST_SEQ=0
      ENDIF
      IF((K76.EQ.5).OR.(PERM.EQ.5).OR.(PERM.EQ.6))THEN
         COST_TRAV=0.0D0
         COST_NODE=0.0d0
      ENDIF        
      DO I=1,NSTEPS
         M(I)=0_8
         IF (SBTR_M.OR.(PERM.EQ.2))  THEN
            IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
               M_TOTAL(I)=0_8
            ENDIF
         ENDIF
      ENDDO
      DO I=1,NSTEPS
         fact(I)=0_8
      ENDDO
        NBROOT = NA(2)
        NBLEAF = NA(1)
        IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        CONTINUE
        fin=NBROOT
        LEAF=NA(1)
        FIRST_LEAF=-9999
        SIZE_SBTR=0
 999    CONTINUE
        INODE=IPOOL(fin)
        IF(INODE.LT.0)THEN
           WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
           CALL MUMPS_ABORT()
        ENDIF
        IF(SIZE_SBTR.NE.0)THEN 
           IF(.NOT.MUMPS_167(PROCNODE(STEP(INODE)),SLAVEF))THEN
              IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
                 IF((SLAVEF.NE.1))THEN
                    MY_FIRST_LEAF(INDICE(ID+1)-1,ID+1)=FIRST_LEAF
                    MY_NB_LEAF(INDICE(ID+1)-1,ID+1)=SIZE_SBTR
                    FIRST_LEAF=-9999
                    SIZE_SBTR=0
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
        IF(MUMPS_283(PROCNODE(STEP(INODE)),SLAVEF))THEN
           ROOT_OF_CUR_SBTR=INODE
        ENDIF
        IF (K76.EQ.4)THEN
           IF(SLAVEF.NE.1)THEN
              WRITE(*,*)'INODE=',INODE,'RANK',RANK_TRAV
              IF(MUMPS_167(PROCNODE(STEP(INODE)),SLAVEF))THEN
                 DEPTH_FIRST_TRAV(STEP(INODE))=DEPTH_FIRST_TRAV(STEP(
     &                ROOT_OF_CUR_SBTR))
              ELSE
                 DEPTH_FIRST_TRAV(STEP(INODE))=RANK_TRAV
              ENDIF
              RANK_TRAV=RANK_TRAV-1
           ENDIF
        ENDIF
        IF (K76.EQ.5)THEN
           IF(SLAVEF.NE.1)THEN
              IF (USE_DAD) THEN
                IFATH=DAD(INODE)
              ELSE
                IN = INODE
 395            IN = FRERE(IN)
                IF (IN.GT.0) GO TO 395
                IFATH = -IN
              ENDIF
              NFR4   = ND(STEP(INODE))
              NFR    = int(NFR4,8)
              NELIM4 = 0
              IN = INODE
 396          NELIM4 = NELIM4 + 1
              IN = FILS(IN)
              IF (IN .GT. 0 ) GOTO 396
              NELIM=int(NELIM4,8)
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=(NFR-NELIM)*(NFR-NELIM)
              ELSE
                 SIZECB=(NFR-NELIM)*(NFR-NELIM+1_8)/2_8
              ENDIF
              CALL MUMPS_511(NFR4,NELIM4,NELIM4,
     &             SYM,1,COST_NODE)
              IF(IFATH.NE.0)THEN
                 IF(MUMPS_167(PROCNODE(STEP(INODE)),SLAVEF))THEN
                    COST_TRAV(STEP(INODE))=COST_TRAV(STEP(
     &                   ROOT_OF_CUR_SBTR))
                 ELSE
                    COST_TRAV(STEP(INODE))=dble(COST_NODE)+
     &                   COST_TRAV(STEP(IFATH))+
     &                   dble(SIZECB*18_8)  
                 ENDIF
              ELSE
                 COST_TRAV(STEP(INODE))=dble(COST_NODE)
              ENDIF
              IF(K76.EQ.5)THEN
                 WRITE(*,*)'INODE=',INODE,'COST=',COST_TRAV(STEP(INODE))
              ENDIF
           ENDIF
        ENDIF
        IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
              IF((SLAVEF.NE.1).AND.
     &          MUMPS_283(PROCNODE(STEP(INODE)),SLAVEF))THEN
                IF (NE(STEP(INODE)).NE.0) THEN
                   ID=MUMPS_275(PROCNODE(STEP(INODE)),SLAVEF)
                   IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
                      MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                     dble(M_TOTAL(STEP(INODE)))
                   ELSE
                      MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                     dble(M(STEP(INODE)))
                   ENDIF
                   MY_ROOT_SBTR(INDICE(ID+1),ID+1)=INODE
                  INDICE(ID+1)=INDICE(ID+1)+1
                ENDIF
              ENDIF
              IF((SLAVEF.EQ.1).AND.FRERE(STEP(INODE)).EQ.0)THEN
                 ID=MUMPS_275(PROCNODE(STEP(INODE)),SLAVEF)
                 IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
                    MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                   dble(M_TOTAL(STEP(INODE)))
                 ELSE
                    MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                   dble(M(STEP(INODE)))
                 ENDIF
                 INDICE(ID+1)=INDICE(ID+1)+1
              ENDIF
        ENDIF
        IN=INODE
 5602   IN = FILS(IN)
        IF (IN .GT. 0 ) THEN
           dernier=IN
           GOTO 5602
        ENDIF
        IN=-IN
        DO I=1,NE(STEP(INODE))
           IPOOL(fin)=IN
           IF(IN.GT.0) IN=FRERE(STEP(IN))
           fin=fin+1
        ENDDO
        IF(NE(STEP(INODE)).EQ.0)THEN
           IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
              IF(SLAVEF.NE.1)THEN
                 IF(MUMPS_167(PROCNODE(STEP(INODE)),SLAVEF))THEN
                    IF(FIRST_LEAF.EQ.-9999)THEN
                       FIRST_LEAF=INODE
                    ENDIF
                    SIZE_SBTR=SIZE_SBTR+1
                 ENDIF
              ENDIF
           ENDIF
           IF(PERM.NE.7)THEN
              NA(LEAF+2)=INODE
           ENDIF
           LEAF=LEAF-1
        ELSE
           fin=fin-1
           GOTO 999
        ENDIF
        fin=fin-1
        IF(fin.EQ.0) THEN
           IF(SIZE_SBTR.NE.0)THEN 
              IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
                 IF((SLAVEF.NE.1))THEN
                    MY_FIRST_LEAF(INDICE(ID+1)-1,ID+1)=FIRST_LEAF
                    MY_NB_LEAF(INDICE(ID+1)-1,ID+1)=SIZE_SBTR
                    FIRST_LEAF=-9999
                    SIZE_SBTR=0
                 ENDIF
              ENDIF
           ENDIF
           GOTO 789
        ENDIF
        GOTO 999
 789    CONTINUE
        IF(K76.EQ.6)THEN
           OOC_CUR_SBTR=1
           DO I=1,NSTEPS
              TNSTK(I) = NE(I)
           ENDDO           
           NBROOT=NA(2)
           NBLEAF=NA(1)
           IPOOL(1:NBLEAF)=NA(3:2+NBLEAF)
           LEAF = NBLEAF + 1
 9100      CONTINUE
           IF (LEAF.NE.1) THEN
              LEAF = LEAF -1
              INODE = IPOOL(LEAF)
           ENDIF
 9600      CONTINUE
           IF(SLAVEF.NE.1)THEN
              ID=MUMPS_275(PROCNODE(STEP(INODE)),SLAVEF)
              DEPTH_FIRST_TRAV(STEP(INODE))=CUR_DEPTH_FIRST_RANK
              DEPTH_FIRST_SEQ(CUR_DEPTH_FIRST_RANK)=INODE
              WRITE(*,*)ID,': INODE -> ',INODE,'DF =',
     &             CUR_DEPTH_FIRST_RANK
              CUR_DEPTH_FIRST_RANK=CUR_DEPTH_FIRST_RANK+1
              IF(MUMPS_170(PROCNODE(STEP(INODE)),
     &             SLAVEF))THEN
                 SBTR_ID(STEP(INODE))=OOC_CUR_SBTR
              ELSE
                 SBTR_ID(STEP(INODE))=-9999
              ENDIF
              IF(MUMPS_283(PROCNODE(STEP(INODE)),
     &             SLAVEF))THEN
                 OOC_CUR_SBTR=OOC_CUR_SBTR+1
              ENDIF
           ENDIF
           IF (USE_DAD) THEN
              IFATH = DAD( STEP(INODE) )
           ELSE
              IN = INODE
 1133         IN = FRERE(IN)
              IF (IN.GT.0) GO TO 1133
              IFATH = -IN
           ENDIF
           IF (IFATH.EQ.0) THEN
              NBROOT = NBROOT - 1
              IF (NBROOT.EQ.0) GOTO 1163
              GOTO 9100
           ENDIF
           TNSTK(STEP(IFATH))=TNSTK(STEP(IFATH))-1
           IF(TNSTK(STEP(IFATH)).EQ.0) THEN
              INODE=IFATH
              GOTO 9600
           ELSE
              GOTO 9100
           ENDIF        
 1163      CONTINUE
        ENDIF
        PEAK=0.0D0
        FACT_SIZE=0_8
        DO I=1,NBROOT
           PEAK=max(PEAK,dble(M(STEP(NA(2+NBLEAF+I)))))
           FACT_SIZE=FACT_SIZE+fact(STEP(NA(2+NBLEAF+I)))
        ENDDO
        CONTINUE
        DEALLOCATE(IPOOL)
        DEALLOCATE(M)
        DEALLOCATE(fact)
        DEALLOCATE(TNSTK)
        DEALLOCATE(SON)
        DEALLOCATE(TAB2)
        DEALLOCATE(TAB1)
        DEALLOCATE(T1)
        DEALLOCATE(T2)
        DEALLOCATE(RESULT)
        DEALLOCATE(TEMP)
        IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
           DEALLOCATE(DEPTH)
        ENDIF
        IF (SBTR_M.OR.(PERM.EQ.2))  THEN
           IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1).OR.(PERM.EQ.2))THEN
              DEALLOCATE(M_TOTAL)
           ENDIF
        ENDIF
      RETURN
      END SUBROUTINE DMUMPS_364
      RECURSIVE SUBROUTINE DMUMPS_462(TAB,DIM,TAB1,TAB2,PERM,
     &  RESULT,TEMP1,TEMP2)
      IMPLICIT NONE
      INTEGER DIM
      INTEGER(8) TAB1(DIM),TAB2(DIM)
      INTEGER(8) TEMP1(DIM),TEMP2(DIM)
      INTEGER TAB(DIM), PERM,RESULT(DIM)
      INTEGER I,J,I1,I2
      IF(DIM.EQ.1) THEN
        RESULT(1)=TAB(1)
        TEMP1(1)=TAB1(1)
        TEMP2(1)=TAB2(1)
        RETURN
      ENDIF
      I=DIM/2
      CALL DMUMPS_462(TAB(1),I,TAB1(1),TAB2(1),PERM,
     &  RESULT(1),TEMP1(1),TEMP2(1))
      CALL DMUMPS_462(TAB(I+1),DIM-I,TAB1(I+1),TAB2(I+1),
     &  PERM,RESULT(I+1),TEMP1(I+1),TEMP2(I+1))
      I1=1
      I2=I+1
      J=1
      DO WHILE ((I1.LE.I).AND.(I2.LE.DIM))
        IF((PERM.EQ.3))THEN
          IF(TEMP1(I1).LE.TEMP1(I2))THEN
            TAB(J)=RESULT(I1)
            TAB1(J)=TEMP1(I1)
            J=J+1
            I1=I1+1
          ELSE
            TAB(J)=RESULT(I2)
            TAB1(J)=TEMP1(I2)
            J=J+1
            I2=I2+1
          ENDIF
          GOTO 3
        ENDIF
        IF((PERM.EQ.4).OR.(PERM.EQ.5))THEN
          IF (TEMP1(I1).GE.TEMP1(I2))THEN
            TAB(J)=RESULT(I1)
            TAB1(J)=TEMP1(I1)
            J=J+1
            I1=I1+1
          ELSE
            TAB(J)=RESULT(I2)
            TAB1(J)=TEMP1(I2)
            J=J+1
            I2=I2+1          
          ENDIF
          GOTO 3
        ENDIF
        IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.(PERM.EQ.2)) THEN
          IF(TEMP1(I1).GT.TEMP1(I2))THEN
            TAB1(J)=TEMP1(I1)
            TAB2(J)=TEMP2(I1)
            TAB(J)=RESULT(I1)
            J=J+1
            I1=I1+1
            GOTO 3
          ENDIF
          IF(TEMP1(I1).LT.TEMP1(I2))THEN
            TAB1(J)=TEMP1(I2)
            TAB2(J)=TEMP2(I2)
            TAB(J)=RESULT(I2)
            J=J+1
            I2=I2+1
            GOTO 3
          ENDIF        
          IF((TEMP1(I1).EQ.TEMP1(I2)))THEN
            IF(TEMP2(I1).LE.TEMP2(I2))THEN
              TAB1(J)=TEMP1(I1)
              TAB2(J)=TEMP2(I1)
              TAB(J)=RESULT(I1)
              J=J+1
              I1=I1+1
            ELSE
              TAB1(J)=TEMP1(I2)
              TAB2(J)=TEMP2(I2)
              TAB(J)=RESULT(I2)
              J=J+1
              I2=I2+1
            ENDIF
          ENDIF
        ENDIF
  3   CONTINUE    
      ENDDO
      IF(I1.GT.I)THEN
        DO WHILE(I2.LE.DIM)
          TAB(J)=RESULT(I2)
          TAB1(J)=TEMP1(I2)
          TAB2(J)=TEMP2(I2)
          J=J+1
          I2=I2+1
        ENDDO
      ELSE
        IF(I2.GT.DIM)THEN
          DO WHILE(I1.LE.I)
            TAB1(J)=TEMP1(I1)
            TAB2(J)=TEMP2(I1)
            TAB(J)=RESULT(I1)
            J=J+1
            I1=I1+1
          ENDDO
        ENDIF
      ENDIF
      DO I=1,DIM
        TEMP1(I)=TAB1(I)
        TEMP2(I)=TAB2(I)
        RESULT(I)=TAB(I)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_462
