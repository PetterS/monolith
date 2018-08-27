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
C $Id$
      SUBROUTINE MUMPS_419 (METRIC, JOB, COMPRESS, N, NBBUCK, 
     &                   IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, 
     &                   WF, 
     &                   NEXT, W, HEAD, AGG4, 
     &                   SIZE_COMPLEM_LIST,
     &                   COMPLEM_LIST)
      IMPLICIT NONE
      INTEGER, intent(in) :: METRIC, JOB, N, NBBUCK
      LOGICAL, intent(in) :: COMPRESS
      INTEGER IWLEN, PFREE, LEN(N), 
     &        ELEN(N), LAST(N), NCMPA, DEGREE(N), NEXT(N),
     &        W(N) 
      INTEGER PE(N), IW(IWLEN), NV(N)
      LOGICAL, intent(in) :: AGG4
      INTEGER, intent(in) :: SIZE_COMPLEM_LIST
      INTEGER, intent(in), optional :: 
     &                     COMPLEM_LIST(max(1,SIZE_COMPLEM_LIST))
      INTEGER HEAD(0:NBBUCK+1), WF(N)
      INTEGER AMD, AMF1, AMF4MA41
      PARAMETER (AMD=1, AMF1=2, AMF4MA41=4)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        LASTD, NELME, N2, PAS
       INTEGER MAXINT_N
       INTEGER WF3, WF4
       INTEGER(8) HASH, HMOD
       DOUBLE PRECISION RMF, RMF1 
       DOUBLE PRECISION dummy
       INTEGER idummy
       LOGICAL SchurON
       LOGICAL NOTDEFINEDAMD  
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod, huge
      INTEGER TOTEL
        NOTDEFINEDAMD = (METRIC.NE.AMD)
        IF (N.EQ.1) THEN
           ELEN(1) = 1
           LAST(1) = 1
           PE(1) = 0
           NV(1) = 1
           RETURN
        ENDIF
        IF (.NOT.present(COMPLEM_LIST)) SchurON=.FALSE.
        IF ( SIZE_COMPLEM_LIST < 0 .OR. SIZE_COMPLEM_LIST > N ) THEN
          WRITE(*,*) 'Internal MUMPS_419 ', SIZE_COMPLEM_LIST,N
          CALL MUMPS_ABORT()
        ENDIF
        IF (JOB.EQ.2) THEN
          SchurON = .FALSE.        
        ENDIF
        IF (JOB.NE.2) THEN
          SchurON   = (SIZE_COMPLEM_LIST > 0)
          IF ((JOB.EQ.1) .AND. (.NOT.SchurON) ) THEN
           WRITE(6,*) ' WARNING MUMPS_419 on Options ', JOB
          ENDIF
        ENDIF
      idummy = huge(idummy) - 1
      dummy = dble(idummy)
      N2 = -NBBUCK-1
      PAS = max((N/8), 1)
      WFLG = 2
      MAXINT_N = huge(MAXINT_N) - N
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, NBBUCK-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DEGREE(1:N)      = LEN(1:N)
      LAST             = 0
      W(1:N)           = 1
      TOTEL            = N 
      IF (.NOT.COMPRESS) THEN
       NV = 1
      ENDIF
      IF (JOB.EQ.2) THEN
        DO I = 1,SIZE_COMPLEM_LIST
             X       = COMPLEM_LIST(I)
             ELEN(X) = -I       
             NV(X)   = LEN(X)+1 
        ENDDO
        NEL = NEL + SIZE_COMPLEM_LIST  
      ELSE
       ELEN(1:N)        = 0
       DO K=1, SIZE_COMPLEM_LIST
        I = COMPLEM_LIST(K)
        DEGREE(I) = N2
        IF ((LEN(I) .EQ.0).OR.(LEN(I).EQ.-N-1)) THEN
         PE (I)     = 0
         LEN(I)     = 0
        ENDIF
        DEG = NBBUCK + 1
        IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
        ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
        ENDIF
       ENDDO
      ENDIF
      IF(COMPRESS) THEN
         TOTEL = 0
         DO I=1,N
            IF (ELEN(I).LT.0) CYCLE   
            IF (DEGREE(I).NE.N2) THEN 
               TOTEL = TOTEL + NV(I)
               DEGREE(I) = ELEN(I)
               DO J= PE(I)+ELEN(I), PE(I)+LEN(I)-1 
                  DEGREE(I) = DEGREE(I) + NV(IW(J))
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      RMF = dble(0)
      DO I = 1, N
        IF (ELEN(I).LT.0) CYCLE   
        DEG = DEGREE (I)
        IF (DEG.EQ.N2) CYCLE
        IF (DEG .GT. 0) THEN
          IF (JOB.EQ.2) THEN
            DEG = DEG - ELEN(I) 
            NVI = NV(I)
            RMF = dble(0)
            IF (ELEN(I).GT.0) THEN
             DO  J= PE(I), PE(I)+ELEN(I)-1
              DEG = DEG + LEN(IW(J)) - NVI  
              IF (NOTDEFINEDAMD) THEN
               RMF1 = dble( LEN(IW(J)))
               RMF1 = (RMF1-dble(NVI))*(RMF1-dble(NVI)-1.0D0)
               RMF = max(RMF, RMF1)
              ENDIF   
             ENDDO
              DEG = min(DEG, TOTEL-NEL-NV(I))  
            ENDIF
          ENDIF
            IF (
     &          ( (JOB.EQ.2).AND.NOTDEFINEDAMD) 
     &         .OR. (METRIC.EQ.AMF4MA41) 
     &         ) THEN
                DEG   = nint ( 
     &              ( (dble(DEG)*dble(DEG-1)) - RMF )  / dble(2) ) 
                DEG = max (DEG,1)
            ENDIF
            IF (NOTDEFINEDAMD) THEN
              WF(I) = DEG
              IF (DEG.GT.N) THEN
               DEG = min(((DEG-N)/PAS) + N , NBBUCK)
              ENDIF
            ELSE 
              DEGREE(I) = DEG
            ENDIF   
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
      ENDDO
      NLEFT = TOTEL-NEL
   30 IF ( ((NEL .LT. TOTEL).AND. (JOB.NE.1)) .OR. 
     &     ((JOB.EQ.1).AND.(NEL.LT.TOTEL-SIZE_COMPLEM_LIST))
     &   ) THEN
        DO 40 DEG = MINDEG, NBBUCK
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (ME.LE.0) THEN
          NCMPA = -N
          CALL MUMPS_ABORT()
        ENDIF
        IF (DEG.GT.N) THEN
      IF (NOTDEFINEDAMD) THEN
         J = NEXT(ME)
         K = WF(ME)
   55    CONTINUE
         IF (J.GT.0) THEN
          IF (WF(J).LT.K) THEN
           ME = J
           K  = WF(ME)
          ENDIF
          J= NEXT(J)
          GOTO 55
         ENDIF
         ILAST = LAST(ME)
         INEXT = NEXT(ME)
         IF (INEXT .NE. 0) LAST (INEXT) = ILAST
         IF (ILAST .NE. 0) THEN
           NEXT (ILAST) = INEXT
         ELSE
           HEAD (DEG) = INEXT
         ENDIF
      ELSE 
           WRITE(6,*) ' Internal error AMD, DEG>N '
           CALL MUMPS_ABORT()
      ENDIF
         ELSE
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).NE.N2) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
               IF (NOTDEFINEDAMD) THEN
                IF (WF(I).GT.N) THEN
                 DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                ELSE
                 DEG = WF(I)
                ENDIF
                HEAD (DEG) = INEXT
               ELSE
                HEAD (DEGREE (I)) = INEXT
               ENDIF
              ENDIF
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                 IF (NOTDEFINEDAMD) THEN
                  IF (WF(I).GT.N) THEN
                   DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
                  HEAD (DEG) = INEXT
                 ELSE
                  HEAD(DEGREE(I)) = INEXT
                 ENDIF
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max(MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).EQ.N2) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
                IF (NOTDEFINEDAMD) WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).EQ.N2) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG  = 0
          IF (NOTDEFINEDAMD) THEN
            WF3  = 0
            WF4  = 0
          ENDIF
          NVI  = -NV(I)
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
             IF (NOTDEFINEDAMD) THEN
              IF ( WF(E) .EQ. 0 ) THEN
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
             ENDIF
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E, kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
             IF (.NOT.AGG4) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
             ELSE
              PE (E) = -ME
              W (E) = 0
             ENDIF
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IF (NOTDEFINEDAMD) WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (DEGREE(I).EQ.N2) DEG = N2
          IF ( (AGG4.AND.(DEG .EQ. 0)).OR.
     &       (ELEN(I).EQ.1 .AND. P3.EQ.PN) ) THEN
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
                 IF ( DEGREE (I).LT.DEG ) THEN
                  IF (NOTDEFINEDAMD) THEN
                   WF4 = 0
                   WF3 = 0
                  ENDIF
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            IF (NOTDEFINEDAMD) THEN
             WF(I)      = WF4 + 2*NVI*WF3
            ENDIF
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                IF (NOTDEFINEDAMD) WF(I)  = max(WF(I),WF(J))
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).NE.N2) THEN
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
       IF (NOTDEFINEDAMD) THEN
            IF(METRIC.EQ.AMF1) THEN
             DEGREE(I) = DEG
             RMF = dble(DEG)*dble(DEG-1) 
     &         -  dble(DEGME-NVI)*dble(DEGME-NVI-1)
            ELSE
             IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
              DEG = DEGREE(I)
              RMF1  = dble(DEG)*dble( (DEG-1) + 2*DEGME )
     &              - dble(WF(I))
              DEGREE(I) = NLEFT - NVI
              DEG       = DEGREE(I) 
              RMF = dble(DEG)*dble(DEG-1) 
     &         -  dble(DEGME-NVI)*dble(DEGME-NVI-1)
              RMF = min(RMF, RMF1)
             ELSE 
              DEG = DEGREE(I)
              DEGREE(I) = DEGREE (I) + DEGME - NVI
              RMF  = dble(DEG)*dble( (DEG-1) + 2*DEGME ) 
     &              - dble(WF(I))
             ENDIF
            ENDIF
            IF (METRIC.EQ.AMF4MA41) THEN
                RMF =  RMF / dble(2*NVI)   
            ELSE
                 RMF =  RMF / dble(NVI+1) 
            ENDIF
            IF (RMF.LT.dummy) THEN
             WF(I) = int ( anint( RMF ))
            ELSEIF (RMF / dble(N) .LT. dummy) THEN 
             WF(I) = int ( anint( RMF/dble(N) ))
            ELSE
             WF(I) = idummy
            ENDIF
            WF(I) = max(1,WF(I))
            DEG = WF(I)
            IF (DEG.GT.N) THEN
              DEG = min(((DEG-N)/PAS) + N , NBBUCK)
            ENDIF
       ELSE  
            DEGREE(I) = DEG
       ENDIF
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = min (MINDEG, DEG)
              ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      IF (NEL.LT.TOTEL) THEN 
         IF (JOB.EQ.1) THEN
            DO I = 1,SIZE_COMPLEM_LIST
             X       = COMPLEM_LIST(I)
             ELEN(X) = -(N-SIZE_COMPLEM_LIST+I)   
             NV(X)   = 1                
             PE(X)   = 0                
            ENDDO
            NEL = NEL+ SIZE_COMPLEM_LIST
         ELSE
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N2) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
           ENDDO
           ELEN(ME) = NELME
           NV(ME)   = SIZE_COMPLEM_LIST
           PE(ME)   = 0
         ENDIF
        IF (NEL.NE.N) THEN
         write(*,*) ' Error 2 in HALO AMD NEL, N=', NEL,N
         NCMPA = -N - 1
         CALL MUMPS_ABORT()   
        ENDIF
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      IF(COMPRESS) THEN
        LAST(1:N) = 0
        DEGREE(1:TOTEL-N)=0
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            DEGREE(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (DEGREE(K-N) .NE. 0) THEN
            LAST(I)=DEGREE(K-N)
            ELEN(DEGREE(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
        DO 300 I = 1, N
           K = abs (ELEN (I))
           LAST (K) = I
           ELEN (I) = K
300     CONTINUE
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_419 
      SUBROUTINE MUMPS_197(N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      INTEGER N, IWLEN, PFREE, NCMPA
      INTEGER NEXT(N), LEN(N), 
     &        ELEN(N), LAST(N), DEGREE(N), HEAD(N), 
     &        W(N)
      INTEGER IW(IWLEN), NV(N), PE(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER MAXINT_N
      INTEGER(8) HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE (I) = LEN (I)
   10 CONTINUE
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG1)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
#if defined (NOAGG1)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE (I) = min (DEGREE (I), DEG)
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = min (MINDEG, DEG)
            DEGREE (I) = DEG
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_197
      SUBROUTINE MUMPS_23(N,IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      INTEGER N, IWLEN, PFREE, NCMPA
      INTEGER PE(N), LEN(N), 
     &        ELEN(N), LAST(N), DEGREE(N), HEAD(N), 
     &        W(N)
      INTEGER IW(IWLEN), NV(N), NEXT(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X, 
     &        NPRINC
      INTEGER MAXINT_N
      INTEGER(8) HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      NPRINC = 0
      DO I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
      ENDDO
      DO I=1, N
        IF (LEN (I).GE.0) THEN
           DEGREE (I) = LEN (I)
           NPRINC = NPRINC + 1
        ELSE
           J        = -LEN (I)
           DEGREE (I) = - 1
           IF ( PE(I) .NE. 0 ) THEN
             LEN (I) = LEN(J)
           ELSE
             LEN (I) = 0
           ENDIF
           PE (I)   = -J
           NV (J)   = NV (J) + NV (I)
           NV (I)   = 0
           ELEN (I) = 0
        ENDIF
      ENDDO
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE IF ( DEG.EQ. 0) THEN
          ELEN (I) = - (NEL + 1)
          NEL = NEL + NV(I)
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG2)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
#if defined (NOAGG2)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE (I) = min (DEGREE (I), DEG)
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = min (MINDEG, DEG)
            DEGREE (I) = DEG
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_23
      SUBROUTINE MUMPS_162(N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     &                   LISTVAR_SCHUR, SIZE_SCHUR)
      INTEGER SIZE_SCHUR
      INTEGER LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER N, IWLEN, PFREE, NCMPA
      INTEGER LEN(N), 
     &        ELEN(N), LAST(N), DEGREE(N), HEAD(N), 
     &        W(N), NEXT(N)
      INTEGER IW(IWLEN), NV(N), PE(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        NBFLAG, NREAL, LASTD, NELME
      INTEGER MAXINT_N
      INTEGER(8) HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      NBFLAG = 0
      LASTD  = 0
      DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE(I) = LEN(I)
   10 CONTINUE
      NBFLAG = SIZE_SCHUR
      DO K=1,SIZE_SCHUR
       I = LISTVAR_SCHUR(K)
       DEGREE(I) = N+1
       IF ((LEN(I) .EQ.0).OR.(LEN(I).EQ.-N-1)) THEN
         PE (I)     = 0
         LEN(I)     = 0
       ENDIF
       DEG = N
       IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
       ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
       ENDIF
      ENDDO
      NREAL = N - NBFLAG
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N+1)  GOTO 20
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
      NLEFT = N-NEL
   30 IF (NEL .LT. NREAL) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (ME.LE.0) THEN
          write (*,*) ' Error 1 in HALO_AMD '
          NCMPA = -N
          GOTO 500
        ENDIF
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.N) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG3)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (DEGREE(I).EQ.N+1) DEG = N+1
#if defined (NOAGG3)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            IF (DEGREE(I).NE.N+1) THEN
                 DEG        = min (DEG, NLEFT)
                 DEGREE (I) = min (DEGREE (I), DEG)
            ENDIF
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            IF (DEG.LE.N) THEN
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH, kind=kind(LAST))
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I) .LE. N) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.N) THEN
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = min (MINDEG, DEG)
            DEGREE (I) = DEG
              ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      IF (NEL.LT.N) THEN 
           DO DEG = MINDEG, N
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           IF (ME.NE.LISTVAR_SCHUR(1)) THEN
             write(6,*) ' error 1 in MUMPS_162 '
             write(6,*) ' wrong principal var for Schur !!'
             NCMPA = -N - 2
             CALL MUMPS_ABORT()
           ENDIF
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
           ENDDO
           ELEN(ME) = NELME
           NV(ME)   = N-NREAL
           PE(ME)   = 0
        IF (NEL.NE.N) THEN
         write(*,*) ' Error 2 in MUMPS_162 NEL, N=', NEL,N
         NCMPA = -N - 1
         CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
 500  PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_162
      SUBROUTINE MUMPS_337(N, NBBUCK, 
     &                   IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, WF, NEXT, W, HEAD)
      IMPLICIT NONE
      INTEGER N, IWLEN, PFREE, LEN(N), 
     &        ELEN(N), LAST(N), NCMPA, DEGREE(N), NEXT(N),
     &        W(N) 
      INTEGER PE(N), IW(IWLEN), NV(N)
      INTEGER NBBUCK
      INTEGER HEAD(0:NBBUCK+1), WF(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        NBFLAG, NREAL, LASTD, NELME, WF3, WF4, N2, PAS
       INTEGER MAXINT_N
       INTEGER(8) HASH, HMOD
       DOUBLE PRECISION RMF, RMF1 
       DOUBLE PRECISION dummy
       INTEGER idummy
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod, huge
      INTEGER TOTEL
      LOGICAL COMPRESS
      idummy = huge(idummy) - 1
      dummy = dble(idummy)
      N2 = -NBBUCK-1
      PAS = max((N/8), 1)
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, NBBUCK-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
      NBFLAG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DO 10 I = 1, N
        LAST(I) = 0
        W(I) = 1
        ELEN (I) = 0
   10 CONTINUE
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF(COMPRESS) THEN
         TOTEL = 0
         DO I=1,N
            IF (LEN(I).LT.0) THEN 
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               IF (LEN(I).EQ.-N-1) THEN
                  LEN (I)    = 0
                  PE (I)     = 0
               ELSE
                  LEN (I)    = - LEN(I)
               ENDIF
            ELSE
               TOTEL = TOTEL + NV(I)
               DEGREE(I) = 0
               DO J= PE(I) , PE(I)+LEN(I)-1 
                  DEGREE(I) = DEGREE(I) + NV(IW(J))
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            IF (LEN(I).LT.0) THEN 
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               IF (LEN(I).EQ.-N-1) THEN
                  LEN (I)    = 0
                  PE (I)     = 0
               ELSE
                  LEN (I)    = - LEN(I)
               ENDIF
            ELSE
               DEGREE (I) = LEN (I)
            ENDIF
         ENDDO
         TOTEL = N - NBFLAG
      ENDIF
      NREAL = N - NBFLAG
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N2) THEN
             DEG = NBBUCK + 1
             IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
         GOTO 20
        ENDIF
        IF (DEG .GT. 0) THEN
          WF(I) = DEG
           IF (DEG.GT.N) THEN
            DEG = min(((DEG-N)/PAS) + N , NBBUCK)
           ENDIF
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
      NLEFT = TOTEL-NEL
   30 IF (NEL .LT. TOTEL) THEN
        DO 40 DEG = MINDEG, NBBUCK
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (ME.LE.0) THEN
          NCMPA = -N
          CALL MUMPS_ABORT()
        ENDIF
       IF (DEG.GT.N) THEN
         J = NEXT(ME)
         K = WF(ME)
   55    CONTINUE
         IF (J.GT.0) THEN
          IF (WF(J).LT.K) THEN
           ME = J
           K  = WF(ME)
          ENDIF
          J= NEXT(J)
          GOTO 55
         ENDIF
         ILAST = LAST(ME)
         INEXT = NEXT(ME)
         IF (INEXT .NE. 0) LAST (INEXT) = ILAST
         IF (ILAST .NE. 0) THEN
           NEXT (ILAST) = INEXT
         ELSE
           HEAD (DEG) = INEXT
         ENDIF
        ELSE
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).NE.N2) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                IF (WF(I).GT.N) THEN
                 DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                ELSE
                 DEG = WF(I)
                ENDIF
                HEAD (DEG) = INEXT
              ENDIF
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  IF (WF(I).GT.N) THEN
                   DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
                  HEAD (DEG) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
                WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG  = 0
          WF3  = 0
          WF4  = 0
          NVI  = -NV(I)
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              IF ( WF(E) .EQ. 0 ) THEN
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E, kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG4)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (DEGREE(I).EQ.N2) DEG = N2
#if defined (NOAGG4)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            IF (DEGREE(I).NE.N2) THEN
                 IF ( DEGREE (I).LT.DEG ) THEN
                   WF4 = 0
                   WF3 = 0
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            ENDIF
            WF(I)      = WF4 + 2*NVI*WF3
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            IF (DEG.NE.N2) THEN
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                WF(I)  = max(WF(I),WF(J))
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).NE.N2) THEN
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
              DEG = DEGREE(I)
              RMF1  = dble(DEG)*dble( (DEG-1) + 2*DEGME )
     &              - dble(WF(I))
              DEGREE(I) = NLEFT - NVI
              DEG       = DEGREE(I) 
              RMF = dble(DEG)*dble(DEG-1) 
     &         -  dble(DEGME-NVI)*dble(DEGME-NVI-1)
              RMF = min(RMF, RMF1)
            ELSE 
              DEG = DEGREE(I)
              DEGREE(I) = DEGREE (I) + DEGME - NVI
              RMF  = dble(DEG)*dble( (DEG-1) + 2*DEGME ) 
     &              - dble(WF(I))
            ENDIF
            RMF =  RMF / dble(NVI+1)
            IF (RMF.LT.dummy) THEN
             WF(I) = int ( anint( RMF ))
            ELSEIF (RMF / dble(N) .LT. dummy) THEN 
             WF(I) = int ( anint( RMF/dble(N) ))
            ELSE
             WF(I) = idummy
            ENDIF
            WF(I) = max(1,WF(I))
            DEG = WF(I)
            IF (DEG.GT.N) THEN
              DEG = min(((DEG-N)/PAS) + N , NBBUCK)
            ENDIF
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = min (MINDEG, DEG)
              ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      IF (NEL.LT.N) THEN 
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N2) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
           ENDDO
           ELEN(ME) = NELME
           NV(ME)   = N-NREAL
           PE(ME)   = 0
        IF (NEL.NE.N) THEN
         NCMPA = -N - 1
         GOTO 500
        ENDIF
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      IF(COMPRESS) THEN
        LAST(1:N) = 0
        DEGREE(1:TOTEL-N)=0
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            DEGREE(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (DEGREE(K-N) .NE. 0) THEN
            LAST(I)=DEGREE(K-N)
            ELEN(DEGREE(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
        DO 300 I = 1, N
           K = abs (ELEN (I))
           LAST (K) = I
           ELEN (I) = K
300     CONTINUE
      ENDIF
 500  PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_337
      SUBROUTINE MUMPS_421 
     &                (TOTEL, IVersion, THRESH, NDENSE, 
     &                 N, IWLEN, PE, PFREE, LEN, IW, NV, 
     &                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W) 
      INTEGER TOTEL
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N), 
     &     ELEN(N), NCMPA, DEGREE(N), 
     &     LAST(TOTEL), HEAD(TOTEL), NEXT(N), 
     &     W(N)
      INTEGER NDENSE(N)
      INTEGER IVersion, THRESH
      INTEGER THRESM, MINDEN, MAXDEN, NDME
      INTEGER NBD,NBED, NBDM, LASTD, NELME
      LOGICAL IDENSE
      DOUBLE PRECISION RELDEN
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER MAXINT_N
      INTEGER(8) HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod
      LOGICAL COMPRESS
      IF (THRESH.GT.0) THEN 
         THRESM  = min(N,THRESH)
         DO I=1,N
             THRESM = max(THRESM, LEN(I))
          ENDDO
           RELDEN = dble(PFREE-1)/dble(N)
         THRESM =  int(RELDEN)*10 + (THRESM-int(RELDEN))/10 + 1
      ELSE
         THRESM = TOTEL
      ENDIF
      IF (THRESM.GE.0) THEN
       IF ((THRESM.GT.TOTEL).OR.(THRESM.LT.2)) THEN 
          THRESM = TOTEL
       ENDIF
      ENDIF
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        NDENSE(I)= 0
        LAST (I) = 0
        HEAD (I) = 0
        W (I) = 1
        ELEN (I) = 0
   10 CONTINUE
      HEAD(N:TOTEL) = 0
      LAST(N:TOTEL) = 0
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF(COMPRESS) THEN
         DO I=1,N
            DEGREE(I) = 0
            DO J= PE(I) , PE(I)+LEN(I)-1 
               DEGREE(I) = DEGREE(I) + NV(IW(J))
            ENDDO
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            DEGREE (I) = LEN (I)
         ENDDO
      ENDIF
      DO 20 I = 1, N
         DEG = DEGREE (I)
         IF (DEG .GT. 0) THEN
            IF ( (THRESM.GE.0) .AND.
     &           (DEG+NV(I).GE.THRESM) ) THEN
               NBD = NBD+1
               IF (DEG+NV(I).NE.TOTEL-NEL) THEN
                  DEGREE(I) = DEGREE(I)+TOTEL+1
                  DEG = TOTEL
                  INEXT = HEAD (DEG)
                  IF (INEXT .NE. 0) LAST (INEXT) = I
                  NEXT (I) = INEXT
                  HEAD (DEG) = I 
                  LAST(I)  = 0
                  IF (LASTD.EQ.0) LASTD=I
               ELSE
                  NBED = NBED+1
                  DEGREE(I) = TOTEL+1
                  DEG = TOTEL
                  IF (LASTD.EQ.0) THEN
                     LASTD     = I 
                     HEAD(DEG) = I
                     NEXT(I)   = 0 
                     LAST(I)   = 0
                  ELSE
                     NEXT(LASTD) = I
                     LAST(I)     = LASTD
                     LASTD       = I
                     NEXT(I)     = 0
                  ENDIF
               ENDIF
            ELSE
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
            ENDIF
         ELSE
            NEL = NEL + NV(I)
            ELEN (I) = -NEL
            PE (I) = 0
            W (I) = 0
         ENDIF
 20   CONTINUE
          IF (NBD.EQ.0) THRESM = TOTEL
          NLEFT = TOTEL - NEL
 30       IF (NEL .LT. TOTEL) THEN
        DO 40 DEG = MINDEG, TOTEL
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (DEG.LT.TOTEL)  THEN
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          NBDM = max(NBDM,NBD)
          IF (DEGREE(ME).GT.TOTEL+1) THEN
            MINDEN = NBD
            MAXDEN = 0
            IF (WFLG .GT. MAXINT_N) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + LEN(ME) -1
            LN       = P1
            ELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0) THEN
                X = E
  53            X = -PE(X)
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(LN) = IW(ELN)
               IW(ELN) = E
               LN  = LN+1
               ELN = ELN + 1
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(LN)=E
               LN = LN+1
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = LN-P1
            ELEN(ME) = ELN- P1
            NDME = NDENSE(ME)+NV(ME)
            MINDEN = min (MINDEN, NDME)
            MAXDEN = max (MAXDEN, NDME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            IF (IVersion.EQ.1) THEN
              DEG = max (DEGREE(ME)-(TOTEL+1), 1)
            ELSE
              DEG = NDENSE(ME)
            ENDIF
            DEGREE(ME) = DEG
            MINDEG = min(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(TOTEL+1) ) GOTO 51
            ENDIF
            HEAD (TOTEL) = ME
            IF (IVersion .EQ.1 ) THEN
             THRESM = TOTEL
            ELSE
             THRESM=max(THRESM*2,MINDEN+(MAXDEN-MINDEN)/2)
             THRESM = min(THRESM,NBD)
             IF (THRESM.GE.NBD) THRESM=TOTEL
            ENDIF
            NBD    = NBED
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.TOTEL+1) THEN
           IF (NBD.NE.NBED) THEN
            write(6,*) ' Internal ERROR quasi dense rows remains'
            CALL MUMPS_ABORT()
           ENDIF
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.TOTEL+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
   59      CONTINUE
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0
           IF (NEL.NE.TOTEL) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL MUMPS_ABORT()
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.TOTEL) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                 NEXT (ILAST) = INEXT
              ELSE
                 HEAD (DEGREE (I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                IF (DEGREE(I).LE.TOTEL) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                   NEXT (ILAST) = INEXT
                ELSE
                   HEAD (DEGREE (I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#if defined (NOAGG5)
            ELSE IF (DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
            ELSE IF ((DEXT .EQ. 0) .AND.
     &                (NDENSE(ME).EQ.NBD)) THEN
                PE (E) = -ME
                W (E)  = 0
            ELSE IF (DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + int(E,kind=8)
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              IF (DEGREE(J).LE.TOTEL) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
#if defined (NOAGG5)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF ((DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD)) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE(I) = min (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.TOTEL) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.TOTEL) THEN
            DEG = min (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
       IF ( (Iversion .NE. 1).AND. (THRESM.GE.0)) THEN
          IF (DEG+NVI .GE. THRESM) THEN
             IF (THRESM.EQ.TOTEL) THEN
                IF ((ELEN(I).LE.2) .AND. ((DEG+NVI).EQ.NLEFT) ) THEN
                   DEGREE(I) = TOTEL+1
                   IDENSE = .TRUE.
                ENDIF
             ELSE
                IDENSE = .TRUE.
                IF ((ELEN(I).LE.2).AND.((DEG+NVI).EQ.NLEFT) ) THEN
                   DEGREE(I) = TOTEL+1
                ELSE
                   DEGREE(I) = TOTEL+1+DEGREE(I)
                ENDIF
             ENDIF
          ENDIF
          IF (IDENSE) THEN
             P1 = PE(I)
             P2 = P1 + ELEN(I) - 1
             IF (P2.GE.P1) THEN
                DO 264 PJ=P1,P2
                   E= IW(PJ)
                   NDENSE (E) = NDENSE(E) + NVI
 264            CONTINUE
             ENDIF
             NBD = NBD+NVI
             DEG = TOTEL
             IF (DEGREE(I).EQ.TOTEL+1) THEN
                NBED = NBED +NVI
                IF (LASTD.EQ.0) THEN
                   LASTD     = I
                   HEAD(DEG) = I
                   NEXT(I)   = 0
                   LAST(I)   = 0
                ELSE
                   NEXT(LASTD) = I
                   LAST(I)     = LASTD
                   LASTD       = I
                   NEXT(I)     = 0
                ENDIF
             ELSE
                INEXT = HEAD(DEG)
                IF (INEXT .NE. 0) LAST (INEXT) = I
                NEXT (I) = INEXT
                HEAD (DEG) = I
                LAST(I)    = 0
                IF (LASTD.EQ.0) LASTD=I
             ENDIF
          ENDIF
       ENDIF
       IF (.NOT.IDENSE) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          LAST (I) = 0
          HEAD (DEG) = I
       ENDIF
       MINDEG = min (MINDEG, DEG)
            ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
  265 CONTINUE
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      IF(COMPRESS) THEN
        LAST(1:N) = 0
        DEGREE(1:TOTEL-N)=0
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            DEGREE(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (DEGREE(K-N) .NE. 0) THEN
            LAST(I)=DEGREE(K-N)
            ELEN(DEGREE(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
         DO 300 I = 1, N
            K = abs (ELEN (I))
            LAST (K) = I
            ELEN (I) = K
 300     CONTINUE
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_421
      SUBROUTINE MUMPS_560(N, NBBUCK, 
     &     IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &     LAST, NCMPA, DEGREE, WF, NEXT, W, HEAD,
     &     CONSTRAINT,THESON)
      IMPLICIT NONE
      INTEGER N, IWLEN, PFREE, LEN(N), 
     &        ELEN(N), LAST(N), NCMPA, DEGREE(N), NEXT(N),
     &        W(N) 
      INTEGER PE(N), IW(IWLEN), NV(N)
      INTEGER NBBUCK
      INTEGER HEAD(0:NBBUCK+1), WF(N)
      INTEGER CONSTRAINT(N),THESON(N)
      INTEGER PREV,TOTO
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        NBFLAG, NREAL, LASTD, NELME, WF3, WF4, N2, PAS
       INTEGER MAXINT_N
       INTEGER(8) HASH, HMOD
       DOUBLE PRECISION    RMF, RMF1 
       DOUBLE PRECISION    dummy
       INTEGER idummy
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod, huge
      INTEGER TOTEL
      idummy = huge(idummy) - 1
      dummy = dble(idummy)
      N2 = -NBBUCK-1
      PAS = max((N/8), 1)
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, NBBUCK-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
      NBFLAG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DO 10 I = 1, N
         THESON(I) = 0
         LAST (I) = 0
         W (I) = 1
         ELEN (I) = 0
   10 CONTINUE
      TOTEL = 0
      DO I=1,N
         IF (LEN(I).LT.0) THEN 
            DEGREE (I) = N2
            NBFLAG     = NBFLAG +1
            IF (LEN(I).EQ.-N-1) THEN
               LEN (I)    = 0
               PE (I)     = 0
            ELSE
               LEN (I)    = - LEN(I)
            ENDIF
         ELSE
            TOTEL = TOTEL + NV(I)
            DEGREE(I) = 0
            DO J= PE(I) , PE(I)+LEN(I)-1 
               DEGREE(I) = DEGREE(I) + NV(IW(J))
            ENDDO
         ENDIF
      ENDDO
      NREAL = N - NBFLAG
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N2) THEN
             DEG = NBBUCK + 1
             IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
         GOTO 20
        ENDIF
        IF (DEG .GT. 0) THEN
           WF(I) = DEG
           IF (DEG.GT.N) THEN
            DEG = min(((DEG-N)/PAS) + N , NBBUCK)
           ENDIF
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
           NEL = NEL + NV(I)
           ELEN (I) = -NEL
           PE (I) = 0
           W (I) = 0
        ENDIF
   20 CONTINUE
      NLEFT = TOTEL-NEL
 30   IF (NEL .LT. TOTEL) THEN
         DO 40 DEG = MINDEG, NBBUCK
            ME = HEAD (DEG)
            IF (ME .GT. 0) GO TO 50
 40      CONTINUE
 50      MINDEG = DEG
         IF (ME.LE.0) THEN
            NCMPA = -N
            CALL MUMPS_ABORT()
         ENDIF
         IF (DEG.GT.N) THEN
            J = NEXT(ME)
            K = WF(ME)
            IF(CONSTRAINT(ME) .LT. 0) THEN
               K = -1
            ENDIF
 55         CONTINUE
            IF (J.GT.0) THEN
               IF(CONSTRAINT(J) .GE. 0) THEN
                  IF (WF(J).LT.K .OR. K .LT. 0) THEN
                     ME = J
                     K  = WF(ME)
                  ENDIF
               ENDIF
               J= NEXT(J)
               GOTO 55
            ENDIF
            ILAST = LAST(ME)
            INEXT = NEXT(ME)
            IF (INEXT .NE. 0) LAST (INEXT) = ILAST
            IF (ILAST .NE. 0) THEN
               NEXT (ILAST) = INEXT
            ELSE
               HEAD (DEG) = INEXT
            ENDIF
         ELSE
            IF(CONSTRAINT(ME) .GE. 0) GOTO 59
 56         CONTINUE
            IF(NEXT(ME) .NE. 0) THEN
               ME = NEXT(ME)
               IF(CONSTRAINT(ME) .GE. 0) THEN
                  GOTO 59
               ELSE
                  GOTO 56
               ENDIF
            ELSE
 57            DEG = DEG+1
               ME = HEAD(DEG) 
               IF(ME .GT. 0) THEN
                  IF(CONSTRAINT(ME) .GE. 0) THEN
                     GOTO 59
                  ELSE
                     GOTO 56
                  ENDIF
               ELSE
                  GOTO 57
               ENDIF
            ENDIF   
 59         PREV = LAST (ME)
            INEXT = NEXT (ME)
            IF(PREV .NE. 0) THEN
               NEXT(PREV) = INEXT
            ELSE
               HEAD (DEG) = INEXT
            ENDIF
            IF (INEXT .NE. 0) LAST (INEXT) = PREV
         ENDIF
         TOTO = ME
 5910    IF(TOTO .NE. 0) THEN
            J = CONSTRAINT(TOTO)
            IF(J .GT. 0) THEN
               CONSTRAINT(J) = 0
            ENDIF
            TOTO = THESON(TOTO)
            GOTO 5910
         ENDIF
            ELENME = ELEN (ME)
            ELEN (ME) = - (NEL + 1)
            NVPIV = NV (ME)
            NEL = NEL + NVPIV
            NV (ME) = -NVPIV
            DEGME = 0
            IF (ELENME .EQ. 0) THEN
               PME1 = PE (ME)
               PME2 = PME1 - 1
               DO 60 P = PME1, PME1 + LEN (ME) - 1
                  I = IW (P)
                  NVI = NV (I)
                  IF (NVI .GT. 0) THEN
                     DEGME = DEGME + NVI
                     NV (I) = -NVI
                     PME2 = PME2 + 1
                     IW (PME2) = I
                     IF (DEGREE(I).NE.N2) THEN
                        ILAST = LAST (I)
                        INEXT = NEXT (I)
                        IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                        IF (ILAST .NE. 0) THEN
                           NEXT (ILAST) = INEXT
                        ELSE
                           IF (WF(I).GT.N) THEN
                              DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                           ELSE
                              DEG = WF(I)
                           ENDIF
                           HEAD (DEG) = INEXT
                        ENDIF
                     ENDIF
                  ENDIF
 60            CONTINUE
               NEWMEM = 0
            ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  IF (WF(I).GT.N) THEN
                   DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
                  HEAD (DEG) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
                WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG  = 0
          WF3  = 0
          WF4  = 0
          NVI  = -NV(I)
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              IF ( WF(E) .EQ. 0 ) THEN
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG4)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (DEGREE(I).EQ.N2) DEG = N2
#if defined (NOAGG4)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            IF (DEGREE(I).NE.N2) THEN
                 IF ( DEGREE (I).LT.DEG ) THEN
                   WF4 = 0
                   WF3 = 0
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            ENDIF
            WF(I)      = WF4 + 2*NVI*WF3
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            IF (DEG.NE.N2) THEN
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                 IF(CONSTRAINT(J) .LT. 0
     &                .AND. CONSTRAINT(I) .LT. 0) THEN
                    GOTO 240
                 ENDIF
                 IF(CONSTRAINT(I) .GE. 0) THEN
                    IF(CONSTRAINT(J) .LT. 0) THEN
                       TOTO = I
 221                   IF(TOTO .NE. 0) THEN
                          IF(CONSTRAINT(TOTO) .EQ. J) THEN
                             GOTO 225
                          ENDIF
                          TOTO =THESON(TOTO)
                          GOTO 221
                       ENDIF
                    ELSE
                       GOTO 225
                    ENDIF
                 ELSE
                    IF(CONSTRAINT(J) .GE. 0) THEN
                       TOTO = J
 222                   IF(TOTO .NE. 0) THEN
                          IF(CONSTRAINT(TOTO) .EQ. I) THEN
                             GOTO 225
                          ENDIF
                          TOTO =THESON(TOTO)
                          GOTO 222
                       ENDIF
                    ENDIF
                 ENDIF
                 GOTO 240
 225             CONTINUE
                 IF (LEN (J) .NE. LN) GO TO 240
                 IF (ELEN (J) .NE. ELN) GO TO 240
                 DO 230 P = PE (J) + 1, PE (J) + LN - 1
                    IF (W (IW (P)) .NE. WFLG) GO TO 240
 230             CONTINUE
                 TOTO = I
 231             IF(THESON(TOTO) .NE. 0) THEN
                    TOTO = THESON(TOTO)
                    GOTO 231
                 ENDIF
                 THESON(TOTO) = J
                 IF(CONSTRAINT(I) .LT. 0) THEN
                    CONSTRAINT(I) = 0
                 ENDIF
                 PE (J) = -I
                 WF(I)  = max(WF(I),WF(J))
                 NV (I) = NV (I) + NV (J)
                 NV (J) = 0
                 ELEN (J) = 0
                 J = NEXT (J)
                 NEXT (JLAST) = J
                 GO TO 220
 240             CONTINUE
                 JLAST = J
                 J = NEXT (J)
                 GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
           ENDIF
          ENDIF
 250   CONTINUE
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
           I = IW (PME)
           NVI = -NV (I)
           IF (NVI .GT. 0) THEN
              NV (I) = NVI
              IF (DEGREE(I).NE.N2) THEN
                 DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
                 IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
                  DEG = DEGREE(I)
                  RMF1  = dble(DEG)*dble( (DEG-1) + 2*DEGME )
     &                 - dble(WF(I))
                  DEGREE(I) = NLEFT - NVI
                  DEG       = DEGREE(I) 
                  RMF = dble(DEG)*dble(DEG-1) 
     &                 -  dble(DEGME-NVI)*dble(DEGME-NVI-1)
                  RMF = min(RMF, RMF1)
               ELSE 
                  DEG = DEGREE(I)
                  DEGREE(I) = DEGREE (I) + DEGME - NVI
                  RMF  = dble(DEG)*dble( (DEG-1) + 2*DEGME ) 
     &                 - dble(WF(I))
               ENDIF
               RMF =  RMF / dble(NVI+1)
               IF (RMF.LT.dummy) THEN
                  WF(I) = int ( anint( RMF ))
               ELSEIF (RMF / dble(N) .LT. dummy) THEN 
                  WF(I) = int ( anint( RMF/dble(N) ))
               ELSE
                  WF(I) = idummy
               ENDIF
               WF(I) = max(1,WF(I))
               DEG = WF(I)
               IF (DEG.GT.N) THEN
                  DEG = min(((DEG-N)/PAS) + N , NBBUCK)
               ENDIF
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               LAST (I) = 0
               HEAD (DEG) = I
               MINDEG = min (MINDEG, DEG)
            ENDIF
            IW (P) = I
            P = P + 1
         ENDIF
 260  CONTINUE
      NV (ME) = NVPIV + DEGME
      LEN (ME) = P - PME1
      IF (LEN (ME) .EQ. 0) THEN
         PE (ME) = 0
         W (ME) = 0
      ENDIF
      IF (NEWMEM .NE. 0) THEN
         PFREE = P
         MEM = MEM - NEWMEM + LEN (ME)
      ENDIF
      GO TO 30
      ENDIF
      IF (NEL.LT.N) THEN 
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N2) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
           ENDDO
           ELEN(ME) = NELME
           NV(ME)   = N-NREAL
           PE(ME)   = 0
        IF (NEL.NE.N) THEN
         NCMPA = -N - 1
         GOTO 500
        ENDIF
      ENDIF
      DO 290 I = 1, N
         IF (ELEN (I) .EQ. 0) THEN
            J = -PE (I)
 270        CONTINUE
            IF (ELEN (J) .GE. 0) THEN
               J = -PE (J)
               GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
 280        CONTINUE
            IF (ELEN (J) .GE. 0) THEN
               JNEXT = -PE (J)
               PE (J) = -E
               IF (ELEN (J) .EQ. 0) THEN
                  ELEN (J) = K
                  K = K + 1
               ENDIF
               J = JNEXT
               GO TO 280
            ENDIF
            ELEN (E) = -K
         ENDIF
 290  CONTINUE
      IF(.TRUE.) THEN
        LAST(1:N) = 0
        DEGREE(1:TOTEL-N)=0
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            DEGREE(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (DEGREE(K-N) .NE. 0) THEN
            LAST(I)=DEGREE(K-N)
            ELEN(DEGREE(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
        DO 300 I = 1, N
           K = abs (ELEN (I))
           LAST (K) = I
           ELEN (I) = K
300     CONTINUE
      ENDIF
 500  PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_560
      SUBROUTINE MUMPS_422 
     &                ( THRESH, NDENSE, 
     &                 N, IWLEN, PE, PFREE, LEN, IW, NV, 
     &                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     &                 PERM, LISTVAR_SCHUR, SIZE_SCHUR, AGG6 ) 
      IMPLICIT NONE
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N), 
     &        ELEN(N), LAST(N), NCMPA, DEGREE(N), HEAD(N), NEXT(N), 
     &        W(N), SIZE_SCHUR
      LOGICAL AGG6
      INTEGER NDENSE(N), LISTVAR_SCHUR(max(1,SIZE_SCHUR))
      INTEGER PERM(N)
      INTEGER THRESH
      INTEGER THRESM, NDME, PERMeqN
      INTEGER NBD,NBED, NBDM, LASTD, NELME
      LOGICAL IDENSE
      INTEGER FDEG, ThresMin, ThresPrev, IBEGSchur, NbSchur, 
     &        ThresMinINIT
      LOGICAL SchurON
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER MAXINT_N
      INTEGER(8) HASH, HMOD 
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod
        IF (N.EQ.1) THEN
           ELEN(1) = 1
           LAST(1) = 1
           PE(1) = 0
           NV(1) = 1
           NCMPA = 0
           RETURN
        ENDIF
        SIZE_SCHUR = min(N,SIZE_SCHUR)
        SIZE_SCHUR = max(0,SIZE_SCHUR)
        SchurON   = (SIZE_SCHUR > 0)
        IBEGSchur = N-SIZE_SCHUR+1
        IF (THRESH.GT.N) THRESH = N
        IF (THRESH.LT.0) THRESH = 0
        IF ( SchurON )  THEN 
           DO I= 1, N
             IF ( PERM(I) .GE. IBEGSchur) THEN 
                 PERM(I) = N + 1
                IF (LEN(I) .EQ.0) THEN
                  PE(I) = 0
                ENDIF
             ENDIF
           ENDDO
        ENDIF
        IF (SchurON) THEN
             THRESM    = N
             ThresMin  = N
             ThresPrev = N
        ELSE
             THRESM    = max(int(31*N/32),THRESH)
             THRESM    = max(THRESM,1)
             ThresMin  = max( 3*THRESM / 4, 1)
             ThresPrev = THRESM
        ENDIF
        ThresMinINIT = ThresMin/4
      IF (THRESM.GT.0) THEN
       IF ((THRESM.GT.N).OR.(THRESM.LT.2)) THEN 
          THRESM = N
       ENDIF
      ENDIF
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        NDENSE(I)= 0
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE (I) = LEN (I)
   10 CONTINUE
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (PERM(I).EQ.N) THEN
           PERMeqN = I
           PERM(I) = N-1
        ENDIF
        FDEG = PERM(I)
        IF ( (DEG .GT. 0).OR.(PERM(I).EQ.N+1) ) THEN
          IF ( (THRESM.GT.0) .AND.
     &         (FDEG .GT.THRESM) ) THEN
            NBD = NBD+1
            IF (FDEG.NE.N+1) THEN
             DEGREE(I) = DEGREE(I)+N+2
             DEG = N
             INEXT = HEAD (DEG)
             IF (INEXT .NE. 0) LAST (INEXT) = I
             NEXT (I) = INEXT
             HEAD (DEG) = I 
             LAST(I)  = 0
             IF (LASTD.EQ.0) LASTD=I
            ELSE
             NBED = NBED+1
             DEGREE(I) = N+1
             DEG = N
             IF (LASTD.EQ.0) THEN
               LASTD     = I 
               HEAD(DEG) = I
               NEXT(I)   = 0 
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
            ENDIF
          ELSE
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            HEAD (FDEG) = I
          ENDIF
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
          IF ((NBD.EQ.0).AND.(THRESM.GT.0)) THRESM = N
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF ( (DEG.NE.N) .AND.
     &    (DEG.GT.THRESM+1) .AND. (NBD.GT.0) ) THEN
           MINDEG = N
           GOTO 30
        ENDIF
        IF (DEGREE(ME).LE.N)  THEN
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          MINDEG = 1
          NBDM = max(NBDM,NBD)
          IF (DEGREE(ME).GT.N+1) THEN
            IF (WFLG .GT. MAXINT_N) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + LEN(ME) -1
            LN       = P1
            ELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0) THEN
                X = E
  53            X = -PE(X)
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(LN) = IW(ELN)
               IW(ELN) = E
               LN  = LN+1
               ELN = ELN + 1
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(LN)=E
               LN = LN+1
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = LN-P1
            ELEN(ME) = ELN- P1
            NDME = NDENSE(ME)+NV(ME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            DEGREE(ME) = NDENSE(ME)
            DEG = PERM(ME)
            MINDEG = min(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(N+1) ) GOTO 51
            ENDIF
            HEAD (N) = ME
            IF (THRESM.LT.N) THEN
             ThresMin  = max(THRESM+ThresMin,ThresPrev+ThresMin/2+1)
             ThresMin  = min(ThresMin, N)
             ThresPrev = ThresPrev+(N-ThresPrev)/2+ThresMinINIT
             THRESM    = max(
     &         THRESM + int(sqrt(dble(ThresMin)))+ ThresMinINIT ,
     &         ThresPrev)
             THRESM    = min(THRESM,N) 
             ThresMin  = min(THRESM, ThresMin)
             ThresPrev = THRESM
            ENDIF
            NBD    = NBED
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.N+1) THEN
           IF (NBD.NE.NBED) THEN
          write(6,*) ' ERROR in MUMPS_422 quasi dense rows remains'
            CALL MUMPS_ABORT()
           ENDIF
           NbSchur = 0   
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
              PE(X) = -LISTVAR_SCHUR(1)
            ELSE IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -LISTVAR_SCHUR(1)
            ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
             NbSchur = NbSchur+ 1
            ENDIF
   59      CONTINUE
           IF (NbSchur.NE.SIZE_SCHUR) then
             write(6,*) ' Internal error 2 in QAMD :',
     &         ' Schur size expected:',SIZE_SCHUR, 'Real:', NbSchur
             CALL MUMPS_ABORT()
           ENDIF
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0
           IF (NEL.NE.N) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL MUMPS_ABORT()
           ENDIF
           IF (ME.NE. LISTVAR_SCHUR(1)) THEN
             DO I=1, SIZE_SCHUR
               PE(LISTVAR_SCHUR(I)) = -LISTVAR_SCHUR(1)
             ENDDO
             PE(LISTVAR_SCHUR(1)) = 0
             NV( LISTVAR_SCHUR(1))= NV(ME)
             NV(ME)               = 0
             ELEN( LISTVAR_SCHUR(1)) = ELEN(ME)
             ELEN(ME)             = 0
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.N) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (PERM(I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (PERM(I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (.NOT. AGG6 .AND. DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (AGG6 .AND. (DEXT .EQ. 0) .AND.
     &            ((NDENSE(ME).EQ.NBD).OR.(NDENSE(E).EQ.0))) THEN
                PE (E) = -ME
                W (E)  = 0
             ELSE IF (AGG6 .AND. DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + int(E,kind=8)
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              IF (DEGREE(J).LE.N) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (((ELEN(I).EQ.1).AND.(P3.EQ.PN))
     &     .OR.
     &         (AGG6.AND.(DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD))
     &       )
     &    THEN
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE(I) = min (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.N) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
             X = I 
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                IF (PERM(J).GT.PERM(X)) THEN
                  PE (J) = -X
                  NV (X) = NV (X) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                ELSE
                  PE (X) = -J
                  NV (J) = NV (X) + NV (J)
                  NV (X) = 0
                  ELEN (X) = 0
                  X = J
                ENDIF
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        IF ( (THRESM .GT. 0).AND.(THRESM.LT.N) ) THEN 
          THRESM = max(ThresMin, THRESM-NVPIV)
        ENDIF
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.N) THEN
            DEG = min (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
            IF (THRESM.GT.0) THEN
             IF (PERM(I) .GT. THRESM) THEN
               IDENSE = .TRUE.
               DEGREE(I) = DEGREE(I)+N+2
             ENDIF
             IF (IDENSE) THEN
               P1 = PE(I)
               P2 = P1 + ELEN(I) - 1
               IF (P2.GE.P1) THEN
               DO 264 PJ=P1,P2
                 E= IW(PJ)
                 NDENSE (E) = NDENSE(E) + NVI
 264           CONTINUE
               ENDIF
               NBD = NBD+NVI
               FDEG = N
               DEG = N
               INEXT = HEAD(DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
               LAST(I)    = 0
               IF (LASTD.EQ.0) LASTD=I
             ENDIF
            ENDIF
            IF (.NOT.IDENSE) THEN
            FDEG = PERM(I)
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (FDEG) = I
            ENDIF
            MINDEG = min (MINDEG, FDEG)
            ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
  265 CONTINUE
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      IF (.NOT.SchurON) THEN
        PERM(PERMeqN) = N
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_422 
      SUBROUTINE MUMPS_276( ICNTL, INFO, COMM, ID )
      INTEGER ICNTL(40), INFO(40), COMM, ID
      INCLUDE 'mpif.h'
      INTEGER IN( 2 ), OUT( 2 )
      INTEGER LP, IERR
      LP      = ICNTL( 1 )
      IN( 1 ) = INFO ( 1 )
      IN( 2 ) = ID
      CALL MPI_ALLREDUCE( IN, OUT, 1, MPI_2INTEGER, MPI_MINLOC,
     &                    COMM, IERR)
      IF ( OUT( 1 ) .LT. 0 .and. INFO(1) .GE. 0 ) THEN
        INFO( 1 ) = -001
        INFO( 2 ) = OUT( 2 )
      END IF
      RETURN
      END SUBROUTINE MUMPS_276
      SUBROUTINE MUMPS_137( INODE, N, PROCNODE_STEPS,
     &           SLAVEF,
     &           ND, FILS, FRERE_STEPS, STEP, PIMASTER,
     &           KEEP28, KEEP50, KEEP253,
     &           FLOP1,
     &           IW, LIW, XSIZE )
      IMPLICIT NONE
      INTEGER INODE, N, KEEP50, LIW, SLAVEF, KEEP28, KEEP253
      INTEGER PROCNODE_STEPS(KEEP28), ND(KEEP28),
     &        FILS(N), FRERE_STEPS(KEEP28),
     &        STEP(N), 
     & PIMASTER(KEEP28),
     &  IW( LIW )
      INTEGER XSIZE
      DOUBLE PRECISION FLOP1
      INTEGER NUMORG, IN, NASS, IFSON, NUMSTK, NFRONT, NPIV, NCB, 
     &        LEVEL, ISON
      LOGICAL MUMPS_170
      INTEGER MUMPS_330
      EXTERNAL MUMPS_170, MUMPS_330
      INCLUDE 'mumps_headers.h'
      FLOP1 = 0.0D0
      IF (MUMPS_170(PROCNODE_STEPS(STEP(INODE)),
     &                SLAVEF) ) RETURN
      IN     = INODE
      NUMORG = 0
   10 NUMORG = NUMORG + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 10
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .EQ. 0) GOTO 30
   20 NUMSTK = NUMSTK + 1
      NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 +XSIZE)
      ISON = FRERE_STEPS(STEP(ISON))
      IF (ISON .GT. 0) GOTO 20
   30 NFRONT = ND(STEP(INODE)) + NASS + KEEP253
      NPIV  = NASS + NUMORG
      NCB   = NFRONT - NPIV
      LEVEL = MUMPS_330(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      CALL MUMPS_511(NFRONT,NPIV,NPIV,KEEP50,LEVEL,FLOP1)
      RETURN
      END SUBROUTINE MUMPS_137
      SUBROUTINE MUMPS_511(NFRONT,NPIV,NASS,
     &                                 KEEP50,LEVEL,COST)
      IMPLICIT NONE
      INTEGER, intent(in) :: NFRONT,NPIV,KEEP50,LEVEL, NASS
      DOUBLE PRECISION, intent(out) :: COST
      IF (KEEP50.EQ.0) THEN
        IF (LEVEL.EQ.1 .OR. LEVEL.EQ.3) THEN
          COST = dble(2) * dble(NFRONT) * dble(NPIV) *
     &      dble(NFRONT - NPIV - 1) +
     &      dble(NPIV) * dble(NPIV + 1) * dble(2 * NPIV + 1)
     &          / dble(3)
          COST = COST + dble(2 * NFRONT - NPIV - 1)
     &      * dble(NPIV) /dble(2)
        ELSEIF (LEVEL.EQ.2) THEN
          COST = dble(2*NASS)*dble(NFRONT) -
     &          dble(NASS+NFRONT)*dble(NPIV+1)
          COST = dble(NPIV)*COST +
     &     dble(2 * NASS - NPIV - 1) * dble(NPIV) / dble(2) +
     &     dble(NPIV) * dble(NPIV + 1) *
     &     dble(2 * NPIV + 1) /dble(3)
        ENDIF
      ELSE
        IF (LEVEL.EQ.1) THEN
          COST = dble(NPIV) * (
     &          dble( NFRONT ) * dble( NFRONT ) +
     &          dble( NFRONT ) - (
     &          dble( NFRONT)*dble(NPIV) + dble(NPIV+1)
     &          )) +( dble(NPIV)*dble(NPIV+1)
     &          *dble(2*NPIV+1))/ dble(6)
        ELSE IF (LEVEL.EQ.3.AND.KEEP50.EQ.2) THEN
          COST = dble(2) * dble(NFRONT) * dble(NPIV) *
     &      dble(NFRONT - NPIV - 1) +
     &      dble(NPIV) * dble(NPIV + 1) *
     &      dble(2 * NPIV + 1) / dble(3)
          COST = COST + dble(2 * NFRONT - NPIV - 1)
     &         * dble(NPIV) / dble(2)
        ELSE
          COST = dble(NPIV) * (
     &          dble( NASS ) * dble( NASS ) + dble( NASS )
     &        - ( dble( NASS) * dble(NPIV) + dble( NPIV + 1 ) ) )
     &        + ( dble(NPIV)*dble(NPIV+1)*dble(2*NPIV+1) )
     &        / dble( 6 )
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_511
      SUBROUTINE MUMPS_81(MYID, INODE, N, IOLDPS,
     &           HF, NFRONT, NFRONT_EFF, DAD,
     &           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     &           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     &           INTARR, ITLOC, RHS_MUMPS, FILS, FRERE, 
     &           SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG,
     &           PROCNODE_STEPS, SLAVEF )
      IMPLICIT NONE
      INTEGER, intent(in) :: INODE, N, IOLDPS, HF, NFRONT, 
     &        NASS1, LIW, NASS,
     &        NUMSTK, NUMORG, IWPOSCB 
      INTEGER, intent(in) :: KEEP(500)
      INTEGER(8) , intent(in) ::KEEP8(150)
      INTEGER STEP(N), 
     &        PIMASTER(KEEP(28)),
     &        PTRAIW(N), IW(LIW), 
     &        ITLOC(N+KEEP(253)), FILS(N), FRERE(KEEP(28))
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER, intent(inout) ::  NBPROCFILS(KEEP(28))
      LOGICAL, intent(in)    :: NIV1
      INTEGER, intent(inout) :: IFLAG
      LOGICAL, intent(out)   :: SON_LEVEL2
      INTEGER, intent(out)   :: NFRONT_EFF
      INTEGER, intent(in) :: PROCNODE_STEPS(KEEP(28)), SLAVEF
      INTEGER, intent(in)    :: DAD (KEEP(28)), IFSON, MYID
      INTEGER NEWEL, INEW, IOLDP2, INEW1,
     &        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     &        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     &        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     &        I, K
      LOGICAL LEVEL1
      INTEGER  MUMPS_810, MUMPS_330
      EXTERNAL MUMPS_810, MUMPS_330
      INTEGER  TYPESPLIT
      INCLUDE 'mumps_headers.h'
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1   
      ICT11  = IOLDP2 + NFRONT
      NTOTFS = 0
      TYPESPLIT  = MUMPS_810 (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
      IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) ) THEN
        J2    = PIMASTER(STEP(IFSON))
        LSTK  = IW(J2    +KEEP(IXSZ))
        NELIM = IW(J2 + 1+KEEP(IXSZ))
        NPIVS  = IW(J2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(J2 + 5+KEEP(IXSZ))
        IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
        LEVEL1    = NSLSON.EQ.0
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF (NIV1) THEN
         write(6,*) MYID, ':',
     &    ' Internal error 2 in MUMPS_BUILD__INDEX ',
     &    ' interior split node of type 1 '
         CALL MUMPS_ABORT()
        ELSE
         I= MUMPS_330(PROCNODE_STEPS(STEP(IFSON)),SLAVEF)
         J= MUMPS_810(PROCNODE_STEPS(STEP(IFSON)), 
     &              SLAVEF)
         IF (LEVEL1.or.J.LT.4) THEN
         write(6,*) MYID, ':',
     &     ' Internal error 3 in MUMPS_81 ',
     &     ' son', IFSON, 
     &     ' of interior split node', INODE, ' of type 1 ', 
     &     ' NSLSON =', NSLSON, ' TYPE_SON=', I, 'TYPESPLIT_SON=', J
          CALL MUMPS_ABORT()
         ELSE
          NBPROCFILS(STEP(IFSON)) = NSLSON
          NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(IFSON))
         ENDIF
        ENDIF
        IF ( J2.GT. IWPOSCB ) THEN
          NROWS = IW(J2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        J1 = J2 + HS + NROWS + NPIVS
        J2 = J1 + LSTK - 1 
        J3 = J1 + NELIM - 1
        IF (NELIM.GT.0) THEN
         DO JJ=J1,J3
          NTOTFS = NTOTFS + 1
          JT1 = IW(JJ)
          IW(ICT11 + NTOTFS) = JT1
          IW(JJ) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
         ENDDO
        ENDIF
        DO JJ =J3+1, J3+NUMORG
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS   
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        DO JJ =J3+NUMORG+1, J2
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS    
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        NFRONT_EFF = NTOTFS
        IBROT = INODE
        DO IORG = 1, NUMORG
          K1 = PTRAIW(IBROT) + 2
          JT1 = INTARR(K1)
          INTARR(K1) = ITLOC(JT1)
          IBROT = FILS(IBROT)
         K2 = K1 + INTARR(K1 - 2) - INTARR(K1 - 1)
         K1 = K1 + 1
         IF (K1 .LE. K2) THEN
          DO JJ = K1, K2
            J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
          ENDDO
         ENDIF
        ENDDO
        K1 = IOLDPS+HF
        DO JJ=K1+NELIM,K1+NFRONT_EFF-1
          ITLOC(IW(JJ)) = 0
        ENDDO
       RETURN
      ENDIF
      NEWEL  = IOLDP2 + NASS1
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF 
      INEW1 = 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW) = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      IF  (TYPESPLIT.EQ.4) THEN
         IBROT = INODE
         DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )  
          IBROT = DAD(STEP(IBROT)) 
          IN = IBROT
          DO WHILE (IN.GT.0)
            NFRONT_EFF = NFRONT_EFF+1
            NEWEL      = NEWEL + 1
            ITLOC(IN)  = NFRONT_EFF
            IW(NEWEL)  = IN
            IN         = FILS( IN )
          ENDDO
         ENDDO
      ENDIF
      IF (NUMSTK .NE. 0) THEN
        NTOTFS = NUMORG
        ISON = IFSON
        DO 100 IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+KEEP(IXSZ))
          NELIM = IW(J2 + 1+KEEP(IXSZ))
          NPIVS = IW(J2 + 3+KEEP(IXSZ))
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON = IW(J2 + 5+KEEP(IXSZ))
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+KEEP(IXSZ))
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1 - KEEP(253)
          J3 = J1 + NELIM - 1
          IF (NELIM .EQ. 0) GOTO 70
          DO 60 JJ = J1, J3
            NTOTFS = NTOTFS + 1
            JT1 = IW(JJ)
            IW(ICT11 + NTOTFS) = JT1
            ITLOC(JT1) = NTOTFS
            IW(JJ) = NTOTFS
            IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
   60     CONTINUE
   70     J1 = J3 + 1
          IF (NASS1 .NE. NFRONT - KEEP(253)) THEN
            DO 80 JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN
                NEWEL = NEWEL + 1
                NFRONT_EFF = NFRONT_EFF + 1
                IW(NEWEL) = J
                IW(JJ) = NFRONT_EFF
                ITLOC(J) = NFRONT_EFF
              ELSE
                IW(JJ) = ITLOC(J)
              ENDIF
   80       CONTINUE
          ELSE
            DO 90 JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
   90       CONTINUE
          ENDIF
          DO JJ=J2+1, J2+KEEP(253)
              IW(JJ)=NFRONT-KEEP(253)+JJ-J2
          ENDDO
          ISON = FRERE(STEP(ISON))
  100   CONTINUE
      ENDIF
      IBROT = INODE
      DO 120 IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO 110 JJ = J1, J2
            J = INTARR(JJ)
            IF (ITLOC(J) .EQ. 0) THEN
              NEWEL = NEWEL + 1
              NFRONT_EFF = NFRONT_EFF + 1
              IW(NEWEL) = J
              INTARR(JJ) = NFRONT_EFF
              ITLOC(J) = NFRONT_EFF
            ELSE
              INTARR(JJ) = ITLOC(J)
            ENDIF
  110     CONTINUE
        ENDIF
  120 CONTINUE
      IF ( (TYPESPLIT.EQ.4).AND.(NFRONT_EFF.LT.NFRONT-KEEP(253)) ) THEN
         IBROT = INODE
         DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )  
          IBROT = DAD(STEP(IBROT)) 
          IN = IBROT
          DO WHILE (IN.GT.0.AND.NFRONT_EFF.LT.NFRONT-KEEP(253))
            J1    = PTRAIW(IN) + 2
            J2    = J1 + INTARR(J1 - 2) - INTARR(J1-1)
            IN = FILS( IN )
            DO JJ = J1+1, J2
              J     = INTARR( JJ )
              IF ( ITLOC( J ) .eq. 0 ) THEN
               NEWEL  = NEWEL + 1
               NFRONT_EFF = NFRONT_EFF + 1
               IW( NEWEL ) = J
               ITLOC( J ) = NFRONT_EFF
              END IF
            ENDDO
          ENDDO
          IF (NFRONT_EFF.EQ.NFRONT-KEEP(253)) EXIT  
         ENDDO
      ENDIF
      IF ( KEEP(253).NE.0) THEN
        IP1 = IOLDPS +  HF + NFRONT_EFF  
        IP2 = IOLDPS + HF + NFRONT + NFRONT_EFF 
        DO I= 1, KEEP(253)
          IW(IP1+I-1) = N+I
          IW(IP2+I-1) = N+I
        ENDDO
        NFRONT_EFF = NFRONT_EFF + KEEP(253)
      ENDIF
      IF (NFRONT.NE.NFRONT_EFF) THEN
       IF (NUMORG.EQ.NASS1) THEN
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NFRONT_EFF - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ELSE
        IP1 = IOLDPS + NFRONT + HF + NUMORG 
        IP2 = IOLDPS + HF + NFRONT_EFF + NUMORG 
        DO I=1,NASS
          IW(IP2+I-1)=IW(IP1+I-1)
        ENDDO
        IP1 = IOLDPS + NASS1 + HF 
        IP2 = IOLDPS + HF + NFRONT - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NUMORG - 1 
        DO I = IP1, IP2
          IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ENDIF
      ELSE
       IP1 = IOLDPS + NASS1 + HF 
       IP2 = IOLDPS + HF + NFRONT - KEEP(253) - 1 
       DO I = IP1, IP2
        IW(I + NFRONT) = IW(I)
       ENDDO
       IP1 = IOLDPS + HF 
       IP2 = IOLDPS + HF + NUMORG - 1 
       DO I = IP1, IP2
         IW(I + NFRONT) = IW(I)
       ENDDO
      ENDIF
      K1 = IOLDPS + HF + NUMORG 
      K2 = K1 + NFRONT_EFF - 1 + NASS
      DO 150 K = K1, K2
        I = IW(K)
        ITLOC(I) = 0
  150 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_81
      SUBROUTINE MUMPS_124(
     &           NUMELT, LIST_ELT,
     &           MYID, INODE, N, IOLDPS,
     &           HF, NFRONT, NFRONT_EFF,
     &           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     &           IFSON, STEP, PIMASTER, PTRAIW, NELT,
     &           IW, LIW, 
     &           INTARR, LINTARR, ITLOC, RHS_MUMPS,
     &           FILS, FRERE_STEPS, 
     &           KEEP,
     &           SON_LEVEL2, NIV1, NBPROCFILS, IFLAG, 
     &           DAD, PROCNODE_STEPS, SLAVEF, 
     &           FRT_PTR, FRT_ELT, Pos_First_NUMORG)
      IMPLICIT NONE
      INTEGER NELT, INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     &        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG,
     &        LINTARR, NUMELT, NFRONT_EFF
      INTEGER KEEP(500)  
      INTEGER LIST_ELT(*)
      INTEGER STEP(N), 
     & PIMASTER(KEEP(28)),
     &  PTRAIW(NELT+1), IW(LIW), 
     &        ITLOC(N+KEEP(253)), FILS(N),
     &        FRERE_STEPS(KEEP(28)),
     &        NBPROCFILS(KEEP(28))
      COMPLEX, POINTER, DIMENSION(:) :: RHS_MUMPS
      INTEGER INTARR(LINTARR)
      LOGICAL SON_LEVEL2, NIV1
      INTEGER, intent(in)    :: DAD (KEEP(28))
      INTEGER, intent(in) :: PROCNODE_STEPS(KEEP(28)), SLAVEF
      INTEGER, intent(in) :: FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER, intent(out) :: Pos_First_NUMORG
      INTEGER NEWEL, INEW, IOLDP2, INEW1,
     &        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     &        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     &        NROWS, HS, IP1, IP2, K1, K2, 
     &        I, K, ELTI
      LOGICAL LEVEL1
      INTEGER  MUMPS_810, MUMPS_330
      EXTERNAL MUMPS_810, MUMPS_330
      INTEGER  TYPESPLIT, NUMELT_IBROT, IBROT
      INCLUDE 'mumps_headers.h'
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1    
      NTOTFS = 0
      ICT11 = IOLDP2 + NFRONT
      TYPESPLIT  = MUMPS_810 (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
      IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) ) THEN
        J2    = PIMASTER(STEP(IFSON))
        LSTK  = IW(J2    +KEEP(IXSZ))
        NELIM = IW(J2 + 1+KEEP(IXSZ))
        NPIVS  = IW(J2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(J2 + 5+KEEP(IXSZ))
        IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
        LEVEL1    = NSLSON.EQ.0
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF (NIV1) THEN
         write(6,*) MYID, ':',
     &    ' Internal error 2 in MUMPS_BUILD__INDEX ',
     &    ' interior split node of type 1 '
         CALL MUMPS_ABORT()
        ELSE
         I= MUMPS_330(PROCNODE_STEPS(STEP(IFSON)),SLAVEF)
         J= MUMPS_810(PROCNODE_STEPS(STEP(IFSON)), 
     &              SLAVEF)
         IF (LEVEL1.or.J.LT.4) THEN
         write(6,*) MYID, ':',
     &     ' Internal error 3 in MUMPS_81 ',
     &     ' son', IFSON, 
     &     ' of interior split node', INODE, ' of type 1 ', 
     &     ' NSLSON =', NSLSON, ' TYPE_SON=', I, 'TYPESPLIT_SON=', J
          CALL MUMPS_ABORT()
         ELSE
          NBPROCFILS(STEP(IFSON)) = NSLSON
          NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(IFSON))
         ENDIF
        ENDIF
        IF ( J2.GT. IWPOSCB ) THEN
          NROWS = IW(J2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        J1 = J2 + HS + NROWS + NPIVS
        J2 = J1 + LSTK - 1
        J3 = J1 + NELIM - 1
        IF (NELIM.GT.0) THEN
         DO JJ=J1,J3
          NTOTFS = NTOTFS + 1
          JT1 = IW(JJ)
          IW(ICT11 + NTOTFS) = JT1
          IW(JJ) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
         ENDDO
        ENDIF
        DO JJ =J3+1, J2
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS    
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        NFRONT_EFF = NTOTFS
        DO IELL=1,NUMELT
         ELTI = LIST_ELT(IELL)
         J1= PTRAIW(ELTI)
         J2= PTRAIW(ELTI+1)-1
         DO JJ=J1,J2
          J = INTARR(JJ)
          INTARR(JJ) = ITLOC(J)
         END DO
        ENDDO
        K1 = IOLDPS+HF
        DO JJ=K1+NELIM,K1+NFRONT_EFF-1
          ITLOC(IW(JJ)) = 0
        ENDDO
        RETURN   
      ENDIF
      NEWEL = IOLDP2 + NASS1
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF 
      INEW1 = 1
      DO WHILE (IN.GT.0)
       ITLOC(IN) = INEW1
       IW(INEW)  = IN
       INEW1     = INEW1 + 1
       INEW      = INEW + 1
       IN = FILS(IN)
      END DO
      IF  (TYPESPLIT.EQ.4) THEN
        IBROT = INODE
         DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )  
          IBROT = DAD(STEP(IBROT)) 
          IN = IBROT
          DO WHILE (IN.GT.0)
            NFRONT_EFF = NFRONT_EFF+1
            NEWEL      = NEWEL + 1
            ITLOC(IN)  = NFRONT_EFF
            IW(NEWEL)  = IN
            IN         = FILS( IN )
          ENDDO
         ENDDO
      ENDIF
      IF (NUMSTK .NE. 0) THEN
        NTOTFS = NUMORG
        ISON = IFSON
        DO 100 IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+KEEP(IXSZ))
          NELIM = IW(J2 + 1+KEEP(IXSZ))
          NPIVS = IW(J2 + 3+KEEP(IXSZ))
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON = IW(J2 + 5+KEEP(IXSZ))
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+KEEP(IXSZ))
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 +KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1 - KEEP(253)
          J3 = J1 + NELIM - 1
          IF (NELIM .EQ. 0) GOTO 70
          DO 60 JJ = J1, J3
            NTOTFS = NTOTFS + 1
            JT1 = IW(JJ)
            IW(ICT11 + NTOTFS) = JT1
            ITLOC(JT1) = NTOTFS
            IW(JJ) = NTOTFS
            IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
   60     CONTINUE
   70     J1 = J3 + 1
          IF (NASS1 .NE. NFRONT) THEN
            DO 80 JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN
                NEWEL = NEWEL + 1
                NFRONT_EFF = NFRONT_EFF + 1
                IW(NEWEL) = J
                IW(JJ) = NFRONT_EFF
                ITLOC(J) = NFRONT_EFF
              ELSE
                IW(JJ) = ITLOC(J)
              ENDIF
   80       CONTINUE
          ELSE
            DO 90 JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
   90       CONTINUE
          ENDIF
          DO JJ=J2+1, J2+KEEP(253)
              IW(JJ)=NFRONT-KEEP(253)+JJ-J2
          ENDDO
          ISON = FRERE_STEPS(STEP(ISON))
  100   CONTINUE
      ENDIF
      DO IELL=1,NUMELT 
       ELTI = LIST_ELT(IELL)
       J1= PTRAIW(ELTI)
       J2= PTRAIW(ELTI+1)-1
       DO JJ=J1,J2
          J = INTARR(JJ)
          IF (ITLOC(J) .EQ. 0) THEN
              NEWEL = NEWEL + 1
              NFRONT_EFF = NFRONT_EFF + 1
              IW(NEWEL)  = J
              INTARR(JJ) = NFRONT_EFF
              ITLOC(J)   = NFRONT_EFF
          ELSE
              INTARR(JJ) = ITLOC(J)
          ENDIF
       END DO
      ENDDO
      IF ( (TYPESPLIT.EQ.4).AND.(NFRONT_EFF.LT.NFRONT-KEEP(253)) ) THEN
         IBROT = INODE
         DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )  
          IBROT = DAD(STEP(IBROT)) 
          NUMELT_IBROT = FRT_PTR(IBROT+1) - FRT_PTR(IBROT)
          IF (NUMELT_IBROT.EQ.0) CYCLE 
          DO IELL = FRT_PTR(IBROT), FRT_PTR(IBROT+1)
            ELTI = FRT_ELT(IELL)
            J1= PTRAIW(ELTI)
            J2= PTRAIW(ELTI+1)-1
            DO JJ=J1,J2
              J     = INTARR( JJ )
              IF ( ITLOC( J ) .eq. 0 ) THEN
               NEWEL  = NEWEL + 1
               NFRONT_EFF = NFRONT_EFF + 1
               IW( NEWEL ) = J
               ITLOC( J ) = NFRONT_EFF
              END IF
            ENDDO
          ENDDO
          IF (NFRONT_EFF.EQ.NFRONT) EXIT  
         ENDDO
      ENDIF  
      IF ( KEEP(253).GT.0) THEN
        IP1 = IOLDPS +  HF + NFRONT_EFF  
        IP2 = IOLDPS + HF + NFRONT + NFRONT_EFF 
        DO I= 1, KEEP(253)
          IW(IP1+I-1) = N+I
          IW(IP2+I-1) = N+I
        ENDDO
        NFRONT_EFF = NFRONT_EFF + KEEP(253)
      ENDIF
      IF (NFRONT.NE.NFRONT_EFF) THEN
       IF (NUMORG.EQ.NASS1) THEN
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NFRONT_EFF - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ELSE
        IP1 = IOLDPS + NFRONT + HF + NUMORG 
        IP2 = IOLDPS + HF + NFRONT_EFF + NUMORG 
        DO I=1,NASS
          IW(IP2+I-1)=IW(IP1+I-1)
        ENDDO
        IP1 = IOLDPS + NASS1 + HF 
        IP2 = IOLDPS + HF + NFRONT - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NUMORG - 1 
        DO I = IP1, IP2
          IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ENDIF
      ELSE
       IP1 = IOLDPS + NASS1 + HF 
       IP2 = IOLDPS + HF + NFRONT - 1 
       DO I = IP1, IP2
        IW(I + NFRONT) = IW(I)
       ENDDO
       IP1 = IOLDPS + HF 
       IP2 = IOLDPS + HF + NUMORG - 1 
       DO I = IP1, IP2
         IW(I + NFRONT) = IW(I)
       ENDDO
      ENDIF
      Pos_First_NUMORG = ITLOC(INODE)
      K1 = IOLDPS + HF + NUMORG 
      K2 = K1 + NFRONT_EFF - 1 + NASS
      DO 150 K = K1, K2
        I = IW(K)
        ITLOC(I) = 0
  150 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_124
      SUBROUTINE MUMPS_86(MYID, INODE, N, IOLDPS,
     &           HF, NFRONT, NFRONT_EFF, PERM, DAD,
     &           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     &           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     &           INTARR, ITLOC, RHS_MUMPS, FILS, FRERE_STEPS, 
     &           SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG,
     &           ISON_IN_PLACE, PROCNODE_STEPS, SLAVEF )
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     &        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID
      INTEGER, intent(in) ::  ISON_IN_PLACE
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N), PIMASTER(KEEP(28)),
     &        PTRAIW(N), IW(LIW), 
     &        ITLOC(N+KEEP(253)), FILS(N), FRERE_STEPS(KEEP(28)),
     &        NBPROCFILS(KEEP(28)), PERM(N)
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      LOGICAL, intent(in)    ::  NIV1
      INTEGER, intent(inout) :: IFLAG
      LOGICAL, intent(out)   :: SON_LEVEL2
      INTEGER, intent(out)   :: NFRONT_EFF
      INTEGER, intent(in)    :: DAD (KEEP(28))
      INTEGER, intent(in) :: PROCNODE_STEPS(KEEP(28)), SLAVEF
      INTEGER NELIM_SON_IN_PLACE 
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     &        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     &        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     &        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     &        I, K, JDEBROW, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     &        LAST_J_ASS, JMIN, MIN_PERM
      LOGICAL LEVEL1
      INTEGER TYPESPLIT
      INCLUDE 'mumps_headers.h'
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      INTEGER  MUMPS_810, MUMPS_330 
      EXTERNAL MUMPS_810, MUMPS_330 
      TYPESPLIT  = MUMPS_810 (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NTOTFS = 0
      NELIM_SON_IN_PLACE = 0
      IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) ) THEN
        J2    = PIMASTER(STEP(IFSON))
        LSTK  = IW(J2    +KEEP(IXSZ))
        NELIM = IW(J2 + 1+KEEP(IXSZ))
        IF ( ISON_IN_PLACE > 0 ) THEN
          IF (ISON_IN_PLACE.NE.IFSON) THEN
         write(6,*) MYID, ':',
     &   ' Internal error 1 in MUMPS_86 ',
     &   ' in place node is not the first son a interior split node '
         CALL MUMPS_ABORT()
          ENDIF
          NELIM_SON_IN_PLACE = NELIM
        ENDIF
        NPIVS  = IW(J2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(J2 + 5+KEEP(IXSZ))
        IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
        LEVEL1    = NSLSON.EQ.0
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF (NIV1) THEN
         write(6,*) MYID, ':',
     &    ' Internal error 2 in MUMPS_86 ',
     &    ' interior split node of type 1 '
         CALL MUMPS_ABORT()
        ELSE
         I= MUMPS_330(PROCNODE_STEPS(STEP(IFSON)),SLAVEF)
         J= MUMPS_810(PROCNODE_STEPS(STEP(IFSON)), 
     &              SLAVEF)
         IF (LEVEL1.or.J.LT.4) THEN
         write(6,*) MYID, ':',
     &     ' Internal error 3 in MUMPS_86 ',
     &     ' son', IFSON, 
     &     ' of interior split node', INODE, ' of type 1 ', 
     &     ' NSLSON =', NSLSON, ' TYPE_SON=', I, 'TYPESPLIT_SON=', J
          CALL MUMPS_ABORT()
         ELSE
          NBPROCFILS(STEP(IFSON)) = NSLSON
          NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(IFSON))
         ENDIF
        ENDIF
        IF ( J2.GT. IWPOSCB ) THEN
          NROWS = IW(J2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        J1 = J2 + HS + NROWS + NPIVS
        J2 = J1 + LSTK - 1
        J3 = J1 + NELIM - 1
        IF (NELIM.GT.0) THEN
         DO JJ=J1,J3
          NTOTFS = NTOTFS + 1
          JT1 = IW(JJ)
          IW(ICT11 + NTOTFS) = JT1
          IW(JJ) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
         ENDDO
        ENDIF
        DO JJ =J3+1, J3+NUMORG  
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS   
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        DO JJ =J3+NUMORG+1, J2
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS    
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        NFRONT_EFF = NTOTFS
        IBROT = INODE
        DO IORG = 1, NUMORG
          K1 = PTRAIW(IBROT) + 2
          JT1 = INTARR(K1)
          INTARR(K1) = ITLOC(JT1)
          IBROT = FILS(IBROT)
         K2 = K1 + INTARR(K1 - 2) - INTARR(K1 - 1)
         K1 = K1 + 1
         IF (K1 .LE. K2) THEN
          DO JJ = K1, K2
            J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
          ENDDO
         ENDIF
        ENDDO
        K1 = IOLDPS+HF
        DO JJ=K1+NELIM,K1+NFRONT_EFF-1
          ITLOC(IW(JJ)) = 0
        ENDDO
        RETURN   
      ENDIF
       ALLOCATE(PTTRI(NUMSTK+1), stat=allocok)
       IF (allocok .GT. 0) THEN
        IFLAG = -13
        GOTO 800
       ENDIF
       ALLOCATE(PTLAST(NUMSTK+1), stat=allocok)
       IF (allocok .GT. 0) THEN
        IFLAG = -13
        GOTO 800
       ENDIF
      NFRONT_EFF = NASS1
      IF ( ISON_IN_PLACE > 0 ) THEN
        ISON  = ISON_IN_PLACE
        J2    = PIMASTER(STEP(ISON))
        LSTK   = IW(J2    +KEEP(IXSZ))
        NELIM = IW(J2 + 1+KEEP(IXSZ))
        NPIVS  = IW(J2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(J2 + 5+KEEP(IXSZ))
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF ( J2.GT. IWPOSCB ) THEN
          NROWS = IW(J2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        J1 = J2 + HS + NROWS + NPIVS
        J2 = J1 + LSTK - 1
        J3 = J1 + NELIM - 1
        DO JJ = J1, J3
          NTOTFS = NTOTFS + 1
          JT1 = IW(JJ)
          IW(ICT11 + NTOTFS) = JT1
          ITLOC(JT1) = NTOTFS
          IW(JJ) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
        ENDDO
        NELIM_SON_IN_PLACE = NTOTFS
      ENDIF
      IN = INODE
      INEW = IOLDPS + HF +  NTOTFS
      INEW1 = NTOTFS + 1
      JDEBROW = PTRAIW(INODE)+3
      PTTRI(NUMSTK+1)  = JDEBROW
      PTLAST(NUMSTK+1) = JDEBROW + INTARR(JDEBROW-3) - 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW)         = JT1
      IW(INEW+NFRONT)  = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      NTOTFS = NTOTFS + NUMORG
      IF (NUMSTK .NE. 0) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK   = IW(J2    +KEEP(IXSZ))
          NELIM  = IW(J2 + 1+KEEP(IXSZ))
          NPIVS  = IW(J2 + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5+KEEP(IXSZ))
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS  = NPIVS + LSTK
          NROWS  = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+KEEP(IXSZ))
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1 - KEEP(253)
          J3 = J1 + NELIM - 1
          IF (NELIM .NE. 0 .AND. ISON.NE.ISON_IN_PLACE) THEN
            DO JJ = J1, J3
              NTOTFS = NTOTFS + 1
              JT1 = IW(JJ)
              IW(ICT11 + NTOTFS) = JT1
              ITLOC(JT1) = NTOTFS
              IW(JJ) = NTOTFS
              IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
            ENDDO
          ENDIF
          PTTRI(IELL)  = J2+1
          PTLAST(IELL) = J2
          J1 = J3 + 1
          IF (NASS1 .NE. NFRONT - KEEP(253)) THEN
            DO JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN 
                PTTRI(IELL) = JJ
                EXIT
              ENDIF
            ENDDO
          ELSE
            DO JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
            ENDDO
            DO JJ=J2+1, J2+KEEP(253)
              IW(JJ)=NFRONT-KEEP(253)+JJ-J2
            ENDDO
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IF (NFRONT-KEEP(253).EQ.NASS1) GOTO 500
 199  CONTINUE
      IF ( PTTRI( NUMSTK + 1 ) .LE. PTLAST( NUMSTK + 1 ) ) THEN
      IF ( ITLOC( INTARR( PTTRI( NUMSTK + 1 ) ) ) .NE. 0 ) THEN
       PTTRI( NUMSTK + 1 ) = PTTRI( NUMSTK + 1 ) + 1
       GOTO 199
      END IF
      END IF
      MIN_PERM = N + 1
      DO IELL = 1, NUMSTK 
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
      IELL = NUMSTK + 1
      ILOC =  PTTRI( IELL )
      IF ( ILOC .LE. PTLAST( IELL ) ) THEN
        IF ( PERM( INTARR( ILOC ) ) .LT. MIN_PERM ) THEN
         JMIN        = INTARR( ILOC )
         MIN_PERM = PERM( JMIN )
        END IF
      END IF
      NEWEL = IOLDP2 + NASS1 + NFRONT
      DO WHILE ( MIN_PERM .NE. N + 1 )
          NEWEL  = NEWEL + 1
          NFRONT_EFF = NFRONT_EFF + 1
          IW( NEWEL ) = JMIN
          ITLOC( JMIN ) = NFRONT_EFF
          LAST_J_ASS = JMIN
          MIN_PERM = N + 1
          DO IELL = 1,  NUMSTK
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
              IF ( IW( PTTRI( IELL ) ) .eq. LAST_J_ASS )
     &        PTTRI( IELL ) = PTTRI( IELL ) + 1
            ENDIF
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
             IF ( PERM(IW( PTTRI( IELL )) ) .LT. MIN_PERM ) THEN
                JMIN        = IW( PTTRI( IELL ) )
                MIN_PERM = PERM( JMIN )
             END IF
            END IF
          END DO
          IELL = NUMSTK + 1
 145      CONTINUE
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
            IF ( INTARR( PTTRI( IELL ) ) .eq. LAST_J_ASS ) THEN
              PTTRI( IELL ) = PTTRI( IELL ) + 1 
              GOTO 145
            END IF
          END IF
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
            IF (PERM(INTARR( PTTRI(IELL) )) .LT. MIN_PERM) THEN
              JMIN        = INTARR( PTTRI(IELL) )
              MIN_PERM = PERM( JMIN )
            END IF
          END IF
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT - KEEP(253)) THEN 
       IBROT = INODE
       DO IORG = 1, NUMORG
         J1    = PTRAIW(IBROT) + 2
         J2    = J1 + INTARR(J1 - 2) - INTARR(J1-1)
         IBROT = FILS( IBROT )
         IF ( IORG.EQ. 1) THEN
           IF ( KEEP(50).NE.0 ) CYCLE
           J1 = J1 + 1 + INTARR(J1-2)
         ELSE
           J1 = J1 + 1
         ENDIF
         DO JJ = J1, J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
       ENDDO
       IF ( (TYPESPLIT.EQ.4).AND.
     &      (NFRONT_EFF.LT.NFRONT-KEEP(253)) ) THEN
         IBROT = INODE
         DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )  
          IBROT = DAD(STEP(IBROT)) 
          IN = IBROT
          DO WHILE (IN.GT.0.AND.NFRONT_EFF.LT.NFRONT-KEEP(253))
            J1    = PTRAIW(IN) + 2
            J2    = J1 + INTARR(J1 - 2) - INTARR(J1-1)
            IN = FILS( IN )
            DO JJ = J1, J2
              J     = INTARR( JJ )
              IF ( ITLOC( J ) .eq. 0 ) THEN
               NEWEL  = NEWEL + 1
               NFRONT_EFF = NFRONT_EFF + 1
               IW( NEWEL ) = J
               ITLOC( J ) = NFRONT_EFF
              END IF
            ENDDO
          ENDDO
          IF (NFRONT_EFF.EQ.NFRONT-KEEP(253)) EXIT  
         ENDDO
       ENDIF
      ENDIF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
         DO JJ=NASS1+1, NFRONT_EFF
           IW( IOLDP2+JJ ) = IW( ICT11+JJ )
         ENDDO
      ELSE
        CALL MUMPS_308( N, PERM, 
     &           IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
        CALL MUMPS_309( N, NASS1, PERM, ITLOC,
     &    IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE, 
     &    IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1, 
     &    IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 ) 
        DO JJ = NASS1+1, NFRONT_EFF
          IW(ICT11 + JJ) = IW(IOLDP2+JJ)
        ENDDO
      END IF
  500 CONTINUE
      IF ( KEEP(253).GT.0) THEN
        IP1 = IOLDPS +  HF + NFRONT_EFF  
        IP2 = IOLDPS + HF + NFRONT + NFRONT_EFF 
        DO I= 1, KEEP(253)
          IW(IP1+I-1) = N+I
          IW(IP2+I-1) = N+I
          ITLOC(N+I)  = NFRONT_EFF + I
        ENDDO
        NFRONT_EFF = NFRONT_EFF + KEEP(253)
      ENDIF
      IF (NFRONT.NE.NFRONT_EFF) THEN
        IP1 = IOLDPS + NFRONT + HF 
        IP2 = IOLDPS + NFRONT_EFF + HF
        DO I=1, NFRONT_EFF
          IW(IP2+I-1)=IW(IP1+I-1)
        ENDDO
      ENDIF
      IF ((NUMSTK .NE. 0).AND.(NFRONT-KEEP(253).GT.NASS1)) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+KEEP(IXSZ))
          NELIM = IW(J2 + 1+KEEP(IXSZ))
          NPIVS = IW(J2 + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5+KEEP(IXSZ))
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+KEEP(IXSZ))
          ENDIF
          HS = NSLSON + 6 +KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          J1 = J3 + 1
          DO JJ = J1, J2
              J = IW(JJ)
                IW(JJ) = ITLOC(J)
          ENDDO
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IBROT = INODE
      DO IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO JJ = J1, J2
            J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
          ENDDO
        ENDIF
      ENDDO
        K1 = IOLDPS + HF
        K2 = K1 + NFRONT_EFF -1
        IF (KEEP(50).EQ.0) K2 = K2 + NELIM_SON_IN_PLACE
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
        IF (KEEP(50).EQ.0) THEN
          K1 = IOLDPS+HF+NFRONT_EFF+NELIM_SON_IN_PLACE+NUMORG
          K2 = K1 + NASS -NELIM_SON_IN_PLACE - 1
          DO K = K1, K2
            I = IW(K)
            ITLOC(I) = 0
          ENDDO
        ENDIF
  800 CONTINUE
      IF (allocated(PTTRI)) DEALLOCATE(PTTRI)
      IF (allocated(PTLAST)) DEALLOCATE(PTLAST)
      RETURN
      END SUBROUTINE MUMPS_86
      SUBROUTINE MUMPS_308( N, PERM, IW, LIW )
      IMPLICIT NONE
      INTEGER N, LIW
      INTEGER PERM( N ), IW( LIW )
      INTEGER I, SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, LIW - 1
          IF ( PERM( IW( I ) ) .GT. PERM( IW( I + 1 ) ) ) THEN
            DONE = .FALSE.
            SWAP  = IW( I + 1 )
            IW( I + 1 ) = IW( I )
            IW( I ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_308
      SUBROUTINE MUMPS_309( N, NASS1, PERM, ITLOC,
     &                             SMALL, LSMALL,
     &                             LARGE, LLARGE,
     &                             MERGE, LMERGE )
      IMPLICIT NONE
      INTEGER N, NASS1, LSMALL, LLARGE, LMERGE
      INTEGER PERM( N ), ITLOC( N ) 
      INTEGER SMALL(LSMALL), LARGE(LLARGE), MERGE(LMERGE)
      INTEGER PSMALL, PLARGE, PMERGE, VSMALL, VLARGE, VMERGE
      PSMALL = 1
      PLARGE = 1
      PMERGE = 1
      DO WHILE ( PSMALL .LE. LSMALL .or. PLARGE.LE. LLARGE )
        IF ( PSMALL .GT. LSMALL ) THEN
          VMERGE = LARGE( PLARGE )
          PLARGE = PLARGE + 1
        ELSE IF ( PLARGE .GT. LLARGE ) THEN
          VMERGE = SMALL( PSMALL )
          PSMALL = PSMALL + 1
        ELSE
          VSMALL = SMALL( PSMALL )
          VLARGE = LARGE( PLARGE )
          IF ( PERM( VSMALL ) .LT. PERM( VLARGE ) ) THEN
            VMERGE = VSMALL
            PSMALL   = PSMALL + 1
          ELSE
            VMERGE = VLARGE
            PLARGE   = PLARGE + 1
          END IF
        END IF
        MERGE( PMERGE ) = VMERGE
        ITLOC( VMERGE ) = PMERGE + NASS1
        PMERGE = PMERGE + 1
      END DO
      PMERGE = PMERGE - 1
      RETURN
      END SUBROUTINE MUMPS_309
      SUBROUTINE MUMPS_125(
     &           NUMELT, LIST_ELT, 
     &           MYID, INODE, N, IOLDPS,
     &           HF, NFRONT, NFRONT_EFF, PERM, 
     &           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     &           IFSON, STEP, PIMASTER, PTRAIW, NELT, 
     &           IW, LIW, 
     &           INTARR, LINTARR, ITLOC, RHS_MUMPS,
     &           FILS, FRERE_STEPS, 
     &           KEEP, SON_LEVEL2, NIV1, NBPROCFILS, IFLAG, 
     &           DAD, PROCNODE_STEPS, SLAVEF, 
     &           FRT_PTR, FRT_ELT, Pos_First_NUMORG)
      IMPLICIT NONE
      INTEGER NELT, INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     &        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG,
     &        LINTARR, NUMELT
      INTEGER KEEP(500)
      INTEGER LIST_ELT(*)
      INTEGER STEP(N), PIMASTER(KEEP(28)),
     &  PTRAIW(NELT+1), IW(LIW), 
     &        ITLOC(N+KEEP(253)), FILS(N), FRERE_STEPS(KEEP(28)),
     &        NBPROCFILS(KEEP(28)), PERM(N)
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(LINTARR)
      LOGICAL, intent(in)    :: NIV1
      LOGICAL, intent(out)   :: SON_LEVEL2
      INTEGER, intent(out)   :: NFRONT_EFF
      INTEGER, intent(in)    :: DAD (KEEP(28))
      INTEGER, intent(in) :: PROCNODE_STEPS(KEEP(28)), SLAVEF
      INTEGER, intent(in) :: FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER, intent(out) :: Pos_First_NUMORG
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     &        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     &        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     &        NROWS, HS, IP1, IP2, K1, K2, IBROT,
     &        I, K, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     &        LAST_J_ASS, JMIN, MIN_PERM
      INTEGER TYPESPLIT, NUMELT_IBROT
      INTEGER ELTI
      INCLUDE 'mumps_headers.h'
      LOGICAL LEVEL1
      INTEGER allocok
      INTEGER , ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      INTEGER  MUMPS_810, MUMPS_330
      EXTERNAL MUMPS_810, MUMPS_330 
      Pos_First_NUMORG = 1
      TYPESPLIT  = MUMPS_810 (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NFRONT_EFF = NASS1
      NTOTFS     = 0
      IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) ) THEN
        J2    = PIMASTER(STEP(IFSON))
        LSTK  = IW(J2    +KEEP(IXSZ))
        NELIM = IW(J2 + 1+KEEP(IXSZ))
        NPIVS  = IW(J2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(J2 + 5+KEEP(IXSZ))
        IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
        LEVEL1    = NSLSON.EQ.0
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF (NIV1) THEN
         write(6,*) MYID, ':',
     &    ' Internal error 2 in MUMPS_86 ',
     &    ' interior split node of type 1 '
         CALL MUMPS_ABORT()
        ELSE
         I= MUMPS_330(PROCNODE_STEPS(STEP(IFSON)),SLAVEF)
         J= MUMPS_810(PROCNODE_STEPS(STEP(IFSON)), 
     &              SLAVEF)
         IF (LEVEL1.or.J.LT.4) THEN
         write(6,*) MYID, ':',
     &     ' Internal error 3 in MUMPS_86 ',
     &     ' son', IFSON, 
     &     ' of interior split node', INODE, ' of type 1 ', 
     &     ' NSLSON =', NSLSON, ' TYPE_SON=', I, 'TYPESPLIT_SON=', J
          CALL MUMPS_ABORT()
         ELSE
          NBPROCFILS(STEP(IFSON)) = NSLSON
          NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(IFSON))
         ENDIF
        ENDIF
        IF ( J2.GT. IWPOSCB ) THEN
          NROWS = IW(J2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        J1 = J2 + HS + NROWS + NPIVS
        J2 = J1 + LSTK - 1
        J3 = J1 + NELIM - 1
        IF (NELIM.GT.0) THEN
         DO JJ=J1,J3
          NTOTFS = NTOTFS + 1
          JT1 = IW(JJ)
          IW(ICT11 + NTOTFS) = JT1
          IW(JJ) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
         ENDDO
        ENDIF
        DO JJ =J3+1, J2
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS    
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        NFRONT_EFF = NTOTFS
        DO IELL=1,NUMELT
         ELTI = LIST_ELT(IELL)
         J1= PTRAIW(ELTI)
         J2= PTRAIW(ELTI+1)-1
         DO JJ=J1,J2
          J = INTARR(JJ)
          INTARR(JJ) = ITLOC(J)
         END DO
        ENDDO
        Pos_First_NUMORG = ITLOC(INODE)
        K1 = IOLDPS+HF
        DO JJ=K1+NELIM,K1+NFRONT_EFF-1
          ITLOC(IW(JJ)) = 0
        ENDDO
        RETURN   
      ENDIF
      IF (NUMSTK.GT.0) THEN
        ALLOCATE(PTTRI(NUMSTK), stat=allocok)
        IF (allocok .GT. 0) THEN
         IFLAG = -13
         GOTO 800
        ENDIF
        ALLOCATE(PTLAST(NUMSTK), stat=allocok)
        IF (allocok .GT. 0) THEN
         IFLAG = -13
         GOTO 800
        ENDIF        
      ENDIF
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      DO WHILE (IN.GT.0)
       ITLOC(IN)        = INEW1
       IW(INEW)         = IN
       IW(INEW+NFRONT)  = IN
       INEW1     = INEW1 + 1
       INEW      = INEW + 1
       IN = FILS(IN)
      END DO
      NTOTFS = NUMORG
      IF (NUMSTK .NE. 0) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK   = IW(J2    +KEEP(IXSZ))
          NELIM  = IW(J2 + 1+KEEP(IXSZ))
          NPIVS  = IW(J2 + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5+KEEP(IXSZ))
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS  = NPIVS + LSTK
          NROWS  = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+KEEP(IXSZ))
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1 - KEEP(253)
          J3 = J1 + NELIM - 1
          IF (NELIM .NE. 0) THEN
            DO JJ = J1, J3
              NTOTFS = NTOTFS + 1
              JT1 = IW(JJ)
              IW(ICT11 + NTOTFS) = JT1
              ITLOC(JT1) = NTOTFS
              IW(JJ) = NTOTFS
              IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
            ENDDO
          ENDIF
          PTTRI(IELL)  = J2+1
          PTLAST(IELL) = J2
          J1 = J3 + 1
          IF (NASS1 .NE. NFRONT - KEEP(253)) THEN
            DO JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN 
                PTTRI(IELL) = JJ
                EXIT
              ENDIF
            ENDDO
          ELSE
            DO JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
            ENDDO
           DO JJ=J2+1, J2+KEEP(253)
             IW(JJ)=NFRONT-KEEP(253)+JJ-J2
           ENDDO
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IF (NFRONT-KEEP(253).EQ.NASS1) GOTO 500
      MIN_PERM = N + 1
      JMIN     = -1
      DO IELL = 1, NUMSTK 
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
      NEWEL = IOLDP2 + NASS1 + NFRONT
      DO WHILE ( MIN_PERM .NE. N + 1 )
          NEWEL  = NEWEL + 1
          NFRONT_EFF = NFRONT_EFF + 1
          IW( NEWEL ) = JMIN
          ITLOC( JMIN ) = NFRONT_EFF
          LAST_J_ASS = JMIN
          MIN_PERM = N + 1
          DO IELL = 1,  NUMSTK
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
              IF ( IW( PTTRI( IELL ) ) .eq. LAST_J_ASS )
     &        PTTRI( IELL ) = PTTRI( IELL ) + 1
            ENDIF
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
             IF ( PERM(IW( PTTRI( IELL )) ) .LT. MIN_PERM ) THEN
                JMIN        = IW( PTTRI( IELL ) )
                MIN_PERM = PERM( JMIN )
             END IF
            END IF
          END DO
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT-KEEP(253)) THEN 
       DO IELL = 1,NUMELT
        ELTI = LIST_ELT(IELL)
         J1= PTRAIW(ELTI)
         J2= PTRAIW(ELTI+1)-1
         DO JJ=J1,J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
       ENDDO
       IF ( (TYPESPLIT.EQ.4).AND.
     &      (NFRONT_EFF.LT.NFRONT-KEEP(253)) ) THEN
         IBROT = INODE
         DO WHILE 
     &      (
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_810 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.6  
     &        ) 
     &      )     
          IBROT = DAD(STEP(IBROT)) 
          NUMELT_IBROT = FRT_PTR(IBROT+1) - FRT_PTR(IBROT)
          IF (NUMELT_IBROT.EQ.0) CYCLE 
          DO IELL = FRT_PTR(IBROT), FRT_PTR(IBROT+1)
            ELTI = FRT_ELT(IELL)
            J1= PTRAIW(ELTI)
            J2= PTRAIW(ELTI+1)-1
            DO JJ=J1,J2
              J     = INTARR( JJ )
              IF ( ITLOC( J ) .eq. 0 ) THEN
               NEWEL  = NEWEL + 1
               NFRONT_EFF = NFRONT_EFF + 1
               IW( NEWEL ) = J
               ITLOC( J ) = NFRONT_EFF
              END IF
            ENDDO
          ENDDO
          IF (NFRONT_EFF.EQ.NFRONT-KEEP(253)) EXIT  
         ENDDO
       ENDIF
      END IF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
         DO JJ=NASS1+1, NFRONT_EFF
           IW( IOLDP2+JJ ) = IW( ICT11+JJ )
         ENDDO
      ELSE
          CALL MUMPS_308( N, PERM, 
     &            IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
          CALL MUMPS_309( N, NASS1, PERM, ITLOC,
     &               IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE,
     &               IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1,
     &               IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 )
        DO JJ = NASS1+1, NFRONT_EFF
          IW(ICT11 + JJ) = IW(IOLDP2+JJ)
        ENDDO
      END IF
  500 CONTINUE
      IF ( KEEP(253).GT.0) THEN
        IP1 = IOLDPS +  HF + NFRONT_EFF  
        IP2 = IOLDPS + HF + NFRONT + NFRONT_EFF 
        DO I= 1, KEEP(253)
          IW(IP1+I-1) = N+I
          IW(IP2+I-1) = N+I
          ITLOC(N+I)  = NFRONT_EFF + I
        ENDDO
        NFRONT_EFF = NFRONT_EFF + KEEP(253)
      ENDIF
      IF (NFRONT.GT.NFRONT_EFF) THEN
        IP1 = IOLDPS + NFRONT + HF
        IP2 = IOLDPS + NFRONT_EFF + HF
        DO I=1,NFRONT_EFF
          IW(IP2+I)=IW(IP1+I)
        ENDDO
      ELSE IF (NFRONT .LT. NFRONT_EFF) THEN
        WRITE(*,*) "Internal error in MUMPS_125",
     &             NFRONT, NFRONT_EFF
        CALL MUMPS_ABORT()
      ENDIF
      IF ((NUMSTK .NE. 0).AND.
     &    (NFRONT-KEEP(253).GT.NASS1)) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+KEEP(IXSZ))
          NELIM = IW(J2 + 1 +KEEP(IXSZ))
          NPIVS = IW(J2 + 3 +KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5 +KEEP(IXSZ))
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+KEEP(IXSZ))
          ENDIF
          HS = NSLSON + 6 +KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          J1 = J3 + 1
          DO JJ = J1, J2
              J = IW(JJ)
                IW(JJ) = ITLOC(J)
          ENDDO
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      DO IELL=1,NUMELT
       ELTI = LIST_ELT(IELL)
       J1= PTRAIW(ELTI)
       J2= PTRAIW(ELTI+1)-1
       DO JJ=J1,J2
          J = INTARR(JJ)
          INTARR(JJ) = ITLOC(J)
       END DO
      ENDDO
        K1 = IOLDPS + HF + NUMORG
        K2 = K1 + NFRONT_EFF - 1 + NASS
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
  800 CONTINUE
      IF (allocated(PTTRI)) DEALLOCATE(PTTRI)
      IF (allocated(PTLAST)) DEALLOCATE(PTLAST)
      RETURN
      END SUBROUTINE MUMPS_125
      INTEGER FUNCTION MUMPS_50 
     &         ( SLAVEF, K48, K821, K50,
     &         NFRONT, NCB)
      IMPLICIT NONE
      INTEGER,    INTENT (IN) :: SLAVEF, K48, K50, NFRONT, NCB
      INTEGER(8), INTENT (IN) :: K821
      INTEGER NSLAVESMIN, NASS, KMAX
      REAL Wmaster, Wtotal, Wmax
      INTEGER ACC,X
      REAL MUMPS_45
      INTEGER MUMPS_497    
      EXTERNAL  MUMPS_45, MUMPS_497
      KMAX = MUMPS_497( K821, NCB )
      NASS = NFRONT - NCB
      NSLAVESMIN  = 1
      IF ( K48 .EQ.0 .OR. (K48.EQ.5 .AND.K50.EQ.0)) THEN 
         NSLAVESMIN = max(NCB/max(1,KMAX),1)
      ELSE IF (K48 .EQ. 3 .OR.(K48.EQ.5 .AND.K50.NE.0) ) THEN
        Wmax    = MUMPS_45(KMAX,NFRONT,NASS)
        Wtotal  = MUMPS_45(NCB,NFRONT,NASS)
        Wmaster = real(NASS*NASS)*real(NASS)/(3.0)
        IF ( Wmaster .GT.  Wmax ) THEN 
         NSLAVESMIN = max ( nint ( Wtotal / Wmaster ), 1 )
        ELSE
         NSLAVESMIN = max ( nint ( Wtotal / Wmax ), 1 )
        ENDIF
        IF (K48 .EQ. 5) THEN
          NSLAVESMIN = max ( NSLAVESMIN/2, 1 )
        END IF
      ELSE IF (K48 .EQ. 4 ) THEN
         IF ( K821 > 0_8 ) THEN
           WRITE(*,*) 'Internal Error 1 in MUMPS_50'
           CALL MUMPS_ABORT()
         ENDIF
         CALL MUMPS_ABORT_ON_OVERFLOW(K821,
     &           "K821 too large in MUMPS_50" )
         KMAX=int(abs(K821))
         IF(K50.EQ.0)THEN
            NSLAVESMIN = max(int(
     &                (int(NCB,8)*int(NCB,8))/int(KMAX,8)
     &                  ),1)
         ELSE
            ACC=0
            NSLAVESMIN=0
            DO WHILE (ACC.NE.NCB)
               X=int((-real(NFRONT-NCB+ACC)
     &              +sqrt(((real(NFRONT-NCB+ACC)*
     &              real(NFRONT-NCB+ACC))+real(4)*
     &              real(KMAX))))/
     &              real(2))
               ACC=ACC+X
               NSLAVESMIN=NSLAVESMIN+1
               IF (((NCB-ACC)*NCB).LT.KMAX)THEN
                  ACC=NCB 
                  NSLAVESMIN=NSLAVESMIN+1
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      NSLAVESMIN = min ( NSLAVESMIN,(SLAVEF-1) )
      MUMPS_50 = 
     &               min ( NSLAVESMIN, NCB )
      RETURN
      END FUNCTION MUMPS_50 
      INTEGER FUNCTION MUMPS_52 
     &        ( SLAVEF, K48, K821, K50,
     &          NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K50,NFRONT, NCB
      INTEGER(8), INTENT(IN) :: K821
      INTEGER NSLAVESMAX, KMAX, KMIN
      INTEGER NSLAVESMIN
      INTEGER MUMPS_497,MUMPS_442,
     &        MUMPS_50,
     &        MUMPS_46
      EXTERNAL MUMPS_497,MUMPS_442,
     &        MUMPS_50,
     &        MUMPS_46
      IF (K48 .eq. 0 .OR. K48.eq.3.OR.K48.EQ.5) THEN
         KMAX = MUMPS_497( K821, NCB )
         KMIN = MUMPS_442( K821, K50, KMAX, NCB)
         NSLAVESMAX = MUMPS_46(
     &                SLAVEF, K48, K50, KMIN, NFRONT, NCB )
      ELSE
         NSLAVESMAX = SLAVEF-1
      ENDIF
      NSLAVESMIN = MUMPS_50(
     &     SLAVEF, K48, K821, K50, NFRONT, NCB )
      NSLAVESMAX = max ( NSLAVESMAX, NSLAVESMIN )
      MUMPS_52 = 
     &               min ( NSLAVESMAX, NCB )
      RETURN
      END FUNCTION MUMPS_52
      SUBROUTINE MUMPS_503( WHAT, KEEP,KEEP8,
     &           NCB, NFR, SLAVEF, NBROWMAX, MAXSURFCB8
     &     )
      IMPLICIT NONE
      INTEGER, intent(in) :: WHAT, NCB, NFR, SLAVEF
      INTEGER, intent(in) :: KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, intent(out) :: NBROWMAX
      INTEGER(8), intent(out) :: MAXSURFCB8
      INTEGER KMAX, KMIN, NSLAVES, SIZEDUMMY, TABDUMMY(1)
      EXTERNAL MUMPS_497, MUMPS_442,
     &         MUMPS_50
      INTEGER MUMPS_497, MUMPS_442,
     &        MUMPS_50
      IF ( WHAT .NE. 1 .and. WHAT .NE. 2 ) THEN
        IF (WHAT .NE. 4 .and. WHAT .NE. 5 .AND.
     &       KEEP(48).NE.5 ) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_503"
        CALL MUMPS_ABORT()
        END IF
      ENDIF
      KMAX    = MUMPS_497( KEEP8(21), NCB )
      IF (WHAT .EQ.1.OR.WHAT.EQ.2) THEN
        NSLAVES = MUMPS_50( SLAVEF, KEEP(48),
     &            KEEP8(21), KEEP(50),
     &            NFR, NCB )
      ELSE
        NSLAVES=SLAVEF
      ENDIF
      IF ( KEEP(48) == 0 .OR. (KEEP(48).EQ.5.AND.KEEP(50).EQ.0)) THEN
        NBROWMAX = NCB / NSLAVES + mod( NCB, NSLAVES )
        IF ( WHAT == 2 .OR. WHAT == 5 )
     &    MAXSURFCB8 = int(NBROWMAX,8) * int(NCB,8)
      ELSE IF (KEEP(48) == 3.OR.(KEEP(48).EQ.5.AND.KEEP(50).NE.0))THEN
        KMIN = MUMPS_442( KEEP8(21), KEEP(50), KMAX, NCB )  
        SIZEDUMMY        = 1
        IF (WHAT.GT.3) THEN
           CALL  MUMPS_440(
     &          WHAT-3, NSLAVES, NFR, NCB, 
     &          KMIN, KMAX, SLAVEF,
     &          NBROWMAX, MAXSURFCB8, TABDUMMY, SIZEDUMMY) 
        ELSE
           CALL  MUMPS_440(
     &          WHAT, NSLAVES, NFR, NCB, 
     &          KMIN, KMAX, SLAVEF,
     &          NBROWMAX, MAXSURFCB8, TABDUMMY, SIZEDUMMY) 
        ENDIF
      ELSE IF ( KEEP(48) == 4 ) THEN
         IF (KEEP8(21) > 0_8) THEN
            WRITE(*,*) "Internal error 2 in MUMPS_503"
            CALL MUMPS_ABORT()
         END IF
         IF(KEEP(50).EQ.0)THEN
            IF ( abs(KEEP8(21)) * int( SLAVEF - 1,8 ) >
     &                            int( NCB,8) * int(NFR,8) ) THEN
              NBROWMAX = (NCB + SLAVEF -2 ) / ( SLAVEF - 1 )
              IF ( WHAT == 2 ) MAXSURFCB8 = int(NBROWMAX,8) *int(NCB,8)
            ELSE
              NBROWMAX=int(
     &                      (abs(KEEP8(21)) + int(NFR - 1,8))
     &                    /  int(NFR,8)
     &                    )
              IF ( WHAT == 2 ) MAXSURFCB8 = abs(KEEP8(21))
            ENDIF
         ELSE
            NBROWMAX=int((-real(NFR-NCB)
     &              +sqrt((real(NFR-NCB)*
     &              real(NFR-NCB))+real(4)*
     &              real(abs(KEEP8(21)))))/
     &              real(2))
            IF ( WHAT == 2 ) MAXSURFCB8 = abs(KEEP8(21))
         ENDIF
      ELSE
        NBROWMAX = NCB
        IF (WHAT == 2) MAXSURFCB8 = int(NCB,8) * int(NCB,8)
      ENDIF 
      NBROWMAX = min ( max(NBROWMAX, 1), NCB)
      RETURN
      END SUBROUTINE MUMPS_503
      INTEGER FUNCTION MUMPS_46( SLAVEF, K48, K50,
     &         BLSIZE, NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K50, BLSIZE, NFRONT, NCB
      INTEGER NSLAVES, NASS
      REAL Wtotal, Wblsize
      REAL MUMPS_45
      EXTERNAL          MUMPS_45
      NASS = NFRONT - NCB
      NSLAVES  = SLAVEF-1
      IF ( K48 .EQ.0 .OR. (K48.EQ.5 .AND. K50.EQ.0)) THEN 
         NSLAVES = max(NCB/max(1,BLSIZE),1)
      ELSE IF (K48.EQ.3 .OR. (K48.EQ.5 .AND. K50.NE.0))THEN
        Wblsize = MUMPS_45(BLSIZE,NFRONT,NASS)
        Wtotal  = MUMPS_45(NCB,NFRONT,NASS)
        NSLAVES = max(nint ( Wtotal / Wblsize ), 1)
      ENDIF
      MUMPS_46 = 
     &               min ( NSLAVES,(SLAVEF-1) )
      RETURN
      END FUNCTION MUMPS_46
      SUBROUTINE  MUMPS_440( 
     &    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     &    KMIN, KMAX, SLAVEF,
     &    NBROWMAX, MAXSURFCB, TABPOS, SIZETABPOS) 
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: GETPOSITIONS, 
     &    NSLAVES, NFRONT, NCB, 
     &    KMIN, KMAX, SLAVEF, SIZETABPOS 
      INTEGER, INTENT (OUT) :: NBROWMAX
      INTEGER(8), INTENT(OUT) :: MAXSURFCB
      INTEGER, INTENT (OUT) :: TABPOS(SIZETABPOS) 
      REAL W, COSTni
      REAL delta
      INTEGER  SumNi, NCOLim1, I, BLSIZE, NASS
      LOGICAL GETROW, GETSURF, GETPOS, GET_AVGROW, GET_AVGSURF
      REAL MUMPS_45
      EXTERNAL          MUMPS_45
      GETROW = (GETPOSITIONS.EQ.1)
      GETSURF= (GETPOSITIONS.EQ.2)
      GETPOS = (GETPOSITIONS.EQ.3)
      GET_AVGROW = (GETPOSITIONS.EQ.4)
      GET_AVGSURF = (GETPOSITIONS.EQ.5)
      NBROWMAX  = 0
      MAXSURFCB = 0_8
      IF (GETPOS) THEN 
        TABPOS (1) = 1
        TABPOS (NSLAVES+1)= NCB+1
        TABPOS (SLAVEF+2) = NSLAVES
      ENDIF
      IF (NSLAVES.EQ.1) THEN
       IF ( GETSURF ) THEN
         NBROWMAX  = NCB
         MAXSURFCB = int(NCB,8)*int(NCB,8)
       ELSEIF ( GETROW ) THEN
         NBROWMAX  = NCB
       ENDIF
      ELSE
        NASS    = NFRONT - NCB
        W       = MUMPS_45(NCB,NFRONT,NASS)
        SumNi   = 0
        NCOLim1 = NASS
        DO I = 1, NSLAVES-1
          delta   = real(2*NCOLim1-NASS+1)**2 + 
     &                  (real(4)*W)/real(NASS*(NSLAVES-I+1))
          delta   = sqrt(delta)
          delta   = (real(-2*NCOLim1+NASS-1) + delta )/real(2)
          BLSIZE  = max(int(delta), 1)
          IF ( (NFRONT-NCOLim1-BLSIZE) .LE. NSLAVES-I ) THEN
            BLSIZE = 1
          ENDIF
          NCOLim1 = NCOLim1+BLSIZE
          COSTni  = MUMPS_45(BLSIZE,NCOLim1,NASS)
          W       = W - COSTni
          IF (GETPOS) TABPOS(I) = SumNi + 1
          IF (GETSURF) THEN
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
            MAXSURFCB = max ( MAXSURFCB, 
     &       int(BLSIZE,8)* int(SumNi+BLSIZE,8) )
          ELSEIF ( GETROW ) THEN         
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
             RETURN
          ELSEIF (GET_AVGSURF) THEN
            NBROWMAX = NBROWMAX + BLSIZE
            MAXSURFCB = MAXSURFCB + int(BLSIZE,8)*int(SumNi+BLSIZE,8) 
          ELSEIF (GET_AVGROW) THEN
             NBROWMAX = NBROWMAX + BLSIZE
          ENDIF
          SumNi   = SumNi + BLSIZE
        ENDDO
        BLSIZE = NCB - SumNi
        IF (BLSIZE.LE.0) THEN
          write(*,*) ' Error in MUMPS_440: ', 
     &     ' size lastbloc ', BLSIZE
          CALL MUMPS_ABORT()
        ENDIF
        if (NCOLim1+BLSIZE.NE.NFRONT) then
          write(*,*) ' Error in MUMPS_440: ', 
     &     ' NCOLim1, BLSIZE, NFRONT=', 
     &       NCOLim1, BLSIZE, NFRONT
          CALL MUMPS_ABORT()
        endif
        IF (GETPOS) TABPOS(NSLAVES) = SumNi + 1
        IF (GETSURF) THEN
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
            MAXSURFCB = max ( MAXSURFCB, 
     &       int(BLSIZE,8)* int(SumNi+BLSIZE,8 ))
        ELSEIF ( GETROW ) THEN         
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
        ELSEIF (GET_AVGSURF) THEN
          NBROWMAX = NBROWMAX + BLSIZE
          MAXSURFCB = MAXSURFCB + int(BLSIZE,8)*int(SumNi+BLSIZE,8)
          NBROWMAX=(NBROWMAX+NSLAVES-1)/NSLAVES
          MAXSURFCB=(MAXSURFCB+int(NSLAVES-1,8))/int(NSLAVES,8)
        ELSEIF (GET_AVGROW) THEN
          NBROWMAX = NBROWMAX + BLSIZE
          NBROWMAX=(NBROWMAX+NSLAVES-1)/NSLAVES
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_440
      SUBROUTINE MUMPS_441(      
     &            KEEP,KEEP8, SLAVEF,
     &            TAB_POS_IN_PERE,
     &            NSLAVES, NFRONT, NCB
     &             )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NSLAVES, SLAVEF, NFRONT,
     &                         KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER TAB_POS_IN_PERE(SLAVEF+2)
      INTEGER :: I, BLSIZE
      INTEGER KMIN, KMAX, NBROWDUMMY,
     &        GETPOSITIONS, SIZECOLTAB
      INTEGER(8) MAXSURFDUMMY8
      INTEGER MUMPS_442, MUMPS_497 
      EXTERNAL MUMPS_442, MUMPS_497, 
     &        MUMPS_440
       IF (KEEP(48).EQ.0) THEN
        BLSIZE = NCB / NSLAVES
        TAB_POS_IN_PERE( 1 ) = 1
        DO I = 1, NSLAVES-1
          TAB_POS_IN_PERE( I+1 ) = TAB_POS_IN_PERE(I) +
     &    BLSIZE
        ENDDO
        TAB_POS_IN_PERE(NSLAVES+1) = NCB+1
        TAB_POS_IN_PERE(SLAVEF+2)  = NSLAVES
        RETURN
      ELSE IF (KEEP(48).EQ.3 ) THEN
        KMAX = MUMPS_497(KEEP8(21), NCB)
        KMIN = MUMPS_442(KEEP8(21), KEEP(50), KMAX, NCB)
        GETPOSITIONS = 3
        SIZECOLTAB       = SLAVEF+2
        CALL  MUMPS_440(
     &    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     &    KMIN, KMAX, SLAVEF,
     &    NBROWDUMMY, MAXSURFDUMMY8,
     &    TAB_POS_IN_PERE(1), SIZECOLTAB)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_441
      SUBROUTINE MUMPS_49(
     &            KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &            ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &            ISLAVE, NCB, NSLAVES, SIZE, FIRST_INDEX )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: ISLAVE, NCB, NSLAVES, SLAVEF,
     &                         KEEP(500), INODE, N 
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT( IN ) :: STEP(N), 
     &          ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER, INTENT( OUT ):: SIZE, FIRST_INDEX
      INTEGER BLSIZE, J
      IF (KEEP(48).EQ.0) THEN
       BLSIZE = NCB / NSLAVES
       IF ( ISLAVE .NE. NSLAVES ) THEN
        SIZE = BLSIZE
       ELSE
        SIZE = BLSIZE + mod( NCB, NSLAVES )
       END IF
       FIRST_INDEX = ( ISLAVE - 1 ) * BLSIZE + 1
      ELSEIF (KEEP(48).EQ.3) THEN
       J = ISTEP_TO_INIV2 ( STEP(INODE) )
       FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
       SIZE        = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSEIF (KEEP(48).EQ.4) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
         SIZE = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSEIF (KEEP(48).EQ.5) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
         SIZE = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSE
       WRITE(*,*) 'Error in MUMPS_BLOC2 undef strat'  
       CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_49
      REAL FUNCTION MUMPS_45(NROW,NCOL,NASS)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: NROW,NCOL,NASS
      MUMPS_45 = real(NASS*NROW)*
     &                 real(2*NCOL - NASS - NROW + 1)
      RETURN
      END FUNCTION MUMPS_45
      INTEGER FUNCTION MUMPS_12 
     &      (K821, K48, K50, SLAVEF, 
     &      NCB, NFRONT, NSLAVES_less, NMB_OF_CAND )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NFRONT, NSLAVES_less, 
     &                  K48, K50, SLAVEF, NMB_OF_CAND
      INTEGER(8), INTENT(IN) :: K821
      INTEGER NSLAVES
      INTEGER KMAX, NPIV,
     &        NSLAVES_ref, NSLAVES_max
      REAL WK_MASTER, WK_SLAVE
      INTEGER  MUMPS_497, MUMPS_50, 
     &         MUMPS_52
      REAL  MUMPS_45
      EXTERNAL MUMPS_497, MUMPS_50, 
     &         MUMPS_52
      EXTERNAL MUMPS_45
      IF (NMB_OF_CAND.LE.0) THEN
      ENDIF
      IF ( (K48.EQ.0).OR. (K48.EQ.3) ) THEN
         KMAX = MUMPS_497( K821, NCB )
         NSLAVES_ref = MUMPS_50(
     &     SLAVEF, K48, K821, K50, NFRONT, NCB )
         NSLAVES = NSLAVES_ref
         IF ( NSLAVES_ref.LT.SLAVEF ) THEN
           NSLAVES_max = MUMPS_52(
     &       SLAVEF, K48, K821, K50, NFRONT, NCB )
           IF ( NSLAVES_max .LT. NSLAVES_less ) THEN
            NSLAVES =  NSLAVES_max
           ELSE 
            NSLAVES =  NSLAVES_less
           ENDIF
           NSLAVES = max(NSLAVES_ref,NSLAVES)
         ENDIF
         NSLAVES = min (NSLAVES, NMB_OF_CAND)
         IF ( NSLAVES.GT.NSLAVES_ref) THEN
          NPIV = NFRONT - NCB
          IF ( K50.EQ.0 ) THEN
           WK_SLAVE = real( NPIV ) * real( NCB ) *
     &         ( 2.0E0 * real(NFRONT) - real(NPIV) )
     &         / real(NSLAVES)
           WK_MASTER = 0.66667E0 *
     &                 real(NPIV)*real(NPIV)*real(NPIV)+
     &                 real(NPIV)*real(NPIV)*real(NCB)
          ELSE
           WK_SLAVE = MUMPS_45(NCB,NFRONT,NPIV)
     &         / real(NSLAVES)
           WK_MASTER =  real(NPIV)*real(NPIV)*real(NPIV)/3.0E0
          ENDIF
          IF ( (WK_MASTER.GT.WK_SLAVE).AND.
     &     (WK_SLAVE.GT.1.0E0) ) THEN
           NSLAVES = 
     &         int( real(NSLAVES) * (WK_SLAVE/WK_MASTER))
           NSLAVES = max(NSLAVES_ref, NSLAVES)
          ENDIF
         ENDIF
      ELSE
       NSLAVES = NSLAVES_less
      ENDIF
      NSLAVES = min (NSLAVES, NCB)
      NSLAVES = min (NSLAVES, NMB_OF_CAND)
      MUMPS_12 = NSLAVES
      RETURN
      END FUNCTION MUMPS_12
      SUBROUTINE MUMPS_47( 
     &   KEEP,KEEP8, INODE, STEP, N, SLAVEF, 
     &   ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &   NASS, NCB, 
     &   NSLAVES, POSITION, ISLAVE, IPOSSLAVE )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: KEEP(500),INODE,N,SLAVEF 
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT( IN ) :: STEP(N), 
     &          ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER, INTENT( IN  ) :: NASS, NCB, 
     &                          NSLAVES, POSITION
      INTEGER, INTENT( OUT ) :: ISLAVE, IPOSSLAVE
      INTEGER BLSIZE, J, ISHIFT
      IF ((NSLAVES.LE.0).OR.(POSITION.LE.NASS)) THEN
       ISLAVE = 0
       IPOSSLAVE = POSITION
       RETURN
      ENDIF
      IF ( KEEP(48).EQ.0) THEN
       BLSIZE = NCB / NSLAVES
       ISLAVE    = min( NSLAVES,
     &               ( POSITION - NASS - 1 ) / BLSIZE + 1 )
       IPOSSLAVE = POSITION - NASS - ( ISLAVE - 1 ) * BLSIZE
      ELSEIF (KEEP(48).EQ.3) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSEIF (KEEP(48).EQ.4) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSEIF (KEEP(48).EQ.5) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSE
       WRITE(*,*) 'Error in MUMPS_47: undef strat'  
       CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_47
      INTEGER FUNCTION MUMPS_442( K821, K50, KMAX, NCB )
      IMPLICIT NONE
      INTEGER, INTENT( IN    )  :: KMAX, NCB, K50
      INTEGER(8), INTENT(IN) :: K821
      INTEGER KMIN, MINGRAN
      INTEGER(8) :: KMINSURF
      IF ( ( NCB .LE.0 ).OR. (KMAX.LE.0) ) THEN 
        MUMPS_442 = 1
        RETURN
      ENDIF
      IF (K50.EQ.0) THEN
       KMINSURF = 60000_8
#if defined(t3e) || defined(sgi)
       MINGRAN = 40 
#else
       MINGRAN = 50 
#endif
      ELSE
       KMINSURF = 30000_8
#if defined(t3e) || defined(sgi)
       MINGRAN = 10 
#else
       MINGRAN = 20 
#endif
      ENDIF
      IF (K821.GT.0_8) THEN
#if defined(t3e) || defined(sgi)
           KMIN = max(MINGRAN,KMAX/10)
#else
           KMIN = max(MINGRAN,KMAX/20)
#endif
      ELSE
           KMINSURF = max( abs(K821)/500_8, KMINSURF )
           KMIN     = max(
     &                     int( KMINSURF / int(max(NCB,1),8) ),
     &                     1
     &                   )
      ENDIF
      KMIN = min(KMIN,KMAX)
      KMIN = max(KMIN,1)  
      MUMPS_442 = KMIN
      RETURN
      END FUNCTION MUMPS_442
      INTEGER FUNCTION MUMPS_497( KEEP821, NCB )
      IMPLICIT NONE
      INTEGER,    intent( in    )  :: NCB
      INTEGER(8), intent( in    )  :: KEEP821
      INTEGER KMAX 
      IF ( NCB .LE.0 ) THEN 
        MUMPS_497 = 1
        RETURN
      ENDIF
      IF ( KEEP821.GT.0_8 ) THEN
       KMAX = int(KEEP821)
      ELSE
       KMAX =  -int(KEEP821/int(NCB,8))
      ENDIF
      KMAX = min (NCB, KMAX)
      MUMPS_497 = max ( KMAX, 1 )
      RETURN
      END FUNCTION MUMPS_497
      SUBROUTINE MUMPS_546( IS, DS )
      INTEGER IS, DS
#if defined(t3e)
      IS = 8
      DS = 16
#else
      IS = 4
      DS = 8
#endif
      END SUBROUTINE MUMPS_546
      SUBROUTINE MUMPS_SET_VERSION( VERSION_STR )
      IMPLICIT NONE
      CHARACTER(LEN=*) :: VERSION_STR
      CHARACTER(LEN=*) :: V;
      PARAMETER (V = "4.10.0" )
      IF ( len(V) .GT. 14 ) THEN
         WRITE(*,*) "Version string too long ( >14 characters )"
         CALL MUMPS_ABORT()
      END IF
      VERSION_STR = V
      RETURN
      END SUBROUTINE MUMPS_SET_VERSION
      SUBROUTINE MUMPS_420
     &                ( JOB, THRESH, NDENSE, 
     &                 N, IWLEN, PE, PFREE, LEN, IW, NV, 
     &                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     &                 PERM, COMPLEM_LIST, SIZE_COMPLEM_LIST, AGG6 ) 
      IMPLICIT NONE
      INTEGER JOB
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N), 
     &        ELEN(N), LAST(N), NCMPA, DEGREE(N), HEAD(N), NEXT(N), 
     &        W(N)
      LOGICAL AGG6
      INTEGER, intent(in) :: SIZE_COMPLEM_LIST
      INTEGER NDENSE(N)
      INTEGER, intent (in) :: COMPLEM_LIST(max(1,SIZE_COMPLEM_LIST))
      INTEGER PERM(N)
      INTEGER THRESH
      INTEGER THRESM, NDME, PERMeqN
      INTEGER NBD,NBED, NBDM, LASTD, NELME
      LOGICAL IDENSE
      INTEGER FDEG, ThresMin, ThresPrev, IBEGSchur, 
     &        ThresMinINIT
      LOGICAL SchurON
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER MAXINT_N
      INTEGER(8) HASH, HMOD 
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC max, min, mod
        IF (N.EQ.1) THEN
           ELEN(1) = 1
           LAST(1) = 1
           PE(1) = 0
           NV(1) = 1
           RETURN
        ENDIF
        IF ( SIZE_COMPLEM_LIST < 0 .OR. SIZE_COMPLEM_LIST > N ) THEN
          WRITE(*,*) "Internal MUMPS_420", SIZE_COMPLEM_LIST,N
          CALL MUMPS_ABORT()
        ENDIF
        IF (JOB.EQ.2) THEN
          SchurON = .FALSE.        
        ENDIF
        IF (JOB.NE.2) THEN
          SchurON   = (SIZE_COMPLEM_LIST > 0)
          IF ((JOB.EQ.1) .AND. (.NOT.SchurON) .AND. (N .GT. 0)) THEN
           WRITE(6,*) ' WARNING MUMPS_420 on Options '
          ENDIF
          IBEGSchur = N-SIZE_COMPLEM_LIST+1
          IF (THRESH.GT.N) THRESH = N
          IF (THRESH.LT.0) THRESH = 0
          IF ( SchurON )  THEN 
           DO I= 1, N
             IF ( PERM(I) .GE. IBEGSchur) THEN 
                 PERM(I) = N + 1
                IF (LEN(I) .EQ.0) THEN
                  PE(I) = 0
                ENDIF
             ENDIF
           ENDDO
          ENDIF
        ENDIF
        IF (SchurON) THEN
             THRESM    = N
             ThresMin  = N
             ThresPrev = N
        ELSE
             THRESM    = max(int(31*N/32),THRESH)
             THRESM    = max(THRESM,1)
             ThresMin  = max( 3*THRESM / 4, 1)
             ThresPrev = THRESM
        ENDIF
        ThresMinINIT = ThresMin/4
      IF (THRESM.GT.0) THEN
       IF ((THRESM.GT.N).OR.(THRESM.LT.2)) THEN 
          THRESM = N
       ENDIF
      ENDIF
      IF (JOB.EQ.2) THEN
      ENDIF
      PERMeqN = 0
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      NEL   = 0
      WFLG   = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA  = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0  
      MEM  = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        NDENSE(I)= 0
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
   10 CONTINUE
      IF (JOB.EQ.2) THEN
        DO I = 1,SIZE_COMPLEM_LIST
             X       = COMPLEM_LIST(I)
             ELEN(X) = -I       
             NV(X)   = LEN(X)+1 
             DMAX = max(DMAX, LEN(X))
        ENDDO
        NEL = NEL + SIZE_COMPLEM_LIST  
        DO I=1,N
          DEGREE (I) = LEN (I)  
        ENDDO
      ELSE
        DO I=1, N
          ELEN (I) = 0
          DEGREE (I) = LEN (I)
        ENDDO
      ENDIF
      DO 20 I = 1, N
        IF (ELEN(I).LT.0) CYCLE   
        DEG = DEGREE (I)
        IF (PERM(I).EQ.N) THEN
           PERMeqN = I
           PERM(I) = N-1
        ENDIF
        FDEG = PERM(I)
        IF ( (DEG .GT. 0).OR.(PERM(I).EQ.N+1) ) THEN
          IF ( (THRESM.GT.0) .AND.
     &         (FDEG .GT.THRESM) ) THEN
            NBD = NBD+1
            IF (FDEG.NE.N+1) THEN
             DEGREE(I) = DEGREE(I)+N+2
             DEG = N
             INEXT = HEAD (DEG)
             IF (INEXT .NE. 0) LAST (INEXT) = I
             NEXT (I) = INEXT
             HEAD (DEG) = I 
             LAST(I)  = 0
             IF (LASTD.EQ.0) LASTD=I
            ELSE
             NBED = NBED+1
             DEGREE(I) = N+1
             DEG = N
             IF (LASTD.EQ.0) THEN
               LASTD     = I 
               HEAD(DEG) = I
               NEXT(I)   = 0 
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
            ENDIF
          ELSE
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            HEAD (FDEG) = I
          ENDIF
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
          IF ((NBD.EQ.0).AND.(THRESM.GT.0)) THRESM = N
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF ( (DEG.NE.N) .AND.
     &    (DEG.GT.THRESM+1) .AND. (NBD.GT.0) ) THEN
           MINDEG = N
           GOTO 30
        ENDIF
        IF (DEGREE(ME).LE.N)  THEN
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          MINDEG = 1
          NBDM = max(NBDM,NBD)
          IF (DEGREE(ME).GT.N+1) THEN
            IF (WFLG .GT. MAXINT_N) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + LEN(ME) -1
            LN       = P1
            ELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0) THEN
                X = E
  53            X = -PE(X)
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(LN) = IW(ELN)
               IW(ELN) = E
               LN  = LN+1
               ELN = ELN + 1
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(LN)=E
               LN = LN+1
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = LN-P1
            ELEN(ME) = ELN- P1
            NDME = NDENSE(ME)+NV(ME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            DEGREE(ME) = NDENSE(ME)
            DEG = PERM(ME)
            MINDEG = min(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(N+1) ) GOTO 51
            ENDIF
            HEAD (N) = ME
            IF (THRESM.LT.N) THEN
             ThresMin  = max(THRESM+ThresMin,ThresPrev+ThresMin/2+1)
             ThresMin  = min(ThresMin, N)
             ThresPrev = ThresPrev+(N-ThresPrev)/2+ThresMinINIT
             THRESM    = max(
     &         THRESM + int(sqrt(dble(ThresMin)))+ ThresMinINIT ,
     &         ThresPrev)
             THRESM    = min(THRESM,N) 
             ThresMin  = min(THRESM, ThresMin)
             ThresPrev = THRESM
            ENDIF
            NBD    = NBED
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.N+1) THEN
           IF (NBD.NE.NBED) THEN
          write(6,*) ' ERROR in MUMPS_420 ',
     &                ' quasi dense rows remains'
            CALL MUMPS_ABORT()
           ENDIF
           IF (JOB.EQ.1) THEN
            DO I = 1,SIZE_COMPLEM_LIST
             X       = COMPLEM_LIST(I)
             ELEN(X) = -(N-SIZE_COMPLEM_LIST+I)   
             NV(X)   = 1                
             PE(X)   = 0                
            ENDDO
            GOTO 265
           ENDIF
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -COMPLEM_LIST(1)
            ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
   59      CONTINUE
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0
           IF (NEL.NE.N) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL MUMPS_ABORT()
           ENDIF
           IF (ME.NE. COMPLEM_LIST(1)) THEN
             DO I=1, SIZE_COMPLEM_LIST
               PE(COMPLEM_LIST(I)) = -COMPLEM_LIST(1)
             ENDDO
             PE(COMPLEM_LIST(1)) = 0
             NV( COMPLEM_LIST(1))= NV(ME)
             NV(ME)               = 0
             ELEN( COMPLEM_LIST(1)) = ELEN(ME)
             ELEN(ME)             = 0
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.N) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (PERM(I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (PERM(I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (.NOT. AGG6 .AND. DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (AGG6 .AND. (DEXT .EQ. 0) .AND.
     &            ((NDENSE(ME).EQ.NBD).OR.(NDENSE(E).EQ.0))) THEN
                PE (E) = -ME
                W (E)  = 0
             ELSE IF (AGG6 .AND. DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + int(E,kind=8)
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              IF (DEGREE(J).LE.N) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (((ELEN(I).EQ.1).AND.(P3.EQ.PN))
     &     .OR.
     &         (AGG6.AND.(DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD))
     &       )
     &    THEN
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE(I) = min (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.N) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
             X = I 
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                IF (PERM(J).GT.PERM(X)) THEN
                  PE (J) = -X
                  NV (X) = NV (X) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                ELSE
                  PE (X) = -J
                  NV (J) = NV (X) + NV (J)
                  NV (X) = 0
                  ELEN (X) = 0
                  X = J
                ENDIF
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        IF ( (THRESM .GT. 0).AND.(THRESM.LT.N) ) THEN 
          THRESM = max(ThresMin, THRESM-NVPIV)
        ENDIF
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.N) THEN
            DEG = min (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
            IF (THRESM.GT.0) THEN
             IF (PERM(I) .GT. THRESM) THEN
               IDENSE = .TRUE.
               DEGREE(I) = DEGREE(I)+N+2
             ENDIF
             IF (IDENSE) THEN
               P1 = PE(I)
               P2 = P1 + ELEN(I) - 1
               IF (P2.GE.P1) THEN
               DO 264 PJ=P1,P2
                 E= IW(PJ)
                 NDENSE (E) = NDENSE(E) + NVI
 264           CONTINUE
               ENDIF
               NBD = NBD+NVI
               FDEG = N
               DEG = N
               INEXT = HEAD(DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
               LAST(I)    = 0
               IF (LASTD.EQ.0) LASTD=I
             ENDIF
            ENDIF
            IF (.NOT.IDENSE) THEN
            FDEG = PERM(I)
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (FDEG) = I
            ENDIF
            MINDEG = min (MINDEG, FDEG)
            ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
  265 CONTINUE
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      IF (.NOT.SchurON) THEN
        IF (PERMeqN.GT.0) PERM(PERMeqN) = N
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_420
      SUBROUTINE MUMPS_209( N, FRERE, FILS, NFSIZ, THEROOT )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: N
      INTEGER, intent( in    )  :: NFSIZ( N )
      INTEGER, intent( inout )  :: FRERE( N ), FILS( N )
      INTEGER, intent( out   )  :: THEROOT
      INTEGER INODE, IROOT, IFILS, IN, IROOTLAST, SIZE
      IROOT = -9999
      SIZE  = 0
      DO INODE = 1, N
        IF ( FRERE( INODE ) .EQ. 0 )  THEN
          IF ( NFSIZ( INODE ) .GT. SIZE ) THEN
            SIZE  = NFSIZ( INODE )
            IROOT = INODE
          END IF
        ENDIF
      END DO
      IN = IROOT
      DO WHILE ( FILS( IN ) .GT. 0 )
        IN = FILS( IN )
      END DO
      IROOTLAST = IN
      IFILS     = - FILS ( IN )
      DO INODE = 1, N
        IF ( FRERE( INODE ) .eq. 0 .and. INODE .ne. IROOT ) THEN
          IF ( IFILS .eq. 0 ) THEN
            FILS( IROOTLAST ) = - INODE
            FRERE( INODE )    = -IROOT
            IFILS             = INODE
          ELSE
            FRERE( INODE ) = -FILS( IROOTLAST )
            FILS( IROOTLAST ) = - INODE
          END IF
        END IF
      END DO
      THEROOT = IROOT
      RETURN
      END SUBROUTINE MUMPS_209
      INTEGER FUNCTION MUMPS_330(PROCINFO_INODE, SLAVEF)
      IMPLICIT NONE
      INTEGER SLAVEF 
      INTEGER PROCINFO_INODE, TPN
      IF (PROCINFO_INODE <= SLAVEF ) THEN
        MUMPS_330 = 1
      ELSE
        TPN = (PROCINFO_INODE-1+2*SLAVEF)/SLAVEF - 1
        IF ( TPN .LT. 1 ) TPN = 1
        IF (TPN.EQ.4.OR.TPN.EQ.5.OR.TPN.EQ.6) TPN = 2
        MUMPS_330 = TPN
      END IF
      RETURN 
      END FUNCTION MUMPS_330
      INTEGER FUNCTION MUMPS_275(PROCINFO_INODE, SLAVEF)
      IMPLICIT NONE
      INTEGER SLAVEF 
      INTEGER PROCINFO_INODE
      IF (SLAVEF == 1) THEN
        MUMPS_275 = 0
      ELSE
        MUMPS_275=mod(2*SLAVEF+PROCINFO_INODE-1,SLAVEF)
      END IF
      RETURN
      END FUNCTION MUMPS_275
      INTEGER FUNCTION MUMPS_810 (PROCINFO_INODE, SLAVEF)
      IMPLICIT NONE
      INTEGER, intent(in) ::  SLAVEF 
      INTEGER PROCINFO_INODE, TPN
      IF (PROCINFO_INODE <= SLAVEF ) THEN
         MUMPS_810 = 1
      ELSE
        TPN = (PROCINFO_INODE-1+2*SLAVEF)/SLAVEF - 1
        IF ( TPN .LT. 1 ) TPN = 1
         MUMPS_810 = TPN
      ENDIF
      RETURN
      END FUNCTION MUMPS_810
      LOGICAL FUNCTION MUMPS_283( PROCINFO_INODE, SLAVEF )
      IMPLICIT NONE
      INTEGER SLAVEF
      INTEGER TPN, PROCINFO_INODE
      TPN = (PROCINFO_INODE-1+2*SLAVEF)/SLAVEF - 1
      MUMPS_283 = ( TPN .eq. 0 )
      RETURN
      END FUNCTION MUMPS_283
      LOGICAL FUNCTION MUMPS_167( PROCINFO_INODE, SLAVEF )
      IMPLICIT NONE
      INTEGER SLAVEF
      INTEGER TPN, PROCINFO_INODE
      TPN = (PROCINFO_INODE-1+SLAVEF+SLAVEF)/SLAVEF - 1
      MUMPS_167 = ( TPN .eq. -1 )
      RETURN 
      END FUNCTION MUMPS_167
      LOGICAL FUNCTION MUMPS_170
     &        ( PROCINFO_INODE, SLAVEF )
      IMPLICIT NONE
      INTEGER SLAVEF
      INTEGER TPN, PROCINFO_INODE
      TPN = (PROCINFO_INODE-1+SLAVEF+SLAVEF)/SLAVEF - 1
      MUMPS_170 =
     &           ( TPN .eq. -1 .OR. TPN .eq. 0 )
      RETURN
      END FUNCTION MUMPS_170
      LOGICAL FUNCTION MUMPS_358( MYID, SLAVEF, INODE,
     &                 NMB_PAR2, ISTEP_TO_INIV2 , K71, STEP, N, 
     &                 CANDIDATES, KEEP24 )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, INODE, NMB_PAR2, KEEP24, I
      INTEGER K71, N
      INTEGER ISTEP_TO_INIV2 ( K71 ), STEP ( N )
      INTEGER CANDIDATES(SLAVEF+1, max(NMB_PAR2,1))
      INTEGER NCAND, POSINODE
      MUMPS_358 = .FALSE.
      IF (KEEP24 .eq. 0) RETURN
      POSINODE = ISTEP_TO_INIV2 ( STEP (INODE) )
      NCAND = CANDIDATES( SLAVEF+1, POSINODE )
      DO I = 1, NCAND
        IF (MYID .EQ. CANDIDATES( I, POSINODE ))
     &     MUMPS_358 = .TRUE.
      END DO
      RETURN
      END FUNCTION MUMPS_358
      SUBROUTINE MUMPS_291(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()
      RETURN
      END SUBROUTINE MUMPS_291
      SUBROUTINE MUMPS_292(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()-T
      RETURN
      END SUBROUTINE MUMPS_292
      SUBROUTINE MUMPS_558( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      DOUBLE PRECISION VAL( N )
      INTEGER I, ISWAP
      DOUBLE PRECISION SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_558
#if defined (PESSL)
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     &                     LLD, INFO )
      INTEGER            ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     &                   LLD_, MB_, M_, NB_, N_, RSRC_
# if defined(DESC8)
      PARAMETER          ( DLEN_ = 8, DTYPE_ = 1,
     &                     CTXT_ = 7, M_ = 1, N_ = 2, MB_ = 3, NB_ = 4,
     &                     RSRC_ = 5, CSRC_ = 6, LLD_ = 8 )
# else
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     &                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     &                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
# endif
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           blacs_gridinfo, PXERBLA
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      INTRINSIC          max, min
      CALL blacs_gridinfo( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( MB.LT.1 ) THEN
         INFO = -4
      ELSE IF( NB.LT.1 ) THEN
         INFO = -5
      ELSE IF( IRSRC.LT.0 .OR. IRSRC.GE.NPROW ) THEN
         INFO = -6
      ELSE IF( ICSRC.LT.0 .OR. ICSRC.GE.NPCOL ) THEN
         INFO = -7
      ELSE IF( NPROW.EQ.-1 ) THEN
         INFO = -8
      ELSE IF( LLD.LT.max( 1, numroc( M, MB, MYROW, IRSRC,
     &                                NPROW ) ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 )
     &   CALL PXERBLA( ICTXT, 'DESCINIT', -INFO )
# ifndef DESC8
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D
# endif
      DESC( M_ )  = max( 0, M )
      DESC( N_ )  = max( 0, N )
      DESC( MB_ ) = max( 1, MB )
      DESC( NB_ ) = max( 1, NB )
      DESC( RSRC_ ) = max( 0, min( IRSRC, NPROW-1 ) )
      DESC( CSRC_ ) = max( 0, min( ICSRC, NPCOL-1 ) )
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ )  = max( LLD, max( 1, numroc( DESC( M_ ), DESC( MB_ ),
     &                              MYROW, DESC( RSRC_ ), NPROW ) ) )
      RETURN
      END SUBROUTINE DESCINIT
      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
      INTEGER            ICTXT, INFO
      CHARACTER*(*)      SRNAME
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           blacs_gridinfo
      CALL blacs_gridinfo( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     &        ' parameter number', I4, ' had an illegal value' )
      END SUBROUTINE PXERBLA
#endif
      SUBROUTINE MUMPS_243(MYID, COMM, INFO, INFOG, IRANK)
      IMPLICIT NONE
      INTEGER MYID, COMM, IRANK, INFO, INFOG(2)
      INCLUDE 'mpif.h'
      INTEGER IERR_MPI, MASTER
      INTEGER TEMP1(2), TEMP2(2)
      PARAMETER( MASTER = 0 )
      CALL MPI_REDUCE( INFO, INFOG(1), 1, MPI_INTEGER,
     &                 MPI_MAX, MASTER, COMM, IERR_MPI )
      CALL MPI_REDUCE( INFO, INFOG(2), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, COMM, IERR_MPI )
      TEMP1(1) = INFO
      TEMP1(2) = MYID
      CALL MPI_REDUCE( TEMP1, TEMP2, 1, MPI_2INTEGER,
     &                 MPI_MAXLOC, MASTER, COMM, IERR_MPI )
      IF ( MYID.eq. MASTER ) THEN
        IF ( INFOG(1) .ne. TEMP2(1) ) THEN
          write(*,*) 'Error in MUMPS_243'
          CALL MUMPS_ABORT()
        END IF
        IRANK    = TEMP2(2)
      ELSE
        IRANK    = -1
      END IF
      RETURN
      END SUBROUTINE MUMPS_243
      SUBROUTINE MUMPS_362(N, LEAF, NBROOT, NROOT_LOC,
     &           MYID_NODES,
     &           SLAVEF, NA, LNA, KEEP,KEEP8, STEP,
     &           PROCNODE_STEPS, IPOOL, LPOOL)
      IMPLICIT NONE
      INTEGER N, LEAF, NROOT_LOC, NBROOT, MYID_NODES,
     &        SLAVEF, LPOOL, LNA
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28)), NA(LNA),
     &        IPOOL(LPOOL)
      INTEGER NBLEAF, INODE, I
      INTEGER MUMPS_275
      EXTERNAL MUMPS_275
      NBLEAF = NA(1)
      NBROOT = NA(2)
      LEAF = 1
      DO I = 1, NBLEAF
        INODE = NA(I+2)
        IF (MUMPS_275(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
     &   .EQ.MYID_NODES) THEN
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
          ENDIF
      ENDDO
      NROOT_LOC = 0
      DO I = 1, NBROOT
        INODE = NA(I+2+NBLEAF)
        IF (MUMPS_275(PROCNODE_STEPS(STEP(INODE)),
     &    SLAVEF).EQ.MYID_NODES) THEN
            NROOT_LOC = NROOT_LOC + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_362
      LOGICAL FUNCTION MUMPS_438(TAB1,TAB2,LEN1,LEN2)
      IMPLICIT NONE
      INTEGER LEN1 , LEN2 ,I
      INTEGER TAB1(LEN1)
      INTEGER TAB2(LEN2)
      MUMPS_438=.FALSE.
      IF(LEN1 .NE. LEN2) THEN
         RETURN
      ENDIF
      DO I=1 , LEN1
         IF(TAB1(I) .NE. TAB2(I)) THEN
            RETURN
         ENDIF
      ENDDO
      MUMPS_438=.TRUE.
      RETURN
      END FUNCTION MUMPS_438
      SUBROUTINE MUMPS_463( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_463
      SUBROUTINE MUMPS_466( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .LT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_466
      SUBROUTINE MUMPS_ABORT()
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, IERRCODE
      IERRCODE = -99
      CALL MPI_ABORT(MPI_COMM_WORLD, IERRCODE, IERR)
      RETURN
      END SUBROUTINE MUMPS_ABORT
      SUBROUTINE MUMPS_633(KEEP12,ICNTL14,
     &     KEEP50,KEEP54,ICNTL6,ICNTL8)
      IMPLICIT NONE
      INTEGER, intent(out)::KEEP12
      INTEGER, intent(in)::ICNTL14,KEEP50,KEEP54,ICNTL6,ICNTL8
      KEEP12 = ICNTL14 
      IF(ICNTL6.EQ.0 .AND. ICNTL8.EQ.0) RETURN
      IF ( (KEEP54.NE.0).AND. (KEEP50.NE.1)
     &     .AND. (KEEP12 .GT. 0) ) KEEP12= KEEP12+5
      RETURN
      END SUBROUTINE MUMPS_633
      SUBROUTINE MUMPS_749( I8_VALUE, ROOT, MYID, COMM, IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ROOT, MYID, COMM, IERR
      INTEGER(8) :: I8_VALUE
      DOUBLE PRECISION :: DBLE_VALUE
      IF (MYID .EQ. ROOT) THEN
        DBLE_VALUE = dble(I8_VALUE)
      ENDIF
      CALL MPI_BCAST( DBLE_VALUE, 1, MPI_DOUBLE_PRECISION,
     &                ROOT,  COMM, IERR )
      I8_VALUE = int( DBLE_VALUE,8)
      RETURN
      END SUBROUTINE MUMPS_749
      SUBROUTINE MUMPS_646( IN, OUT, MPI_OP, ROOT, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ROOT, COMM, MPI_OP
      INTEGER(8) IN, OUT
      INTEGER IERR
      DOUBLE PRECISION DIN, DOUT
      DIN =dble(IN)
      DOUT=0.0D0
      CALL MPI_REDUCE(DIN, DOUT, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_OP, ROOT, COMM, IERR)
      OUT=int(DOUT,kind=8)
      RETURN
      END SUBROUTINE MUMPS_646
      SUBROUTINE MUMPS_736( IN, OUT, MPI_OP, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, MPI_OP
      INTEGER(8) IN, OUT
      INTEGER IERR
      DOUBLE PRECISION DIN, DOUT
      DIN =dble(IN)
      DOUT=0.0D0
      CALL MPI_ALLREDUCE(DIN, DOUT, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_OP, COMM, IERR)
      OUT=int(DOUT,kind=8)
      RETURN
      END SUBROUTINE MUMPS_736
      SUBROUTINE MUMPS_754(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      INTEGER, POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      INTEGER, POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER                    :: FMTA*60, FMTD*60
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL = (/ERRCODE, MINSIZE/)
      ELSE
         ERRTPL = (/-13, MINSIZE/)
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_754
      SUBROUTINE MUMPS_750(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(kind(1.E0)), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(kind(1.E0)), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER                    :: FMTA*60, FMTD*60
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL = (/ERRCODE, MINSIZE/)
      ELSE
         ERRTPL = (/-13, MINSIZE/)
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_750
      SUBROUTINE MUMPS_752(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(kind(1.D0)), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(kind(1.D0)), POINTER    :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER                    :: FMTA*60, FMTD*60
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL = (/ERRCODE, MINSIZE/)
      ELSE
         ERRTPL = (/-13, MINSIZE/)
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_752
      SUBROUTINE MUMPS_751(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(kind((1.E0,1.E0))), POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(kind((1.E0,1.E0))), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER                    :: FMTA*60, FMTD*60
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL = (/ERRCODE, MINSIZE/)
      ELSE
         ERRTPL = (/-13, MINSIZE/)
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_751
      SUBROUTINE MUMPS_753(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(kind((1.D0,1.D0))), POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(kind((1.D0,1.D0))), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER                    :: FMTA*60, FMTD*60
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL = (/ERRCODE, MINSIZE/)
      ELSE
         ERRTPL = (/-13, MINSIZE/)
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_753
      SUBROUTINE MUMPS_735(I8, I4)
      IMPLICIT NONE
      INTEGER   , INTENT(OUT) :: I4
      INTEGER(8), INTENT(IN)  :: I8
      IF ( I8 .GT. int(huge(I4),8) ) THEN
        I4 = -int(I8/1000000_8,kind(I4))
      ELSE
        I4 = int(I8,kind(I4))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_735
      SUBROUTINE MUMPS_ABORT_ON_OVERFLOW(I8, STRING)
      IMPLICIT NONE
      INTEGER(8), INTENT(IN) :: I8
      CHARACTER(*), INTENT(IN) :: STRING
      INTEGER I4
      IF ( I8 .GT. int(huge(I4),8)) THEN
        WRITE(*,*) STRING
        CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_ABORT_ON_OVERFLOW
      SUBROUTINE MUMPS_731( SIZE8, IERROR  )
      INTEGER(8), INTENT(IN) :: SIZE8
      INTEGER, INTENT(OUT) :: IERROR
      CALL MUMPS_735(SIZE8, IERROR)
      RETURN
      END SUBROUTINE MUMPS_731
      SUBROUTINE MUMPS_730(I8, INT_ARRAY)
      IMPLICIT NONE
      INTEGER(8), intent(in)  :: I8
      INTEGER,    intent(out) :: INT_ARRAY(2)
      INTEGER(kind(0_4)) :: I32
      INTEGER(8) :: IDIV, IPAR
      PARAMETER (IPAR=int(huge(I32),8))
      PARAMETER (IDIV=IPAR+1_8)
      IF ( I8 .LT. IDIV ) THEN
        INT_ARRAY(1) = 0
        INT_ARRAY(2) = int(I8)
      ELSE
        INT_ARRAY(1) = int(I8 / IDIV)
        INT_ARRAY(2) = int(mod(I8,IDIV))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_730
      SUBROUTINE MUMPS_729(I8, INT_ARRAY)
      IMPLICIT NONE
      INTEGER(8), intent(out)  :: I8
      INTEGER,    intent(in)  :: INT_ARRAY(2)
      INTEGER(kind(0_4)) :: I32
      INTEGER(8) :: IDIV, IPAR
      PARAMETER (IPAR=int(huge(I32),8))
      PARAMETER (IDIV=IPAR+1_8)
      IF ( INT_ARRAY(1) .EQ. 0 ) THEN
        I8=int(INT_ARRAY(2),8)
      ELSE
        I8=int(INT_ARRAY(1),8)*IDIV+int(INT_ARRAY(2),8)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_729
      SUBROUTINE MUMPS_723( INT_ARRAY, I8 )
      IMPLICIT NONE
      INTEGER(8), intent(in) :: I8
      INTEGER, intent(inout) :: INT_ARRAY(2)
      INTEGER(8) :: I8TMP
      CALL MUMPS_729(I8TMP, INT_ARRAY)
      I8TMP = I8TMP + I8
      CALL MUMPS_730(I8TMP, INT_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_723
      SUBROUTINE MUMPS_724( INT_ARRAY, I8 )
      IMPLICIT NONE
      INTEGER(8), intent(in) :: I8
      INTEGER, intent(inout) :: INT_ARRAY(2)
      INTEGER(8) :: I8TMP
      CALL MUMPS_729(I8TMP, INT_ARRAY)
      I8TMP = I8TMP - I8
      CALL MUMPS_730(I8TMP, INT_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_724
      FUNCTION MUMPS_815(WHICH)
      LOGICAL :: MUMPS_815
      CHARACTER :: WHICH*(*)
      LOGICAL :: PTSCOTCH=.FALSE., PARMETIS=.FALSE.
#if defined(ptscotch)
      PTSCOTCH = .TRUE.
#endif
#if defined(parmetis)
      PARMETIS = .TRUE.
#endif
      SELECT CASE(WHICH)
      CASE('ptscotch','PTSCOTCH')
         MUMPS_815 = PTSCOTCH
      CASE('parmetis','PARMETIS')
         MUMPS_815 = PARMETIS
      CASE('both','BOTH')
         MUMPS_815 = PTSCOTCH .AND. PARMETIS
      CASE('any','ANY')
         MUMPS_815 = PTSCOTCH .OR. PARMETIS
      CASE default
         write(*,'("Invalid input in MUMPS_815")')
      END SELECT
      RETURN
      END FUNCTION MUMPS_815
