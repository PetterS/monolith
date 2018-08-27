!
!  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
!
!
!  This version of MUMPS is provided to you free of charge. It is public
!  domain, based on public domain software developed during the Esprit IV
!  European project PARASOL (1996-1999). Since this first public domain
!  version in 1999, research and developments have been supported by the
!  following institutions: CERFACS, CNRS, ENS Lyon, INPT(ENSEEIHT)-IRIT,
!  INRIA, and University of Bordeaux.
!
!  The MUMPS team at the moment of releasing this version includes
!  Patrick Amestoy, Maurice Bremond, Alfredo Buttari, Abdou Guermouche,
!  Guillaume Joslin, Jean-Yves L'Excellent, Francois-Henry Rouet, Bora
!  Ucar and Clement Weisbecker.
!
!  We are also grateful to Emmanuel Agullo, Caroline Bousquet, Indranil
!  Chowdhury, Philippe Combes, Christophe Daniel, Iain Duff, Vincent Espirat,
!  Aurelia Fevre, Jacko Koster, Stephane Pralet, Chiara Puglisi, Gregoire
!  Richard, Tzvetomila Slavova, Miroslav Tuma and Christophe Voemel who
!  have been contributing to this project.
!
!  Up-to-date copies of the MUMPS package can be obtained
!  from the Web pages:
!  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
!
!
!   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
!   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
!
!
!  User documentation of any code that uses this software can
!  include this complete notice. You can acknowledge (using
!  references [1] and [2]) the contribution of this package
!  in any scientific publication dependent upon the use of the
!  package. You shall use reasonable endeavours to notify
!  the authors of the package of this publication.
!
!   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
!   A fully asynchronous multifrontal solver using distributed dynamic
!   scheduling, SIAM Journal of Matrix Analysis and Applications,
!   Vol 23, No 1, pp 15-41 (2001).
!
!   [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
!   S. Pralet, Hybrid scheduling for the parallel solution of linear
!   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
!
!
!      Dummy mpif.h file including symbols used by MUMPS.
!
      INTEGER MPI_2DOUBLE_PRECISION
      INTEGER MPI_2INTEGER
      INTEGER MPI_2REAL
      INTEGER MPI_ANY_SOURCE
      INTEGER MPI_ANY_TAG
      INTEGER MPI_BYTE
      INTEGER MPI_CHARACTER
      INTEGER MPI_COMM_NULL
      INTEGER MPI_COMM_WORLD
      INTEGER MPI_COMPLEX
      INTEGER MPI_DOUBLE_COMPLEX
      INTEGER MPI_DOUBLE_PRECISION
      INTEGER MPI_INTEGER
      INTEGER MPI_LOGICAL
      INTEGER MPI_MAX
      INTEGER MPI_MAX_PROCESSOR_NAME
      INTEGER MPI_MAXLOC
      INTEGER MPI_MIN
      INTEGER MPI_MINLOC
      INTEGER MPI_PACKED
      INTEGER MPI_PROD
      INTEGER MPI_REAL
      INTEGER MPI_REPLACE
      INTEGER MPI_REQUEST_NULL
      INTEGER MPI_SOURCE
      INTEGER MPI_STATUS_SIZE
      INTEGER MPI_SUM
      INTEGER MPI_TAG
      INTEGER MPI_UNDEFINED
      INTEGER MPI_WTIME_IS_GLOBAL
      INTEGER MPI_LOR
      INTEGER MPI_LAND
      INTEGER MPI_INTEGER8
      INTEGER MPI_REAL8
      INTEGER MPI_BSEND_OVERHEAD
      PARAMETER (MPI_2DOUBLE_PRECISION=1)
      PARAMETER (MPI_2INTEGER=2)
      PARAMETER (MPI_2REAL=3)
      PARAMETER (MPI_ANY_SOURCE=4)
      PARAMETER (MPI_ANY_TAG=5)
      PARAMETER (MPI_BYTE=6)
      PARAMETER (MPI_CHARACTER=7)
      PARAMETER (MPI_COMM_NULL=8)
      PARAMETER (MPI_COMM_WORLD=9)
      PARAMETER (MPI_COMPLEX=10)
      PARAMETER (MPI_DOUBLE_COMPLEX=11)
      PARAMETER (MPI_DOUBLE_PRECISION=12)
      PARAMETER (MPI_INTEGER=13)
      PARAMETER (MPI_LOGICAL=14)
      PARAMETER (MPI_MAX=15)
      PARAMETER (MPI_MAX_PROCESSOR_NAME=31)
      PARAMETER (MPI_MAXLOC=16)
      PARAMETER (MPI_MIN=17)
      PARAMETER (MPI_MINLOC=18)
      PARAMETER (MPI_PACKED=19)
      PARAMETER (MPI_PROD=20)
      PARAMETER (MPI_REAL=21)
      PARAMETER (MPI_REPLACE=22)
      PARAMETER (MPI_REQUEST_NULL=23)
      PARAMETER (MPI_SOURCE=1)
      PARAMETER (MPI_STATUS_SIZE=2)
      PARAMETER (MPI_SUM=26)
      PARAMETER (MPI_TAG=2)
      PARAMETER (MPI_UNDEFINED=28)
      PARAMETER (MPI_WTIME_IS_GLOBAL=30)
      PARAMETER (MPI_LOR=31)
      PARAMETER (MPI_LAND=32)
      PARAMETER (MPI_INTEGER8=33)
      PARAMETER (MPI_REAL8=34)

      PARAMETER (MPI_BSEND_OVERHEAD=0)
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
