/*
 *
 *  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
 *
 *
 *  This version of MUMPS is provided to you free of charge. It is public
 *  domain, based on public domain software developed during the Esprit IV
 *  European project PARASOL (1996-1999). Since this first public domain
 *  version in 1999, research and developments have been supported by the
 *  following institutions: CERFACS, CNRS, ENS Lyon, INPT(ENSEEIHT)-IRIT,
 *  INRIA, and University of Bordeaux.
 *
 *  The MUMPS team at the moment of releasing this version includes
 *  Patrick Amestoy, Maurice Bremond, Alfredo Buttari, Abdou Guermouche,
 *  Guillaume Joslin, Jean-Yves L'Excellent, Francois-Henry Rouet, Bora
 *  Ucar and Clement Weisbecker.
 *
 *  We are also grateful to Emmanuel Agullo, Caroline Bousquet, Indranil
 *  Chowdhury, Philippe Combes, Christophe Daniel, Iain Duff, Vincent Espirat,
 *  Aurelia Fevre, Jacko Koster, Stephane Pralet, Chiara Puglisi, Gregoire
 *  Richard, Tzvetomila Slavova, Miroslav Tuma and Christophe Voemel who
 *  have been contributing to this project.
 *
 *  Up-to-date copies of the MUMPS package can be obtained
 *  from the Web pages:
 *  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
 *
 *
 *   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 *   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
 *
 *
 *  User documentation of any code that uses this software can
 *  include this complete notice. You can acknowledge (using
 *  references [1] and [2]) the contribution of this package
 *  in any scientific publication dependent upon the use of the
 *  package. You shall use reasonable endeavours to notify
 *  the authors of the package of this publication.
 *
 *   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
 *   A fully asynchronous multifrontal solver using distributed dynamic
 *   scheduling, SIAM Journal of Matrix Analysis and Applications,
 *   Vol 23, No 1, pp 15-41 (2001).
 *
 *   [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
 *   S. Pralet, Hybrid scheduling for the parallel solution of linear
 *   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
 *
 */
#ifndef MUMPS_ORDERINGS_H
#define MUMPS_ORDERINGS_H
#include "mumps_common.h"
#if defined(pord)
#include <space.h>
int mumps_pord( int, int, int *, int *, int * );
#define MUMPS_PORDF \
    F_SYMBOL(pordf,PORDF)
void MUMPS_CALL
MUMPS_PORDF( int *nvtx, int *nedges,
             int *xadj, int *adjncy,
             int *nv, int *ncmpa );
int mumps_pord_wnd( int, int, int *, int *, int *, int * );
#define MUMPS_PORDF_WND          \
    F_SYMBOL(pordf_wnd,PORDF_WND)
void MUMPS_CALL
MUMPS_PORDF_WND( int *nvtx, int *nedges,
                 int *xadj, int *adjncy,
                 int *nv, int *ncmpa, int *totw );
#endif /*PORD*/
#if defined(scotch) || defined(ptscotch)
int esmumps( const int n, const int iwlen, int * const pe, const int pfree,
             int * const len, int * const iw, int * const nv, int * const elen,
             int * const last);
#define MUMPS_SCOTCH        \
    F_SYMBOL(scotch,SCOTCH)
void MUMPS_CALL
MUMPS_SCOTCH( const int * const  n,
              const int * const  iwlen,
              int * const        petab,
              const int * const  pfree,
              int * const        lentab,
              int * const        iwtab,
              int * const        nvtab,
              int * const        elentab,
              int * const        lasttab,
              int * const        ncmpa );
#endif /*scotch or ptscotch*/
#if defined(ptscotch)
#include "mpi.h"
#include <stdio.h>
#include "ptscotch.h"
int mumps_dgraphinit( SCOTCH_Dgraph *, MPI_Fint *, MPI_Fint *);
#define MUMPS_DGRAPHINIT \
  F_SYMBOL(dgraphinit,DGRAPHINIT)
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr);
#endif /*ptscotch*/
#if defined(parmetis)
#include "mpi.h"
#include "parmetis.h"
void mumps_parmetis(int *first,      int *vertloctab, 
                   int *edgeloctab, int *numflag, 
                   int *options,    int *order, 
                   int *sizes,      int *comm);
#define MUMPS_PARMETIS \
  F_SYMBOL(parmetis,PARMETIS)
void MUMPS_CALL
MUMPS_PARMETIS(int *first,      int *vertloctab, 
               int *edgeloctab, int *numflag, 
               int *options,    int *order, 
               int *sizes,      int *comm);
#endif /*PARMETIS*/
#endif /* MUMPS_ORDERINGS_H */
