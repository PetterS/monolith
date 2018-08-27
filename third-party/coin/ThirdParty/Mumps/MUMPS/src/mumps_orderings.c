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
/*
 * This file contains interfaces to external ordering packages.
 * At the moment, PORD (J. Schulze) and SCOTCH are interfaced.
 */
#include "mumps_orderings.h"
#if defined(pord)
/* Interface to PORD */
/*int mumps_pord( int, int, int *, int *, int * );
#define MUMPS_PORDF   \
F_SYMBOL(pordf,PORDF)*/
void MUMPS_CALL
MUMPS_PORDF( int *nvtx, int *nedges,
             int *xadj, int *adjncy,
             int *nv, int *ncmpa )
{
    *ncmpa = mumps_pord( *nvtx, *nedges, xadj, adjncy, nv );
}
/* Interface to PORD with weighted graph*/
/*int mumps_pord_wnd( int, int, int *, int *, int *, int * );
#define MUMPS_PORDF_WND           \
    F_SYMBOL(pordf_wnd,PORDF_WND)*/
void MUMPS_CALL
MUMPS_PORDF_WND( int *nvtx, int *nedges,
                 int *xadj, int *adjncy,
                 int *nv, int *ncmpa, int *totw )
{
    *ncmpa = mumps_pord_wnd( *nvtx, *nedges, xadj, adjncy, nv, totw );
}
/************************************************************
 mumps_pord is used in ana_aux.F
        permutation and inverse permutation not set in output,
        but may be printed in default file: "perm_pord" and "iperm_pord",
        if associated part uncommneted.
        But, if uncommetnted a bug occurs in psl_ma41_analysi.F
******************************************************************/
/*********************************************************/
int mumps_pord
(
   int nvtx,
   int nedges,
   int *xadj_pe,
   int *adjncy,
   int *nv
)
{
/**********************************
Argument Comments:
input:
-----
- nvtx          : dimension of the Problem (N)
- nedges        : number of entries (NZ)
- adjncy        : non-zeros entries (IW input)
input/output:
-------------
- xadj_pe       : pointer through beginning of column non-zeros entries (PTRAR)
- on exit, "father array" (PE)
ouput:
------
- nv            : "nfront array" (NV)
*************************************/
  graph_t    *G;
  elimtree_t *T;
  timings_t  cpus[12];
  options_t  options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                    SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                    SPACE_DOMAIN_SIZE, 0 };
  int *ncolfactor, *ncolupdate, *parent, *vtx2front;
  int *first, *link, nfronts, J, K, u, vertex, vertex_root, count;
      /**************************************************
       declaration to uncomment if printing ordering
      ***************************************************
         FILE *fp1, *fp2;
         int  *perm,  *iperm;
      */
/*** decalage des indices couteux dans un premier temps:
****  A modifier dans une version ulterieure de MA41GD  */
  for (u = nvtx; u >= 0; u--)
   {
     xadj_pe[u] = xadj_pe[u] - 1;
   }
   for (K = nedges-1; K >= 0; K--)
   {
      adjncy[K] = adjncy[K] - 1;
   }
 /* initialization of the graph */
   mymalloc(G, 1, graph_t);
   G->xadj   = xadj_pe;
   G->adjncy = adjncy;
   mymalloc(G->vwght, nvtx, int);
   G->nvtx = nvtx;
   G->nedges = nedges;
   G->type = UNWEIGHTED;
   G->totvwght = nvtx;
   for (u = 0; u < nvtx; u++)
     G->vwght[u] = 1;
  /* main function of the Ordering */
   T = SPACE_ordering(G, options, cpus);
   nfronts = T->nfronts;
   ncolfactor = T->ncolfactor;
   ncolupdate = T->ncolupdate;
   parent = T->parent;
  /*    firstchild = T->firstchild; */
   vtx2front = T->vtx2front;
    /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
   mymalloc(first, nfronts, int);
   mymalloc(link, nvtx, int);
   for (J = 0; J < nfronts; J++)
      first[J] = -1;
   for (u = nvtx-1; u >= 0; u--)
      {
        J = vtx2front[u];
        link[u] = first[J];
        first[J] = u;
      }
  /* -----------------------------------------------------------
     fill the two arrays corresponding to the MUMPS tree structure
     ----------------------------------------------------------- */
  count = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
     {
       vertex_root = first[K];
       if (vertex_root == -1)
          {
            /* JY: I think this cannot happen */
            printf(" Internal error in mumps_pord (cf JY), %d\n",K);
            exit(-1);
          }
       /* for the principal column of the supervariable */
       if (parent[K] == -1)
          xadj_pe[vertex_root] = 0; /* root of the tree */
       else
          xadj_pe[vertex_root] = - (first[parent[K]]+1);
          nv[vertex_root] = ncolfactor[K] + ncolupdate[K];
          count++;
       for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
        /* for the secondary columns of the supervariable */
       {
         xadj_pe[vertex] = - (vertex_root+1);
         nv[vertex] = 0;
         count++;
        }
  }
  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
  free(G->vwght);
  free(G);
  freeElimTree(T);
  return (0);
}
/*********************************************************/
int mumps_pord_wnd
(
        int nvtx,
        int nedges,
        int *xadj_pe,
        int *adjncy,
        int *nv,
        int *totw
)
{
/**********************************
Argument Comments:
input:
-----
- nvtx   : dimension of the Problem (N)
- nedges : number of entries (NZ)
- adjncy : non-zeros entries (IW input)
- totw   : sum of the weigth of the vertices
input/output:
-------------
- xadj_pe : pointer through beginning of column non-zeros entries (PTRAR)
- on exit, "father array" (PE)
ouput:
------
- nv      : weight of the vertices
- on exit "nfront array" (NV)
*************************************/
        graph_t    *G;
        elimtree_t *T;
        timings_t  cpus[12];
        options_t  options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                    SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                    SPACE_DOMAIN_SIZE, 0 };
        int *ncolfactor, *ncolupdate, *parent, *vtx2front;
        int *first, *link, nfronts, J, K, u, vertex, vertex_root, count;
      /**************************************************
       declaration to uncomment if printing ordering
      ***************************************************
         FILE *fp1, *fp2;
         int  *perm,  *iperm;
      */
/*** decalage des indices couteux dans un premier temps:
****  A modifier dans une version ulterieure de MA41GD  */
        for (u = nvtx; u >= 0; u--)
        {
          xadj_pe[u] = xadj_pe[u] - 1;
        }
        for (K = nedges-1; K >= 0; K--)
        {
          adjncy[K] = adjncy[K] - 1;
        }
 /* initialization of the graph */
        mymalloc(G, 1, graph_t);
        G->xadj  = xadj_pe;
        G->adjncy= adjncy;
        mymalloc(G->vwght, nvtx, int);
        G->nvtx = nvtx;
        G->nedges = nedges;
        G->type = WEIGHTED;
        G->totvwght = (*totw);
        for (u = 0; u < nvtx; u++)
          G->vwght[u] = nv[u];
  /* main function of the Ordering */
        T = SPACE_ordering(G, options, cpus);
        nfronts = T->nfronts;
        ncolfactor = T->ncolfactor;
        ncolupdate = T->ncolupdate;
        parent = T->parent;
  /*    firstchild = T->firstchild; */
        vtx2front = T->vtx2front;
    /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
        mymalloc(first, nfronts, int);
        mymalloc(link, nvtx, int);
        for (J = 0; J < nfronts; J++)
          first[J] = -1;
        for (u = nvtx-1; u >= 0; u--)
        {
          J = vtx2front[u];
          link[u] = first[J];
          first[J] = u;
        }
  /* -----------------------------------------------------------
     fill the two arrays corresponding to the MUMPS tree structure
     ----------------------------------------------------------- */
  count = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
     {
        vertex_root = first[K];
        if (vertex_root == -1)
          {
            /* JY: I think this cannot happen */
            printf(" Internal error in mumps_pord (cf JY), %d\n",K);
            exit(-1);
          }
         /* for the principal column of the supervariable */
        if (parent[K] == -1)
          xadj_pe[vertex_root] = 0; /* root of the tree */
        else
          xadj_pe[vertex_root] = - (first[parent[K]]+1);
          nv[vertex_root] = ncolfactor[K] + ncolupdate[K];
          count++;
          for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
          /* for the secondary columns of the supervariable */
            {
              xadj_pe[vertex] = - (vertex_root+1);
              nv[vertex] = 0;
              count++;
        }
  }
  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
  free(G->vwght);
  free(G);
  freeElimTree(T);
  return (0);
}
#endif /* pord */
/************************************************************/
#if defined(scotch) || defined(ptscotch)
/*int esmumps( const int n, const int iwlen, int * const pe, const int pfree,
             int * const len, int * const iw, int * const nv, int * const elen,
             int * const last);*/
/* Fortran interface to SCOTCH */
/*#define MUMPS_SCOTCH    \
  F_SYMBOL(scotch,SCOTCH)*/
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
              int * const        ncmpa )
{
     *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                       lentab, iwtab, nvtab, elentab, lasttab );
}
#endif /* scotch */
#if defined(ptscotch)
/*#include "mpi.h"
#include <stdio.h>
#include "ptscotch.h"
int mumps_dgraphinit( SCOTCH_Dgraph *, MPI_Fint *, MPI_Fint *);
#define MUMPS_DGRAPHINIT        \
F_SYMBOL(dgraphinit,DGRAPHINIT)*/
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr)
{
  MPI_Comm  int_comm;
  int_comm = MPI_Comm_f2c(*comm);
  *ierr = SCOTCH_dgraphInit(graphptr, int_comm);
  return;
}
#endif
#if defined(parmetis)
void MUMPS_CALL
MUMPS_PARMETIS(int *first,      int *vertloctab, 
               int *edgeloctab, int *numflag, 
               int *options,    int *order, 
               int *sizes,      int *comm)
{
  MPI_Comm  int_comm;
  int_comm = MPI_Comm_f2c(*comm);
  ParMETIS_V3_NodeND(first, vertloctab, edgeloctab, numflag, options, order, sizes, &int_comm);
  return;
}
#endif
