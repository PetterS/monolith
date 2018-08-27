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
/* Written by JYL, march 2002 */
/* This file groups so far all C functions and symbols that vary with the
   arithmetic */
/* Header used for debug purpose only
#include <stdio.h>
*/
#include <string.h>
#include "mumps_common.h"
#if MUMPS_ARITH == MUMPS_ARITH_s
# include "smumps_c.h"
# define MUMPS_REAL    SMUMPS_REAL
# define MUMPS_COMPLEX SMUMPS_COMPLEX
#elif MUMPS_ARITH == MUMPS_ARITH_d
# include "dmumps_c.h"
# define MUMPS_REAL    DMUMPS_REAL
# define MUMPS_COMPLEX DMUMPS_COMPLEX
#elif MUMPS_ARITH == MUMPS_ARITH_c
# include "cmumps_c.h"
# define MUMPS_REAL    CMUMPS_REAL
# define MUMPS_COMPLEX CMUMPS_COMPLEX
#elif MUMPS_ARITH == MUMPS_ARITH_z
# include "zmumps_c.h"
# define MUMPS_REAL    ZMUMPS_REAL
# define MUMPS_COMPLEX ZMUMPS_COMPLEX
#endif
/**
 * F_SYM_ARITH is the same as F_SYMBOL (see mumps_commn.h) for the symbols
 * that depend on the arithmetic.
 * Example: For CMUMPS_XXX, first define
 *   #define CMUMPS_XXX F_SYM_ARITH(xxx,XXX) and then use
 *   CMUMPS_XXX in the code to get rid of any symbol convention annoyance.
 */
#if MUMPS_ARITH == MUMPS_ARITH_s
# if defined(UPPER) || defined(MUMPS_WIN32)
#  define F_SYM_ARITH(lower_case,upper_case) SMUMPS_##upper_case
# elif defined(Add_)
#  define F_SYM_ARITH(lower_case,upper_case) smumps_##lower_case##_
# elif defined(Add__)
#  define F_SYM_ARITH(lower_case,upper_case) smumps_##lower_case##__
# else
#  define F_SYM_ARITH(lower_case,upper_case) smumps_##lower_case
# endif
#elif MUMPS_ARITH == MUMPS_ARITH_d
# if defined(UPPER) || defined(MUMPS_WIN32)
#  define F_SYM_ARITH(lower_case,upper_case) DMUMPS_##upper_case
# elif defined(Add_)
#  define F_SYM_ARITH(lower_case,upper_case) dmumps_##lower_case##_
# elif defined(Add__)
#  define F_SYM_ARITH(lower_case,upper_case) dmumps_##lower_case##__
# else
#  define F_SYM_ARITH(lower_case,upper_case) dmumps_##lower_case
# endif
#elif MUMPS_ARITH == MUMPS_ARITH_c
# if defined(UPPER) || defined(MUMPS_WIN32)
#  define F_SYM_ARITH(lower_case,upper_case) CMUMPS_##upper_case
# elif defined(Add_)
#  define F_SYM_ARITH(lower_case,upper_case) cmumps_##lower_case##_
# elif defined(Add__)
#  define F_SYM_ARITH(lower_case,upper_case) cmumps_##lower_case##__
# else
#  define F_SYM_ARITH(lower_case,upper_case) cmumps_##lower_case
# endif
#elif MUMPS_ARITH == MUMPS_ARITH_z
# if defined(UPPER) || defined(MUMPS_WIN32)
#  define F_SYM_ARITH(lower_case,upper_case) ZMUMPS_##upper_case
# elif defined(Add_)
#  define F_SYM_ARITH(lower_case,upper_case) zmumps_##lower_case##_
# elif defined(Add__)
#  define F_SYM_ARITH(lower_case,upper_case) zmumps_##lower_case##__
# else
#  define F_SYM_ARITH(lower_case,upper_case) zmumps_##lower_case
# endif
#endif
#define MUMPS_F77       \
    F_SYM_ARITH(f77,F77)
void MUMPS_CALL
MUMPS_F77( MUMPS_INT      *job,
           MUMPS_INT      *sym,
           MUMPS_INT      *par,
           MUMPS_INT      *comm_fortran,
           MUMPS_INT      *n,
           MUMPS_INT      *icntl,
           MUMPS_REAL     *cntl,
           MUMPS_INT      *nz,
           MUMPS_INT      *irn,
           MUMPS_INT      *irn_avail,
           MUMPS_INT      *jcn,
           MUMPS_INT      *jcn_avail,
           MUMPS_COMPLEX  *a,
           MUMPS_INT      *a_avail,
           MUMPS_INT      *nz_loc,
           MUMPS_INT      *irn_loc,
           MUMPS_INT      *irn_loc_avail,
           MUMPS_INT      *jcn_loc,
           MUMPS_INT      *jcn_loc_avail,
           MUMPS_COMPLEX  *a_loc,
           MUMPS_INT      *a_loc_avail,
           MUMPS_INT      *nelt,
           MUMPS_INT      *eltptr,
           MUMPS_INT      *eltptr_avail,
           MUMPS_INT      *eltvar,
           MUMPS_INT      *eltvar_avail,
           MUMPS_COMPLEX  *a_elt,
           MUMPS_INT      *a_elt_avail,
           MUMPS_INT      *perm_in,
           MUMPS_INT      *perm_in_avail,
           MUMPS_COMPLEX  *rhs,
           MUMPS_INT      *rhs_avail,
           MUMPS_COMPLEX  *redrhs,
           MUMPS_INT      *redrhs_avail,
           MUMPS_INT      *info,
           MUMPS_REAL     *rinfo,
           MUMPS_INT      *infog,
           MUMPS_REAL     *rinfog,
           MUMPS_INT      *deficiency,
           MUMPS_INT      *lwk_user,
           MUMPS_INT      *size_schur,
           MUMPS_INT      *listvar_schur,
           MUMPS_INT      *listvar_schur_avail,
           MUMPS_COMPLEX  *schur,
           MUMPS_INT      *schur_avail,
           MUMPS_COMPLEX  *wk_user,
           MUMPS_INT      *wk_user_avail,
           MUMPS_REAL     *colsca,
           MUMPS_INT      *colsca_avail,
           MUMPS_REAL     *rowsca,
           MUMPS_INT      *rowsca_avail,
           MUMPS_INT      *instance_number,
           MUMPS_INT      *nrhs,
           MUMPS_INT      *lrhs,
           MUMPS_INT      *lredrhs,
           MUMPS_COMPLEX  *rhs_sparse,
           MUMPS_INT      *rhs_sparse_avail,
           MUMPS_COMPLEX  *sol_loc,
           MUMPS_INT      *sol_loc_avail,
           MUMPS_INT      *irhs_sparse,
           MUMPS_INT      *irhs_sparse_avail,
           MUMPS_INT      *irhs_ptr,
           MUMPS_INT      *irhs_ptr_avail,
           MUMPS_INT      *isol_loc,
           MUMPS_INT      *isol_loc_avail,
           MUMPS_INT      *nz_rhs,
           MUMPS_INT      *lsol_loc,
           MUMPS_INT      *schur_mloc,
           MUMPS_INT      *schur_nloc,
           MUMPS_INT      *schur_lld,
           MUMPS_INT      *schur_mblock,
           MUMPS_INT      *schur_nblock,
           MUMPS_INT      *schur_nprow,
           MUMPS_INT      *schur_npcol,
           MUMPS_INT      *ooc_tmpdir,
           MUMPS_INT      *ooc_prefix,
           MUMPS_INT      *write_problem,
           MUMPS_INT      *ooc_tmpdirlen,
           MUMPS_INT      *ooc_prefixlen,
           MUMPS_INT      *write_problemlen
           );
#ifdef return_scaling
/*
 * Those two are static. They are passed inside cmumps_f77 but
 * might also be changed on return by MUMPS_AFFECT_COLSCA/ROWSCA
 * NB: They are put here because they use MUMPS_REAL and need thus
 * one symbol per arithmetic.
 */
#if MUMPS_ARITH == MUMPS_ARITH_s
# define MUMPS_COLSCA_STATIC SMUMPS_COLSCA_STATIC
# define MUMPS_ROWSCA_STATIC SMUMPS_ROWSCA_STATIC
#elif MUMPS_ARITH == MUMPS_ARITH_d
# define MUMPS_COLSCA_STATIC SMUMPS_COLSCA_STATIC
# define MUMPS_ROWSCA_STATIC SMUMPS_ROWSCA_STATIC
#elif MUMPS_ARITH == MUMPS_ARITH_c
# define MUMPS_COLSCA_STATIC CMUMPS_COLSCA_STATIC
# define MUMPS_ROWSCA_STATIC CMUMPS_ROWSCA_STATIC
#elif MUMPS_ARITH == MUMPS_ARITH_z
# define MUMPS_COLSCA_STATIC ZMUMPS_COLSCA_STATIC
# define MUMPS_ROWSCA_STATIC ZMUMPS_ROWSCA_STATIC
#endif
static MUMPS_REAL * MUMPS_COLSCA_STATIC;
static MUMPS_REAL * MUMPS_ROWSCA_STATIC;
#define MUMPS_AFFECT_COLSCA \
    F_SYM_ARITH(affect_colsca,AFFECT_COLSCA)
void MUMPS_CALL
MUMPS_AFFECT_COLSCA(MUMPS_REAL * f77colsca)
{
  MUMPS_COLSCA_STATIC = f77colsca;
}
#define MUMPS_NULLIFY_C_COLSCA \
    F_SYM_ARITH(nullify_c_colsca,NULLIFY_C_COLSCA)
void MUMPS_CALL
MUMPS_NULLIFY_C_COLSCA()
{
  MUMPS_COLSCA_STATIC = 0;
}
#define MUMPS_AFFECT_ROWSCA \
    F_SYM_ARITH(affect_rowsca,AFFECT_ROWSCA)
void MUMPS_CALL
MUMPS_AFFECT_ROWSCA(MUMPS_REAL * f77rowsca)
{
  MUMPS_ROWSCA_STATIC = f77rowsca;
}
#define MUMPS_NULLIFY_C_ROWSCA \
    F_SYM_ARITH(nullify_c_rowsca,NULLIFY_C_ROWSCA)
void MUMPS_CALL
MUMPS_NULLIFY_C_ROWSCA()
{
  MUMPS_ROWSCA_STATIC = 0;
}
#endif /* return_scaling */
#if MUMPS_ARITH == MUMPS_ARITH_s
# define mumps_c       smumps_c
# define MUMPS_STRUC_C SMUMPS_STRUC_C
#elif MUMPS_ARITH == MUMPS_ARITH_d
# define mumps_c       dmumps_c
# define MUMPS_STRUC_C DMUMPS_STRUC_C
#elif MUMPS_ARITH == MUMPS_ARITH_c
# define mumps_c       cmumps_c
# define MUMPS_STRUC_C CMUMPS_STRUC_C
#elif MUMPS_ARITH == MUMPS_ARITH_z
# define mumps_c       zmumps_c
# define MUMPS_STRUC_C ZMUMPS_STRUC_C
#endif
void MUMPS_CALL
mumps_c(MUMPS_STRUC_C * mumps_par)
{
    /*
     * The following local variables will 
     *  be passed to the F77 interface.
     */
    MUMPS_INT *icntl;
    MUMPS_REAL *cntl;
    MUMPS_INT *irn; MUMPS_INT *jcn; MUMPS_COMPLEX *a;
    MUMPS_INT *irn_loc; MUMPS_INT *jcn_loc; MUMPS_COMPLEX *a_loc;
    MUMPS_INT *eltptr, *eltvar; MUMPS_COMPLEX *a_elt;
    MUMPS_INT *perm_in; MUMPS_INT perm_in_avail;
    MUMPS_INT *listvar_schur; MUMPS_INT listvar_schur_avail;
    MUMPS_COMPLEX *schur; MUMPS_INT schur_avail;
    MUMPS_COMPLEX *rhs; MUMPS_COMPLEX *redrhs;
    MUMPS_COMPLEX *wk_user; MUMPS_INT wk_user_avail;
    MUMPS_REAL *colsca; MUMPS_REAL *rowsca;
    MUMPS_COMPLEX *rhs_sparse, *sol_loc;
    MUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    MUMPS_INT irn_avail, jcn_avail, a_avail, rhs_avail, redrhs_avail;
    /* These are actually used
     * as booleans, but we stick
     * to simple types for the
     * C-F77 interface */
    MUMPS_INT irn_loc_avail, jcn_loc_avail, a_loc_avail;
    MUMPS_INT eltptr_avail, eltvar_avail, a_elt_avail;
    MUMPS_INT colsca_avail, rowsca_avail;
    MUMPS_INT irhs_ptr_avail, rhs_sparse_avail, sol_loc_avail;
    MUMPS_INT irhs_sparse_avail, isol_loc_avail;
    MUMPS_INT *info; MUMPS_INT *infog;
    MUMPS_REAL *rinfo; MUMPS_REAL *rinfog;
    MUMPS_INT ooc_tmpdir[255]; MUMPS_INT ooc_prefix[63];
    MUMPS_INT write_problem[255];
    /* Other local variables */
    MUMPS_INT idummy; MUMPS_INT *idummyp;
    MUMPS_REAL rdummy; MUMPS_REAL *rdummyp;
    MUMPS_COMPLEX cdummy; MUMPS_COMPLEX *cdummyp;
    /* String lengths to be passed to Fortran by address */
    int ooc_tmpdirlen;
    int ooc_prefixlen;
    int write_problemlen;
    int i;
    static const MUMPS_INT no = 0;
    static const MUMPS_INT yes = 1;
    idummyp = &idummy;
    cdummyp = &cdummy;
    rdummyp = &rdummy;
#ifdef return_scaling
    /* Don't forget to initialize those two before
     * each call to mumps as we may copy values from
     * old instances otherwise ! */
    MUMPS_COLSCA_STATIC=0;
    MUMPS_ROWSCA_STATIC=0;
#endif
    /* Initialize pointers to zero for job == -1 */
    if ( mumps_par->job == -1 )
      { /* job = -1: we just reset all pointers to 0 */
        mumps_par->irn=0; mumps_par->jcn=0; mumps_par->a=0; mumps_par->rhs=0; mumps_par->wk_user=0;
        mumps_par->redrhs=0;
        mumps_par->eltptr=0; mumps_par->eltvar=0; mumps_par->a_elt=0; mumps_par->perm_in=0; mumps_par->sym_perm=0; mumps_par->uns_perm=0; mumps_par->irn_loc=0;mumps_par->jcn_loc=0;mumps_par->a_loc=0; mumps_par->listvar_schur=0;mumps_par->schur=0;mumps_par->mapping=0;mumps_par->pivnul_list=0;mumps_par->colsca=0;mumps_par->rowsca=0; mumps_par->rhs_sparse=0; mumps_par->irhs_sparse=0; mumps_par->sol_loc=0; mumps_par->irhs_ptr=0; mumps_par->isol_loc=0;
        strcpy(mumps_par->ooc_tmpdir,"NAME_NOT_INITIALIZED");
        strcpy(mumps_par->ooc_prefix,"NAME_NOT_INITIALIZED");
        strcpy(mumps_par->write_problem,"NAME_NOT_INITIALIZED");
        strncpy(mumps_par->version_number,MUMPS_VERSION,MUMPS_VERSION_MAX_LEN);
        mumps_par->version_number[MUMPS_VERSION_MAX_LEN+1] = '\0';
        /* Next line initializes scalars to arbitrary values.
         * Some of those will anyway be overwritten during the
         * call to Fortran routine [SDCZ]MUMPS_INIT_PHASE */
        mumps_par->n=0; mumps_par->nz=0; mumps_par->nz_loc=0; mumps_par->nelt=0;mumps_par->instance_number=0;mumps_par->deficiency=0;mumps_par->lwk_user=0;mumps_par->size_schur=0;mumps_par->lrhs=0; mumps_par->lredrhs=0; mumps_par->nrhs=0; mumps_par->nz_rhs=0; mumps_par->lsol_loc=0;
 mumps_par->schur_mloc=0; mumps_par->schur_nloc=0; mumps_par->schur_lld=0; mumps_par->mblock=0; mumps_par->nblock=0; mumps_par->nprow=0; mumps_par->npcol=0;
      }
     ooc_tmpdirlen=(int)strlen(mumps_par->ooc_tmpdir);
     ooc_prefixlen=(int)strlen(mumps_par->ooc_prefix);
     write_problemlen=(int)strlen(mumps_par->write_problem);
    /* Avoid the use of strnlen which may not be
     * available on all systems. Allow strings without
     * \0 at the end, if the file is not found, the
     * Fortran layer is responsible for raising an
     * error.  */
    if(ooc_tmpdirlen > 255){
        ooc_tmpdirlen=255;
      }
    if(ooc_prefixlen > 63){
        ooc_prefixlen=63;
      }
    if(write_problemlen > 255){
        write_problemlen=255;
      }
    /*
     * Extract info from the C structure to call the F77 interface. The
     * following macro avoids repeating the same code with risks of errors.
     */
#define EXTRACT_POINTERS(component,dummypointer) \
    if ( mumps_par-> component == 0) \
      { component = dummypointer; \
        component ## _avail = no; }  \
    else  \
      { component = mumps_par-> component; \
        component ## _avail = yes; }
    /*
     * For example, EXTRACT_POINTERS(irn,idummyp) produces the following line of code:
       if (mumps_par->irn== 0) {irn= idummyp;irn_avail = no; } else {  irn  = mumps_par->irn;irn_avail = yes; } ;
     * which says that irn is set to mumps_par->irn except if
     * mumps_par->irn is 0, which means that it is not available.
     */
    EXTRACT_POINTERS(irn,idummyp);
    EXTRACT_POINTERS(jcn,idummyp);
    EXTRACT_POINTERS(rhs,cdummyp);
    EXTRACT_POINTERS(wk_user,cdummyp);
    EXTRACT_POINTERS(redrhs,cdummyp);
    EXTRACT_POINTERS(irn_loc,idummyp);
    EXTRACT_POINTERS(jcn_loc,idummyp);
    EXTRACT_POINTERS(a_loc,cdummyp);
    EXTRACT_POINTERS(a,cdummyp);
    EXTRACT_POINTERS(eltptr,idummyp);
    EXTRACT_POINTERS(eltvar,idummyp);
    EXTRACT_POINTERS(a_elt,cdummyp);
    EXTRACT_POINTERS(perm_in,idummyp);
    EXTRACT_POINTERS(listvar_schur,idummyp);
    EXTRACT_POINTERS(schur,cdummyp);
    EXTRACT_POINTERS(colsca,rdummyp);
    EXTRACT_POINTERS(rowsca,rdummyp);
    EXTRACT_POINTERS(rhs_sparse,cdummyp);
    EXTRACT_POINTERS(sol_loc,cdummyp);
    EXTRACT_POINTERS(irhs_sparse,idummyp);
    EXTRACT_POINTERS(isol_loc,idummyp);
    EXTRACT_POINTERS(irhs_ptr,idummyp);
    /* printf("irn_avail,jcn_avail, rhs_avail, a_avail, eltptr_avail, eltvar_avail,a_elt_avail,perm_in_avail= %d %d %d %d %d %d %d \n", irn_avail,jcn_avail, rhs_avail, a_avail, eltptr_avail, eltvar_avail, a_elt_avail, perm_in_avail);*/
    /*
     * Extract integers (input) or pointers that are
     * always allocated (such as ICNTL, INFO, ...)
     */
    /* size_schur = mumps_par->size_schur; */
    /* instance_number = mumps_par->instance_number; */
    icntl = mumps_par->icntl;
    cntl = mumps_par->cntl;
    info = mumps_par->info;
    infog = mumps_par->infog;
    rinfo = mumps_par->rinfo;
    rinfog = mumps_par->rinfog;
    for(i=0;i<ooc_tmpdirlen;i++){
      ooc_tmpdir[i]=(int)mumps_par->ooc_tmpdir[i];
    }
    for(i=0;i<ooc_prefixlen;i++){
      ooc_prefix[i]=(int)mumps_par->ooc_prefix[i];
    }
    for(i=0;i<write_problemlen;i++){
      write_problem[i]=(int)mumps_par->write_problem[i];
    }
    /* Call F77 interface */
    MUMPS_F77(&(mumps_par->job), &(mumps_par->sym), &(mumps_par->par), &(mumps_par->comm_fortran),
          &(mumps_par->n), icntl, cntl,
          &(mumps_par->nz), irn, &irn_avail, jcn, &jcn_avail, a, &a_avail,
          &(mumps_par->nz_loc), irn_loc, &irn_loc_avail, jcn_loc, &jcn_loc_avail,
          a_loc, &a_loc_avail,
          &(mumps_par->nelt), eltptr, &eltptr_avail, eltvar, &eltvar_avail, a_elt, &a_elt_avail,
          perm_in, &perm_in_avail,
          rhs, &rhs_avail, redrhs, &redrhs_avail, info, rinfo, infog, rinfog,
          &(mumps_par->deficiency), &(mumps_par->lwk_user), &(mumps_par->size_schur), listvar_schur, &listvar_schur_avail, schur,
          &schur_avail, wk_user, &wk_user_avail, colsca, &colsca_avail, rowsca, &rowsca_avail,
          &(mumps_par->instance_number), &(mumps_par->nrhs), &(mumps_par->lrhs),
          &(mumps_par->lredrhs),
          rhs_sparse, &rhs_sparse_avail, sol_loc, &sol_loc_avail, irhs_sparse,
          &irhs_sparse_avail, irhs_ptr, &irhs_ptr_avail, isol_loc,
          &isol_loc_avail, &(mumps_par->nz_rhs), &(mumps_par->lsol_loc)
          , &(mumps_par->schur_mloc)
          , &(mumps_par->schur_nloc)
          , &(mumps_par->schur_lld)
          , &(mumps_par->mblock)
          , &(mumps_par->nblock)
          , &(mumps_par->nprow)
          , &(mumps_par->npcol)
          , ooc_tmpdir
          , ooc_prefix
          , write_problem
          , &ooc_tmpdirlen
          , &ooc_prefixlen
          , &write_problemlen
    );
    /*
     * mapping and pivnul_list are usually 0 except if
     * MUMPS_AFFECT_MAPPING/MUMPS_AFFECT_PIVNUL_LIST was called.
     */
    mumps_par->mapping=mumps_get_mapping();
    mumps_par->pivnul_list=mumps_get_pivnul_list();
    /* to get permutations computed during analysis */
    mumps_par->sym_perm=mumps_get_sym_perm();
    mumps_par->uns_perm=mumps_get_uns_perm();
#ifdef return_scaling
    /*
     * colsca/rowsca can either be user data or have been
     * modified within mumps by calls to MUMPS_AFFECT_COLSCA/ROWSCA.
     */
    if (colsca_avail == no) mumps_par->colsca = MUMPS_COLSCA_STATIC;
    if (rowsca_avail == no) mumps_par->rowsca = MUMPS_ROWSCA_STATIC;
#endif
}
