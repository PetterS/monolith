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
#include <errno.h>
#include "mumps_common.h"
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
# include <pthread.h>
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
extern pthread_mutex_t err_mutex;
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
/* Exported functions */
#define MUMPS_LOW_LEVEL_INIT_ERR_STR \
    F_SYMBOL(low_level_init_err_str,LOW_LEVEL_INIT_ERR_STR)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_ERR_STR( MUMPS_INT *dim, char *err_str, mumps_ftnlen l1 );
/* Export an error to the Fortran layer
   Returns mumps_errno for convenience */
int mumps_io_error(int mumps_errno, const char* desc);
/* Export a system error to the Fortran layer (errno must be set)
   Returns mumps_errno for convenience */
int mumps_io_sys_error(int mumps_errno, const char* desc);
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
int mumps_io_init_err_lock();
int mumps_io_destroy_err_lock();
int mumps_check_error_th();
MUMPS_INLINE int mumps_io_protect_err();
MUMPS_INLINE int mumps_io_unprotect_err();
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
