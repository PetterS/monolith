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
#include "mumps_io_err.h"
#include "mumps_io_basic.h"
#if defined( MUMPS_WIN32 )
# include <string.h>
#endif
/* Exported global variables */
char* mumps_err;
MUMPS_INT* dim_mumps_err;
int mumps_err_max_len;
int err_flag;
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
pthread_mutex_t err_mutex;
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
/* Functions */
/* Keeps a C pointer to store error description string that will be
   displayed by the Fortran layers.
   * dim contains the size of the Fortran character array to store the
   description.
*/
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_ERR_STR(MUMPS_INT *dim, char* err_str, mumps_ftnlen l1){
  mumps_err = err_str;
  dim_mumps_err = (MUMPS_INT *) dim;
  mumps_err_max_len = (int) *dim;
  err_flag = 0;
  return;
}
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
MUMPS_INLINE int
mumps_io_protect_err()
{
  if(mumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_lock(&err_mutex);
  }
  return 0;
}
MUMPS_INLINE int
mumps_io_unprotect_err()
{
  if(mumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_unlock(&err_mutex);
  }
  return 0;
}
int
mumps_io_init_err_lock()
{
  pthread_mutex_init(&err_mutex,NULL);
  return 0;
}
int
mumps_io_destroy_err_lock()
{
  pthread_mutex_destroy(&err_mutex);
  return 0;
}
int
mumps_check_error_th()
{
  /* If err_flag != 0, then error_str is set */
  return err_flag;
}
#endif /* MUMPS_WIN32 && WITHOUT_PTHREAD */
int
mumps_io_error(int mumps_errno, const char* desc)
{
    int len;
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_protect_err();
#endif
  if(err_flag == 0){
    strncpy(mumps_err, desc, mumps_err_max_len);
    /* mumps_err is a FORTRAN string, we do not care about adding a final 0 */
    len = (int) strlen(desc);
    *dim_mumps_err = (len <= mumps_err_max_len ) ? len : mumps_err_max_len;
    err_flag = mumps_errno;
  }
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_unprotect_err();
#endif
  return mumps_errno;
}
int
mumps_io_sys_error(int mumps_errno, const char* desc)
{
  int len = 2; /* length of ": " */
  const char* _desc;
  char* _err;
#if defined( MUMPS_WIN32 )
  int _err_len;
#endif
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_protect_err();
#endif
  if(err_flag==0){
    if(desc == NULL) {
      _desc = "";
    } else {
        len += (int) strlen(desc);
      _desc = desc;
    }
#if ! defined( MUMPS_WIN32 )
    _err = strerror(errno);
    len += (int) strlen(_err);
    snprintf(mumps_err, mumps_err_max_len, "%s: %s", _desc, _err);
    /* mumps_err is a FORTRAN string, we do not care about adding a final 0 */
#else
    /* This a VERY UGLY workaround for snprintf: this function has been
     * integrated quite lately into the ANSI stdio: some windows compilers are
     * not up-to-date yet. */
    if( len >= mumps_err_max_len - 1 ) { /* then do not print sys error msg at all */
      len -= 2;
      len = (len >= mumps_err_max_len ) ? mumps_err_max_len - 1 : len;
      _err = strdup( _desc );
      _err[len] = '\0';
      sprintf(mumps_err, "%s", _err);
    } else {
      _err = strdup(strerror(errno));
      _err_len = (int) strlen(_err);
      /* We will use sprintf, so make space for the final '\0' ! */
      if((len + _err_len) >= mumps_err_max_len) {
        /* truncate _err, not to overtake mumps_err_max_len at the end. */
        _err[mumps_err_max_len - len - 1] = '\0';
        len = mumps_err_max_len - 1;
      } else {
        len += _err_len;
      }
      sprintf(mumps_err, "%s: %s", _desc, _err);
    }
    free(_err);
#endif
    *dim_mumps_err = (len <= mumps_err_max_len ) ? len : mumps_err_max_len;
    err_flag = mumps_errno;
  }
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_unprotect_err();
#endif
  return mumps_errno;
}
