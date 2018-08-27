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
#ifndef MUMPS_IO_THREAD_H
#define MUMPS_IO_THREAD_H
#include "mumps_compat.h"
#if ! defined (MUMPS_WIN32) && ! defined (WITHOUT_PTHREAD)
# include <unistd.h>
# include <pthread.h>
# include <sys/types.h>
# include <sys/time.h>
# include <time.h>
# define MAX_IO 20
# define MAX_FINISH_REQ 40
# define IO_FLAG_STOP 1
# define IO_FLAG_RUN 0
# define IO_READ 1
# define IO_WRITE 0
struct request_io{
  int inode;
  int req_num; /*request number*/
  void* addr;  /*memory address (either source or dest)*/
  long long size;    /* size of the requested io (unit=size of elementary mumps data)*/
  long long vaddr; /* virtual address for file management */
  int io_type; /*read or write*/
  int file_type; /* cb or lu or ... */
  pthread_cond_t local_cond;
  int int_local_cond;
};
/* Exported global variables */
extern int io_flag_stop,current_req_num;
extern pthread_t io_thread,main_thread;
extern pthread_mutex_t io_mutex;
extern pthread_cond_t cond_io,cond_nb_free_finished_requests,cond_nb_free_active_requests,cond_stop;
extern pthread_mutex_t io_mutex_cond;
extern int int_sem_io,int_sem_nb_free_finished_requests,int_sem_nb_free_active_requests,int_sem_stop;
extern int with_sem;
extern struct request_io *io_queue;
extern int first_active,last_active,nb_active;
extern int *finished_requests_inode,*finished_requests_id,first_finished_requests,
  last_finished_requests,nb_finished_requests,smallest_request_id;
extern int mumps_owns_mutex;
extern int test_request_called_from_mumps;
/* Exported functions */
void* mumps_async_thread_function_with_sem (void* arg);
int   mumps_is_there_finished_request_th(int* flag);
int   mumps_clean_request_th(int* request_id);
int   mumps_wait_req_sem_th(int *request_id);
int   mumps_test_request_th(int* request_id,int *flag);
int   mumps_wait_request_th(int *request_id);
int   mumps_low_level_init_ooc_c_th(int* async, int* ierr);
int   mumps_async_write_th(const int * strat_IO,void * address_block,long long block_size,
                           int * inode,int * request_arg,int * type,long long vaddr,int * ierr);
int   mumps_async_read_th(const int * strat_IO,void * address_block,long long block_size,int * inode,int * request_arg,
                           int * type,long long vaddr,int * ierr);
int mumps_clean_io_data_c_th(int *myid);
int mumps_get_sem(void *arg,int *value);
int mumps_wait_sem(void *arg,pthread_cond_t *cond);
int mumps_post_sem(void *arg,pthread_cond_t *cond);
int mumps_clean_finished_queue_th();
#endif /*_WIN32 && WITHOUT_PTHREAD*/
#endif /* MUMPS_IO_THREAD_H */
