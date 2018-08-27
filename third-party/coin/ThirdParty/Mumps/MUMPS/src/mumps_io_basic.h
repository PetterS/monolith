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
#ifndef MUMPS_IO_BASIC_H
#define MUMPS_IO_BASIC_H
#include "mumps_compat.h"
#if ! defined(WITHOUT_PTHREAD) && defined(MUMPS_WIN32)
# define WITHOUT_PTHREAD 1
#endif
#if defined(_AIX)
# if ! defined(_ALL_SOURCE)
/* Macro needed for direct I/O on IBM AIX */
#  define _ALL_SOURCE 1
# endif
#endif
#if ! defined (MUMPS_WIN32)
# if ! defined(_XOPEN_SOURCE)
/* Setting this macro avoids the warnings ("missing
 * prototype") related to the use of pread /pwrite */
#  define _XOPEN_SOURCE 500
# endif
#endif
#define MAX_FILE_SIZE 1879048192 /* (2^31)-1-(2^27) */
/* #define MAX_FILE_SIZE 1000000 */  /* (2^31)-1-(2^27) */
/*                                                      */
/* Important Note :                                     */
/* ================                                     */
/* On GNU systems, __USE_GNU must be defined to have    */
/* access to the O_DIRECT I/O flag.                     */
/*                                                      */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if ! defined (MUMPS_WIN32)
# include <unistd.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <fcntl.h>
# include <sys/time.h>
# include <time.h>
#endif
#if ! defined (MUMPS_WIN32)
# define MUMPS_IO_FLAG_O_DIRECT 0
#endif
/* Force WITH_PFUNC on architectures where we know that it should work */
#if (defined (sgi) || defined (__sgi)) || defined(_AIX) || (defined(sun) || defined(__sun)) || defined(_GNU_SOURCE)
# undef WITH_PFUNC
# define WITH_PFUNC
#endif
#define IO_SYNC      0
#define IO_ASYNC_TH  1
#define IO_ASYNC_AIO 2
#define IO_READ 1
#define IO_WRITE 0
#define UNITIALIZED "NAME_NOT_INITIALIZED"
#define MUMPS_OOC_DEFAULT_DIR "/tmp"
#ifdef MUMPS_WIN32
# define SEPARATOR "\\"
#else
# define SEPARATOR "/"
#endif
/* #define NB_FILE_TYPE_FACTO 1 */
/* #define NB_FILE_TYPE_SOLVE 1 */
#define my_max(x,y) ( (x) > (y) ? (x) : (y) ) 
#define my_ceil(x) ( (int)(x) >= (x) ? (int)(x) : ( (int)(x) + 1 ) )
typedef struct __mumps_file_struct{
  int write_pos;
  int current_pos;
  int is_opened;
#if ! defined (MUMPS_WIN32)
  int file;
#else
  FILE* file;
#endif
  char name[351]; /* Should be large enough to hold tmpdir, prefix, suffix */
}mumps_file_struct;
typedef struct __mumps_file_type{
#if ! defined (MUMPS_WIN32)
  int mumps_flag_open;
#else
  char mumps_flag_open[6];
#endif
  int mumps_io_current_file_number;
  int mumps_io_last_file_opened;
  int mumps_io_nb_file_opened;
  int mumps_io_nb_file;
  mumps_file_struct* mumps_io_pfile_pointer_array;
  mumps_file_struct* mumps_io_current_file;
}mumps_file_type;
/* Exported global variables */
#if ! defined (MUMPS_WIN32)
# if defined (WITH_PFUNC) && ! defined (WITHOUT_PTHREAD)
#  include <pthread.h>
extern pthread_mutex_t mumps_io_pwrite_mutex;
# endif
/* extern int* mumps_io_pfile_pointer_array; */
/* extern int* mumps_io_current_file; */
/* #else /\*_WIN32*\/ */
/* extern FILE** mumps_io_current_file; */
/* extern FILE** mumps_io_pfile_pointer_array; */
#endif /* MUMPS_WIN32 */
/*extern mumps_file_struct* mumps_io_pfile_pointer_array;
  extern mumps_file_struct* mumps_io_current_file;*/
extern mumps_file_type* mumps_files;
/* extern int mumps_io_current_file_number; */
extern char* mumps_ooc_file_prefix;
/* extern char** mumps_io_pfile_name; */
/* extern int mumps_io_current_file_position; */
/* extern int mumps_io_write_pos; */
/* extern int mumps_io_last_file_opened; */
extern int mumps_elementary_data_size;
extern int mumps_io_is_init_called;
extern int mumps_io_myid;
extern int mumps_io_max_file_size;
/* extern int mumps_io_nb_file; */
extern int mumps_io_flag_async;
extern int mumps_io_k211;
/* extern int mumps_flag_open; */
extern int directio_flag;
extern int mumps_directio_flag;
extern int mumps_io_nb_file_type;
/* Exported functions */
int mumps_set_file(int type,int file_number_arg);
void mumps_update_current_file_position(mumps_file_struct* file_arg);
int mumps_compute_where_to_write(const double to_be_written,const int type,long long vaddr,size_t already_written);
int mumps_prepare_pointers_for_write(double to_be_written,int * pos_in_file, int * file_number,const int type,long long vaddr,size_t already_written);
int mumps_io_do_write_block(void * address_block,long long block_size,int * type,long long vaddr,int * ierr);
int mumps_io_do_read_block(void * address_block,long long block_size,int * type,long long vaddr,int * ierr);
int mumps_compute_nb_concerned_files(long long block_size,int * nb_concerned_files,long long vaddr);
MUMPS_INLINE int mumps_gen_file_info(long long vaddr, int * pos, int * file);
int mumps_free_file_pointers(int* step);
int mumps_init_file_structure(int *_myid, long long *total_size_io,int *size_element,int *nb_file_type,int *flag_tab);
int mumps_init_file_name(char* mumps_dir,char* mumps_file,int* mumps_dim_dir,int* mumps_dim_file,int* _myid);
void mumps_io_init_file_struct(int* nb,int which);
int mumps_io_alloc_file_struct(int* nb,int which);
int mumps_io_get_nb_files(int* nb_files, const int* type);
int mumps_io_get_file_name(int* indice,char* name,int* length,int* type);
int mumps_io_alloc_pointers(int * nb_file_type, int * dim);
int mumps_io_init_vars(int* myid_arg,int* size_element,int* async_arg);
int mumps_io_set_file_name(int* indice,char* name,int* length,int* type);
int mumps_io_open_files_for_read();
int mumps_io_set_last_file(int* dim,int* type);
int mumps_io_write__(void *file, void *loc_add, size_t write_size, int where,int type);
#if ! defined (MUMPS_WIN32)
int mumps_io_write_os_buff__(void *file, void *loc_add, size_t write_size, int where);
int mumps_io_write_direct_io__(void *file, void *loc_addr, size_t write_size, int where,int type);
int mumps_io_flush_write__(int type);
#else
int mumps_io_write_win32__(void *file, void *loc_add, size_t write_size, int where);
#endif
int mumps_io_read__(void * file,void * loc_addr,size_t size,int local_offset,int type);
#if ! defined (MUMPS_WIN32)
int mumps_io_read_os_buff__(void * file,void * loc_addr,size_t size,int local_offset);
int mumps_io_read_direct_io__(void * file,void * loc_addr,size_t size,int local_offset,int type);
#else
int mumps_io_read_win32__(void * file,void * loc_addr,size_t size,int local_offset);
#endif
int mumps_compute_file_size(void *file,size_t *size);
#if ! defined (MUMPS_WIN32) && ! defined (WITHOUT_PTHREAD)
# ifdef WITH_PFUNC
int mumps_io_protect_pointers();
int mumps_io_unprotect_pointers();
int mumps_io_init_pointers_lock();
int mumps_io_destroy_pointers_lock();
# endif /* WITH_PFUNC */
#endif /* MUMPS_WIN32 */
#endif /* MUMPS_IO_BASIC_H */
