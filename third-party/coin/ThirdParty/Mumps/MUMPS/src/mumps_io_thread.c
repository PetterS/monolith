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
#include "mumps_io_basic.h"
#include "mumps_io_err.h"
#include "mumps_io_thread.h"
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
/* Exported global variables */
int io_flag_stop,current_req_num;
pthread_t io_thread,main_thread;
pthread_mutex_t io_mutex;
pthread_cond_t cond_io,cond_nb_free_finished_requests,cond_nb_free_active_requests,cond_stop;
pthread_mutex_t io_mutex_cond;
int int_sem_io,int_sem_nb_free_finished_requests,int_sem_nb_free_active_requests,int_sem_stop;
int with_sem;
struct request_io *io_queue;
int first_active,last_active,nb_active;
int *finished_requests_inode,*finished_requests_id,first_finished_requests,
  last_finished_requests,nb_finished_requests,smallest_request_id;
int mumps_owns_mutex;
int test_request_called_from_mumps;
/* Other global variables */
double inactive_time_io_thread;
int time_flag_io_thread;
struct timeval origin_time_io_thread;
/**
 * Main function of the io thread when semaphores are used.
 */
void*  mumps_async_thread_function_with_sem (void* arg){
   struct request_io *current_io_request;
   int ierr,_sem_stop;
   struct timeval start_time,end_time;
   int ret_code;
   for (;;){      
     gettimeofday(&start_time,NULL);
       if(with_sem==2){
         mumps_wait_sem(&int_sem_io,&cond_io);
       }
     /*     sem_wait(&sem_io);  */
     gettimeofday(&end_time,NULL);
     if(time_flag_io_thread){
       inactive_time_io_thread=inactive_time_io_thread+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
     }else{
       inactive_time_io_thread=((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)origin_time_io_thread.tv_sec+((double)origin_time_io_thread.tv_usec/1000000));
     }
     if(!time_flag_io_thread){
       time_flag_io_thread=1;
     }
     /* Check if the main thread ordered to stop this slave thread */
     /*     sem_getvalue(&sem_stop,&_sem_stop); */
       if(with_sem==2){
         mumps_get_sem(&int_sem_stop,&_sem_stop);
       }
     if(_sem_stop==IO_FLAG_STOP){
       /* The thread must stop */
       break; /* Breaks the while loop. */
     }
      current_io_request=&io_queue[first_active];
      switch(current_io_request->io_type)
        {
         case IO_WRITE:
           ret_code=mumps_io_do_write_block(current_io_request->addr,
                                     current_io_request->size,
                                     &(current_io_request->file_type),
                                     current_io_request->vaddr,
                                     &ierr);
           if(ret_code<0){
             goto end;
           }
           break;
         case IO_READ:
           ret_code=mumps_io_do_read_block(current_io_request->addr,
           current_io_request->size,
           &(current_io_request->file_type),
           current_io_request->vaddr,
           &ierr);
           if(ret_code<0){
             goto end;
           }
           break;
         default:
           printf("Error : Mumps_IO : Operation %d is neither READ nor WRITE\n",current_io_request->io_type);
           exit (-3);
         }
      /* Notify that the IO was performed */
      /* Wait that finished_requests queue could register 
       the notification */
       if(with_sem==2){
         mumps_wait_sem(&int_sem_nb_free_finished_requests,&cond_nb_free_finished_requests);
       }
      pthread_mutex_lock(&io_mutex);
      /* Updates active queue bounds */
      /* Register the notification in finished_requests queue
       and updates its bounds. */
      finished_requests_id[last_finished_requests]=current_io_request->req_num;
      finished_requests_inode[last_finished_requests]=current_io_request->inode;
      last_finished_requests=(last_finished_requests+1)%(MAX_FINISH_REQ); /* ??? */
      nb_finished_requests++;      
      /* Realeases the lock : ***UNLOCK*** */
      nb_active--;      
      if(first_active<MAX_IO-1){
        first_active++; 
      }
      else{
        first_active=0;
      }
      if(with_sem==2){
        mumps_post_sem(&(current_io_request->int_local_cond),&(current_io_request->local_cond));
      }
      pthread_mutex_unlock(&io_mutex);  
      /* Finally increases the number of free active requests.*/
      /*      sem_post(&sem_nb_free_active_requests); */
          mumps_post_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
   }
 end:
   /* The main thread ordered the end of the IO thread (it changed sem_stop).
    We exit. */
   pthread_exit(NULL);
/* Not reached */
   return NULL;
}
int mumps_test_request_th(int* request_id,int *flag){
  /* Tests if the request "request_id" has finished. It sets the flag  */
  /* argument to 1 if the request has finished (0 otherwise)           */
  int request_pos;
  int i;
  i=mumps_check_error_th();
  if(i!=0){
    return i;
  }
  pthread_mutex_lock(&io_mutex);
  /*  printf("entering test !!! \n"); */
  if(*request_id < smallest_request_id){
    *flag=1;
    /*    exit (-2); */
  }else{
    if(nb_finished_requests==0){
      *flag=0;
    }else{
      request_pos=(first_finished_requests+nb_finished_requests-1)%(MAX_IO*2);    
      if(*request_id > finished_requests_id[request_pos]){
        /*the request has not been treated yet since it is not in 
        the list of treated requests*/
        i=0;
        /*this loop is only for checking (no special treatment is done*/
        while(i<nb_active){
          request_pos=(first_active+i)%(MAX_IO);
          if(io_queue[request_pos].req_num==*request_id){
            break;
          }
          i++;
        }
        if(i==nb_active){
          return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_test_request_th (1))\n");
        }
        *flag=0;
      }else{
        /*the request has been treated yet since it is the list of
          treated requests*/
        i=0;
        while(i<nb_finished_requests){
          request_pos=(first_finished_requests+i)%(MAX_IO*2);
          if(finished_requests_id[request_pos]==*request_id){
            break;
          }
          i++;
        }
        if(i==nb_finished_requests){
          return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_test_request_th (2))\n");
        }      
        *flag=1;
      }
    }
  }
  /*  printf("je vais essayer de nettoyer %d\n",nb_finished_requests);
    pthread_mutex_unlock(&io_mutex);
    sleep (10); */
  mumps_owns_mutex=1;
  /* printf("ici nb_finished_requests=%d\n",nb_finished_requests);*/
  mumps_clean_finished_queue_th();
  mumps_owns_mutex=0;
  pthread_mutex_unlock(&io_mutex);
  return 0;
}
int mumps_wait_req_sem_th(int *request_id){
  int i,j;
  j=first_active;
  for(i=0;i<nb_active;i++){
    if(io_queue[j].req_num==*request_id) break;
    j=(j+1)%MAX_IO;
  }
  if(i<nb_active){
    mumps_wait_sem(&(io_queue[j].int_local_cond),&(io_queue[j].local_cond));
  }
  return 0;
}
int mumps_wait_request_th(int *request_id){
  /* waits for the termination of the request "request_id" */
  int flag=0,ierr;
  if(with_sem!=2){
    while(!flag){
      ierr=mumps_test_request_th(request_id,&flag);
      if(ierr!=0)return ierr;
    }
  }else{
    ierr=mumps_test_request_th(request_id,&flag);
    if(ierr!=0)return ierr;
    if(!flag){
      mumps_wait_req_sem_th(request_id);
      ierr=mumps_test_request_th(request_id,&flag);
      if(ierr!=0)return ierr;
    }
  }
  return 0;
}
int mumps_is_there_finished_request_th(int* flag){
  if(!mumps_owns_mutex) pthread_mutex_lock(&io_mutex);
  if(nb_finished_requests==0){
    *flag=0;
  }else{
    *flag=1;
    /*    printf("finished : %d\n",nb_finished_requests);     */
  }
  if(!mumps_owns_mutex) pthread_mutex_unlock(&io_mutex);
  return 0;
}
int mumps_clean_finished_queue_th(){
   /* Cleans the finished request queue. On exit, the queue is empty.*/
   int local_flag;
   int cur_req;
   int loc_owned_mutex=0,ierr;
   if(!mumps_owns_mutex){
      pthread_mutex_lock(&io_mutex);
      mumps_owns_mutex=1;
      loc_owned_mutex=1;
  }
  /* this block of code is designed for avoiding deadlocks between
     the two threads*/
   mumps_is_there_finished_request_th(&local_flag);
   while(local_flag){
     ierr=mumps_clean_request_th(&cur_req);
     if(ierr!=0){
       return ierr;
     }
     mumps_is_there_finished_request_th(&local_flag);
   }
   if((!mumps_owns_mutex)||(loc_owned_mutex)){
      pthread_mutex_unlock(&io_mutex);
      mumps_owns_mutex=0;
   }
   return 0;
}
int mumps_clean_request_th(int* request_id){
  int ierr;
  ierr=mumps_check_error_th();
  if(ierr!=0){
    return ierr;
  }
  if(!mumps_owns_mutex)pthread_mutex_lock(&io_mutex);
  *request_id=finished_requests_id[first_finished_requests];
  if(smallest_request_id!=finished_requests_id[first_finished_requests]){
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_clean_request_th)\n");
  }
  finished_requests_id[first_finished_requests]=-9999;
  first_finished_requests=(first_finished_requests+1)%(MAX_FINISH_REQ);
  nb_finished_requests--;
  /*we treat the io requests in their arrival order => we just have to
    increase smallest_request_id*/
  smallest_request_id++;
  if(!mumps_owns_mutex) pthread_mutex_unlock(&io_mutex);
  if(with_sem) {
      if(with_sem==2){
        mumps_post_sem(&int_sem_nb_free_finished_requests,&cond_nb_free_finished_requests);
      }
  }
  return 0;
}
int mumps_low_level_init_ooc_c_th(int* async, int* ierr){
  int i, ret_code;    
  char buf[64];
  /* Computes the number of files needed. Uses ceil value. */
  *ierr=0;
  current_req_num=0;
  with_sem=2;
  first_active=0;
  last_active=0;
  nb_active=0;
  first_finished_requests=0;
  last_finished_requests=0;
  nb_finished_requests=0;
  smallest_request_id=0;
  mumps_owns_mutex=0;
  inactive_time_io_thread=0;
  time_flag_io_thread=0;
  gettimeofday(&origin_time_io_thread,NULL);
  /*  mumps_io_flag_async=*async; */
  if(*async!=IO_ASYNC_TH){
    *ierr = -91;
    sprintf(buf,"Internal error: mumps_low_level_init_ooc_c_th should not to be called with strat_IO=%d\n",*async);
    return mumps_io_error(*ierr,buf);
  }
  if(*async){
    pthread_mutex_init(&io_mutex,NULL);
    mumps_io_init_err_lock();
#ifdef WITH_PFUNC
    mumps_io_init_pointers_lock();
#endif
    io_queue=(struct request_io *)malloc(MAX_IO*sizeof(struct request_io));
    if(with_sem==2){
      for(i=0;i<MAX_IO;i++){
        pthread_cond_init(&(io_queue[i].local_cond),NULL);
        io_queue[i].int_local_cond=0;
      }
    }
    finished_requests_id=(int *)malloc(MAX_IO*2*sizeof(int));
    finished_requests_inode=(int *)malloc(MAX_IO*2*sizeof(int));
    for(i=0;i<MAX_IO*2;i++){
      finished_requests_id[i]=-9999;
      finished_requests_inode[i]=-9999;
    }
    if(with_sem){
      switch(with_sem){
      case 2:
        int_sem_io=0;
        int_sem_stop=0;
        int_sem_nb_free_finished_requests=MAX_FINISH_REQ;
        int_sem_nb_free_active_requests=MAX_IO;
        pthread_cond_init(&cond_stop,NULL);
        pthread_cond_init(&cond_io,NULL);
        pthread_cond_init(&cond_nb_free_active_requests,NULL);
        pthread_cond_init(&cond_nb_free_finished_requests,NULL);
        pthread_mutex_init(&io_mutex_cond,NULL);
        break;
      default:
        *ierr = -92;
        sprintf(buf,"Internal error: mumps_low_level_init_ooc_c_th should not to be called with strat_IO=%d\n",*async);
        return mumps_io_error(*ierr,buf);
      }
      ret_code=pthread_create(&io_thread,NULL,mumps_async_thread_function_with_sem,NULL);
    }
    if(ret_code!=0){
      errno = ret_code;
      return mumps_io_sys_error(-92,"Unable to create I/O thread");
    }
    main_thread=pthread_self();
  }
  return 0;
}
int mumps_async_write_th(const int * strat_IO, 
                        void * address_block,
                        long long block_size,
                        int * inode,
                        int * request_arg,
                        int * type,
                        long long vaddr,
                        int * ierr){
  int cur_req;
  *ierr=mumps_check_error_th();
  if(*ierr!=0){
    return *ierr;
  }
  if(with_sem){ 
    mumps_clean_finished_queue_th();
      if(with_sem==2){
        mumps_wait_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
      }
    /*    sem_wait(&sem_nb_free_active_requests); */
    pthread_mutex_lock(&io_mutex);
  }
  if(nb_active<=MAX_IO){
    if(nb_active==0){
      first_active=last_active;
    }
    else{
      last_active=(last_active+1)%MAX_IO;
    }
    cur_req=last_active;
    nb_active++;
    io_queue[cur_req].inode=*inode; 
    io_queue[cur_req].req_num=current_req_num; 
    io_queue[cur_req].addr=address_block; 
    io_queue[cur_req].size=block_size;  
    io_queue[cur_req].vaddr=vaddr;
    io_queue[cur_req].io_type=0;
    io_queue[cur_req].file_type=*type;
    if(with_sem==2){
      io_queue[cur_req].int_local_cond=0;      
    }
    *request_arg=current_req_num; 
    current_req_num++;
  }else{
    *ierr = -91;
    return mumps_io_error(*ierr,"Internal error in OOC Management layer (mumps_async_write_th)\n");
    /*    exit(-3);*/
  }
  pthread_mutex_unlock(&io_mutex);
  if(with_sem){
    /*    sem_post(&sem_io); */
      if(with_sem==2){
        mumps_post_sem(&int_sem_io,&cond_io);
      }
  }
  return 0;
}
int mumps_async_read_th(const int * strat_IO, 
                       void * address_block,
                       long long  block_size,
                       int * inode,
                       int * request_arg,
                        int * type,
                       long long vaddr,
                       int * ierr){
  int cur_req;  
  *ierr=mumps_check_error_th();
  if(*ierr!=0){
    return *ierr;
  }
  if(with_sem){
    mumps_clean_finished_queue_th();
    /* end of the block*/
      if(with_sem==2){
        mumps_wait_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
      }
    /*    sem_wait(&sem_nb_free_active_requests); */
    pthread_mutex_lock(&io_mutex);
  }
  if(nb_active<MAX_IO){
    if(nb_active==0){
      first_active=last_active;
    }else{
      last_active=(last_active+1)%MAX_IO;
      /*        if(last_active<MAX_IO-1){
                cur_req=last_active+1;
                last_active++;
                }else{
                cur_req=0;
                last_active=0;
                }*/
    }
    cur_req=last_active;
    nb_active++;
    io_queue[cur_req].inode=*inode;
    io_queue[cur_req].req_num=current_req_num; 
    io_queue[cur_req].addr=address_block; 
    io_queue[cur_req].size=block_size;  
    io_queue[cur_req].vaddr=vaddr;
    io_queue[cur_req].io_type=1;
    io_queue[cur_req].file_type=*type;
    if(with_sem==2){
      io_queue[cur_req].int_local_cond=0;      
    }
    *request_arg=current_req_num; 
    current_req_num++;
  }else{
    *ierr = -91;
    return mumps_io_error(*ierr,"Internal error in OOC Management layer (mumps_async_read_th)\n");
  }
  if(with_sem){
    /*    sem_post(&sem_io); */
      if(with_sem==2){
        mumps_post_sem(&int_sem_io,&cond_io);
      }
  }
  pthread_mutex_unlock(&io_mutex);
  return 0;
}
int mumps_clean_io_data_c_th(int *myid){
  int i;
  /* cleans the thread/io management data*/
  if(mumps_io_flag_async){
    /*we can use either signals or mutexes for this step */
    if(with_sem){
        if(with_sem==2){
          mumps_post_sem(&int_sem_stop,&cond_stop);
          mumps_post_sem(&int_sem_io,&cond_io);
        }
    }else{
      pthread_mutex_lock(&io_mutex);
      io_flag_stop=1;
      pthread_mutex_unlock(&io_mutex);
    }
    pthread_join(io_thread,NULL);
    pthread_mutex_destroy(&io_mutex);
    mumps_io_destroy_err_lock();
#ifdef WITH_PFUNC
    mumps_io_destroy_pointers_lock();
#endif
    if(with_sem){
        if(with_sem==2){
          pthread_cond_destroy(&cond_stop);
          pthread_cond_destroy(&cond_io);
          pthread_cond_destroy(&cond_nb_free_active_requests);
          pthread_cond_destroy(&cond_nb_free_finished_requests);
          pthread_mutex_destroy(&io_mutex_cond);
        }
    }
  }
  if(with_sem==2){
    for(i=0;i<MAX_IO;i++){
      pthread_cond_destroy(&(io_queue[i].local_cond));
    }
  }
  free(io_queue);
  free(finished_requests_id);
  free(finished_requests_inode);
  return 0;
}
int mumps_get_sem(void *arg,int *value){
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    *value=*((int *)arg);
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps__get_sem)\n");
  }
  return 0;
}
int mumps_wait_sem(void *arg,pthread_cond_t *cond){
  int *tmp_pointer;
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    tmp_pointer=(int *)arg;
    while(*tmp_pointer==0){
      pthread_cond_wait(cond, &io_mutex_cond);
    }
    (*tmp_pointer)--;
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_wait_sem)\n");
  }
  return 0;  
}
int mumps_post_sem(void *arg,pthread_cond_t *cond){
  int *tmp_pointer;
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    tmp_pointer=(int *)arg;
    (*tmp_pointer)++;
    if(*tmp_pointer==1){
      pthread_cond_broadcast(cond);
    }
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_post_sem)\n");
  }
  return 0;  
}
#endif /* MUMPS_WIN32 && WITHOUT_PTHREAD */
