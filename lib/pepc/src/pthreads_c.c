/*
* This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
*
* Copyright (C) 2002-2014 Juelich Supercomputing Centre,
*                         Forschungszentrum Juelich GmbH,
*                         Germany
*
* PEPC is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* PEPC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
*/


/*************************************************************************
>
>  pthreads interface library
>
*************************************************************************/

#ifndef __APPLE__
  #define _GNU_SOURCE
  #include <features.h>
  #include <sys/prctl.h>
#endif

#if defined(__TOS_BGQ__)
#define _GNU_SOURCE
#include <stdint.h>
#include <spi/include/kernel/location.h>
#include <spi/include/kernel/thread.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <pthread.h>

#ifdef TAU_PAPI
#include <stdbool.h>
#endif

pthread_rwlock_t *my_rwlocks;
pthread_attr_t thread_attr;

const int THREAD_TYPE_DEFAULT = 0;
const int THREAD_TYPE_MAIN = 1;
const int THREAD_TYPE_COMMUNICATOR = 2;
const int THREAD_TYPE_WORKER = 3;

const char threadnames[4][20] = {"Tdefault", "Tmain", "Tcomm", "Twork"};

typedef struct {
  pthread_t* thread;
  void* (*start_routine)(void*);
  void* arg;
  int thread_type;
  int counter;
} pthread_with_type_t;

#define CHECKRES do {if (iret != 0) return iret;} while(0);

//////////////// PThreads //////////////////////

int pthreads_init()
{
    int iret = 0;

    iret = pthread_attr_init(&thread_attr);
    CHECKRES;

    iret = pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
    CHECKRES;

    iret = pthread_attr_setscope(&thread_attr, PTHREAD_SCOPE_SYSTEM);
    CHECKRES;

    return 0;
}


int pthreads_uninit()
{
    int iret = 0;

    iret = pthread_attr_destroy(&thread_attr);
    CHECKRES;

    return 0;
}


#if defined(__TOS_BGQ__)
  #if defined(__GNUC__)
    #define cntlz8(x) (__builtin_ctzl(x))
    #define popcnt8(x) (__builtin_popcountl(x))
  #elif defined(__IBMC__) || defined(__IBMCPP__)
    #define cntlz8(x) (__cntlz8(x))
    #define popcnt8(x) (__popcnt8(x))
  #endif


void place_thread(int thread_type, int counter)
{
    const int reserved = 4;
    pthread_t tid = pthread_self();
    cpu_set_t cpumask;
    // accessible is a 64-bit mask identifying all hardware threads
    // allocated to this process. Most significant bit corresponds
    // to hardware thread 0.
    uint64_t accessible = Kernel_ThreadMask(Kernel_MyTcoord());
    unsigned int first;
    int count, selected = -1;

    if (thread_type == THREAD_TYPE_MAIN) {
      first = cntlz8(accessible); // identify first accessible hardware thread
      selected = first;
    } else if (thread_type == THREAD_TYPE_COMMUNICATOR) {
      first = cntlz8(accessible); // identify first accessible hardware thread

      if (counter == 1) {
        // first communicator goes on the free hardware thread of the first core
        selected = first + 2;
      } else if (counter == 2) {
        // second communicator shares a hardware thread with the main thread
        selected = first;
      } else {
        printf("WARNING: No placement policy for communicator thread no. %d implemented!\n", counter);
      }
    } else if (thread_type == THREAD_TYPE_WORKER) {
      first = cntlz8(accessible); // identify first accessible hardware thread
      count = popcnt8(accessible) - reserved; // number of accessible threads minus reserved first core

      if (count > 0) {
        selected = first + reserved  // skip first core
          // put the first count / 2 workers on even numbered hardware threads
          // (these are not shared with MPI threads)
          + ((2 * (counter - 1)) % count)
          // put additional threads on odd numbered hardware threads
          + (((2 * (counter - 1)) / count) % 2);
          // strictly speaking, this will (intentionally) wrap around to even
          // numbered hardware threads again after placing count workers.

      } else {
        // only one core is available, put workers on odd numbered hardware threads,
        // away from main and communicator threads on the even numbered cores
        selected = first + 1 + 2 * ((counter - 1) % 2);
      }
    }

    if (selected > -1) {
      // settle down on a hardware thread according to policy
      CPU_ZERO(&cpumask);
      CPU_SET(selected, &cpumask);
      pthread_setaffinity_np(tid, sizeof(cpumask), &cpumask);
    }
}
#else
void place_thread(int thread_type, int counter) { }
#endif

void rename_thread(pthread_with_type_t* thread)
{
    char threadname[16];

    sprintf(threadname, "%s [%d]", threadnames[thread->thread_type], thread->counter);
    // see http://stackoverflow.com/questions/2369738/can-i-set-the-name-of-a-thread-in-pthreads-linux for details
    #ifdef __APPLE__
    pthread_setname_np(threadname);
    #else
      #if ((__GLIBC__ == 2) && (__GLIBC_MINOR__ >= 12)) || (__GLIBC__ > 2)
        // try the recommended way first
        if (0!=pthread_setname_np(*thread->thread, threadname))
      #endif // __GLIBC >=2.12
      // then with some more brute force
      prctl(PR_SET_NAME,threadname,0,0,0);
    #endif // __APPLE__
}


#ifdef TAU_PAPI
typedef struct {
  void* (*start_routine)(void*);
  void* arg;
} tau_pthread_pack;

extern pthread_key_t wrapper_flags_key;

extern void* tau_pthread_function(void*);

void* thread_helper(pthread_with_type_t* thread)
{
    tau_pthread_pack pack;
    pack.start_routine = thread->start_routine;
    pack.arg = thread->arg;

    rename_thread(thread);
    place_thread(thread->thread_type, thread->counter);

    return tau_pthread_function((void*)&pack);
}


int pthreads_createthread_c(pthread_with_type_t* thread)
{
    int ret;
    bool* wrapped;
    thread->thread = (pthread_t*) malloc(sizeof (pthread_t));

    wrapped = pthread_getspecific(wrapper_flags_key);
    if (!wrapped) {
      wrapped = (bool*) malloc(sizeof(bool));
      pthread_setspecific(wrapper_flags_key, (void*) wrapped);
    }

    *wrapped = true;
    ret = pthread_create(thread->thread, &thread_attr, (void* (*)(void*)) thread_helper, thread);
    *wrapped = false;
    return ret;
}

#else

void* thread_helper(pthread_with_type_t* thread)
{
    rename_thread(thread);
    place_thread(thread->thread_type, thread->counter);

    return (thread->start_routine)(thread->arg);
}


int pthreads_createthread_c(pthread_with_type_t* thread)
{
    thread->thread = (pthread_t*) malloc(sizeof (pthread_t));

    return pthread_create(thread->thread, &thread_attr, (void* (*)(void*)) thread_helper, thread);
}
#endif


int pthreads_jointhread(pthread_with_type_t thread)
{
    void *retval; // for convenience we do not pass it to fortran
    int ret;

    ret = pthread_join(*(thread.thread), &retval);
    free(thread.thread);
    return ret;
}


int pthreads_exitthread()
{
    pthread_exit(NULL);
    return 0;
}


int pthreads_sched_yield()
{
    return sched_yield();
}

///////////////// Utils //////////////////////

int get_my_tid()
{
    return (int)syscall(SYS_gettid);
}

int get_my_pid()
{
    return (int)getpid();
}
