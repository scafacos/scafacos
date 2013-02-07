/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <dcmfdimpl.h>

_BGP_Atomic global_atomic = _BGP_ATOMIC_INIT(0);
volatile uint32_t global_lock __attribute__((__aligned__(16)));

static __inline__ uint32_t testandset(volatile uint32_t *global_lock)
{
  uint32_t ret;
  uint32_t val = 1;  

  __asm__ volatile(
                   "loop:	lwarx 	%0,0,%1   \n"
                   "            stwcx.	%2,0,%1   \n"
                   "            bne-    loop      "
                   : "=r" (ret)
                   : "r" (global_lock)
                   : "r" (val)
                  );  
   return ret;
}

static __inline__ void reset(volatile uint32_t *global_lock)
{
  uint32_t val = 0;

  __asm__ volatile(
                   "            mr      %0,%1"
                   : "=r" (global_lock)
                   : "r" (val)
                  );
   return;
}

void A1DI_Global_lock_acquire()
{
   while(testandset(&global_lock)); 
   return;
}

void A1DI_Global_lock_release()
{
   reset(&global_lock);
   return;
}
