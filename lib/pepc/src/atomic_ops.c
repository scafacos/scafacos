/*
* This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
* 
* Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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
>  atomic operations library
>
*************************************************************************/

#include <stdlib.h>

#include "opa_primitives.h"

OPA_int_t* _atomic_alloc_int()
{
  void* ret;
  if (0 != posix_memalign(&ret, (size_t)8, sizeof(OPA_int_t)))
  {
    ret = NULL;
  }
  
  return (OPA_int_t*)ret;
}

void _atomic_free_int(OPA_int_t* storage)
{
  free(storage);
}

int _atomic_load_int(OPA_int_t* storage)
{
  return OPA_load_int(storage);
}

void _atomic_store_int(OPA_int_t* storage, int val)
{
  OPA_store_int(storage, val);
}

int _atomic_fetch_and_increment_int(OPA_int_t* storage)
{
  return OPA_fetch_and_incr_int(storage);
}

int _atomic_mod_increment_and_fetch_int(OPA_int_t* storage, int mod)
{
  int cmp;
  int prev = OPA_load_int(storage);

  do {
      cmp = prev;
      prev = OPA_cas_int(storage, cmp, prev % mod + 1);
  } while (prev != cmp);

  return prev % mod + 1;
}

int _atomic_compare_and_swap_int(OPA_int_t* storage, int oldval, int newval)
{
  return OPA_cas_int(storage, oldval, newval);
}

void _atomic_write_barrier()
{
  OPA_write_barrier();
}

void _atomic_read_barrier()
{
  OPA_read_barrier();
}

void _atomic_read_write_barrier()
{
  OPA_read_write_barrier();
}

void _critical_section_enter(OPA_int_t* storage)
{
  do {
    /* wait until it is open and immediately block it */
  } while(0 != OPA_cas_int(storage, 0, 1));
}

void _critical_section_leave(OPA_int_t* storage)
{
  /* open it again */
  OPA_store_int(storage, 0);
}
