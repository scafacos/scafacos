/*
 * Copyright (c) 2011-2013 Michael Pippig
 *
 * This file is part of PNFFT.
 *
 * PNFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PNFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PNFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "pnfft.h"
#include "ipnfft.h"


void *PNX(malloc)(size_t n){
  return PX(malloc)(n);  
}

R *PNX(alloc_real)(size_t n){
  return PX(alloc_real)(n);
}

C *PNX(alloc_complex)(size_t n){
  return PX(alloc_complex)(n);
}
  
void PNX(free)(void *p){
  PX(free)(p);
}



void PNX(die)(
    const char *s, MPI_Comm comm
    )
{
  fflush(stdout);
  PX(fprintf)(comm, stderr, s);
  exit(EXIT_FAILURE);
}
