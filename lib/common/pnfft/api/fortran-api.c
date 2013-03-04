/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 * Copyright (c) 2010-2013 Michael Pippig
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

#include <stddef.h> /* for ptrdiff_t */

#include <pnfft.h>
#include <inpfft.h>
#include <util.h>
#include "fortran-mangling.h"

/* if FC_FUNC is not defined, then we don't know how to mangle identifiers
   for the Fortran linker, and we must omit the Fortran API. */
#if defined(FC_FUNC) || defined(WINDOWS_FC_MANGLING)

/*-----------------------------------------------------------------------*/
/* some internal functions used by the Fortran api */

/* in Fortran, the natural array ordering is column-major, which
   corresponds to reversing the dimensions relative to C's row-major */
// static void revert_int(int rnk, const int *n, int *n_rev)
// {
//   for (int i = 0; i < rnk; i++)
//     n_rev[rnk-1-i] = n[i];
// }

static void revert_int(int rnk, const int *n, int *n_rev)
{
  for (int i = 0; i < rnk; i++)
    n_rev[rnk-1-i] = n[i];
}

static void revert_INT(int rnk, const INT *n, INT *n_rev)
{
  for (int i = 0; i < rnk; i++)
    n_rev[rnk-1-i] = n[i];
}

static void revert_R(int rnk, const R *n, R *n_rev)
{
  for (int i = 0; i < rnk; i++)
    n_rev[rnk-1-i] = n[i];
}

static void revert_and_sub_ones_INT(int rnk, const INT *n, INT *n_rev)
{
  revert_INT(rnk, n, n_rev);
  for (int i = 0; i < rnk; i++)
    n_rev[i]--;
}

static void revert_and_add_ones_INT(int rnk, const INT *n, INT *n_rev)
{
  revert_INT(rnk, n, n_rev);
  for (int i = 0; i < rnk; i++)
    n_rev[i]++;
}

static int *malloc_and_revert_int(int rnk, const int *n)
{
  int *n_rev;
  
  if(n == NULL)
    return NULL;
  
  n_rev = PNX(malloc_int)(rnk);
  revert_int(rnk, n, n_rev);
  
  return n_rev;
}

static INT *malloc_and_revert_INT(int rnk, const INT *n)
{
  INT *n_rev;
  
  if(n == NULL)
    return NULL;
  
  n_rev = PNX(malloc_INT)(rnk);
  revert_INT(rnk, n, n_rev);
  
  return n_rev;
}

static INT *malloc_and_revert_R(int rnk, const R *n)
{
  R *n_rev;
  
  if(n == NULL)
    return NULL;
  
  n_rev = PNX(malloc_R)(rnk);
  revert_R(rnk, n, n_rev);
  
  return n_rev;
}




/*-----------------------------------------------------------------------*/

#  define FORT(a, A) FCpx(pxf(a), PXF(A))

#  ifndef WINDOWS_FC_MANGLING

#    if defined(FC_FUNC)
#      define FCpx(a, A) FC_FUNC(a, A)
#      include "fortran-wrappers.h"
#    endif

/* If identifiers with underscores are mangled differently than those
   without underscores, then we include *both* mangling versions.  The
   reason is that the only Fortran compiler that does such differing
   mangling is currently g77 (which adds an extra underscore to names
   with underscores), whereas other compilers running on the same
   machine are likely to use non-underscored mangling.  (I'm sick
   of users complaining that FFTW works with g77 but not with e.g.
   pgf77 or ifc on the same machine.)  Note that all FFTW identifiers
   contain underscores, and configure picks g77 by default. */
#    if defined(FC_FUNC_) && !defined(FC_FUNC_EQUIV)
#      undef FCpx
#      define FCpx(a, A) FC_FUNC_(a, A)
#      include "fortran-wrappers.h"
#    endif

#  else /* WINDOWS_FC_MANGLING */

/* Various mangling conventions common (?) under Windows. */

/* g77 */
#    define WINDOWS_FC_FUNC(a, A) a ## __
#    define FCpx(a, A) WINDOWS_FC_FUNC(a, A)
#    include "fortran-wrappers.h"

/* Intel, etc. */
#    undef WINDOWS_FC_FUNC
#    define WINDOWS_FC_FUNC(a, A) a ## _
#    include "fortran-wrappers.h"

/* Digital/Compaq/HP Visual Fortran, Intel Fortran.  stdcall attribute
   is apparently required to adjust for calling conventions (callee
   pops stack in stdcall).  See also:
       http://msdn.microsoft.com/library/en-us/vccore98/html/_core_mixed.2d.language_programming.3a_.overview.asp
*/
#    undef WINDOWS_FC_FUNC
#    if defined(__GNUC__)
#      define WINDOWS_FC_FUNC(a, A) __attribute__((stdcall)) A
#    elif defined(_MSC_VER) || defined(_ICC) || defined(_STDCALL_SUPPORTED)
#      define WINDOWS_FC_FUNC(a, A) __stdcall A
#    else
#      define WINDOWS_FC_FUNC(a, A) A /* oh well */
#    endif
#    include "fortran-wrappers.h"

#  endif /* WINDOWS_FC_MANGLING */

#endif /* FC_FUNC || WINDOWS_FC_MANGLING */
