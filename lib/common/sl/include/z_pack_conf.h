/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */


#ifndef __Z_PACK_CONF_H__
#define __Z_PACK_CONF_H__


#include <string.h>

#include "sl_rename.h"

#include "sl_config.h"
#include "sl_config_intern.h"

#ifdef SL_USE_MPI
 #define SL_PROC_RANK  SL_DEFCON(mpi.rank)
#else
 #define SL_PROC_RANK  SL_DEFCON(sl.dummy_rank)
#endif

#include "sl_environment.h"
#include "sl_environment_intern.h"

#include "sl_elements.h"

#include "sl_pelem.h"

#include "sl_rti_tids.h"
#include "sl_rti_intern.h"

#include "spec_public_conf.h"
#include "spec_public.h"

#include "sl_types.h"
#include "sl_types_intern.h"

#include "sl_globals.h"

#include "sl_adds.h"
#include "sl_adds_intern.h"

#ifdef SL_PREFIX
# define Z_PACK_RENAME
# define Z_PREFIX  SL_PREFIX
#endif

typedef sl_int_type_c z_int_t;
#define z_int_fmt  sl_int_type_fmt

#ifdef SL_USE_MPI
# define Z_PACK_MPI
#endif
#define Z_PACK_MPI_RANK  SL_PROC_RANK

#define Z_PACK_NUMERIC

#define Z_PACK_DEBUG

#ifdef SLDEBUG_OUTPUT
# define Z_DEBUG_LEVEL  SLDEBUG_OUTPUT
#endif

#define Z_PACK_ALLOC

#ifdef SLDEBUG_ALLOC
# define Z_ALLOC_DEBUG  SLDEBUG_ALLOC
#endif

#define Z_PACK_TIME

#if defined(SL_TIMING) && (SL_TIMING >= 2)
# define Z_PACK_TIMING
#endif
#ifdef SL_TIMING_PRINT_PREFIX
# define Z_TIMING_PRINT_PREFIX  SL_TIMING_PRINT_PREFIX
#endif
#ifdef SL_TIMING_PRINT
# define Z_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)  SL_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)
#endif

#define Z_PACK_RANDOM

#define Z_PACK_DIGEST

#define Z_PACK_CRC32

typedef unsigned int z_crc32_t;
#define z_crc32_fmt  "X"

#define Z_PACK_GMP


#endif /* __Z_PACK_CONF_H__ */
