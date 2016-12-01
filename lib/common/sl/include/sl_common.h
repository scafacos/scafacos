/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
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


#ifndef __SL_COMMON_H__
#define __SL_COMMON_H__


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>


#ifdef SLDEBUG
# ifndef SLDEBUG_OUTPUT
#  define SLDEBUG_OUTPUT  SLDEBUG
# endif
# ifndef SLDEBUG_ALLOC
#  define SLDEBUG_ALLOC   SLDEBUG
# endif
#endif

#ifdef SLDEBUG_OUTPUT_NOT
# undef SLDEBUG_OUTPUT
#endif

#ifdef SLDEBUG_ALLOC_NOT
# undef SLDEBUG_ALLOC
#endif

#ifndef SL_TIMING_PRINT_DEFAULT
# define SL_TIMING_PRINT_DEFAULT  z_timing_print_default
#endif


#include "sl_rename.h"

#ifndef __Z_PACK_H__
# define NO__Z_PACK_H__
#endif

#include "sl_config.h"
#include "sl_config_intern.h"

#ifdef SL_USE_MPI
 #include <mpi.h>
 #define SL_PROC_RANK  SL_DEFCON(mpi.rank)
#else
 #define SL_PROC_RANK  SL_DEFCON(sl.dummy_rank)
#endif

#ifdef SL_USE_OMP
# ifdef _OPENMP
#  include <omp.h>
# else
#  error _OPENMP not defined!
# endif
#endif

#include "z_config.h"

#include "sl_tune.h"
#include "sl_tune_auto.h"
#include "sl_tune_intern.h"

#include "sl_deprecated.h"

#include "sl_environment.h"
#include "sl_environment_intern.h"

#ifdef SL_EXTRA
# include "sl_extra.h"
#endif

#if defined(NO__Z_PACK_H__) && defined(__Z_PACK_H__)
# error "z_pack.h" is not allowed to be used in config/tune/environment (the sorting library already uses z-pack!)
#endif

#include "sl_elements.h"

#include "sl_pelem.h"

#include "sl_rti_tids.h"
#include "sl_rti_intern.h"

#include "spec_public_conf.h"
#include "spec_public.h"

#include "sl_types.h"
#include "sl_types_intern.h"

#include "sl_macros.h"

#include "sl_globals.h"

#include "sl_adds.h"
#include "sl_adds_intern.h"


#define Z_IGNORE_CONFIG_H    /* do not include <config.h> */
#define Z_IGNORE_Z_CONFIG_H  /* do not include "z_config.h" */
#include "z_pack.h"

#define SL_PROTO(_f_)  _f_
#include "sl_protos.h"
#include "sl_protos_intern.h"
#undef SL_PROTO
#undef __SL_PROTOS_H__
#undef __SL_PROTOS_INTERN_H__
#define SL_PROTO(_f_)  _f_##_di
#include "sl_protos.h"
#include "sl_protos_intern.h"
#undef SL_PROTO

#ifdef SL_USE_MPI
# define SL_PROTO(_f_)  _f_
#  include "sl_protos_mpi.h"
# undef SL_PROTO
# undef __SL_MPI_PROTOS_H__
# define SL_PROTO(_f_)  _f_##_di
#  include "sl_protos_mpi.h"
# undef SL_PROTO
#endif


#endif /* __SL_COMMON_H__ */
